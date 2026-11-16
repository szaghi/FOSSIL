!< FOSSIL, field transfer between the wall and its sample points (#109 WP2 step 04).
module fossil_field_transfer
!< FOSSIL, field transfer between the wall and its sample points (#109 WP2 step 04).
!<
!< A solver keeps fields on the wall (a wall temperature, a wall velocity, a deformation, per vertex or per facet) and
!< needs their values where its closure and its loads look at the wall; and it computes fields at those points (the wall
!< stress at the load samples) that it wants back on the wall. Two directions:
!<
!< **Wall to points** (`transfer_field_3d` / `transfer_field_2d`): for points on (or projected onto) an element given by
!< its id (3-D a facet, 2-D a segment = its first vertex) the value of
!<
!< - a per-vertex field `vertex_field(ncomp, nv)` (3-D: the vertices of the surface's vertex pool, `facet_at(f)%vertex_id`,
!<   nv its size; 2-D: the curve's vertices) interpolated linearly on the element with the weights of the point's closest
!<   point on it: 3-D the barycentric coordinates (1 - s - t, s, t) of `triangle_closest_st` (the solve of every distance
!<   query; a point equal to a vertex takes that vertex's value exactly), 2-D (1 - s, s) of the clamped projection on the
!<   segment. A linear field is reproduced at any point of the element to round-off;
!< - a per-element field `facet_field(ncomp, nf)` (2-D `segment_field(ncomp, nv)`): piecewise constant, copied.
!<
!< ncomp components (a scalar ncomp = 1, a vector 2 or 3). An id outside 1..nf (0: an ILW slot without a crossing, a ghost
!< without a foot) gives zeros. `pose=`: the points are in the world, the body at the pose (the frame of a posed closure
!< table): each point is taken to the body frame (`to_body`) before its weights. Wrappers: `transfer_to_closure` (the feet
!< of an `ib_closure_t` and, optionally, its ILW crossings: the ghost's coordinates with the slot's axis set to `ilw_t`, on
!< `ilw_element`) and `transfer_to_samples` (the points of `load_samples_t`, interpolated with the rule's own barycentric
!< coordinates, not the closest point: no pose needed).
!<
!< **Points to wall** (`samples_to_vertices`): values at the load samples' quadrature points to a per-vertex field, the
!< projection a solver needs for a wall stress known at its samples (van Noordt et al. 2022 average theirs to the
!< vertices). Per element, a linear fit of its q samples in its barycentric coordinates (weighted least squares with the
!< rule's weights; order 1: the constant; `linear=.false.`: the element's weighted mean), evaluated at each of its
!< vertices; per vertex, the average of those values over its incident elements weighted by their areas (2-D lengths:
!< the sums of their samples' weights), in increasing element order with Neumaier sums. A smoothing projection: a
!< constant is reproduced exactly (the sums are of deviations from the first incident value), values linear on every
!< element to round-off at every vertex from order 2 (`linear=.true.`): any field linear in space, on any mesh, and any
!< vertex field brought to the samples by `transfer_to_samples` (the projection is its left inverse); anything else is
!< fitted per element and averaged, and noise at the samples is smoothed over the vertex's ring. A vertex without an
!< incident element of positive area gets zeros.
!<
!< Two paths, the same values bit for bit: the host and the device (`device=.true.`), one set of loops with `if(dev)` and
!< the shared kernels (`transfer_weights_3d_oac`, `transfer_weights_2d_oac`, the rule kernels of `fossil_surface_loads`,
!< `to_body`). The inputs are read on the device where they are mapped (a resident closure table or load samples, a field
!< mapped by the caller with `!$acc enter data`), else copied in; the output is written on the device where the caller
!< mapped it (then `!$acc update host` to read it), else copied out. The host path reads and writes the host arrays
!< (`update_host` resident tables first). The facet vertex ids (3-D) and the curve's vertices (2-D) are gathered on the
!< host per call and copied in. Built without OpenACC the device loops run on the host in order.
use fossil_curve_object,       only : curve_object
use fossil_fd_closure,         only : ib_closure_t, to_body_3d_oac, to_body_2d_oac
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t
use fossil_surface_loads,      only : load_samples_t, rule_3d_oac, load_point_2d_oac
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use fossil_triangle_closest,   only : triangle_closest_st
use fossil_vertex_pool_object, only : vertex_pool_object
use iso_fortran_env,           only : stderr => error_unit
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
private
public :: transfer_field_3d, transfer_field_2d, transfer_to_closure, transfer_to_samples, samples_to_vertices
! the per-point kernels, shared by the two paths (public for the tests)
public :: transfer_weights_3d_oac, transfer_weights_2d_oac

interface transfer_to_closure
   !< A wall field at the feet (and the ILW crossings) of a closure table: 3-D (surface) or 2-D (curve).
   module procedure transfer_to_closure_3d, transfer_to_closure_2d
endinterface transfer_to_closure

interface transfer_to_samples
   !< A wall field at the points of load samples: 3-D (surface) or 2-D (curve).
   module procedure transfer_to_samples_3d, transfer_to_samples_2d
endinterface transfer_to_samples

interface samples_to_vertices
   !< Values at the load samples to a per-vertex field: 3-D (surface) or 2-D (curve).
   module procedure samples_to_vertices_3d, samples_to_vertices_2d
endinterface samples_to_vertices

integer(I4P), parameter :: MODE_FACET  = 1_I4P !< Per-element field, copied.
integer(I4P), parameter :: MODE_VERTEX = 2_I4P !< Per-vertex field, closest-point weights.
integer(I4P), parameter :: MODE_RULE   = 3_I4P !< Per-vertex field, the load rule's barycentric coordinates.

contains
   ! wall to points

   subroutine transfer_field_3d(surface, facets, points, values, status, vertex_field, facet_field, pose, device)
   !< The values of a wall field at points on (or projected onto) the facets `facets(n)` (see the module header):
   !< `points(3, n)`, `values(ncomp, n)`; exactly one of `vertex_field(ncomp, nv)` (nv = the vertex pool's size) and
   !< `facet_field(ncomp, nf)`. `pose`: the points in the world, the body at the pose. `device` (default .false.): the
   !< device path. `status = STATUS_INVALID_INPUT` (stderr message, values zero) on no field or both, a size that differs,
   !< or a surface with facets not analyzed.
   type(surface_stl_object), intent(in), target    :: surface             !< Surface (analyzed).
   integer(I4P),             intent(in)            :: facets(:)           !< Facet of each point [n].
   real(R8P),                intent(in)            :: points(:,:)         !< Points [3, n].
   real(R8P),                intent(out)           :: values(:,:)         !< Values [ncomp, n].
   integer(I4P),             intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                intent(in),  optional :: vertex_field(:,:)   !< Per-vertex field [ncomp, nv].
   real(R8P),                intent(in),  optional :: facet_field(:,:)    !< Per-facet field [ncomp, nf].
   type(rigid_pose_t),       intent(in),  optional :: pose                !< Pose of the body.
   logical,                  intent(in),  optional :: device              !< Device path (default .false.).
   real(R8P)    :: pr(3, 3), pt(3) !< Pose.
   integer(I4P) :: st, n, nc       !< Status, points, components.
   logical      :: dev             !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = size(facets, kind=I4P) ; nc = size(values, 1, kind=I4P)
   st = check_3d(surface, nc, present(vertex_field), present(facet_field), vertex_field, facet_field)
   if (size(points, 1) /= 3 .or. size(points, 2, kind=I4P) /= n .or. size(values, 2, kind=I4P) /= n) &
      st = STATUS_INVALID_INPUT
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_field_3d: refused -- one field, sizes of the points and values'
      values = 0._R8P
      return
   endif
   call pose_3d(pose, pr, pt)
   if (present(vertex_field)) then
      call run_3d(surface, MODE_VERTEX, 0_I4P, 0_I4P, n, n, nc, size(vertex_field, 2, kind=I4P), facets, points, &
                  vertex_field, present(pose), pr, pt, dev, values)
   else
      call run_3d(surface, MODE_FACET, 0_I4P, 0_I4P, n, n, nc, size(facet_field, 2, kind=I4P), facets, points, &
                  facet_field, present(pose), pr, pt, dev, values)
   endif
   endsubroutine transfer_field_3d

   subroutine transfer_field_2d(curve, segments, points, values, status, vertex_field, segment_field, pose, device)
   !< 2-D twin of `transfer_field_3d` in the curve's (u, v) plane: points on (or projected onto) the segments
   !< `segments(n)` (segment i -> next(i), by its first vertex), `points(2, n)`; exactly one of `vertex_field(ncomp, nv)`
   !< and `segment_field(ncomp, nv)` (nv the curve's vertices = segments). The missing segment at the end of an open loop
   !< gives the value of its vertex.
   type(curve_object),    intent(in)            :: curve               !< Curve.
   integer(I4P),          intent(in)            :: segments(:)         !< Segment of each point [n].
   real(R8P),             intent(in)            :: points(:,:)         !< Points [2, n].
   real(R8P),             intent(out)           :: values(:,:)         !< Values [ncomp, n].
   integer(I4P),          intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),             intent(in),  optional :: vertex_field(:,:)   !< Per-vertex field [ncomp, nv].
   real(R8P),             intent(in),  optional :: segment_field(:,:)  !< Per-segment field [ncomp, nv].
   type(rigid_pose_2d_t), intent(in),  optional :: pose                !< Pose of the curve.
   logical,               intent(in),  optional :: device              !< Device path (default .false.).
   real(R8P)    :: pr(2, 2), pt(2) !< Pose.
   integer(I4P) :: st, n, nc       !< Status, points, components.
   logical      :: dev             !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = size(segments, kind=I4P) ; nc = size(values, 1, kind=I4P)
   st = check_2d(curve, nc, present(vertex_field), present(segment_field), vertex_field, segment_field)
   if (size(points, 1) /= 2 .or. size(points, 2, kind=I4P) /= n .or. size(values, 2, kind=I4P) /= n) &
      st = STATUS_INVALID_INPUT
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_field_2d: refused -- one field, sizes of the points and values'
      values = 0._R8P
      return
   endif
   call pose_2d_arrays(pose, pr, pt)
   if (present(vertex_field)) then
      call run_2d(curve, MODE_VERTEX, 0_I4P, 0_I4P, n, n, nc, segments, points, vertex_field, present(pose), pr, pt, &
                  dev, values)
   else
      call run_2d(curve, MODE_FACET, 0_I4P, 0_I4P, n, n, nc, segments, points, segment_field, present(pose), pr, pt, &
                  dev, values)
   endif
   endsubroutine transfer_field_2d

   subroutine transfer_to_closure_3d(surface, closure, values, status, vertex_field, facet_field, ilw_values, pose, device)
   !< A wall field at the feet of a 3-D closure table, `values(ncomp, closure%n)` (on `closure%element`), and with
   !< `ilw_values(ncomp, 6, closure%n)` at its ILW crossings (the table built with `ilw=.true.`; a slot without a crossing
   !< gives zeros). `pose`: the pose the table was built with. The rest as `transfer_field_3d`; a resident table is read on
   !< the device.
   type(surface_stl_object), intent(in), target    :: surface              !< Surface (analyzed).
   type(ib_closure_t),       intent(in)            :: closure              !< Closure table (3-D).
   real(R8P),                intent(out)           :: values(:,:)          !< Values at the feet [ncomp, n].
   integer(I4P),             intent(out), optional :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                intent(in),  optional :: vertex_field(:,:)    !< Per-vertex field [ncomp, nv].
   real(R8P),                intent(in),  optional :: facet_field(:,:)     !< Per-facet field [ncomp, nf].
   real(R8P),                intent(out), optional :: ilw_values(:,:,:)    !< Values at the ILW crossings [ncomp, 6, n].
   type(rigid_pose_t),       intent(in),  optional :: pose                 !< Pose of the body (of the table).
   logical,                  intent(in),  optional :: device               !< Device path (default .false.).
   real(R8P), allocatable :: ip(:,:,:)     !< ILW crossing points [3, 6, n].
   integer(I4P), allocatable :: ifc(:,:)   !< ILW crossing facets [6, n].
   real(R8P)    :: pr(3, 3), pt(3)         !< Pose.
   integer(I4P) :: st, n, nc, mode, nfld   !< Status, ghosts, components, mode, field size.
   logical      :: dev                     !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = closure%n ; nc = size(values, 1, kind=I4P)
   st = check_3d(surface, nc, present(vertex_field), present(facet_field), vertex_field, facet_field)
   if (closure%dims /= 3_I4P .or. size(values, 2, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   if (present(ilw_values)) then
      if (.not. closure%with_ilw .or. size(ilw_values, 1, kind=I4P) /= nc .or. size(ilw_values, 2) /= 6 .or. &
          size(ilw_values, 3, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_to_closure: refused -- a 3-D table (with ILW for ilw_values), one field, sizes'
      values = 0._R8P
      if (present(ilw_values)) ilw_values = 0._R8P
      return
   endif
   call pose_3d(pose, pr, pt)
   mode = MODE_FACET ; if (present(vertex_field)) mode = MODE_VERTEX
   if (present(vertex_field)) then
      nfld = size(vertex_field, 2, kind=I4P)
      call run_3d(surface, mode, 0_I4P, 0_I4P, n, n, nc, nfld, closure%element, closure%foot, vertex_field, &
                  present(pose), pr, pt, dev, values)
   else
      nfld = size(facet_field, 2, kind=I4P)
      call run_3d(surface, mode, 0_I4P, 0_I4P, n, n, nc, nfld, closure%element, closure%foot, facet_field, &
                  present(pose), pr, pt, dev, values)
   endif
   if (.not. present(ilw_values)) return
   if (n < 1_I4P) return
   allocate(ip(3, 6, n), ifc(6, n))
   !$acc data create(ip, ifc) if(dev)
   call ilw_gather(3_I4P, n, closure%points, closure%ilw_t, closure%ilw_element, dev, ip, ifc)
   if (present(vertex_field)) then
      call run_3d(surface, mode, 0_I4P, 0_I4P, 6_I4P * n, 6_I4P * n, nc, nfld, ifc, ip, vertex_field, present(pose), &
                  pr, pt, dev, ilw_values)
   else
      call run_3d(surface, mode, 0_I4P, 0_I4P, 6_I4P * n, 6_I4P * n, nc, nfld, ifc, ip, facet_field, present(pose), &
                  pr, pt, dev, ilw_values)
   endif
   !$acc end data
   endsubroutine transfer_to_closure_3d

   subroutine transfer_to_closure_2d(curve, closure, values, status, vertex_field, segment_field, ilw_values, pose, device)
   !< 2-D twin of `transfer_to_closure_3d`: a 2-D table, `ilw_values(ncomp, 4, closure%n)`, the fields of
   !< `transfer_field_2d`.
   type(curve_object),    intent(in)            :: curve                !< Curve.
   type(ib_closure_t),    intent(in)            :: closure              !< Closure table (2-D).
   real(R8P),             intent(out)           :: values(:,:)          !< Values at the feet [ncomp, n].
   integer(I4P),          intent(out), optional :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),             intent(in),  optional :: vertex_field(:,:)    !< Per-vertex field [ncomp, nv].
   real(R8P),             intent(in),  optional :: segment_field(:,:)   !< Per-segment field [ncomp, nv].
   real(R8P),             intent(out), optional :: ilw_values(:,:,:)    !< Values at the ILW crossings [ncomp, 4, n].
   type(rigid_pose_2d_t), intent(in),  optional :: pose                 !< Pose of the curve (of the table).
   logical,               intent(in),  optional :: device               !< Device path (default .false.).
   real(R8P), allocatable :: ip(:,:,:)     !< ILW crossing points [2, 4, n].
   integer(I4P), allocatable :: ifc(:,:)   !< ILW crossing segments [4, n].
   real(R8P)    :: pr(2, 2), pt(2)         !< Pose.
   integer(I4P) :: st, n, nc, mode         !< Status, ghosts, components, mode.
   logical      :: dev                     !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = closure%n ; nc = size(values, 1, kind=I4P)
   st = check_2d(curve, nc, present(vertex_field), present(segment_field), vertex_field, segment_field)
   if (closure%dims /= 2_I4P .or. size(values, 2, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   if (present(ilw_values)) then
      if (.not. closure%with_ilw .or. size(ilw_values, 1, kind=I4P) /= nc .or. size(ilw_values, 2) /= 4 .or. &
          size(ilw_values, 3, kind=I4P) /= n) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_to_closure: refused -- a 2-D table (with ILW for ilw_values), one field, sizes'
      values = 0._R8P
      if (present(ilw_values)) ilw_values = 0._R8P
      return
   endif
   call pose_2d_arrays(pose, pr, pt)
   mode = MODE_FACET ; if (present(vertex_field)) mode = MODE_VERTEX
   if (present(vertex_field)) then
      call run_2d(curve, mode, 0_I4P, 0_I4P, n, n, nc, closure%element, closure%foot, vertex_field, present(pose), pr, &
                  pt, dev, values)
   else
      call run_2d(curve, mode, 0_I4P, 0_I4P, n, n, nc, closure%element, closure%foot, segment_field, present(pose), pr, &
                  pt, dev, values)
   endif
   if (.not. present(ilw_values)) return
   if (n < 1_I4P) return
   allocate(ip(2, 4, n), ifc(4, n))
   !$acc data create(ip, ifc) if(dev)
   call ilw_gather(2_I4P, n, closure%points, closure%ilw_t, closure%ilw_element, dev, ip, ifc)
   if (present(vertex_field)) then
      call run_2d(curve, mode, 0_I4P, 0_I4P, 4_I4P * n, 4_I4P * n, nc, ifc, ip, vertex_field, present(pose), pr, pt, &
                  dev, ilw_values)
   else
      call run_2d(curve, mode, 0_I4P, 0_I4P, 4_I4P * n, 4_I4P * n, nc, ifc, ip, segment_field, present(pose), pr, pt, &
                  dev, ilw_values)
   endif
   !$acc end data
   endsubroutine transfer_to_closure_2d

   subroutine transfer_to_samples_3d(surface, samples, values, status, vertex_field, facet_field, device)
   !< A wall field at the points of 3-D load samples, `values(ncomp, samples%n)`: a vertex field interpolated with the
   !< rule's barycentric coordinates of each sample (the coordinates that placed it: no closest point, pose-free), a facet
   !< field copied. The samples must be of this surface (as many facets). Resident samples are read on the device.
   type(surface_stl_object), intent(in), target    :: surface            !< Surface (analyzed).
   type(load_samples_t),     intent(in)            :: samples            !< Load samples (3-D) of the surface.
   real(R8P),                intent(out)           :: values(:,:)        !< Values [ncomp, n].
   integer(I4P),             intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                intent(in),  optional :: vertex_field(:,:)  !< Per-vertex field [ncomp, nv].
   real(R8P),                intent(in),  optional :: facet_field(:,:)   !< Per-facet field [ncomp, nf].
   logical,                  intent(in),  optional :: device             !< Device path (default .false.).
   real(R8P)    :: pr(3, 3), pt(3), pz(3, 1) !< Pose (unused), placeholder points.
   integer(I4P) :: st, n, nc                 !< Status, samples, components.
   logical      :: dev                       !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = samples%n ; nc = size(values, 1, kind=I4P)
   st = check_3d(surface, nc, present(vertex_field), present(facet_field), vertex_field, facet_field)
   if (samples%dims /= 3_I4P .or. samples%ne /= surface%get_facets_number() .or. size(values, 2, kind=I4P) /= n) &
      st = STATUS_INVALID_INPUT
   if (st == STATUS_OK .and. n > 0_I4P) then
      if (.not. allocated(samples%element)) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_to_samples: refused -- 3-D samples of this surface, one field, sizes'
      values = 0._R8P
      return
   endif
   call pose_3d(rp=pr, tp=pt) ; pz = 0._R8P
   if (present(vertex_field)) then
      call run_3d(surface, MODE_RULE, samples%order, samples%q, n, 1_I4P, nc, size(vertex_field, 2, kind=I4P), &
                  samples%element, pz, vertex_field, .false., pr, pt, dev, values)
   else
      call run_3d(surface, MODE_FACET, samples%order, samples%q, n, 1_I4P, nc, size(facet_field, 2, kind=I4P), &
                  samples%element, pz, facet_field, .false., pr, pt, dev, values)
   endif
   endsubroutine transfer_to_samples_3d

   subroutine transfer_to_samples_2d(curve, samples, values, status, vertex_field, segment_field, device)
   !< 2-D twin of `transfer_to_samples_3d`: (1 - s, s) with s the Gauss point of the sample on its segment.
   type(curve_object),   intent(in)            :: curve               !< Curve.
   type(load_samples_t), intent(in)            :: samples             !< Load samples (2-D) of the curve.
   real(R8P),            intent(out)           :: values(:,:)         !< Values [ncomp, n].
   integer(I4P),         intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),            intent(in),  optional :: vertex_field(:,:)   !< Per-vertex field [ncomp, nv].
   real(R8P),            intent(in),  optional :: segment_field(:,:)  !< Per-segment field [ncomp, nv].
   logical,              intent(in),  optional :: device              !< Device path (default .false.).
   real(R8P)    :: pr(2, 2), pt(2), pz(2, 1) !< Pose (unused), placeholder points.
   integer(I4P) :: st, n, nc                 !< Status, samples, components.
   logical      :: dev                       !< Device path.

   dev = .false. ; if (present(device)) dev = device
   n = samples%n ; nc = size(values, 1, kind=I4P)
   st = check_2d(curve, nc, present(vertex_field), present(segment_field), vertex_field, segment_field)
   if (samples%dims /= 2_I4P .or. samples%ne /= curve%get_vertices_number() .or. size(values, 2, kind=I4P) /= n) &
      st = STATUS_INVALID_INPUT
   if (st == STATUS_OK .and. n > 0_I4P) then
      if (.not. allocated(samples%element)) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil transfer_to_samples: refused -- 2-D samples of this curve, one field, sizes'
      values = 0._R8P
      return
   endif
   call pose_2d_arrays(rp=pr, tp=pt) ; pz = 0._R8P
   if (present(vertex_field)) then
      call run_2d(curve, MODE_RULE, samples%order, samples%q, n, 1_I4P, nc, samples%element, pz, vertex_field, .false., &
                  pr, pt, dev, values)
   else
      call run_2d(curve, MODE_FACET, samples%order, samples%q, n, 1_I4P, nc, samples%element, pz, segment_field, &
                  .false., pr, pt, dev, values)
   endif
   endsubroutine transfer_to_samples_2d

   ! points to wall

   subroutine samples_to_vertices_3d(surface, samples, values, field, status, linear, device)
   !< Values at the 3-D load samples `values(ncomp, samples%n)` to a per-vertex field `field(ncomp, nv)` (nv = the vertex
   !< pool's size; see the module header: per-facet linear fit, area-weighted vertex average; a smoothing projection).
   !< `linear` (default .true.): the per-facet linear fit (order >= 2); .false.: the facet's weighted mean. `device`
   !< (default .false.): the device path (resident samples, mapped values read there; a mapped field written there).
   !< `status = STATUS_INVALID_INPUT` (stderr message, field zero) on samples not 3-D or not of this surface, or sizes.
   type(surface_stl_object), intent(in), target    :: surface     !< Surface (analyzed).
   type(load_samples_t),     intent(in)            :: samples     !< Load samples (3-D) of the surface.
   real(R8P),                intent(in)            :: values(:,:) !< Values at the samples [ncomp, n].
   real(R8P),                intent(out)           :: field(:,:)  !< Per-vertex field [ncomp, nv].
   integer(I4P),             intent(out), optional :: status      !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,                  intent(in),  optional :: linear      !< Per-facet linear fit (default .true.).
   logical,                  intent(in),  optional :: device      !< Device path (default .false.).
   type(vertex_pool_object), pointer :: pool             !< Vertex pool.
   real(R8P),    allocatable :: cm(:,:)                  !< Fit coefficients [3, q].
   integer(I4P), allocatable :: fvid(:,:), vstart(:), vent(:) !< Facet vertices, vertex -> (facet, corner) CSR.
   integer(I4P) :: st, n, nc, nf, nv                     !< Status, samples, components, facets, vertices.
   logical      :: dev, lin                              !< Options.

   dev = .false. ; if (present(device)) dev = device
   lin = .true. ; if (present(linear)) lin = linear
   n = samples%n ; nc = size(values, 1, kind=I4P) ; nf = surface%get_facets_number()
   pool => surface%get_vertex_pool() ; nv = 0_I4P
   if (pool%get_is_initialized()) nv = pool%vertex_count()
   st = STATUS_OK
   if (samples%dims /= 3_I4P .or. samples%ne /= nf .or. size(values, 2, kind=I4P) /= n .or. &
       size(field, 1, kind=I4P) /= nc .or. size(field, 2, kind=I4P) /= nv) st = STATUS_INVALID_INPUT
   if (st == STATUS_OK .and. n > 0_I4P) then
      if (.not. allocated(samples%weight) .or. samples%q < 1_I4P) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil samples_to_vertices: refused -- 3-D samples of this surface, sizes of values and field'
      field = 0._R8P
      return
   endif
   if (n < 1_I4P .or. nv < 1_I4P .or. nc < 1_I4P) then
      field = 0._R8P
      return
   endif
   call facet_vertices(surface, fvid)
   call vertex_csr(3_I4P, nf, nv, fvid, vstart, vent)
   call fit_matrix(3_I4P, samples%order, samples%q, lin, cm)
   call s2v_run(nc, 3_I4P, samples%q, nf, n, nv, size(vent, kind=I4P), cm, samples%weight, values, vstart, vent, dev, &
                field)
   endsubroutine samples_to_vertices_3d

   subroutine samples_to_vertices_2d(curve, samples, values, field, status, linear, device)
   !< 2-D twin of `samples_to_vertices_3d`: per-segment linear fit in (1 - s, s), length-weighted vertex average over the
   !< (at most two) segments of each vertex; `field(ncomp, nv)` on the curve's vertices.
   type(curve_object),   intent(in)            :: curve       !< Curve.
   type(load_samples_t), intent(in)            :: samples     !< Load samples (2-D) of the curve.
   real(R8P),            intent(in)            :: values(:,:) !< Values at the samples [ncomp, n].
   real(R8P),            intent(out)           :: field(:,:)  !< Per-vertex field [ncomp, nv].
   integer(I4P),         intent(out), optional :: status      !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,              intent(in),  optional :: linear      !< Per-segment linear fit (default .true.).
   logical,              intent(in),  optional :: device      !< Device path (default .false.).
   real(R8P),    allocatable :: cm(:,:)                  !< Fit coefficients [2, q].
   integer(I4P), allocatable :: svid(:,:), vstart(:), vent(:) !< Segment vertices, vertex -> (segment, end) CSR.
   integer(I4P) :: st, n, nc, nv, i                      !< Status, samples, components, vertices, counter.
   logical      :: dev, lin                              !< Options.

   dev = .false. ; if (present(device)) dev = device
   lin = .true. ; if (present(linear)) lin = linear
   n = samples%n ; nc = size(values, 1, kind=I4P) ; nv = curve%get_vertices_number()
   st = STATUS_OK
   if (samples%dims /= 2_I4P .or. samples%ne /= nv .or. size(values, 2, kind=I4P) /= n .or. &
       size(field, 1, kind=I4P) /= nc .or. size(field, 2, kind=I4P) /= nv) st = STATUS_INVALID_INPUT
   if (st == STATUS_OK .and. n > 0_I4P) then
      if (.not. allocated(samples%weight) .or. samples%q < 1_I4P) st = STATUS_INVALID_INPUT
   endif
   if (present(status)) status = st
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil samples_to_vertices: refused -- 2-D samples of this curve, sizes of values and field'
      field = 0._R8P
      return
   endif
   if (n < 1_I4P .or. nv < 1_I4P .or. nc < 1_I4P) then
      field = 0._R8P
      return
   endif
   allocate(svid(2, nv))
   do i = 1_I4P, nv
      svid(1, i) = i ; svid(2, i) = curve%get_next(i) ! 0 at the end of an open loop: no entry
   enddo
   call vertex_csr(2_I4P, nv, nv, svid, vstart, vent)
   call fit_matrix(2_I4P, samples%order, samples%q, lin, cm)
   call s2v_run(nc, 2_I4P, samples%q, nv, n, nv, size(vent, kind=I4P), cm, samples%weight, values, vstart, vent, dev, &
                field)
   endsubroutine samples_to_vertices_2d

   ! private: checks and inputs

   integer(I4P) function check_3d(surface, nc, hv, hf, vertex_field, facet_field) result(st)
   !< The checks common to the 3-D wall-to-points routines: exactly one field, its sizes (ncomp, nv or nf), an analyzed
   !< surface.
   type(surface_stl_object), intent(in), target    :: surface            !< Surface.
   integer(I4P),             intent(in)            :: nc                 !< Components.
   logical,                  intent(in)            :: hv, hf             !< Fields given.
   real(R8P),                intent(in),  optional :: vertex_field(:,:)  !< Per-vertex field.
   real(R8P),                intent(in),  optional :: facet_field(:,:)   !< Per-facet field.
   type(vertex_pool_object), pointer :: pool !< Vertex pool.
   integer(I4P) :: nf, nv                    !< Facets, vertices.

   st = STATUS_OK
   nf = surface%get_facets_number()
   pool => surface%get_vertex_pool() ; nv = 0_I4P
   if (pool%get_is_initialized()) nv = pool%vertex_count()
   if (hv .eqv. hf) then
      st = STATUS_INVALID_INPUT
   elseif (hv) then
      if (size(vertex_field, 1, kind=I4P) /= nc .or. size(vertex_field, 2, kind=I4P) /= nv) st = STATUS_INVALID_INPUT
      if (nf > 0_I4P .and. .not. allocated(surface%tri_vertices)) st = STATUS_INVALID_INPUT
   else
      if (size(facet_field, 1, kind=I4P) /= nc .or. size(facet_field, 2, kind=I4P) /= nf) st = STATUS_INVALID_INPUT
   endif
   endfunction check_3d

   integer(I4P) function check_2d(curve, nc, hv, hs, vertex_field, segment_field) result(st)
   !< The checks common to the 2-D wall-to-points routines: exactly one field, its sizes (ncomp, nv).
   type(curve_object), intent(in)           :: curve               !< Curve.
   integer(I4P),       intent(in)           :: nc                  !< Components.
   logical,            intent(in)           :: hv, hs              !< Fields given.
   real(R8P),          intent(in), optional :: vertex_field(:,:)   !< Per-vertex field.
   real(R8P),          intent(in), optional :: segment_field(:,:)  !< Per-segment field.
   integer(I4P) :: nv                                              !< Vertices.

   st = STATUS_OK
   nv = curve%get_vertices_number()
   if (hv .eqv. hs) then
      st = STATUS_INVALID_INPUT
   elseif (hv) then
      if (size(vertex_field, 1, kind=I4P) /= nc .or. size(vertex_field, 2, kind=I4P) /= nv) st = STATUS_INVALID_INPUT
   else
      if (size(segment_field, 1, kind=I4P) /= nc .or. size(segment_field, 2, kind=I4P) /= nv) st = STATUS_INVALID_INPUT
   endif
   endfunction check_2d

   subroutine pose_3d(pose, rp, tp)
   !< The pose's rotation and translation (the identity without a pose).
   type(rigid_pose_t), intent(in), optional :: pose     !< Pose.
   real(R8P),          intent(out)          :: rp(3, 3) !< Rotation.
   real(R8P),          intent(out)          :: tp(3)    !< Translation.

   rp = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ; tp = 0._R8P
   if (present(pose)) then
      rp = pose%r ; tp = pose%t
   endif
   endsubroutine pose_3d

   subroutine pose_2d_arrays(pose, rp, tp)
   !< The 2-D pose's rotation and translation (the identity without a pose).
   type(rigid_pose_2d_t), intent(in), optional :: pose     !< Pose.
   real(R8P),             intent(out)          :: rp(2, 2) !< Rotation.
   real(R8P),             intent(out)          :: tp(2)    !< Translation.

   rp = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) ; tp = 0._R8P
   if (present(pose)) then
      rp = pose%r ; tp = pose%t
   endif
   endsubroutine pose_2d_arrays

   subroutine facet_vertices(surface, fvid)
   !< The vertex-pool ids of the corners of every facet [3, nf].
   type(surface_stl_object),  intent(in), target :: surface   !< Surface (analyzed).
   integer(I4P), allocatable, intent(out)        :: fvid(:,:) !< Facet vertices.
   type(vertex_pool_object), pointer :: pool !< Vertex pool.
   integer(I4P) :: nf, f, k                  !< Facets, counters.

   nf = surface%get_facets_number()
   allocate(fvid(3, nf))
   pool => surface%get_vertex_pool()
   do f = 1_I4P, nf
      do k = 1_I4P, 3_I4P
         fvid(k, f) = pool%facet_vid(f, k)
      enddo
   enddo
   endsubroutine facet_vertices

   subroutine vertex_csr(kc, ne, nv, evid, vstart, vent)
   !< The vertex -> (element, corner) incidence as CSR: row v lists the entries (e - 1) kc + k with evid(k, e) = v, in
   !< increasing e then k (corners with an id outside 1..nv are skipped). vent has at least one slot.
   integer(I4P),              intent(in)  :: kc, ne, nv     !< Corners per element, elements, vertices.
   integer(I4P),              intent(in)  :: evid(kc, ne)   !< Element corners.
   integer(I4P), allocatable, intent(out) :: vstart(:)      !< Row starts [nv + 1].
   integer(I4P), allocatable, intent(out) :: vent(:)        !< Entries.
   integer(I4P), allocatable :: pos(:)                      !< Fill positions.
   integer(I4P) :: e, k, v                                  !< Counters.

   allocate(vstart(nv + 1), pos(nv)) ; pos = 0_I4P
   do e = 1_I4P, ne
      do k = 1_I4P, kc
         v = evid(k, e)
         if (v >= 1_I4P .and. v <= nv) pos(v) = pos(v) + 1_I4P
      enddo
   enddo
   vstart(1) = 1_I4P
   do v = 1_I4P, nv
      vstart(v + 1) = vstart(v) + pos(v)
   enddo
   allocate(vent(max(1_I4P, vstart(nv + 1) - 1_I4P))) ; vent = 0_I4P
   pos = vstart(1:nv)
   do e = 1_I4P, ne
      do k = 1_I4P, kc
         v = evid(k, e)
         if (v >= 1_I4P .and. v <= nv) then
            vent(pos(v)) = (e - 1_I4P) * kc + k ; pos(v) = pos(v) + 1_I4P
         endif
      enddo
   enddo
   endsubroutine vertex_csr

   subroutine fit_matrix(d, ord, q, lin, cm)
   !< The coefficients cm(k, iq) of the per-element fit: the value at corner k = sum_iq cm(k, iq) (value_iq - value_1) +
   !< value_1. Linear: weighted least squares of the rule's barycentric coordinates lambda_iq with its weights w_iq,
   !< cm = M^-1 [w_iq lambda_iq], M = sum_iq w_iq lambda_iq lambda_iq^T (q >= kc; one point: the constant); mean: cm(k, iq)
   !< = w_iq. Host code, once per call; the same coefficients on both paths.
   integer(I4P),           intent(in)  :: d, ord, q !< Dimensions (corners), rule order, points.
   logical,                intent(in)  :: lin       !< Linear fit.
   real(R8P), allocatable, intent(out) :: cm(:,:)   !< Coefficients [d, q].
   real(R8P)    :: la(d, q), wr(q), m(d, d), mi(d, d), pu, pv, ww, nu, nw, det !< Work.
   integer(I4P) :: iq, i, j                         !< Counters.

   allocate(cm(d, q))
   do iq = 1_I4P, q
      if (d == 3_I4P) then
         call rule_3d_oac(ord, iq, la(1, iq), la(2, iq), la(3, iq), wr(iq))
      else
         call load_point_2d_oac(ord, iq, 0._R8P, 0._R8P, 1._R8P, 0._R8P, pu, pv, ww, nu, nw)
         la(1, iq) = 1._R8P - pu ; la(2, iq) = pu ; wr(iq) = ww
      endif
   enddo
   if (q == 1_I4P) then
      cm = 1._R8P
   elseif (.not. lin) then
      do iq = 1_I4P, q
         cm(:, iq) = wr(iq)
      enddo
   else
      m = 0._R8P
      do iq = 1_I4P, q
         do j = 1_I4P, d
            do i = 1_I4P, d
               m(i, j) = m(i, j) + wr(iq) * la(i, iq) * la(j, iq)
            enddo
         enddo
      enddo
      if (d == 3_I4P) then
         mi(1, 1) = m(2, 2) * m(3, 3) - m(2, 3) * m(3, 2)
         mi(1, 2) = m(1, 3) * m(3, 2) - m(1, 2) * m(3, 3)
         mi(1, 3) = m(1, 2) * m(2, 3) - m(1, 3) * m(2, 2)
         mi(2, 1) = m(2, 3) * m(3, 1) - m(2, 1) * m(3, 3)
         mi(2, 2) = m(1, 1) * m(3, 3) - m(1, 3) * m(3, 1)
         mi(2, 3) = m(1, 3) * m(2, 1) - m(1, 1) * m(2, 3)
         mi(3, 1) = m(2, 1) * m(3, 2) - m(2, 2) * m(3, 1)
         mi(3, 2) = m(1, 2) * m(3, 1) - m(1, 1) * m(3, 2)
         mi(3, 3) = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
         det = m(1, 1) * mi(1, 1) + m(1, 2) * mi(2, 1) + m(1, 3) * mi(3, 1)
      else
         mi(1, 1) = m(2, 2) ; mi(1, 2) = -m(1, 2) ; mi(2, 1) = -m(2, 1) ; mi(2, 2) = m(1, 1)
         det = m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)
      endif
      mi = mi / det
      do iq = 1_I4P, q
         do i = 1_I4P, d
            cm(i, iq) = wr(iq) * sum(mi(i, :) * la(:, iq))
         enddo
      enddo
   endif
   endsubroutine fit_matrix

   ! private: the runs (host and device in one set of loops)

   subroutine run_3d(surface, mode, ord, q, n, npt, nc, nfld, fac, pts, fld, posed, pr, pt, dev, val)
   !< The 3-D wall-to-points run after the checks: the facet vertex ids (vertex modes) gathered on the host, the facet
   !< vertices passed only to the closest-point mode, then one stage.
   type(surface_stl_object), intent(in), target :: surface          !< Surface.
   integer(I4P),             intent(in)         :: mode, ord, q     !< Mode, rule order and points (rule mode).
   integer(I4P),             intent(in)         :: n, npt, nc, nfld !< Points, size of pts, components, field size.
   integer(I4P),             intent(in)         :: fac(n)           !< Facets.
   real(R8P),                intent(in)         :: pts(3, npt)      !< Points (closest-point mode).
   real(R8P),                intent(in)         :: fld(nc, nfld)    !< Field.
   logical,                  intent(in)         :: posed, dev       !< A pose, device.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   real(R8P),                intent(inout)      :: val(nc, n)       !< Values.
   integer(I4P), allocatable :: fvid(:,:)                           !< Facet vertices.
   real(R8P)                 :: tvz(9, 1)                           !< Placeholder vertices.
   integer(I4P)              :: nf                                  !< Facets.

   if (n < 1_I4P .or. nc < 1_I4P) return
   nf = surface%get_facets_number()
   if (nf < 1_I4P) then
      val = 0._R8P
      return
   endif
   if (mode == MODE_FACET) then
      allocate(fvid(3, 1)) ; fvid = 1_I4P
   else
      call facet_vertices(surface, fvid)
   endif
   tvz = 0._R8P
   if (mode == MODE_VERTEX) then
      call stage_3d(mode, ord, q, n, npt, nc, nfld, nf, nf, size(fvid, 2, kind=I4P), fac, pts, surface%tri_vertices, &
                    fvid, fld, posed, pr, pt, dev, val)
   else
      call stage_3d(mode, ord, q, n, npt, nc, nfld, nf, 1_I4P, size(fvid, 2, kind=I4P), fac, pts, tvz, fvid, fld, &
                    posed, pr, pt, dev, val)
   endif
   endsubroutine run_3d

   subroutine stage_3d(mode, ord, q, n, npt, nc, nfld, nf, ntv, nfv, fac, pts, tv, fvid, fld, posed, pr, pt, dev, val)
   !< One thread per point: zero for a facet outside 1..nf; the facet's value (facet mode); the vertex field with the
   !< rule's barycentric coordinates (rule mode: point i is sample iq = i - (f - 1) q of facet f) or the closest point's
   !< (`transfer_weights_3d_oac`, the point taken to the body frame with a pose). The inputs copied in on the device path
   !< unless mapped (read there), the values copied out unless mapped (written there).
   integer(I4P), intent(in)    :: mode, ord, q           !< Mode, rule order and points.
   integer(I4P), intent(in)    :: n, npt, nc, nfld       !< Points, size of pts, components, field size.
   integer(I4P), intent(in)    :: nf, ntv, nfv           !< Facets, sizes of tv and fvid.
   integer(I4P), intent(in)    :: fac(n)                 !< Facets.
   real(R8P),    intent(in)    :: pts(3, npt)            !< Points.
   real(R8P),    intent(in)    :: tv(9, ntv)             !< Facet vertices (closest-point mode).
   integer(I4P), intent(in)    :: fvid(3, nfv)           !< Facet vertex ids (vertex modes).
   real(R8P),    intent(in)    :: fld(nc, nfld)          !< Field.
   logical,      intent(in)    :: posed, dev             !< A pose, device.
   real(R8P),    intent(in)    :: pr(3, 3), pt(3)        !< Pose.
   real(R8P),    intent(inout) :: val(nc, n)             !< Values.
   real(R8P)    :: r(3, 3), t(3), bx, by, bz, l1, l2, l3, wr !< Work.
   integer(I4P) :: i, f, c, iq, i1, i2, i3               !< Work.

   r = pr ; t = pt
   !$acc data copyin(fac, pts, tv, fvid, fld) copyout(val) if(dev)
   !$acc parallel loop gang vector default(present) present(fac, pts, tv, fvid, fld, val) copyin(r, t) &
   !$acc&  private(f, c, iq, i1, i2, i3, bx, by, bz, l1, l2, l3, wr) if(dev)
   do i = 1_I4P, n
      f = fac(i)
      if (f < 1_I4P .or. f > nf) then
         do c = 1_I4P, nc
            val(c, i) = 0._R8P
         enddo
      elseif (mode == MODE_FACET) then
         do c = 1_I4P, nc
            val(c, i) = fld(c, f)
         enddo
      else
         if (mode == MODE_RULE) then
            iq = i - (f - 1_I4P) * q
            call rule_3d_oac(ord, iq, l1, l2, l3, wr)
         else
            if (posed) then
               call to_body_3d_oac(r, t, pts(1, i), pts(2, i), pts(3, i), bx, by, bz)
            else
               bx = pts(1, i) ; by = pts(2, i) ; bz = pts(3, i)
            endif
            call transfer_weights_3d_oac(bx, by, bz, tv(1, f), tv(2, f), tv(3, f), tv(4, f), tv(5, f), tv(6, f), &
                                         tv(7, f), tv(8, f), tv(9, f), l1, l2, l3)
         endif
         i1 = fvid(1, f) ; i2 = fvid(2, f) ; i3 = fvid(3, f)
         do c = 1_I4P, nc
            val(c, i) = (l1 * fld(c, i1) + l2 * fld(c, i2)) + l3 * fld(c, i3)
         enddo
      endif
   enddo
   !$acc end data
   endsubroutine stage_3d

   subroutine run_2d(curve, mode, ord, q, n, npt, nc, fac, pts, fld, posed, pr, pt, dev, val)
   !< The 2-D wall-to-points run after the checks: the curve's vertices and next gathered on the host (vertex modes),
   !< then one stage.
   type(curve_object), intent(in)    :: curve            !< Curve.
   integer(I4P),       intent(in)    :: mode, ord, q     !< Mode, rule order and points (rule mode).
   integer(I4P),       intent(in)    :: n, npt, nc       !< Points, size of pts, components.
   integer(I4P),       intent(in)    :: fac(n)           !< Segments.
   real(R8P),          intent(in)    :: pts(2, npt)      !< Points (closest-point mode).
   real(R8P),          intent(in)    :: fld(:,:)         !< Field [nc, nv].
   logical,            intent(in)    :: posed, dev       !< A pose, device.
   real(R8P),          intent(in)    :: pr(2, 2), pt(2)  !< Pose.
   real(R8P),          intent(inout) :: val(nc, n)       !< Values.
   real(R8P),    allocatable :: u(:), v(:)               !< Vertices.
   integer(I4P), allocatable :: nxt(:)                   !< Next vertex.
   real(R8P)                 :: uv(2)                    !< A vertex.
   integer(I4P)              :: nv, i                    !< Vertices, counter.

   if (n < 1_I4P .or. nc < 1_I4P) return
   nv = curve%get_vertices_number()
   if (nv < 1_I4P) then
      val = 0._R8P
      return
   endif
   allocate(u(nv), v(nv), nxt(nv))
   do i = 1_I4P, nv
      uv = curve%get_vertex(i) ; u(i) = uv(1) ; v(i) = uv(2) ; nxt(i) = curve%get_next(i)
   enddo
   call stage_2d(mode, ord, q, n, npt, nc, nv, fac, pts, u, v, nxt, fld, posed, pr, pt, dev, val)
   endsubroutine run_2d

   subroutine stage_2d(mode, ord, q, n, npt, nc, nv, fac, pts, u, v, nxt, fld, posed, pr, pt, dev, val)
   !< 2-D twin of `stage_3d`: segment e -> nxt(e) (its first vertex's value where nxt(e) < 1), weights (1 - s, s) from the
   !< rule's Gauss point (rule mode) or the clamped projection (`transfer_weights_2d_oac`).
   integer(I4P), intent(in)    :: mode, ord, q           !< Mode, rule order and points.
   integer(I4P), intent(in)    :: n, npt, nc, nv         !< Points, size of pts, components, vertices.
   integer(I4P), intent(in)    :: fac(n)                 !< Segments.
   real(R8P),    intent(in)    :: pts(2, npt)            !< Points.
   real(R8P),    intent(in)    :: u(nv), v(nv)           !< Vertices.
   integer(I4P), intent(in)    :: nxt(nv)                !< Next vertex.
   real(R8P),    intent(in)    :: fld(nc, nv)            !< Field.
   logical,      intent(in)    :: posed, dev             !< A pose, device.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2)        !< Pose.
   real(R8P),    intent(inout) :: val(nc, n)             !< Values.
   real(R8P)    :: r(2, 2), t(2), bu, bv, s, su, sv, wr, nu, nw !< Work.
   integer(I4P) :: i, e, j, c, iq                        !< Work.

   r = pr ; t = pt
   !$acc data copyin(fac, pts, u, v, nxt, fld) copyout(val) if(dev)
   !$acc parallel loop gang vector default(present) present(fac, pts, u, v, nxt, fld, val) copyin(r, t) &
   !$acc&  private(e, j, c, iq, bu, bv, s, su, sv, wr, nu, nw) if(dev)
   do i = 1_I4P, n
      e = fac(i)
      if (e < 1_I4P .or. e > nv) then
         do c = 1_I4P, nc
            val(c, i) = 0._R8P
         enddo
      elseif (mode == MODE_FACET) then
         do c = 1_I4P, nc
            val(c, i) = fld(c, e)
         enddo
      else
         j = nxt(e)
         if (j < 1_I4P) then
            s = 0._R8P ; j = e
         elseif (mode == MODE_RULE) then
            iq = i - (e - 1_I4P) * q
            call load_point_2d_oac(ord, iq, 0._R8P, 0._R8P, 1._R8P, 0._R8P, su, sv, wr, nu, nw)
            s = su
         else
            if (posed) then
               call to_body_2d_oac(r, t, pts(1, i), pts(2, i), bu, bv)
            else
               bu = pts(1, i) ; bv = pts(2, i)
            endif
            call transfer_weights_2d_oac(bu, bv, u(e), v(e), u(j), v(j), s)
         endif
         do c = 1_I4P, nc
            val(c, i) = (1._R8P - s) * fld(c, e) + s * fld(c, j)
         enddo
      endif
   enddo
   !$acc end data
   endsubroutine stage_2d

   subroutine ilw_gather(d, n, gp, ilt, ile, dev, ip, ifc)
   !< One thread per (ghost, slot): the ILW crossing point (the ghost's coordinates, the slot's axis (sl + 1) / 2 set to
   !< ilw_t) and its element. The table's arrays read on the device where mapped (a resident table), else copied in;
   !< ip, ifc present.
   integer(I4P), intent(in)    :: d, n                !< Dimensions, ghosts.
   real(R8P),    intent(in)    :: gp(d, n)            !< Ghost points.
   real(R8P),    intent(in)    :: ilt(2 * d, n)       !< ILW crossing coordinates.
   integer(I4P), intent(in)    :: ile(2 * d, n)       !< ILW crossing elements.
   logical,      intent(in)    :: dev                 !< Device.
   real(R8P),    intent(inout) :: ip(d, 2 * d, n)     !< Crossing points.
   integer(I4P), intent(inout) :: ifc(2 * d, n)       !< Crossing elements.
   integer(I4P) :: g, sl, c, a                        !< Counters.

   !$acc data copyin(gp, ilt, ile) if(dev)
   !$acc parallel loop gang vector collapse(2) default(present) present(gp, ilt, ile, ip, ifc) private(c, a) if(dev)
   do g = 1_I4P, n
      do sl = 1_I4P, 2_I4P * d
         a = (sl + 1_I4P) / 2_I4P
         do c = 1_I4P, d
            ip(c, sl, g) = gp(c, g)
         enddo
         ip(a, sl, g) = ilt(sl, g)
         ifc(sl, g) = ile(sl, g)
      enddo
   enddo
   !$acc end data
   endsubroutine ilw_gather

   subroutine s2v_run(nc, kc, q, ne, n, nv, nent, cm, wt, val, vstart, vent, dev, field)
   !< `samples_to_vertices` after the checks: the per-element fits and areas (one thread per element), then the vertex
   !< averages (one thread per vertex, its row in order). The samples' weights and the values read on the device where
   !< mapped, else copied in; the field written on the device where mapped, else copied out.
   integer(I4P), intent(in)    :: nc, kc, q, ne, n, nv, nent !< Components, corners, points, elements, samples, vertices.
   real(R8P),    intent(in)    :: cm(kc, q)          !< Fit coefficients.
   real(R8P),    intent(in)    :: wt(n)              !< Sample weights.
   real(R8P),    intent(in)    :: val(nc, n)         !< Values at the samples.
   integer(I4P), intent(in)    :: vstart(nv + 1)     !< CSR row starts.
   integer(I4P), intent(in)    :: vent(nent)         !< CSR entries.
   logical,      intent(in)    :: dev                !< Device.
   real(R8P),    intent(inout) :: field(nc, nv)      !< Per-vertex field.
   real(R8P), allocatable :: ar(:), fe(:,:,:)        !< Element areas, corner values [nc, kc, ne].

   allocate(ar(ne), fe(nc, kc, ne))
   !$acc data copyin(cm, wt, val, vstart, vent) create(ar, fe) copyout(field) if(dev)
   call s2v_elements(nc, kc, q, ne, n, cm, wt, val, dev, ar, fe)
   call s2v_vertices(nc, kc, ne, nv, nent, vstart, vent, ar, fe, dev, field)
   !$acc end data
   endsubroutine s2v_run

   ! nvfortran 26.9 vectorizes the host copy of the inner sums at -O2 ("vector simd code for the loop containing
   ! reductions"): partial sums in another order, 1-ulp host/device differences (1142 of 24330 dragon-fine values).
   ! novector for these two routines (the device code is unchanged); a comment to every other compiler.
!pgi$r novector
   subroutine s2v_elements(nc, kc, q, ne, n, cm, wt, val, dev, ar, fe)
   !< One thread per element: its area (the sum of its samples' weights in order) and the fit's value at each corner,
   !< value_1 + sum_iq cm(k, iq) (value_iq - value_1).
   integer(I4P), intent(in)    :: nc, kc, q, ne, n   !< Components, corners, points, elements, samples.
   real(R8P),    intent(in)    :: cm(kc, q)          !< Fit coefficients (present).
   real(R8P),    intent(in)    :: wt(n)              !< Sample weights (present).
   real(R8P),    intent(in)    :: val(nc, n)         !< Values (present).
   logical,      intent(in)    :: dev                !< Device.
   real(R8P),    intent(inout) :: ar(ne)             !< Areas (present).
   real(R8P),    intent(inout) :: fe(nc, kc, ne)     !< Corner values (present).
   real(R8P)    :: a, ref, acc                       !< Work.
   integer(I4P) :: e, s0, iq, c, k                   !< Work.

   !$acc parallel loop gang vector default(present) present(cm, wt, val, ar, fe) private(s0, a, ref, acc, iq, c, k) if(dev)
   do e = 1_I4P, ne
      s0 = (e - 1_I4P) * q
      a = 0._R8P
      do iq = 1_I4P, q
         a = a + wt(s0 + iq)
      enddo
      ar(e) = a
      do c = 1_I4P, nc
         ref = val(c, s0 + 1_I4P)
         do k = 1_I4P, kc
            acc = 0._R8P
            do iq = 1_I4P, q
               acc = acc + cm(k, iq) * (val(c, s0 + iq) - ref)
            enddo
            fe(c, k, e) = ref + acc
         enddo
      enddo
   enddo
   endsubroutine s2v_elements

!pgi$r novector
   subroutine s2v_vertices(nc, kc, ne, nv, nent, vstart, vent, ar, fe, dev, field)
   !< One thread per vertex: over its incident elements of positive area in row order (increasing element), the
   !< area-weighted mean of their corner values, as the first one plus the weighted mean of the deviations from it
   !< (Neumaier sums); zeros without such an element.
   integer(I4P), intent(in)    :: nc, kc, ne, nv, nent !< Components, corners, elements, vertices, entries.
   integer(I4P), intent(in)    :: vstart(nv + 1)     !< CSR row starts (present).
   integer(I4P), intent(in)    :: vent(nent)         !< CSR entries (present).
   real(R8P),    intent(in)    :: ar(ne)             !< Areas (present).
   real(R8P),    intent(in)    :: fe(nc, kc, ne)     !< Corner values (present).
   logical,      intent(in)    :: dev                !< Device.
   real(R8P),    intent(inout) :: field(nc, nv)      !< Field (present).
   real(R8P)    :: den, cden, sm, cs, ref, a         !< Work.
   integer(I4P) :: v, j, jf, en, e, k, c             !< Work.

   !$acc parallel loop gang vector default(present) present(vstart, vent, ar, fe, field) &
   !$acc&  private(den, cden, sm, cs, ref, a, j, jf, en, e, k, c) if(dev)
   do v = 1_I4P, nv
      den = 0._R8P ; cden = 0._R8P ; jf = 0_I4P
      do j = vstart(v), vstart(v + 1) - 1_I4P
         e = (vent(j) - 1_I4P) / kc + 1_I4P
         a = ar(e)
         if (a > 0._R8P) then
            if (jf == 0_I4P) jf = j
            call neumaier_oac(den, cden, a)
         endif
      enddo
      if (jf == 0_I4P) then
         do c = 1_I4P, nc
            field(c, v) = 0._R8P
         enddo
      else
         en = vent(jf) ; e = (en - 1_I4P) / kc + 1_I4P ; k = en - (e - 1_I4P) * kc
         den = den + cden
         do c = 1_I4P, nc
            ref = fe(c, k, e) ; sm = 0._R8P ; cs = 0._R8P
            do j = jf, vstart(v + 1) - 1_I4P
               en = vent(j) ; e = (en - 1_I4P) / kc + 1_I4P
               a = ar(e)
               if (a > 0._R8P) call neumaier_oac(sm, cs, a * (fe(c, en - (e - 1_I4P) * kc, e) - ref))
            enddo
            field(c, v) = ref + (sm + cs) / den
         enddo
      endif
   enddo
   endsubroutine s2v_vertices

   ! the per-point kernels (`!$acc routine seq`, shared by the two paths)

   pure subroutine transfer_weights_3d_oac(px, py, pz, x1, y1, z1, x2, y2, z2, x3, y3, z3, l1, l2, l3)
   !< Kernel: the barycentric weights (l1, l2, l3) = (1 - s - t, s, t) of the closest point to p on the triangle
   !< (v1, v2, v3), (s, t) from `triangle_closest_st` (the solve of the distance queries); a point equal to a vertex
   !< takes (1, 0, 0) on it exactly; a triangle with all three vertices equal (or NaN edges) gives (1, 0, 0). The ONE
   !< routine of this module with `vector_R8P` locals (nvfortran duplicate-type bug, CLAUDE.md "Device-code rules").
   !$acc routine seq
   real(R8P), intent(in)  :: px, py, pz                         !< Point (body frame).
   real(R8P), intent(in)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3 !< Vertices.
   real(R8P), intent(out) :: l1, l2, l3                         !< Weights.
   type(vector_R8P) :: v1, e12, e13, p                          !< Triangle, point.
   real(R8P)        :: a, b, c, det, sq, tq, d2                 !< Gram entries, coordinates, distance.

   l1 = 0._R8P ; l2 = 0._R8P ; l3 = 0._R8P
   if (px == x1 .and. py == y1 .and. pz == z1) then
      l1 = 1._R8P
   elseif (px == x2 .and. py == y2 .and. pz == z2) then
      l2 = 1._R8P
   elseif (px == x3 .and. py == y3 .and. pz == z3) then
      l3 = 1._R8P
   else
      v1%x = x1 ; v1%y = y1 ; v1%z = z1
      e12%x = x2 - x1 ; e12%y = y2 - y1 ; e12%z = z2 - z1
      e13%x = x3 - x1 ; e13%y = y3 - y1 ; e13%z = z3 - z1
      p%x = px ; p%y = py ; p%z = pz
      a = (e12%x * e12%x + e12%y * e12%y) + e12%z * e12%z
      b = (e12%x * e13%x + e12%y * e13%y) + e12%z * e13%z
      c = (e13%x * e13%x + e13%y * e13%y) + e13%z * e13%z
      det = a * c - b * b
      if (a > 0._R8P .or. c > 0._R8P) then
         call triangle_closest_st(v1, e12, e13, a, b, c, det, p, sq, tq, d2)
         l1 = (1._R8P - sq) - tq ; l2 = sq ; l3 = tq
      else
         l1 = 1._R8P
      endif
   endif
   endsubroutine transfer_weights_3d_oac

   pure subroutine transfer_weights_2d_oac(pu, pv, u1, v1, u2, v2, s)
   !< Kernel: the parameter s in [0, 1] of the closest point to p on the segment (u1, v1) -> (u2, v2), the clamped
   !< projection ((p - a) . d) / (d . d) (the weights (1 - s, s)); a point equal to an end gives 0 or 1 exactly; a zero
   !< segment gives 0.
   !$acc routine seq
   real(R8P), intent(in)  :: pu, pv         !< Point (body frame).
   real(R8P), intent(in)  :: u1, v1, u2, v2 !< Segment ends.
   real(R8P), intent(out) :: s              !< Parameter.
   real(R8P) :: du, dv, ll                  !< Direction, its squared length.

   du = u2 - u1 ; dv = v2 - v1
   ll = du * du + dv * dv
   s = 0._R8P
   if (ll > 0._R8P) then
      s = ((pu - u1) * du + (pv - v1) * dv) / ll
      s = min(max(s, 0._R8P), 1._R8P)
   endif
   endsubroutine transfer_weights_2d_oac

   pure subroutine neumaier_oac(sm, cm, x)
   !< Kernel: add x to the running sum sm with the compensation cm (Neumaier): the total is sm + cm.
   !$acc routine seq
   real(R8P), intent(inout) :: sm, cm !< Sum, compensation.
   real(R8P), intent(in)    :: x      !< Term.
   real(R8P) :: t                     !< New sum.

   t = sm + x
   if (abs(sm) >= abs(x)) then
      cm = cm + ((sm - t) + x)
   else
      cm = cm + ((x - t) + sm)
   endif
   sm = t
   endsubroutine neumaier_oac
endmodule fossil_field_transfer
