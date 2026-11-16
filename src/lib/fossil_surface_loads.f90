!< FOSSIL, surface loads: the sample geometry of the wall integrals of a solver (#109 WP2 step 01).
module fossil_surface_loads
!< FOSSIL, surface loads: the sample geometry of the wall integrals of a solver (#109 WP2 step 01).
!<
!< The force on a body is the integral over its wall of (-p n + tau) dA, n the outward normal (into the fluid); FOSSIL
!< gives the geometry of that integral, the solver evaluates p and tau. `load_samples_3d` / `load_samples_2d` give, for
!< every facet (2-D: segment), the q points of a quadrature rule of order `order` on the element itself and per point:
!<
!< - the point, its weight (the rule's weight times the facet area or the segment length: the weights of an element sum to
!<   its area to round-off), the element's unit normal from its vertices (3-D: (v2 - v1) x (v3 - v1) / |...|, the vertex
!<   winding, outward on a sanitized surface; 2-D: (dv, -du) / |...| of the segment i -> next(i), the solid on its left),
!<   the element (3-D the facet, 2-D the segment = its first vertex), its patch and its body;
!< - m probes along the normal at the caller's signed lengths (+ into the fluid, absolute units, the same for every
!<   point; a zero length is the point itself): probe(:, k, s) = point + lengths(k) normal; visible(k, s) when the open
!<   segment from the point to the probe crosses no element (the kernels of `probe_many`, `fossil_segment_visibility`:
!<   3-D a crossing at a parameter above 1e-9 of the surface's box diagonal, the element of the point crossed at 0 and
!<   ignored; 2-D exact, the segment of the point excluded); probe_cell(a, k, s) the enclosing cell of the probe on the
!<   solver's lines x, y[, z], low-side rule (x(i) < p <= x(i+1) gives i, 0 outside), as `fd_closure` gives it. The lines
!<   are optional: without them there are no cells.
!<
!< Rules (`load_rule_size`): 3-D order 1 the centroid (exact for degree 1), order 2 the 3-point rule (2/3, 1/6, 1/6)
!< (degree 2), order 3 the 6-point rule of Dunavant (1985, degree 4; positive weights, interior points; its constants
!< solved to 40 digits); 2-D Gauss-Legendre with 1, 2, 3 points (degrees 1, 3, 5). A degenerate element (zero area or
!< length, or the missing segment at the end of an open loop) gets zero weights and a zero normal; its probes sit at the
!< point and are not visible. The samples of element e are s = (e - 1) q + 1 ... e q.
!<
!< Poses (`pose=`): the body at the pose, every output in the world frame. As the probes of `fd_closure`: the samples,
!< probes and visibility are computed in the body frame on the given geometry (its tree), the points and probes mapped
!< by `to_world`, the normals by `rotate`; no moved copy and no refitted tree (no decision is taken on the moved
!< geometry); the cells are those of the world probes.
!<
!< Two paths, the same samples bit for bit: the host and the device (`device=.true.`; the surface or curve mapped by
!< `enter_device`), one set of loops with `if(dev)` and the shared kernels (`load_point_3d_oac`, `load_point_2d_oac`, the
!< visibility kernels, `fd_cell_index_oac`; 2-D the device visibility through `curve%probes_present`, the host through
!< `curve%path_visible`, the same kernel). `resident=.true.` (implies the device path) leaves the samples mapped: the
!< resident pattern of the closure table (`on_device`, `device_key`, `enter_device` / `update_host` / `exit_device`;
!< `load_samples_reset`; every writer takes the samples intent(inout) and releases them first). Built without OpenACC
!< the device loops run on the host in order.
!<
!< Loads (#109 WP2 step 02): `integrate_loads(samples, force, moment, ...)` sums, with the solver's pressure p(n) and/or
!< viscous traction tau(d, n) at the samples, F = sum_s w_s (-p_s n_s + tau_s) and M = sum_s (x_s - c) x f_s (c = `about`,
!< default the origin; 2-D the z component, per unit depth), the totals per patch and per body and the wetted area
!< sum_s w_s. Each sample's contribution comes from one kernel (`load_force_3d_oac` / `load_force_2d_oac`), on the device
!< with `device=.true.` (the samples and p, tau read where they are mapped, else copied in); the contributions are then
!< summed on the host in sample order with Neumaier compensation, the same code for both paths: the totals agree bit for
!< bit. A linear p = p0 + g . x is integrated exactly by every rule for the force (F = -g V); the moment's integrand is
!< then quadratic, exact from order 2 (M = V (x_c - c) x (-g)).
!<
!< Several bodies (#109 WP3 step 02): `load_samples_3d_bodies(surfaces, ...)` / `load_samples_2d_bodies(curves, ...)`
!< give one sample set of an array of bodies, each at its pose: the samples of every body as `load_samples_3d/_2d` gives
!< them (the same kernels), concatenated in body order (element = the facet / segment of its own surface), `body` = the
!< surface (curve) index and `sub_body` = the element's facet_body (loop body); a probe is visible when its segment
!< crosses no element of ANY body (its own as the single-body set decides it, every other body by the closed segment
!< point - probe in that body's frame: `blocked_3d` / `blocked_2d` of the closure of several bodies, the 3-D tolerance
!< the sample's own body's); `clearance` the unsigned distance from the sample point to the nearest OTHER surface
!< (huge if none; the shells of one surface are not separated), `second_body` that surface, `second_foot` its closest
!< point (world). Ties keep the lower surface index. With one body the set equals `load_samples_3d/_2d` bit for bit in
!< every single-body field when the surface has one facet_body (body = 1 = its facet_body). `integrate_loads` sums
!< `per_body` by the `body` field: the facet_body of a single-body set, the surface of a bodies set (the per-shell
!< totals of one surface: a single-body set of it). Both paths bit for bit (the merge kernels are copies and
!< comparisons; the device runs one launch set per body, then the merges in body order).
use fossil_aabb_node_object,   only : aabb_node_object
use fossil_aabb_tree_object,   only : facet_distance_payload, AABB_TREE_SAH_BVH
use fossil_curve_object,       only : curve_object
use fossil_fd_closure,         only : fd_cell_index_oac, to_world_3d_oac, rotate_3d_oac, to_world_2d_oac, rotate_2d_oac, &
                                      blocked_3d, blocked_2d, merge_clear, body_points_3d_dev, body_points_2d_dev, &
                                      world_points_3d_dev, world_points_2d_dev
use fossil_rigid_pose,         only : rigid_pose_t, rigid_pose_2d_t
use fossil_segment_visibility, only : bvh_segment_blocked_oac
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use iso_c_binding,             only : c_intptr_t, c_loc
use iso_fortran_env,           only : stderr => error_unit, int64
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
private
public :: load_samples_t, load_samples_3d, load_samples_2d, load_samples_reset, load_samples_mismatches, load_rule_size
public :: load_samples_3d_bodies, load_samples_2d_bodies
public :: integrate_loads
! the per-point kernels, shared by the two paths (public for the tests)
public :: load_point_3d_oac, load_point_2d_oac, load_force_3d_oac, load_force_2d_oac
public :: rule_3d_oac ! the 3-D rules, for `fossil_field_transfer`

interface integrate_loads
   !< Force, moment, per patch and per body totals of the solver's p and tau at the load samples (3-D: moment(3); 2-D: a
   !< scalar moment).
   module procedure integrate_loads_3d, integrate_loads_2d
endinterface integrate_loads

type :: load_samples_t
   !< The load samples of a surface or curve (see the module header). dims = 3 or 2 (2-D: the curve's (u, v) plane). Every
   !< array is allocated after a call (zero-size when not asked: the probes with m = 0, probe_cell without lines).
   integer(I4P)              :: n = 0_I4P             !< Samples (q per element).
   integer(I4P)              :: dims = 0_I4P          !< Dimensions (2 or 3).
   integer(I4P)              :: order = 0_I4P         !< Rule order (1, 2, 3).
   integer(I4P)              :: q = 0_I4P             !< Points per element.
   integer(I4P)              :: ne = 0_I4P            !< Elements (facets, segments = vertices).
   integer(I4P)              :: m = 0_I4P             !< Probes per sample.
   integer(I4P)              :: ni = 0_I4P            !< Points of the x (u) line (0: no lines).
   integer(I4P)              :: nj = 0_I4P            !< Points of the y (v) line.
   integer(I4P)              :: nk = 0_I4P            !< Points of the z line (0 in 2-D).
   logical                   :: with_cells = .false.  !< probe_cell filled (lines given).
   logical                   :: with_bodies = .false. !< The several-bodies fields filled (`load_samples_*_bodies`).
   integer(I4P)              :: bodies = 0_I4P        !< Surfaces (curves) of a bodies set (0: single body).
   real(R8P),    allocatable :: points(:,:)           !< Sample points [dims, n].
   real(R8P),    allocatable :: weight(:)             !< Weights [n] (rule weight times area / length).
   real(R8P),    allocatable :: normal(:,:)           !< Element unit normal [dims, n] (outward; 0 if degenerate).
   integer(I4P), allocatable :: element(:)            !< Facet / segment [n].
   integer(I4P), allocatable :: patch(:)              !< Patch of the element [n].
   integer(I4P), allocatable :: body(:)               !< Body of the element [n] (bodies set: the surface index).
   real(R8P),    allocatable :: lengths(:)            !< Probe lengths [m] (+ outward, absolute).
   real(R8P),    allocatable :: probe(:,:,:)          !< Probe points [dims, m, n].
   logical,      allocatable :: visible(:,:)          !< The segment point-probe crosses no element [m, n].
   integer(I4P), allocatable :: probe_cell(:,:,:)     !< Enclosing cell per axis [dims, m, n] (0: outside).
   ! several bodies (`load_samples_*_bodies`; zero-size in a single-body set): `body` is then the surface (curve) index
   integer(I4P), allocatable :: sub_body(:)           !< The element's facet_body (2-D: loop body) [n].
   real(R8P),    allocatable :: clearance(:)          !< Unsigned distance point - nearest other surface [n] (huge: none).
   integer(I4P), allocatable :: second_body(:)        !< That nearest other surface [n] (0: none).
   real(R8P),    allocatable :: second_foot(:,:)      !< Its closest point to the sample point [dims, n].
   logical                   :: on_device = .false.   !< The arrays are mapped on the device (resident).
   integer(c_intptr_t)       :: device_key = 0        !< Host address of points when mapped (a copy has another one).
   contains
      procedure, pass(self) :: enter_device !< Map the samples to the device (copyin).
      procedure, pass(self) :: update_host  !< Copy the mapped samples back to the host.
      procedure, pass(self) :: exit_device  !< Release the device copies.
endtype load_samples_t

! 3-D rules: barycentric points and weights (summing to 1)
real(R8P), parameter :: T1 = 0.33333333333333333333_R8P  !< Order 1: the centroid.
real(R8P), parameter :: T2A = 0.66666666666666666667_R8P !< Order 2: 2/3.
real(R8P), parameter :: T2B = 0.16666666666666666667_R8P !< Order 2: 1/6.
real(R8P), parameter :: T3A1 = 0.44594849091596488632_R8P  !< Order 3 (Dunavant degree 4): first orbit, a.
real(R8P), parameter :: T3B1 = 0.10810301816807022736_R8P  !< First orbit, 1 - 2a.
real(R8P), parameter :: T3W1 = 0.22338158967801146570_R8P  !< First orbit, weight.
real(R8P), parameter :: T3A2 = 0.091576213509770743460_R8P !< Second orbit, a.
real(R8P), parameter :: T3B2 = 0.81684757298045851308_R8P  !< Second orbit, 1 - 2a.
real(R8P), parameter :: T3W2 = 0.10995174365532186764_R8P  !< Second orbit, weight.
! 2-D rules: Gauss-Legendre on [0, 1]
real(R8P), parameter :: G2A = 0.21132486540518711775_R8P !< 2 points: 1/2 - sqrt(3)/6.
real(R8P), parameter :: G2B = 0.78867513459481288225_R8P !< 2 points: 1/2 + sqrt(3)/6.
real(R8P), parameter :: G3A = 0.11270166537925831148_R8P !< 3 points: 1/2 - sqrt(15)/10.
real(R8P), parameter :: G3B = 0.88729833462074168852_R8P !< 3 points: 1/2 + sqrt(15)/10.
real(R8P), parameter :: G3W1 = 0.27777777777777777778_R8P !< 3 points: 5/18.
real(R8P), parameter :: G3W2 = 0.44444444444444444444_R8P !< 3 points: 8/18.

interface mis_i
   module procedure mis_i1, mis_i3
endinterface mis_i
interface mis_r
   module procedure mis_r1, mis_r2, mis_r3
endinterface mis_r

contains
   pure integer(I4P) function load_rule_size(order, dims) result(q)
   !< Points per element of the rule of order `order` (1, 2, 3) in `dims` dimensions: 3-D 1, 3, 6; 2-D 1, 2, 3; 0 if the
   !< order or the dimensions are not supported.
   integer(I4P), intent(in) :: order !< Rule order.
   integer(I4P), intent(in) :: dims  !< Dimensions (2 or 3).

   q = 0_I4P
   if (order < 1_I4P .or. order > 3_I4P) return
   if (dims == 3_I4P) then
      q = 1_I4P ; if (order == 2_I4P) q = 3_I4P ; if (order == 3_I4P) q = 6_I4P
   elseif (dims == 2_I4P) then
      q = order
   endif
   endfunction load_rule_size

   subroutine load_samples_3d(surface, x, y, z, samples, status, order, lengths, pose, device, resident)
   !< The load samples of the surface (see the module header). `x, y, z` (optional, all or none): the solver's lines, for
   !< the cells of the probes. `order` (default 2): the rule, 1, 2 or 3. `lengths(m)`: probe lengths (absolute, + outward;
   !< default none; they need the SAH BVH). `pose`: the body at the pose. `device` (default .false.): the device path
   !< (`enter_device` first); `resident` (default .false., implies device): leave the samples mapped.
   !< `status = STATUS_INVALID_INPUT` on a bad order, lines given in part, with fewer than two points or not increasing,
   !< probes on an octree surface, or (device) a surface with facets not mapped (a stderr message); the samples are then
   !< reset (n = 0). A surface without facets gives no sample.
   type(surface_stl_object), intent(in), target    :: surface          !< Surface.
   real(R8P),                intent(in),  optional :: x(:), y(:), z(:) !< The solver's lines (strictly increasing).
   type(load_samples_t),     intent(inout)         :: samples          !< Samples (released and reset first).
   integer(I4P),             intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),             intent(in),  optional :: order            !< Rule order (default 2).
   real(R8P),                intent(in),  optional :: lengths(:)       !< Probe lengths [m] (+ outward, absolute).
   type(rigid_pose_t),       intent(in),  optional :: pose             !< Pose of the body.
   logical,                  intent(in),  optional :: device           !< Device path (default .false.).
   logical,                  intent(in),  optional :: resident         !< Leave the samples mapped (default .false.).
   real(R8P), allocatable :: ll(:), xl(:), yl(:), zl(:) !< Probe lengths, lines.
   real(R8P)              :: pr(3, 3), pt(3)            !< Pose.
   integer(I4P)           :: st, ord, nf                !< Status, order, facets.
   logical                :: dev, res, wl               !< Options.

   call load_samples_reset(samples)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   ord = 2_I4P ; if (present(order)) ord = order
   call resolve_inputs(3_I4P, ord, lengths, x, y, z, st, ll, wl, xl, yl, zl)
   nf = surface%get_facets_number()
   if (st == STATUS_OK .and. nf > 0_I4P) then
      if (.not. (allocated(surface%tri_vertices) .and. allocated(surface%body_payload))) then
         write(stderr, '(A)') 'fossil load_samples_3d: refused -- the surface is not analyzed'
         st = STATUS_INVALID_INPUT
      elseif (size(ll) > 0 .and. .not. (surface%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. &
                                         allocated(surface%aabb%node) .and. allocated(surface%aabb%payload))) then
         write(stderr, '(A)') 'fossil load_samples_3d: refused -- the probes need the SAH BVH'
         st = STATUS_INVALID_INPUT
      elseif (dev .and. .not. surface%is_on_device) then
         write(stderr, '(A)') 'fossil load_samples_3d: device path refused -- call surface%enter_device first'
         st = STATUS_INVALID_INPUT
      endif
   endif
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   pr = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ; pt = 0._R8P
   if (present(pose)) then
      pr = pose%r ; pt = pose%t
   endif
   call samples_3d_run(surface, ord, ll, present(pose), pr, pt, wl, xl, yl, zl, dev .and. nf > 0_I4P, res, samples)
   endsubroutine load_samples_3d

   subroutine load_samples_2d(curve, x, y, samples, status, order, lengths, pose, device, resident)
   !< The load samples of the curve in its (u, v) plane (see the module header): x (u), y (v) the solver's lines (optional,
   !< both or none); one sample set per segment i -> next(i), i = 1 ... vertices (the missing segment at the end of an open
   !< loop is degenerate). The other arguments and the statuses as `load_samples_3d` (device: the curve mapped by
   !< `enter_device`).
   type(curve_object),    intent(in)            :: curve      !< Curve.
   real(R8P),             intent(in),  optional :: x(:), y(:) !< The solver's lines (u, v; strictly increasing).
   type(load_samples_t),  intent(inout)         :: samples    !< Samples (released and reset first).
   integer(I4P),          intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),          intent(in),  optional :: order      !< Rule order (default 2).
   real(R8P),             intent(in),  optional :: lengths(:) !< Probe lengths [m] (+ outward, absolute).
   type(rigid_pose_2d_t), intent(in),  optional :: pose       !< Pose of the curve.
   logical,               intent(in),  optional :: device     !< Device path (default .false.).
   logical,               intent(in),  optional :: resident   !< Leave the samples mapped (default .false.).
   real(R8P), allocatable :: ll(:), xl(:), yl(:), zl(:) !< Probe lengths, lines.
   real(R8P)              :: pr(2, 2), pt(2)            !< Pose.
   integer(I4P)           :: st, ord, nv                !< Status, order, vertices.
   logical                :: dev, res, wl               !< Options.

   call load_samples_reset(samples)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   ord = 2_I4P ; if (present(order)) ord = order
   call resolve_inputs(2_I4P, ord, lengths, x, y, st=st, ll=ll, wl=wl, xl=xl, yl=yl, zl=zl)
   nv = curve%get_vertices_number()
   if (st == STATUS_OK .and. dev .and. nv > 0_I4P .and. .not. curve%get_is_on_device()) then
      write(stderr, '(A)') 'fossil load_samples_2d: device path refused -- call curve%enter_device first'
      st = STATUS_INVALID_INPUT
   endif
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   pr = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) ; pt = 0._R8P
   if (present(pose)) then
      pr = pose%r ; pt = pose%t
   endif
   call samples_2d_run(curve, ord, ll, present(pose), pr, pt, wl, xl, yl, dev .and. nv > 0_I4P, res, samples)
   endsubroutine load_samples_2d

   subroutine load_samples_3d_bodies(surfaces, x, y, z, samples, status, poses, order, lengths, device, resident)
   !< The load samples of several bodies (#109 WP3 step 02; see the module header, "Several bodies"): surfaces(b) at
   !< poses(b) (absent: no body posed), the samples of every body concatenated in body order, body = b, sub_body = the
   !< facet_body; the probes visible against every body; clearance, second_body, second_foot. The other arguments as
   !< `load_samples_3d`. Every non-empty surface needs its SAH BVH (the default tree: the visibility and the clearance
   !< against the other bodies traverse it) and (device path) to be mapped by `enter_device`. `status =
   !< STATUS_INVALID_INPUT` as `load_samples_3d`, and on no surface, size(poses) /= size(surfaces), a surface with an
   !< octree (a stderr message); the samples are then reset (n = 0). Empty surfaces give no sample.
   type(surface_stl_object), intent(in), target    :: surfaces(:)      !< Surfaces (the bodies).
   real(R8P),                intent(in),  optional :: x(:), y(:), z(:) !< The solver's lines (strictly increasing).
   type(load_samples_t),     intent(inout)         :: samples          !< Samples (released and reset first).
   integer(I4P),             intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_t),       intent(in),  optional :: poses(:)         !< Poses of the bodies [nb].
   integer(I4P),             intent(in),  optional :: order            !< Rule order (default 2).
   real(R8P),                intent(in),  optional :: lengths(:)       !< Probe lengths [m] (+ outward, absolute).
   logical,                  intent(in),  optional :: device           !< Device path (default .false.).
   logical,                  intent(in),  optional :: resident         !< Leave the samples mapped (default .false.).
   real(R8P),    allocatable :: ll(:), xl(:), yl(:), zl(:)      !< Probe lengths, lines.
   real(R8P),    allocatable :: prs(:,:,:), pts(:,:), tols(:)   !< Poses, visibility tolerances per body.
   integer(I4P), allocatable :: nfs(:)                          !< Facets per body.
   real(R8P)                 :: dbx, dby, dbz                   !< Box.
   type(vector_R8P)          :: bmin, bmax                      !< Box.
   integer(I4P)              :: st, ord, nb, b                  !< Work.
   logical                   :: dev, res, wl, posed, some       !< Options.

   call load_samples_reset(samples)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   ord = 2_I4P ; if (present(order)) ord = order
   posed = present(poses)
   nb = size(surfaces, kind=I4P)
   call resolve_inputs(3_I4P, ord, lengths, x, y, z, st, ll, wl, xl, yl, zl)
   if (nb < 1_I4P) st = STATUS_INVALID_INPUT
   if (posed) then
      if (size(poses) /= nb) st = STATUS_INVALID_INPUT
   endif
   allocate(prs(3, 3, max(nb, 1_I4P)), pts(3, max(nb, 1_I4P)), tols(max(nb, 1_I4P)), nfs(max(nb, 1_I4P)))
   some = .false.
   do b = 1_I4P, nb
      prs(:, :, b) = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3])
      pts(:, b) = 0._R8P
      if (posed .and. st == STATUS_OK) then
         prs(:, :, b) = poses(b)%r ; pts(:, b) = poses(b)%t
      endif
      nfs(b) = surfaces(b)%get_facets_number()
      bmin = surfaces(b)%get_bmin() ; bmax = surfaces(b)%get_bmax()
      dbx = bmax%x - bmin%x ; dby = bmax%y - bmin%y ; dbz = bmax%z - bmin%z
      tols(b) = 1.e-9_R8P * sqrt(dbx * dbx + dby * dby + dbz * dbz) ! the tolerance of probe_many
      if (nfs(b) < 1_I4P .or. st /= STATUS_OK) cycle
      some = .true.
      if (.not. (allocated(surfaces(b)%tri_vertices) .and. allocated(surfaces(b)%body_payload))) then
         write(stderr, '(A,I0,A)') 'fossil load_samples_3d_bodies: refused -- surface ', b, ' is not analyzed'
         st = STATUS_INVALID_INPUT
      elseif (.not. (surfaces(b)%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(surfaces(b)%aabb%node) .and. &
                     allocated(surfaces(b)%aabb%payload))) then
         write(stderr, '(A,I0,A)') 'fossil load_samples_3d_bodies: refused -- surface ', b, ' has no SAH BVH'
         st = STATUS_INVALID_INPUT
      elseif (dev .and. .not. surfaces(b)%is_on_device) then
         write(stderr, '(A,I0,A)') 'fossil load_samples_3d_bodies: device path refused -- call surface%enter_device '// &
                                   'first (surface ', b, ')'
         st = STATUS_INVALID_INPUT
      endif
   enddo
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (posed) then
      call bodies_3d_run(surfaces, nfs, ord, ll, posed, prs, pts, tols, wl, xl, yl, zl, dev .and. some, res, samples, poses)
   else
      call bodies_3d_run(surfaces, nfs, ord, ll, posed, prs, pts, tols, wl, xl, yl, zl, dev .and. some, res, samples)
   endif
   endsubroutine load_samples_3d_bodies

   subroutine load_samples_2d_bodies(curves, x, y, samples, status, poses, order, lengths, device, resident)
   !< 2-D twin of `load_samples_3d_bodies` (#109 WP3 step 02): curves(b) at poses(b) in their common (u, v) plane, the
   !< samples of `load_samples_2d` per curve, sub_body = the loop body. The device path needs every non-empty curve mapped
   !< by `enter_device`. The other arguments and the statuses as `load_samples_3d_bodies` (no tree requirement: every
   !< curve has its segment BVH).
   type(curve_object),    intent(in), target    :: curves(:)  !< Curves (the bodies).
   real(R8P),             intent(in),  optional :: x(:), y(:) !< The solver's lines (u, v; strictly increasing).
   type(load_samples_t),  intent(inout)         :: samples    !< Samples (released and reset first).
   integer(I4P),          intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_2d_t), intent(in),  optional :: poses(:)   !< Poses of the curves [nb].
   integer(I4P),          intent(in),  optional :: order      !< Rule order (default 2).
   real(R8P),             intent(in),  optional :: lengths(:) !< Probe lengths [m] (+ outward, absolute).
   logical,               intent(in),  optional :: device     !< Device path (default .false.).
   logical,               intent(in),  optional :: resident   !< Leave the samples mapped (default .false.).
   real(R8P),    allocatable :: ll(:), xl(:), yl(:), zl(:)    !< Probe lengths, lines.
   real(R8P),    allocatable :: prs(:,:,:), pts(:,:)          !< Poses per curve.
   integer(I4P), allocatable :: nvs(:)                        !< Vertices per curve.
   integer(I4P)              :: st, ord, nb, b                !< Work.
   logical                   :: dev, res, wl, posed, some     !< Options.

   call load_samples_reset(samples)
   if (present(status)) status = STATUS_OK
   res = .false. ; if (present(resident)) res = resident
   dev = res ; if (present(device)) dev = dev .or. device
   ord = 2_I4P ; if (present(order)) ord = order
   posed = present(poses)
   nb = size(curves, kind=I4P)
   call resolve_inputs(2_I4P, ord, lengths, x, y, st=st, ll=ll, wl=wl, xl=xl, yl=yl, zl=zl)
   if (nb < 1_I4P) st = STATUS_INVALID_INPUT
   if (posed) then
      if (size(poses) /= nb) st = STATUS_INVALID_INPUT
   endif
   allocate(prs(2, 2, max(nb, 1_I4P)), pts(2, max(nb, 1_I4P)), nvs(max(nb, 1_I4P)))
   some = .false.
   do b = 1_I4P, nb
      prs(:, :, b) = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P], [2, 2]) ; pts(:, b) = 0._R8P
      if (posed .and. st == STATUS_OK) then
         prs(:, :, b) = poses(b)%r ; pts(:, b) = poses(b)%t
      endif
      nvs(b) = curves(b)%get_vertices_number()
      if (nvs(b) < 1_I4P .or. st /= STATUS_OK) cycle
      some = .true.
      if (dev .and. .not. curves(b)%get_is_on_device()) then
         write(stderr, '(A,I0,A)') 'fossil load_samples_2d_bodies: device path refused -- call curve%enter_device '// &
                                   'first (curve ', b, ')'
         st = STATUS_INVALID_INPUT
      endif
   enddo
   if (st /= STATUS_OK) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (posed) then
      call bodies_2d_run(curves, nvs, ord, ll, posed, prs, pts, wl, xl, yl, dev .and. some, res, samples, poses)
   else
      call bodies_2d_run(curves, nvs, ord, ll, posed, prs, pts, wl, xl, yl, dev .and. some, res, samples)
   endif
   endsubroutine load_samples_2d_bodies

   subroutine load_samples_reset(samples)
   !< Release resident samples from the device (`exit_device`), then reset them to the default (no sample, arrays
   !< released): what `intent(out)` did, without freeing arrays that are still mapped.
   type(load_samples_t), intent(inout) :: samples !< Samples.
   type(load_samples_t)                :: fresh   !< Default.

   if (samples%on_device) call samples%exit_device
   samples = fresh
   endsubroutine load_samples_reset

   function load_samples_mismatches(a, b) result(nmis)
   !< Number of differing entries of two sample sets, every field, reals and logicals by bits (the twin checks of the tests
   !< and benches). Arrays of different shapes count their larger size. Read on the host: `update_host` resident samples
   !< first.
   type(load_samples_t), intent(in) :: a, b !< Samples.
   integer(I4P)                     :: nmis !< Mismatches.

   nmis = 0_I4P
   if (a%n /= b%n) nmis = nmis + 1_I4P
   if (a%dims /= b%dims .or. a%order /= b%order .or. a%q /= b%q .or. a%ne /= b%ne) nmis = nmis + 1_I4P
   if (a%m /= b%m) nmis = nmis + 1_I4P
   if (a%ni /= b%ni .or. a%nj /= b%nj .or. a%nk /= b%nk) nmis = nmis + 1_I4P
   if (a%with_cells .neqv. b%with_cells) nmis = nmis + 1_I4P
   nmis = nmis + mis_r(a%points, b%points) + mis_r(a%weight, b%weight) + mis_r(a%normal, b%normal)
   nmis = nmis + mis_i(a%element, b%element) + mis_i(a%patch, b%patch) + mis_i(a%body, b%body)
   nmis = nmis + mis_r(a%lengths, b%lengths) + mis_r(a%probe, b%probe) + mis_l(a%visible, b%visible)
   nmis = nmis + mis_i(a%probe_cell, b%probe_cell)
   if ((a%with_bodies .neqv. b%with_bodies) .or. a%bodies /= b%bodies) nmis = nmis + 1_I4P
   nmis = nmis + mis_i(a%sub_body, b%sub_body) + mis_r(a%clearance, b%clearance) + mis_i(a%second_body, b%second_body)
   nmis = nmis + mis_r(a%second_foot, b%second_foot)
   endfunction load_samples_mismatches

   ! the loads

   subroutine integrate_loads_3d(samples, force, moment, status, pressure, shear, about, per_patch_force, &
                                 per_patch_moment, per_body_force, per_body_moment, area, device)
   !< Force and moment on the wall of 3-D load samples (see the module header): with the pressure p(n) and/or the viscous
   !< traction tau(3, n) at the samples (at least one), F = sum_s w_s (-p_s n_s + tau_s), M = sum_s (x_s - about) x f_s
   !< (`about` default the origin), in sample order with Neumaier compensation. Optional: the totals per patch
   !< `per_patch_force(3, np)`, `per_patch_moment(3, np)` and per body `per_body_force(3, nb)`, `per_body_moment(3, nb)`
   !< (by the samples' `body`: the facet_body of `load_samples_3d`, the surface index of `load_samples_3d_bodies`; np, nb
   !< the largest patch and body ids of the samples; a sample with an id below 1 enters the totals only), the
   !< wetted area sum_s w_s. `device` (default .false.): the contributions computed on the device, the samples and p, tau
   !< read on the device where they are mapped (resident samples, a pressure mapped by the caller), else copied in; the
   !< host path reads the host arrays (`update_host` resident samples first). Both paths give the same bits.
   !< `status = STATUS_INVALID_INPUT` (stderr message, outputs zero, per-patch / per-body arrays of size 0) if the samples
   !< are not 3-D, neither p nor tau is given, or a size differs (p(n), tau(3, n), force(3), moment(3), about(3)).
   type(load_samples_t),   intent(in)            :: samples                !< Load samples (3-D).
   real(R8P),              intent(out)           :: force(:)               !< Force [3].
   real(R8P),              intent(out)           :: moment(:)              !< Moment about `about` [3].
   integer(I4P),           intent(out), optional :: status                 !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),              intent(in),  optional :: pressure(:)            !< Pressure at the samples [n].
   real(R8P),              intent(in),  optional :: shear(:,:)             !< Viscous traction at the samples [3, n].
   real(R8P),              intent(in),  optional :: about(:)               !< Point of the moments [3] (default origin).
   real(R8P), allocatable, intent(out), optional :: per_patch_force(:,:)   !< Force per patch [3, np].
   real(R8P), allocatable, intent(out), optional :: per_patch_moment(:,:)  !< Moment per patch [3, np].
   real(R8P), allocatable, intent(out), optional :: per_body_force(:,:)    !< Force per body [3, nb].
   real(R8P), allocatable, intent(out), optional :: per_body_moment(:,:)   !< Moment per body [3, nb].
   real(R8P),              intent(out), optional :: area                   !< Wetted area sum_s w_s.
   logical,                intent(in),  optional :: device                 !< Device path (default .false.).
   real(R8P), allocatable :: tot(:), pp(:,:), pb(:,:) !< Totals: all, per patch, per body [7, ...].
   real(R8P)              :: ab(3)                    !< About.
   integer(I4P)           :: st                       !< Status.
   logical                :: dev                      !< Device path.

   dev = .false. ; if (present(device)) dev = device
   ab = 0._R8P
   st = loads_check(samples, 3_I4P, present(pressure), present(shear), size(force, kind=I4P), &
                    size(moment, kind=I4P) == 3_I4P)
   if (st == STATUS_OK .and. present(pressure)) then
      if (size(pressure, kind=I4P) /= samples%n) st = STATUS_INVALID_INPUT
   endif
   if (st == STATUS_OK .and. present(shear)) then
      if (size(shear, 1, kind=I4P) /= 3_I4P .or. size(shear, 2, kind=I4P) /= samples%n) st = STATUS_INVALID_INPUT
   endif
   if (st == STATUS_OK .and. present(about)) then
      if (size(about, kind=I4P) /= 3_I4P) then
         st = STATUS_INVALID_INPUT
      else
         ab = about
      endif
   endif
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil integrate_loads: refused -- 3-D samples, pressure and/or shear of the samples'' sizes'
      allocate(tot(7), pp(7, 0), pb(7, 0)) ; tot = 0._R8P
   else
      call loads_run(samples, 3_I4P, 7_I4P, pressure, shear, ab, dev, tot, pp, pb)
   endif
   force = 0._R8P ; moment = 0._R8P
   if (size(force) == 3) force = tot(1:3)
   if (size(moment) == 3) moment = tot(4:6)
   if (present(per_patch_force)) per_patch_force = pp(1:3, :)
   if (present(per_patch_moment)) per_patch_moment = pp(4:6, :)
   if (present(per_body_force)) per_body_force = pb(1:3, :)
   if (present(per_body_moment)) per_body_moment = pb(4:6, :)
   if (present(area)) area = tot(7)
   if (present(status)) status = st
   endsubroutine integrate_loads_3d

   subroutine integrate_loads_2d(samples, force, moment, status, pressure, shear, about, per_patch_force, &
                                 per_patch_moment, per_body_force, per_body_moment, area, device)
   !< 2-D twin of `integrate_loads_3d` per unit depth in the curve's (u, v) plane: tau(2, n), force(2), about(2), the
   !< scalar moment M = sum_s (r_u f_v - r_v f_u), r = x_s - about (its component along the plane's normal axis);
   !< `per_patch_moment(np)`, `per_body_moment(nb)`; the wetted length as `area`.
   type(load_samples_t),   intent(in)            :: samples                !< Load samples (2-D).
   real(R8P),              intent(out)           :: force(:)               !< Force [2].
   real(R8P),              intent(out)           :: moment                 !< Moment about `about`.
   integer(I4P),           intent(out), optional :: status                 !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),              intent(in),  optional :: pressure(:)            !< Pressure at the samples [n].
   real(R8P),              intent(in),  optional :: shear(:,:)             !< Viscous traction at the samples [2, n].
   real(R8P),              intent(in),  optional :: about(:)               !< Point of the moments [2] (default origin).
   real(R8P), allocatable, intent(out), optional :: per_patch_force(:,:)   !< Force per patch [2, np].
   real(R8P), allocatable, intent(out), optional :: per_patch_moment(:)    !< Moment per patch [np].
   real(R8P), allocatable, intent(out), optional :: per_body_force(:,:)    !< Force per body [2, nb].
   real(R8P), allocatable, intent(out), optional :: per_body_moment(:)     !< Moment per body [nb].
   real(R8P),              intent(out), optional :: area                   !< Wetted length sum_s w_s.
   logical,                intent(in),  optional :: device                 !< Device path (default .false.).
   real(R8P), allocatable :: tot(:), pp(:,:), pb(:,:) !< Totals: all, per patch, per body [4, ...].
   real(R8P)              :: ab(2)                    !< About.
   integer(I4P)           :: st                       !< Status.
   logical                :: dev                      !< Device path.

   dev = .false. ; if (present(device)) dev = device
   ab = 0._R8P
   st = loads_check(samples, 2_I4P, present(pressure), present(shear), size(force, kind=I4P), .true.)
   if (st == STATUS_OK .and. present(pressure)) then
      if (size(pressure, kind=I4P) /= samples%n) st = STATUS_INVALID_INPUT
   endif
   if (st == STATUS_OK .and. present(shear)) then
      if (size(shear, 1, kind=I4P) /= 2_I4P .or. size(shear, 2, kind=I4P) /= samples%n) st = STATUS_INVALID_INPUT
   endif
   if (st == STATUS_OK .and. present(about)) then
      if (size(about, kind=I4P) /= 2_I4P) then
         st = STATUS_INVALID_INPUT
      else
         ab = about
      endif
   endif
   if (st /= STATUS_OK) then
      write(stderr, '(A)') 'fossil integrate_loads: refused -- 2-D samples, pressure and/or shear of the samples'' sizes'
      allocate(tot(4), pp(4, 0), pb(4, 0)) ; tot = 0._R8P
   else
      call loads_run(samples, 2_I4P, 4_I4P, pressure, shear, ab, dev, tot, pp, pb)
   endif
   force = 0._R8P
   if (size(force) == 2) force = tot(1:2)
   moment = tot(3)
   if (present(per_patch_force)) per_patch_force = pp(1:2, :)
   if (present(per_patch_moment)) per_patch_moment = pp(3, :)
   if (present(per_body_force)) per_body_force = pb(1:2, :)
   if (present(per_body_moment)) per_body_moment = pb(3, :)
   if (present(area)) area = tot(4)
   if (present(status)) status = st
   endsubroutine integrate_loads_2d

   ! the resident pattern

   subroutine enter_device(self)
   !< Map host samples to the device (copyin). Member arrays only, never `self` (#20). No-op if already mapped or if there
   !< is no sample; a copy of resident samples is mapped as new samples.
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc, gd                    !< Non-empty groups.

   if (is_mapped(self)) return
   self%on_device = .false. ; self%device_key = 0
   if (self%n < 1_I4P .or. .not. allocated(self%points)) return
   call groups(self, ga, gb, gc, gd)
   !$acc enter data copyin(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
   !$acc enter data copyin(self%lengths, self%probe, self%visible) if(gb)
   !$acc enter data copyin(self%probe_cell) if(gc)
   !$acc enter data copyin(self%sub_body, self%clearance, self%second_body, self%second_foot) if(gd)
   call mark_resident(self)
   endsubroutine enter_device

   subroutine update_host(self)
   !< Copy the mapped samples back to the host; they stay mapped. No-op if not mapped; on a copy of resident samples a
   !< warning on stderr (its arrays are not on the device: update the original before copying it).
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc, gd                    !< Non-empty groups.

   if (.not. self%on_device) return
   if (.not. is_mapped(self)) then
      write(stderr, '(A)') 'fossil load_samples_t%update_host: a copy of resident samples is not on the device; '// &
                           'call update_host on the original before copying it'
      self%on_device = .false. ; self%device_key = 0
      return
   endif
   call groups(self, ga, gb, gc, gd)
   !$acc update host(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
   !$acc update host(self%lengths, self%probe, self%visible) if(gb)
   !$acc update host(self%probe_cell) if(gc)
   !$acc update host(self%sub_body, self%clearance, self%second_body, self%second_foot) if(gd)
   endsubroutine update_host

   subroutine exit_device(self)
   !< Release the device copies; the host arrays are kept as they are (call `update_host` first to keep the device
   !< values). No-op if not mapped; a copy of resident samples only drops the flag.
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc, gd                    !< Non-empty groups.

   if (.not. self%on_device) return
   if (is_mapped(self)) then
      call groups(self, ga, gb, gc, gd)
      !$acc exit data delete(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
      !$acc exit data delete(self%lengths, self%probe, self%visible) if(gb)
      !$acc exit data delete(self%probe_cell) if(gc)
      !$acc exit data delete(self%sub_body, self%clearance, self%second_body, self%second_foot) if(gd)
   endif
   self%on_device = .false. ; self%device_key = 0
   endsubroutine exit_device

   ! private: the two runs (host and device in one set of loops)

   subroutine samples_3d_run(surface, ord, ll, posed, pr, pt, wl, xl, yl, zl, dev, res, s)
   !< `load_samples_3d` after the checks: the samples built in a local set (mapped on the device path), then moved into s.
   type(surface_stl_object), intent(in), target :: surface          !< Surface (analyzed; mapped on the device path).
   integer(I4P),             intent(in)         :: ord              !< Rule order.
   real(R8P),                intent(in)         :: ll(:)            !< Probe lengths.
   logical,                  intent(in)         :: posed            !< A pose.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3)  !< Pose.
   logical,                  intent(in)         :: wl               !< Lines given.
   real(R8P),                intent(in)         :: xl(:), yl(:), zl(:) !< Lines (one dummy point each without).
   logical,                  intent(in)         :: dev, res         !< Device path, resident.
   type(load_samples_t),     intent(inout)      :: s                !< Samples.
   type(load_samples_t)   :: w                    !< The samples being built.
   integer(I4P), allocatable :: fpatch(:)         !< Facet patches.
   real(R8P)              :: tol, dbx, dby, dbz   !< Visibility tolerance.
   type(vector_R8P)       :: bmin, bmax           !< Surface box.
   integer(I4P)           :: nf, q, n, m, f, ni, nj, nk !< Work.

   nf = surface%get_facets_number() ; q = load_rule_size(ord, 3_I4P) ; n = nf * q ; m = size(ll, kind=I4P)
   ni = 0_I4P ; nj = 0_I4P ; nk = 0_I4P
   if (wl) then
      ni = size(xl, kind=I4P) ; nj = size(yl, kind=I4P) ; nk = size(zl, kind=I4P)
   endif
   call samples_arrays(w, 3_I4P, ord, q, nf, m, wl, ni, nj, nk)
   w%lengths = ll
   if (dev) call map_create(w)
   if (n > 0_I4P) then
      allocate(fpatch(nf))
      do f = 1_I4P, nf
         fpatch(f) = surface%facet_patch(f)
      enddo
      bmin = surface%get_bmin() ; bmax = surface%get_bmax()
      dbx = bmax%x - bmin%x ; dby = bmax%y - bmin%y ; dbz = bmax%z - bmin%z
      tol = 1.e-9_R8P * sqrt(dbx * dbx + dby * dby + dbz * dbz) ! the tolerance of probe_many
      !$acc data copyin(fpatch, xl, yl, zl) if(dev)
      call points_3d_stage(surface, n, q, ord, nf, fpatch, posed, pr, pt, dev, w%points, w%weight, w%normal, w%element, &
                           w%patch, w%body)
      if (m > 0_I4P) call probes_3d_stage(surface, n, q, ord, m, w%lengths, tol, posed, pr, pt, dev, w%probe, w%visible)
      if (m > 0_I4P .and. wl) call cells_stage(3_I4P, n, m, ni, nj, nk, xl, yl, zl, dev, w%probe, w%probe_cell)
      !$acc end data
   endif
   call finish(w, dev, res, s)
   endsubroutine samples_3d_run

   subroutine samples_2d_run(curve, ord, ll, posed, pr, pt, wl, xl, yl, dev, res, s)
   !< `load_samples_2d` after the checks: the curve's vertices, next, segment bodies and patches gathered on the host, the
   !< samples built in a local set (mapped on the device path), then moved into s.
   type(curve_object),   intent(in)    :: curve           !< Curve (mapped on the device path).
   integer(I4P),         intent(in)    :: ord             !< Rule order.
   real(R8P),            intent(in)    :: ll(:)           !< Probe lengths.
   logical,              intent(in)    :: posed           !< A pose.
   real(R8P),            intent(in)    :: pr(2, 2), pt(2) !< Pose.
   logical,              intent(in)    :: wl              !< Lines given.
   real(R8P),            intent(in)    :: xl(:), yl(:)    !< Lines (one dummy point each without).
   logical,              intent(in)    :: dev, res        !< Device path, resident.
   type(load_samples_t), intent(inout) :: s               !< Samples.
   type(load_samples_t)      :: w                         !< The samples being built.
   real(R8P),    allocatable :: u(:), v(:), qb(:,:), nb(:,:), ppb(:,:,:), zl(:) !< Vertices, body-frame samples, probes.
   integer(I4P), allocatable :: nxt(:), sb(:), sp(:), sgl(:), vxl(:) !< Next, segment bodies and patches, segments.
   real(R8P)                 :: uv(2)                     !< A vertex.
   integer(I4P)              :: nv, q, n, m, i, l, i0, i1, ni, nj, k, sm !< Work.

   nv = curve%get_vertices_number() ; q = load_rule_size(ord, 2_I4P) ; n = nv * q ; m = size(ll, kind=I4P)
   ni = 0_I4P ; nj = 0_I4P
   if (wl) then
      ni = size(xl, kind=I4P) ; nj = size(yl, kind=I4P)
   endif
   call samples_arrays(w, 2_I4P, ord, q, nv, m, wl, ni, nj, 0_I4P)
   w%lengths = ll
   if (dev) call map_create(w)
   if (n > 0_I4P) then
      allocate(u(nv), v(nv), nxt(nv), sb(nv), sp(nv), zl(1)) ; zl = 0._R8P ; sb = 0_I4P
      do i = 1_I4P, nv
         uv = curve%get_vertex(i) ; u(i) = uv(1) ; v(i) = uv(2) ; nxt(i) = curve%get_next(i) ; sp(i) = curve%segment_patch(i)
      enddo
      do l = 1_I4P, curve%get_loops_number()
         call curve%get_loop_range(l, i0, i1)
         if (i1 >= i0) sb(i0:i1) = curve%loop_body(l)
      enddo
      allocate(qb(2, n), nb(2, n), sgl(n), vxl(n))
      !$acc data copyin(u, v, nxt, sb, sp, xl, yl, zl) create(qb, nb, sgl, vxl) if(dev)
      call points_2d_stage(n, q, ord, nv, u, v, nxt, sb, sp, posed, pr, pt, dev, w%points, w%weight, w%normal, w%element, &
                           w%patch, w%body, qb, nb, sgl, vxl)
      if (m > 0_I4P) then
         allocate(ppb(2, m, n))
         !$acc data create(ppb) if(dev)
         if (dev) then
            call curve%probes_present(n, m, qb, nb, sgl, vxl, w%lengths, ppb, w%visible)
         else
            do sm = 1_I4P, n ! the host visibility: the arithmetic and the kernel of `probes_present`
               do k = 1_I4P, m
                  ppb(1, k, sm) = qb(1, sm) + w%lengths(k) * nb(1, sm) ; ppb(2, k, sm) = qb(2, sm) + w%lengths(k) * nb(2, sm)
                  w%visible(k, sm) = curve%path_visible(qb(1, sm), qb(2, sm), ppb(1, k, sm), ppb(2, k, sm), sgl(sm), &
                                                        vxl(sm))
               enddo
            enddo
         endif
         call probes_2d_stage(n, m, posed, pr, pt, dev, nb, ppb, w%probe, w%visible)
         !$acc end data
         if (wl) call cells_stage(2_I4P, n, m, ni, nj, 1_I4P, xl, yl, zl, dev, w%probe, w%probe_cell)
      endif
      !$acc end data
   endif
   call finish(w, dev, res, s)
   endsubroutine samples_2d_run

   ! private: several bodies (#109 WP3 step 02)

   subroutine bodies_3d_run(surfaces, nfs, ord, ll, posed, prs, pts, tols, wl, xl, yl, zl, dev, res, s, poses)
   !< `load_samples_3d_bodies` after the checks: per body the single-body samples (`samples_3d_run`, resident on the
   !< device path), concatenated in body order (`concat_stage`); with two or more bodies per body b the clearance of
   !< every sample point to b (host `query_many`, device `query_points_present` in b's frame; `merge_clear`, the samples
   !< of b skipped) and the probes' visibility against b (`blocked_3d`). One set of loops for both paths.
   type(surface_stl_object), intent(in), target   :: surfaces(:)          !< Surfaces (checked; mapped on the device path).
   integer(I4P),             intent(in)           :: nfs(:)               !< Facets per body.
   integer(I4P),             intent(in)           :: ord                  !< Rule order.
   real(R8P),                intent(in)           :: ll(:)                !< Probe lengths.
   logical,                  intent(in)           :: posed                !< Poses given.
   real(R8P),                intent(in)           :: prs(:,:,:), pts(:,:) !< Poses [3, 3, nb], [3, nb].
   real(R8P),                intent(in)           :: tols(:)              !< Visibility tolerances [nb].
   logical,                  intent(in)           :: wl                   !< Lines given.
   real(R8P),                intent(in)           :: xl(:), yl(:), zl(:)  !< Lines (one dummy point each without).
   logical,                  intent(in)           :: dev, res             !< Device path, resident.
   type(load_samples_t),     intent(inout)        :: s                    !< Samples.
   type(rigid_pose_t),       intent(in), optional :: poses(:)             !< Poses (host queries).
   type(load_samples_t), allocatable :: sb(:)                             !< Samples per body.
   type(load_samples_t)              :: w                                 !< The samples being built.
   type(vector_R8P),     allocatable :: fp(:)                             !< Sample points (host queries).
   real(R8P),    allocatable :: tl(:), qp(:,:), qd(:), qc(:,:), qn(:,:), qw(:,:) !< Tolerances, queries.
   integer(I4P), allocatable :: qf(:), qft(:), qfi(:), qb(:)              !< Queries.
   integer(I4P)              :: nb, b, q, ne, n, m, off, ni, nj, nk, g    !< Work.

   nb = size(nfs, kind=I4P) ; q = load_rule_size(ord, 3_I4P) ; m = size(ll, kind=I4P)
   ni = 0_I4P ; nj = 0_I4P ; nk = 0_I4P
   if (wl) then
      ni = size(xl, kind=I4P) ; nj = size(yl, kind=I4P) ; nk = size(zl, kind=I4P)
   endif
   allocate(sb(nb), tl(nb)) ; tl = tols(1:nb)
   ne = 0_I4P
   do b = 1_I4P, nb
      if (nfs(b) < 1_I4P) cycle
      call samples_3d_run(surfaces(b), ord, ll, posed, prs(:, :, b), pts(:, b), wl, xl, yl, zl, dev, dev, sb(b))
      ne = ne + nfs(b)
   enddo
   call samples_arrays(w, 3_I4P, ord, q, ne, m, wl, ni, nj, nk) ; call samples_arrays_bodies(w, nb)
   w%lengths = ll ; n = w%n
   if (dev) call map_create(w)
   off = 0_I4P
   do b = 1_I4P, nb
      if (sb(b)%n < 1_I4P) cycle
      call concat_stage(dev, 3_I4P, m, sb(b)%n, size(sb(b)%probe_cell, 3, kind=I4P), n, size(w%probe_cell, 3, kind=I4P), &
                        off, b, sb(b)%points, sb(b)%weight, sb(b)%normal, sb(b)%element, sb(b)%patch, sb(b)%body, &
                        sb(b)%probe, sb(b)%visible, sb(b)%probe_cell, w%points, w%weight, w%normal, w%element, w%patch, &
                        w%body, w%sub_body, w%probe, w%visible, w%probe_cell)
      off = off + sb(b)%n
      call load_samples_reset(sb(b))
   enddo
   if (n > 0_I4P) then
      !$acc data copyin(tl) if(dev)
      call bodies_init_stage(dev, 3_I4P, n, w%clearance, w%second_body, w%second_foot)
      if (nb > 1_I4P) then
         allocate(qp(3, n), qd(n), qc(3, n), qn(3, n), qw(3, n), qf(n), qft(n), qfi(n), qb(n))
         if (.not. dev) then
            allocate(fp(n))
            do g = 1_I4P, n
               fp(g)%x = w%points(1, g) ; fp(g)%y = w%points(2, g) ; fp(g)%z = w%points(3, g)
            enddo
         endif
         do b = 1_I4P, nb
            if (nfs(b) < 1_I4P) cycle
            if (dev) then
               !$acc data create(qp, qd, qc, qn, qw, qf, qft, qfi, qb)
               call body_points_3d_dev(n, w%points, posed, prs(:, :, b), pts(:, b), qp)
               call surfaces(b)%query_points_present(n, qp, .true., .false., qd, qc, qn, qf, qft, qfi, qb)
               call world_points_3d_dev(n, posed, prs(:, :, b), pts(:, b), qc, qw)
               call merge_clear(.true., 3_I4P, n, b, qd, qw, w%body, w%clearance, w%second_body, w%second_foot)
               !$acc end data
            else
               if (posed) then
                  call surfaces(b)%query_many(points=fp, distances=qd, closest_points=qc, pose=poses(b))
               else
                  call surfaces(b)%query_many(points=fp, distances=qd, closest_points=qc)
               endif
               call merge_clear(.false., 3_I4P, n, b, qd, qc, w%body, w%clearance, w%second_body, w%second_foot)
            endif
            if (m > 0_I4P) call blocked_3d(dev, surfaces(b), b, n, m, posed, prs(:, :, b), pts(:, b), nb, tl, w%body, &
                                           w%points, w%probe, w%visible)
         enddo
      endif
      !$acc end data
   endif
   call finish(w, dev, res, s)
   endsubroutine bodies_3d_run

   subroutine bodies_2d_run(curves, nvs, ord, ll, posed, prs, pts, wl, xl, yl, dev, res, s, poses)
   !< 2-D twin of `bodies_3d_run`: `samples_2d_run` per curve, the clearance by `query_many` / `query_present` in the
   !< curve's frame, the visibility by `blocked_2d`.
   type(curve_object),    intent(in), target   :: curves(:)            !< Curves (mapped on the device path).
   integer(I4P),          intent(in)           :: nvs(:)               !< Vertices per curve.
   integer(I4P),          intent(in)           :: ord                  !< Rule order.
   real(R8P),             intent(in)           :: ll(:)                !< Probe lengths.
   logical,               intent(in)           :: posed                !< Poses given.
   real(R8P),             intent(in)           :: prs(:,:,:), pts(:,:) !< Poses [2, 2, nb], [2, nb].
   logical,               intent(in)           :: wl                   !< Lines given.
   real(R8P),             intent(in)           :: xl(:), yl(:)         !< Lines (one dummy point each without).
   logical,               intent(in)           :: dev, res             !< Device path, resident.
   type(load_samples_t),  intent(inout)        :: s                    !< Samples.
   type(rigid_pose_2d_t), intent(in), optional :: poses(:)             !< Poses (host queries).
   type(load_samples_t), allocatable :: sb(:)                          !< Samples per curve.
   type(load_samples_t)              :: w                              !< The samples being built.
   real(R8P),    allocatable :: qp(:,:), qd(:), qc(:,:), qn(:,:), qw(:,:) !< Queries.
   integer(I4P), allocatable :: qs(:), qv(:), qb(:)                    !< Queries.
   integer(I4P)              :: nb, b, q, ne, n, m, off, ni, nj         !< Work.

   nb = size(nvs, kind=I4P) ; q = load_rule_size(ord, 2_I4P) ; m = size(ll, kind=I4P)
   ni = 0_I4P ; nj = 0_I4P
   if (wl) then
      ni = size(xl, kind=I4P) ; nj = size(yl, kind=I4P)
   endif
   allocate(sb(nb))
   ne = 0_I4P
   do b = 1_I4P, nb
      if (nvs(b) < 1_I4P) cycle
      call samples_2d_run(curves(b), ord, ll, posed, prs(:, :, b), pts(:, b), wl, xl, yl, dev, dev, sb(b))
      ne = ne + nvs(b)
   enddo
   call samples_arrays(w, 2_I4P, ord, q, ne, m, wl, ni, nj, 0_I4P) ; call samples_arrays_bodies(w, nb)
   w%lengths = ll ; n = w%n
   if (dev) call map_create(w)
   off = 0_I4P
   do b = 1_I4P, nb
      if (sb(b)%n < 1_I4P) cycle
      call concat_stage(dev, 2_I4P, m, sb(b)%n, size(sb(b)%probe_cell, 3, kind=I4P), n, size(w%probe_cell, 3, kind=I4P), &
                        off, b, sb(b)%points, sb(b)%weight, sb(b)%normal, sb(b)%element, sb(b)%patch, sb(b)%body, &
                        sb(b)%probe, sb(b)%visible, sb(b)%probe_cell, w%points, w%weight, w%normal, w%element, w%patch, &
                        w%body, w%sub_body, w%probe, w%visible, w%probe_cell)
      off = off + sb(b)%n
      call load_samples_reset(sb(b))
   enddo
   if (n > 0_I4P) then
      call bodies_init_stage(dev, 2_I4P, n, w%clearance, w%second_body, w%second_foot)
      if (nb > 1_I4P) then
         allocate(qp(2, n), qd(n), qc(2, n), qn(2, n), qw(2, n), qs(n), qv(n), qb(n))
         do b = 1_I4P, nb
            if (nvs(b) < 1_I4P) cycle
            if (dev) then
               !$acc data create(qp, qd, qc, qn, qw, qs, qv, qb)
               call body_points_2d_dev(n, w%points, posed, prs(:, :, b), pts(:, b), qp)
               call curves(b)%query_present(n, qp, .true., .false., qd, qc, qs, qv, qn, qb)
               call world_points_2d_dev(n, posed, prs(:, :, b), pts(:, b), qc, qw)
               call merge_clear(.true., 2_I4P, n, b, qd, qw, w%body, w%clearance, w%second_body, w%second_foot)
               !$acc end data
            else
               if (posed) then
                  call curves(b)%query_many(points=w%points, distances=qd, closest_points=qc, pose=poses(b))
               else
                  call curves(b)%query_many(points=w%points, distances=qd, closest_points=qc)
               endif
               call merge_clear(.false., 2_I4P, n, b, qd, qc, w%body, w%clearance, w%second_body, w%second_foot)
            endif
            if (m > 0_I4P) call blocked_2d(dev, curves(b), b, n, m, posed, prs(:, :, b), pts(:, b), w%body, w%points, &
                                           w%probe, w%visible)
         enddo
      endif
   endif
   call finish(w, dev, res, s)
   endsubroutine bodies_2d_run

   subroutine concat_stage(dev, d, m, ns, ncs, n, nct, off, b, spts, swt, snrm, sel, spa, sbo, sprb, svis, spc, pts, wt, &
                           nrm, el, pa, bo, sub, prb, vis, pc)
   !< Copy one body's samples (ns, mapped on the device path) into the slots off+1 ... off+ns of the bodies set: every
   !< field as it is, body = b, sub_body = the single-body set's body. Copies only (the same bits on both paths); one
   !< kernel per non-empty group (zero-size arrays are never mapped).
   logical,      intent(in)    :: dev                         !< On the device (arrays present).
   integer(I4P), intent(in)    :: d, m, ns, ncs, n, nct, off, b !< Dims, probes, samples, cells extents, offset, body.
   real(R8P),    intent(in)    :: spts(d, ns), swt(ns), snrm(d, ns) !< Body b: points, weights, normals.
   integer(I4P), intent(in)    :: sel(ns), spa(ns), sbo(ns)   !< Body b: elements, patches, facet_bodies.
   real(R8P),    intent(in)    :: sprb(d, m, ns)              !< Body b: probes.
   logical,      intent(in)    :: svis(m, ns)                 !< Body b: visibility.
   integer(I4P), intent(in)    :: spc(d, m, ncs)              !< Body b: probe cells.
   real(R8P),    intent(inout) :: pts(d, n), wt(n), nrm(d, n) !< Points, weights, normals.
   integer(I4P), intent(inout) :: el(n), pa(n), bo(n), sub(n) !< Elements, patches, bodies, sub-bodies.
   real(R8P),    intent(inout) :: prb(d, m, n)                !< Probes.
   logical,      intent(inout) :: vis(m, n)                   !< Visibility.
   integer(I4P), intent(inout) :: pc(d, m, nct)               !< Probe cells.
   integer(I4P) :: s, k, a, g                                 !< Counters.

   !$acc parallel loop gang vector default(present) present(spts, swt, snrm, sel, spa, sbo, pts, wt, nrm, el, pa, bo, sub) &
   !$acc&  private(a, g) if(dev)
   do s = 1_I4P, ns
      g = off + s
      wt(g) = swt(s) ; el(g) = sel(s) ; pa(g) = spa(s) ; bo(g) = b ; sub(g) = sbo(s)
      do a = 1_I4P, d
         pts(a, g) = spts(a, s) ; nrm(a, g) = snrm(a, s)
      enddo
   enddo
   if (m < 1_I4P) return
   !$acc parallel loop gang vector collapse(2) default(present) present(sprb, svis, prb, vis) private(a, g) if(dev)
   do s = 1_I4P, ns
      do k = 1_I4P, m
         g = off + s
         vis(k, g) = svis(k, s)
         do a = 1_I4P, d
            prb(a, k, g) = sprb(a, k, s)
         enddo
      enddo
   enddo
   if (ncs < 1_I4P .or. nct < 1_I4P) return
   !$acc parallel loop gang vector collapse(2) default(present) present(spc, pc) private(a, g) if(dev)
   do s = 1_I4P, ns
      do k = 1_I4P, m
         g = off + s
         do a = 1_I4P, d
            pc(a, k, g) = spc(a, k, s)
         enddo
      enddo
   enddo
   endsubroutine concat_stage

   subroutine bodies_init_stage(dev, d, n, clr, sb, sf)
   !< Before the clearance merges: clearance huge, no second body, second foot 0 (arrays present if dev).
   logical,      intent(in)    :: dev         !< On the device.
   integer(I4P), intent(in)    :: d, n        !< Dimensions, samples.
   real(R8P),    intent(inout) :: clr(n)      !< Clearances.
   integer(I4P), intent(inout) :: sb(n)       !< Second bodies.
   real(R8P),    intent(inout) :: sf(d, n)    !< Second feet.
   integer(I4P) :: s, a                       !< Counters.

   !$acc parallel loop gang vector default(present) present(clr, sb, sf) private(a) if(dev)
   do s = 1_I4P, n
      clr(s) = huge(1._R8P) ; sb(s) = 0_I4P
      do a = 1_I4P, d
         sf(a, s) = 0._R8P
      enddo
   enddo
   endsubroutine bodies_init_stage

   ! private: the stages (every array present on the device path; on the host the loops run in order)

   subroutine points_3d_stage(surface, n, q, ord, nf, fpatch, posed, pr, pt, dev, pts, wt, nrm, el, pa, bo)
   !< One thread per sample: the rule's point on its facet (`load_point_3d_oac`), its weight and the facet normal, mapped
   !< to the world with a pose; the facet, its patch and body.
   type(surface_stl_object), intent(in), target :: surface         !< Surface.
   integer(I4P),             intent(in)         :: n, q, ord, nf   !< Samples, points per facet, order, facets.
   integer(I4P),             intent(in)         :: fpatch(nf)      !< Facet patches (present).
   logical,                  intent(in)         :: posed, dev      !< A pose, device.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3) !< Pose.
   real(R8P),                intent(inout)      :: pts(3, n)       !< Points (present).
   real(R8P),                intent(inout)      :: wt(n)           !< Weights (present).
   real(R8P),                intent(inout)      :: nrm(3, n)       !< Normals (present).
   integer(I4P),             intent(inout)      :: el(n), pa(n), bo(n) !< Facets, patches, bodies (present).
   real(R8P),    pointer :: tv_p(:,:)                              !< Alias of surface%tri_vertices.
   integer(I4P), pointer :: body_p(:)                              !< Alias of surface%body_payload.
   real(R8P)             :: r(3, 3), t(3), px, py, pz, ww, nx, ny, nz !< Work.
   integer(I4P)          :: s, f, iq                               !< Work.

   r = pr ; t = pt
   tv_p => surface%tri_vertices ; body_p => surface%body_payload
   !$acc parallel loop gang vector default(present) present(tv_p, body_p, fpatch, pts, wt, nrm, el, pa, bo) copyin(r, t) &
   !$acc&  private(f, iq, px, py, pz, ww, nx, ny, nz) if(dev)
   do s = 1_I4P, n
      f = (s - 1_I4P) / q + 1_I4P ; iq = s - (f - 1_I4P) * q
      call load_point_3d_oac(ord, iq, tv_p(1, f), tv_p(2, f), tv_p(3, f), tv_p(4, f), tv_p(5, f), tv_p(6, f), tv_p(7, f), &
                             tv_p(8, f), tv_p(9, f), px, py, pz, ww, nx, ny, nz)
      if (posed) then
         call to_world_3d_oac(r, t, px, py, pz, pts(1, s), pts(2, s), pts(3, s))
         call rotate_3d_oac(r, nx, ny, nz, nrm(1, s), nrm(2, s), nrm(3, s))
      else
         pts(1, s) = px ; pts(2, s) = py ; pts(3, s) = pz
         nrm(1, s) = nx ; nrm(2, s) = ny ; nrm(3, s) = nz
      endif
      wt(s) = ww ; el(s) = f ; pa(s) = fpatch(f) ; bo(s) = body_p(f)
   enddo
   endsubroutine points_3d_stage

   subroutine probes_3d_stage(surface, n, q, ord, m, ll, tol, posed, pr, pt, dev, probe, vis)
   !< One thread per (sample, probe): the body-frame sample again (`load_point_3d_oac`), the probe point + l n, its
   !< visibility (`bvh_segment_blocked_oac` from the point, the kernel of `probe_many`; a degenerate facet: not visible),
   !< the probe in the world.
   type(surface_stl_object), intent(in), target :: surface         !< Surface (SAH BVH).
   integer(I4P),             intent(in)         :: n, q, ord, m    !< Samples, points per facet, order, probes.
   real(R8P),                intent(in)         :: ll(m)           !< Lengths (present).
   real(R8P),                intent(in)         :: tol             !< Visibility tolerance.
   logical,                  intent(in)         :: posed, dev      !< A pose, device.
   real(R8P),                intent(in)         :: pr(3, 3), pt(3) !< Pose.
   real(R8P),                intent(inout)      :: probe(3, m, n)  !< Probes (present).
   logical,                  intent(inout)      :: vis(m, n)       !< Visibility (present).
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias of surface%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias of surface%aabb%payload.
   real(R8P),                    pointer :: tv_p(:,:)    !< Alias of surface%tri_vertices.
   real(R8P)    :: r(3, 3), t(3), px, py, pz, ww, nx, ny, nz, bx, by, bz, lk, wx, wy, wz !< Work.
   integer(I4P) :: s, k, f, iq                           !< Work.

   r = pr ; t = pt
   node_p => surface%aabb%node ; payload_p => surface%aabb%payload ; tv_p => surface%tri_vertices
   !$acc parallel loop gang vector collapse(2) default(present) present(node_p, payload_p, tv_p, ll, probe, vis) &
   !$acc&  copyin(r, t) private(f, iq, px, py, pz, ww, nx, ny, nz, bx, by, bz, lk, wx, wy, wz) if(dev)
   do s = 1_I4P, n
      do k = 1_I4P, m
         f = (s - 1_I4P) / q + 1_I4P ; iq = s - (f - 1_I4P) * q
         call load_point_3d_oac(ord, iq, tv_p(1, f), tv_p(2, f), tv_p(3, f), tv_p(4, f), tv_p(5, f), tv_p(6, f), &
                                tv_p(7, f), tv_p(8, f), tv_p(9, f), px, py, pz, ww, nx, ny, nz)
         lk = ll(k)
         bx = px + lk * nx ; by = py + lk * ny ; bz = pz + lk * nz
         if (nx == 0._R8P .and. ny == 0._R8P .and. nz == 0._R8P) then
            vis(k, s) = .false.
         else
            vis(k, s) = .not. bvh_segment_blocked_oac(node_p, payload_p, tv_p, px, py, pz, bx, by, bz, tol)
         endif
         if (posed) then
            call to_world_3d_oac(r, t, bx, by, bz, wx, wy, wz)
         else
            wx = bx ; wy = by ; wz = bz
         endif
         probe(1, k, s) = wx ; probe(2, k, s) = wy ; probe(3, k, s) = wz
      enddo
   enddo
   endsubroutine probes_3d_stage

   subroutine points_2d_stage(n, q, ord, nv, u, v, nxt, sb, sp, posed, pr, pt, dev, pts, wt, nrm, el, pa, bo, qb, nb, sgl, &
                              vxl)
   !< 2-D twin of `points_3d_stage` (`load_point_2d_oac` on segment i -> nxt(i); none at the end of an open loop: a
   !< degenerate sample at the vertex), with the body-frame points and normals and the segments of the probes.
   integer(I4P), intent(in)    :: n, q, ord, nv   !< Samples, points per segment, order, vertices.
   real(R8P),    intent(in)    :: u(nv), v(nv)    !< Vertices, body frame (present).
   integer(I4P), intent(in)    :: nxt(nv)         !< Next vertex (present).
   integer(I4P), intent(in)    :: sb(nv), sp(nv)  !< Segment bodies, patches (present).
   logical,      intent(in)    :: posed, dev      !< A pose, device.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),    intent(inout) :: pts(2, n)       !< Points (present).
   real(R8P),    intent(inout) :: wt(n)           !< Weights (present).
   real(R8P),    intent(inout) :: nrm(2, n)       !< Normals (present).
   integer(I4P), intent(inout) :: el(n), pa(n), bo(n) !< Segments, patches, bodies (present).
   real(R8P),    intent(inout) :: qb(2, n)        !< Body-frame points (present).
   real(R8P),    intent(inout) :: nb(2, n)        !< Body-frame normals (present).
   integer(I4P), intent(inout) :: sgl(n), vxl(n)  !< Segments and vertices of the probes (present).
   real(R8P)    :: r(2, 2), t(2), pu, pv, ww, nu, nw !< Work.
   integer(I4P) :: s, e, iq, j                    !< Work.

   r = pr ; t = pt
   !$acc parallel loop gang vector default(present) present(u, v, nxt, sb, sp, pts, wt, nrm, el, pa, bo, qb, nb, sgl, vxl) &
   !$acc&  copyin(r, t) private(e, iq, j, pu, pv, ww, nu, nw) if(dev)
   do s = 1_I4P, n
      e = (s - 1_I4P) / q + 1_I4P ; iq = s - (e - 1_I4P) * q ; j = nxt(e)
      if (j < 1_I4P) then
         pu = u(e) ; pv = v(e) ; ww = 0._R8P ; nu = 0._R8P ; nw = 0._R8P
      else
         call load_point_2d_oac(ord, iq, u(e), v(e), u(j), v(j), pu, pv, ww, nu, nw)
      endif
      qb(1, s) = pu ; qb(2, s) = pv ; nb(1, s) = nu ; nb(2, s) = nw ; sgl(s) = e ; vxl(s) = 0_I4P
      if (posed) then
         call to_world_2d_oac(r, t, pu, pv, pts(1, s), pts(2, s))
         call rotate_2d_oac(r, nu, nw, nrm(1, s), nrm(2, s))
      else
         pts(1, s) = pu ; pts(2, s) = pv ; nrm(1, s) = nu ; nrm(2, s) = nw
      endif
      wt(s) = ww ; el(s) = e ; pa(s) = sp(e) ; bo(s) = sb(e)
   enddo
   endsubroutine points_2d_stage

   subroutine probes_2d_stage(n, m, posed, pr, pt, dev, nb, ppb, probe, vis)
   !< One thread per (sample, probe): a degenerate sample's probes not visible, the body-frame probes mapped to the world.
   integer(I4P), intent(in)    :: n, m            !< Samples, probes.
   logical,      intent(in)    :: posed, dev      !< A pose, device.
   real(R8P),    intent(in)    :: pr(2, 2), pt(2) !< Pose.
   real(R8P),    intent(in)    :: nb(2, n)        !< Body-frame normals (present).
   real(R8P),    intent(in)    :: ppb(2, m, n)    !< Body-frame probes (present).
   real(R8P),    intent(inout) :: probe(2, m, n)  !< Probes (present).
   logical,      intent(inout) :: vis(m, n)       !< Visibility (present).
   real(R8P)    :: r(2, 2), t(2)                  !< Pose copies.
   integer(I4P) :: s, k                           !< Counters.

   r = pr ; t = pt
   !$acc parallel loop gang vector collapse(2) default(present) present(nb, ppb, probe, vis) copyin(r, t) if(dev)
   do s = 1_I4P, n
      do k = 1_I4P, m
         if (nb(1, s) == 0._R8P .and. nb(2, s) == 0._R8P) vis(k, s) = .false.
         if (posed) then
            call to_world_2d_oac(r, t, ppb(1, k, s), ppb(2, k, s), probe(1, k, s), probe(2, k, s))
         else
            probe(1, k, s) = ppb(1, k, s) ; probe(2, k, s) = ppb(2, k, s)
         endif
      enddo
   enddo
   endsubroutine probes_2d_stage

   subroutine cells_stage(nd, n, m, ni, nj, nk, x, y, z, dev, probe, pcell)
   !< One thread per (sample, probe): the enclosing cell of the world probe on the lines (`fd_cell_index_oac`).
   integer(I4P), intent(in)    :: nd, n, m        !< Dimensions, samples, probes.
   integer(I4P), intent(in)    :: ni, nj, nk      !< Line points (nk unused in 2-D).
   real(R8P),    intent(in)    :: x(ni), y(nj), z(nk) !< Lines (present).
   logical,      intent(in)    :: dev             !< Device.
   real(R8P),    intent(in)    :: probe(nd, m, n) !< Probes (present).
   integer(I4P), intent(inout) :: pcell(nd, m, n) !< Cells (present).
   integer(I4P) :: s, k                           !< Counters.

   !$acc parallel loop gang vector collapse(2) default(present) present(x, y, z, probe, pcell) if(dev)
   do s = 1_I4P, n
      do k = 1_I4P, m
         pcell(1, k, s) = fd_cell_index_oac(ni, x, probe(1, k, s))
         pcell(2, k, s) = fd_cell_index_oac(nj, y, probe(2, k, s))
         if (nd == 3_I4P) pcell(3, k, s) = fd_cell_index_oac(nk, z, probe(3, k, s))
      enddo
   enddo
   endsubroutine cells_stage

   ! private: the loads

   integer(I4P) function loads_check(samples, d, hp, ht, nforce, moment_ok) result(st)
   !< The checks of `integrate_loads` common to both dimensions: the samples' dimensions, p and/or tau, the force's size,
   !< the moment's.
   type(load_samples_t), intent(in) :: samples   !< Samples.
   integer(I4P),         intent(in) :: d         !< Dimensions of the routine.
   logical,              intent(in) :: hp, ht    !< Pressure, shear given.
   integer(I4P),         intent(in) :: nforce    !< Size of force.
   logical,              intent(in) :: moment_ok !< The moment's size is right.

   st = STATUS_OK
   if (samples%dims /= d .or. .not. (hp .or. ht) .or. nforce /= d .or. .not. moment_ok) st = STATUS_INVALID_INPUT
   if (st == STATUS_OK .and. samples%n > 0_I4P) then
      if (.not. (allocated(samples%points) .and. allocated(samples%weight) .and. allocated(samples%normal))) &
         st = STATUS_INVALID_INPUT
   endif
   endfunction loads_check

   subroutine loads_run(samples, d, nc, pressure, shear, ab, dev, tot, pp, pb)
   !< `integrate_loads` after the checks: the contributions of every sample (one kernel per sample, `if(dev)` on the
   !< device), then their compensated sums in sample order on the host (`loads_reduce`), the same code for both paths.
   type(load_samples_t),   intent(in)           :: samples          !< Samples (checked).
   integer(I4P),           intent(in)           :: d, nc            !< Dimensions, components (3-D 7, 2-D 4).
   real(R8P),              intent(in), optional :: pressure(:)      !< Pressure [n].
   real(R8P),              intent(in), optional :: shear(:,:)       !< Shear [d, n].
   real(R8P),              intent(in)           :: ab(d)            !< About.
   logical,                intent(in)           :: dev              !< Device path.
   real(R8P), allocatable, intent(out)          :: tot(:), pp(:,:), pb(:,:) !< Totals: all, per patch, per body.
   real(R8P),    allocatable :: ctr(:,:), p1(:), t1(:,:) !< Contributions [nc, n]; placeholders of an absent p or tau.
   integer(I4P), allocatable :: idp(:), idb(:)           !< Patch and body of every sample.
   integer(I4P)              :: n                        !< Samples.
   logical                   :: hp, ht                   !< Pressure, shear given.

   n = samples%n ; hp = present(pressure) ; ht = present(shear)
   allocate(ctr(nc, n), idp(n), idb(n), p1(1), t1(d, 1)) ; p1 = 0._R8P ; t1 = 0._R8P
   if (n > 0_I4P) then
      if (hp .and. ht) then
         call loads_stage(d, nc, n, hp, ht, samples%points, samples%weight, samples%normal, samples%patch, samples%body, &
                          n, pressure, n, shear, ab, dev, ctr, idp, idb)
      elseif (hp) then
         call loads_stage(d, nc, n, hp, ht, samples%points, samples%weight, samples%normal, samples%patch, samples%body, &
                          n, pressure, 1_I4P, t1, ab, dev, ctr, idp, idb)
      else
         call loads_stage(d, nc, n, hp, ht, samples%points, samples%weight, samples%normal, samples%patch, samples%body, &
                          1_I4P, p1, n, shear, ab, dev, ctr, idp, idb)
      endif
   endif
   call loads_reduce(nc, n, ctr, idp, idb, tot, pp, pb)
   endsubroutine loads_run

   subroutine loads_stage(d, nc, n, hp, ht, pts, wt, nrm, pa, bo, npr, pres, nsh, shr, ab, dev, ctr, idp, idb)
   !< One thread per sample: its contribution (`load_force_3d_oac` / `load_force_2d_oac`) and its patch and body ids. The
   !< samples, p and tau copied in on the device path unless already mapped (resident samples, a pressure mapped by the
   !< caller: read there); the contributions copied out.
   integer(I4P), intent(in)    :: d, nc, n         !< Dimensions, components, samples.
   logical,      intent(in)    :: hp, ht           !< Pressure, shear given.
   real(R8P),    intent(in)    :: pts(d, n)        !< Points.
   real(R8P),    intent(in)    :: wt(n)            !< Weights.
   real(R8P),    intent(in)    :: nrm(d, n)        !< Normals.
   integer(I4P), intent(in)    :: pa(n), bo(n)     !< Patches, bodies.
   integer(I4P), intent(in)    :: npr, nsh         !< Sizes of pres and shr (1 if absent).
   real(R8P),    intent(in)    :: pres(npr)        !< Pressure.
   real(R8P),    intent(in)    :: shr(d, nsh)      !< Shear.
   real(R8P),    intent(in)    :: ab(d)            !< About.
   logical,      intent(in)    :: dev              !< Device.
   real(R8P),    intent(inout) :: ctr(nc, n)       !< Contributions.
   integer(I4P), intent(inout) :: idp(n), idb(n)   !< Patch and body ids.
   real(R8P)    :: a(3), p, t1, t2, t3             !< About copy, the sample's p and tau.
   integer(I4P) :: s                               !< Counter.

   a = 0._R8P ; a(1:d) = ab
   !$acc data copyin(pts, wt, nrm, pa, bo, pres, shr) copyout(ctr, idp, idb) if(dev)
   !$acc parallel loop gang vector default(present) present(pts, wt, nrm, pa, bo, pres, shr, ctr, idp, idb) copyin(a) &
   !$acc&  private(p, t1, t2, t3) if(dev)
   do s = 1_I4P, n
      p = 0._R8P ; t1 = 0._R8P ; t2 = 0._R8P ; t3 = 0._R8P
      if (hp) p = pres(s)
      if (ht) then
         t1 = shr(1, s) ; t2 = shr(2, s)
         if (d == 3_I4P) t3 = shr(3, s)
      endif
      if (d == 3_I4P) then
         call load_force_3d_oac(hp, ht, p, t1, t2, t3, wt(s), nrm(1, s), nrm(2, s), nrm(3, s), pts(1, s), pts(2, s), &
                                pts(3, s), a(1), a(2), a(3), ctr(1, s), ctr(2, s), ctr(3, s), ctr(4, s), ctr(5, s), &
                                ctr(6, s), ctr(7, s))
      else
         call load_force_2d_oac(hp, ht, p, t1, t2, wt(s), nrm(1, s), nrm(2, s), pts(1, s), pts(2, s), a(1), a(2), &
                                ctr(1, s), ctr(2, s), ctr(3, s), ctr(4, s))
      endif
      idp(s) = pa(s) ; idb(s) = bo(s)
   enddo
   !$acc end data
   endsubroutine loads_stage

   subroutine loads_reduce(nc, n, ctr, idp, idb, tot, pp, pb)
   !< The compensated (Neumaier) sums of the contributions in sample order: all, per patch, per body (ids below 1 enter
   !< the totals only). Host code, the same for both paths.
   integer(I4P),           intent(in)  :: nc, n          !< Components, samples.
   real(R8P),              intent(in)  :: ctr(nc, n)     !< Contributions.
   integer(I4P),           intent(in)  :: idp(n), idb(n) !< Patch and body ids.
   real(R8P), allocatable, intent(out) :: tot(:), pp(:,:), pb(:,:) !< Totals.
   real(R8P), allocatable :: ct(:), cp(:,:), cb(:,:)     !< Compensations.
   integer(I4P)           :: np, nb, s, c, ip, ib        !< Work.

   np = 0_I4P ; nb = 0_I4P
   if (n > 0_I4P) then
      np = max(0_I4P, maxval(idp)) ; nb = max(0_I4P, maxval(idb))
   endif
   allocate(tot(nc), pp(nc, np), pb(nc, nb), ct(nc), cp(nc, np), cb(nc, nb))
   tot = 0._R8P ; pp = 0._R8P ; pb = 0._R8P ; ct = 0._R8P ; cp = 0._R8P ; cb = 0._R8P
   do s = 1_I4P, n
      ip = idp(s) ; ib = idb(s)
      do c = 1_I4P, nc
         call neumaier_add(tot(c), ct(c), ctr(c, s))
         if (ip >= 1_I4P) call neumaier_add(pp(c, ip), cp(c, ip), ctr(c, s))
         if (ib >= 1_I4P) call neumaier_add(pb(c, ib), cb(c, ib), ctr(c, s))
      enddo
   enddo
   tot = tot + ct ; pp = pp + cp ; pb = pb + cb
   endsubroutine loads_reduce

   pure subroutine neumaier_add(sm, cm, x)
   !< Add x to the running sum sm with the compensation cm (Neumaier): the total is sm + cm.
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
   endsubroutine neumaier_add

   ! private: the samples' arrays

   subroutine resolve_inputs(dims, ord, lengths, x, y, z, st, ll, wl, xl, yl, zl)
   !< Check the order and the lines (all or none, >= 2 points, strictly increasing), copy the lengths and the lines (one
   !< dummy point each without lines).
   integer(I4P),           intent(in)           :: dims             !< Dimensions.
   integer(I4P),           intent(in)           :: ord              !< Order.
   real(R8P),              intent(in), optional :: lengths(:)       !< Probe lengths.
   real(R8P),              intent(in), optional :: x(:), y(:), z(:) !< Lines (z ignored in 2-D).
   integer(I4P),           intent(out)          :: st               !< Status.
   real(R8P), allocatable, intent(out)          :: ll(:), xl(:), yl(:), zl(:) !< Copies.
   logical,                intent(out)          :: wl               !< Lines given.
   logical :: lx, ly, lz                                            !< Lines present.

   st = STATUS_OK
   if (load_rule_size(ord, dims) < 1_I4P) st = STATUS_INVALID_INPUT
   if (present(lengths)) then
      allocate(ll(size(lengths))) ; ll = lengths
   else
      allocate(ll(0))
   endif
   lx = present(x) ; ly = present(y) ; lz = present(z) .or. dims == 2_I4P
   wl = present(x) .or. present(y) .or. (dims == 3_I4P .and. present(z))
   allocate(xl(1), yl(1), zl(1)) ; xl = 0._R8P ; yl = 0._R8P ; zl = 0._R8P
   if (.not. wl) return
   if (.not. (lx .and. ly .and. lz)) then
      st = STATUS_INVALID_INPUT
      return
   endif
   if (.not. (line_ok(x) .and. line_ok(y))) st = STATUS_INVALID_INPUT
   deallocate(xl, yl) ; allocate(xl(size(x)), yl(size(y))) ; xl = x ; yl = y
   if (dims == 3_I4P) then
      if (.not. line_ok(z)) st = STATUS_INVALID_INPUT
      deallocate(zl) ; allocate(zl(size(z))) ; zl = z
   endif
   endsubroutine resolve_inputs

   pure logical function line_ok(e)
   !< At least two points, strictly increasing.
   real(R8P), intent(in) :: e(:) !< Line.
   integer :: i                  !< Counter.

   line_ok = size(e) >= 2
   do i = 2, size(e)
      if (.not. (e(i) > e(i - 1))) line_ok = .false.
   enddo
   endfunction line_ok

   subroutine samples_arrays(w, dims, ord, q, ne, m, wl, ni, nj, nk)
   !< Allocate every array of the samples (zero-size the ones not asked) and set their sizes and flags.
   type(load_samples_t), intent(inout) :: w                   !< Samples.
   integer(I4P),         intent(in)    :: dims, ord, q, ne, m !< Dimensions, order, points per element, elements, probes.
   logical,              intent(in)    :: wl                  !< Lines given.
   integer(I4P),         intent(in)    :: ni, nj, nk          !< Line points.
   integer(I4P) :: n, nc                                      !< Samples, cells extent.

   n = ne * q
   w%n = n ; w%dims = dims ; w%order = ord ; w%q = q ; w%ne = ne ; w%m = m
   w%ni = ni ; w%nj = nj ; w%nk = nk ; w%with_cells = wl
   nc = 0_I4P ; if (wl) nc = n
   allocate(w%points(dims, n), w%weight(n), w%normal(dims, n), w%element(n), w%patch(n), w%body(n), w%lengths(m), &
            w%probe(dims, m, n), w%visible(m, n), w%probe_cell(dims, m, nc))
   allocate(w%sub_body(0), w%clearance(0), w%second_body(0), w%second_foot(dims, 0))
   endsubroutine samples_arrays

   subroutine samples_arrays_bodies(w, nb)
   !< The several-bodies arrays of samples allocated by `samples_arrays` (sized n) and their flags.
   type(load_samples_t), intent(inout) :: w  !< Samples.
   integer(I4P),         intent(in)    :: nb !< Bodies.

   w%with_bodies = .true. ; w%bodies = nb
   deallocate(w%sub_body, w%clearance, w%second_body, w%second_foot)
   allocate(w%sub_body(w%n), w%clearance(w%n), w%second_body(w%n), w%second_foot(w%dims, w%n))
   endsubroutine samples_arrays_bodies

   subroutine map_create(w)
   !< Create the device copies of the arrays of samples built on the device (unstructured data on a local set: its arrays
   !< are later moved into the caller's, the present table keyed on the unchanged host addresses).
   type(load_samples_t), intent(inout) :: w !< Samples.
   logical :: ga, gb, gc, gd                !< Non-empty groups (zero-size arrays are never mapped).

   call groups(w, ga, gb, gc, gd)
   !$acc enter data create(w%points, w%weight, w%normal, w%element, w%patch, w%body) if(ga)
   !$acc enter data create(w%lengths, w%probe, w%visible) if(gb)
   !$acc update device(w%lengths) if(gb)
   !$acc enter data create(w%probe_cell) if(gc)
   !$acc enter data create(w%sub_body, w%clearance, w%second_body, w%second_foot) if(gd)
   endsubroutine map_create

   subroutine finish(w, dev, res, s)
   !< End of a run: on the device path resident, the arrays move into the caller's samples, flagged; else they are copied
   !< back and released first. The host path moves them.
   type(load_samples_t), intent(inout) :: w        !< Samples built.
   logical,              intent(in)    :: dev, res !< Device path, resident.
   type(load_samples_t), intent(inout) :: s        !< Caller's samples.
   logical :: ga, gb, gc, gd                       !< Non-empty groups.

   if (dev .and. (.not. res .or. w%n < 1_I4P)) then
      call groups(w, ga, gb, gc, gd)
      !$acc exit data copyout(w%points, w%weight, w%normal, w%element, w%patch, w%body) if(ga)
      !$acc exit data copyout(w%lengths, w%probe, w%visible) if(gb)
      !$acc exit data copyout(w%probe_cell) if(gc)
      !$acc exit data copyout(w%sub_body, w%clearance, w%second_body, w%second_foot) if(gd)
   endif
   s%n = w%n ; s%dims = w%dims ; s%order = w%order ; s%q = w%q ; s%ne = w%ne ; s%m = w%m
   s%ni = w%ni ; s%nj = w%nj ; s%nk = w%nk ; s%with_cells = w%with_cells
   s%with_bodies = w%with_bodies ; s%bodies = w%bodies
   call move_alloc(w%sub_body, s%sub_body) ; call move_alloc(w%clearance, s%clearance)
   call move_alloc(w%second_body, s%second_body) ; call move_alloc(w%second_foot, s%second_foot)
   call move_alloc(w%points, s%points) ; call move_alloc(w%weight, s%weight) ; call move_alloc(w%normal, s%normal)
   call move_alloc(w%element, s%element) ; call move_alloc(w%patch, s%patch) ; call move_alloc(w%body, s%body)
   call move_alloc(w%lengths, s%lengths) ; call move_alloc(w%probe, s%probe) ; call move_alloc(w%visible, s%visible)
   call move_alloc(w%probe_cell, s%probe_cell)
   if (dev .and. res .and. s%n > 0_I4P) call mark_resident(s)
   endsubroutine finish

   subroutine groups(c, ga, gb, gc, gd)
   !< The non-empty array groups of samples (the device copies exist for these only): samples, probes, cells, bodies.
   type(load_samples_t), intent(in)  :: c              !< Samples.
   logical,              intent(out) :: ga, gb, gc, gd !< Groups.

   ga = c%n > 0_I4P
   gb = ga .and. c%m > 0_I4P
   gc = gb .and. c%with_cells
   gd = ga .and. c%with_bodies
   endsubroutine groups

   function host_key(c) result(key)
   !< Host address of points (the key of mapped samples; 0 if points is not allocated or empty).
   type(load_samples_t), intent(in), target :: c   !< Samples.
   integer(c_intptr_t)                      :: key !< Address.

   key = 0
   if (allocated(c%points)) then
      if (size(c%points) > 0) key = transfer(c_loc(c%points), key)
   endif
   endfunction host_key

   subroutine mark_resident(c)
   !< Flag the arrays of `c` as mapped, keyed by the host address of points.
   type(load_samples_t), intent(inout) :: c !< Samples.

   c%on_device = .true. ; c%device_key = host_key(c)
   endsubroutine mark_resident

   logical function is_mapped(c)
   !< The flag is set and the arrays are the ones mapped (not a copy of resident samples).
   type(load_samples_t), intent(in) :: c !< Samples.

   is_mapped = .false.
   if (c%on_device) is_mapped = c%device_key == host_key(c)
   endfunction is_mapped

   integer(I4P) function mis_i1(a, b) result(k)
   !< Differing entries of two integer rank-1 arrays (shapes differ: the larger size).
   integer(I4P), intent(in) :: a(:), b(:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a /= b)
   endif
   endfunction mis_i1

   integer(I4P) function mis_i3(a, b) result(k)
   !< Differing entries of two integer rank-3 arrays (shapes differ: the larger size).
   integer(I4P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a /= b)
   endif
   endfunction mis_i3

   integer(I4P) function mis_r1(a, b) result(k)
   !< Differing entries of two real rank-1 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:), b(:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(a) /= bits(b))
   endif
   endfunction mis_r1

   integer(I4P) function mis_r2(a, b) result(k)
   !< Differing entries of two real rank-2 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:,:), b(:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(reshape(a, [size(a)])) /= bits(reshape(b, [size(b)])))
   endif
   endfunction mis_r2

   integer(I4P) function mis_r3(a, b) result(k)
   !< Differing entries of two real rank-3 arrays by bits (shapes differ: the larger size).
   real(R8P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(bits(reshape(a, [size(a)])) /= bits(reshape(b, [size(b)])))
   endif
   endfunction mis_r3

   integer(I4P) function mis_l(a, b) result(k)
   !< Differing entries of two logical rank-2 arrays (shapes differ: the larger size).
   logical, intent(in) :: a(:,:), b(:,:) !< Arrays.

   if (any(shape(a) /= shape(b))) then
      k = max(size(a), size(b))
   else
      k = count(a .neqv. b)
   endif
   endfunction mis_l

   pure function bits(a) result(b)
   !< The bit patterns of reals.
   real(R8P), intent(in) :: a(:)       !< Reals.
   integer(int64)        :: b(size(a)) !< Bits.

   b = transfer(a, 0_int64, size(a))
   endfunction bits

   ! the per-point kernels (`!$acc routine seq`, shared by the two paths; plain reals)

   pure subroutine load_point_3d_oac(order, iq, x1, y1, z1, x2, y2, z2, x3, y3, z3, px, py, pz, w, nx, ny, nz)
   !< Kernel: point iq of the rule of order `order` on the triangle (v1, v2, v3): p = l1 v1 + l2 v2 + l3 v3, its weight
   !< w = rule weight x area, the unit normal (v2 - v1) x (v3 - v1) / |...|; a zero (or NaN) cross product gives w = 0 and
   !< a zero normal. Explicit parentheses: the same operations in the same order on the host and the device.
   !$acc routine seq
   integer(I4P), intent(in)  :: order, iq                            !< Rule order, point.
   real(R8P),    intent(in)  :: x1, y1, z1, x2, y2, z2, x3, y3, z3   !< Vertices.
   real(R8P),    intent(out) :: px, py, pz                           !< Point.
   real(R8P),    intent(out) :: w                                    !< Weight.
   real(R8P),    intent(out) :: nx, ny, nz                           !< Unit normal.
   real(R8P) :: l1, l2, l3, wr, e1x, e1y, e1z, e2x, e2y, e2z, cx, cy, cz, a2 !< Work.

   call rule_3d_oac(order, iq, l1, l2, l3, wr)
   px = (l1 * x1 + l2 * x2) + l3 * x3
   py = (l1 * y1 + l2 * y2) + l3 * y3
   pz = (l1 * z1 + l2 * z2) + l3 * z3
   e1x = x2 - x1 ; e1y = y2 - y1 ; e1z = z2 - z1
   e2x = x3 - x1 ; e2y = y3 - y1 ; e2z = z3 - z1
   cx = e1y * e2z - e1z * e2y
   cy = e1z * e2x - e1x * e2z
   cz = e1x * e2y - e1y * e2x
   a2 = sqrt((cx * cx + cy * cy) + cz * cz)
   if (a2 > 0._R8P) then
      w = wr * (0.5_R8P * a2)
      nx = cx / a2 ; ny = cy / a2 ; nz = cz / a2
   else
      w = 0._R8P ; nx = 0._R8P ; ny = 0._R8P ; nz = 0._R8P
   endif
   endsubroutine load_point_3d_oac

   pure subroutine rule_3d_oac(order, iq, l1, l2, l3, wr)
   !< Kernel: barycentric coordinates and weight (the weights of a rule sum to 1) of point iq of the 3-D rule of order
   !< `order` (see the module header).
   !$acc routine seq
   integer(I4P), intent(in)  :: order, iq      !< Rule order, point.
   real(R8P),    intent(out) :: l1, l2, l3, wr !< Barycentric coordinates, weight.

   select case (order)
   case (1_I4P)
      l1 = T1 ; l2 = T1 ; l3 = T1 ; wr = 1._R8P
   case (2_I4P)
      l1 = T2B ; l2 = T2B ; l3 = T2B ; wr = T1
      select case (iq)
      case (1_I4P) ; l1 = T2A
      case (2_I4P) ; l2 = T2A
      case default ; l3 = T2A
      endselect
   case default
      if (iq <= 3_I4P) then
         l1 = T3A1 ; l2 = T3A1 ; l3 = T3A1 ; wr = T3W1
         select case (iq)
         case (1_I4P) ; l3 = T3B1
         case (2_I4P) ; l2 = T3B1
         case default ; l1 = T3B1
         endselect
      else
         l1 = T3A2 ; l2 = T3A2 ; l3 = T3A2 ; wr = T3W2
         select case (iq)
         case (4_I4P) ; l3 = T3B2
         case (5_I4P) ; l2 = T3B2
         case default ; l1 = T3B2
         endselect
      endif
   endselect
   endsubroutine rule_3d_oac

   pure subroutine load_point_2d_oac(order, iq, u1, v1, u2, v2, pu, pv, w, nu, nv)
   !< Kernel: point iq of the Gauss-Legendre rule of `order` points on the segment (u1, v1) -> (u2, v2): p = a + s (b - a),
   !< its weight w = rule weight x length, the outward unit normal (dv, -du) / |...| (the solid on the left); a zero (or
   !< NaN) length gives w = 0 and a zero normal.
   !$acc routine seq
   integer(I4P), intent(in)  :: order, iq      !< Rule order (points), point.
   real(R8P),    intent(in)  :: u1, v1, u2, v2 !< Segment ends.
   real(R8P),    intent(out) :: pu, pv         !< Point.
   real(R8P),    intent(out) :: w              !< Weight.
   real(R8P),    intent(out) :: nu, nv         !< Unit normal.
   real(R8P) :: sp, wr, du, dv, ln             !< Work.

   select case (order)
   case (1_I4P)
      sp = 0.5_R8P ; wr = 1._R8P
   case (2_I4P)
      wr = 0.5_R8P ; sp = G2A ; if (iq == 2_I4P) sp = G2B
   case default
      select case (iq)
      case (1_I4P) ; sp = G3A ; wr = G3W1
      case (2_I4P) ; sp = 0.5_R8P ; wr = G3W2
      case default ; sp = G3B ; wr = G3W1
      endselect
   endselect
   du = u2 - u1 ; dv = v2 - v1
   pu = u1 + sp * du ; pv = v1 + sp * dv
   ln = sqrt(du * du + dv * dv)
   if (ln > 0._R8P) then
      w = wr * ln
      nu = dv / ln ; nv = -du / ln
   else
      w = 0._R8P ; nu = 0._R8P ; nv = 0._R8P
   endif
   endsubroutine load_point_2d_oac

   pure subroutine load_force_3d_oac(hp, ht, p, tx, ty, tz, w, nx, ny, nz, px, py, pz, ax, ay, az, fx, fy, fz, mx, my, &
                                     mz, wa)
   !< Kernel: the contribution of one sample to the loads: the traction g = -p n + tau (only the given terms: hp, ht),
   !< f = w g, m = (x - a) x f, and the weight wa = w (the area). The same operations in the same order on the host and
   !< the device.
   !$acc routine seq
   logical,   intent(in)  :: hp, ht                 !< Pressure, shear given.
   real(R8P), intent(in)  :: p, tx, ty, tz          !< Pressure, shear.
   real(R8P), intent(in)  :: w, nx, ny, nz          !< Weight, unit normal.
   real(R8P), intent(in)  :: px, py, pz, ax, ay, az !< Point, point of the moments.
   real(R8P), intent(out) :: fx, fy, fz, mx, my, mz !< Force and moment contributions.
   real(R8P), intent(out) :: wa                     !< Area contribution.
   real(R8P) :: gx, gy, gz, rx, ry, rz              !< Traction, arm.

   if (hp .and. ht) then
      gx = tx - p * nx ; gy = ty - p * ny ; gz = tz - p * nz
   elseif (hp) then
      gx = -(p * nx) ; gy = -(p * ny) ; gz = -(p * nz)
   else
      gx = tx ; gy = ty ; gz = tz
   endif
   fx = w * gx ; fy = w * gy ; fz = w * gz
   rx = px - ax ; ry = py - ay ; rz = pz - az
   mx = ry * fz - rz * fy
   my = rz * fx - rx * fz
   mz = rx * fy - ry * fx
   wa = w
   endsubroutine load_force_3d_oac

   pure subroutine load_force_2d_oac(hp, ht, p, tu, tv, w, nu, nv, pu, pv, au, av, fu, fv, m, wa)
   !< Kernel: 2-D twin of `load_force_3d_oac`: f = w (-p n + tau), the scalar moment m = r_u f_v - r_v f_u, r = x - a.
   !$acc routine seq
   logical,   intent(in)  :: hp, ht         !< Pressure, shear given.
   real(R8P), intent(in)  :: p, tu, tv      !< Pressure, shear.
   real(R8P), intent(in)  :: w, nu, nv      !< Weight, unit normal.
   real(R8P), intent(in)  :: pu, pv, au, av !< Point, point of the moments.
   real(R8P), intent(out) :: fu, fv, m      !< Force and moment contributions.
   real(R8P), intent(out) :: wa             !< Length contribution.
   real(R8P) :: gu, gv, ru, rv              !< Traction, arm.

   if (hp .and. ht) then
      gu = tu - p * nu ; gv = tv - p * nv
   elseif (hp) then
      gu = -(p * nu) ; gv = -(p * nv)
   else
      gu = tu ; gv = tv
   endif
   fu = w * gu ; fv = w * gv
   ru = pu - au ; rv = pv - av
   m = ru * fv - rv * fu
   wa = w
   endsubroutine load_force_2d_oac
endmodule fossil_surface_loads
