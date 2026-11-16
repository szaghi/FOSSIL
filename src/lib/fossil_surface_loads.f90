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
use fossil_aabb_node_object,   only : aabb_node_object
use fossil_aabb_tree_object,   only : facet_distance_payload, AABB_TREE_SAH_BVH
use fossil_curve_object,       only : curve_object
use fossil_fd_closure,         only : fd_cell_index_oac, to_world_3d_oac, rotate_3d_oac, to_world_2d_oac, rotate_2d_oac
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
! the per-point kernels, shared by the two paths (public for the tests)
public :: load_point_3d_oac, load_point_2d_oac

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
   real(R8P),    allocatable :: points(:,:)           !< Sample points [dims, n].
   real(R8P),    allocatable :: weight(:)             !< Weights [n] (rule weight times area / length).
   real(R8P),    allocatable :: normal(:,:)           !< Element unit normal [dims, n] (outward; 0 if degenerate).
   integer(I4P), allocatable :: element(:)            !< Facet / segment [n].
   integer(I4P), allocatable :: patch(:)              !< Patch of the element [n].
   integer(I4P), allocatable :: body(:)               !< Body of the element [n].
   real(R8P),    allocatable :: lengths(:)            !< Probe lengths [m] (+ outward, absolute).
   real(R8P),    allocatable :: probe(:,:,:)          !< Probe points [dims, m, n].
   logical,      allocatable :: visible(:,:)          !< The segment point-probe crosses no element [m, n].
   integer(I4P), allocatable :: probe_cell(:,:,:)     !< Enclosing cell per axis [dims, m, n] (0: outside).
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
   endfunction load_samples_mismatches

   ! the resident pattern

   subroutine enter_device(self)
   !< Map host samples to the device (copyin). Member arrays only, never `self` (#20). No-op if already mapped or if there
   !< is no sample; a copy of resident samples is mapped as new samples.
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc                        !< Non-empty groups.

   if (is_mapped(self)) return
   self%on_device = .false. ; self%device_key = 0
   if (self%n < 1_I4P .or. .not. allocated(self%points)) return
   call groups(self, ga, gb, gc)
   !$acc enter data copyin(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
   !$acc enter data copyin(self%lengths, self%probe, self%visible) if(gb)
   !$acc enter data copyin(self%probe_cell) if(gc)
   call mark_resident(self)
   endsubroutine enter_device

   subroutine update_host(self)
   !< Copy the mapped samples back to the host; they stay mapped. No-op if not mapped; on a copy of resident samples a
   !< warning on stderr (its arrays are not on the device: update the original before copying it).
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc                        !< Non-empty groups.

   if (.not. self%on_device) return
   if (.not. is_mapped(self)) then
      write(stderr, '(A)') 'fossil load_samples_t%update_host: a copy of resident samples is not on the device; '// &
                           'call update_host on the original before copying it'
      self%on_device = .false. ; self%device_key = 0
      return
   endif
   call groups(self, ga, gb, gc)
   !$acc update host(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
   !$acc update host(self%lengths, self%probe, self%visible) if(gb)
   !$acc update host(self%probe_cell) if(gc)
   endsubroutine update_host

   subroutine exit_device(self)
   !< Release the device copies; the host arrays are kept as they are (call `update_host` first to keep the device
   !< values). No-op if not mapped; a copy of resident samples only drops the flag.
   class(load_samples_t), intent(inout) :: self !< Samples.
   logical :: ga, gb, gc                        !< Non-empty groups.

   if (.not. self%on_device) return
   if (is_mapped(self)) then
      call groups(self, ga, gb, gc)
      !$acc exit data delete(self%points, self%weight, self%normal, self%element, self%patch, self%body) if(ga)
      !$acc exit data delete(self%lengths, self%probe, self%visible) if(gb)
      !$acc exit data delete(self%probe_cell) if(gc)
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
   endsubroutine samples_arrays

   subroutine map_create(w)
   !< Create the device copies of the arrays of samples built on the device (unstructured data on a local set: its arrays
   !< are later moved into the caller's, the present table keyed on the unchanged host addresses).
   type(load_samples_t), intent(inout) :: w !< Samples.
   logical :: ga, gb, gc                    !< Non-empty groups (zero-size arrays are never mapped).

   call groups(w, ga, gb, gc)
   !$acc enter data create(w%points, w%weight, w%normal, w%element, w%patch, w%body) if(ga)
   !$acc enter data create(w%lengths, w%probe, w%visible) if(gb)
   !$acc update device(w%lengths) if(gb)
   !$acc enter data create(w%probe_cell) if(gc)
   endsubroutine map_create

   subroutine finish(w, dev, res, s)
   !< End of a run: on the device path resident, the arrays move into the caller's samples, flagged; else they are copied
   !< back and released first. The host path moves them.
   type(load_samples_t), intent(inout) :: w        !< Samples built.
   logical,              intent(in)    :: dev, res !< Device path, resident.
   type(load_samples_t), intent(inout) :: s        !< Caller's samples.
   logical :: ga, gb, gc                           !< Non-empty groups.

   if (dev .and. (.not. res .or. w%n < 1_I4P)) then
      call groups(w, ga, gb, gc)
      !$acc exit data copyout(w%points, w%weight, w%normal, w%element, w%patch, w%body) if(ga)
      !$acc exit data copyout(w%lengths, w%probe, w%visible) if(gb)
      !$acc exit data copyout(w%probe_cell) if(gc)
   endif
   s%n = w%n ; s%dims = w%dims ; s%order = w%order ; s%q = w%q ; s%ne = w%ne ; s%m = w%m
   s%ni = w%ni ; s%nj = w%nj ; s%nk = w%nk ; s%with_cells = w%with_cells
   call move_alloc(w%points, s%points) ; call move_alloc(w%weight, s%weight) ; call move_alloc(w%normal, s%normal)
   call move_alloc(w%element, s%element) ; call move_alloc(w%patch, s%patch) ; call move_alloc(w%body, s%body)
   call move_alloc(w%lengths, s%lengths) ; call move_alloc(w%probe, s%probe) ; call move_alloc(w%visible, s%visible)
   call move_alloc(w%probe_cell, s%probe_cell)
   if (dev .and. res .and. s%n > 0_I4P) call mark_resident(s)
   endsubroutine finish

   subroutine groups(c, ga, gb, gc)
   !< The non-empty array groups of samples (the device copies exist for these only): samples, probes, cells.
   type(load_samples_t), intent(in)  :: c          !< Samples.
   logical,              intent(out) :: ga, gb, gc !< Groups.

   ga = c%n > 0_I4P
   gb = ga .and. c%m > 0_I4P
   gc = gb .and. c%with_cells
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
endmodule fossil_surface_loads
