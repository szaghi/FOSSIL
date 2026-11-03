!< FOSSIL, closed planar curves: the 2-D geometry for quadtree AMR (issue #36).
module fossil_curve_object
!< FOSSIL, closed planar curves: the 2-D geometry for quadtree AMR (issue #36).
!<
!< A `curve_object` is a set of closed polylines (loops) in an axis-aligned plane, the 2-D
!< analogue of `surface_stl_object`.
!<
!< Coordinates. The plane is named by `normal_axis` (1, 2 or 3; default 3, the xy plane); the
!< in-plane coordinates (u, v) are the two remaining axes in cyclic order: (y, z) for
!< normal_axis 1, (z, x) for 2, (x, y) for 3, so that u x v points along +normal_axis.
!< Vertices are stored as structure-of-arrays `u(:)`, `v(:)` (device-mappable, like the 3-D
!< distance payload); a third coordinate would waste memory and carry no information.
!<
!< Topology. Loop l owns the vertices `loop_start(l) : loop_start(l+1) - 1`; a closed loop joins
!< its last vertex to the first implicitly. Segment i runs from vertex i to the next vertex of its
!< loop, so segments and vertices share indices; the last vertex of an *open* loop (a polyline,
!< e.g. from slicing a mesh with holes) starts no segment (`next = 0`). Open loops are invalid
!< (`check_validity`), have no area, and are skipped by the exact (integer) winding number, but
!< the distance and the generalized winding number (`winding_number`) handle them, which is what
!< `SIGN_FAST_WINDING_NUMBER` is for.
!<
!< Orientation. The solid is on the left of every segment: outer loops run counter-clockwise
!< (positive signed area), holes clockwise, islands inside holes counter-clockwise again. The
!< outward normal of segment i is `n = (dv, -du) / length`. The signed distance (step 3) is
!< negative inside, as in 3-D.
!<
!< Exactness. Metrics are plain floating point. Nesting (and so orientation validity) uses the
!< exact `orient2d_xy` crossing test, so every loop gets the same answer from every query.
use fossil_box_query,          only : curve_box_contact_oac, cell_bands, cell_components, BOX_NONE, BOX_CUT, BOX_INSIDE, BOX_OUTSIDE
use fossil_predicates,         only : orient2d_xy
use fossil_segment_bvh,        only : segment_bvh, FEATURE_VERTEX_START, FEATURE_VERTEX_END, segment_point_distance2, &
                                     segment_closest_oac, segment_winding_number_oac
use fossil_surface_stl_object, only : STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL, &
                                     SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use penf,                      only : I4P, R8P
use, intrinsic :: iso_fortran_env, only : stderr => error_unit

implicit none
private
public :: curve_object
public :: CURVE_VALID, CURVE_TOO_FEW_VERTICES, CURVE_ZERO_LENGTH_SEGMENT, CURVE_ZERO_AREA_LOOP, CURVE_WRONG_ORIENTATION
public :: CURVE_OPEN_LOOP

! Validity codes returned by `check_validity` (the first defect found, loop by loop).
integer(I4P), parameter :: CURVE_VALID               = 0_I4P !< Valid curve.
integer(I4P), parameter :: CURVE_TOO_FEW_VERTICES    = 1_I4P !< A loop has fewer than 3 vertices.
integer(I4P), parameter :: CURVE_ZERO_LENGTH_SEGMENT = 2_I4P !< Two consecutive vertices coincide.
integer(I4P), parameter :: CURVE_ZERO_AREA_LOOP      = 3_I4P !< A loop has zero signed area.
integer(I4P), parameter :: CURVE_WRONG_ORIENTATION   = 4_I4P !< A loop's orientation does not match its nesting depth.
integer(I4P), parameter :: CURVE_OPEN_LOOP           = 5_I4P !< A loop is open (a polyline).

type :: curve_object
   !< Closed planar curve: loops of vertices in an axis-aligned plane.
   private
   integer(I4P)              :: normal_axis = 3_I4P !< Plane normal axis (1, 2, 3).
   integer(I4P)              :: nv = 0_I4P          !< Number of vertices (= number of segments).
   integer(I4P)              :: nl = 0_I4P          !< Number of loops.
   real(R8P),    allocatable :: u(:)                !< Vertex u coordinates [nv].
   real(R8P),    allocatable :: v(:)                !< Vertex v coordinates [nv].
   integer(I4P), allocatable :: loop_start(:)       !< Loop offsets [nl+1]; loop l = loop_start(l):loop_start(l+1)-1.
   integer(I4P), allocatable :: loop_of(:)          !< Loop of each vertex [nv].
   logical,      allocatable :: loop_closed(:)      !< Loop is closed [nl].
   integer(I4P), allocatable :: next(:)             !< Next vertex in the loop [nv] (segment i = i -> next(i); 0: none).
   integer(I4P), allocatable :: prev(:)             !< Previous vertex in the loop [nv] (0: none).
   real(R8P),    allocatable :: seg_len(:)          !< Segment lengths [nv].
   real(R8P),    allocatable :: seg_nu(:)           !< Segment outward unit normals, u [nv].
   real(R8P),    allocatable :: seg_nv(:)           !< Segment outward unit normals, v [nv].
   real(R8P),    allocatable :: pn_u(:)             !< Vertex pseudo-normals, u [nv].
   real(R8P),    allocatable :: pn_v(:)             !< Vertex pseudo-normals, v [nv].
   real(R8P),    allocatable :: loop_area(:)        !< Signed area of each loop [nl].
   integer(I4P), allocatable :: loop_body_(:)       !< Body of each loop [nl] (issue #37).
   integer(I4P), allocatable :: vertex_body(:)      !< Body of each vertex/segment [nv] (device-mappable, issue #37).
   integer(I4P)              :: bodies_number = 0_I4P !< Number of bodies.
   real(R8P)                 :: bmin(2) = 0._R8P    !< Bounding box minimum (u, v).
   real(R8P)                 :: bmax(2) = 0._R8P    !< Bounding box maximum (u, v).
   real(R8P)                 :: area = 0._R8P       !< Total signed area (solid area for a valid curve).
   real(R8P)                 :: perimeter = 0._R8P  !< Total length.
   real(R8P)                 :: centroid(2) = 0._R8P !< Centroid of the solid region (u, v).
   type(segment_bvh)         :: bvh                 !< SAH BVH over the segments.
   logical                   :: is_on_device = .false. !< Arrays mapped by enter_device.
   contains
      ! public methods
      procedure, pass(self) :: adopt_loops          !< Build the curve from vertex arrays and loop sizes.
      procedure, pass(self) :: box_classify_many    !< Classify boxes: cut, near, inside, outside (issue #38).
      procedure, pass(self) :: box_intersects_curve !< True if a box meets the curve (issue #38).
      procedure, pass(self) :: box_classify_many_device !< box_classify_many on the device (issue #38).
      procedure, pass(self) :: tag_block            !< Signed tag bands on a block of cells (issue #38).
      procedure, pass(self) :: check_validity       !< First validity defect (CURVE_* code).
      procedure, pass(self) :: compute_distance     !< Distance with closest point, segment and vertex.
      procedure, pass(self) :: distance             !< (Signed) distance of a point.
      procedure, pass(self) :: destroy              !< Reset to an empty curve.
      procedure, pass(self) :: distance_many        !< Batched (signed) distances, CPU (OpenMP).
      procedure, pass(self) :: distance_many_device !< Batched (signed) distances on the device (OpenACC).
      procedure, pass(self) :: enter_device         !< Map the query arrays to the device.
      procedure, pass(self) :: query_many           !< Batched rich query (issue #37), CPU (OpenMP).
      procedure, pass(self) :: query_many_device    !< Batched rich query on the device (OpenACC).
      procedure, pass(self) :: probe_many           !< Probe points along the normal, with exact visibility (issue #37).
      procedure, pass(self) :: image_points         !< Mirror points through the closest point, with visibility (issue #37).
      procedure, pass(self) :: exit_device          !< Release the device arrays.
      procedure, pass(self) :: get_area             !< Total signed area.
      procedure, pass(self) :: get_bmax             !< Bounding box maximum.
      procedure, pass(self) :: get_bmin             !< Bounding box minimum.
      procedure, pass(self) :: get_bodies_number    !< Number of bodies (outer loops with their holes).
      procedure, pass(self) :: get_centroid         !< Centroid of the solid region.
      procedure, pass(self) :: loop_body            !< Body of a loop.
      procedure, pass(self) :: get_is_on_device     !< True if enter_device mapped the arrays.
      procedure, pass(self) :: get_loop_area        !< Signed area of a loop.
      procedure, pass(self) :: get_loop_range       !< First and last vertex of a loop.
      procedure, pass(self) :: get_loops_number     !< Number of loops.
      procedure, pass(self) :: get_normal_axis      !< Plane normal axis.
      procedure, pass(self) :: get_perimeter        !< Total length.
      procedure, pass(self) :: get_vertex           !< Coordinates of a vertex.
      procedure, pass(self) :: get_vertices_number  !< Number of vertices (= segments).
      procedure, pass(self) :: is_loop_closed       !< True if a loop is closed.
      procedure, pass(self) :: is_point_inside      !< Inside test.
      procedure, pass(self) :: load_from_file       !< Load from the ASCII loop format.
      procedure, pass(self) :: nesting_depth        !< Nesting depth of a loop (exact).
      procedure, pass(self) :: repair               !< Join/close chains, merge vertices, drop degenerate loops, orient.
      procedure, pass(self) :: save_into_file       !< Save into the ASCII loop format.
      procedure, pass(self) :: self_intersections   !< Pairs of segments that cross, touch or overlap (exact).
      procedure, pass(self) :: winding_number       !< Generalized winding number (fast or brute force).
      procedure, pass(self) :: winding_number_exact !< Exact integer winding number of a point (brute force).
      ! private methods
      procedure, pass(self), private :: analyze     !< Derived data: topology, metrics, normals.
endtype curve_object

contains
   ! public methods
   subroutine adopt_loops(self, u, v, loop_sizes, normal_axis, closed, status)
   !< Build the curve from vertex coordinates and loop sizes.
   !<
   !< `u(:)`, `v(:)` hold the vertices of loop 1, then loop 2, ... ; `loop_sizes(l)` is the number
   !< of vertices of loop l, closed implicitly unless `closed(l)` is false (an open polyline). The
   !< orientation is taken as given:
   !< use `check_validity` to verify it. `status = STATUS_INVALID_INPUT` (curve left empty) if the
   !< sizes do not add up, a size is negative, an axis is not 1-3, or a coordinate is not finite.
   class(curve_object), intent(inout)         :: self          !< Curve.
   real(R8P),           intent(in)            :: u(:)          !< Vertex u coordinates.
   real(R8P),           intent(in)            :: v(:)          !< Vertex v coordinates.
   integer(I4P),        intent(in)            :: loop_sizes(:) !< Vertices per loop.
   integer(I4P),        intent(in),  optional :: normal_axis   !< Plane normal axis (default 3).
   logical,             intent(in),  optional :: closed(:)     !< Loop closed (default all .true.).
   integer(I4P),        intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                               :: axis, l       !< Resolved axis, counter.

   call self%destroy
   if (present(status)) status = STATUS_OK
   axis = 3_I4P ; if (present(normal_axis)) axis = normal_axis
   if (axis < 1_I4P .or. axis > 3_I4P .or. size(u) /= size(v) .or. any(loop_sizes < 0_I4P) .or. &
       sum(loop_sizes) /= size(u) .or. .not. (all(abs(u) <= huge(1._R8P)) .and. all(abs(v) <= huge(1._R8P)))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(closed)) then
      if (size(closed) /= size(loop_sizes)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   self%normal_axis = axis
   self%nv = size(u, kind=I4P)
   self%nl = size(loop_sizes, kind=I4P)
   allocate(self%u, source=u)
   allocate(self%v, source=v)
   allocate(self%loop_start(self%nl + 1), self%loop_closed(self%nl))
   self%loop_closed = .true. ; if (present(closed)) self%loop_closed = closed
   self%loop_start(1) = 1_I4P
   do l = 1, self%nl
      self%loop_start(l + 1) = self%loop_start(l) + loop_sizes(l)
   enddo
   call self%analyze
   endsubroutine adopt_loops

   subroutine box_classify_many(self, lo, hi, codes, radius, sign_algorithm, inside, status)
   !< Classify boxes (u, v) against the curve (issue #38), OpenMP over the boxes.
   !<
   !< `codes(i)` is BOX_CUT if box i meets a segment (exact, `orient2d`), BOX_NEAR if not but a
   !< segment is within `radius` of it (default 0: no BOX_NEAR; floating-point segment-box distance),
   !< else BOX_INSIDE or BOX_OUTSIDE from `is_point_inside` at the box centre with `sign_algorithm`
   !< (default SIGN_PSEUDO_NORMAL). `inside(i)` (optional): the centre sign of every box that is not
   !< cut (false for cut boxes). Segment BVH traversal (`curve_box_contact_oac`, shared with the
   !< device). `lo(2, n)`, `hi(2, n)`; `status = STATUS_INVALID_INPUT` on a size mismatch or a
   !< negative radius. An empty curve: every box outside.
   class(curve_object), intent(in)            :: self           !< Curve.
   real(R8P),           intent(in)            :: lo(:,:)        !< Box minima [2, n].
   real(R8P),           intent(in)            :: hi(:,:)        !< Box maxima [2, n].
   integer(I4P),        intent(out)           :: codes(:)       !< BOX_* [n].
   real(R8P),           intent(in),  optional :: radius         !< Near radius (default 0).
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_* for the centre query (default SIGN_PSEUDO_NORMAL).
   logical,             intent(out), optional :: inside(:)      !< Centre inside, for boxes not cut [n].
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                  :: r              !< Radius.
   integer(I4P)                               :: algo, i, n, c  !< Work.
   logical                                    :: bad, ins       !< Work.

   if (present(status)) status = STATUS_OK
   n = size(codes, kind=I4P)
   r = 0._R8P ; if (present(radius)) r = radius
   bad = size(lo, 1) /= 2 .or. size(lo, 2) /= n .or. size(hi, 1) /= 2 .or. size(hi, 2) /= n .or. r < 0._R8P
   if (present(inside)) bad = bad .or. size(inside) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   if (self%bvh%nodes_number == 0_I4P) then
      codes = BOX_OUTSIDE
      if (present(inside)) inside = .false.
      return
   endif
   !$omp parallel do default(shared) schedule(dynamic, 16) private(i, c, ins)
   do i = 1_I4P, n
      c = curve_box_contact_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, &
                                self%bvh%right, self%bvh%first, self%bvh%last, self%bvh%seg_id, self%next, self%u, &
                                self%v, lo(:, i), hi(:, i), r)
      ins = .false.
      if (c /= BOX_CUT .and. (c == BOX_NONE .or. present(inside))) &
         ins = self%is_point_inside(0.5_R8P * (lo(:, i) + hi(:, i)), sign_algorithm=algo)
      if (c == BOX_NONE) then
         c = BOX_OUTSIDE ; if (ins) c = BOX_INSIDE
      endif
      codes(i) = c
      if (present(inside)) inside(i) = ins
   enddo
   !$omp end parallel do
   endsubroutine box_classify_many

   pure function box_intersects_curve(self, lo, hi) result(cut)
   !< True if the closed box [lo, hi] (u, v) meets the curve (issue #38, exact).
   class(curve_object), intent(in) :: self         !< Curve.
   real(R8P),           intent(in) :: lo(2), hi(2) !< Box.
   logical                         :: cut          !< Result.

   cut = .false.
   if (self%bvh%nodes_number == 0_I4P) return
   cut = curve_box_contact_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, self%bvh%right, &
                               self%bvh%first, self%bvh%last, self%bvh%seg_id, self%next, self%u, self%v, lo, hi, &
                               0._R8P) == BOX_CUT
   endfunction box_intersects_curve

   subroutine box_classify_many_device(self, lo, hi, codes, radius, sign_algorithm, inside, status)
   !< `box_classify_many` on the device (OpenACC, issue #38): same arguments, codes and defaults.
   !<
   !< Requires `enter_device` (else STATUS_INVALID_INPUT and a stderr message). Each box runs
   !< `curve_box_contact_oac` and, when not cut, `curve_signed_distance_oac` at the centre (sign < 0:
   !< inside), the kernels of the CPU path (`is_point_inside` takes the same sign). SIGN_PSEUDO_NORMAL
   !< or SIGN_FAST_WINDING_NUMBER only. #20 rules: top-level pointer aliases of the mapped arrays.
   !< Validated on the host (tests-gnu, directives inert); the GPU validation is pending.
   class(curve_object), intent(inout), target :: self           !< Curve (target: device aliases).
   real(R8P),           intent(in)            :: lo(:,:)        !< Box minima [2, n].
   real(R8P),           intent(in)            :: hi(:,:)        !< Box maxima [2, n].
   integer(I4P),        intent(out)           :: codes(:)       !< BOX_* [n].
   real(R8P),           intent(in),  optional :: radius         !< Near radius (default 0).
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   logical,             intent(out), optional :: inside(:)      !< Centre inside, for boxes not cut [n].
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                  :: r, blo(2), bhi(2) !< Radius, per-box work.
   integer(I4P)                               :: algo, i, n, c  !< Work.
   logical                                    :: bad, fwn, want_in, ins !< Work.
   integer(I4P), allocatable                  :: cl(:)          !< Device-side codes.
   logical,      allocatable                  :: il(:)          !< Device-side inside flags.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Aliases.
   real(R8P),    pointer :: u(:), v(:), pn_u(:), pn_v(:), wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_diag2(:) !< Aliases.
   real(R8P),    pointer :: wn_c(:,:) !< Alias.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:) !< Aliases.

   if (present(status)) status = STATUS_OK
   n = size(codes, kind=I4P)
   r = 0._R8P ; if (present(radius)) r = radius
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   bad = size(lo, 1) /= 2 .or. size(lo, 2) /= n .or. size(hi, 1) /= 2 .or. size(hi, 2) /= n .or. r < 0._R8P .or. &
         (algo /= SIGN_PSEUDO_NORMAL .and. algo /= SIGN_FAST_WINDING_NUMBER)
   if (present(inside)) bad = bad .or. size(inside) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (self%bvh%nodes_number == 0_I4P) then
      codes = BOX_OUTSIDE
      if (present(inside)) inside = .false.
      return
   endif
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_curve_object%box_classify_many_device: curve not on device -- call enter_device first.'
      return
   endif
   if (n == 0_I4P) return
   fwn = algo == SIGN_FAST_WINDING_NUMBER
   want_in = present(inside)
   allocate(cl(n), il(n))
   lo_u => self%bvh%lo_u ; lo_v => self%bvh%lo_v ; hi_u => self%bvh%hi_u ; hi_v => self%bvh%hi_v
   left => self%bvh%left ; right => self%bvh%right ; first => self%bvh%first ; last => self%bvh%last
   a_u => self%bvh%a_u ; a_v => self%bvh%a_v ; d_u => self%bvh%d_u ; d_v => self%bvh%d_v ; len2 => self%bvh%len2
   seg_id => self%bvh%seg_id ; next => self%next ; u => self%u ; v => self%v ; pn_u => self%pn_u ; pn_v => self%pn_v
   wn_cu => self%bvh%wn_cu ; wn_cv => self%bvh%wn_cv ; wn_du => self%bvh%wn_du ; wn_dv => self%bvh%wn_dv
   wn_c => self%bvh%wn_c ; wn_diag2 => self%bvh%wn_diag2
   !$acc parallel loop present(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, seg_id, &
   !$acc&                      next, u, v, pn_u, pn_v, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2) &
   !$acc&                      copyin(lo, hi) copyout(cl, il) private(blo, bhi, c, ins)
   do i = 1_I4P, n
      blo(1) = lo(1, i) ; blo(2) = lo(2, i) ; bhi(1) = hi(1, i) ; bhi(2) = hi(2, i)
      c = curve_box_contact_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, blo, bhi, r)
      ins = .false.
      if (c /= BOX_CUT .and. (c == BOX_NONE .or. want_in)) &
         ins = curve_signed_distance_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, &
                                         seg_id, next, u, v, pn_u, pn_v, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2, &
                                         0.5_R8P * (blo(1) + bhi(1)), 0.5_R8P * (blo(2) + bhi(2)), .true., fwn) < 0._R8P
      if (c == BOX_NONE) then
         c = BOX_OUTSIDE ; if (ins) c = BOX_INSIDE
      endif
      cl(i) = c ; il(i) = ins
   enddo
   codes = cl
   if (present(inside)) inside = il
   endsubroutine box_classify_many_device

   subroutine tag_block(self, x, y, band, tags, sign_algorithm, status)
   !< Signed tag bands on a rectilinear block of cells in the (u, v) plane (issue #38): the 2-D twin
   !< of `surface_stl_object%tag_block`. `x(nx+1)`, `y(ny+1)` are the u and v cell edges; `tags(nx,
   !< ny)`: 0 on a cut cell (exact), else the Chebyshev cell distance to the nearest cut cell of the
   !< block capped at band + 1, negative inside; one `is_point_inside` per 8-connected component of
   !< uncut cells (one for the whole block when it meets no segment). Distances see this block only.
   !< `status = STATUS_INVALID_INPUT` on bad edges, shape or band.
   class(curve_object), intent(in)            :: self           !< Curve.
   real(R8P),           intent(in)            :: x(:), y(:)     !< Cell edges (u, v).
   integer(I4P),        intent(in)            :: band           !< Band width (>= 0).
   integer(I4P),        intent(out)           :: tags(:,:)      !< Tags [nx, ny].
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_* (default SIGN_PSEUDO_NORMAL).
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,      allocatable :: cut(:,:,:)  !< Cut cells.
   integer(I4P), allocatable :: dist(:,:,:), comp(:,:,:), seed(:,:) !< Bands, components, seeds.
   integer(I4P), allocatable :: sgn(:)      !< Sign of each component.
   integer(I4P)              :: nx, ny, i, j, m, ncomp, algo !< Work.
   logical                   :: bad         !< Size check.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   bad = nx < 1 .or. ny < 1 .or. band < 0_I4P
   if (.not. bad) bad = any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny)) .or. size(tags, 1) /= nx .or. size(tags, 2) /= ny
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   if (.not. self%box_intersects_curve([x(1), y(1)], [x(nx+1), y(ny+1)])) then
      tags = band + 1_I4P
      if (self%bvh%nodes_number > 0_I4P) then
         if (self%is_point_inside([0.5_R8P * (x(1) + x(nx+1)), 0.5_R8P * (y(1) + y(ny+1))], sign_algorithm=algo)) tags = -tags
      endif
      return
   endif
   allocate(cut(nx, ny, 1), dist(nx, ny, 1), comp(nx, ny, 1))
   !$omp parallel do default(shared) schedule(dynamic) private(i, j)
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         cut(i, j, 1) = self%box_intersects_curve([x(i), y(j)], [x(i+1), y(j+1)])
      enddo
   enddo
   !$omp end parallel do
   call cell_bands(cut, band, dist)
   call cell_components(cut, comp, ncomp, seed)
   allocate(sgn(0:ncomp))
   sgn(0) = 0_I4P
   !$omp parallel do default(shared) schedule(dynamic) private(m, i, j)
   do m = 1_I4P, ncomp
      i = seed(1, m) ; j = seed(2, m)
      sgn(m) = 1_I4P
      if (self%is_point_inside([0.5_R8P * (x(i) + x(i+1)), 0.5_R8P * (y(j) + y(j+1))], sign_algorithm=algo)) sgn(m) = -1_I4P
   enddo
   !$omp end parallel do
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         tags(i, j) = sgn(comp(i, j, 1)) * dist(i, j, 1)
      enddo
   enddo
   endsubroutine tag_block

   subroutine distance_many(self, points, distances, is_signed, sign_algorithm, status)
   !< Batched distances of the points (u, v) = points(:, i), CPU, OpenMP over the points.
   !<
   !< Each point runs `curve_signed_distance_oac`, the same kernel as `distance_many_device`, and
   !< matches `compute_distance` bit for bit. Signed by default (negative inside), sign by the
   !< pseudo-normal (default) or the fast winding number. `status = STATUS_INVALID_INPUT` if the
   !< sizes do not match or the sign algorithm is neither of the two.
   class(curve_object), intent(in)            :: self           !< Curve.
   real(R8P),           intent(in)            :: points(:,:)    !< Points [2, n] (u, v).
   real(R8P),           intent(out)           :: distances(:)   !< Distances [n].
   logical,             intent(in),  optional :: is_signed      !< Signed (default .true.).
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   logical                                    :: signed_, fwn   !< Resolved options.
   integer(I4P)                               :: i, n           !< Work.

   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   n = size(distances, kind=I4P)
   if (self%bvh%nodes_number == 0_I4P) then
      distances = huge(1._R8P)
      return
   endif
   !$omp parallel do default(shared) private(i) schedule(static)
   do i = 1, n
      distances(i) = curve_signed_distance_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, &
                     self%bvh%left, self%bvh%right, self%bvh%first, self%bvh%last, self%bvh%a_u, self%bvh%a_v, &
                     self%bvh%d_u, self%bvh%d_v, self%bvh%len2, self%bvh%seg_id, self%next, self%u, self%v, &
                     self%pn_u, self%pn_v, self%bvh%wn_cu, self%bvh%wn_cv, self%bvh%wn_du, self%bvh%wn_dv, &
                     self%bvh%wn_c, self%bvh%wn_diag2, points(1, i), points(2, i), signed_, fwn)
   enddo
   !$omp end parallel do
   endsubroutine distance_many

   subroutine distance_many_device(self, points, distances, is_signed, sign_algorithm, status)
   !< Batched distances on the device (OpenACC): same contract and kernel as `distance_many`.
   !<
   !< Requires `enter_device` (else STATUS_INVALID_INPUT and a stderr message). The compute region
   !< references top-level pointer aliases of the mapped member arrays, never `self%...`
   !< components (#20: an aggregate without a data clause is copied with host addresses).
   !< Validation on the GPU is pending (issue #36 step 7: the GPU was busy); the CPU path runs the
   !< same kernel.
   class(curve_object), intent(inout), target :: self           !< Curve (target: device aliases).
   real(R8P),           intent(in)            :: points(:,:)    !< Points [2, n] (u, v), host.
   real(R8P),           intent(out)           :: distances(:)   !< Distances [n], host.
   logical,             intent(in),  optional :: is_signed      !< Signed (default .true.).
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   logical                                    :: signed_, fwn   !< Resolved options.
   integer(I4P)                               :: i, n           !< Work.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Aliases.
   real(R8P),    pointer :: u(:), v(:), pn_u(:), pn_v(:), wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_diag2(:) !< Aliases.
   real(R8P),    pointer :: wn_c(:,:) !< Alias.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:) !< Aliases.

   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   n = size(distances, kind=I4P)
   if (self%bvh%nodes_number == 0_I4P) then
      distances = huge(1._R8P)
      return
   endif
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_curve_object%distance_many_device: curve not on device -- call enter_device first.'
      return
   endif
   lo_u => self%bvh%lo_u ; lo_v => self%bvh%lo_v ; hi_u => self%bvh%hi_u ; hi_v => self%bvh%hi_v
   left => self%bvh%left ; right => self%bvh%right ; first => self%bvh%first ; last => self%bvh%last
   a_u => self%bvh%a_u ; a_v => self%bvh%a_v ; d_u => self%bvh%d_u ; d_v => self%bvh%d_v ; len2 => self%bvh%len2
   seg_id => self%bvh%seg_id ; next => self%next ; u => self%u ; v => self%v ; pn_u => self%pn_u ; pn_v => self%pn_v
   wn_cu => self%bvh%wn_cu ; wn_cv => self%bvh%wn_cv ; wn_du => self%bvh%wn_du ; wn_dv => self%bvh%wn_dv
   wn_c => self%bvh%wn_c ; wn_diag2 => self%bvh%wn_diag2
   !$acc parallel loop present(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, seg_id, &
   !$acc&                      next, u, v, pn_u, pn_v, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2) &
   !$acc&                      copyin(points) copyout(distances)
   do i = 1, n
      distances(i) = curve_signed_distance_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, &
                     len2, seg_id, next, u, v, pn_u, pn_v, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2, &
                     points(1, i), points(2, i), signed_, fwn)
   enddo
   endsubroutine distance_many_device

   subroutine query_many(self, points, distances, closest_points, segments, vertices, normals, bodies, is_signed, &
                         sign_algorithm, ambiguous, second_foot_points, rel_tolerance, foot_tolerance, status)
   !< Batched rich query (issue #37), CPU, OpenMP over the points: distance (signed by default),
   !< closest point, closest segment, closest vertex (0 in a segment interior), pseudo-normal of the
   !< closest feature (segment normal or vertex bisector), body. Every output after `distances` is
   !< optional; `closest_points(2, n)`, `normals(2, n)`. Each point runs `curve_query_oac`, the kernel
   !< of `query_many_device`, and matches `compute_distance` bit for bit. Refusals as in
   !< `distance_many`.
   !< Multiple foot points (issue #37 step 7, CPU only): `ambiguous(i)` when another segment has its
   !< closest point farther than `foot_tolerance` (default 1e-9 of the box diagonal) from the foot
   !< point within (1 + `rel_tolerance`) (default 1e-6) of the distance; `second_foot_points(:, i)` is
   !< the nearest such point (the foot point itself otherwise).
   class(curve_object), intent(in)            :: self               !< Curve.
   real(R8P),           intent(in)            :: points(:,:)        !< Points [2, n].
   real(R8P),           intent(out)           :: distances(:)       !< Distances [n].
   real(R8P),           intent(out), optional :: closest_points(:,:) !< Closest points [2, n].
   integer(I4P),        intent(out), optional :: segments(:)        !< Closest segments [n].
   integer(I4P),        intent(out), optional :: vertices(:)        !< Closest vertices [n] (0: segment interior).
   real(R8P),           intent(out), optional :: normals(:,:)       !< Pseudo-normals [2, n].
   integer(I4P),        intent(out), optional :: bodies(:)          !< Bodies [n].
   logical,             intent(in),  optional :: is_signed          !< Signed (default .true.).
   integer(I4P),        intent(in),  optional :: sign_algorithm     !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   logical,             intent(out), optional :: ambiguous(:)       !< Another foot point within tolerance [n].
   real(R8P),           intent(out), optional :: second_foot_points(:,:) !< Second foot points [2, n].
   real(R8P),           intent(in),  optional :: rel_tolerance      !< Relative tolerance (default 1e-6).
   real(R8P),           intent(in),  optional :: foot_tolerance     !< Foot separation (default 1e-9 box diagonal).
   integer(I4P),        intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                  :: cp(2), nn(2), cp2(2), rel_, foot_ !< Per-point outputs, tolerances.
   logical                                    :: amb, want_amb      !< Ambiguity.
   integer(I4P)                               :: sg, vx, b, i, n    !< Per-point outputs, counters.
   logical                                    :: signed_, fwn       !< Resolved options.

   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   if (.not. rich_sizes_ok(size(distances, kind=I4P), closest_points, segments, vertices, normals, bodies) .or. &
       .not. rich_sizes_ok(size(distances, kind=I4P), closest_points=second_foot_points)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(ambiguous)) then
      if (size(ambiguous) /= size(distances)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   n = size(distances, kind=I4P)
   if (self%bvh%nodes_number == 0_I4P) then
      distances = huge(1._R8P)
      return
   endif
   want_amb = present(ambiguous) .or. present(second_foot_points)
   rel_ = 1.e-6_R8P ; if (present(rel_tolerance)) rel_ = rel_tolerance
   foot_ = 1.e-9_R8P * hypot(self%bmax(1) - self%bmin(1), self%bmax(2) - self%bmin(2))
   if (present(foot_tolerance)) foot_ = foot_tolerance
   !$omp parallel do default(shared) private(i, cp, nn, sg, vx, b, cp2, amb) schedule(static)
   do i = 1, n
      call curve_query_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, self%bvh%right, &
                           self%bvh%first, self%bvh%last, self%bvh%a_u, self%bvh%a_v, self%bvh%d_u, self%bvh%d_v, &
                           self%bvh%len2, self%bvh%seg_id, self%next, self%u, self%v, self%pn_u, self%pn_v, &
                           self%seg_nu, self%seg_nv, self%vertex_body, self%bvh%wn_cu, self%bvh%wn_cv, self%bvh%wn_du, &
                           self%bvh%wn_dv, self%bvh%wn_c, self%bvh%wn_diag2, points(1, i), points(2, i), signed_, fwn, &
                           distances(i), cp(1), cp(2), sg, vx, nn(1), nn(2), b)
      if (present(closest_points)) closest_points(:, i) = cp
      if (present(segments)) segments(i) = sg
      if (present(vertices)) vertices(i) = vx
      if (present(normals)) normals(:, i) = nn
      if (present(bodies)) bodies(i) = b
      if (want_amb) then
         call second_foot_2d(self, points(:, i), cp, abs(distances(i)), rel_, foot_, amb, cp2)
         if (present(ambiguous)) ambiguous(i) = amb
         if (present(second_foot_points)) second_foot_points(:, i) = cp2
      endif
   enddo
   !$omp end parallel do
   endsubroutine query_many

   subroutine query_many_device(self, points, distances, closest_points, segments, vertices, normals, bodies, is_signed, &
                                sign_algorithm, status)
   !< `query_many` on the device (issue #37): same contract and kernel. Requires `enter_device`.
   !< Every output is computed into local arrays (no absent optional in device code), the
   !< requested ones copied out on the host; #20 pointer aliases. GPU validation pending.
   class(curve_object), intent(inout), target :: self               !< Curve.
   real(R8P),           intent(in)            :: points(:,:)        !< Points [2, n].
   real(R8P),           intent(out)           :: distances(:)       !< Distances [n].
   real(R8P),           intent(out), optional :: closest_points(:,:) !< Closest points [2, n].
   integer(I4P),        intent(out), optional :: segments(:)        !< Closest segments [n].
   integer(I4P),        intent(out), optional :: vertices(:)        !< Closest vertices [n].
   real(R8P),           intent(out), optional :: normals(:,:)       !< Pseudo-normals [2, n].
   integer(I4P),        intent(out), optional :: bodies(:)          !< Bodies [n].
   logical,             intent(in),  optional :: is_signed          !< Signed (default .true.).
   integer(I4P),        intent(in),  optional :: sign_algorithm     !< Sign algorithm.
   integer(I4P),        intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable :: dd(:), cpl(:,:), nml(:,:) !< Device-side outputs.
   integer(I4P), allocatable :: sgl(:), vxl(:), bdl(:)    !< Device-side outputs.
   logical                   :: signed_, fwn              !< Resolved options.
   integer(I4P)              :: i, n                      !< Work.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Aliases.
   real(R8P),    pointer :: u(:), v(:), pn_u(:), pn_v(:), seg_nu(:), seg_nv(:) !< Aliases.
   real(R8P),    pointer :: wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_diag2(:), wn_c(:,:) !< Aliases.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:), vbody(:) !< Aliases.

   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   if (.not. rich_sizes_ok(size(distances, kind=I4P), closest_points, segments, vertices, normals, bodies)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   n = size(distances, kind=I4P)
   if (self%bvh%nodes_number == 0_I4P) then
      distances = huge(1._R8P)
      return
   endif
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_curve_object%query_many_device: curve not on device -- call enter_device first.'
      return
   endif
   allocate(dd(n), cpl(2, n), nml(2, n), sgl(n), vxl(n), bdl(n))
   lo_u => self%bvh%lo_u ; lo_v => self%bvh%lo_v ; hi_u => self%bvh%hi_u ; hi_v => self%bvh%hi_v
   left => self%bvh%left ; right => self%bvh%right ; first => self%bvh%first ; last => self%bvh%last
   a_u => self%bvh%a_u ; a_v => self%bvh%a_v ; d_u => self%bvh%d_u ; d_v => self%bvh%d_v ; len2 => self%bvh%len2
   seg_id => self%bvh%seg_id ; next => self%next ; u => self%u ; v => self%v ; pn_u => self%pn_u ; pn_v => self%pn_v
   seg_nu => self%seg_nu ; seg_nv => self%seg_nv ; vbody => self%vertex_body
   wn_cu => self%bvh%wn_cu ; wn_cv => self%bvh%wn_cv ; wn_du => self%bvh%wn_du ; wn_dv => self%bvh%wn_dv
   wn_c => self%bvh%wn_c ; wn_diag2 => self%bvh%wn_diag2
   !$acc parallel loop present(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, seg_id, &
   !$acc&                      next, u, v, pn_u, pn_v, seg_nu, seg_nv, vbody, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2) &
   !$acc&                      copyin(points) copyout(dd, cpl, nml, sgl, vxl, bdl)
   do i = 1, n
      call curve_query_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, seg_id, next, u, v, &
                           pn_u, pn_v, seg_nu, seg_nv, vbody, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2, &
                           points(1, i), points(2, i), signed_, fwn, dd(i), cpl(1, i), cpl(2, i), sgl(i), vxl(i), &
                           nml(1, i), nml(2, i), bdl(i))
   enddo
   distances = dd
   if (present(closest_points)) closest_points = cpl
   if (present(segments)) segments = sgl
   if (present(vertices)) vertices = vxl
   if (present(normals)) normals = nml
   if (present(bodies)) bodies = bdl
   endsubroutine query_many_device

   subroutine probe_many(self, points, lengths, probe_points, visible, foot_points, normals, status)
   !< Probe points (issue #37): for each point (u, v), the point at distance lengths(i) along the
   !< pseudo-normal of its closest feature from its foot point (+ outward, - inward), and whether
   !< the segment from the foot point to the probe point crosses no segment of the curve other than
   !< those containing the foot point (the closest segment, or the two at the closest vertex). The
   !< crossing test is exact (`orient2d_xy`, candidates from the BVH). OpenMP over the points.
   class(curve_object), intent(in)            :: self              !< Curve.
   real(R8P),           intent(in)            :: points(:,:)       !< Points [2, n].
   real(R8P),           intent(in)            :: lengths(:)        !< Probe lengths [n].
   real(R8P),           intent(out)           :: probe_points(:,:) !< Probe points [2, n].
   logical,             intent(out)           :: visible(:)        !< Segment foot-probe crosses nothing [n].
   real(R8P),           intent(out), optional :: foot_points(:,:)  !< Foot points [2, n].
   real(R8P),           intent(out), optional :: normals(:,:)      !< Pseudo-normals [2, n].
   integer(I4P),        intent(out), optional :: status            !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                  :: d, cp(2), nn(2), tip(2) !< Per-point work.
   integer(I4P)                               :: i, n, sg, vx      !< Work.

   if (present(status)) status = STATUS_OK
   n = size(points, 2, kind=I4P)
   if (size(points, 1) /= 2 .or. size(lengths) /= n .or. size(probe_points, 1) /= 2 .or. size(probe_points, 2) /= n .or. &
       size(visible) /= n .or. .not. rich_sizes_ok(n, closest_points=foot_points, normals=normals)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   !$omp parallel do default(shared) schedule(static) private(i, d, cp, nn, tip, sg, vx)
   do i = 1, n
      call self%compute_distance(p=points(:, i), distance=d, closest_point=cp, closest_segment=sg, closest_vertex=vx, &
                                 normal=nn)
      tip = cp + lengths(i) * nn
      probe_points(:, i) = tip
      visible(i) = path_clear(self, cp, tip, sg, vx)
      if (present(foot_points)) foot_points(:, i) = cp
      if (present(normals)) normals(:, i) = nn
   enddo
   !$omp end parallel do
   endsubroutine probe_many

   subroutine image_points(self, points, images, visible, foot_points, status)
   !< Image points (issue #37): the mirror of each point through its closest point, 2 foot - x, and
   !< whether the segment foot-image crosses no other segment (exact, see `probe_many`).
   class(curve_object), intent(in)            :: self             !< Curve.
   real(R8P),           intent(in)            :: points(:,:)      !< Points [2, n].
   real(R8P),           intent(out)           :: images(:,:)      !< Image points [2, n].
   logical,             intent(out), optional :: visible(:)       !< Segment foot-image crosses nothing [n].
   real(R8P),           intent(out), optional :: foot_points(:,:) !< Foot points [2, n].
   integer(I4P),        intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                  :: d, cp(2), img(2) !< Per-point work.
   integer(I4P)                               :: i, n, sg, vx     !< Work.
   logical                                    :: bad              !< Size check.

   if (present(status)) status = STATUS_OK
   n = size(points, 2, kind=I4P)
   bad = size(points, 1) /= 2 .or. size(images, 1) /= 2 .or. size(images, 2) /= n .or. &
         .not. rich_sizes_ok(n, closest_points=foot_points)
   if (present(visible)) bad = bad .or. size(visible) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   !$omp parallel do default(shared) schedule(static) private(i, d, cp, img, sg, vx)
   do i = 1, n
      call self%compute_distance(p=points(:, i), distance=d, closest_point=cp, closest_segment=sg, closest_vertex=vx)
      img = 2._R8P * cp - points(:, i)
      images(:, i) = img
      if (present(visible)) visible(i) = path_clear(self, cp, img, sg, vx)
      if (present(foot_points)) foot_points(:, i) = cp
   enddo
   !$omp end parallel do
   endsubroutine image_points

   subroutine enter_device(self)
   !< Map the query arrays (tree, payload, winding cache, topology, pseudo-normals) to the device.
   !< Member arrays only, never `self` (#20). No-op if already mapped or the curve is empty.
   class(curve_object), intent(inout) :: self !< Curve.

   if (self%is_on_device .or. self%bvh%nodes_number == 0_I4P) return
   !$acc enter data copyin(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, self%bvh%right, &
   !$acc&                  self%bvh%first, self%bvh%last, self%bvh%a_u, self%bvh%a_v, self%bvh%d_u, self%bvh%d_v, &
   !$acc&                  self%bvh%len2, self%bvh%seg_id, self%next, self%u, self%v, self%pn_u, self%pn_v, &
   !$acc&                  self%bvh%wn_cu, self%bvh%wn_cv, self%bvh%wn_du, self%bvh%wn_dv, self%bvh%wn_c, self%bvh%wn_diag2, &
   !$acc&                  self%seg_nu, self%seg_nv, self%vertex_body)
   self%is_on_device = .true.
   endsubroutine enter_device

   subroutine exit_device(self)
   !< Release the device arrays mapped by `enter_device`. No-op if not mapped.
   class(curve_object), intent(inout) :: self !< Curve.

   if (.not. self%is_on_device) return
   !$acc exit data delete(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, self%bvh%right, &
   !$acc&                 self%bvh%first, self%bvh%last, self%bvh%a_u, self%bvh%a_v, self%bvh%d_u, self%bvh%d_v, &
   !$acc&                 self%bvh%len2, self%bvh%seg_id, self%next, self%u, self%v, self%pn_u, self%pn_v, &
   !$acc&                 self%bvh%wn_cu, self%bvh%wn_cv, self%bvh%wn_du, self%bvh%wn_dv, self%bvh%wn_c, self%bvh%wn_diag2, &
   !$acc&                 self%seg_nu, self%seg_nv, self%vertex_body)
   self%is_on_device = .false.
   endsubroutine exit_device

   function check_validity(self, loop) result(defect)
   !< First validity defect of the curve (CURVE_VALID if none), and the loop it was found in.
   !<
   !< Checked loop by loop: closed; at least 3 vertices; no zero-length segment (exactly coincident
   !< consecutive vertices); nonzero signed area; orientation matching the nesting depth (depth
   !< even = outer or island: counter-clockwise, positive area; odd = hole: clockwise). The
   !< nesting depth is exact (`nesting_depth`). Self-intersections are checked separately, by
   !< `self_intersections`.
   class(curve_object), intent(in)            :: self   !< Curve.
   integer(I4P),        intent(out), optional :: loop   !< Loop of the defect (0 if valid).
   integer(I4P)                               :: defect !< CURVE_* code.
   integer(I4P)                               :: l, i   !< Counters.

   defect = CURVE_VALID
   if (present(loop)) loop = 0_I4P
   do l = 1, self%nl
      if (.not. self%loop_closed(l)) then
         defect = CURVE_OPEN_LOOP
      elseif (self%loop_start(l + 1) - self%loop_start(l) < 3_I4P) then
         defect = CURVE_TOO_FEW_VERTICES
      else
         do i = self%loop_start(l), self%loop_start(l + 1) - 1
            if (self%u(i) == self%u(self%next(i)) .and. self%v(i) == self%v(self%next(i))) then
               defect = CURVE_ZERO_LENGTH_SEGMENT
               exit
            endif
         enddo
         if (defect == CURVE_VALID) then
            if (self%loop_area(l) == 0._R8P) then
               defect = CURVE_ZERO_AREA_LOOP
            elseif ((mod(self%nesting_depth(l), 2_I4P) == 0_I4P) .neqv. (self%loop_area(l) > 0._R8P)) then
               defect = CURVE_WRONG_ORIENTATION
            endif
         endif
      endif
      if (defect /= CURVE_VALID) then
         if (present(loop)) loop = l
         return
      endif
   enddo
   endfunction check_validity

   elemental subroutine destroy(self)
   !< Reset to an empty curve.
   class(curve_object), intent(inout) :: self  !< Curve.
   type(curve_object)                 :: fresh !< Fresh instance.

   select type(self)
   type is(curve_object)
      self = fresh
   endselect
   endsubroutine destroy

   pure function get_area(self) result(a)
   !< Total signed area (the solid area for a valid curve).
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P)                       :: a    !< Area.

   a = self%area
   endfunction get_area

   pure function get_bmax(self) result(b)
   !< Bounding box maximum (u, v).
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P)                       :: b(2) !< Maximum.

   b = self%bmax
   endfunction get_bmax

   pure function get_bmin(self) result(b)
   !< Bounding box minimum (u, v).
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P)                       :: b(2) !< Minimum.

   b = self%bmin
   endfunction get_bmin

   pure function get_centroid(self) result(c)
   !< Centroid (u, v) of the solid region (0 for a zero-area curve).
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P)                       :: c(2) !< Centroid.

   c = self%centroid
   endfunction get_centroid

   pure function get_bodies_number(self) result(n)
   !< Number of bodies: closed loops at even nesting depth (with the holes they contain) and open
   !< loops (see `analyze`).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P)                    :: n    !< Bodies.

   n = self%bodies_number
   endfunction get_bodies_number

   pure function loop_body(self, l) result(b)
   !< Body of loop l (0 if out of range).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P),        intent(in) :: l    !< Loop index.
   integer(I4P)                    :: b    !< Body.

   b = 0_I4P
   if (l >= 1_I4P .and. l <= self%nl) b = self%loop_body_(l)
   endfunction loop_body

   pure function get_is_on_device(self) result(yes)
   !< True if `enter_device` mapped the query arrays.
   class(curve_object), intent(in) :: self !< Curve.
   logical                         :: yes  !< Result.

   yes = self%is_on_device
   endfunction get_is_on_device

   pure function get_loop_area(self, l) result(a)
   !< Signed area of loop l (0 if out of range).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P),        intent(in) :: l    !< Loop index.
   real(R8P)                       :: a    !< Signed area.

   a = 0._R8P
   if (l >= 1_I4P .and. l <= self%nl) a = self%loop_area(l)
   endfunction get_loop_area

   pure subroutine get_loop_range(self, l, first, last)
   !< First and last vertex of loop l (first > last if out of range or empty).
   class(curve_object), intent(in)  :: self  !< Curve.
   integer(I4P),        intent(in)  :: l     !< Loop index.
   integer(I4P),        intent(out) :: first !< First vertex.
   integer(I4P),        intent(out) :: last  !< Last vertex.

   first = 1_I4P ; last = 0_I4P
   if (l < 1_I4P .or. l > self%nl) return
   first = self%loop_start(l)
   last  = self%loop_start(l + 1) - 1_I4P
   endsubroutine get_loop_range

   pure function get_loops_number(self) result(n)
   !< Number of loops.
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P)                    :: n    !< Loops.

   n = self%nl
   endfunction get_loops_number

   pure function get_normal_axis(self) result(axis)
   !< Plane normal axis (1, 2, 3).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P)                    :: axis !< Axis.

   axis = self%normal_axis
   endfunction get_normal_axis

   pure function get_perimeter(self) result(p)
   !< Total length of the curve.
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P)                       :: p    !< Perimeter.

   p = self%perimeter
   endfunction get_perimeter

   pure function get_vertex(self, i) result(uv)
   !< Coordinates (u, v) of vertex i (0 if out of range).
   class(curve_object), intent(in) :: self  !< Curve.
   integer(I4P),        intent(in) :: i     !< Vertex index.
   real(R8P)                       :: uv(2) !< Coordinates.

   uv = 0._R8P
   if (i >= 1_I4P .and. i <= self%nv) uv = [self%u(i), self%v(i)]
   endfunction get_vertex

   pure function get_vertices_number(self) result(n)
   !< Number of vertices (= number of segments).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P)                    :: n    !< Vertices.

   n = self%nv
   endfunction get_vertices_number

   pure subroutine compute_distance(self, p, distance, closest_point, closest_segment, closest_vertex, is_signed, &
                                    use_brute_force, sign_algorithm, normal, body)
   !< Distance of the point p = (u, v) from the curve, with its closest point, segment and vertex.
   !<
   !< `closest_vertex` is the curve vertex when the closest point is a vertex, else 0. Signed
   !< (default): negative inside the solid. The sign comes from the closest feature: the side of
   !< the segment, or the vertex pseudo-normal (a point on the curve has distance 0, so the side
   !< test needs no exact predicate); with `sign_algorithm = SIGN_FAST_WINDING_NUMBER`, from the
   !< generalized winding number (inside iff w > 0.5), robust on open loops and gaps. `use_brute_force` scans every segment (reference). An empty
   !< curve gives distance = huge and zero indices.
   class(curve_object), intent(in)            :: self             !< Curve.
   real(R8P),           intent(in)            :: p(2)             !< Point (u, v).
   real(R8P),           intent(out)           :: distance         !< Distance.
   real(R8P),           intent(out), optional :: closest_point(2) !< Closest point (u, v).
   integer(I4P),        intent(out), optional :: closest_segment  !< Closest segment.
   integer(I4P),        intent(out), optional :: closest_vertex   !< Closest vertex (0 in a segment interior).
   logical,             intent(in),  optional :: is_signed        !< Signed distance (default .true.).
   logical,             intent(in),  optional :: use_brute_force  !< Brute force (default .false.).
   integer(I4P),        intent(in),  optional :: sign_algorithm   !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   real(R8P),           intent(out), optional :: normal(2)        !< Pseudo-normal of the closest feature (issue #37).
   integer(I4P),        intent(out), optional :: body             !< Body of the closest segment (issue #37).
   real(R8P)                                  :: d2, t, side      !< Work.
   integer(I4P)                               :: k, feature, sg, vx !< Work.
   logical                                    :: brute, signed_   !< Resolved options.

   brute = .false. ; if (present(use_brute_force)) brute = use_brute_force
   signed_ = .true. ; if (present(is_signed)) signed_ = is_signed
   if (brute) then
      call self%bvh%closest_brute(p(1), p(2), d2, k, t, feature)
   else
      call self%bvh%closest(p(1), p(2), d2, k, t, feature)
   endif
   if (k == 0_I4P) then
      distance = huge(1._R8P)
      if (present(closest_point)) closest_point = 0._R8P
      if (present(closest_segment)) closest_segment = 0_I4P
      if (present(closest_vertex)) closest_vertex = 0_I4P
      if (present(normal)) normal = 0._R8P
      if (present(body)) body = 0_I4P
      return
   endif
   distance = sqrt(d2)
   sg = self%bvh%seg_id(k)
   select case (feature)
   case (FEATURE_VERTEX_START)
      vx = sg
   case (FEATURE_VERTEX_END)
      vx = self%next(sg)
   case default
      vx = 0_I4P
   endselect
   if (signed_) then
      if (present(sign_algorithm)) then
         if (sign_algorithm == SIGN_FAST_WINDING_NUMBER) then
            side = 0.5_R8P - self%winding_number(p, use_brute_force=brute)
            if (side < 0._R8P) distance = -distance
            signed_ = .false. ! sign done
         endif
      endif
   endif
   if (signed_) then
      if (vx > 0_I4P) then
         side = (p(1) - self%u(vx)) * self%pn_u(vx) + (p(2) - self%v(vx)) * self%pn_v(vx)
      else ! outward normal (dv, -du): side = -cross(d, p - a)
         side = -(self%bvh%d_u(k) * (p(2) - self%bvh%a_v(k)) - self%bvh%d_v(k) * (p(1) - self%bvh%a_u(k)))
      endif
      if (side < 0._R8P) distance = -distance
   endif
   if (present(closest_point)) then
      if (vx > 0_I4P) then
         closest_point = [self%u(vx), self%v(vx)]
      else
         closest_point = [self%bvh%a_u(k) + t * self%bvh%d_u(k), self%bvh%a_v(k) + t * self%bvh%d_v(k)]
      endif
   endif
   if (present(closest_segment)) closest_segment = sg
   if (present(closest_vertex)) closest_vertex = vx
   if (present(normal)) then
      if (vx > 0_I4P) then
         normal = [self%pn_u(vx), self%pn_v(vx)]
      else
         normal = [self%seg_nu(sg), self%seg_nv(sg)]
      endif
   endif
   if (present(body)) body = self%vertex_body(sg)
   endsubroutine compute_distance

   pure function distance(self, p, is_signed, sign_algorithm) result(d)
   !< (Signed) distance of the point p = (u, v) from the curve; negative inside (default signed).
   class(curve_object), intent(in)           :: self      !< Curve.
   real(R8P),           intent(in)           :: p(2)      !< Point (u, v).
   logical,             intent(in), optional :: is_signed !< Signed distance (default .true.).
   integer(I4P),        intent(in), optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   real(R8P)                                 :: d         !< Distance.

   call self%compute_distance(p=p, distance=d, is_signed=is_signed, sign_algorithm=sign_algorithm)
   endfunction distance

   pure function is_loop_closed(self, l) result(closed)
   !< True if loop l is closed (.false. if open or out of range).
   class(curve_object), intent(in) :: self   !< Curve.
   integer(I4P),        intent(in) :: l      !< Loop index.
   logical                         :: closed !< Result.

   closed = .false.
   if (l >= 1_I4P .and. l <= self%nl) closed = self%loop_closed(l)
   endfunction is_loop_closed

   pure function is_point_inside(self, p, sign_algorithm) result(inside)
   !< True if p = (u, v) is inside the solid: signed distance < 0 with the pseudo-normal sign
   !< (default), or generalized winding number > 0.5 with SIGN_FAST_WINDING_NUMBER.
   class(curve_object), intent(in)           :: self           !< Curve.
   real(R8P),           intent(in)           :: p(2)           !< Point (u, v).
   integer(I4P),        intent(in), optional :: sign_algorithm !< Sign algorithm.
   logical                                   :: inside         !< Result.

   if (present(sign_algorithm)) then
      if (sign_algorithm == SIGN_FAST_WINDING_NUMBER) then
         inside = self%winding_number(p) > 0.5_R8P
         return
      endif
   endif
   inside = self%distance(p) < 0._R8P
   endfunction is_point_inside

   pure function winding_number(self, p, use_brute_force, beta) result(w)
   !< Generalized winding number of p = (u, v): (1/2pi) times the sum of the angles the segments
   !< subtend at p. 1 inside and 0 outside closed loops (0 in a hole, 2 where two outer loops
   !< overlap); fractional near open ends. Fast hierarchical evaluation on the BVH (default) or
   !< the exact sum (`use_brute_force`); `beta` is the admissibility parameter (default 2).
   class(curve_object), intent(in)           :: self            !< Curve.
   real(R8P),           intent(in)           :: p(2)            !< Point (u, v).
   logical,             intent(in), optional :: use_brute_force !< Exact sum (default .false.).
   real(R8P),           intent(in), optional :: beta            !< Admissibility (default 2).
   real(R8P)                                 :: w               !< Winding number.
   logical                                   :: brute           !< Resolved option.

   brute = .false. ; if (present(use_brute_force)) brute = use_brute_force
   if (brute) then
      w = self%bvh%winding_number_brute(p(1), p(2))
   else
      w = self%bvh%winding_number(p(1), p(2), beta)
   endif
   endfunction winding_number

   subroutine load_from_file(self, file_name, normal_axis, status)
   !< Load the curve from the ASCII loop format.
   !<
   !< One vertex `u v` per line; loops are separated by one or more blank lines; text after `#`
   !< is a comment, except the directive `# open` before a loop's first vertex, which marks the loop
   !< open (a polyline). A closed loop whose last vertex repeats its first exactly is closed
   !< explicitly: the repeated vertex is dropped. Gnuplot plots the file as it is (`plot 'file' w l`).
   !< `status`: STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL, STATUS_INVALID_INPUT (unreadable
   !< line or invalid data), else STATUS_OK.
   class(curve_object), intent(inout)         :: self        !< Curve.
   character(*),        intent(in)            :: file_name   !< File name.
   integer(I4P),        intent(in),  optional :: normal_axis !< Plane normal axis (default 3).
   integer(I4P),        intent(out), optional :: status      !< Status.
   real(R8P),    allocatable                  :: u(:), v(:)  !< Vertices read.
   integer(I4P), allocatable                  :: sizes(:)    !< Loop sizes read.
   logical,      allocatable                  :: closed(:)   !< Loop closed flags read.
   logical                                    :: open_next   !< Next loop is open.
   character(1024)                            :: line        !< Line buffer.
   real(R8P)                                  :: a, b        !< Pair read.
   integer(I4P)                               :: unit, ios, n, nl, cur, h, status_ !< Work.
   logical                                    :: exists      !< File exists.

   call self%destroy
   status_ = STATUS_OK
   inquire(file=file_name, exist=exists)
   if (.not. exists) then
      if (present(status)) status = STATUS_FILE_NOT_FOUND
      return
   endif
   open(newunit=unit, file=file_name, status='old', action='read', iostat=ios)
   if (ios /= 0) then
      if (present(status)) status = STATUS_FILE_OPEN_FAIL
      return
   endif
   allocate(u(64), v(64), sizes(8), closed(8))
   n = 0_I4P ; nl = 0_I4P ; cur = 0_I4P ; open_next = .false.
   do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      h = index(line, '#')
      if (h > 0) then
         if (adjustl(line(h+1:)) == 'open' .and. len_trim(line(1:h-1)) == 0) open_next = .true.
         line(h:) = ' '
      endif
      if (len_trim(line) == 0) then
         call close_loop
         cycle
      endif
      read(line, *, iostat=ios) a, b
      if (ios /= 0) then
         status_ = STATUS_INVALID_INPUT
         exit
      endif
      if (n == size(u)) then
         u = [u, u] ; v = [v, v]
      endif
      n = n + 1_I4P ; u(n) = a ; v(n) = b ; cur = cur + 1_I4P
   enddo
   close(unit)
   if (status_ == STATUS_OK) then
      call close_loop
      call self%adopt_loops(u=u(1:n), v=v(1:n), loop_sizes=sizes(1:nl), normal_axis=normal_axis, closed=closed(1:nl), &
                            status=status_)
   endif
   if (present(status)) status = status_
   contains
      subroutine close_loop
      !< End the current loop: drop an explicit closing vertex, record the size.
      if (cur == 0_I4P) return
      if (cur > 1_I4P .and. .not. open_next) then
         if (u(n) == u(n - cur + 1) .and. v(n) == v(n - cur + 1)) then
            n = n - 1_I4P ; cur = cur - 1_I4P
         endif
      endif
      if (nl == size(sizes)) then
         sizes = [sizes, sizes] ; closed = [closed, closed]
      endif
      nl = nl + 1_I4P ; sizes(nl) = cur ; closed(nl) = .not. open_next ; cur = 0_I4P ; open_next = .false.
      endsubroutine close_loop
   endsubroutine load_from_file

   pure function nesting_depth(self, l) result(depth)
   !< Nesting depth of loop l: the number of other loops that contain it (exact).
   !<
   !< Loop l is represented by its first vertex; a loop m contains it if the exact winding number
   !< of that vertex with respect to loop m is nonzero. Exact for loops that do not touch (a valid
   !< curve). Open loops contain nothing. 0 if l is out of range.
   class(curve_object), intent(in) :: self  !< Curve.
   integer(I4P),        intent(in) :: l     !< Loop index.
   integer(I4P)                    :: depth !< Nesting depth.
   integer(I4P)                    :: m, i  !< Counters.

   depth = 0_I4P
   if (l < 1_I4P .or. l > self%nl) return
   i = self%loop_start(l)
   do m = 1, self%nl
      if (m == l .or. .not. self%loop_closed(m)) cycle
      if (loop_winding(self, m, self%u(i), self%v(i)) /= 0_I4P) depth = depth + 1_I4P
   enddo
   endfunction nesting_depth

   subroutine repair(self, tolerance, gaps_closed, vertices_merged, loops_removed, loops_reversed, loops_open)
   !< Repair the curve: join and close open chains, merge vertices, drop degenerate loops, orient.
   !<
   !< 1. Open chains: the end of an open chain is joined to the nearest start or end of another
   !<    open chain within `tolerance` (that chain reversed if matched by its end), or closed onto
   !<    its own start, whichever is nearest (ties: own start, then lower loop index, start before
   !<    end); repeated until nothing changes. The gap becomes a segment of length <= tolerance.
   !< 2. Consecutive vertices of a loop within `tolerance` of the last kept one are merged (the
   !<    later dropped); a closed loop also merges its last vertex into its first.
   !< 3. Closed loops with fewer than 3 vertices or zero signed area are removed; open loops with
   !<    fewer than 2 vertices too.
   !< 4. Closed loops are reversed where their orientation does not match their exact nesting
   !<    depth (outer and islands counter-clockwise, holes clockwise).
   !< `tolerance = 0` merges exact duplicates only. Each count is optional.
   class(curve_object), intent(inout)         :: self            !< Curve.
   real(R8P),           intent(in)            :: tolerance       !< Distance tolerance (>= 0).
   integer(I4P),        intent(out), optional :: gaps_closed     !< Gaps closed (joins and self-closures).
   integer(I4P),        intent(out), optional :: vertices_merged !< Vertices merged.
   integer(I4P),        intent(out), optional :: loops_removed   !< Degenerate loops removed.
   integer(I4P),        intent(out), optional :: loops_reversed  !< Loops reversed.
   integer(I4P),        intent(out), optional :: loops_open      !< Loops still open after repair.
   type :: chain_t
      real(R8P), allocatable :: u(:), v(:) !< Vertices.
      logical                :: closed     !< Closed.
      logical                :: alive      !< Not merged into another chain.
   endtype chain_t
   type(chain_t), allocatable :: ch(:)        !< Loops as chains.
   real(R8P),     allocatable :: nu(:), nv(:) !< Rebuilt vertices.
   integer(I4P),  allocatable :: sizes(:)     !< Rebuilt loop sizes.
   logical,       allocatable :: closed(:)    !< Rebuilt closed flags.
   real(R8P),     allocatable :: tu(:), tv(:) !< Work vertices.
   type(curve_object)         :: tmp          !< Curve after steps 1-3, for the nesting depths.
   integer(I4P)               :: l, m, i, k, best, best_end, gaps, merged, removed, reversed, first, last, n, axis !< Work.
   real(R8P)                  :: d, dbest     !< Distances.
   logical                    :: changed      !< Join loop flag.

   gaps = 0_I4P ; merged = 0_I4P ; removed = 0_I4P ; reversed = 0_I4P
   axis = self%normal_axis
   allocate(ch(self%nl))
   do l = 1, self%nl
      ch(l)%u = self%u(self%loop_start(l):self%loop_start(l + 1) - 1)
      ch(l)%v = self%v(self%loop_start(l):self%loop_start(l + 1) - 1)
      ch(l)%closed = self%loop_closed(l)
      ch(l)%alive = .true.
   enddo
   ! 1. join and close open chains
   changed = .true.
   do while (changed)
      changed = .false.
      do l = 1, size(ch)
         if (.not. ch(l)%alive .or. ch(l)%closed .or. size(ch(l)%u) == 0) cycle
         n = size(ch(l)%u)
         best = 0_I4P ; best_end = 0_I4P ; dbest = huge(1._R8P)
         if (n > 2) then ! closing onto itself
            d = hypot(ch(l)%u(n) - ch(l)%u(1), ch(l)%v(n) - ch(l)%v(1))
            if (d <= tolerance) then
               best = l ; dbest = d
            endif
         endif
         do m = 1, size(ch)
            if (m == l .or. .not. ch(m)%alive .or. ch(m)%closed .or. size(ch(m)%u) == 0) cycle
            d = hypot(ch(l)%u(n) - ch(m)%u(1), ch(l)%v(n) - ch(m)%v(1))
            if (d <= tolerance .and. d < dbest) then
               best = m ; best_end = 0_I4P ; dbest = d
            endif
            k = size(ch(m)%u)
            d = hypot(ch(l)%u(n) - ch(m)%u(k), ch(l)%v(n) - ch(m)%v(k))
            if (d <= tolerance .and. d < dbest) then
               best = m ; best_end = 1_I4P ; dbest = d
            endif
         enddo
         if (best == 0_I4P) cycle
         gaps = gaps + 1_I4P ; changed = .true.
         if (best == l) then
            ch(l)%closed = .true.
         else
            if (best_end == 1_I4P) then
               ch(best)%u = ch(best)%u(size(ch(best)%u):1:-1) ; ch(best)%v = ch(best)%v(size(ch(best)%v):1:-1)
            endif
            ch(l)%u = [ch(l)%u, ch(best)%u] ; ch(l)%v = [ch(l)%v, ch(best)%v]
            ch(best)%alive = .false.
         endif
         exit ! restart the scan: chain indices and ends changed
      enddo
   enddo
   ! 2. merge consecutive vertices
   do l = 1, size(ch)
      if (.not. ch(l)%alive) cycle
      n = size(ch(l)%u)
      if (n == 0) cycle
      allocate(tu(n), tv(n))
      k = 1_I4P ; tu(1) = ch(l)%u(1) ; tv(1) = ch(l)%v(1)
      do i = 2, n
         if (hypot(ch(l)%u(i) - tu(k), ch(l)%v(i) - tv(k)) <= tolerance) then
            merged = merged + 1_I4P
         else
            k = k + 1_I4P ; tu(k) = ch(l)%u(i) ; tv(k) = ch(l)%v(i)
         endif
      enddo
      if (ch(l)%closed .and. k > 1_I4P) then
         if (hypot(tu(k) - tu(1), tv(k) - tv(1)) <= tolerance) then
            merged = merged + 1_I4P ; k = k - 1_I4P
         endif
      endif
      ch(l)%u = tu(1:k) ; ch(l)%v = tv(1:k)
      deallocate(tu, tv)
   enddo
   ! 3. drop degenerate loops
   do l = 1, size(ch)
      if (.not. ch(l)%alive) cycle
      n = size(ch(l)%u)
      if (ch(l)%closed) then
         if (n < 3) then
            ch(l)%alive = .false.
         elseif (shoelace(ch(l)%u, ch(l)%v) == 0._R8P) then
            ch(l)%alive = .false.
         endif
      elseif (n < 2) then
         ch(l)%alive = .false.
      endif
      if (.not. ch(l)%alive) removed = removed + 1_I4P
   enddo
   call pack_chains
   ! 4. orient closed loops by nesting depth
   call tmp%adopt_loops(u=nu, v=nv, loop_sizes=sizes, normal_axis=axis, closed=closed)
   first = 1_I4P
   do l = 1, size(sizes)
      last = first + sizes(l) - 1_I4P
      if (closed(l)) then
         if ((mod(tmp%nesting_depth(l), 2_I4P) == 0_I4P) .neqv. (tmp%loop_area(l) > 0._R8P)) then
            nu(first:last) = nu(last:first:-1) ; nv(first:last) = nv(last:first:-1)
            reversed = reversed + 1_I4P
         endif
      endif
      first = last + 1_I4P
   enddo
   call self%adopt_loops(u=nu, v=nv, loop_sizes=sizes, normal_axis=axis, closed=closed)
   if (present(gaps_closed)) gaps_closed = gaps
   if (present(vertices_merged)) vertices_merged = merged
   if (present(loops_removed)) loops_removed = removed
   if (present(loops_reversed)) loops_reversed = reversed
   if (present(loops_open)) loops_open = count(.not. closed)
   contains
      subroutine pack_chains
      !< Rebuild the vertex arrays, sizes and closed flags from the alive chains.
      integer(I4P) :: q, nt, nl2 !< Work.

      nt = 0_I4P ; nl2 = 0_I4P
      do q = 1, size(ch)
         if (.not. ch(q)%alive) cycle
         nt = nt + size(ch(q)%u, kind=I4P) ; nl2 = nl2 + 1_I4P
      enddo
      allocate(nu(nt), nv(nt), sizes(nl2), closed(nl2))
      nt = 0_I4P ; nl2 = 0_I4P
      do q = 1, size(ch)
         if (.not. ch(q)%alive) cycle
         nl2 = nl2 + 1_I4P
         sizes(nl2) = size(ch(q)%u, kind=I4P) ; closed(nl2) = ch(q)%closed
         nu(nt + 1:nt + sizes(nl2)) = ch(q)%u ; nv(nt + 1:nt + sizes(nl2)) = ch(q)%v
         nt = nt + sizes(nl2)
      enddo
      endsubroutine pack_chains
   endsubroutine repair

   subroutine save_into_file(self, file_name, status)
   !< Save the curve into the ASCII loop format (see `load_from_file`), bit-exact on reload.
   class(curve_object), intent(in)            :: self      !< Curve.
   character(*),        intent(in)            :: file_name !< File name.
   integer(I4P),        intent(out), optional :: status    !< STATUS_OK or STATUS_FILE_OPEN_FAIL.
   integer(I4P)                               :: unit, ios, l, i !< Work.

   if (present(status)) status = STATUS_OK
   open(newunit=unit, file=file_name, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      if (present(status)) status = STATUS_FILE_OPEN_FAIL
      return
   endif
   write(unit, '(A,I0,A,I0)') '# FOSSIL curve: loops ', self%nl, ', normal_axis ', self%normal_axis
   do l = 1, self%nl
      if (l > 1_I4P) write(unit, '(A)') ''
      if (.not. self%loop_closed(l)) write(unit, '(A)') '# open'
      do i = self%loop_start(l), self%loop_start(l + 1) - 1
         write(unit, '(2ES26.17E3)') self%u(i), self%v(i)
      enddo
   enddo
   close(unit)
   endsubroutine save_into_file

   pure subroutine self_intersections(self, pairs, n)
   !< Pairs of segments (i < j) that cross, touch or overlap, decided exactly (`orient2d_xy`).
   !<
   !< Adjacent segments (sharing a vertex) are reported only if they overlap beyond the shared
   !< vertex (a collinear fold-back). Candidates come from the BVH; `pairs(1:2, 1:n)` is
   !< (re)allocated, in increasing i, then traversal order.
   class(curve_object),       intent(in)    :: self       !< Curve.
   integer(I4P), allocatable, intent(inout) :: pairs(:,:) !< Segment pairs [2, n].
   integer(I4P),              intent(out)   :: n          !< Number of pairs.
   integer(I4P), allocatable                :: ids(:), buf(:,:), tmp(:,:) !< Candidates, buffers.
   integer(I4P)                             :: i, j, c, m, q1, q2 !< Work.
   logical                                  :: hit        !< Pair intersects.

   n = 0_I4P
   allocate(buf(2, 16))
   do i = 1, self%nv
      j = self%next(i)
      if (j == 0_I4P) cycle ! last vertex of an open loop
      call self%bvh%segments_in_box(min(self%u(i), self%u(j)), min(self%v(i), self%v(j)), &
                                    max(self%u(i), self%u(j)), max(self%v(i), self%v(j)), ids, m)
      do c = 1, m
         q1 = self%bvh%seg_id(ids(c))
         if (q1 <= i) cycle
         q2 = self%next(q1)
         if (q1 == j .and. q2 == i) then ! a 2-vertex loop: the two segments coincide
            hit = .true.
         elseif (q1 == j) then         ! adjacent, sharing j: fold-back only
            hit = folds_back(self%u(i), self%v(i), self%u(j), self%v(j), self%u(q2), self%v(q2))
         elseif (q2 == i) then         ! adjacent, sharing i
            hit = folds_back(self%u(j), self%v(j), self%u(i), self%v(i), self%u(q1), self%v(q1))
         else
            hit = segments_intersect(self%u(i), self%v(i), self%u(j), self%v(j), &
                                     self%u(q1), self%v(q1), self%u(q2), self%v(q2))
         endif
         if (hit) then
            if (n == size(buf, dim=2)) then
               allocate(tmp(2, 2 * n))
               tmp(:, 1:n) = buf(:, 1:n)
               call move_alloc(from=tmp, to=buf)
            endif
            n = n + 1_I4P ; buf(:, n) = [i, q1]
         endif
      enddo
   enddo
   if (allocated(pairs)) deallocate(pairs)
   allocate(pairs(2, n))
   pairs(:, 1:n) = buf(:, 1:n)
   endsubroutine self_intersections

   pure function winding_number_exact(self, p) result(w)
   !< Exact integer winding number of the point p = (u, v) with respect to the whole curve.
   !<
   !< Brute force over all segments, with the exact `orient2d_xy` crossing rule (half-open in v):
   !< 1 inside an outer loop, 0 outside and inside a hole. Undefined (but deterministic) for a
   !< point exactly on the curve. Closed loops only (open loops are skipped). Reference for tests
   !< and for repair; `winding_number` is the generalized (and fast) one.
   class(curve_object), intent(in) :: self !< Curve.
   real(R8P),           intent(in) :: p(2) !< Point (u, v).
   integer(I4P)                    :: w    !< Winding number.
   integer(I4P)                    :: l    !< Counter.

   w = 0_I4P
   do l = 1, self%nl
      if (self%loop_closed(l)) w = w + loop_winding(self, l, p(1), p(2))
   enddo
   endfunction winding_number_exact

   ! private methods
   pure subroutine analyze(self)
   !< Derived data: loop topology, segment lengths and normals, vertex pseudo-normals, metrics.
   !<
   !< Vertex pseudo-normal: the normalized sum of the incident segment normals (in 2-D the
   !< angle-weighted pseudo-normal of Baerentzen and Aanaes reduces to this bisector, both
   !< incident "angles" being equal; at the end of an open loop, the one segment normal). With
   !< it, the sign of (p - x) . n at the closest feature x is the inside/outside sign. Open loops
   !< have no closing segment and contribute no area.
   class(curve_object), intent(inout) :: self !< Curve.
   integer(I4P), allocatable          :: ids(:) !< Vertices that start a segment.
   integer(I4P)                       :: l, i, j, first, last !< Counters.
   real(R8P)                          :: du, dv, cross, cu, cv, s !< Work.

   allocate(self%loop_of(self%nv), self%next(self%nv), self%prev(self%nv))
   allocate(self%seg_len(self%nv), self%seg_nu(self%nv), self%seg_nv(self%nv))
   allocate(self%pn_u(self%nv), self%pn_v(self%nv), self%loop_area(self%nl))
   do l = 1, self%nl
      first = self%loop_start(l) ; last = self%loop_start(l + 1) - 1_I4P
      do i = first, last
         self%loop_of(i) = l
         self%next(i) = i + 1_I4P
         self%prev(i) = i - 1_I4P
      enddo
      if (last >= first) then
         if (self%loop_closed(l)) then
            self%next(last) = first
            self%prev(first) = last
         else
            self%next(last) = 0_I4P
            self%prev(first) = 0_I4P
         endif
      endif
   enddo
   self%perimeter = 0._R8P
   cu = 0._R8P ; cv = 0._R8P
   self%loop_area = 0._R8P
   self%seg_len = 0._R8P ; self%seg_nu = 0._R8P ; self%seg_nv = 0._R8P
   do i = 1, self%nv
      j = self%next(i)
      if (j == 0_I4P) cycle
      du = self%u(j) - self%u(i)
      dv = self%v(j) - self%v(i)
      self%seg_len(i) = sqrt(du * du + dv * dv)
      if (self%seg_len(i) > 0._R8P) then
         self%seg_nu(i) =  dv / self%seg_len(i)
         self%seg_nv(i) = -du / self%seg_len(i)
      endif
      self%perimeter = self%perimeter + self%seg_len(i)
      if (.not. self%loop_closed(self%loop_of(i))) cycle
      cross = self%u(i) * self%v(j) - self%u(j) * self%v(i) ! shoelace
      self%loop_area(self%loop_of(i)) = self%loop_area(self%loop_of(i)) + 0.5_R8P * cross
      cu = cu + (self%u(i) + self%u(j)) * cross
      cv = cv + (self%v(i) + self%v(j)) * cross
   enddo
   do i = 1, self%nv
      j = self%prev(i)
      self%pn_u(i) = self%seg_nu(i) ; self%pn_v(i) = self%seg_nv(i) ! zero if i starts no segment
      if (j > 0_I4P) then
         self%pn_u(i) = self%pn_u(i) + self%seg_nu(j)
         self%pn_v(i) = self%pn_v(i) + self%seg_nv(j)
      endif
      s = sqrt(self%pn_u(i)**2 + self%pn_v(i)**2)
      if (s > 0._R8P) then
         self%pn_u(i) = self%pn_u(i) / s
         self%pn_v(i) = self%pn_v(i) / s
      endif
   enddo
   self%area = sum(self%loop_area)
   self%centroid = 0._R8P
   if (self%area /= 0._R8P) self%centroid = [cu, cv] / (6._R8P * self%area)
   if (self%nv > 0_I4P) then
      self%bmin = [minval(self%u), minval(self%v)]
      self%bmax = [maxval(self%u), maxval(self%v)]
   else
      self%bmin = 0._R8P ; self%bmax = 0._R8P
   endif
   call compute_bodies(self)
   allocate(self%vertex_body(self%nv))
   do i = 1, self%nv
      self%vertex_body(i) = self%loop_body_(self%loop_of(i))
   enddo
   ids = pack([(i, i = 1, self%nv)], self%next > 0_I4P)
   call self%bvh%build(self%u(ids), self%v(ids), self%u(self%next(ids)), self%v(self%next(ids)), ids=ids)
   endsubroutine analyze

   ! private non type-bound procedures
   pure function loop_winding(self, l, pu, pv) result(w)
   !< Exact integer winding number of (pu, pv) with respect to loop l.
   !<
   !< Crossing rule half-open in v: an upward segment (v_i <= pv < v_j) counts +1 if the point is
   !< strictly left of it, a downward one (v_j <= pv < v_i) counts -1 if strictly right. The side
   !< test is the exact `orient2d_xy`.
   type(curve_object), intent(in) :: self   !< Curve.
   integer(I4P),       intent(in) :: l      !< Loop.
   real(R8P),          intent(in) :: pu, pv !< Point.
   integer(I4P)                   :: w      !< Winding number.
   integer(I4P)                   :: i, j   !< Counters.

   w = 0_I4P
   do i = self%loop_start(l), self%loop_start(l + 1) - 1
      j = self%next(i)
      if (self%v(i) <= pv) then
         if (self%v(j) > pv) then
            if (orient2d_xy(self%u(i), self%v(i), self%u(j), self%v(j), pu, pv) > 0._R8P) w = w + 1_I4P
         endif
      else
         if (self%v(j) <= pv) then
            if (orient2d_xy(self%u(i), self%v(i), self%u(j), self%v(j), pu, pv) < 0._R8P) w = w - 1_I4P
         endif
      endif
   enddo
   endfunction loop_winding

   pure subroutine compute_bodies(self)
   !< Bodies (issue #37): a closed loop at even nesting depth starts a body, one at odd depth (a
   !< hole) joins its immediate container, an open loop is a body of its own; an island inside a
   !< hole is a body of its own. Containment is tested (exact winding of the loop's first vertex)
   !< only when the container's box contains the loop's box.
   type(curve_object), intent(inout) :: self !< Curve.
   real(R8P),    allocatable         :: lo(:,:), hi(:,:) !< Loop boxes.
   integer(I4P), allocatable         :: depth(:), container(:), body_of(:) !< Nesting.
   integer(I4P)                      :: l, m, i, nb, f, la !< Work.

   allocate(lo(2, self%nl), hi(2, self%nl), depth(self%nl), container(self%nl), body_of(self%nl))
   allocate(self%loop_body_(self%nl))
   do l = 1, self%nl
      f = self%loop_start(l) ; la = self%loop_start(l + 1) - 1_I4P
      if (la < f) then
         lo(:, l) = 0._R8P ; hi(:, l) = 0._R8P
      else
         lo(:, l) = [minval(self%u(f:la)), minval(self%v(f:la))]
         hi(:, l) = [maxval(self%u(f:la)), maxval(self%v(f:la))]
      endif
   enddo
   depth = 0_I4P ; container = 0_I4P
   do l = 1, self%nl
      if (.not. self%loop_closed(l) .or. self%loop_start(l + 1) <= self%loop_start(l)) cycle
      i = self%loop_start(l)
      do m = 1, self%nl
         if (m == l .or. .not. self%loop_closed(m)) cycle
         if (any(lo(:, m) > lo(:, l)) .or. any(hi(:, m) < hi(:, l))) cycle
         if (loop_winding(self, m, self%u(i), self%v(i)) == 0_I4P) cycle
         depth(l) = depth(l) + 1_I4P
         if (container(l) == 0_I4P) then
            container(l) = m
         elseif (product(hi(:, m) - lo(:, m)) < product(hi(:, container(l)) - lo(:, container(l)))) then
            container(l) = m
         endif
      enddo
   enddo
   body_of = 0_I4P ; nb = 0_I4P
   do l = 1, self%nl
      if (mod(depth(l), 2_I4P) == 0_I4P) then
         nb = nb + 1_I4P ; body_of(l) = nb
      endif
   enddo
   do l = 1, self%nl
      if (body_of(l) > 0_I4P) then
         self%loop_body_(l) = body_of(l)
      else
         m = container(l)
         if (m > 0_I4P) then
            self%loop_body_(l) = body_of(m)
         else
            self%loop_body_(l) = 0_I4P
         endif
         if (self%loop_body_(l) == 0_I4P) then ! defensive: inconsistent nesting
            nb = nb + 1_I4P ; self%loop_body_(l) = nb
         endif
      endif
   enddo
   self%bodies_number = nb
   endsubroutine compute_bodies

   pure function shoelace(u, v) result(a)
   !< Signed area of the closed polygon (u, v).
   real(R8P), intent(in) :: u(:), v(:) !< Vertices.
   real(R8P)             :: a          !< Signed area.
   integer(I4P)          :: i, j       !< Counters.

   a = 0._R8P
   do i = 1, size(u)
      j = mod(i, size(u)) + 1
      a = a + 0.5_R8P * (u(i) * v(j) - u(j) * v(i))
   enddo
   endfunction shoelace

   pure function curve_signed_distance_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, &
                                           seg_id, next, u, v, pn_u, pn_v, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2, &
                                           pu, pv, signed, fwn) result(d)
   !< Distance of (pu, pv) on the flat arrays of a non-empty curve, signed by the pseudo-normal of
   !< the closest feature or (fwn) by the fast winding number with beta = 2. The same expressions
   !< as `compute_distance`. CPU and device kernel.
   !$acc routine seq
   real(R8P),    intent(in) :: lo_u(:), lo_v(:), hi_u(:), hi_v(:)  !< Node boxes.
   integer(I4P), intent(in) :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   real(R8P),    intent(in) :: a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Segment payload.
   integer(I4P), intent(in) :: seg_id(:), next(:)                   !< Payload -> segment, next vertex.
   real(R8P),    intent(in) :: u(:), v(:), pn_u(:), pn_v(:)         !< Vertices, vertex pseudo-normals.
   real(R8P),    intent(in) :: wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_c(:,:), wn_diag2(:) !< Winding cache.
   real(R8P),    intent(in) :: pu, pv                               !< Point.
   logical,      intent(in) :: signed                               !< Signed distance.
   logical,      intent(in) :: fwn                                  !< Sign by the fast winding number.
   real(R8P)                :: d                                    !< Distance.
   real(R8P)                :: d2, t, side                          !< Work.
   integer(I4P)             :: k, feature, sg, vx                   !< Work.

   call segment_closest_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, pu, pv, &
                            d2, k, t, feature)
   d = sqrt(d2)
   if (.not. signed) return
   if (fwn) then
      side = 0.5_R8P - segment_winding_number_oac(left, right, first, last, a_u, a_v, d_u, d_v, wn_cu, wn_cv, wn_du, wn_dv, &
                                                  wn_c, wn_diag2, pu, pv, 4._R8P)
   else
      sg = seg_id(k)
      vx = 0_I4P
      if (feature == FEATURE_VERTEX_START) vx = sg
      if (feature == FEATURE_VERTEX_END) vx = next(sg)
      if (vx > 0_I4P) then
         side = (pu - u(vx)) * pn_u(vx) + (pv - v(vx)) * pn_v(vx)
      else
         side = -(d_u(k) * (pv - a_v(k)) - d_v(k) * (pu - a_u(k)))
      endif
   endif
   if (side < 0._R8P) d = -d
   endfunction curve_signed_distance_oac

   pure subroutine curve_query_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, seg_id, &
                                   next, u, v, pn_u, pn_v, seg_nu, seg_nv, vbody, wn_cu, wn_cv, wn_du, wn_dv, wn_c, wn_diag2, &
                                   pu, pv, signed, fwn, d, cpu, cpv, sg, vx, nu, nv, b)
   !< Rich query of (pu, pv) on the flat arrays of a non-empty curve (issue #37): the expressions of
   !< `compute_distance`. CPU and device kernel.
   !$acc routine seq
   real(R8P),    intent(in)  :: lo_u(:), lo_v(:), hi_u(:), hi_v(:)  !< Node boxes.
   integer(I4P), intent(in)  :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   real(R8P),    intent(in)  :: a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Segment payload.
   integer(I4P), intent(in)  :: seg_id(:), next(:)                   !< Payload -> segment, next vertex.
   real(R8P),    intent(in)  :: u(:), v(:), pn_u(:), pn_v(:)         !< Vertices, vertex pseudo-normals.
   real(R8P),    intent(in)  :: seg_nu(:), seg_nv(:)                 !< Segment normals.
   integer(I4P), intent(in)  :: vbody(:)                             !< Body of each vertex/segment.
   real(R8P),    intent(in)  :: wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_c(:,:), wn_diag2(:) !< Winding cache.
   real(R8P),    intent(in)  :: pu, pv                               !< Point.
   logical,      intent(in)  :: signed, fwn                          !< Options.
   real(R8P),    intent(out) :: d, cpu, cpv, nu, nv                  !< Distance, closest point, normal.
   integer(I4P), intent(out) :: sg, vx, b                            !< Segment, vertex, body.
   real(R8P)                 :: d2, t, side                          !< Work.
   integer(I4P)              :: k, feature                           !< Work.

   call segment_closest_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, pu, pv, &
                            d2, k, t, feature)
   d = sqrt(d2)
   sg = seg_id(k)
   vx = 0_I4P
   if (feature == FEATURE_VERTEX_START) vx = sg
   if (feature == FEATURE_VERTEX_END) vx = next(sg)
   if (vx > 0_I4P) then
      cpu = u(vx) ; cpv = v(vx) ; nu = pn_u(vx) ; nv = pn_v(vx)
   else
      cpu = a_u(k) + t * d_u(k) ; cpv = a_v(k) + t * d_v(k) ; nu = seg_nu(sg) ; nv = seg_nv(sg)
   endif
   b = vbody(sg)
   if (.not. signed) return
   if (fwn) then
      side = 0.5_R8P - segment_winding_number_oac(left, right, first, last, a_u, a_v, d_u, d_v, wn_cu, wn_cv, wn_du, wn_dv, &
                                                  wn_c, wn_diag2, pu, pv, 4._R8P)
   elseif (vx > 0_I4P) then
      side = (pu - u(vx)) * pn_u(vx) + (pv - v(vx)) * pn_v(vx)
   else
      side = -(d_u(k) * (pv - a_v(k)) - d_v(k) * (pu - a_u(k)))
   endif
   if (side < 0._R8P) d = -d
   endsubroutine curve_query_oac

   pure subroutine second_foot_2d(self, p, cp1, d1, rel, foot_tol, found, cp2)
   !< Nearest segment whose closest point is farther than foot_tol from cp1 within (1 + rel) d1
   !< (issue #37): candidates from a box query of half-width (1 + rel) d1 around p.
   type(curve_object), intent(in)  :: self     !< Curve.
   real(R8P),          intent(in)  :: p(2)     !< Point.
   real(R8P),          intent(in)  :: cp1(2)   !< First foot point.
   real(R8P),          intent(in)  :: d1       !< Minimum distance.
   real(R8P),          intent(in)  :: rel      !< Relative tolerance.
   real(R8P),          intent(in)  :: foot_tol !< Foot separation.
   logical,            intent(out) :: found    !< Second foot point found.
   real(R8P),          intent(out) :: cp2(2)   !< Second foot point (cp1 if none).
   integer(I4P), allocatable       :: ids(:)   !< Candidates.
   real(R8P)                       :: r, bound, d2, t, c(2) !< Work.
   integer(I4P)                    :: m, j, k, f !< Work.

   found = .false. ; cp2 = cp1
   r = (1._R8P + rel) * d1
   bound = r * r
   call self%bvh%segments_in_box(p(1) - r, p(2) - r, p(1) + r, p(2) + r, ids, m)
   do j = 1, m
      k = ids(j)
      call segment_point_distance2(self%bvh%a_u(k), self%bvh%a_v(k), self%bvh%d_u(k), self%bvh%d_v(k), self%bvh%len2(k), &
                                   p(1), p(2), d2, t, f)
      if (d2 > bound) cycle
      c = [self%bvh%a_u(k) + t * self%bvh%d_u(k), self%bvh%a_v(k) + t * self%bvh%d_v(k)]
      if (hypot(c(1) - cp1(1), c(2) - cp1(2)) <= foot_tol) cycle
      found = .true. ; cp2 = c ; bound = d2
   enddo
   endsubroutine second_foot_2d

   pure function path_clear(self, a, b, sg, vx) result(clear)
   !< True if the segment a-b crosses no segment of the curve except those containing the foot point
   !< a: the closest segment sg (vx = 0) or the segments at the closest vertex vx. Exact.
   type(curve_object), intent(in) :: self  !< Curve.
   real(R8P),          intent(in) :: a(2)  !< Foot point.
   real(R8P),          intent(in) :: b(2)  !< Probe / image point.
   integer(I4P),       intent(in) :: sg    !< Closest segment.
   integer(I4P),       intent(in) :: vx    !< Closest vertex (0: segment interior).
   logical                        :: clear !< Result.
   integer(I4P), allocatable      :: ids(:) !< Candidates.
   integer(I4P)                   :: m, c, q1, q2 !< Work.

   clear = .true.
   call self%bvh%segments_in_box(min(a(1), b(1)), min(a(2), b(2)), max(a(1), b(1)), max(a(2), b(2)), ids, m)
   do c = 1, m
      q1 = self%bvh%seg_id(ids(c)) ; q2 = self%next(q1)
      if (vx == 0_I4P .and. q1 == sg) cycle
      if (vx > 0_I4P .and. (q1 == vx .or. q2 == vx)) cycle
      if (segments_intersect(a(1), a(2), b(1), b(2), self%u(q1), self%v(q1), self%u(q2), self%v(q2))) then
         clear = .false.
         return
      endif
   enddo
   endfunction path_clear

   pure function rich_sizes_ok(n, closest_points, segments, vertices, normals, bodies) result(good)
   !< True if every present rich output has the right size for n points.
   integer(I4P),           intent(in) :: n                   !< Points.
   real(R8P),    optional, intent(in) :: closest_points(:,:) !< [2, n].
   integer(I4P), optional, intent(in) :: segments(:)         !< [n].
   integer(I4P), optional, intent(in) :: vertices(:)         !< [n].
   real(R8P),    optional, intent(in) :: normals(:,:)        !< [2, n].
   integer(I4P), optional, intent(in) :: bodies(:)           !< [n].
   logical                            :: good                !< Result.

   good = .true.
   if (present(closest_points)) good = good .and. size(closest_points, 1) == 2 .and. size(closest_points, 2) == n
   if (present(normals)) good = good .and. size(normals, 1) == 2 .and. size(normals, 2) == n
   if (present(segments)) good = good .and. size(segments) == n
   if (present(vertices)) good = good .and. size(vertices) == n
   if (present(bodies)) good = good .and. size(bodies) == n
   endfunction rich_sizes_ok

   function resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed, fwn) result(good)
   !< Check the batch arguments and resolve the options; false (status set) if invalid.
   type(curve_object), intent(in)            :: self           !< Curve (unused: kept for the call symmetry).
   real(R8P),          intent(in)            :: points(:,:)    !< Points [2, n].
   real(R8P),          intent(in)            :: distances(:)   !< Distances [n].
   logical,            intent(in),  optional :: is_signed      !< Signed.
   integer(I4P),       intent(in),  optional :: sign_algorithm !< Sign algorithm.
   integer(I4P),       intent(out), optional :: status         !< Status.
   logical,            intent(out)           :: signed         !< Resolved signed.
   logical,            intent(out)           :: fwn            !< Resolved FWN sign.
   logical                                   :: good           !< Arguments valid.
   integer(I4P)                              :: algo           !< Resolved algorithm.

   if (present(status)) status = STATUS_OK
   signed = .true. ; if (present(is_signed)) signed = is_signed
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   fwn = algo == SIGN_FAST_WINDING_NUMBER
   good = size(points, dim=1) == 2 .and. size(points, dim=2) == size(distances) .and. &
          (algo == SIGN_PSEUDO_NORMAL .or. algo == SIGN_FAST_WINDING_NUMBER)
   if (.not. good .and. present(status)) status = STATUS_INVALID_INPUT
   endfunction resolve_batch_options

   pure function segments_intersect(p1u, p1v, p2u, p2v, q1u, q1v, q2u, q2v) result(hit)
   !< True if the closed segments p1p2 and q1q2 share a point (exact).
   real(R8P), intent(in) :: p1u, p1v, p2u, p2v, q1u, q1v, q2u, q2v !< Endpoints.
   logical               :: hit                                     !< Result.
   real(R8P)             :: o1, o2, o3, o4                           !< Orientations.

   o1 = orient2d_xy(p1u, p1v, p2u, p2v, q1u, q1v)
   o2 = orient2d_xy(p1u, p1v, p2u, p2v, q2u, q2v)
   o3 = orient2d_xy(q1u, q1v, q2u, q2v, p1u, p1v)
   o4 = orient2d_xy(q1u, q1v, q2u, q2v, p2u, p2v)
   hit = opposite(o1, o2) .and. opposite(o3, o4)
   if (hit) return
   hit = (o1 == 0._R8P .and. on_box(q1u, q1v, p1u, p1v, p2u, p2v)) .or. &
         (o2 == 0._R8P .and. on_box(q2u, q2v, p1u, p1v, p2u, p2v)) .or. &
         (o3 == 0._R8P .and. on_box(p1u, p1v, q1u, q1v, q2u, q2v)) .or. &
         (o4 == 0._R8P .and. on_box(p2u, p2v, q1u, q1v, q2u, q2v))
   endfunction segments_intersect

   pure function folds_back(au, av, bu, bv, cu, cv) result(hit)
   !< For adjacent segments ab and bc: true if c is collinear with ab and bc runs back along ab.
   real(R8P), intent(in) :: au, av, bu, bv, cu, cv !< Points.
   logical               :: hit                    !< Result.

   hit = .false.
   if (orient2d_xy(au, av, bu, bv, cu, cv) /= 0._R8P) return
   hit = (cu - bu) * (au - bu) + (cv - bv) * (av - bv) > 0._R8P
   endfunction folds_back

   pure function opposite(a, b) result(yes)
   !< True if a and b are nonzero with opposite signs.
   real(R8P), intent(in) :: a, b !< Values.
   logical               :: yes  !< Result.

   yes = (a > 0._R8P .and. b < 0._R8P) .or. (a < 0._R8P .and. b > 0._R8P)
   endfunction opposite

   pure function on_box(pu, pv, au, av, bu, bv) result(yes)
   !< True if p lies in the closed bounding box of segment ab (p collinear with ab: then on ab).
   real(R8P), intent(in) :: pu, pv, au, av, bu, bv !< Points.
   logical               :: yes                    !< Result.

   yes = pu >= min(au, bu) .and. pu <= max(au, bu) .and. pv >= min(av, bv) .and. pv <= max(av, bv)
   endfunction on_box
endmodule fossil_curve_object
