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
use fossil_box_query,          only : curve_box_contact_oac, curve_box_thin_oac, cell_bands, cell_components, BOX_NONE, &
                                     BOX_CUT, BOX_INSIDE, BOX_OUTSIDE
use fossil_line_query,         only : curve_line_crossings_oac
use fossil_cut_cells_2d,       only : cut_cells_2d_t, compute_cut_cells_2d
use fossil_rigid_pose,         only : rigid_pose_2d_t, to_body, to_world, rotate
use fossil_predicates,         only : orient2d_xy
use fossil_mass_properties,    only : moments_accumulator_t, accumulate_segment, moments_to_properties_2d, symmetric2_eigen
use fossil_feature_edges,      only : FEATURE_EDGE_DEFAULT_ANGLE, feature_cos
use fossil_thickness,          only : THICKNESS_RAY, THICKNESS_MEDIAL, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS, &
                                     THICKNESS_DEFAULT_SEPARATION, lower_median
use fossil_deviation,          only : deviation_evaluator_2d_t, deviation_sampled_2d, deviation_certified_2d, &
                                     DEV_STATUS_OK, DEV_STATUS_BAD_INPUT
use fossil_segment_bvh,        only : segment_bvh, FEATURE_VERTEX_START, FEATURE_VERTEX_END, segment_point_distance2, &
                                     segment_closest_oac, segment_winding_number_oac
use fossil_surface_stl_object, only : STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL, &
                                     SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use penf,                      only : I4P, R8P
use, intrinsic :: iso_fortran_env, only : stderr => error_unit

implicit none
private
public :: curve_object
public :: cut_cells_2d_t
public :: CURVE_VALID, CURVE_TOO_FEW_VERTICES, CURVE_ZERO_LENGTH_SEGMENT, CURVE_ZERO_AREA_LOOP, CURVE_WRONG_ORIENTATION
public :: CURVE_OPEN_LOOP, CURVE_SELF_INTERSECTING

! Validity codes returned by `check_validity` (the first defect found, loop by loop).
integer(I4P), parameter :: CURVE_VALID               = 0_I4P !< Valid curve.
integer(I4P), parameter :: CURVE_TOO_FEW_VERTICES    = 1_I4P !< A loop has fewer than 3 vertices.
integer(I4P), parameter :: CURVE_ZERO_LENGTH_SEGMENT = 2_I4P !< Two consecutive vertices coincide.
integer(I4P), parameter :: CURVE_ZERO_AREA_LOOP      = 3_I4P !< A loop has zero signed area.
integer(I4P), parameter :: CURVE_WRONG_ORIENTATION   = 4_I4P !< A loop's orientation does not match its nesting depth.
integer(I4P), parameter :: CURVE_OPEN_LOOP           = 5_I4P !< A loop is open (a polyline).
integer(I4P), parameter :: CURVE_SELF_INTERSECTING   = 6_I4P !< Two segments cross, touch or overlap (only with
                                                             !< `check_self_intersections`, issue #44).

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
   integer(I4P), allocatable :: seg_patch(:)         !< Patch of each segment [nv] (issue #66; segment i = i -> next(i)).
   character(64), allocatable :: patch_names(:)      !< Patch names (default one unnamed patch).
   real(R8P)                 :: bmin(2) = 0._R8P    !< Bounding box minimum (u, v).
   real(R8P)                 :: bmax(2) = 0._R8P    !< Bounding box maximum (u, v).
   real(R8P)                 :: area = 0._R8P       !< Total signed area (solid area for a valid curve).
   real(R8P)                 :: perimeter = 0._R8P  !< Total length.
   real(R8P)                 :: centroid(2) = 0._R8P !< Centroid of the solid region (u, v).
   type(segment_bvh)         :: bvh                 !< SAH BVH over the segments.
   logical                   :: is_on_device = .false. !< Arrays mapped by enter_device.
   real(R8P),    allocatable :: thin_pay(:)         !< Local size per segment-tree payload entry (issue #62).
   real(R8P),    allocatable :: thin_node(:)        !< Minimum local size per tree node.
   logical                   :: thin_built = .false. !< build_thin ran on the current curve.
   logical                   :: thin_on_device = .false. !< The thin arrays were mapped by enter_device.
   contains
      ! public methods
      procedure, pass(self) :: adopt_loops          !< Build the curve from vertex arrays and loop sizes.
      procedure, pass(self) :: get_patches_number   !< Number of patches (issue #66).
      procedure, pass(self) :: patch_name           !< Name of a patch.
      procedure, pass(self) :: segment_patch        !< Patch of a segment.
      procedure, pass(self) :: set_patch            !< Assign segments to a named patch.
      procedure, pass(self) :: distance_to          !< Distance to another curve, with the closest pair (issue #67).
      procedure, pass(self) :: box_classify_many    !< Classify boxes: cut, near, inside, outside (issue #38).
      procedure, pass(self) :: box_intersects_curve !< True if a box meets the curve (issue #38).
      procedure, pass(self) :: box_classify_many_device !< box_classify_many on the device (issue #38).
      procedure, pass(self) :: tag_block            !< Signed tag bands on a block of cells (issue #38).
      procedure, pass(self) :: line_crossings       !< Exact crossings of an axis-aligned line (issue #40).
      procedure, pass(self) :: cell_wall_fractions  !< Wall fractions along the grid lines of a block (issue #40).
      procedure, pass(self) :: line_crossings_many  !< Crossings of many parallel lines, fixed buffers (issue #40).
      procedure, pass(self) :: line_crossings_many_device !< line_crossings_many on the device (issue #40).
      procedure, pass(self) :: cut_cells            !< Exact cut cells of a block (issue #40).
      procedure, pass(self) :: move                 !< Move rigidly by an in-plane pose, refit the tree (issue #41).
      procedure, pass(self) :: state_changes        !< Fresh and dead points between two poses (issue #41).
      procedure, pass(self) :: wall_distance_many   !< Wall distance, body id, cutoff (issue #42).
      procedure, pass(self) :: curve_markers        !< Exact arc-length resampling: markers, lengths, normals (issue #42).
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
      procedure, pass(self) :: get_next             !< End vertex of a segment (issue #41).
      procedure, pass(self) :: get_vertices_number  !< Number of vertices (= segments).
      procedure, pass(self) :: is_loop_closed       !< True if a loop is closed.
      procedure, pass(self) :: is_point_inside      !< Inside test.
      procedure, pass(self) :: load_from_file       !< Load from the ASCII loop format.
      procedure, pass(self) :: mass_properties      !< Mass, centroid, second moments, principal axes (issue #59).
      procedure, pass(self) :: deviation            !< One-sided deviation (Hausdorff) from another curve (issue #60).
      procedure, pass(self) :: hausdorff            !< Symmetric Hausdorff distance to another curve (issue #60).
      procedure, pass(self) :: feature_corners      !< Corners: turning above an angle, ends of open loops (issue #61).
      procedure, pass(self) :: thickness => curve_thickness !< Local thickness and gap per segment (issue #62).
      procedure, pass(self) :: thin_regions => curve_thin_regions !< Segments thinner than factor*h, in runs (issue #62).
      procedure, pass(self) :: build_thin => curve_build_thin !< Local sizes for the thin tagging (issue #62).
      procedure, pass(self) :: box_thin_many => curve_box_thin_many !< Boxes meeting a thin segment (issue #62).
      procedure, pass(self) :: box_thin_many_device => curve_box_thin_many_device !< box_thin_many on the device (#62).
      procedure, pass(self) :: nesting_depth        !< Nesting depth of a loop (exact).
      procedure, pass(self) :: repair               !< Join/close chains, merge vertices, drop degenerate loops, orient.
      procedure, pass(self) :: save_into_file       !< Save into the ASCII loop format.
      procedure, pass(self) :: self_intersections   !< Pairs of segments that cross, touch or overlap (exact).
      procedure, pass(self) :: winding_number       !< Generalized winding number (fast or brute force).
      procedure, pass(self) :: winding_number_exact !< Exact integer winding number of a point (brute force).
      ! private methods
      procedure, pass(self), private :: analyze     !< Derived data: topology, metrics, normals.
endtype curve_object

type, extends(deviation_evaluator_2d_t) :: curve_evaluator_t
   !< A curve seen through its batched closest-point query, for `fossil_deviation` (issue #60).
   type(curve_object), pointer :: curve => null() !< The curve B.
   logical                     :: device = .false. !< Query on the device (query_many_device).
   contains
      procedure, pass(self) :: eval => curve_evaluator_eval !< Unsigned distances, closest points and segments.
endtype curve_evaluator_t

contains
   ! public methods
   pure function get_patches_number(self) result(n)
   !< Number of patches (issue #66): one per `# patch <name>` directive met, 1 by default.
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P)                    :: n    !< Patches.

   n = 1_I4P
   if (allocated(self%patch_names)) n = max(1_I4P, int(size(self%patch_names), I4P))
   endfunction get_patches_number

   pure function patch_name(self, i) result(name)
   !< Name of patch i ('' if unnamed or out of range).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P),        intent(in) :: i    !< Patch.
   character(len=:), allocatable   :: name !< Name.

   name = ''
   if (allocated(self%patch_names)) then
      if (i >= 1 .and. i <= size(self%patch_names)) name = trim(self%patch_names(i))
   endif
   endfunction patch_name

   elemental function segment_patch(self, i) result(p)
   !< Patch of segment i (i -> next(i)); 0 out of range.
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P),        intent(in) :: i    !< Segment.
   integer(I4P)                    :: p    !< Patch.

   p = 0_I4P
   if (allocated(self%seg_patch)) then
      if (i >= 1 .and. i <= size(self%seg_patch)) p = self%seg_patch(i)
   endif
   endfunction segment_patch

   subroutine distance_to(self, other, distance, point_self, point_other, segment_self, segment_other, intersecting, &
                          status, pose_self, pose_other)
   !< Distance between this curve and `other` (issue #67, the 2-D twin of `surface%distance_to`): the minimum over their
   !< segment pairs, 0 with `intersecting` when two segments cross or touch (exact orient2d), else the minimum of the
   !< four endpoint-segment distances; the closest segments (ties: the smallest pair) and the witness points. Every pair
   !< is evaluated in order behind a box pre-filter (curves are small: no tree descent). Poses (`rigid_pose_2d_t`) place
   !< the curves as `move` would (to_world of each vertex).
   class(curve_object),   intent(in)            :: self, other    !< Curves.
   real(R8P),             intent(out)           :: distance       !< Distance.
   real(R8P),             intent(out), optional :: point_self(2), point_other(2) !< Witness points.
   integer(I4P),          intent(out), optional :: segment_self, segment_other   !< Closest segments.
   logical,               intent(out), optional :: intersecting   !< The curves touch or cross.
   integer(I4P),          intent(out), optional :: status         !< STATUS_OK, STATUS_INVALID_INPUT (empty curve).
   type(rigid_pose_2d_t), intent(in),  optional :: pose_self, pose_other !< Poses.
   real(R8P)                                    :: d2, best, a0(2), a1(2), b0(2), b1(2), pa(2), pb(2), qa(2), qb(2)
   integer(I4P)                                 :: i, j, ia, ib
   logical                                      :: meet, hit

   if (present(status)) status = STATUS_OK
   distance = huge(1._R8P)
   if (self%nv == 0 .or. other%nv == 0) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   best = huge(1._R8P) ; ia = 0 ; ib = 0 ; meet = .false. ; pa = 0._R8P ; pb = 0._R8P
   do i = 1, self%nv
      if (self%next(i) == 0) cycle
      a0 = place(self, i, pose_self) ; a1 = place(self, self%next(i), pose_self)
      do j = 1, other%nv
         if (other%next(j) == 0) cycle
         b0 = place(other, j, pose_other) ; b1 = place(other, other%next(j), pose_other)
         if (max(0._R8P, min(b0(1), b1(1)) - max(a0(1), a1(1)), min(a0(1), a1(1)) - max(b0(1), b1(1)))**2 + &
             max(0._R8P, min(b0(2), b1(2)) - max(a0(2), a1(2)), min(a0(2), a1(2)) - max(b0(2), b1(2)))**2 > &
             best * (1._R8P + 1.e-10_R8P)) cycle
         call segment_pair(a0, a1, b0, b1, d2, qa, qb, hit)
         if (d2 < best .or. (d2 == best .and. (i < ia .or. (i == ia .and. j < ib)))) then
            best = d2 ; ia = i ; ib = j ; pa = qa ; pb = qb ; meet = hit
         endif
      enddo
   enddo
   distance = sqrt(best)
   if (present(point_self)) point_self = pa
   if (present(point_other)) point_other = pb
   if (present(segment_self)) segment_self = ia
   if (present(segment_other)) segment_other = ib
   if (present(intersecting)) intersecting = meet
   contains
      function place(c, k, pose) result(x)
      class(curve_object),             intent(in) :: c
      integer(I4P),                    intent(in) :: k
      type(rigid_pose_2d_t), optional, intent(in) :: pose
      real(R8P)                                   :: x(2)

      x = [c%u(k), c%v(k)]
      if (present(pose)) x = to_world(pose, x)
      endfunction place
   endsubroutine distance_to

   subroutine segment_pair(a0, a1, b0, b1, d2, pa, pb, meet)
   !< Squared distance between segments (a0, a1) and (b0, b1) with witness points; 0 and `meet` when they cross or touch
   !< (exact orient2d, collinear overlaps included).
   real(R8P), intent(in)  :: a0(2), a1(2), b0(2), b1(2)
   real(R8P), intent(out) :: d2, pa(2), pb(2)
   logical,   intent(out) :: meet
   integer(I4P)           :: s1, s2, s3, s4
   real(R8P)              :: c(2), dd, x(2)
   integer(I4P)           :: k

   s1 = sgn(orient2d_xy(b0(1), b0(2), b1(1), b1(2), a0(1), a0(2)))
   s2 = sgn(orient2d_xy(b0(1), b0(2), b1(1), b1(2), a1(1), a1(2)))
   s3 = sgn(orient2d_xy(a0(1), a0(2), a1(1), a1(2), b0(1), b0(2)))
   s4 = sgn(orient2d_xy(a0(1), a0(2), a1(1), a1(2), b1(1), b1(2)))
   meet = (s1 * s2 < 0 .and. s3 * s4 < 0) .or. (s1 == 0 .and. within(b0, b1, a0)) .or. (s2 == 0 .and. within(b0, b1, a1)) &
          .or. (s3 == 0 .and. within(a0, a1, b0)) .or. (s4 == 0 .and. within(a0, a1, b1))
   if (meet) then
      d2 = 0._R8P
      if (s1 == 0 .and. within(b0, b1, a0)) then
         x = a0
      elseif (s2 == 0 .and. within(b0, b1, a1)) then
         x = a1
      elseif (s3 == 0 .and. within(a0, a1, b0)) then
         x = b0
      elseif (s4 == 0 .and. within(a0, a1, b1)) then
         x = b1
      else ! proper crossing: the point on (a0, a1)
         x = a0 + (orient2d_xy(b0(1), b0(2), b1(1), b1(2), a0(1), a0(2)) / &
             (orient2d_xy(b0(1), b0(2), b1(1), b1(2), a0(1), a0(2)) - orient2d_xy(b0(1), b0(2), b1(1), b1(2), a1(1), a1(2)))) * &
             (a1 - a0)
      endif
      pa = x ; pb = x
      return
   endif
   d2 = huge(1._R8P)
   do k = 1, 4
      select case (k)
      case (1) ; c = closest(b0, b1, a0) ; dd = sum((a0 - c)**2) ; if (dd < d2) then ; d2 = dd ; pa = a0 ; pb = c ; endif
      case (2) ; c = closest(b0, b1, a1) ; dd = sum((a1 - c)**2) ; if (dd < d2) then ; d2 = dd ; pa = a1 ; pb = c ; endif
      case (3) ; c = closest(a0, a1, b0) ; dd = sum((b0 - c)**2) ; if (dd < d2) then ; d2 = dd ; pa = c ; pb = b0 ; endif
      case (4) ; c = closest(a0, a1, b1) ; dd = sum((b1 - c)**2) ; if (dd < d2) then ; d2 = dd ; pa = c ; pb = b1 ; endif
      endselect
   enddo
   contains
      pure function sgn(x) result(s)
      real(R8P), intent(in) :: x
      integer(I4P)          :: s

      s = 0_I4P ; if (x > 0._R8P) s = 1_I4P ; if (x < 0._R8P) s = -1_I4P
      endfunction sgn

      pure function within(p, q, r) result(yes)
      real(R8P), intent(in) :: p(2), q(2), r(2)
      logical               :: yes

      yes = r(1) >= min(p(1), q(1)) .and. r(1) <= max(p(1), q(1)) .and. r(2) >= min(p(2), q(2)) .and. r(2) <= max(p(2), q(2))
      endfunction within

      pure function closest(p, q, r) result(x)
      real(R8P), intent(in) :: p(2), q(2), r(2)
      real(R8P)             :: x(2), d(2), l2, t

      d = q - p ; l2 = sum(d * d)
      t = 0._R8P ; if (l2 > 0._R8P) t = min(1._R8P, max(0._R8P, sum((r - p) * d) / l2))
      x = p + t * d
      endfunction closest
   endsubroutine segment_pair

   subroutine set_patch(self, segments, name)
   !< Assign the segments `segments` to the patch named `name`, created if new (issue #66).
   class(curve_object), intent(inout) :: self        !< Curve.
   integer(I4P),        intent(in)    :: segments(:) !< Segment ids.
   character(*),        intent(in)    :: name        !< Patch name.
   character(64)                      :: nm          !< Padded name.
   integer(I4P)                       :: p, k        !< Patch, counter.

   if (.not. allocated(self%seg_patch)) allocate(self%seg_patch(self%nv), source=1_I4P)
   if (.not. allocated(self%patch_names)) then
      allocate(self%patch_names(1)) ; self%patch_names = ''
   endif
   nm = name
   p = 0_I4P
   do k = 1, size(self%patch_names)
      if (self%patch_names(k) == nm) p = k
   enddo
   if (p == 0_I4P) then
      self%patch_names = [self%patch_names, nm] ; p = size(self%patch_names)
   endif
   do k = 1, size(segments)
      if (segments(k) >= 1 .and. segments(k) <= self%nv) self%seg_patch(segments(k)) = p
   enddo
   endsubroutine set_patch

   subroutine adopt_loops(self, u, v, loop_sizes, normal_axis, closed, status, patches, patch_names)
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
   integer(I4P),        intent(in),  optional :: patches(:)    !< Patch of each segment (by its first vertex), issue #66.
   character(*),        intent(in),  optional :: patch_names(:) !< Patch names.
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
   ! patches (issue #66): one unnamed patch unless given
   allocate(self%seg_patch(self%nv), source=1_I4P)
   if (present(patches)) then
      if (size(patches) == self%nv) self%seg_patch = max(1_I4P, patches)
   endif
   if (present(patch_names)) then
      allocate(self%patch_names(max(1, size(patch_names))))
      self%patch_names = ''
      if (size(patch_names) > 0) self%patch_names(1:size(patch_names)) = patch_names
   else
      allocate(self%patch_names(1)) ; self%patch_names = ''
   endif
   do while (size(self%patch_names) < maxval([1_I4P, self%seg_patch]))
      self%patch_names = [self%patch_names, repeat(' ', 64)]
   enddo
   endsubroutine adopt_loops

   subroutine box_classify_many(self, lo, hi, codes, radius, sign_algorithm, inside, status)
   !< Classify boxes (u, v) against the curve (issue #38), OpenMP over the boxes.
   !<
   !< `codes(i)` is BOX_CUT if box i meets a segment (exact, `orient2d`), BOX_NEAR if not but a
   !< segment is within `radius` of it (default 0: no BOX_NEAR; floating-point segment-box distance),
   !< else BOX_INSIDE or BOX_OUTSIDE from `is_point_inside` at the box centre with `sign_algorithm`
   !< (default SIGN_FAST_WINDING_NUMBER). `inside(i)` (optional): the centre sign of every box that is not
   !< cut (false for cut boxes). Segment BVH traversal (`curve_box_contact_oac`, shared with the
   !< device). `lo(2, n)`, `hi(2, n)`; `status = STATUS_INVALID_INPUT` on a size mismatch or a
   !< negative radius. An empty curve: every box outside.
   class(curve_object), intent(in)            :: self           !< Curve.
   real(R8P),           intent(in)            :: lo(:,:)        !< Box minima [2, n].
   real(R8P),           intent(in)            :: hi(:,:)        !< Box maxima [2, n].
   integer(I4P),        intent(out)           :: codes(:)       !< BOX_* [n].
   real(R8P),           intent(in),  optional :: radius         !< Near radius (default 0).
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_* for the centre query (default SIGN_FAST_WINDING_NUMBER).
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
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
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
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_FAST_WINDING_NUMBER (default) or SIGN_PSEUDO_NORMAL.
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
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
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

   subroutine tag_block(self, x, y, band, tags, sign_algorithm, status, feature_band, feature_tags, feature_angle, &
                        thin_factor, thin_band, thin_tags)
   !< Corner tags (issue #61): with `feature_band` and `feature_tags(nx, ny)`, the Chebyshev cell distance to the nearest
   !< cell holding a corner (`feature_corners`, `feature_angle`), capped at feature_band + 1.
   !< Thin tags (issue #62): with `thin_factor`, `thin_band` and `thin_tags(nx, ny)`, the Chebyshev cell distance to the
   !< nearest cell meeting a segment thinner than thin_factor times the cell's larger side (`box_thin_many`), capped at
   !< thin_band + 1; requires `build_thin` (else STATUS_INVALID_INPUT).
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
   integer(I4P),        intent(in),  optional :: sign_algorithm !< SIGN_* (default SIGN_FAST_WINDING_NUMBER).
   integer(I4P),        intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),        intent(in),  optional :: feature_band   !< Corner band width (>= 0, issue #61).
   integer(I4P),        intent(out), optional :: feature_tags(:,:) !< Chebyshev cells to the nearest corner cell [nx, ny].
   real(R8P),           intent(in),  optional :: feature_angle  !< Corner threshold [degrees] (default 30).
   real(R8P),           intent(in),  optional :: thin_factor    !< Thin threshold factor (issue #62).
   integer(I4P),        intent(in),  optional :: thin_band      !< Thin band width (>= 0).
   integer(I4P),        intent(out), optional :: thin_tags(:,:) !< Thin tags [nx, ny].
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
   if (present(feature_tags)) then ! cells to the nearest corner (issue #61)
      bad = .not. present(feature_band)
      if (.not. bad) bad = feature_band < 0_I4P .or. size(feature_tags, 1) /= nx .or. size(feature_tags, 2) /= ny
      if (bad) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      block
         integer(I4P), allocatable :: cn(:), fd(:,:,:)
         logical,      allocatable :: hit(:,:,:)
         real(R8P)                 :: t_
         integer(I4P)              :: c_, ii, jj
         call self%feature_corners(cn, angle=feature_angle)
         allocate(hit(nx, ny, 1), fd(nx, ny, 1))
         hit = .false.
         t_ = 1.e-12_R8P * max(norm2(self%bmax - self%bmin), maxval(abs(self%bmin)), maxval(abs(self%bmax)))
         do c_ = 1, size(cn)                  ! every cell whose closed box (inflated by t_) holds the corner
            do jj = 1_I4P, ny
               if (self%v(cn(c_)) < y(jj) - t_ .or. self%v(cn(c_)) > y(jj+1) + t_) cycle
               do ii = 1_I4P, nx
                  if (self%u(cn(c_)) < x(ii) - t_ .or. self%u(cn(c_)) > x(ii+1) + t_) cycle
                  hit(ii, jj, 1) = .true.
               enddo
            enddo
         enddo
         call cell_bands(hit, feature_band, fd)
         feature_tags = fd(:, :, 1)
      endblock
   endif
   if (present(thin_tags)) then ! cells meeting a thin segment (issue #62)
      bad = .not. (present(thin_band) .and. present(thin_factor) .and. self%thin_built)
      if (.not. bad) bad = thin_band < 0_I4P .or. thin_factor <= 0._R8P .or. size(thin_tags, 1) /= nx .or. &
                           size(thin_tags, 2) /= ny
      if (bad) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      block
         real(R8P),    allocatable :: blo(:,:), bhi(:,:)
         logical,      allocatable :: hit(:)
         integer(I4P), allocatable :: td(:,:,:)
         integer(I4P)              :: c_, ii, jj
         allocate(blo(2, nx * ny), bhi(2, nx * ny), hit(nx * ny), td(nx, ny, 1))
         c_ = 0
         do jj = 1_I4P, ny
            do ii = 1_I4P, nx
               c_ = c_ + 1 ; blo(:, c_) = [x(ii), y(jj)] ; bhi(:, c_) = [x(ii+1), y(jj+1)]
            enddo
         enddo
         call self%box_thin_many(blo, bhi, thin_factor, hit)
         call cell_bands(reshape(hit, [nx, ny, 1]), thin_band, td)
         thin_tags = td(:, :, 1)
      endblock
   endif
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
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

   subroutine curve_thickness(self, thickness, gap, method, cone_angle, rays, invalid, status, separation_angle)
   !< Local thickness and gap of every segment (issue #62), the 2-D twin of `surface_stl_object%thickness`, measured from
   !< the segment midpoint p along its outward normal n; `huge` where undefined, where nothing is in front, and on the
   !< missing closing segment of an open loop.
   !< - THICKNESS_RAY: `rays` rays in a fan of full aperture `cone_angle` [degrees] around -n (thickness) or +n (gap):
   !<   ray 1 on the axis, then pairs at +-theta_j, theta_j = (cone_angle/2) j / (rays/2). First hit by the segment tree
   !<   (the source segment excluded); counted as for surfaces (leaving the solid for the thickness, entering or missing
   !<   for the gap), lower median, huge with fewer than half counted (`invalid`).
   !< - THICKNESS_MEDIAL: the shrinking disc, started at min(R0, t_axis/2), shrunk through the closest point (batched
   !<   `query_many`), stopped below `separation_angle` (default 30 degrees) as for surfaces.
   class(curve_object),    intent(in)            :: self             !< Curve.
   real(R8P), allocatable, intent(out), optional :: thickness(:)     !< Thickness [nv].
   real(R8P), allocatable, intent(out), optional :: gap(:)           !< Gap [nv].
   integer(I4P),           intent(in),  optional :: method           !< THICKNESS_RAY (default) or THICKNESS_MEDIAL.
   real(R8P),              intent(in),  optional :: cone_angle       !< Full aperture [degrees] (default 60).
   integer(I4P),           intent(in),  optional :: rays             !< Rays per segment (default 31).
   integer(I4P),           intent(out), optional :: invalid          !< Segments with an undefined value.
   integer(I4P),           intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),              intent(in),  optional :: separation_angle !< Medial denoising angle [degrees] (default 30).
   real(R8P), allocatable                        :: th(:), gp(:)     !< Work values.
   logical,   allocatable                        :: bad(:)           !< Undefined segments.
   real(R8P)                                     :: cone, sep        !< Options.
   integer(I4P)                                  :: meth, m, i       !< Options, counter.
   logical                                       :: want_t, want_g   !< Requested outputs.

   if (present(invalid)) invalid = 0
   if (present(status)) status = STATUS_OK
   meth = THICKNESS_RAY ; if (present(method)) meth = method
   cone = THICKNESS_DEFAULT_CONE ; if (present(cone_angle)) cone = cone_angle
   m = THICKNESS_DEFAULT_RAYS ; if (present(rays)) m = rays
   sep = THICKNESS_DEFAULT_SEPARATION ; if (present(separation_angle)) sep = separation_angle
   want_t = present(thickness) ; want_g = present(gap)
   if (self%nv == 0 .or. self%bvh%nodes_number == 0 .or. (meth /= THICKNESS_RAY .and. meth /= THICKNESS_MEDIAL) .or. &
       m < 1 .or. cone < 0._R8P .or. cone >= 180._R8P .or. sep < 0._R8P .or. sep >= 180._R8P .or. &
       .not. (want_t .or. want_g)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(th(self%nv), gp(self%nv), bad(self%nv))
   th = huge(1._R8P) ; gp = huge(1._R8P) ; bad = .false.
   if (meth == THICKNESS_RAY) then
      !$omp parallel do default(shared) schedule(dynamic, 16) private(i)
      do i = 1, self%nv
         if (self%next(i) > 0) call ray_segment(i)
      enddo
      !$omp end parallel do
   else
      if (want_t) call medial(-1._R8P, th)
      if (want_g) call medial(1._R8P, gp)
   endif
   if (present(invalid)) invalid = count(bad)
   if (present(thickness)) call move_alloc(th, thickness)
   if (present(gap)) call move_alloc(gp, gap)
   contains
      subroutine ray_segment(s)
      !< Thickness and gap of segment s by the ray method.
      integer(I4P), intent(in) :: s                             !< Segment.
      real(R8P)                :: p(2), n(2), d(2), buf(m), t, th0 !< Work.
      integer(I4P)             :: k, nv_, nmiss, hs             !< Counters, hit segment.
      real(R8P)                :: side                          !< -1 thickness, +1 gap.
      integer(I4P)             :: pass_                         !< Pass.

      p = 0.5_R8P * [self%u(s) + self%u(self%next(s)), self%v(s) + self%v(self%next(s))]
      n = [self%seg_nu(s), self%seg_nv(s)]
      if (n(1) ** 2 + n(2) ** 2 < 0.25_R8P) then
         bad(s) = .true. ; return
      endif
      do pass_ = 1, 2
         if (pass_ == 1 .and. .not. want_t) cycle
         if (pass_ == 2 .and. .not. want_g) cycle
         side = merge(-1._R8P, 1._R8P, pass_ == 1)
         nv_ = 0 ; nmiss = 0
         do k = 1, m
            th0 = fan_angle(k)
            d = side * [cos(th0) * n(1) - sin(th0) * n(2), sin(th0) * n(1) + cos(th0) * n(2)]
            call ray_first_2d(self, p, d, s, t, hs)
            if (hs == 0) then
               nmiss = nmiss + 1 ; cycle
            endif
            if (side < 0._R8P) then
               if (self%seg_nu(hs) * d(1) + self%seg_nv(hs) * d(2) <= 0._R8P) cycle
            else
               if (self%seg_nu(hs) * d(1) + self%seg_nv(hs) * d(2) >= 0._R8P) cycle
            endif
            nv_ = nv_ + 1 ; buf(nv_) = t
         enddo
         if (side < 0._R8P) then
            if (2 * nv_ < m) then
               bad(s) = .true.
            else
               th(s) = lower_median(buf(1:nv_))
            endif
         else
            if (2 * (nv_ + nmiss) < m) then
               bad(s) = .true.
            else
               buf(nv_ + 1:nv_ + nmiss) = huge(1._R8P)
               gp(s) = lower_median(buf(1:nv_ + nmiss))
            endif
         endif
      enddo
      endsubroutine ray_segment

      pure function fan_angle(k) result(a)
      !< Angle of ray k of the fan [rad]: 0, then +-theta_j.
      integer(I4P), intent(in) :: k !< Ray.
      real(R8P)                :: a !< Angle.

      a = 0._R8P
      if (k > 1) a = merge(1._R8P, -1._R8P, mod(k, 2) == 0) * 0.5_R8P * cone * (acos(-1._R8P) / 180._R8P) * &
                     real(k / 2, R8P) / real(m / 2, R8P)
      endfunction fan_angle

      subroutine medial(side, val)
      !< Shrinking disc on the side `side` n of every segment (-1 thickness, +1 gap); val = 2r.
      real(R8P), intent(in)    :: side   !< -1 or +1.
      real(R8P), intent(inout) :: val(:) !< Diameters [nv] (huge on entry).
      integer(I4P), parameter  :: MAXIT = 64 !< Iteration cap.
      real(R8P),    allocatable :: r(:), pc(:,:), dir(:,:), dq(:), q(:,:), c(:,:) !< Work.
      integer(I4P), allocatable :: act(:)   !< Segments not converged.
      real(R8P)                 :: r0, rn, cos_sep, w(2), t !< Work.
      integer(I4P)              :: na, it, k, s, nk, hs !< Counters.

      cos_sep = cos(sep * acos(-1._R8P) / 180._R8P)
      r0 = norm2(self%bmax - self%bmin)
      allocate(r(self%nv), pc(2, self%nv), dir(2, self%nv), act(self%nv))
      na = 0
      do s = 1, self%nv
         r(s) = r0 ; pc(:, s) = 0._R8P ; dir(:, s) = 0._R8P
         if (self%next(s) == 0) cycle
         pc(:, s) = 0.5_R8P * [self%u(s) + self%u(self%next(s)), self%v(s) + self%v(self%next(s))]
         dir(:, s) = side * [self%seg_nu(s), self%seg_nv(s)]
         if (dir(1, s) ** 2 + dir(2, s) ** 2 < 0.25_R8P) then
            bad(s) = .true.
         else
            na = na + 1 ; act(na) = s
         endif
      enddo
      !$omp parallel do default(shared) schedule(dynamic, 16) private(k, s, t, hs)
      do k = 1, na                                         ! start: half the axial chord
         s = act(k)
         call ray_first_2d(self, pc(:, s), dir(:, s), s, t, hs)
         if (hs > 0) r(s) = min(r0, 0.5_R8P * t)
      enddo
      !$omp end parallel do
      do it = 1, MAXIT
         if (na == 0) exit
         allocate(c(2, na), dq(na), q(2, na))
         do k = 1, na
            s = act(k) ; c(:, k) = pc(:, s) + r(s) * dir(:, s)
         enddo
         call self%query_many(points=c, distances=dq, closest_points=q, is_signed=.false.)
         nk = 0
         do k = 1, na
            s = act(k)
            if (dq(k) >= r(s) * (1._R8P - 1e-12_R8P)) cycle
            w = q(:, k) - pc(:, s)
            t = 2._R8P * (dir(1, s) * w(1) + dir(2, s) * w(2))
            if (t <= 0._R8P) cycle
            rn = (w(1) * w(1) + w(2) * w(2)) / t
            if (rn >= r(s)) cycle
            w = q(:, k) - pc(:, s) - rn * dir(:, s)
            if (-(dir(1, s) * w(1) + dir(2, s) * w(2)) > cos_sep * rn) cycle
            r(s) = rn ; nk = nk + 1 ; act(nk) = s
         enddo
         na = nk
         deallocate(c, dq, q)
      enddo
      do k = 1, na
         bad(act(k)) = .true. ; r(act(k)) = huge(1._R8P)
      enddo
      do s = 1, self%nv
         if (self%next(s) == 0 .or. bad(s) .or. r(s) >= huge(1._R8P)) cycle
         if (side > 0._R8P .and. r(s) >= r0) cycle
         val(s) = 2._R8P * r(s)
      enddo
      endsubroutine medial
   endsubroutine curve_thickness

   subroutine ray_first_2d(self, o, d, excl, t, hs)
   !< First segment hit by the ray o + t d (t >= 0), the segment `excl` skipped: hs (0: none) and t. Segment tree,
   !< slab test on the node boxes, nearest first.
   type(curve_object), intent(in)  :: self   !< Curve.
   real(R8P),          intent(in)  :: o(2)   !< Origin.
   real(R8P),          intent(in)  :: d(2)   !< Direction (unit).
   integer(I4P),       intent(in)  :: excl   !< Segment to skip.
   real(R8P),          intent(out) :: t      !< Parameter of the hit.
   integer(I4P),       intent(out) :: hs     !< Hit segment (0: none).
   integer(I4P)                    :: stk(128), sp, nd, j, s !< Traversal.
   real(R8P)                       :: e(2), ao(2), den, tt, w, t0, t1, ta, tb !< Work.

   t = huge(1._R8P) ; hs = 0
   if (self%bvh%nodes_number == 0) return
   sp = 1 ; stk(1) = 1
   do while (sp > 0)
      nd = stk(sp) ; sp = sp - 1
      t0 = 0._R8P ; t1 = t                                  ! slab test against [0, t)
      if (d(1) /= 0._R8P) then
         ta = (self%bvh%lo_u(nd) - o(1)) / d(1) ; tb = (self%bvh%hi_u(nd) - o(1)) / d(1)
         t0 = max(t0, min(ta, tb)) ; t1 = min(t1, max(ta, tb))
      elseif (o(1) < self%bvh%lo_u(nd) .or. o(1) > self%bvh%hi_u(nd)) then
         cycle
      endif
      if (d(2) /= 0._R8P) then
         ta = (self%bvh%lo_v(nd) - o(2)) / d(2) ; tb = (self%bvh%hi_v(nd) - o(2)) / d(2)
         t0 = max(t0, min(ta, tb)) ; t1 = min(t1, max(ta, tb))
      elseif (o(2) < self%bvh%lo_v(nd) .or. o(2) > self%bvh%hi_v(nd)) then
         cycle
      endif
      if (t0 > t1) cycle
      if (self%bvh%left(nd) == 0) then
         do j = self%bvh%first(nd), self%bvh%last(nd)
            s = self%bvh%seg_id(j)
            if (s == excl) cycle
            e = [self%u(self%next(s)) - self%u(s), self%v(self%next(s)) - self%v(s)]
            ao = [self%u(s) - o(1), self%v(s) - o(2)]
            den = d(1) * e(2) - d(2) * e(1)
            if (den == 0._R8P) cycle
            tt = (ao(1) * e(2) - ao(2) * e(1)) / den
            w = (ao(1) * d(2) - ao(2) * d(1)) / den
            if (tt < 0._R8P .or. w < 0._R8P .or. w > 1._R8P .or. tt >= t) cycle
            t = tt ; hs = s
         enddo
      else
         sp = sp + 1 ; stk(sp) = self%bvh%right(nd)
         sp = sp + 1 ; stk(sp) = self%bvh%left(nd)
      endif
   enddo
   endsubroutine ray_first_2d

   subroutine curve_thin_regions(self, thickness, gap, h, factor, region, region_length, status)
   !< The segments thinner than factor * h (local size min(thickness, gap) below it, issue #62), grouped into runs of
   !< consecutive thin segments along each loop (a closed loop wraps around), sorted by decreasing length (ties: the
   !< lowest segment first): region(nv) is 0 or the region id, region_length(nr) the lengths.
   class(curve_object),       intent(in)            :: self             !< Curve.
   real(R8P),                 intent(in),  optional :: thickness(:)     !< Thickness per segment.
   real(R8P),                 intent(in),  optional :: gap(:)           !< Gap per segment.
   real(R8P),                 intent(in)            :: h                !< Cell size.
   real(R8P),                 intent(in)            :: factor           !< Threshold factor.
   integer(I4P), allocatable, intent(out)           :: region(:)        !< Region of each segment (0: not thin) [nv].
   real(R8P),    allocatable, intent(out)           :: region_length(:) !< Length of each region [nr].
   integer(I4P),              intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable                        :: ls(:), length(:) !< Local sizes, root lengths.
   integer(I4P), allocatable                        :: parent(:), roots(:), rank(:), root_id(:) !< Union-find, order.
   logical,      allocatable                        :: thin(:)          !< Thin segments.
   integer(I4P)                                     :: nv, s, a, b, nr, i, j, key !< Counters.
   logical                                          :: bad              !< Input check.

   nv = self%nv
   allocate(region(nv), region_length(0)) ; region = 0
   if (present(status)) status = STATUS_INVALID_INPUT
   bad = .not. (present(thickness) .or. present(gap)) .or. nv == 0 .or. h <= 0._R8P .or. factor <= 0._R8P
   if (.not. bad .and. present(thickness)) bad = size(thickness) /= nv
   if (.not. bad .and. present(gap)) bad = size(gap) /= nv
   if (bad) return
   allocate(ls(nv)) ; ls = huge(1._R8P)
   if (present(thickness)) ls = min(ls, thickness)
   if (present(gap)) ls = min(ls, gap)
   thin = ls < factor * h .and. self%next > 0
   allocate(parent(nv)) ; parent = [(s, s = 1, nv)]
   do s = 1, nv
      if (.not. thin(s)) cycle
      if (.not. thin(self%next(s))) cycle                 ! segment next(s) follows segment s along the loop
      a = find(s) ; b = find(self%next(s))
      if (a /= b) parent(max(a, b)) = min(a, b)
   enddo
   allocate(length(nv), root_id(nv), roots(nv)) ; length = 0._R8P ; root_id = 0
   do s = 1, nv
      if (thin(s)) then
         a = find(s) ; length(a) = length(a) + self%seg_len(s)
      endif
   enddo
   nr = 0
   do s = 1, nv
      if (thin(s) .and. find(s) == s) then
         nr = nr + 1 ; roots(nr) = s
      endif
   enddo
   allocate(rank(nr)) ; rank = roots(1:nr)
   do i = 2, nr
      key = rank(i) ; j = i - 1
      do while (j >= 1)
         if (length(rank(j)) >= length(key)) exit
         rank(j + 1) = rank(j) ; j = j - 1
      enddo
      rank(j + 1) = key
   enddo
   deallocate(region_length) ; allocate(region_length(nr))
   do i = 1, nr
      root_id(rank(i)) = i ; region_length(i) = length(rank(i))
   enddo
   do s = 1, nv
      if (thin(s)) region(s) = root_id(find(s))
   enddo
   if (present(status)) status = STATUS_OK
   contains
      function find(x) result(r)
      !< Root of x, with path halving.
      integer(I4P), intent(in) :: x !< Segment.
      integer(I4P)             :: r !< Root.

      r = x
      do while (parent(r) /= r)
         parent(r) = parent(parent(r)) ; r = parent(r)
      enddo
      endfunction find
   endsubroutine curve_thin_regions

   subroutine curve_build_thin(self, thickness, gap, status)
   !< Store min(thickness, gap) per segment tree payload entry, with its minimum per node (issue #62), for
   !< `box_thin_many(_device)` and `tag_block(thin_factor=)`. `adopt_loops` drops it; `move` keeps it. Call it before
   !< `enter_device` for the device twin.
   class(curve_object),    intent(inout)         :: self         !< Curve.
   real(R8P),              intent(in),  optional :: thickness(:) !< Thickness per segment.
   real(R8P),              intent(in),  optional :: gap(:)       !< Gap per segment.
   integer(I4P),           intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable                     :: ls(:)        !< Local sizes.
   integer(I4P), allocatable                     :: order(:), stk(:) !< Pre-order, stack.
   integer(I4P)                                  :: nn, no, sp, nd, j, k !< Counters.
   logical                                       :: bad          !< Input check.

   call curve_drop_thin(self)
   if (present(status)) status = STATUS_INVALID_INPUT
   bad = .not. (present(thickness) .or. present(gap)) .or. self%bvh%nodes_number == 0
   if (.not. bad .and. present(thickness)) bad = size(thickness) /= self%nv
   if (.not. bad .and. present(gap)) bad = size(gap) /= self%nv
   if (bad) return
   allocate(ls(self%nv)) ; ls = huge(1._R8P)
   if (present(thickness)) ls = min(ls, thickness)
   if (present(gap)) ls = min(ls, gap)
   allocate(self%thin_pay(self%bvh%segments_number))
   do j = 1, self%bvh%segments_number
      self%thin_pay(j) = ls(self%bvh%seg_id(j))
   enddo
   nn = self%bvh%nodes_number
   allocate(self%thin_node(nn), order(nn), stk(nn)) ; self%thin_node = huge(1._R8P)
   no = 0 ; sp = 1 ; stk(1) = 1
   do while (sp > 0)
      nd = stk(sp) ; sp = sp - 1
      no = no + 1 ; order(no) = nd
      if (self%bvh%left(nd) > 0) then
         sp = sp + 1 ; stk(sp) = self%bvh%left(nd) ; sp = sp + 1 ; stk(sp) = self%bvh%right(nd)
      endif
   enddo
   do k = no, 1, -1                                       ! reversed pre-order: children before parents
      nd = order(k)
      if (self%bvh%left(nd) == 0) then
         do j = self%bvh%first(nd), self%bvh%last(nd)
            self%thin_node(nd) = min(self%thin_node(nd), self%thin_pay(j))
         enddo
      else
         self%thin_node(nd) = min(self%thin_node(self%bvh%left(nd)), self%thin_node(self%bvh%right(nd)))
      endif
   enddo
   self%thin_built = .true.
   if (present(status)) status = STATUS_OK
   endsubroutine curve_build_thin

   subroutine curve_drop_thin(self)
   !< Forget the stored local sizes. The device copy, if any, is released by exit_device.
   type(curve_object), intent(inout) :: self !< Curve.

   if (self%thin_on_device) return
   if (allocated(self%thin_pay)) deallocate(self%thin_pay)
   if (allocated(self%thin_node)) deallocate(self%thin_node)
   self%thin_built = .false.
   endsubroutine curve_drop_thin

   subroutine curve_box_thin_many(self, lo, hi, factor, touched, status)
   !< Which boxes [lo, hi] (u, v) meet (exactly) a segment thinner than factor * h, h the larger side of each box
   !< (issue #62; `build_thin` first, else STATUS_INVALID_INPUT). OpenMP over the boxes; the device twin runs the same
   !< routine (`curve_box_thin_oac`).
   class(curve_object), intent(in)            :: self       !< Curve.
   real(R8P),           intent(in)            :: lo(:,:)    !< Box minima [2, n].
   real(R8P),           intent(in)            :: hi(:,:)    !< Box maxima [2, n].
   real(R8P),           intent(in)            :: factor     !< Threshold factor.
   logical,             intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),        intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                               :: i          !< Counter.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. self%thin_built .or. size(lo, 1) /= 2 .or. size(hi, 1) /= 2 .or. size(lo, 2) /= size(hi, 2) .or. &
       size(touched) /= size(lo, 2)) return
   !$omp parallel do default(shared) schedule(dynamic, 64) private(i)
   do i = 1, size(lo, 2)
      touched(i) = curve_box_thin_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, &
                                      self%bvh%right, self%bvh%first, self%bvh%last, self%bvh%seg_id, self%next, self%u, &
                                      self%v, self%thin_pay, self%thin_node, lo(:, i), hi(:, i), factor)
   enddo
   !$omp end parallel do
   if (present(status)) status = STATUS_OK
   endsubroutine curve_box_thin_many

   subroutine curve_box_thin_many_device(self, lo, hi, factor, touched, status)
   !< `box_thin_many` on the device (issue #62): requires `build_thin` before `enter_device` (else STATUS_INVALID_INPUT).
   !< #20 rules: top-level pointer aliases of the mapped arrays.
   class(curve_object), intent(inout), target :: self       !< Curve (target: device aliases).
   real(R8P),           intent(in)            :: lo(:,:)    !< Box minima [2, n].
   real(R8P),           intent(in)            :: hi(:,:)    !< Box maxima [2, n].
   real(R8P),           intent(in)            :: factor     !< Threshold factor.
   logical,             intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),        intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), u(:), v(:), tp(:), tn(:) !< Aliases.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:)  !< Aliases.
   logical, allocatable  :: tl(:)          !< Device-side result.
   real(R8P)             :: blo(2), bhi(2) !< Per-box work.
   integer(I4P)          :: i, n           !< Counters.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. (self%is_on_device .and. self%thin_on_device) .or. size(lo, 1) /= 2 .or. size(hi, 1) /= 2 .or. &
       size(lo, 2) /= size(hi, 2) .or. size(touched) /= size(lo, 2)) then
      write(stderr, '(A)') 'fossil_curve_object%box_thin_many_device: build_thin, then enter_device, first'
      return
   endif
   n = size(lo, 2)
   allocate(tl(n))
   lo_u => self%bvh%lo_u ; lo_v => self%bvh%lo_v ; hi_u => self%bvh%hi_u ; hi_v => self%bvh%hi_v
   left => self%bvh%left ; right => self%bvh%right ; first => self%bvh%first ; last => self%bvh%last
   seg_id => self%bvh%seg_id ; next => self%next ; u => self%u ; v => self%v ; tp => self%thin_pay ; tn => self%thin_node
   !$acc parallel loop present(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, tp, tn) &
   !$acc&              copyin(lo, hi) copyout(tl) private(blo, bhi)
   do i = 1, n
      blo(1) = lo(1, i) ; blo(2) = lo(2, i) ; bhi(1) = hi(1, i) ; bhi(2) = hi(2, i)
      tl(i) = curve_box_thin_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, tp, tn, blo, bhi, &
                                 factor)
   enddo
   touched = tl
   if (present(status)) status = STATUS_OK
   endsubroutine curve_box_thin_many_device

   subroutine line_crossings(self, axis, a, t0, t1, t, segments, entering, n, status)
   !< Exact crossings of the line along `axis` (1: along u at v = a; 2: along v at u = a) with the
   !< curve (issue #40): every crossing with parameter in [t0, t1], sorted, with its segment and
   !< whether the line (in +axis) enters the solid. Exact side test with a symbolic tie-break
   !< (`fossil_line_query`), so on closed loops every crossing is counted once; the parameters are
   !< rounded. `status = STATUS_INVALID_INPUT` if axis is not 1 or 2, or t1 < t0.
   class(curve_object),       intent(in)            :: self        !< Curve.
   integer(I4P),              intent(in)            :: axis        !< 1 (u) or 2 (v).
   real(R8P),                 intent(in)            :: a           !< Line position.
   real(R8P),                 intent(in)            :: t0, t1      !< Parameter range.
   real(R8P),    allocatable, intent(out)           :: t(:)        !< Crossing parameters [n].
   integer(I4P), allocatable, intent(out)           :: segments(:) !< Segments [n].
   logical,      allocatable, intent(out)           :: entering(:) !< Entering flags [n].
   integer(I4P),              intent(out)           :: n           !< Number of crossings.
   integer(I4P),              intent(out), optional :: status      !< STATUS_OK or STATUS_INVALID_INPUT.

   if (present(status)) status = STATUS_OK
   n = 0_I4P
   if (axis < 1_I4P .or. axis > 2_I4P .or. t1 < t0) then
      if (present(status)) status = STATUS_INVALID_INPUT
      allocate(t(0), segments(0), entering(0))
      return
   endif
   call curve_crossings_core(self, axis, a, t0, t1, t, segments, entering, n)
   endsubroutine line_crossings

   subroutine cell_wall_fractions(self, x, y, frac, inside, bad_lines, status)
   !< Wall fractions along the grid lines through the cell centres x(nx), y(ny) (u, v) of a block
   !< (issue #40): `frac(d, i, j)` for d = 1..4 (-u, +u, -v, +v) as in the surface twin (-1: no
   !< crossing before the neighbouring centre, or a block end); `inside(i, j)` by the parity along the
   !< u line from outside the curve box; `bad_lines` counts lines whose flags do not alternate.
   class(curve_object), intent(in)            :: self          !< Curve.
   real(R8P),           intent(in)            :: x(:), y(:)    !< Cell centres (u, v).
   real(R8P),           intent(out)           :: frac(:,:,:)   !< Fractions [4, nx, ny].
   logical,             intent(out)           :: inside(:,:)   !< Centres inside [nx, ny].
   integer(I4P),        intent(out), optional :: bad_lines     !< Lines with non-alternating flags.
   integer(I4P),        intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: nx, ny, i, j, nbad !< Work.
   logical      :: bad                !< Size check.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) ; ny = size(y, kind=I4P)
   bad = nx < 1 .or. ny < 1
   if (.not. bad) bad = any(x(2:) <= x(:nx-1)) .or. any(y(2:) <= y(:ny-1)) .or. size(frac, 1) /= 4 .or. &
                        size(frac, 2) /= nx .or. size(frac, 3) /= ny .or. size(inside, 1) /= nx .or. size(inside, 2) /= ny
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   frac = -1._R8P ; inside = .false. ; nbad = 0_I4P
   if (self%bvh%nodes_number == 0_I4P) then
      if (present(bad_lines)) bad_lines = 0_I4P
      return
   endif
   !$omp parallel default(shared) private(i, j) reduction(+:nbad)
   !$omp do schedule(dynamic)
   do j = 1_I4P, ny
      call sweep(1_I4P, y(j), x, self%bmin(1), frac(1:2, :, j), nbad, inside(:, j))
   enddo
   !$omp end do nowait
   !$omp do schedule(dynamic)
   do i = 1_I4P, nx
      call sweep(2_I4P, x(i), y, self%bmin(2), frac(3:4, i, :), nbad)
   enddo
   !$omp end do
   !$omp end parallel
   if (present(bad_lines)) bad_lines = nbad
   contains
      subroutine sweep(ax, a, c, cmin, fr, nb, ins)
      !< One grid line: fractions toward both neighbours, parity.
      integer(I4P), intent(in)              :: ax      !< Axis.
      real(R8P),    intent(in)              :: a       !< Line position.
      real(R8P),    intent(in)              :: c(:)    !< Centres along the line.
      real(R8P),    intent(in)              :: cmin    !< Curve box minimum along the axis.
      real(R8P),    intent(inout)           :: fr(:,:) !< Fractions [2, centres].
      integer(I4P), intent(inout)           :: nb      !< Bad lines.
      logical,      intent(inout), optional :: ins(:)  !< Parity inside flags.
      real(R8P),    allocatable :: t(:)   !< Crossings.
      integer(I4P), allocatable :: g(:)   !< Segments.
      logical,      allocatable :: e(:)   !< Entering flags.
      integer(I4P)              :: n, m, q, nc !< Work.

      nc = size(c, kind=I4P)
      call curve_crossings_core(self, ax, a, min(c(1), cmin) - 1._R8P, c(nc), t, g, e, n)
      do q = 1_I4P, n
         if (e(q) .neqv. (mod(q, 2_I4P) == 1_I4P)) then
            nb = nb + 1_I4P
            exit
         endif
      enddo
      m = 0_I4P
      do q = 1_I4P, nc
         do while (m < n)
            if (t(m + 1) >= c(q)) exit
            m = m + 1_I4P
         enddo
         if (present(ins)) ins(q) = mod(m, 2_I4P) == 1_I4P
         if (q > 1_I4P .and. m > 0_I4P) then
            if (t(m) >= c(q - 1)) fr(1, q) = (c(q) - t(m)) / (c(q) - c(q - 1))
         endif
         if (q < nc .and. m < n) then
            if (t(m + 1) <= c(q + 1)) fr(2, q) = (t(m + 1) - c(q)) / (c(q + 1) - c(q))
         endif
      enddo
      endsubroutine sweep
   endsubroutine cell_wall_fractions

   subroutine cut_cells(self, x, y, cut, status)
   !< Exact cut cells of the block with cell edges x(nx+1), y(ny+1) (u, v) (issue #40): fluid area
   !< fraction, edge apertures, fluid centroid, wetted length, centroid and normal integral, and the
   !< fluid pieces of every cell (`fossil_cut_cells_2d`; the fluid is outside the solid). Exact
   !< topology with the symbolic tie-break of the line crossings; the measures are exact to round-off.
   !< Cells whose boundary does not close (open loops) are flagged in `cut%closed` with no measures.
   !< `status = STATUS_INVALID_INPUT` on bad edges.
   class(curve_object),  intent(in)            :: self   !< Curve.
   real(R8P),            intent(in)            :: x(:), y(:) !< Cell edges (u, v), strictly increasing.
   type(cut_cells_2d_t), intent(out)           :: cut    !< Cut cells.
   integer(I4P),         intent(out), optional :: status !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)                                :: st     !< Module status.

   if (present(status)) status = STATUS_OK
   if (self%nv > 0_I4P) then
      call compute_cut_cells_2d(self%u, self%v, self%next, x, y, cut, st)
   else
      block
         real(R8P)    :: uu(0), vv(0)
         integer(I4P) :: nn(0)
         call compute_cut_cells_2d(uu, vv, nn, x, y, cut, st)
      end block
   endif
   if (st /= 0_I4P .and. present(status)) status = STATUS_INVALID_INPUT
   endsubroutine cut_cells

   subroutine line_crossings_many(self, axis, a, t0, t1, t, segments, entering, counts, overflow, status)
   !< Crossings of the lines along `axis` (1: u, 2: v) at a(i), i = 1..n, with parameter in [t0, t1]
   !< (issue #40), into fixed buffers t(nmax, n), segments(nmax, n), entering(nmax, n), sorted;
   !< counts(i), overflow(i) as in the surface twin (on overflow, re-run with a larger buffer). OpenMP over the lines; the kernel of the device
   !< twin. `status = STATUS_INVALID_INPUT` on a bad axis, range or shape, or an empty curve.
   class(curve_object), intent(in)            :: self          !< Curve.
   integer(I4P),        intent(in)            :: axis          !< 1 or 2.
   real(R8P),           intent(in)            :: a(:)          !< Line positions [n].
   real(R8P),           intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),           intent(out)           :: t(:,:)        !< Parameters [nmax, n].
   integer(I4P),        intent(out)           :: segments(:,:) !< Segments [nmax, n].
   logical,             intent(out)           :: entering(:,:) !< Entering flags [nmax, n].
   integer(I4P),        intent(out)           :: counts(:)     !< Crossings per line [n].
   logical,             intent(out)           :: overflow(:)   !< Buffer overflow per line [n].
   integer(I4P),        intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: i !< Counter.

   if (.not. curve_many_ok(self, axis, a, t0, t1, t, segments, entering, counts, overflow, status)) return
   !$omp parallel do default(shared) schedule(dynamic)
   do i = 1_I4P, size(a, kind=I4P)
      call curve_line_crossings_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, &
                                    self%bvh%right, self%bvh%first, self%bvh%last, self%bvh%seg_id, self%next, self%u, &
                                    self%v, axis, a(i), t0, t1, t(:, i), segments(:, i), entering(:, i), counts(i), &
                                    overflow(i))
   enddo
   !$omp end parallel do
   endsubroutine line_crossings_many

   subroutine line_crossings_many_device(self, axis, a, t0, t1, t, segments, entering, counts, overflow, status)
   !< `line_crossings_many` on the device (OpenACC, issue #40): same results bit for bit. Requires
   !< `enter_device`. #20 rules (pointer aliases). Validated on the host; the GPU gate is pending.
   class(curve_object), intent(inout), target :: self          !< Curve (target: device aliases).
   integer(I4P),        intent(in)            :: axis          !< 1 or 2.
   real(R8P),           intent(in)            :: a(:)          !< Line positions [n].
   real(R8P),           intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),           intent(out)           :: t(:,:)        !< Parameters [nmax, n].
   integer(I4P),        intent(out)           :: segments(:,:) !< Segments [nmax, n].
   logical,             intent(out)           :: entering(:,:) !< Entering flags [nmax, n].
   integer(I4P),        intent(out)           :: counts(:)     !< Crossings per line [n].
   logical,             intent(out)           :: overflow(:)   !< Buffer overflow per line [n].
   integer(I4P),        intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), u(:), v(:) !< Aliases.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:) !< Aliases.
   real(R8P),    allocatable :: tl(:,:) !< Device-side parameters.
   integer(I4P), allocatable :: sl(:,:) !< Device-side segments.
   logical,      allocatable :: el(:,:) !< Device-side flags.
   integer(I4P), allocatable :: cl(:)   !< Device-side counts.
   logical,      allocatable :: ol(:)   !< Device-side overflow.
   integer(I4P)              :: i, n, nmax !< Work.

   if (.not. curve_many_ok(self, axis, a, t0, t1, t, segments, entering, counts, overflow, status)) return
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_curve_object%line_crossings_many_device: curve not on device -- call enter_device first.'
      return
   endif
   n = size(a, kind=I4P) ; nmax = size(t, 1, kind=I4P)
   if (n == 0_I4P) return
   allocate(tl(nmax, n), sl(nmax, n), el(nmax, n), cl(n), ol(n))
   tl = 0._R8P ; sl = 0_I4P ; el = .false.
   lo_u => self%bvh%lo_u ; lo_v => self%bvh%lo_v ; hi_u => self%bvh%hi_u ; hi_v => self%bvh%hi_v
   left => self%bvh%left ; right => self%bvh%right ; first => self%bvh%first ; last => self%bvh%last
   seg_id => self%bvh%seg_id ; next => self%next ; u => self%u ; v => self%v
   !$acc parallel loop present(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v) copyin(a) &
   !$acc&              copy(tl, sl, el) copyout(cl, ol)
   do i = 1_I4P, n
      call curve_line_crossings_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, axis, a(i), &
                                    t0, t1, tl(:, i), sl(:, i), el(:, i), cl(i), ol(i))
   enddo
   t = tl ; segments = sl ; entering = el ; counts = cl ; overflow = ol
   endsubroutine line_crossings_many_device

   recursive subroutine distance_many(self, points, distances, is_signed, sign_algorithm, status, pose)
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
   type(rigid_pose_2d_t), intent(in), optional :: pose          !< Pose of the curve (issue #41): world-frame points.
   logical                                    :: signed_, fwn   !< Resolved options.
   integer(I4P)                               :: i, n           !< Work.

   if (present(pose)) then ! body-frame query (issue #41): distances are invariant
      call self%distance_many(points=body_points2(pose, points), distances=distances, is_signed=is_signed, &
                              sign_algorithm=sign_algorithm, status=status)
      return
   endif
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

   recursive subroutine distance_many_device(self, points, distances, is_signed, sign_algorithm, status, pose)
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
   type(rigid_pose_2d_t), intent(in), optional :: pose          !< Pose of the curve (issue #41): world-frame points.
   logical                                    :: signed_, fwn   !< Resolved options.
   integer(I4P)                               :: i, n           !< Work.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Aliases.
   real(R8P),    pointer :: u(:), v(:), pn_u(:), pn_v(:), wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_diag2(:) !< Aliases.
   real(R8P),    pointer :: wn_c(:,:) !< Alias.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:) !< Aliases.

   if (present(pose)) then ! body-frame query (issue #41): the transform of the CPU path, on the host
      call self%distance_many_device(points=body_points2(pose, points), distances=distances, is_signed=is_signed, &
                                     sign_algorithm=sign_algorithm, status=status)
      return
   endif
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

   recursive subroutine query_many(self, points, distances, closest_points, segments, vertices, normals, bodies, is_signed, &
                                   sign_algorithm, ambiguous, second_foot_points, rel_tolerance, foot_tolerance, status, pose, &
                                   curvatures, feature_distances, feature_angle, patches)
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
   !< Curvature at the foot point (issue #42 step 3, CPU only): `curvatures(i)` is the vertex curvature
   !< (turning angle over the mean length of the two segments) interpolated linearly along the closest
   !< segment. Positive where the curve is convex seen from outside the solid (an outer loop turning
   !< left), negative at concave walls and on holes; 0 at the ends of an open loop. A regular N-gon of
   !< circumradius R gives (1/R) (pi/N) / sin(pi/N) = (1/R) (1 + O(h^2)).
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
   type(rigid_pose_2d_t), intent(in), optional :: pose              !< Pose of the curve (issue #41): world points, outputs.
   real(R8P),           intent(out), optional :: feature_distances(:) !< Distance to the nearest corner [n] (issue #61).
   real(R8P),           intent(in),  optional :: feature_angle      !< Corner threshold [degrees] (default 30).
   real(R8P),           intent(out), optional :: curvatures(:)      !< Curvature at the foot point [n].
   integer(I4P),        intent(out), optional :: patches(:)         !< Patch of the closest segment [n] (issue #66).
   real(R8P)                                  :: cp(2), nn(2), cp2(2), rel_, foot_ !< Per-point outputs, tolerances.
   logical                                    :: amb, want_amb      !< Ambiguity.
   integer(I4P)                               :: sg, vx, b, i, n    !< Per-point outputs, counters.
   logical                                    :: signed_, fwn       !< Resolved options.

   if (present(pose)) then ! body-frame query (issue #41), outputs mapped back to the world frame
      call self%query_many(points=body_points2(pose, points), distances=distances, closest_points=closest_points, &
                           segments=segments, vertices=vertices, normals=normals, bodies=bodies, is_signed=is_signed, &
                           patches=patches, &
                           sign_algorithm=sign_algorithm, ambiguous=ambiguous, second_foot_points=second_foot_points, &
                           rel_tolerance=rel_tolerance, foot_tolerance=foot_tolerance, status=status, curvatures=curvatures, &
                           feature_distances=feature_distances, feature_angle=feature_angle)
      if (present(status)) then
         if (status /= STATUS_OK) return
      endif
      do i = 1_I4P, size(points, 2, kind=I4P)
         if (present(closest_points)) closest_points(:, i) = to_world(pose, closest_points(:, i))
         if (present(second_foot_points)) second_foot_points(:, i) = to_world(pose, second_foot_points(:, i))
         if (present(normals)) normals(:, i) = rotate(pose, normals(:, i))
      enddo
      return
   endif
   if (present(feature_distances)) then ! distance to the nearest corner (issue #61)
      if (size(feature_distances) /= size(points, 2)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      call corner_distances(self, points, feature_angle, feature_distances)
   endif
   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   if (.not. rich_sizes_ok(size(distances, kind=I4P), closest_points, segments, vertices, normals, bodies) .or. &
       .not. rich_sizes_ok(size(distances, kind=I4P), closest_points=second_foot_points)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (present(patches)) then
      if (size(patches) /= size(distances)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   if (present(ambiguous)) then
      if (size(ambiguous) /= size(distances)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
   endif
   if (present(curvatures)) then
      if (size(curvatures) /= size(distances)) then
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
      if (present(patches)) patches(i) = self%segment_patch(sg)
      if (present(vertices)) vertices(i) = vx
      if (present(normals)) normals(:, i) = nn
      if (present(bodies)) bodies(i) = b
      if (want_amb) then
         call second_foot_2d(self, points(:, i), cp, abs(distances(i)), rel_, foot_, amb, cp2)
         if (present(ambiguous)) ambiguous(i) = amb
         if (present(second_foot_points)) second_foot_points(:, i) = cp2
      endif
      if (present(curvatures)) curvatures(i) = foot_curvature_2d(self, sg, cp)
   enddo
   !$omp end parallel do
   endsubroutine query_many

   pure function vertex_curvature(self, i) result(kappa)
   !< Discrete curvature at vertex i (issue #42 step 3): the turning angle between the incoming and
   !< outgoing segments over the mean of their lengths; positive for a left turn (convex seen from
   !< outside the solid). 0 at the ends of an open loop and at a degenerate vertex.
   class(curve_object), intent(in) :: self  !< Curve.
   integer(I4P),        intent(in) :: i     !< Vertex.
   real(R8P)                       :: kappa !< Curvature.
   real(R8P)                       :: au, av, bu, bv, la, lb !< Incoming, outgoing segments.
   integer(I4P)                    :: ip, ix !< Neighbours.

   kappa = 0._R8P
   ip = self%prev(i) ; ix = self%next(i)
   if (ip <= 0_I4P .or. ix <= 0_I4P) return
   au = self%u(i) - self%u(ip) ; av = self%v(i) - self%v(ip)
   bu = self%u(ix) - self%u(i) ; bv = self%v(ix) - self%v(i)
   la = hypot(au, av) ; lb = hypot(bu, bv)
   if (la <= 0._R8P .or. lb <= 0._R8P) return
   kappa = atan2(au * bv - av * bu, au * bu + av * bv) / (0.5_R8P * (la + lb))
   endfunction vertex_curvature

   pure function foot_curvature_2d(self, sg, cp) result(kappa)
   !< Curvature at a foot point on segment sg (issue #42 step 3): the two vertex curvatures
   !< interpolated linearly in the foot's parameter along the segment (clamped to [0, 1]).
   class(curve_object), intent(in) :: self  !< Curve.
   integer(I4P),        intent(in) :: sg    !< Segment (its start vertex; 0: none).
   real(R8P),           intent(in) :: cp(2) !< Foot point.
   real(R8P)                       :: kappa !< Curvature.
   real(R8P)                       :: du, dv, l2, t !< Segment, parameter.
   integer(I4P)                    :: j     !< End vertex.

   kappa = 0._R8P
   if (sg <= 0_I4P) return
   j = self%next(sg)
   if (j <= 0_I4P) return
   du = self%u(j) - self%u(sg) ; dv = self%v(j) - self%v(sg)
   l2 = du * du + dv * dv
   t = 0._R8P
   if (l2 > 0._R8P) t = min(max(((cp(1) - self%u(sg)) * du + (cp(2) - self%v(sg)) * dv) / l2, 0._R8P), 1._R8P)
   kappa = (1._R8P - t) * vertex_curvature(self, sg) + t * vertex_curvature(self, j)
   endfunction foot_curvature_2d

   subroutine curve_markers(self, dx, centroids, lengths, normals, spacing_fraction, spacings, status)
   !< Grid-adapted curve markers (issue #42 step 4), the 2-D twin of `surface%surface_markers`.
   !<
   !< Each loop of length L is split exactly along the polyline into n = max(nint(L / dx), 1) pieces of
   !< arc length h = L / n. Marker k sits at arc length (k - 1/2) h on the polyline (no chord error),
   !< with length h and the outward normal of the segment it lies on. `spacings` is the chord distance
   !< from each marker to the next one of its loop (to the first for a closed loop; 0 for the last
   !< marker of an open loop), and `spacing_fraction` the fraction of those (excluding the zeros) in
   !< [0.8, 1.2] dx. `status = STATUS_INVALID_INPUT` for dx <= 0 or an empty curve (empty outputs).
   class(curve_object),    intent(in)            :: self             !< Curve.
   real(R8P),              intent(in)            :: dx               !< Grid spacing.
   real(R8P), allocatable, intent(out)           :: centroids(:,:)   !< Marker positions [2, m].
   real(R8P), allocatable, intent(out)           :: lengths(:)       !< Marker lengths [m].
   real(R8P), allocatable, intent(out)           :: normals(:,:)     !< Marker unit normals [2, m].
   real(R8P),              intent(out), optional :: spacing_fraction !< Fraction of spacings in [0.8, 1.2] dx.
   real(R8P), allocatable, intent(out), optional :: spacings(:)      !< Chord spacings [m].
   integer(I4P),           intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P), allocatable                     :: nm(:)            !< Markers per loop.
   real(R8P), allocatable                        :: ll(:), sp(:)     !< Loop lengths, spacings.
   real(R8P)                                     :: h, sk, s0, t     !< Piece length, marker arc length, segment start, parameter.
   integer(I4P)                                  :: l, i, j, k, m, m0, ns, nsp !< Counters.

   if (present(status)) status = STATUS_OK
   allocate(centroids(2, 0), lengths(0), normals(2, 0))
   if (present(spacings)) allocate(spacings(0))
   if (present(spacing_fraction)) spacing_fraction = 0._R8P
   if (.not. (dx > 0._R8P) .or. self%nl == 0_I4P) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(nm(self%nl), ll(self%nl))
   do l = 1_I4P, self%nl
      ll(l) = 0._R8P
      do i = self%loop_start(l), self%loop_start(l + 1) - 1_I4P
         if (self%next(i) > 0_I4P) ll(l) = ll(l) + self%seg_len(i)
      enddo
      nm(l) = 0_I4P ; if (ll(l) > 0._R8P) nm(l) = max(nint(ll(l) / dx, I4P), 1_I4P)
   enddo
   m = sum(nm)
   deallocate(centroids, lengths, normals)
   allocate(centroids(2, m), lengths(m), normals(2, m), sp(m))
   m0 = 0_I4P ; ns = 0_I4P ; nsp = 0_I4P
   do l = 1_I4P, self%nl
      if (nm(l) == 0_I4P) cycle
      h = ll(l) / nm(l)
      i = self%loop_start(l) ; s0 = 0._R8P
      do k = 1_I4P, nm(l)
         sk = (k - 0.5_R8P) * h
         ! advance to the segment holding arc length sk (the last segment takes any round-off excess)
         do while (s0 + self%seg_len(i) < sk .and. self%next(self%next(i)) > 0_I4P .and. self%next(i) /= self%loop_start(l))
            s0 = s0 + self%seg_len(i) ; i = self%next(i)
         enddo
         j = self%next(i)
         t = 0._R8P ; if (self%seg_len(i) > 0._R8P) t = min(max((sk - s0) / self%seg_len(i), 0._R8P), 1._R8P)
         centroids(:, m0 + k) = [self%u(i) + t * (self%u(j) - self%u(i)), self%v(i) + t * (self%v(j) - self%v(i))]
         normals(:, m0 + k) = [self%seg_nu(i), self%seg_nv(i)]
         lengths(m0 + k) = h
      enddo
      do k = 1_I4P, nm(l)
         j = k + 1_I4P
         if (j > nm(l)) then
            if (.not. self%loop_closed(l) .or. nm(l) == 1_I4P) then
               sp(m0 + k) = 0._R8P
               cycle
            endif
            j = 1_I4P
         endif
         sp(m0 + k) = hypot(centroids(1, m0 + j) - centroids(1, m0 + k), centroids(2, m0 + j) - centroids(2, m0 + k))
         nsp = nsp + 1_I4P
         if (sp(m0 + k) >= 0.8_R8P * dx .and. sp(m0 + k) <= 1.2_R8P * dx) ns = ns + 1_I4P
      enddo
      m0 = m0 + nm(l)
   enddo
   if (present(spacing_fraction) .and. nsp > 0_I4P) spacing_fraction = ns / real(nsp, R8P)
   if (present(spacings)) call move_alloc(from=sp, to=spacings)
   endsubroutine curve_markers

   recursive subroutine query_many_device(self, points, distances, closest_points, segments, vertices, normals, bodies, &
                                          is_signed, sign_algorithm, status, pose, feature_distances, feature_angle, patches)
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
   type(rigid_pose_2d_t), intent(in), optional :: pose              !< Pose of the curve (issue #41).
   real(R8P),           intent(out), optional :: feature_distances(:) !< Distance to the nearest corner [n] (issue #61).
   real(R8P),           intent(in),  optional :: feature_angle      !< Corner threshold [degrees] (default 30).
   integer(I4P),        intent(out), optional :: patches(:)         !< Patch of the closest segment [n] (issue #66), host lookup.
   real(R8P),    allocatable :: dd(:), cpl(:,:), nml(:,:) !< Device-side outputs.
   integer(I4P), allocatable :: sgl(:), vxl(:), bdl(:)    !< Device-side outputs.
   logical                   :: signed_, fwn              !< Resolved options.
   integer(I4P)              :: i, n                      !< Work.
   real(R8P),    pointer :: lo_u(:), lo_v(:), hi_u(:), hi_v(:), a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Aliases.
   real(R8P),    pointer :: u(:), v(:), pn_u(:), pn_v(:), seg_nu(:), seg_nv(:) !< Aliases.
   real(R8P),    pointer :: wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:), wn_diag2(:), wn_c(:,:) !< Aliases.
   integer(I4P), pointer :: left(:), right(:), first(:), last(:), seg_id(:), next(:), vbody(:) !< Aliases.

   if (present(pose)) then ! body-frame query (issue #41), outputs mapped back on the host
      call self%query_many_device(points=body_points2(pose, points), distances=distances, closest_points=closest_points, &
                                  segments=segments, vertices=vertices, normals=normals, bodies=bodies, is_signed=is_signed, &
                                  sign_algorithm=sign_algorithm, status=status, feature_distances=feature_distances, &
                                  feature_angle=feature_angle, patches=patches)
      if (present(status)) then
         if (status /= STATUS_OK) return
      endif
      do i = 1_I4P, size(points, 2, kind=I4P)
         if (present(closest_points)) closest_points(:, i) = to_world(pose, closest_points(:, i))
         if (present(normals)) normals(:, i) = rotate(pose, normals(:, i))
      enddo
      return
   endif
   if (present(feature_distances)) then ! distance to the nearest corner (issue #61), on the device
      if (size(feature_distances) /= size(points, 2)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      call corner_distances(self, points, feature_angle, feature_distances, device=.true.)
   endif
   if (.not. resolve_batch_options(self, points, distances, is_signed, sign_algorithm, status, signed_, fwn)) return
   if (.not. rich_sizes_ok(size(distances, kind=I4P), closest_points, segments, vertices, normals, bodies, patches)) then
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
   if (present(patches)) patches = self%segment_patch(sgl)
   endsubroutine query_many_device

   subroutine probe_many(self, points, lengths, probe_points, visible, foot_points, normals, status, patches)
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
   integer(I4P),        intent(out), optional :: patches(:)        !< Patch of the closest segment [n] (issue #66).
   real(R8P)                                  :: d, cp(2), nn(2), tip(2) !< Per-point work.
   integer(I4P)                               :: i, n, sg, vx      !< Work.

   if (present(status)) status = STATUS_OK
   n = size(points, 2, kind=I4P)
   if (size(points, 1) /= 2 .or. size(lengths) /= n .or. size(probe_points, 1) /= 2 .or. size(probe_points, 2) /= n .or. &
       size(visible) /= n .or. .not. rich_sizes_ok(n, closest_points=foot_points, normals=normals, patches=patches)) then
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
      if (present(patches)) patches(i) = self%segment_patch(sg)
   enddo
   !$omp end parallel do
   endsubroutine probe_many

   subroutine image_points(self, points, images, visible, foot_points, status, patches)
   !< Image points (issue #37): the mirror of each point through its closest point, 2 foot - x, and
   !< whether the segment foot-image crosses no other segment (exact, see `probe_many`).
   class(curve_object), intent(in)            :: self             !< Curve.
   real(R8P),           intent(in)            :: points(:,:)      !< Points [2, n].
   real(R8P),           intent(out)           :: images(:,:)      !< Image points [2, n].
   logical,             intent(out), optional :: visible(:)       !< Segment foot-image crosses nothing [n].
   real(R8P),           intent(out), optional :: foot_points(:,:) !< Foot points [2, n].
   integer(I4P),        intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),        intent(out), optional :: patches(:)       !< Patch of the closest segment [n] (issue #66).
   real(R8P)                                  :: d, cp(2), img(2) !< Per-point work.
   integer(I4P)                               :: i, n, sg, vx     !< Work.
   logical                                    :: bad              !< Size check.

   if (present(status)) status = STATUS_OK
   n = size(points, 2, kind=I4P)
   bad = size(points, 1) /= 2 .or. size(images, 1) /= 2 .or. size(images, 2) /= n .or. &
         .not. rich_sizes_ok(n, closest_points=foot_points, patches=patches)
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
      if (present(patches)) patches(i) = self%segment_patch(sg)
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
   if (self%thin_built) then                       ! the local sizes, if built (issue #62)
      !$acc enter data copyin(self%thin_pay, self%thin_node)
      self%thin_on_device = .true.
   endif
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
   if (self%thin_on_device) then
      !$acc exit data delete(self%thin_pay, self%thin_node)
      self%thin_on_device = .false.
   endif
   self%is_on_device = .false.
   endsubroutine exit_device

   function check_validity(self, loop, check_self_intersections) result(defect)
   !< First validity defect of the curve (CURVE_VALID if none), and the loop it was found in.
   !<
   !< Checked loop by loop: closed; at least 3 vertices; no zero-length segment (exactly coincident
   !< consecutive vertices); nonzero signed area; orientation matching the nesting depth (depth
   !< even = outer or island: counter-clockwise, positive area; odd = hole: clockwise). The
   !< nesting depth is exact (`nesting_depth`). Self-intersections are checked only with
   !< `check_self_intersections` (default .false., the meaning of CURVE_VALID before issue #44): after
   !< every loop passed, CURVE_SELF_INTERSECTING if `self_intersections` finds a pair, `loop` the loop of
   !< its first segment. A curve that crosses itself breaks the cut cells and the swept geometry (the
   !< regions are counted by winding number there, by parity in the line crossings).
   class(curve_object), intent(in)            :: self   !< Curve.
   integer(I4P),        intent(out), optional :: loop   !< Loop of the defect (0 if valid).
   logical,             intent(in),  optional :: check_self_intersections !< Also look for crossings.
   integer(I4P)                               :: defect !< CURVE_* code.
   integer(I4P)                               :: l, i, n !< Counters, crossing pairs.
   integer(I4P), allocatable                  :: pairs(:,:) !< Crossing pairs.

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
   if (.not. present(check_self_intersections)) return
   if (.not. check_self_intersections) return
   call self%self_intersections(pairs, n)
   if (n > 0_I4P) then
      defect = CURVE_SELF_INTERSECTING
      if (present(loop)) then
         do l = 1, self%nl
            if (pairs(1, 1) >= self%loop_start(l) .and. pairs(1, 1) < self%loop_start(l + 1)) loop = l
         enddo
      endif
   endif
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

   pure function get_next(self, i) result(j)
   !< End vertex of segment i (segment i runs from vertex i to vertex j; 0 at the end of an open loop
   !< or out of range).
   class(curve_object), intent(in) :: self !< Curve.
   integer(I4P),        intent(in) :: i    !< Segment (= its start vertex).
   integer(I4P)                    :: j    !< End vertex.

   j = 0_I4P
   if (i >= 1_I4P .and. i <= self%nv) j = self%next(i)
   endfunction get_next

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

   pure recursive function is_point_inside(self, p, sign_algorithm, pose) result(inside)
   !< True if p = (u, v) is inside the solid: signed distance < 0 with the pseudo-normal sign
   !< (default), or generalized winding number > 0.5 with SIGN_FAST_WINDING_NUMBER.
   class(curve_object), intent(in)           :: self           !< Curve.
   real(R8P),           intent(in)           :: p(2)           !< Point (u, v).
   integer(I4P),        intent(in), optional :: sign_algorithm !< Sign algorithm.
   type(rigid_pose_2d_t), intent(in), optional :: pose         !< Pose of the curve (issue #41): world-frame point.
   logical                                   :: inside         !< Result.

   if (present(pose)) then
      inside = self%is_point_inside(to_body(pose, p), sign_algorithm=sign_algorithm)
      return
   endif
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

   subroutine deviation(self, other, lower, upper, mean, rms, per_vertex, tolerance, samples_per_length, max_segments, &
                        samples, created, device, status)
   !< One-sided deviation of this curve from `other` (issue #60): h(self, other) = sup over x on self of d(x, other),
   !< bracketed as lower <= h <= upper; the 2-D twin of `surface%deviation`. `lower` from deterministic samples (every
   !< vertex, every segment midpoint, the midpoints of sub-segments, `samples_per_length` per unit length, default 4 per
   !< mean segment length); `upper` certified by halving segments (Lipschitz, convexity and one-segment bounds) until
   !< upper - lower <= `tolerance` (default 1e-3 of the larger bounding-box diagonal) or `max_segments` (default 2e6,
   !< then DEV_STATUS_NOT_CONVERGED). `mean`, `rms`: length-weighted estimates; `per_vertex`: d at every vertex. Open
   !< loops are allowed. `device=.true.`: the queries run through `query_many_device` (other entered).
   class(curve_object), intent(in)            :: self               !< Curve A.
   type(curve_object),  intent(in),   target  :: other              !< Curve B.
   real(R8P),           intent(out), optional :: lower              !< Lower bound of h(self, other).
   real(R8P),           intent(out), optional :: upper              !< Certified upper bound.
   real(R8P),           intent(out), optional :: mean, rms          !< Estimates of d over self.
   real(R8P), allocatable, intent(out), optional :: per_vertex(:)   !< d at the vertices of self.
   real(R8P),           intent(in),  optional :: tolerance          !< Target upper - lower.
   real(R8P),           intent(in),  optional :: samples_per_length !< Interior samples per unit length.
   integer(I4P),        intent(in),  optional :: max_segments       !< Work cap.
   integer(I4P),        intent(out), optional :: samples            !< Points sampled.
   integer(I4P),        intent(out), optional :: created            !< Segments created by the halving.
   logical,             intent(in),  optional :: device             !< Query other on the device.
   integer(I4P),        intent(out), optional :: status             !< DEV_STATUS_*.
   type(curve_evaluator_t)                    :: eval_other         !< other's distance query.
   real(R8P),    allocatable                  :: v(:,:), pv(:), feet(:,:) !< Vertices, distances, feet.
   integer(I4P), allocatable                  :: sg(:,:), fs(:)     !< Segments, feet's segments.
   real(R8P)                                  :: lower_, upper_, mean_, rms_, spl, tol !< Results, parameters.
   integer(I4P)                               :: samples_, created_, status_, cap, i, ns !< Results, counters.

   if (present(lower)) lower = 0._R8P
   if (present(upper)) upper = 0._R8P
   if (present(mean)) mean = 0._R8P
   if (present(rms)) rms = 0._R8P
   if (present(samples)) samples = 0_I4P
   if (present(created)) created = 0_I4P
   if (present(status)) status = DEV_STATUS_BAD_INPUT
   if (self%nv == 0_I4P .or. other%nv == 0_I4P) return
   eval_other%curve => other
   if (present(device)) eval_other%device = device
   if (eval_other%device .and. .not. other%is_on_device) return
   ns = count(self%next(1:self%nv) /= 0_I4P)
   if (ns == 0_I4P) return
   allocate(v(2, self%nv), sg(2, ns), pv(self%nv), feet(2, self%nv), fs(self%nv))
   v(1, :) = self%u(1:self%nv) ; v(2, :) = self%v(1:self%nv)
   ns = 0
   do i = 1, self%nv
      if (self%next(i) == 0_I4P) cycle
      ns = ns + 1 ; sg(:, ns) = [i, self%next(i)]
   enddo
   spl = 4._R8P * ns / max(sum(self%seg_len(1:self%nv), mask=self%next(1:self%nv) /= 0_I4P), tiny(1._R8P))
   if (present(samples_per_length)) spl = samples_per_length
   call deviation_sampled_2d(v, sg, eval_other, spl, lower_, mean_, rms_, pv, samples_, status_, vertex_feet=feet, &
                             vertex_segments=fs)
   if (present(status)) status = status_
   if (status_ /= DEV_STATUS_OK) return
   created_ = 0_I4P
   if (present(upper)) then
      tol = 1.e-3_R8P * max(norm2(self%bmax - self%bmin), norm2(other%bmax - other%bmin))
      if (present(tolerance)) tol = tolerance
      cap = 2000000_I4P ; if (present(max_segments)) cap = max_segments
      call deviation_certified_2d(v, sg, pv, feet, fs, eval_other, max(tol, tiny(1._R8P)), cap, lower_, upper_, &
                                  created_, status_)
      upper = upper_
      if (present(status)) status = status_
   endif
   if (present(lower)) lower = lower_
   if (present(mean)) mean = mean_
   if (present(rms)) rms = rms_
   if (present(per_vertex)) call move_alloc(pv, per_vertex)
   if (present(samples)) samples = samples_
   if (present(created)) created = created_
   endsubroutine deviation

   subroutine hausdorff(self, other, lower, upper, tolerance, samples_per_length, max_segments, device, status)
   !< Symmetric Hausdorff distance max(h(self, other), h(other, self)) of two curves (issue #60), bracketed.
   class(curve_object), intent(in)            :: self               !< Curve A.
   type(curve_object),  intent(in)            :: other              !< Curve B.
   real(R8P),           intent(out), optional :: lower, upper       !< Bounds of H.
   real(R8P),           intent(in),  optional :: tolerance          !< Target upper - lower of each direction.
   real(R8P),           intent(in),  optional :: samples_per_length !< Interior samples per unit length.
   integer(I4P),        intent(in),  optional :: max_segments       !< Work cap of each direction.
   logical,             intent(in),  optional :: device             !< Query on the device (both entered).
   integer(I4P),        intent(out), optional :: status             !< DEV_STATUS_*.
   real(R8P)                                  :: l1, l2, u1, u2     !< One-sided bounds.
   integer(I4P)                               :: s1, s2             !< Statuses.

   if (present(lower)) lower = 0._R8P
   if (present(upper)) upper = 0._R8P
   u1 = 0._R8P ; u2 = 0._R8P
   if (present(upper)) then
      call self%deviation(other, lower=l1, upper=u1, tolerance=tolerance, samples_per_length=samples_per_length, &
                          max_segments=max_segments, device=device, status=s1)
      call other%deviation(self, lower=l2, upper=u2, tolerance=tolerance, samples_per_length=samples_per_length, &
                           max_segments=max_segments, device=device, status=s2)
   else
      call self%deviation(other, lower=l1, samples_per_length=samples_per_length, device=device, status=s1)
      call other%deviation(self, lower=l2, samples_per_length=samples_per_length, device=device, status=s2)
   endif
   if (present(status)) status = max(s1, s2)
   if (s1 == DEV_STATUS_BAD_INPUT .or. s2 == DEV_STATUS_BAD_INPUT) return
   if (present(lower)) lower = max(l1, l2)
   if (present(upper)) upper = max(u1, u2)
   endsubroutine hausdorff

   subroutine feature_corners(self, corners, angle, status)
   !< The corners of the curve (issue #61), the 2-D feature edges: vertices of closed loops where the curve turns by more
   !< than `angle` (degrees, default 30: the angle between the incoming and the outgoing segment), and both ends of every
   !< open loop. Vertex ids as `get_vertex`.
   class(curve_object),       intent(in)            :: self    !< Curve.
   integer(I4P), allocatable, intent(out)           :: corners(:) !< Corner vertices.
   real(R8P),                 intent(in),  optional :: angle   !< Threshold [degrees].
   integer(I4P),              intent(out), optional :: status  !< STATUS_OK or STATUS_INVALID_INPUT.
   logical,      allocatable                        :: is_c(:) !< Corner flags.
   real(R8P)                                        :: a(2), b(2), c_thr, la, lb !< Work.
   integer(I4P)                                     :: i       !< Counter.

   if (present(status)) status = STATUS_INVALID_INPUT
   c_thr = feature_cos(FEATURE_EDGE_DEFAULT_ANGLE)
   if (present(angle)) then
      if (.not. (angle >= 0._R8P .and. angle <= 180._R8P)) then
         allocate(corners(0)) ; return
      endif
      c_thr = feature_cos(angle)
   endif
   allocate(is_c(self%nv))
   is_c = .false.
   do i = 1, self%nv
      if (self%prev(i) == 0_I4P .or. self%next(i) == 0_I4P) then
         is_c(i) = .true. ; cycle
      endif
      a = [self%u(i) - self%u(self%prev(i)), self%v(i) - self%v(self%prev(i))]
      b = [self%u(self%next(i)) - self%u(i), self%v(self%next(i)) - self%v(i)]
      la = norm2(a) ; lb = norm2(b)
      if (la > 0._R8P .and. lb > 0._R8P) is_c(i) = dot_product(a, b) / (la * lb) < c_thr
   enddo
   corners = pack([(i, i = 1_I4P, self%nv)], is_c)
   if (present(status)) status = STATUS_OK
   endsubroutine feature_corners

   subroutine corner_distances(self, points, angle, d, device)
   !< Distance from each point (u, v) to the nearest corner of the curve (huge without corners); with device=.true. an
   !< OpenACC loop over the points (the same arithmetic, so the same bits).
   type(curve_object),  intent(in)           :: self        !< Curve.
   real(R8P),           intent(in)           :: points(:,:) !< Points [2, n].
   real(R8P),           intent(in), optional :: angle       !< Corner threshold [degrees].
   real(R8P),           intent(out)          :: d(:)        !< Distances [n].
   logical,             intent(in), optional :: device      !< On the device.
   integer(I4P), allocatable                 :: c(:)        !< Corners.
   real(R8P),    allocatable                 :: cu(:), cv(:) !< Corner coordinates.
   real(R8P)                                 :: best        !< Work.
   integer(I4P)                              :: i, k, n, nc !< Counters.
   logical                                   :: dev         !< On the device.

   call self%feature_corners(c, angle=angle)
   nc = size(c, kind=I4P) ; n = size(points, 2, kind=I4P)
   allocate(cu(max(1_I4P, nc)), cv(max(1_I4P, nc)))
   cu = 0._R8P ; cv = 0._R8P
   if (nc > 0) then
      cu(1:nc) = self%u(c) ; cv(1:nc) = self%v(c)
   endif
   dev = .false. ; if (present(device)) dev = device
   if (dev) then
      !$acc parallel loop copyin(points, cu, cv) copyout(d) private(best, k)
      do i = 1_I4P, n
         best = huge(1._R8P)
         do k = 1_I4P, nc
            best = min(best, (points(1, i) - cu(k)) ** 2 + (points(2, i) - cv(k)) ** 2)
         enddo
         d(i) = huge(1._R8P)
         if (best < huge(1._R8P)) d(i) = sqrt(best)
      enddo
   else
      !$omp parallel do default(shared) private(i, k, best)
      do i = 1_I4P, n
         best = huge(1._R8P)
         do k = 1_I4P, nc
            best = min(best, (points(1, i) - cu(k)) ** 2 + (points(2, i) - cv(k)) ** 2)
         enddo
         d(i) = huge(1._R8P)
         if (best < huge(1._R8P)) d(i) = sqrt(best)
      enddo
      !$omp end parallel do
   endif
   endsubroutine corner_distances

   subroutine curve_evaluator_eval(self, points, distances, closest, segments)
   !< Unsigned distances of the points to the evaluator's curve, the closest points and their segments (`query_many`).
   class(curve_evaluator_t), intent(in)  :: self         !< Evaluator.
   real(R8P),                intent(in)  :: points(:,:)  !< Points [2, n].
   real(R8P),                intent(out) :: distances(:) !< Distances [n].
   real(R8P),                intent(out) :: closest(:,:) !< Closest points [2, n].
   integer(I4P),             intent(out) :: segments(:)  !< Closest segments [n].

   if (self%device) then
      call self%curve%query_many_device(points=points, distances=distances, closest_points=closest, segments=segments, &
                                        is_signed=.false.)
   else
      call self%curve%query_many(points=points, distances=distances, closest_points=closest, segments=segments, &
                                 is_signed=.false.)
   endif
   endsubroutine curve_evaluator_eval

   subroutine load_from_file(self, file_name, normal_axis, status)
   !< Load the curve from the ASCII loop format.
   !<
   !< One vertex `u v` per line; loops are separated by one or more blank lines; text after `#`
   !< is a comment, except the directive `# open` before a loop's first vertex, which marks the loop
   !< open (a polyline), and `# patch <name>` (issue #66), which puts the segments starting at the
   !< vertices that follow in the named patch (it does not end the loop). A closed loop whose last vertex repeats its first exactly is closed
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
   integer(I4P), allocatable                  :: vp(:)       !< Patch of each vertex's segment (issue #66).
   character(64), allocatable                 :: names(:)    !< Patch names met.
   integer(I4P)                               :: p           !< Current patch (0: none yet).
   character(1024)                            :: directive   !< Text after '#'.
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
   allocate(u(64), v(64), sizes(8), closed(8), vp(64), names(0))
   n = 0_I4P ; nl = 0_I4P ; cur = 0_I4P ; open_next = .false. ; p = 0_I4P
   do
      read(unit, '(A)', iostat=ios) line
      if (ios /= 0) exit
      h = index(line, '#')
      if (h > 0) then
         if (adjustl(line(h+1:)) == 'open' .and. len_trim(line(1:h-1)) == 0) open_next = .true.
         directive = adjustl(line(h+1:))
         if (directive(1:6) == 'patch ' .and. len_trim(line(1:h-1)) == 0) then
            p = name_index(trim(adjustl(directive(7:))))
            cycle ! a patch directive does not end the loop (other comment-only lines still do)
         endif
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
         u = [u, u] ; v = [v, v] ; vp = [vp, vp]
      endif
      if (p == 0_I4P) p = name_index('')
      n = n + 1_I4P ; u(n) = a ; v(n) = b ; vp(n) = p ; cur = cur + 1_I4P
   enddo
   close(unit)
   if (status_ == STATUS_OK) then
      call close_loop
      if (size(names) == 0) names = [repeat(' ', 64)]
      call self%adopt_loops(u=u(1:n), v=v(1:n), loop_sizes=sizes(1:nl), normal_axis=normal_axis, closed=closed(1:nl), &
                            status=status_, patches=vp(1:n), patch_names=names)
   endif
   if (present(status)) status = status_
   contains
      function name_index(nm) result(q)
      !< Index of patch `nm`, added if new.
      character(*), intent(in) :: nm
      integer(I4P)             :: q
      character(64)            :: pad

      pad = nm
      do q = 1, size(names)
         if (names(q) == pad) return
      enddo
      names = [names, pad]
      q = size(names)
      endfunction name_index

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

   subroutine mass_properties(self, mass, center, second_moments, polar_moment, principal_moments, principal_angle, &
                              density, about, area, status)
   !< Mass properties per unit depth of the region bounded by the curve (issue #59): mass, centroid, second moments
   !< [I_uu, I_vv, I_uv] = density * [int (u - p_u)^2, int (v - p_v)^2, int (u - p_u)(v - p_v)] dA about p = the
   !< centroid (or `about`), the polar moment I_uu + I_vv (the inertia for a rotation in the plane), the principal values
   !< of the second-moment matrix (ascending) and the angle (radians, in (-pi/2, pi/2], from the u axis) of the axis of
   !< the larger one.
   !<
   !< Exact polygon integrals by Green's theorem (`fossil_mass_properties`), summed about the bounding-box centre with
   !< compensation: the only error is round-off. Holes, wound clockwise, subtract. Every loop must be closed and the
   !< enclosed area positive: otherwise `status` is `STATUS_INVALID_INPUT` and the outputs are zero.
   class(curve_object), intent(in)            :: self                 !< Curve.
   real(R8P),           intent(out), optional :: mass                 !< Mass per unit depth (density * area).
   real(R8P),           intent(out), optional :: center(2)            !< Centroid (u, v).
   real(R8P),           intent(out), optional :: second_moments(3)    !< [I_uu, I_vv, I_uv].
   real(R8P),           intent(out), optional :: polar_moment         !< I_uu + I_vv.
   real(R8P),           intent(out), optional :: principal_moments(2) !< Principal values, ascending.
   real(R8P),           intent(out), optional :: principal_angle      !< Angle of the axis of the larger value.
   real(R8P),           intent(in),  optional :: density              !< Uniform density per unit area (default 1).
   real(R8P),           intent(in),  optional :: about(2)             !< Point of the moments (default: the centroid).
   real(R8P),           intent(out), optional :: area                 !< Enclosed area.
   integer(I4P),        intent(out), optional :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   type(moments_accumulator_t)                :: acc                  !< Raw moments.
   real(R8P)                                  :: o(2)                 !< Reference point (bounding-box centre).
   real(R8P)                                  :: density_             !< Density, local.
   real(R8P)                                  :: area_, mass_         !< Area, mass, local.
   real(R8P)                                  :: center_(2)           !< Centroid, local.
   real(R8P)                                  :: sm(3)                !< Second moments, local.
   real(R8P)                                  :: w(2), angle          !< Principal values and angle, local.
   integer(I4P)                               :: i                    !< Counter.

   if (present(mass)) mass = 0._R8P
   if (present(center)) center = 0._R8P
   if (present(second_moments)) second_moments = 0._R8P
   if (present(polar_moment)) polar_moment = 0._R8P
   if (present(principal_moments)) principal_moments = 0._R8P
   if (present(principal_angle)) principal_angle = 0._R8P
   if (present(area)) area = 0._R8P
   if (present(status)) status = STATUS_INVALID_INPUT
   if (self%nv == 0_I4P .or. self%nl == 0_I4P) return
   if (.not. all(self%loop_closed)) return
   density_ = 1._R8P ; if (present(density)) density_ = density
   o = 0.5_R8P * (self%bmin + self%bmax)
   do i = 1, self%nv
      call accumulate_segment(acc, [self%u(i), self%v(i)] - o, [self%u(self%next(i)), self%v(self%next(i))] - o)
   enddo
   if (.not. (acc%s(1) + acc%c(1) > 0._R8P)) return
   call moments_to_properties_2d(acc, o, density_, about, area_, mass_, center_, sm)
   if (present(mass)) mass = mass_
   if (present(center)) center = center_
   if (present(second_moments)) second_moments = sm
   if (present(polar_moment)) polar_moment = sm(1) + sm(2)
   if (present(area)) area = area_
   if (present(principal_moments) .or. present(principal_angle)) then
      if (present(about)) call moments_to_properties_2d(acc, o, density_, area=area_, mass=mass_, center=center_, &
                                                        second_moments=sm)
      call symmetric2_eigen(sm, w, angle)
      if (present(principal_moments)) principal_moments = w
      if (present(principal_angle)) principal_angle = angle
   endif
   if (present(status)) status = STATUS_OK
   endsubroutine mass_properties

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
      integer(I4P), allocatable :: p(:)    !< Patch of the segment starting at each vertex (issue #66).
      logical                :: closed     !< Closed.
      logical                :: alive      !< Not merged into another chain.
   endtype chain_t
   type(chain_t), allocatable :: ch(:)        !< Loops as chains.
   real(R8P),     allocatable :: nu(:), nv(:) !< Rebuilt vertices.
   integer(I4P),  allocatable :: np_(:), tp(:) !< Rebuilt and work patches.
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
      ch(l)%p = self%seg_patch(self%loop_start(l):self%loop_start(l + 1) - 1)
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
               ch(best)%p = reversed_patches(ch(best)%p)
            endif
            ch(l)%u = [ch(l)%u, ch(best)%u] ; ch(l)%v = [ch(l)%v, ch(best)%v]
            ! the joining segment (end of l -> start of best) takes the patch of l's last segment
            if (size(ch(l)%p) > 1) ch(l)%p(size(ch(l)%p)) = ch(l)%p(size(ch(l)%p) - 1)
            ch(l)%p = [ch(l)%p, ch(best)%p]
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
      allocate(tu(n), tv(n), tp(n))
      k = 1_I4P ; tu(1) = ch(l)%u(1) ; tv(1) = ch(l)%v(1) ; tp(1) = ch(l)%p(1)
      do i = 2, n
         if (hypot(ch(l)%u(i) - tu(k), ch(l)%v(i) - tv(k)) <= tolerance) then
            merged = merged + 1_I4P
            tp(k) = ch(l)%p(i) ! the segment leaving the merged vertex survives
         else
            k = k + 1_I4P ; tu(k) = ch(l)%u(i) ; tv(k) = ch(l)%v(i) ; tp(k) = ch(l)%p(i)
         endif
      enddo
      if (ch(l)%closed .and. k > 1_I4P) then
         if (hypot(tu(k) - tu(1), tv(k) - tv(1)) <= tolerance) then
            merged = merged + 1_I4P ; k = k - 1_I4P
         endif
      endif
      ch(l)%u = tu(1:k) ; ch(l)%v = tv(1:k) ; ch(l)%p = tp(1:k)
      deallocate(tu, tv, tp)
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
            np_(first:last) = reversed_patches(np_(first:last))
            reversed = reversed + 1_I4P
         endif
      endif
      first = last + 1_I4P
   enddo
   block
      character(64), allocatable :: names(:)
      names = self%patch_names
      call self%adopt_loops(u=nu, v=nv, loop_sizes=sizes, normal_axis=axis, closed=closed, patches=np_, patch_names=names)
   endblock
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
      allocate(nu(nt), nv(nt), sizes(nl2), closed(nl2), np_(nt))
      nt = 0_I4P ; nl2 = 0_I4P
      do q = 1, size(ch)
         if (.not. ch(q)%alive) cycle
         nl2 = nl2 + 1_I4P
         sizes(nl2) = size(ch(q)%u, kind=I4P) ; closed(nl2) = ch(q)%closed
         nu(nt + 1:nt + sizes(nl2)) = ch(q)%u ; nv(nt + 1:nt + sizes(nl2)) = ch(q)%v
         np_(nt + 1:nt + sizes(nl2)) = ch(q)%p
         nt = nt + sizes(nl2)
      enddo
      endsubroutine pack_chains

      pure function reversed_patches(pp) result(rp)
      !< Segment patches of a reversed vertex sequence: segment j of the reversed sequence is segment n-j of the
      !< original (the last entry, an open chain's unused one or a closed loop's closing segment, stays).
      integer(I4P), intent(in)  :: pp(:)
      integer(I4P)              :: rp(size(pp))
      integer(I4P)              :: m

      m = size(pp)
      rp = pp
      if (m > 1) rp(1:m - 1) = pp(m - 1:1:-1)
      endfunction reversed_patches
   endsubroutine repair

   subroutine save_into_file(self, file_name, status)
   !< Save the curve into the ASCII loop format (see `load_from_file`), bit-exact on reload.
   class(curve_object), intent(in)            :: self      !< Curve.
   character(*),        intent(in)            :: file_name !< File name.
   integer(I4P),        intent(out), optional :: status    !< STATUS_OK or STATUS_FILE_OPEN_FAIL.
   integer(I4P)                               :: unit, ios, l, i, last !< Work.
   logical                                    :: named   !< Write patch directives (issue #66).

   if (present(status)) status = STATUS_OK
   named = self%get_patches_number() > 1_I4P
   if (.not. named .and. allocated(self%patch_names)) named = len_trim(self%patch_names(1)) > 0
   open(newunit=unit, file=file_name, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      if (present(status)) status = STATUS_FILE_OPEN_FAIL
      return
   endif
   write(unit, '(A,I0,A,I0)') '# FOSSIL curve: loops ', self%nl, ', normal_axis ', self%normal_axis
   do l = 1, self%nl
      if (l > 1_I4P) write(unit, '(A)') ''
      if (.not. self%loop_closed(l)) write(unit, '(A)') '# open'
      last = 0_I4P
      do i = self%loop_start(l), self%loop_start(l + 1) - 1
         if (named .and. self%segment_patch(i) /= last) then
            write(unit, '(A)') '# patch '//self%patch_name(self%segment_patch(i))
            last = self%segment_patch(i)
         endif
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
   subroutine wall_distance_many(self, points, distances, bodies, cutoff, status)
   !< Unsigned distance to the nearest wall and its body, with a cutoff (issue #42): the surface twin
   !< (distance = cutoff, body 0 beyond it), from `query_many` (no early exit in 2-D: cheap).
   class(curve_object), intent(in)            :: self         !< Curve.
   real(R8P),           intent(in)            :: points(:,:)  !< Points [2, n].
   real(R8P),           intent(out)           :: distances(:) !< Distances [n].
   integer(I4P),        intent(out)           :: bodies(:)    !< Bodies [n].
   real(R8P),           intent(in),  optional :: cutoff       !< Cutoff (default: none).
   integer(I4P),        intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P) :: co !< Cutoff.

   co = huge(1._R8P) ; if (present(cutoff)) co = cutoff
   call self%query_many(points=points, distances=distances, bodies=bodies, is_signed=.false., status=status)
   if (present(status)) then
      if (status /= STATUS_OK) return
   endif
   where (distances >= co)
      distances = co ; bodies = 0_I4P
   endwhere
   endsubroutine wall_distance_many

   subroutine state_changes(self, pose0, pose1, points, fresh, dead, sign_algorithm, checked, status)
   !< Fresh (solid -> fluid) and dead (fluid -> solid) points (u, v) between two in-plane poses (issue
   !< #41): the surface twin, with the same band filter (|d(b0)| > |b0 - b1|: no change).
   class(curve_object),   intent(in)            :: self           !< Curve (body frame).
   type(rigid_pose_2d_t), intent(in)            :: pose0, pose1   !< Poses.
   real(R8P),             intent(in)            :: points(:,:)    !< World points [2, n].
   logical,               intent(out)           :: fresh(:)       !< Solid -> fluid [n].
   logical,               intent(out)           :: dead(:)        !< Fluid -> solid [n].
   integer(I4P),          intent(in),  optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   integer(I4P),          intent(out), optional :: checked        !< Points that needed the second query.
   integer(I4P),          intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: i, n, nchk !< Work.
   real(R8P)    :: b0(2), b1(2), d0 !< Body points, signed distance.
   logical      :: in0, in1 !< States.

   if (present(status)) status = STATUS_OK
   n = size(points, 2, kind=I4P)
   if (size(points, 1) /= 2 .or. size(fresh) /= n .or. size(dead) /= n) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   nchk = 0_I4P
   !$omp parallel do default(shared) schedule(dynamic, 64) private(i, b0, b1, d0, in0, in1) reduction(+:nchk)
   do i = 1_I4P, n
      b0 = to_body(pose0, points(:, i)) ; b1 = to_body(pose1, points(:, i))
      d0 = self%distance(b0, sign_algorithm=sign_algorithm)
      in0 = d0 < 0._R8P
      fresh(i) = .false. ; dead(i) = .false.
      if (abs(d0) > norm2(b1 - b0)) cycle
      nchk = nchk + 1_I4P
      in1 = self%distance(b1, sign_algorithm=sign_algorithm) < 0._R8P
      fresh(i) = in0 .and. .not. in1
      dead(i) = in1 .and. .not. in0
   enddo
   !$omp end parallel do
   if (present(checked)) checked = nchk
   endsubroutine state_changes

   subroutine move(self, pose)
   !< Move the curve rigidly by the in-plane `pose` ((u, v) <- R (u, v) + t) and refit its segment
   !< tree instead of rebuilding it (issue #41): the coordinate-derived data are recomputed as
   !< `analyze` does (`curve_metrics`), the topology and the bodies are kept, the tree is refitted
   !< (`refit_vertices`, O(N)). Queries then equal those on the same curve moved and rebuilt. A device
   !< copy is refreshed.
   class(curve_object),   intent(inout) :: self !< Curve.
   type(rigid_pose_2d_t), intent(in)    :: pose !< Pose.
   integer(I4P)                         :: i    !< Counter.
   real(R8P)                            :: w(2) !< Moved vertex.
   logical                              :: was_on_device !< Device copy to refresh.

   if (self%nv == 0_I4P) return
   was_on_device = self%is_on_device
   if (was_on_device) call self%exit_device
   do i = 1_I4P, self%nv
      w = to_world(pose, [self%u(i), self%v(i)]) ; self%u(i) = w(1) ; self%v(i) = w(2)
   enddo
   call curve_metrics(self)
   call self%bvh%refit_vertices(self%u, self%v, self%next)
   if (was_on_device) call self%enter_device
   endsubroutine move

   pure subroutine curve_metrics(self)
   !< Coordinate-derived data of `analyze`: segment lengths and normals, vertex pseudo-normals,
   !< perimeter, loop areas, centroid, bounding box (also used by `move`, issue #41).
   type(curve_object), intent(inout) :: self !< Curve.
   integer(I4P)                      :: i, j !< Counters.
   real(R8P)                         :: du, dv, cross, cu, cv, s !< Work.

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
   endsubroutine curve_metrics

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
   integer(I4P)                       :: l, i, first, last !< Counters.

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
   call curve_metrics(self)
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

   pure function rich_sizes_ok(n, closest_points, segments, vertices, normals, bodies, patches) result(good)
   !< True if every present rich output has the right size for n points.
   integer(I4P),           intent(in) :: n                   !< Points.
   real(R8P),    optional, intent(in) :: closest_points(:,:) !< [2, n].
   integer(I4P), optional, intent(in) :: segments(:)         !< [n].
   integer(I4P), optional, intent(in) :: vertices(:)         !< [n].
   real(R8P),    optional, intent(in) :: normals(:,:)        !< [2, n].
   integer(I4P), optional, intent(in) :: bodies(:)           !< [n].
   integer(I4P), optional, intent(in) :: patches(:)          !< [n].
   logical                            :: good                !< Result.

   good = .true.
   if (present(closest_points)) good = good .and. size(closest_points, 1) == 2 .and. size(closest_points, 2) == n
   if (present(normals)) good = good .and. size(normals, 1) == 2 .and. size(normals, 2) == n
   if (present(segments)) good = good .and. size(segments) == n
   if (present(vertices)) good = good .and. size(vertices) == n
   if (present(patches)) good = good .and. size(patches) == n
   if (present(bodies)) good = good .and. size(bodies) == n
   endfunction rich_sizes_ok

   function curve_many_ok(self, axis, a, t0, t1, t, segments, entering, counts, overflow, status) result(good)
   !< Argument check of `line_crossings_many(_device)` (issue #40); outputs zeroed when good.
   type(curve_object), intent(in)            :: self          !< Curve.
   integer(I4P),       intent(in)            :: axis          !< 1 or 2.
   real(R8P),          intent(in)            :: a(:)          !< Line positions.
   real(R8P),          intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),          intent(inout)         :: t(:,:)        !< Parameters.
   integer(I4P),       intent(inout)         :: segments(:,:) !< Segments.
   logical,            intent(inout)         :: entering(:,:) !< Entering flags.
   integer(I4P),       intent(inout)         :: counts(:)     !< Counts.
   logical,            intent(inout)         :: overflow(:)   !< Overflow.
   integer(I4P),       intent(out), optional :: status        !< Status.
   logical                                   :: good          !< Arguments valid.
   integer(I4P)                              :: n             !< Lines.

   if (present(status)) status = STATUS_OK
   n = size(a, kind=I4P)
   good = (axis == 1_I4P .or. axis == 2_I4P) .and. t1 >= t0 .and. size(t, 2) == n .and. size(segments, 2) == n .and. &
          size(entering, 2) == n .and. size(counts) == n .and. size(overflow) == n .and. size(segments, 1) == size(t, 1) &
          .and. size(entering, 1) == size(t, 1) .and. size(t, 1) > 0 .and. self%bvh%nodes_number > 0_I4P
   if (.not. good) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   t = 0._R8P ; segments = 0_I4P ; entering = .false. ; counts = 0_I4P ; overflow = .false.
   endfunction curve_many_ok

   subroutine curve_crossings_core(self, axis, a, t0, t1, t, ids, es, n)
   !< Sorted crossings of an axis-aligned line with the curve (issue #40): segment-BVH kernel with a
   !< growing buffer.
   type(curve_object),        intent(in)  :: self   !< Curve.
   integer(I4P),              intent(in)  :: axis   !< 1 or 2.
   real(R8P),                 intent(in)  :: a      !< Line position.
   real(R8P),                 intent(in)  :: t0, t1 !< Parameter range.
   real(R8P),    allocatable, intent(out) :: t(:)   !< Parameters.
   integer(I4P), allocatable, intent(out) :: ids(:) !< Segments.
   logical,      allocatable, intent(out) :: es(:)  !< Entering flags.
   integer(I4P),              intent(out) :: n      !< Count.
   integer(I4P) :: nmax     !< Buffer size.
   logical      :: overflow !< Buffer too small.

   n = 0_I4P
   if (self%bvh%nodes_number == 0_I4P) then
      allocate(t(0), ids(0), es(0))
      return
   endif
   nmax = 64_I4P
   do
      if (allocated(t)) deallocate(t, ids, es)
      allocate(t(nmax), ids(nmax), es(nmax))
      call curve_line_crossings_oac(self%bvh%lo_u, self%bvh%lo_v, self%bvh%hi_u, self%bvh%hi_v, self%bvh%left, &
                                    self%bvh%right, self%bvh%first, self%bvh%last, self%bvh%seg_id, self%next, self%u, &
                                    self%v, axis, a, t0, t1, t, ids, es, n, overflow)
      if (.not. overflow) exit
      nmax = 4_I4P * nmax
   enddo
   t = t(1:n) ; ids = ids(1:n) ; es = es(1:n)
   endsubroutine curve_crossings_core

   pure function body_points2(pose, points) result(bp)
   !< World-frame points (u, v) mapped to the body frame of `pose` (issue #41).
   type(rigid_pose_2d_t), intent(in) :: pose        !< Pose.
   real(R8P),             intent(in) :: points(:,:) !< World points [2, n].
   real(R8P)                         :: bp(2, size(points, 2)) !< Body points.
   integer(I4P)                      :: i           !< Counter.

   do i = 1_I4P, size(points, 2, kind=I4P)
      bp(:, i) = to_body(pose, points(:, i))
   enddo
   endfunction body_points2

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
