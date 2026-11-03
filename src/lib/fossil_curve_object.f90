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
use fossil_predicates,         only : orient2d_xy
use fossil_segment_bvh,        only : segment_bvh, FEATURE_VERTEX_START, FEATURE_VERTEX_END
use fossil_surface_stl_object, only : STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL, &
                                     SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use penf,                      only : I4P, R8P

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
   real(R8P)                 :: bmin(2) = 0._R8P    !< Bounding box minimum (u, v).
   real(R8P)                 :: bmax(2) = 0._R8P    !< Bounding box maximum (u, v).
   real(R8P)                 :: area = 0._R8P       !< Total signed area (solid area for a valid curve).
   real(R8P)                 :: perimeter = 0._R8P  !< Total length.
   real(R8P)                 :: centroid(2) = 0._R8P !< Centroid of the solid region (u, v).
   type(segment_bvh)         :: bvh                 !< SAH BVH over the segments.
   contains
      ! public methods
      procedure, pass(self) :: adopt_loops          !< Build the curve from vertex arrays and loop sizes.
      procedure, pass(self) :: check_validity       !< First validity defect (CURVE_* code).
      procedure, pass(self) :: compute_distance     !< Distance with closest point, segment and vertex.
      procedure, pass(self) :: distance             !< (Signed) distance of a point.
      procedure, pass(self) :: destroy              !< Reset to an empty curve.
      procedure, pass(self) :: get_area             !< Total signed area.
      procedure, pass(self) :: get_bmax             !< Bounding box maximum.
      procedure, pass(self) :: get_bmin             !< Bounding box minimum.
      procedure, pass(self) :: get_centroid         !< Centroid of the solid region.
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
                                    use_brute_force, sign_algorithm)
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
