!< FOSSIL, bounding volume hierarchy over 2-D segments (issue #36 step 3).
module fossil_segment_bvh
!< FOSSIL, bounding volume hierarchy over 2-D segments (issue #36 step 3).
!<
!< A binary BVH built top-down with the bucketed surface-area heuristic (16 buckets per axis;
!< in 2-D the "surface area" of a box is its half-perimeter), the 2-D twin of the SAH BVH of
!< `fossil_aabb_tree_object`. Everything is stored in flat arrays (public, so a later step can
!< map them to the device): node boxes, children, leaf ranges, and the segment payload in leaf
!< order (start point, direction, squared length, id of the segment in the curve).
!<
!< Queries are iterative (explicit stack, no recursion) and deterministic: nearer child first,
!< pruning by squared distance, a candidate replaces the best only if strictly closer. The tree
!< depth is capped (`MAX_SAH_DEPTH`, then median splits), so the fixed traversal stack
!< (`STACK_SIZE`) never overflows.
!<
!< Winding number (issue #36 step 4). The generalized winding number of a point p is
!< w = (1/2pi) sum of the angles the segments subtend at p (exactly 1 inside and 0 outside a
!< closed counter-clockwise loop; fractional near an open end). `winding_number` evaluates it
!< hierarchically, the 2-D analogue of the fast winding number (Barill et al. 2018): a node far
!< enough from p (|c - p| > beta * diagonal) contributes the second-order expansion of
!< (1/2pi) integral n . (x - p) / |x - p|^2 ds about its box centre c, from its dipole moment
!< D = sum length * n and second moment C(i,k) = sum length * n(i) * (m - c)(k) (m segment
!< midpoint; both exact for straight segments); other nodes recurse, leaves sum exact angles.
!<
!< Device (issue #36 step 7). The traversals are free `_oac` procedures (`!$acc routine seq`) on
!< the flat arrays, and the type-bound `closest` and `winding_number` call them: the CPU and the
!< device run the same code, as the 3-D `bvh_winding_number_oac` does (#20 rules: no type-bound
!< calls or `class` dummies in device code).
use penf, only : I4P, R8P, MaxR8P

implicit none
private
public :: segment_bvh
public :: FEATURE_NONE, FEATURE_SEGMENT, FEATURE_VERTEX_START, FEATURE_VERTEX_END
public :: segment_point_distance2, segment_winding
public :: segment_closest_oac, segment_winding_number_oac, STACK_SIZE

integer(I4P), parameter :: LEAF_SIZE     = 4_I4P   !< Maximum segments per leaf.
integer(I4P), parameter :: BUCKETS       = 16_I4P  !< SAH buckets per axis.
integer(I4P), parameter :: MAX_SAH_DEPTH = 48_I4P  !< Below this depth, median splits only.
integer(I4P), parameter :: STACK_SIZE    = 128_I4P !< Traversal stack (> MAX_SAH_DEPTH + log2(huge)).

! Closest feature of a segment: its interior, or one of its end vertices.
integer(I4P), parameter :: FEATURE_NONE         = 0_I4P !< No feature (empty tree).
integer(I4P), parameter :: FEATURE_SEGMENT      = 1_I4P !< Segment interior.
integer(I4P), parameter :: FEATURE_VERTEX_START = 2_I4P !< Start vertex of the segment.
integer(I4P), parameter :: FEATURE_VERTEX_END   = 3_I4P !< End vertex of the segment.

type :: segment_bvh
   !< Flat SAH BVH over 2-D segments. Node 1 is the root; a leaf has `left == 0`.
   integer(I4P)              :: nodes_number = 0_I4P !< Number of nodes.
   integer(I4P)              :: segments_number = 0_I4P !< Number of segments.
   integer(I4P)              :: depth = 0_I4P        !< Tree depth (root = 1).
   real(R8P),    allocatable :: lo_u(:), lo_v(:)     !< Node box minimum [nodes].
   real(R8P),    allocatable :: hi_u(:), hi_v(:)     !< Node box maximum [nodes].
   integer(I4P), allocatable :: left(:), right(:)    !< Children [nodes] (0 for a leaf).
   integer(I4P), allocatable :: first(:), last(:)    !< Leaf range in payload order [nodes].
   real(R8P),    allocatable :: a_u(:), a_v(:)       !< Segment start [segments, payload order].
   real(R8P),    allocatable :: d_u(:), d_v(:)       !< Segment direction (end - start) [segments].
   real(R8P),    allocatable :: len2(:)              !< Squared length [segments].
   integer(I4P), allocatable :: seg_id(:)            !< Segment index in the curve [segments].
   real(R8P),    allocatable :: wn_cu(:), wn_cv(:)   !< Winding cache: expansion centre [nodes].
   real(R8P),    allocatable :: wn_du(:), wn_dv(:)   !< Winding cache: dipole moment [nodes].
   real(R8P),    allocatable :: wn_c(:,:)            !< Winding cache: second moment [4, nodes] (C11, C21, C12, C22).
   real(R8P),    allocatable :: wn_diag2(:)          !< Winding cache: squared box diagonal [nodes].
   contains
      procedure, pass(self) :: build             !< Build over segments.
      procedure, pass(self) :: closest           !< Closest segment, point and feature.
      procedure, pass(self) :: closest_brute     !< Same, brute force (reference).
      procedure, pass(self) :: destroy           !< Reset.
      procedure, pass(self) :: segments_in_box   !< Segments whose box overlaps a box.
      procedure, pass(self) :: refit_vertices    !< Refit to moved vertices, same topology (issue #41).
      procedure, pass(self) :: winding_number    !< Fast (hierarchical) winding number.
      procedure, pass(self) :: winding_number_brute !< Exact sum of the subtended angles.
endtype segment_bvh

contains
   pure subroutine build(self, au, av, bu, bv, ids)
   !< Build the tree over the segments (au(i), av(i)) -> (bu(i), bv(i)), i = 1..n.
   !<
   !< `seg_id` stores ids(i) (default i) for segment i: the caller's numbering. The winding cache
   !< is built too.
   class(segment_bvh), intent(inout) :: self           !< Tree.
   real(R8P),          intent(in)    :: au(:), av(:)   !< Segment starts.
   real(R8P),          intent(in)    :: bu(:), bv(:)   !< Segment ends.
   integer(I4P),       intent(in), optional :: ids(:) !< Caller's segment ids (default 1..n).
   integer(I4P), allocatable         :: perm(:)        !< Segments in leaf order.
   real(R8P),    allocatable         :: cu(:), cv(:)   !< Segment box centres.
   real(R8P),    allocatable         :: sb(:,:)        !< Segment boxes [lo_u, lo_v, hi_u, hi_v; n].
   integer(I4P), allocatable         :: st_node(:), st_lo(:), st_hi(:), st_depth(:) !< Build stack.
   integer(I4P)                      :: n, top, node, lo, hi, dep, mid, i, k !< Work.

   call self%destroy
   n = size(au, kind=I4P)
   self%segments_number = n
   if (n == 0_I4P) return
   allocate(perm(n), cu(n), cv(n), sb(4, n))
   do i = 1, n
      perm(i) = i
      cu(i) = 0.5_R8P * (min(au(i), bu(i)) + max(au(i), bu(i)))
      cv(i) = 0.5_R8P * (min(av(i), bv(i)) + max(av(i), bv(i)))
      sb(:, i) = [min(au(i), bu(i)), min(av(i), bv(i)), max(au(i), bu(i)), max(av(i), bv(i))]
   enddo
   allocate(self%lo_u(2 * n), self%lo_v(2 * n), self%hi_u(2 * n), self%hi_v(2 * n))
   allocate(self%left(2 * n), self%right(2 * n), self%first(2 * n), self%last(2 * n))
   allocate(st_node(2 * n), st_lo(2 * n), st_hi(2 * n), st_depth(2 * n))
   self%nodes_number = 1_I4P
   top = 1_I4P ; st_node(1) = 1_I4P ; st_lo(1) = 1_I4P ; st_hi(1) = n ; st_depth(1) = 1_I4P
   do while (top > 0_I4P)
      node = st_node(top) ; lo = st_lo(top) ; hi = st_hi(top) ; dep = st_depth(top) ; top = top - 1_I4P
      self%depth = max(self%depth, dep)
      self%lo_u(node) = minval(min(au(perm(lo:hi)), bu(perm(lo:hi))))
      self%lo_v(node) = minval(min(av(perm(lo:hi)), bv(perm(lo:hi))))
      self%hi_u(node) = maxval(max(au(perm(lo:hi)), bu(perm(lo:hi))))
      self%hi_v(node) = maxval(max(av(perm(lo:hi)), bv(perm(lo:hi))))
      self%first(node) = lo ; self%last(node) = hi
      self%left(node) = 0_I4P ; self%right(node) = 0_I4P
      if (hi - lo + 1_I4P <= LEAF_SIZE) cycle
      call split(perm(lo:hi), cu, cv, sb, dep >= MAX_SAH_DEPTH, mid)
      mid = lo + mid - 1_I4P
      self%left(node)  = self%nodes_number + 1_I4P
      self%right(node) = self%nodes_number + 2_I4P
      self%nodes_number = self%nodes_number + 2_I4P
      top = top + 1_I4P ; st_node(top) = self%right(node) ; st_lo(top) = mid + 1_I4P ; st_hi(top) = hi
      st_depth(top) = dep + 1_I4P
      top = top + 1_I4P ; st_node(top) = self%left(node) ; st_lo(top) = lo ; st_hi(top) = mid
      st_depth(top) = dep + 1_I4P
   enddo
   allocate(self%a_u(n), self%a_v(n), self%d_u(n), self%d_v(n), self%len2(n), self%seg_id(n))
   do k = 1, n
      i = perm(k)
      self%seg_id(k) = i
      if (present(ids)) self%seg_id(k) = ids(i)
      self%a_u(k) = au(i) ; self%a_v(k) = av(i)
      self%d_u(k) = bu(i) - au(i) ; self%d_v(k) = bv(i) - av(i)
      self%len2(k) = self%d_u(k)**2 + self%d_v(k)**2
   enddo
   call build_winding_cache(self)
   endsubroutine build

   pure subroutine closest(self, pu, pv, d2, k, t, feature)
   !< Closest segment to the point (pu, pv): squared distance, payload index, parameter, feature.
   !<
   !< `k` indexes the payload (`seg_id(k)` is the curve segment); `t` in [0, 1] is the parameter
   !< of the closest point a + t d (exactly 0 or 1 at a vertex). With an empty tree: d2 = huge,
   !< k = 0, feature = FEATURE_NONE. Runs `segment_closest_oac`, the device kernel.
   class(segment_bvh), intent(in)  :: self    !< Tree.
   real(R8P),          intent(in)  :: pu, pv  !< Query point.
   real(R8P),          intent(out) :: d2      !< Squared distance.
   integer(I4P),       intent(out) :: k       !< Payload index of the closest segment.
   real(R8P),          intent(out) :: t       !< Parameter of the closest point.
   integer(I4P),       intent(out) :: feature !< FEATURE_* of the closest point.

   if (self%nodes_number == 0_I4P) then
      d2 = MaxR8P ; k = 0_I4P ; t = 0._R8P ; feature = FEATURE_NONE
      return
   endif
   call segment_closest_oac(self%lo_u, self%lo_v, self%hi_u, self%hi_v, self%left, self%right, self%first, self%last, &
                            self%a_u, self%a_v, self%d_u, self%d_v, self%len2, pu, pv, d2, k, t, feature)
   endsubroutine closest

   pure subroutine segment_closest_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, len2, &
                                       pu, pv, d2, k, t, feature)
   !< Closest segment to (pu, pv) on the flat tree arrays (root = node 1, non-empty tree).
   !<
   !< Nearer child first, pruning on squared box distance, a candidate wins only if strictly
   !< closer. CPU and device kernel.
   !$acc routine seq
   real(R8P),    intent(in)  :: lo_u(:), lo_v(:), hi_u(:), hi_v(:) !< Node boxes.
   integer(I4P), intent(in)  :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   real(R8P),    intent(in)  :: a_u(:), a_v(:), d_u(:), d_v(:), len2(:) !< Segment payload.
   real(R8P),    intent(in)  :: pu, pv  !< Query point.
   real(R8P),    intent(out) :: d2      !< Squared distance.
   integer(I4P), intent(out) :: k       !< Payload index.
   real(R8P),    intent(out) :: t       !< Parameter.
   integer(I4P), intent(out) :: feature !< FEATURE_*.
   integer(I4P)              :: stack(STACK_SIZE), top, node, c1, c2, j, f !< Work.
   real(R8P)                 :: b1, b2, dj, tj !< Work.

   d2 = MaxR8P ; k = 0_I4P ; t = 0._R8P ; feature = FEATURE_NONE
   top = 1_I4P ; stack(1) = 1_I4P
   do while (top > 0_I4P)
      node = stack(top) ; top = top - 1_I4P
      if (box_d2_oac(lo_u(node), lo_v(node), hi_u(node), hi_v(node), pu, pv) >= d2) cycle
      if (left(node) == 0_I4P) then
         do j = first(node), last(node)
            call segment_point_distance2(a_u(j), a_v(j), d_u(j), d_v(j), len2(j), pu, pv, dj, tj, f)
            if (dj < d2) then
               d2 = dj ; k = j ; t = tj ; feature = f
            endif
         enddo
      else
         c1 = left(node) ; c2 = right(node)
         b1 = box_d2_oac(lo_u(c1), lo_v(c1), hi_u(c1), hi_v(c1), pu, pv)
         b2 = box_d2_oac(lo_u(c2), lo_v(c2), hi_u(c2), hi_v(c2), pu, pv)
         if (b1 <= b2) then ! push the farther child first, visit the nearer one next
            top = top + 1_I4P ; stack(top) = c2
            top = top + 1_I4P ; stack(top) = c1
         else
            top = top + 1_I4P ; stack(top) = c1
            top = top + 1_I4P ; stack(top) = c2
         endif
      endif
   enddo
   endsubroutine segment_closest_oac

   pure subroutine closest_brute(self, pu, pv, d2, k, t, feature)
   !< Same as `closest`, by brute force over the payload in payload order (reference for tests).
   class(segment_bvh), intent(in)  :: self    !< Tree.
   real(R8P),          intent(in)  :: pu, pv  !< Query point.
   real(R8P),          intent(out) :: d2      !< Squared distance.
   integer(I4P),       intent(out) :: k       !< Payload index.
   real(R8P),          intent(out) :: t       !< Parameter.
   integer(I4P),       intent(out) :: feature !< FEATURE_*.
   integer(I4P)                    :: j, f    !< Work.
   real(R8P)                       :: dj, tj  !< Work.

   d2 = MaxR8P ; k = 0_I4P ; t = 0._R8P ; feature = FEATURE_NONE
   do j = 1, self%segments_number
      call segment_point_distance2(self%a_u(j), self%a_v(j), self%d_u(j), self%d_v(j), self%len2(j), pu, pv, dj, tj, f)
      if (dj < d2) then
         d2 = dj ; k = j ; t = tj ; feature = f
      endif
   enddo
   endsubroutine closest_brute

   elemental subroutine destroy(self)
   !< Reset to an empty tree.
   class(segment_bvh), intent(inout) :: self  !< Tree.
   type(segment_bvh)                 :: fresh !< Fresh instance.

   select type(self)
   type is(segment_bvh)
      self = fresh
   endselect
   endsubroutine destroy

   pure subroutine refit_vertices(self, u, v, next)
   !< Refit the tree to moved vertices (issue #41): segment seg_id(k) runs from vertex seg_id(k) to
   !< next(seg_id(k)) (the numbering `curve_object` builds with). The payload is recomputed in the same
   !< leaf order, the node boxes bottom-up (the builder's min/max over the leaf ranges, so the same
   !< values a rebuild gives its leaves), and the winding cache rebuilt.
   class(segment_bvh), intent(inout) :: self       !< Tree.
   real(R8P),          intent(in)    :: u(:), v(:) !< Moved vertices.
   integer(I4P),       intent(in)    :: next(:)    !< Next vertex.
   integer(I4P)                      :: k, s, node, c1, c2 !< Work.

   if (self%nodes_number == 0_I4P) return
   do k = 1, self%segments_number
      s = self%seg_id(k)
      self%a_u(k) = u(s) ; self%a_v(k) = v(s)
      self%d_u(k) = u(next(s)) - u(s) ; self%d_v(k) = v(next(s)) - v(s)
      self%len2(k) = self%d_u(k)**2 + self%d_v(k)**2
   enddo
   do node = self%nodes_number, 1, -1
      if (self%left(node) == 0_I4P) then
         self%lo_u(node) = huge(1._R8P) ; self%lo_v(node) = huge(1._R8P)
         self%hi_u(node) = -huge(1._R8P) ; self%hi_v(node) = -huge(1._R8P)
         do k = self%first(node), self%last(node)
            s = self%seg_id(k)
            self%lo_u(node) = min(self%lo_u(node), u(s), u(next(s))) ; self%hi_u(node) = max(self%hi_u(node), u(s), u(next(s)))
            self%lo_v(node) = min(self%lo_v(node), v(s), v(next(s))) ; self%hi_v(node) = max(self%hi_v(node), v(s), v(next(s)))
         enddo
      else
         c1 = self%left(node) ; c2 = self%right(node)
         self%lo_u(node) = min(self%lo_u(c1), self%lo_u(c2)) ; self%hi_u(node) = max(self%hi_u(c1), self%hi_u(c2))
         self%lo_v(node) = min(self%lo_v(c1), self%lo_v(c2)) ; self%hi_v(node) = max(self%hi_v(c1), self%hi_v(c2))
      endif
   enddo
   deallocate(self%wn_cu, self%wn_cv, self%wn_du, self%wn_dv, self%wn_c, self%wn_diag2)
   call build_winding_cache(self)
   endsubroutine refit_vertices

   pure subroutine segments_in_box(self, lo_u, lo_v, hi_u, hi_v, ids, n)
   !< Payload indices of the segments whose bounding box overlaps the box [lo, hi] (closed).
   !<
   !< `ids` is (re)allocated to the number found, n, in deterministic traversal order.
   class(segment_bvh),        intent(in)    :: self       !< Tree.
   real(R8P),                 intent(in)    :: lo_u, lo_v !< Box minimum.
   real(R8P),                 intent(in)    :: hi_u, hi_v !< Box maximum.
   integer(I4P), allocatable, intent(inout) :: ids(:)     !< Payload indices found.
   integer(I4P),              intent(out)   :: n          !< Number found.
   integer(I4P)                             :: stack(STACK_SIZE), top, node, j !< Work.
   integer(I4P), allocatable                :: buf(:)     !< Growing buffer.

   n = 0_I4P
   allocate(buf(16))
   if (self%nodes_number > 0_I4P) then
      top = 1_I4P ; stack(1) = 1_I4P
      do while (top > 0_I4P)
         node = stack(top) ; top = top - 1_I4P
         if (self%lo_u(node) > hi_u .or. self%hi_u(node) < lo_u .or. self%lo_v(node) > hi_v .or. self%hi_v(node) < lo_v) cycle
         if (self%left(node) == 0_I4P) then
            do j = self%first(node), self%last(node)
               if (min(self%a_u(j), self%a_u(j) + self%d_u(j)) > hi_u) cycle
               if (max(self%a_u(j), self%a_u(j) + self%d_u(j)) < lo_u) cycle
               if (min(self%a_v(j), self%a_v(j) + self%d_v(j)) > hi_v) cycle
               if (max(self%a_v(j), self%a_v(j) + self%d_v(j)) < lo_v) cycle
               if (n == size(buf)) buf = [buf, buf]
               n = n + 1_I4P ; buf(n) = j
            enddo
         else
            top = top + 1_I4P ; stack(top) = self%right(node)
            top = top + 1_I4P ; stack(top) = self%left(node)
         endif
      enddo
   endif
   if (allocated(ids)) deallocate(ids)
   allocate(ids(n))
   ids(1:n) = buf(1:n)
   endsubroutine segments_in_box

   pure function winding_number(self, pu, pv, beta) result(w)
   !< Fast generalized winding number of (pu, pv) (see the module notes); beta default 2.
   !< Runs `segment_winding_number_oac`, the device kernel.
   class(segment_bvh), intent(in)           :: self   !< Tree.
   real(R8P),          intent(in)           :: pu, pv !< Point.
   real(R8P),          intent(in), optional :: beta   !< Admissibility (default 2, as in 3-D).
   real(R8P)                                :: w      !< Winding number.
   real(R8P)                                :: b2     !< Squared admissibility.

   w = 0._R8P
   if (self%nodes_number == 0_I4P) return
   b2 = 4._R8P ; if (present(beta)) b2 = beta * beta
   w = segment_winding_number_oac(self%left, self%right, self%first, self%last, self%a_u, self%a_v, self%d_u, self%d_v, &
                                  self%wn_cu, self%wn_cv, self%wn_du, self%wn_dv, self%wn_c, self%wn_diag2, pu, pv, b2)
   endfunction winding_number

   pure function segment_winding_number_oac(left, right, first, last, a_u, a_v, d_u, d_v, wn_cu, wn_cv, wn_du, wn_dv, &
                                            wn_c, wn_diag2, pu, pv, b2) result(w)
   !< Fast generalized winding number on the flat tree arrays (non-empty tree). CPU and device.
   !$acc routine seq
   integer(I4P), intent(in) :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   real(R8P),    intent(in) :: a_u(:), a_v(:), d_u(:), d_v(:)       !< Segment payload.
   real(R8P),    intent(in) :: wn_cu(:), wn_cv(:), wn_du(:), wn_dv(:) !< Winding cache: centre, dipole.
   real(R8P),    intent(in) :: wn_c(:,:), wn_diag2(:)               !< Winding cache: second moment, diagonal.
   real(R8P),    intent(in) :: pu, pv                               !< Point.
   real(R8P),    intent(in) :: b2                                   !< Squared admissibility.
   real(R8P)                :: w                                    !< Winding number.
   real(R8P), parameter     :: TWO_PI = 2._R8P * 3.14159265358979323846264338327950288_R8P !< 2 pi.
   integer(I4P)             :: stack(STACK_SIZE), top, node, j       !< Work.
   real(R8P)                :: ru, rv, r2, omega, tr, rcr            !< Work.

   omega = 0._R8P
   top = 1_I4P ; stack(1) = 1_I4P
   do while (top > 0_I4P)
      node = stack(top) ; top = top - 1_I4P
      ru = wn_cu(node) - pu ; rv = wn_cv(node) - pv
      r2 = ru * ru + rv * rv
      if (r2 > b2 * wn_diag2(node) .and. wn_diag2(node) > 0._R8P) then
         ! n . r / r^2 + trace(C) / r^2 - 2 r^T C r / r^4
         tr  = wn_c(1, node) + wn_c(4, node)
         rcr = ru * (wn_c(1, node) * ru + wn_c(3, node) * rv) + rv * (wn_c(2, node) * ru + wn_c(4, node) * rv)
         omega = omega + (wn_du(node) * ru + wn_dv(node) * rv) / r2 + tr / r2 - 2._R8P * rcr / (r2 * r2)
      elseif (left(node) == 0_I4P) then
         do j = first(node), last(node)
            omega = omega + segment_winding(a_u(j), a_v(j), d_u(j), d_v(j), pu, pv)
         enddo
      else
         top = top + 1_I4P ; stack(top) = right(node)
         top = top + 1_I4P ; stack(top) = left(node)
      endif
   enddo
   w = omega / TWO_PI
   endfunction segment_winding_number_oac

   pure function winding_number_brute(self, pu, pv) result(w)
   !< Generalized winding number of (pu, pv): exact sum of the subtended angles / 2 pi (reference).
   class(segment_bvh), intent(in) :: self   !< Tree.
   real(R8P),          intent(in) :: pu, pv !< Point.
   real(R8P)                      :: w      !< Winding number.
   real(R8P), parameter :: TWO_PI = 2._R8P * 3.14159265358979323846264338327950288_R8P !< 2 pi.
   integer(I4P)                   :: j      !< Counter.

   w = 0._R8P
   do j = 1, self%segments_number
      w = w + segment_winding(self%a_u(j), self%a_v(j), self%d_u(j), self%d_v(j), pu, pv)
   enddo
   w = w / TWO_PI
   endfunction winding_number_brute

   pure function segment_winding(au, av, du, dv, pu, pv) result(theta)
   !< Signed angle subtended at p by the segment a -> a + d (positive counter-clockwise).
   !$acc routine seq
   real(R8P), intent(in) :: au, av, du, dv !< Segment.
   real(R8P), intent(in) :: pu, pv         !< Point.
   real(R8P)             :: theta          !< Angle in (-pi, pi].
   real(R8P)             :: xu, xv, yu, yv !< a - p, b - p.

   xu = au - pu ; xv = av - pv
   yu = au + du - pu ; yv = av + dv - pv
   theta = atan2(xu * yv - xv * yu, xu * yu + xv * yv)
   endfunction segment_winding

   pure subroutine segment_point_distance2(au, av, du, dv, len2, pu, pv, d2, t, feature)
   !< Squared distance from (pu, pv) to the segment a + t d, t in [0, 1], its parameter and feature.
   !<
   !< At a vertex the distance is computed from the vertex itself (not from a + t d), so the two
   !< segments sharing a vertex give the same bits there. In the interior, d2 = cross^2 / len2,
   !< accurate close to the segment.
   !$acc routine seq
   real(R8P),    intent(in)  :: au, av   !< Segment start.
   real(R8P),    intent(in)  :: du, dv   !< Segment direction.
   real(R8P),    intent(in)  :: len2     !< Squared length.
   real(R8P),    intent(in)  :: pu, pv   !< Point.
   real(R8P),    intent(out) :: d2       !< Squared distance.
   real(R8P),    intent(out) :: t        !< Parameter of the closest point.
   integer(I4P), intent(out) :: feature  !< FEATURE_*.
   real(R8P)                 :: wu, wv, dot, cr !< Work.

   wu = pu - au ; wv = pv - av
   dot = wu * du + wv * dv
   if (dot <= 0._R8P .or. len2 == 0._R8P) then
      d2 = wu * wu + wv * wv ; t = 0._R8P ; feature = FEATURE_VERTEX_START
   elseif (dot >= len2) then
      wu = pu - (au + du) ; wv = pv - (av + dv)
      d2 = wu * wu + wv * wv ; t = 1._R8P ; feature = FEATURE_VERTEX_END
   else
      cr = du * wv - dv * wu
      d2 = cr * cr / len2 ; t = dot / len2 ; feature = FEATURE_SEGMENT
   endif
   endsubroutine segment_point_distance2

   ! private
   pure function box_d2_oac(lo_u, lo_v, hi_u, hi_v, pu, pv) result(d2)
   !< Squared distance from (pu, pv) to the box [lo, hi] (0 inside). CPU and device.
   !$acc routine seq
   real(R8P), intent(in) :: lo_u, lo_v, hi_u, hi_v !< Box.
   real(R8P), intent(in) :: pu, pv                 !< Point.
   real(R8P)             :: d2                     !< Squared distance.
   real(R8P)             :: eu, ev                 !< Excess per axis.

   eu = max(lo_u - pu, 0._R8P, pu - hi_u)
   ev = max(lo_v - pv, 0._R8P, pv - hi_v)
   d2 = eu * eu + ev * ev
   endfunction box_d2_oac

   pure subroutine build_winding_cache(self)
   !< Per-node expansion centre (box centre), dipole and second moments, squared box diagonal.
   !<
   !< Accumulated bottom-up: children always have larger indices than their parent, so a reverse
   !< sweep over the nodes sees the children first.
   type(segment_bvh), intent(inout) :: self   !< Tree.
   integer(I4P)                     :: node, j, c, h !< Work.
   real(R8P)                        :: lnu, lnv, mu, mv !< length * normal, midpoint offset.

   allocate(self%wn_cu(self%nodes_number), self%wn_cv(self%nodes_number))
   allocate(self%wn_du(self%nodes_number), self%wn_dv(self%nodes_number))
   allocate(self%wn_c(4, self%nodes_number), self%wn_diag2(self%nodes_number))
   do node = self%nodes_number, 1, -1
      self%wn_cu(node) = 0.5_R8P * (self%lo_u(node) + self%hi_u(node))
      self%wn_cv(node) = 0.5_R8P * (self%lo_v(node) + self%hi_v(node))
      self%wn_diag2(node) = (self%hi_u(node) - self%lo_u(node))**2 + (self%hi_v(node) - self%lo_v(node))**2
      self%wn_du(node) = 0._R8P ; self%wn_dv(node) = 0._R8P ; self%wn_c(:, node) = 0._R8P
      if (self%left(node) == 0_I4P) then
         do j = self%first(node), self%last(node)
            lnu = self%d_v(j) ; lnv = -self%d_u(j) ! length * outward normal
            mu = self%a_u(j) + 0.5_R8P * self%d_u(j) - self%wn_cu(node)
            mv = self%a_v(j) + 0.5_R8P * self%d_v(j) - self%wn_cv(node)
            self%wn_du(node) = self%wn_du(node) + lnu
            self%wn_dv(node) = self%wn_dv(node) + lnv
            self%wn_c(:, node) = self%wn_c(:, node) + [lnu * mu, lnv * mu, lnu * mv, lnv * mv]
         enddo
      else
         do h = 1, 2
            c = self%left(node) ; if (h == 2) c = self%right(node)
            ! shift the child's moments to this centre: C += D_child (c_child - c)^T
            mu = self%wn_cu(c) - self%wn_cu(node) ; mv = self%wn_cv(c) - self%wn_cv(node)
            self%wn_du(node) = self%wn_du(node) + self%wn_du(c)
            self%wn_dv(node) = self%wn_dv(node) + self%wn_dv(c)
            self%wn_c(:, node) = self%wn_c(:, node) + self%wn_c(:, c) + &
                                 [self%wn_du(c) * mu, self%wn_dv(c) * mu, self%wn_du(c) * mv, self%wn_dv(c) * mv]
         enddo
      endif
   enddo
   endsubroutine build_winding_cache

   pure subroutine split(perm, cu, cv, sb, median_only, mid)
   !< Partition perm (segment ids) in place; return the size of the left part (1 <= mid < n).
   !<
   !< SAH over 16 centroid buckets per axis, cost = nL * hp(L) + nR * hp(R) with hp the box
   !< half-perimeter. Falls back to a median split by centroid on the widest axis when no
   !< bucket split separates the segments (coincident centroids) or when `median_only`.
   integer(I4P), intent(inout) :: perm(:)     !< Segment ids of the node.
   real(R8P),    intent(in)    :: cu(:), cv(:) !< Centroids of all segments.
   real(R8P),    intent(in)    :: sb(:,:)     !< Boxes of all segments [lo_u, lo_v, hi_u, hi_v; n].
   logical,      intent(in)    :: median_only !< Skip the SAH.
   integer(I4P), intent(out)   :: mid         !< Left part size.
   real(R8P)                   :: cmin(2), cmax(2), c, best_cost, cost !< Work.
   real(R8P)                   :: blo(2, BUCKETS), bhi(2, BUCKETS), llo(2), lhi(2), rlo(2), rhi(2) !< Bucket boxes.
   integer(I4P)                :: cnt(BUCKETS), n, axis, best_axis, best_b, b, i, j, nl, tmp !< Work.

   n = size(perm, kind=I4P)
   cmin = [minval(cu(perm)), minval(cv(perm))]
   cmax = [maxval(cu(perm)), maxval(cv(perm))]
   best_cost = MaxR8P ; best_axis = 0_I4P ; best_b = 0_I4P
   if (.not. median_only) then
      do axis = 1, 2
         if (cmax(axis) <= cmin(axis)) cycle
         cnt = 0_I4P
         blo = MaxR8P ; bhi = -MaxR8P
         do i = 1, n
            c = centroid(perm(i), axis)
            b = bucket_of(c, cmin(axis), cmax(axis))
            cnt(b) = cnt(b) + 1_I4P
            blo(:, b) = min(blo(:, b), sb(1:2, perm(i)))
            bhi(:, b) = max(bhi(:, b), sb(3:4, perm(i)))
         enddo
         do b = 1, BUCKETS - 1
            nl = sum(cnt(1:b))
            if (nl == 0_I4P .or. nl == n) cycle
            llo = MaxR8P ; lhi = -MaxR8P ; rlo = MaxR8P ; rhi = -MaxR8P
            do j = 1, b
               if (cnt(j) > 0_I4P) then
                  llo = min(llo, blo(:, j)) ; lhi = max(lhi, bhi(:, j))
               endif
            enddo
            do j = b + 1, BUCKETS
               if (cnt(j) > 0_I4P) then
                  rlo = min(rlo, blo(:, j)) ; rhi = max(rhi, bhi(:, j))
               endif
            enddo
            cost = nl * sum(lhi - llo) + (n - nl) * sum(rhi - rlo)
            if (cost < best_cost) then
               best_cost = cost ; best_axis = axis ; best_b = b
            endif
         enddo
      enddo
   endif
   if (best_axis > 0_I4P) then
      ! partition: bucket <= best_b to the left (stable two-pointer, deterministic)
      mid = 0_I4P
      do i = 1, n
         if (bucket_of(centroid(perm(i), best_axis), cmin(best_axis), cmax(best_axis)) <= best_b) then
            mid = mid + 1_I4P
            tmp = perm(mid) ; perm(mid) = perm(i) ; perm(i) = tmp
         endif
      enddo
   else
      axis = 1_I4P ; if (cmax(2) - cmin(2) > cmax(1) - cmin(1)) axis = 2_I4P
      call sort_by_centroid(perm, axis)
      mid = n / 2_I4P
   endif
   contains
      pure function centroid(s, ax) result(cc)
      !< Centroid coordinate of segment s on axis ax.
      integer(I4P), intent(in) :: s, ax !< Segment, axis.
      real(R8P)                :: cc    !< Coordinate.

      if (ax == 1_I4P) then
         cc = cu(s)
      else
         cc = cv(s)
      endif
      endfunction centroid

      pure subroutine sort_by_centroid(p, ax)
      !< Insertion sort of p by centroid on axis ax, ties by id (deterministic; nodes are small
      !< here or degenerate).
      integer(I4P), intent(inout) :: p(:) !< Ids.
      integer(I4P), intent(in)    :: ax   !< Axis.
      integer(I4P)                :: q, r, key !< Work.

      do q = 2, size(p)
         key = p(q) ; r = q - 1_I4P
         do while (r >= 1_I4P)
            if (centroid(p(r), ax) < centroid(key, ax) .or. &
                (centroid(p(r), ax) == centroid(key, ax) .and. p(r) < key)) exit
            p(r + 1) = p(r) ; r = r - 1_I4P
         enddo
         p(r + 1) = key
      enddo
      endsubroutine sort_by_centroid
   endsubroutine split

   pure function bucket_of(c, cmin, cmax) result(b)
   !< SAH bucket (1..BUCKETS) of the centroid coordinate c in [cmin, cmax].
   real(R8P), intent(in) :: c, cmin, cmax !< Coordinate and range.
   integer(I4P)          :: b             !< Bucket.

   b = 1_I4P + int(BUCKETS * (c - cmin) / (cmax - cmin), I4P)
   b = min(max(b, 1_I4P), BUCKETS)
   endfunction bucket_of
endmodule fossil_segment_bvh
