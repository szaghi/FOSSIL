!< FOSSIL, a bounding-volume hierarchy over 3-D segments (the feature edges, issue #61 step 2).
module fossil_feature_bvh
!< FOSSIL, a bounding-volume hierarchy over 3-D segments (the feature edges, issue #61 step 2).
!<
!< Plain arrays, so the traversals run on the device as they are (issue #20 rules; no vector_R8P here):
!< - seg(6, n): the segments, endpoints (x1, y1, z1, x2, y2, z2), reordered by the build;
!< - box(6, nn): node boxes (lo, hi); link(4, nn): left child, right child (0 for a leaf), first segment, count.
!< Built top-down by median split of the segment midpoints on the longest axis of their box, at most LEAF segments per
!< leaf. `segment_box_overlap` is the separating-axis test of a segment and an axis-aligned box inflated by `tol`
!< (conservative: it never misses a contact); `feature_box_touch_oac` traverses the tree with it.
use penf, only : I4P, R8P

implicit none
private
public :: build_segment_bvh, segment_box_overlap, feature_box_touch_oac, feature_point_distance_oac

integer(I4P), parameter :: LEAF  = 4_I4P   !< Segments per leaf (at most).
integer(I4P), parameter :: STACK = 128_I4P !< Traversal stack (depth of a median split of 2^60 segments).

contains
   subroutine build_segment_bvh(seg, box, link)
   !< Build the tree over seg(6, n), reordering seg; box(6, nn), link(4, nn) allocated to the nodes used.
   real(R8P),                 intent(inout) :: seg(:,:)   !< Segments [6, n] (reordered).
   real(R8P),    allocatable, intent(out)   :: box(:,:)   !< Node boxes [6, nn].
   integer(I4P), allocatable, intent(out)   :: link(:,:)  !< Node links [4, nn].
   real(R8P),    allocatable                :: b(:,:)     !< Work boxes.
   integer(I4P), allocatable                :: l(:,:)     !< Work links.
   integer(I4P), allocatable                :: todo(:)    !< Nodes to split.
   integer(I4P)                             :: n, nn, top, node, first, cnt, axis, mid, k !< Counters.

   n = size(seg, 2)
   allocate(b(6, max(1, 2 * n)), l(4, max(1, 2 * n)), todo(max(1, 2 * n)))
   nn = 0
   if (n == 0) then
      allocate(box(6, 0), link(4, 0)) ; return
   endif
   nn = 1 ; l(:, 1) = [0, 0, 1, n] ; top = 1 ; todo(1) = 1
   do while (top > 0)
      node = todo(top) ; top = top - 1
      first = l(3, node) ; cnt = l(4, node)
      do k = 1, 3
         b(k, node) = minval(min(seg(k, first:first + cnt - 1), seg(k + 3, first:first + cnt - 1)))
         b(k + 3, node) = maxval(max(seg(k, first:first + cnt - 1), seg(k + 3, first:first + cnt - 1)))
      enddo
      if (cnt <= LEAF) cycle
      axis = maxloc(b(4:6, node) - b(1:3, node), dim=1)
      mid = cnt / 2
      call select_by_midpoint(seg(:, first:first + cnt - 1), axis, mid)
      l(1, node) = nn + 1 ; l(2, node) = nn + 2
      l(:, nn + 1) = [0, 0, first, mid] ; l(:, nn + 2) = [0, 0, first + mid, cnt - mid]
      top = top + 1 ; todo(top) = nn + 1 ; top = top + 1 ; todo(top) = nn + 2
      nn = nn + 2
   enddo
   box = b(:, 1:nn) ; link = l(:, 1:nn)
   endsubroutine build_segment_bvh

   subroutine select_by_midpoint(s, axis, k)
   !< Partial sort (quickselect, deterministic middle pivot) of the columns of s by their midpoint on `axis`, so the
   !< first k columns have midpoints <= those of the others.
   real(R8P),    intent(inout) :: s(:,:) !< Segments [6, m].
   integer(I4P), intent(in)    :: axis   !< Axis (1..3).
   integer(I4P), intent(in)    :: k      !< Split position.
   real(R8P)                   :: pivot, t(6) !< Work.
   integer(I4P)                :: lo, hi, i, j !< Counters.

   lo = 1 ; hi = size(s, 2)
   do while (lo < hi)
      pivot = key((lo + hi) / 2)
      i = lo ; j = hi
      do while (i <= j)
         do while (key(i) < pivot)
            i = i + 1
         enddo
         do while (key(j) > pivot)
            j = j - 1
         enddo
         if (i <= j) then
            t = s(:, i) ; s(:, i) = s(:, j) ; s(:, j) = t
            i = i + 1 ; j = j - 1
         endif
      enddo
      if (k <= j) then
         hi = j
      elseif (k >= i) then
         lo = i
      else
         exit
      endif
   enddo
   contains
      pure function key(c) result(m)
      !< Midpoint of column c on `axis`.
      integer(I4P), intent(in) :: c
      real(R8P)                :: m

      m = s(axis, c) + s(axis + 3, c)
      endfunction key
   endsubroutine select_by_midpoint

   pure function segment_box_overlap(lo, hi, p, q, tol) result(overlap)
   !$acc routine seq
   !< Separating-axis test of the segment [p, q] and the box [lo - tol, hi + tol]: the 3 box axes and the 3 cross products
   !< of the segment direction with them. Conservative for tol >= the round-off of the projections (box_tolerance).
   real(R8P), intent(in) :: lo(3), hi(3) !< Box.
   real(R8P), intent(in) :: p(3), q(3)   !< Segment.
   real(R8P), intent(in) :: tol          !< Inflation.
   logical               :: overlap      !< The segment meets the inflated box.
   real(R8P)             :: c(3), h(3), d(3), m(3), a(3), r, s !< Work.
   integer               :: k            !< Counter.

   overlap = .false.
   do k = 1, 3
      if (max(p(k), q(k)) < lo(k) - tol .or. min(p(k), q(k)) > hi(k) + tol) return
   enddo
   c = 0.5_R8P * (lo + hi) ; h = 0.5_R8P * (hi - lo) + tol
   d = q - p ; m = 0.5_R8P * (p + q) - c
   do k = 1, 3
      if (k == 1) a = [0._R8P, d(3), -d(2)]
      if (k == 2) a = [-d(3), 0._R8P, d(1)]
      if (k == 3) a = [d(2), -d(1), 0._R8P]
      r = h(1) * abs(a(1)) + h(2) * abs(a(2)) + h(3) * abs(a(3))
      s = a(1) * m(1) + a(2) * m(2) + a(3) * m(3)
      if (abs(s) > r + tol * (abs(a(1)) + abs(a(2)) + abs(a(3)))) return
   enddo
   overlap = .true.
   endfunction segment_box_overlap

   pure function feature_box_touch_oac(seg, box, link, lo, hi, tol) result(touch)
   !$acc routine seq
   !< True if a segment of the tree meets the box [lo, hi] inflated by tol (stops at the first). The CPU and the device
   !< run this same routine.
   real(R8P),    intent(in) :: seg(:,:)     !< Segments [6, n].
   real(R8P),    intent(in) :: box(:,:)     !< Node boxes [6, nn].
   integer(I4P), intent(in) :: link(:,:)    !< Node links [4, nn].
   real(R8P),    intent(in) :: lo(3), hi(3) !< Query box.
   real(R8P),    intent(in) :: tol          !< Inflation.
   logical                  :: touch        !< Contact.
   integer(I4P)             :: stk(STACK), top, node, k !< Traversal.

   touch = .false.
   if (size(link, 2) == 0) return
   top = 1 ; stk(1) = 1
   do while (top > 0)
      node = stk(top) ; top = top - 1
      if (box(1, node) > hi(1) + tol .or. box(4, node) < lo(1) - tol .or. &
          box(2, node) > hi(2) + tol .or. box(5, node) < lo(2) - tol .or. &
          box(3, node) > hi(3) + tol .or. box(6, node) < lo(3) - tol) cycle
      if (link(1, node) == 0) then
         do k = link(3, node), link(3, node) + link(4, node) - 1
            if (segment_box_overlap(lo, hi, seg(1:3, k), seg(4:6, k), tol)) then
               touch = .true. ; return
            endif
         enddo
      else
         top = top + 1 ; stk(top) = link(2, node)
         top = top + 1 ; stk(top) = link(1, node)
      endif
   enddo
   endfunction feature_box_touch_oac

   pure function feature_point_distance_oac(seg, box, link, p) result(d)
   !$acc routine seq
   !< Distance from p to the nearest segment of the tree (huge if the tree is empty): depth-first, nodes pruned when
   !< their box is farther than the best so far. Exact point-segment distance (clamped projection); the CPU and the
   !< device run this same routine, so their results are bit-identical.
   real(R8P),    intent(in) :: seg(:,:)  !< Segments [6, n].
   real(R8P),    intent(in) :: box(:,:)  !< Node boxes [6, nn].
   integer(I4P), intent(in) :: link(:,:) !< Node links [4, nn].
   real(R8P),    intent(in) :: p(3)      !< Point.
   real(R8P)                :: d         !< Distance.
   real(R8P)                :: best, e(3), w(3), t, dd, g(3) !< Work.
   integer(I4P)             :: stk(STACK), top, node, k      !< Traversal.

   best = huge(1._R8P)
   if (size(link, 2) > 0) then
      top = 1 ; stk(1) = 1
      do while (top > 0)
         node = stk(top) ; top = top - 1
         g = max(box(1:3, node) - p, 0._R8P) + max(p - box(4:6, node), 0._R8P)
         if (g(1) * g(1) + g(2) * g(2) + g(3) * g(3) >= best) cycle
         if (link(1, node) == 0) then
            do k = link(3, node), link(3, node) + link(4, node) - 1
               e = seg(4:6, k) - seg(1:3, k) ; w = p - seg(1:3, k)
               t = e(1) * e(1) + e(2) * e(2) + e(3) * e(3)
               if (t > 0._R8P) then
                  t = max(0._R8P, min(1._R8P, (w(1) * e(1) + w(2) * e(2) + w(3) * e(3)) / t))
               else
                  t = 0._R8P
               endif
               w = w - t * e
               dd = w(1) * w(1) + w(2) * w(2) + w(3) * w(3)
               if (dd < best) best = dd
            enddo
         else
            top = top + 1 ; stk(top) = link(2, node)
            top = top + 1 ; stk(top) = link(1, node)
         endif
      enddo
   endif
   d = huge(1._R8P)
   if (best < huge(1._R8P)) d = sqrt(best)
   endfunction feature_point_distance_oac
endmodule fossil_feature_bvh
