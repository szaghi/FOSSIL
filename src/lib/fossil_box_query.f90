!< FOSSIL, box queries on the trees for AMR tagging (issue #38).
module fossil_box_query
!< FOSSIL, box queries on the trees for AMR tagging (issue #38).
!<
!< Classification codes of a box against the geometry:
!<
!< - `BOX_CUT`: the box meets the geometry (3-D: conservative by the tolerance; 2-D: exact);
!< - `BOX_NEAR`: not cut, but the geometry comes within the radius (3-D: within the box inflated by
!<   the radius, so up to sqrt(3) radius at the corners; 2-D: segment-box distance <= radius, the
!<   distance in floating point, so exact up to its round-off);
!< - `BOX_INSIDE`, `BOX_OUTSIDE`: neither, the box centre is inside or outside the solid.
!<
!< The traversals are free `_oac` procedures (`!$acc routine seq`) on the flat tree arrays, shared
!< by the CPU and the device paths (#20 rules: no type-bound calls, no `class` dummies). They prune
!< a node when its box does not meet the query box inflated by the reach, and stop at the first cut.
!<
!< Tag bands (issue #38 step 4) on a block of cells: `cell_bands` is the Chebyshev cell distance to
!< the nearest cut cell (breadth-first search over the 26 neighbours, capped), `cell_components`
!< labels the 26-connected components of the cells that are not cut. Two closed cells that share a
!< point and are both uncut lie on the same side of a closed surface (a path between their centres
!< through the shared point stays in the two cells, which the surface does not meet), so one point
!< query per component signs it. A 2-D block is a block with one cell in z.
use fossil_aabb_node_object, only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, bvh_bmax_x_oac, &
                                    bvh_bmax_y_oac, bvh_bmax_z_oac, bvh_left_child_oac, bvh_right_child_oac, &
                                    bvh_payload_first_oac, bvh_payload_count_oac
use fossil_aabb_tree_object, only : facet_distance_payload
use fossil_box_predicates,   only : box_triangle_overlap, box_segment_overlap_2d, box_segment_distance2_2d
use penf,                    only : I4P, R8P

implicit none
private
public :: BOX_NONE, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
public :: bvh_box_contact_oac, curve_box_contact_oac, box_tolerance, bvh_box_thin_oac
public :: cell_bands, cell_components

integer(I4P), parameter :: BOX_NONE    = 0_I4P !< Contact result: neither cut nor near (inside/outside not decided).
integer(I4P), parameter :: BOX_CUT     = 1_I4P !< The box meets the geometry.
integer(I4P), parameter :: BOX_NEAR    = 2_I4P !< Not cut, the geometry within the radius.
integer(I4P), parameter :: BOX_INSIDE  = 3_I4P !< Neither, box inside the solid.
integer(I4P), parameter :: BOX_OUTSIDE = 4_I4P !< Neither, box outside the solid.

integer(I4P), parameter :: STACK = 256_I4P !< Traversal stack (> BVH depth cap 127 + 2; 2-D cap 128).

contains
   pure function box_tolerance(bmin, bmax) result(tol)
   !< Conservative inflation of the 3-D box-triangle test for a geometry with bounding box
   !< [bmin, bmax]: 1e-12 of the larger of the box diagonal and the largest coordinate magnitude
   !< (the round-off of a separating-axis projection scales with both).
   real(R8P), intent(in) :: bmin(3), bmax(3) !< Geometry bounding box.
   real(R8P)             :: tol              !< Tolerance.

   tol = 1.e-12_R8P * max(norm2(bmax - bmin), maxval(abs(bmin)), maxval(abs(bmax)))
   endfunction box_tolerance

   pure function bvh_box_contact_oac(node, payload, lo, hi, tol, reach) result(code)
   !< Contact of the box [lo, hi] with the triangles of a SAH BVH (0-based nodes, root 0).
   !<
   !< BOX_CUT if a triangle meets the box inflated by `tol`; else BOX_NEAR if one meets the box
   !< inflated by `reach` (>= tol); else BOX_NONE. Nodes whose box misses the box inflated by
   !< `reach` are pruned; the traversal stops at the first cut.
   !$acc routine seq
   type(aabb_node_object),       intent(in) :: node(0:)   !< Tree nodes.
   type(facet_distance_payload), intent(in) :: payload(:) !< Triangle payload.
   real(R8P),                    intent(in) :: lo(3), hi(3) !< Query box.
   real(R8P),                    intent(in) :: tol        !< Cut inflation.
   real(R8P),                    intent(in) :: reach      !< Near inflation (>= tol).
   integer(I4P)                             :: code       !< BOX_CUT, BOX_NEAR or BOX_NONE.
   integer(I4P)                             :: stk(STACK), sp, cur, lc, rc, k, idx !< Work.
   real(R8P)                                :: v1(3), v2(3), v3(3) !< Triangle.

   code = BOX_NONE
   sp = 1_I4P ; stk(1) = 0_I4P
   do while (sp > 0_I4P)
      cur = stk(sp) ; sp = sp - 1_I4P
      if (bvh_bmin_x_oac(node(cur)) > hi(1) + reach .or. bvh_bmax_x_oac(node(cur)) < lo(1) - reach .or. &
          bvh_bmin_y_oac(node(cur)) > hi(2) + reach .or. bvh_bmax_y_oac(node(cur)) < lo(2) - reach .or. &
          bvh_bmin_z_oac(node(cur)) > hi(3) + reach .or. bvh_bmax_z_oac(node(cur)) < lo(3) - reach) cycle
      do k = 0_I4P, bvh_payload_count_oac(node(cur)) - 1_I4P
         idx = bvh_payload_first_oac(node(cur)) + k
         v1 = payload(idx)%v1 ; v2 = v1 + payload(idx)%e12 ; v3 = v1 + payload(idx)%e13
         if (code == BOX_NONE) then
            if (.not. box_triangle_overlap(lo, hi, v1, v2, v3, reach)) cycle
            code = BOX_NEAR
         endif
         if (box_triangle_overlap(lo, hi, v1, v2, v3, tol)) then
            code = BOX_CUT
            return
         endif
      enddo
      lc = bvh_left_child_oac(node(cur)) ; rc = bvh_right_child_oac(node(cur))
      if (rc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = rc
      endif
      if (lc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = lc
      endif
   enddo
   endfunction bvh_box_contact_oac

   pure function bvh_box_thin_oac(node, payload, thin_pay, thin_node, lo, hi, tol, factor) result(touch)
   !< True if the box [lo, hi] (inflated by `tol`, conservative) meets a triangle of the SAH BVH whose local size (min of
   !< thickness and gap, issue #62) is below factor * h, h the largest side of the box. `thin_pay(idx)` is the local size
   !< of payload triangle idx, `thin_node(cur)` the minimum over the subtree of node cur: a subtree with no triangle that
   !< thin is pruned whole. Stops at the first such triangle. The CPU and the device run this same routine.
   !$acc routine seq
   type(aabb_node_object),       intent(in) :: node(0:)     !< Tree nodes.
   type(facet_distance_payload), intent(in) :: payload(:)   !< Triangle payload.
   real(R8P),                    intent(in) :: thin_pay(:)  !< Local size per payload triangle.
   real(R8P),                    intent(in) :: thin_node(0:) !< Minimum local size per node.
   real(R8P),                    intent(in) :: lo(3), hi(3) !< Query box.
   real(R8P),                    intent(in) :: tol          !< Inflation.
   real(R8P),                    intent(in) :: factor       !< Threshold factor.
   logical                                  :: touch        !< Contact with a thin triangle.
   integer(I4P)                             :: stk(STACK), sp, cur, lc, rc, k, idx !< Work.
   real(R8P)                                :: v1(3), v2(3), v3(3), thr !< Triangle, threshold.

   touch = .false.
   thr = factor * max(hi(1) - lo(1), hi(2) - lo(2), hi(3) - lo(3))
   sp = 1_I4P ; stk(1) = 0_I4P
   do while (sp > 0_I4P)
      cur = stk(sp) ; sp = sp - 1_I4P
      if (thin_node(cur) >= thr) cycle
      if (bvh_bmin_x_oac(node(cur)) > hi(1) + tol .or. bvh_bmax_x_oac(node(cur)) < lo(1) - tol .or. &
          bvh_bmin_y_oac(node(cur)) > hi(2) + tol .or. bvh_bmax_y_oac(node(cur)) < lo(2) - tol .or. &
          bvh_bmin_z_oac(node(cur)) > hi(3) + tol .or. bvh_bmax_z_oac(node(cur)) < lo(3) - tol) cycle
      do k = 0_I4P, bvh_payload_count_oac(node(cur)) - 1_I4P
         idx = bvh_payload_first_oac(node(cur)) + k
         if (thin_pay(idx) >= thr) cycle
         v1 = payload(idx)%v1 ; v2 = v1 + payload(idx)%e12 ; v3 = v1 + payload(idx)%e13
         if (box_triangle_overlap(lo, hi, v1, v2, v3, tol)) then
            touch = .true.
            return
         endif
      enddo
      lc = bvh_left_child_oac(node(cur)) ; rc = bvh_right_child_oac(node(cur))
      if (rc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = rc
      endif
      if (lc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = lc
      endif
   enddo
   endfunction bvh_box_thin_oac

   pure function curve_box_contact_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, lo, hi, &
                                       radius) result(code)
   !< Contact of the box [lo, hi] (u, v) with the segments of a curve on its segment BVH (root 1).
   !<
   !< BOX_CUT if a segment meets the closed box (exact); else BOX_NEAR if one is within `radius`
   !< (floating-point segment-box distance <= radius); else BOX_NONE. Node pruning widens the
   !< radius by a few ulps, so round-off never prunes a segment the distance test would accept. The segments are taken from the curve's own
   !< vertices (`u`, `v`, `next`), not from the payload's start + direction, so no rounding enters.
   !$acc routine seq
   real(R8P),    intent(in) :: lo_u(:), lo_v(:), hi_u(:), hi_v(:)  !< Node boxes.
   integer(I4P), intent(in) :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   integer(I4P), intent(in) :: seg_id(:)    !< Curve segment of each payload entry.
   integer(I4P), intent(in) :: next(:)      !< Next vertex of each vertex (segment i = i -> next(i)).
   real(R8P),    intent(in) :: u(:), v(:)   !< Curve vertices.
   real(R8P),    intent(in) :: lo(2), hi(2) !< Query box.
   real(R8P),    intent(in) :: radius       !< Near radius (>= 0).
   integer(I4P)             :: code         !< BOX_CUT, BOX_NEAR or BOX_NONE.
   integer(I4P)             :: stk(STACK), sp, nd, j, s !< Work.
   real(R8P)                :: a(2), b(2), r2, mu, mv !< Work.

   code = BOX_NONE
   r2 = radius * radius
   mu = 0._R8P ; mv = 0._R8P ! node pruning margin (node boxes are exact; with radius 0 the test is exact)
   if (radius > 0._R8P) then
      mu = radius + 4._R8P * spacing(max(abs(lo(1)), abs(hi(1)), radius))
      mv = radius + 4._R8P * spacing(max(abs(lo(2)), abs(hi(2)), radius))
   endif
   sp = 1_I4P ; stk(1) = 1_I4P
   do while (sp > 0_I4P)
      nd = stk(sp) ; sp = sp - 1_I4P
      if (lo_u(nd) > hi(1) + mu .or. hi_u(nd) < lo(1) - mu .or. lo_v(nd) > hi(2) + mv .or. hi_v(nd) < lo(2) - mv) cycle
      if (left(nd) == 0_I4P) then
         do j = first(nd), last(nd)
            s = seg_id(j)
            a = [u(s), v(s)] ; b = [u(next(s)), v(next(s))]
            if (box_segment_overlap_2d(lo, hi, a, b)) then
               code = BOX_CUT
               return
            endif
            if (code == BOX_NONE .and. radius > 0._R8P) then
               if (box_segment_distance2_2d(lo, hi, a, b) <= r2) code = BOX_NEAR
            endif
         enddo
      else
         sp = sp + 1_I4P ; stk(sp) = right(nd)
         sp = sp + 1_I4P ; stk(sp) = left(nd)
      endif
   enddo
   endfunction curve_box_contact_oac

   pure subroutine cell_bands(cut, band, dist)
   !< Chebyshev cell distance to the nearest cut cell, capped at band + 1 (also where no cell is cut).
   logical,      intent(in)  :: cut(:,:,:)  !< Cut cells.
   integer(I4P), intent(in)  :: band        !< Band width (>= 0).
   integer(I4P), intent(out) :: dist(:,:,:) !< Distance [same shape], 0 on cut cells.
   integer(I4P), allocatable :: q(:,:)      !< Queue of cells [3, cells].
   integer(I4P)              :: n(3), head, tail, i, j, k, di, dj, dk, a, b, c, d !< Work.

   n = shape(cut)
   dist = band + 1_I4P
   allocate(q(3, max(1, size(cut))))
   tail = 0_I4P
   do k = 1, n(3)
      do j = 1, n(2)
         do i = 1, n(1)
            if (cut(i, j, k)) then
               dist(i, j, k) = 0_I4P
               tail = tail + 1_I4P ; q(:, tail) = [i, j, k]
            endif
         enddo
      enddo
   enddo
   head = 0_I4P
   do while (head < tail)
      head = head + 1_I4P
      i = q(1, head) ; j = q(2, head) ; k = q(3, head)
      d = dist(i, j, k) + 1_I4P
      if (d > band) cycle
      do dk = -1, 1
         c = k + dk ; if (c < 1 .or. c > n(3)) cycle
         do dj = -1, 1
            b = j + dj ; if (b < 1 .or. b > n(2)) cycle
            do di = -1, 1
               a = i + di ; if (a < 1 .or. a > n(1)) cycle
               if (dist(a, b, c) > d) then
                  dist(a, b, c) = d
                  tail = tail + 1_I4P ; q(:, tail) = [a, b, c]
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine cell_bands

   pure subroutine cell_components(cut, comp, ncomp, seed)
   !< 26-connected components of the cells that are not cut: comp = 0 on cut cells, 1..ncomp else;
   !< seed(:, m) is the first cell (in storage order) of component m.
   logical,      intent(in)               :: cut(:,:,:)  !< Cut cells.
   integer(I4P), intent(out)              :: comp(:,:,:) !< Component labels [same shape].
   integer(I4P), intent(out)              :: ncomp       !< Number of components.
   integer(I4P), allocatable, intent(out) :: seed(:,:)   !< First cell of each component [3, ncomp].
   integer(I4P), allocatable              :: q(:,:), sd(:,:) !< Queue, seeds buffer.
   integer(I4P)                           :: n(3), head, tail, i, j, k, ii, jj, kk, di, dj, dk, a, b, c !< Work.

   n = shape(cut)
   comp = 0_I4P ; ncomp = 0_I4P
   allocate(q(3, max(1, size(cut))), sd(3, max(1, size(cut))))
   do k = 1, n(3)
      do j = 1, n(2)
         do i = 1, n(1)
            if (cut(i, j, k) .or. comp(i, j, k) > 0_I4P) cycle
            ncomp = ncomp + 1_I4P ; sd(:, ncomp) = [i, j, k]
            comp(i, j, k) = ncomp
            head = 0_I4P ; tail = 1_I4P ; q(:, 1) = [i, j, k]
            do while (head < tail)
               head = head + 1_I4P
               ii = q(1, head) ; jj = q(2, head) ; kk = q(3, head)
               do dk = -1, 1
                  c = kk + dk ; if (c < 1 .or. c > n(3)) cycle
                  do dj = -1, 1
                     b = jj + dj ; if (b < 1 .or. b > n(2)) cycle
                     do di = -1, 1
                        a = ii + di ; if (a < 1 .or. a > n(1)) cycle
                        if (cut(a, b, c) .or. comp(a, b, c) > 0_I4P) cycle
                        comp(a, b, c) = ncomp
                        tail = tail + 1_I4P ; q(:, tail) = [a, b, c]
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   enddo
   allocate(seed(3, ncomp))
   seed = sd(:, 1:ncomp)
   endsubroutine cell_components
endmodule fossil_box_query
