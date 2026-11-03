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
use fossil_aabb_node_object, only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, bvh_bmax_x_oac, &
                                    bvh_bmax_y_oac, bvh_bmax_z_oac, bvh_left_child_oac, bvh_right_child_oac, &
                                    bvh_payload_first_oac, bvh_payload_count_oac
use fossil_aabb_tree_object, only : facet_distance_payload
use fossil_box_predicates,   only : box_triangle_overlap, box_segment_overlap_2d, box_segment_distance2_2d
use penf,                    only : I4P, R8P

implicit none
private
public :: BOX_NONE, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
public :: bvh_box_contact_oac, curve_box_contact_oac, box_tolerance

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
endmodule fossil_box_query
