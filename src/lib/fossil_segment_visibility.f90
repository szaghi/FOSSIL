!< FOSSIL, segment visibility against the geometry (WP1 step 01): does a segment from a foot point cross the surface or
!< the curve? The any-hit kernels of `probe_many` / `image_points`, shared by the CPU and the device paths.
module fossil_segment_visibility
!< FOSSIL, segment visibility against the geometry (WP1 step 01): does a segment from a foot point cross the surface or
!< the curve? The any-hit kernels of `probe_many` / `image_points`, shared by the CPU and the device paths.
!<
!< 3-D, `bvh_segment_blocked_oac`: the open segment from a (on the surface) to b is blocked when a facet is crossed at a
!< parameter tol < t <= |b - a| of the ray from a along (b - a) / |b - a|, the semantics of issue #37 (the facets
!< containing the foot are crossed at t ~ 0 and ignored through tol). Every decision is the exact ray-facet crossing
!< of issue #82 (`fossil_ray_kernel`, the very routines of the host `intersect_ray_*`), on the ray prepared against the
!< root box of the SAH BVH as `intersect_ray_all` prepares it, so the result equals the host's former
!< `intersect_ray_all` + filter bit for bit. The traversal prunes a node when the ray misses its box padded by the
!< ray's pad (the slab test of the host traversal, the same arithmetic) or when every point of the padded box projects
!< beyond |b - a| + pad on the ray (a crossing point is a convex combination of the facet's vertices, so it lies in
!< the box up to round-off, far below the pad, and its t is its projection): neither prunes a facet the filter would
!< keep. It stops at the first blocking facet.
!<
!< 2-D, `curve_segment_blocked_oac`: the closed segment a-b is blocked when it shares a point (exact `orient2d_xy`,
!< `segments_intersect_2d`) with a segment of the curve other than those containing the foot point a (the closest
!< segment sg when vx = 0, else the two segments at the closest vertex vx); candidates are the segments whose payload
!< box meets the box of a-b, as the host `segments_in_box` lists them.
!<
!< Free `!$acc routine seq` procedures on flat arrays (#20 rules: no type-bound calls, no `class` dummies, no
!< vector_R8P, scalar arguments).
use fossil_aabb_node_object, only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, bvh_bmax_x_oac, &
                                    bvh_bmax_y_oac, bvh_bmax_z_oac, bvh_left_child_oac, bvh_right_child_oac, &
                                    bvh_payload_first_oac, bvh_payload_count_oac
use fossil_aabb_tree_object, only : facet_distance_payload
use fossil_box_predicates,   only : segments_intersect_2d
use fossil_ray_kernel,       only : ray_prepare_oac, ray_crossing_oac
use penf,                    only : I4P, R8P

implicit none
private
public :: bvh_segment_blocked_oac, curve_segment_blocked_oac

integer(I4P), parameter :: STACK = 256_I4P !< Traversal stack (> BVH depth cap 127 + 2; 2-D cap 128).

contains
   pure function bvh_segment_blocked_oac(node, payload, tv, ax, ay, az, bx, by, bz, tol) result(blocked)
   !< The open segment from a (the foot point) to b crosses a facet of the SAH BVH (0-based nodes, root 0; exact
   !< vertices `tv(9, facets)` indexed by facet id) at tol < t <= |b - a| (see the module header).
   !$acc routine seq
   type(aabb_node_object),       intent(in) :: node(0:)       !< Tree nodes.
   type(facet_distance_payload), intent(in) :: payload(:)     !< Payload (facet ids).
   real(R8P),                    intent(in) :: tv(:,:)        !< Exact vertices [9, facets].
   real(R8P),                    intent(in) :: ax, ay, az     !< Foot point a.
   real(R8P),                    intent(in) :: bx, by, bz     !< End point b.
   real(R8P),                    intent(in) :: tol            !< Crossings with t <= tol (the facets at a) are ignored.
   logical                                  :: blocked        !< A facet is crossed in (tol, |b - a|].
   integer(I4P)                             :: stk(STACK), sp, cur, lc, rc, k, f, axis !< Work.
   real(R8P)                                :: o(3), d(3), q(3), lo(3), hi(3), ex, ey, ez, len, pad, dd, tlim !< Ray.
   real(R8P)                                :: tmin, tmax, t1, t2, tmp, ood, tp, t !< Slab test, crossing.
   real(R8P), parameter                     :: PTOL = 4._R8P * tiny(1._R8P) !< Parallel-ray cutoff of the padded slab test.
   logical                                  :: hit, miss      !< Crossing, slab test.

   blocked = .false.
   if (size(node) == 0 .or. size(payload) == 0) return
   ex = bx - ax ; ey = by - ay ; ez = bz - az
   len = sqrt(ex * ex + ey * ey + ez * ez)
   if (len <= tol) return
   o(1) = ax ; o(2) = ay ; o(3) = az
   d(1) = ex / len ; d(2) = ey / len ; d(3) = ez / len
   lo(1) = bvh_bmin_x_oac(node(0)) ; lo(2) = bvh_bmin_y_oac(node(0)) ; lo(3) = bvh_bmin_z_oac(node(0))
   hi(1) = bvh_bmax_x_oac(node(0)) ; hi(2) = bvh_bmax_y_oac(node(0)) ; hi(3) = bvh_bmax_z_oac(node(0))
   call ray_prepare_oac(o, d, lo, hi, 0._R8P, q, axis, pad)
   dd = ((0._R8P + d(1) * d(1)) + d(2) * d(2)) + d(3) * d(3)
   tlim = len + pad
   sp = 1_I4P ; stk(1) = 0_I4P
   do while (sp > 0_I4P)
      cur = stk(sp) ; sp = sp - 1_I4P
      lo(1) = bvh_bmin_x_oac(node(cur)) - pad ; hi(1) = bvh_bmax_x_oac(node(cur)) + pad
      lo(2) = bvh_bmin_y_oac(node(cur)) - pad ; hi(2) = bvh_bmax_y_oac(node(cur)) + pad
      lo(3) = bvh_bmin_z_oac(node(cur)) - pad ; hi(3) = bvh_bmax_z_oac(node(cur)) + pad
      ! slab test of the host traversal (`aabb_object%ray_slab_interval` with a pad), the same arithmetic
      tmin = 0._R8P ; tmax = huge(1._R8P) ; miss = .false.
      do k = 1, 3
         if (abs(d(k)) < PTOL) then
            if (o(k) < lo(k) .or. o(k) > hi(k)) then
               miss = .true. ; exit
            endif
         else
            ood = 1._R8P / d(k)
            t1 = (lo(k) - o(k)) * ood
            t2 = (hi(k) - o(k)) * ood
            if (t1 > t2) then
               tmp = t1 ; t1 = t2 ; t2 = tmp
            endif
            if (t1 > tmin) tmin = t1
            if (t2 < tmax) tmax = t2
            if (tmin > tmax) then
               miss = .true. ; exit
            endif
         endif
      enddo
      if (miss) cycle
      ! least projection on the ray of a point of the padded box: beyond |b - a| + pad, no crossing there is kept
      tp = 0._R8P
      do k = 1, 3
         if (d(k) >= 0._R8P) then
            tp = tp + (lo(k) - o(k)) * d(k)
         else
            tp = tp + (hi(k) - o(k)) * d(k)
         endif
      enddo
      if (tp > tlim * dd) cycle
      do k = 0_I4P, bvh_payload_count_oac(node(cur)) - 1_I4P
         f = payload(bvh_payload_first_oac(node(cur)) + k)%facet_id
         call ray_crossing_oac(o, q, d, axis, tv(:, f), hit, t)
         if (hit .and. t > tol .and. t <= len) then
            blocked = .true.
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
   endfunction bvh_segment_blocked_oac

   pure function curve_segment_blocked_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, a_u, a_v, d_u, d_v, seg_id, &
                                           next, u, v, au, av, bu, bv, sg, vx) result(blocked)
   !< The closed segment a-b (a the foot point) shares a point with a segment of the curve (segment BVH, root 1) other
   !< than the closest segment sg (vx = 0) or the two segments at the closest vertex vx (see the module header).
   !$acc routine seq
   real(R8P),    intent(in) :: lo_u(:), lo_v(:), hi_u(:), hi_v(:)  !< Node boxes.
   integer(I4P), intent(in) :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   real(R8P),    intent(in) :: a_u(:), a_v(:), d_u(:), d_v(:)      !< Payload start points and directions.
   integer(I4P), intent(in) :: seg_id(:)    !< Curve segment of each payload entry.
   integer(I4P), intent(in) :: next(:)      !< Next vertex of each vertex (segment i = i -> next(i)).
   real(R8P),    intent(in) :: u(:), v(:)   !< Curve vertices.
   real(R8P),    intent(in) :: au, av       !< Foot point a.
   real(R8P),    intent(in) :: bu, bv       !< End point b.
   integer(I4P), intent(in) :: sg           !< Closest segment (its start vertex).
   integer(I4P), intent(in) :: vx           !< Closest vertex (0: segment interior).
   logical                  :: blocked      !< The segment meets another segment.
   integer(I4P)             :: stk(STACK), sp, nd, j, q1, q2 !< Work.
   real(R8P)                :: bl_u, bl_v, bh_u, bh_v !< Box of a-b.

   blocked = .false.
   if (size(lo_u) == 0) return
   bl_u = min(au, bu) ; bl_v = min(av, bv) ; bh_u = max(au, bu) ; bh_v = max(av, bv)
   sp = 1_I4P ; stk(1) = 1_I4P
   do while (sp > 0_I4P)
      nd = stk(sp) ; sp = sp - 1_I4P
      if (lo_u(nd) > bh_u .or. hi_u(nd) < bl_u .or. lo_v(nd) > bh_v .or. hi_v(nd) < bl_v) cycle
      if (left(nd) == 0_I4P) then
         do j = first(nd), last(nd)
            if (min(a_u(j), a_u(j) + d_u(j)) > bh_u) cycle
            if (max(a_u(j), a_u(j) + d_u(j)) < bl_u) cycle
            if (min(a_v(j), a_v(j) + d_v(j)) > bh_v) cycle
            if (max(a_v(j), a_v(j) + d_v(j)) < bl_v) cycle
            q1 = seg_id(j) ; q2 = next(q1)
            if (vx == 0_I4P .and. q1 == sg) cycle
            if (vx > 0_I4P .and. (q1 == vx .or. q2 == vx)) cycle
            if (segments_intersect_2d(au, av, bu, bv, u(q1), v(q1), u(q2), v(q2))) then
               blocked = .true.
               return
            endif
         enddo
      else
         sp = sp + 1_I4P ; stk(sp) = right(nd)
         sp = sp + 1_I4P ; stk(sp) = left(nd)
      endif
   enddo
   endfunction curve_segment_blocked_oac
endmodule fossil_segment_visibility
