!< FOSSIL, exact crossings of axis-aligned lines with triangles and segments (issue #40).
module fossil_line_query
!< FOSSIL, exact crossings of axis-aligned lines with triangles and segments (issue #40).
!<
!< A line along axis k (1, 2, 3) is given by its two other coordinates (a, b) in cyclic order: for
!< k = 1 the line is y = a, z = b; for k = 2, z = a, x = b; for k = 3, x = a, y = b. It crosses a
!< triangle exactly when the point (a, b) lies in the triangle projected on the plane of those two
!< coordinates. That is decided with the exact `orient2d_xy`, and a point on a projected edge or
!< vertex is resolved by one global symbolic perturbation of the grid (issue #69), x + eps1,
!< y + eps2, z + eps3 with eps1 >> eps2 >> eps3: (a, b) + (eps, eps^2) on x and z lines, (a, b) +
!< (eps^2, eps) on y lines ((a, b) = (z, x)). The sign of orient2d(A, B, P + delta) is the sign of
!< orient2d(A, B, P) if it is not zero, else of its derivative along the first perturbed coordinate
!< (A_b - B_b for a, B_a - A_a for b), else along the second (comparisons, so exact). Before #69
!< each axis perturbed its own pair cyclically (y >> z, z >> x, x >> y: no global order), so a
!< point on the surface could be inside along one axis and outside along another. The three perturbed signs agree exactly when the
!< line crosses the triangle, so on a watertight surface every crossing is counted once (Woop et al.
!< 2013, made exact); a triangle whose projection is degenerate (parallel to the line) is never
!< crossed. The common sign is the sign of the triangle normal along the axis: positive, the line
!< (going in +k) leaves the solid; negative, it enters.
!<
!< 2-D: a line along u (v = a) crosses the segment pq when p and q lie on opposite sides of
!< v = a + eps (a vertex with v = a counts as below); along v (u = a), likewise in u. With the solid on
!< the left of every segment, a line in +u leaves the solid through an upward segment and enters
!< through a downward one; a line in +v enters through a segment pointing in +u.
!<
!< The crossing parameter (the coordinate along the axis) is constructed in floating point from the
!< barycentric weights: the topology is exact, the positions are rounded (decision 1 of #40).
!<
!< The traversals are free `_oac` procedures (`!$acc routine seq`) on flat arrays, shared by the CPU
!< and the device paths: they gather the crossings with the parameter in [t0, t1] into fixed buffers
!< (`overflow` if the buffer is too small) and sort them by parameter, then by id.
use fossil_aabb_node_object, only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, bvh_bmax_x_oac, &
                                    bvh_bmax_y_oac, bvh_bmax_z_oac, bvh_left_child_oac, bvh_right_child_oac, &
                                    bvh_payload_first_oac, bvh_payload_count_oac
use fossil_aabb_tree_object, only : facet_distance_payload
use fossil_predicates,       only : orient2d_xy, orient3d_xyz
use penf,                    only : I4P, R8P

implicit none
private
public :: line_crossing_before, line_crossing_before_2d
public :: tri_line_hit_oac, seg_line_hit_oac, bvh_line_crossings_oac, curve_line_crossings_oac, sort_crossings_oac

integer(I4P), parameter :: STACK = 256_I4P !< Traversal stack (> BVH depth cap + 2).

contains
   function line_crossing_before(tri, axis, a, b, c, entering) result(before)
   !< The crossing of triangle tri with the line along `axis` at (a, b) lies before the point of the line at coordinate
   !< c, exactly (issue #69): host only. With n = (B - A) x (C - A) and P the point, orient3d(A, B, C, P) =
   !< -n_axis (c - t), t the exact crossing parameter, and the crossing enters (n_axis < 0) or leaves (n_axis > 0):
   !< the sign of c - t is one exact orient3d. A point on the triangle's plane takes the global perturbation P +
   !< (eps1, eps2, eps3) (issue #69, as the line kernel): orient3d(A, B, C, P + delta) = -n . delta, the sign of
   !< minus the first nonzero component of n in the order x, y, z.
   real(R8P),    intent(in) :: tri(3, 3) !< Triangle vertices (columns), as the crossing kernel saw them.
   integer(I4P), intent(in) :: axis     !< Line axis.
   real(R8P),    intent(in) :: a, b     !< Line position (cyclic order, see the module header).
   real(R8P),    intent(in) :: c        !< Point coordinate along the axis.
   logical,      intent(in) :: entering !< The crossing enters the solid.
   logical                  :: before   !< Crossing before the (perturbed) point.
   real(R8P)                :: p(3), o

   select case (axis)
   case (1_I4P) ; p = [c, a, b]
   case (2_I4P) ; p = [b, c, a]
   case default ; p = [a, b, c]
   endselect
   o = orient3d_xyz(tri(:, 1), tri(:, 2), tri(:, 3), p)
   if (o == 0._R8P) o = -first_normal_sign(tri)
   if (o == 0._R8P) then
      before = .true. ! degenerate triangle: never crossed by the kernel
   elseif (entering) then
      before = o > 0._R8P
   else
      before = o < 0._R8P
   endif
   endfunction line_crossing_before

   function first_normal_sign(tri) result(sg)
   !< Sign (-1, 0, +1) of the first nonzero component of n = (B - A) x (C - A) in the order x, y, z, exactly (the
   !< components are the orientations of the projections on (y, z), (z, x), (x, y)).
   real(R8P), intent(in) :: tri(3, 3) !< Triangle vertices (columns).
   real(R8P)             :: sg
   integer(I4P)          :: k, i1, i2

   sg = 0._R8P
   do k = 1, 3
      i1 = mod(k, 3) + 1 ; i2 = mod(k + 1, 3) + 1
      sg = orient2d_xy(tri(i1, 1), tri(i2, 1), tri(i1, 2), tri(i2, 2), tri(i1, 3), tri(i2, 3))
      if (sg /= 0._R8P) then
         sg = sign(1._R8P, sg) ; return
      endif
   enddo
   endfunction first_normal_sign

   function line_crossing_before_2d(pu, pv, qu, qv, axis, a, c, entering) result(before)
   !< 2-D twin of `line_crossing_before`: the crossing of the segment p -> q with the line along `axis` (1: v = a, 2:
   !< u = a) lies before the point of the line at coordinate c, exactly. For P on the line, orient2d(p, q, P) =
   !< -(q_v - p_v) (c - t) along u and (q_u - p_u) (c - t) along v, and with the module's entering convention (the
   !< solid on the left) both give: entering, before iff orient2d > 0; leaving, iff < 0. A point on the segment takes
   !< the global perturbation P + (eps1, eps2) (u first, as the 2-D node test): orient2d(p, q, P + delta) has the sign
   !< of -(q_v - p_v), else of q_u - p_u.
   real(R8P),    intent(in) :: pu, pv, qu, qv !< Segment.
   integer(I4P), intent(in) :: axis           !< Line axis.
   real(R8P),    intent(in) :: a              !< Line position.
   real(R8P),    intent(in) :: c              !< Point coordinate along the axis.
   logical,      intent(in) :: entering       !< The crossing enters the solid.
   logical                  :: before         !< Crossing before the (perturbed) point.
   real(R8P)                :: o

   if (axis == 1_I4P) then
      o = orient2d_xy(pu, pv, qu, qv, c, a)
   else
      o = orient2d_xy(pu, pv, qu, qv, a, c)
   endif
   if (o == 0._R8P) then
      o = -(qv - pv) ; if (o == 0._R8P) o = qu - pu
   endif
   if (o == 0._R8P) then
      before = .true. ! degenerate segment
   elseif (entering) then
      before = o > 0._R8P
   else
      before = o < 0._R8P
   endif
   endfunction line_crossing_before_2d

   pure subroutine tri_line_hit_oac(tv, axis, a, b, hit, t, enter)
   !< Crossing of the line along `axis` at (a, b) with the triangle tv = [A, B, C] (9 coordinates).
   !$acc routine seq
   real(R8P),    intent(in)  :: tv(9)  !< Triangle vertices (x, y, z of A, then B, then C).
   integer(I4P), intent(in)  :: axis   !< Line axis (1, 2, 3).
   real(R8P),    intent(in)  :: a, b   !< Line position (the other two coordinates, cyclic).
   logical,      intent(out) :: hit    !< The line crosses the triangle.
   real(R8P),    intent(out) :: t      !< Crossing coordinate along the axis.
   logical,      intent(out) :: enter  !< The line (in +axis) enters the solid there.
   integer(I4P)              :: i1, i2 !< Projected coordinate indices.
   real(R8P)                 :: wa, wb, wc, sa, sb, sc !< Weights and perturbed signs.
   logical                   :: bfirst !< Perturb b before a.

   i1 = mod(axis, 3_I4P) + 1_I4P ; i2 = mod(axis + 1_I4P, 3_I4P) + 1_I4P
   hit = .false. ; t = 0._R8P ; enter = .false.
   bfirst = axis == 2_I4P ! global perturbation x + eps1, y + eps2, z + eps3: on a y line (a, b) = (z, x), x first
   call psign(tv(3 + i1), tv(3 + i2), tv(6 + i1), tv(6 + i2), a, b, bfirst, wa, sa) ! edge B -> C: weight of A
   call psign(tv(6 + i1), tv(6 + i2), tv(i1), tv(i2), a, b, bfirst, wb, sb)         ! edge C -> A: weight of B
   call psign(tv(i1), tv(i2), tv(3 + i1), tv(3 + i2), a, b, bfirst, wc, sc)         ! edge A -> B: weight of C
   if (sa == 0._R8P .or. sa /= sb .or. sa /= sc) return
   hit = .true.
   enter = sa < 0._R8P
   if (wa + wb + wc /= 0._R8P) then
      t = (wa * tv(axis) + wb * tv(3 + axis) + wc * tv(6 + axis)) / (wa + wb + wc)
   else
      t = (tv(axis) + tv(3 + axis) + tv(6 + axis)) / 3._R8P
   endif
   endsubroutine tri_line_hit_oac

   pure subroutine psign(p1, p2, q1, q2, a, b, bfirst, w, s)
   !< Weight orient2d(p, q, (a, b)) (approximate magnitude, exact sign) and its perturbed sign: (a, b) + (eps, eps^2),
   !< or (a, b) + (eps^2, eps) with `bfirst` (issue #69: the two orders make one global perturbation of the grid).
   !$acc routine seq
   real(R8P), intent(in)  :: p1, p2, q1, q2 !< Edge end points (projected).
   real(R8P), intent(in)  :: a, b           !< Line position.
   logical,   intent(in)  :: bfirst         !< Perturb b first.
   real(R8P), intent(out) :: w, s           !< Weight, perturbed sign (-1, 0, +1).
   real(R8P)              :: da, db         !< d/da = p2 - q2, d/db = q1 - p1 of orient2d(p, q, P).

   w = orient2d_xy(p1, p2, q1, q2, a, b)
   da = 0._R8P ; db = 0._R8P
   if (p2 > q2) da = 1._R8P
   if (p2 < q2) da = -1._R8P
   if (q1 > p1) db = 1._R8P
   if (q1 < p1) db = -1._R8P
   if (w > 0._R8P) then
      s = 1._R8P
   elseif (w < 0._R8P) then
      s = -1._R8P
   elseif (bfirst) then
      s = db ; if (s == 0._R8P) s = da
   else
      s = da ; if (s == 0._R8P) s = db
   endif
   endsubroutine psign

   pure subroutine seg_line_hit_oac(pu, pv, qu, qv, axis, a, hit, t, enter)
   !< Crossing of the line along `axis` (1: along u at v = a; 2: along v at u = a) with the segment pq.
   !$acc routine seq
   real(R8P),    intent(in)  :: pu, pv, qu, qv !< Segment end points.
   integer(I4P), intent(in)  :: axis           !< 1 (u) or 2 (v).
   real(R8P),    intent(in)  :: a              !< Line position.
   logical,      intent(out) :: hit            !< The line crosses the segment.
   real(R8P),    intent(out) :: t              !< Crossing coordinate along the line.
   logical,      intent(out) :: enter          !< The line (in +axis) enters the solid there.
   logical                   :: above_p, above_q !< Sides of v = a + eps (u = a + eps).

   hit = .false. ; t = 0._R8P ; enter = .false.
   if (axis == 1_I4P) then
      above_p = pv > a ; above_q = qv > a
      if (above_p .eqv. above_q) return
      hit = .true. ; enter = qv < pv
      t = pu + (a - pv) * (qu - pu) / (qv - pv)
   else
      above_p = pu > a ; above_q = qu > a
      if (above_p .eqv. above_q) return
      hit = .true. ; enter = qu > pu
      t = pv + (a - pu) * (qv - pv) / (qu - pu)
   endif
   endsubroutine seg_line_hit_oac

   pure subroutine bvh_line_crossings_oac(node, payload, tv, axis, a, b, t0, t1, ts, ids, es, n, overflow)
   !< Crossings with parameter in [t0, t1] of the line along `axis` at (a, b) with the triangles of a
   !< SAH BVH (0-based nodes, root 0; exact vertices `tv(9, facets)` indexed by facet id), sorted.
   !$acc routine seq
   type(aabb_node_object),       intent(in)  :: node(0:)   !< Tree nodes.
   type(facet_distance_payload), intent(in)  :: payload(:) !< Payload (facet ids).
   real(R8P),                    intent(in)  :: tv(:,:)    !< Exact vertices [9, facets].
   integer(I4P),                 intent(in)  :: axis       !< Line axis.
   real(R8P),                    intent(in)  :: a, b       !< Line position.
   real(R8P),                    intent(in)  :: t0, t1     !< Parameter range.
   real(R8P),                    intent(out) :: ts(:)      !< Crossing parameters [nmax].
   integer(I4P),                 intent(out) :: ids(:)     !< Facets [nmax].
   logical,                      intent(out) :: es(:)      !< Entering flags [nmax].
   integer(I4P),                 intent(out) :: n          !< Crossings found.
   logical,                      intent(out) :: overflow   !< More crossings than the buffer holds.
   integer(I4P)                              :: stk(STACK), sp, cur, lc, rc, k, f, i1, i2 !< Work.
   real(R8P)                                 :: lo(3), hi(3), t !< Node box, crossing.
   logical                                   :: hit, enter      !< Crossing.

   n = 0_I4P ; overflow = .false.
   i1 = mod(axis, 3_I4P) + 1_I4P ; i2 = mod(axis + 1_I4P, 3_I4P) + 1_I4P
   sp = 1_I4P ; stk(1) = 0_I4P
   do while (sp > 0_I4P)
      cur = stk(sp) ; sp = sp - 1_I4P
      lo = [bvh_bmin_x_oac(node(cur)), bvh_bmin_y_oac(node(cur)), bvh_bmin_z_oac(node(cur))]
      hi = [bvh_bmax_x_oac(node(cur)), bvh_bmax_y_oac(node(cur)), bvh_bmax_z_oac(node(cur))]
      if (a < lo(i1) .or. a > hi(i1) .or. b < lo(i2) .or. b > hi(i2) .or. hi(axis) < t0 .or. lo(axis) > t1) cycle
      do k = 0_I4P, bvh_payload_count_oac(node(cur)) - 1_I4P
         f = payload(bvh_payload_first_oac(node(cur)) + k)%facet_id
         call tri_line_hit_oac(tv(:, f), axis, a, b, hit, t, enter)
         if (.not. hit .or. t < t0 .or. t > t1) cycle
         if (n == size(ts)) then
            overflow = .true.
            cycle
         endif
         n = n + 1_I4P ; ts(n) = t ; ids(n) = f ; es(n) = enter
      enddo
      lc = bvh_left_child_oac(node(cur)) ; rc = bvh_right_child_oac(node(cur))
      if (rc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = rc
      endif
      if (lc > 0_I4P) then
         sp = sp + 1_I4P ; stk(sp) = lc
      endif
   enddo
   call sort_crossings_oac(ts, ids, es, n)
   endsubroutine bvh_line_crossings_oac

   pure subroutine curve_line_crossings_oac(lo_u, lo_v, hi_u, hi_v, left, right, first, last, seg_id, next, u, v, axis, a, &
                                            t0, t1, ts, ids, es, n, overflow)
   !< Crossings with parameter in [t0, t1] of the line along `axis` (1: u, 2: v) at a with the segments
   !< of a curve (segment BVH, root 1; segments from the curve's own vertices), sorted.
   !$acc routine seq
   real(R8P),    intent(in)  :: lo_u(:), lo_v(:), hi_u(:), hi_v(:)  !< Node boxes.
   integer(I4P), intent(in)  :: left(:), right(:), first(:), last(:) !< Children, leaf ranges.
   integer(I4P), intent(in)  :: seg_id(:), next(:) !< Payload -> segment, next vertex.
   real(R8P),    intent(in)  :: u(:), v(:)         !< Vertices.
   integer(I4P), intent(in)  :: axis               !< 1 or 2.
   real(R8P),    intent(in)  :: a, t0, t1          !< Line position, parameter range.
   real(R8P),    intent(out) :: ts(:)              !< Crossing parameters [nmax].
   integer(I4P), intent(out) :: ids(:)             !< Segments [nmax].
   logical,      intent(out) :: es(:)              !< Entering flags [nmax].
   integer(I4P), intent(out) :: n                  !< Crossings found.
   logical,      intent(out) :: overflow           !< More crossings than the buffer holds.
   integer(I4P)              :: stk(STACK), sp, nd, j, s !< Work.
   real(R8P)                 :: t                  !< Crossing.
   logical                   :: hit, enter, skip   !< Crossing, pruning.

   n = 0_I4P ; overflow = .false.
   sp = 1_I4P ; stk(1) = 1_I4P
   do while (sp > 0_I4P)
      nd = stk(sp) ; sp = sp - 1_I4P
      if (axis == 1_I4P) then
         skip = a < lo_v(nd) .or. a > hi_v(nd) .or. hi_u(nd) < t0 .or. lo_u(nd) > t1
      else
         skip = a < lo_u(nd) .or. a > hi_u(nd) .or. hi_v(nd) < t0 .or. lo_v(nd) > t1
      endif
      if (skip) cycle
      if (left(nd) == 0_I4P) then
         do j = first(nd), last(nd)
            s = seg_id(j)
            call seg_line_hit_oac(u(s), v(s), u(next(s)), v(next(s)), axis, a, hit, t, enter)
            if (.not. hit .or. t < t0 .or. t > t1) cycle
            if (n == size(ts)) then
               overflow = .true.
               cycle
            endif
            n = n + 1_I4P ; ts(n) = t ; ids(n) = s ; es(n) = enter
         enddo
      else
         sp = sp + 1_I4P ; stk(sp) = right(nd)
         sp = sp + 1_I4P ; stk(sp) = left(nd)
      endif
   enddo
   call sort_crossings_oac(ts, ids, es, n)
   endsubroutine curve_line_crossings_oac

   pure subroutine sort_crossings_oac(ts, ids, es, n)
   !< Insertion sort of the first n crossings by parameter, then by id (deterministic).
   !$acc routine seq
   real(R8P),    intent(inout) :: ts(:)  !< Parameters.
   integer(I4P), intent(inout) :: ids(:) !< Ids.
   logical,      intent(inout) :: es(:)  !< Entering flags.
   integer(I4P), intent(in)    :: n      !< Count.
   integer(I4P)                :: i, j, ki !< Work.
   real(R8P)                   :: kt     !< Work.
   logical                     :: ke     !< Work.

   do i = 2_I4P, n
      kt = ts(i) ; ki = ids(i) ; ke = es(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (ts(j) < kt .or. (ts(j) == kt .and. ids(j) <= ki)) exit
         ts(j + 1) = ts(j) ; ids(j + 1) = ids(j) ; es(j + 1) = es(j)
         j = j - 1_I4P
      enddo
      ts(j + 1) = kt ; ids(j + 1) = ki ; es(j + 1) = ke
   enddo
   endsubroutine sort_crossings_oac
endmodule fossil_line_query
