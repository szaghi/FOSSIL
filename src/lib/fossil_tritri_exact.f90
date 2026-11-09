!< FOSSIL, triangle-triangle intersection with exact predicates and canonical points (issue #33).
module fossil_tritri_exact
!< FOSSIL, triangle-triangle intersection with exact predicates and canonical points (issue #33).
!<
!< Two families. `pair_cut` / `coplanar_cut` (issue #65, used by the boolean arrangement) are exact with no snap and
!< return implicit points (`fossil_implicit_points`): see their comments. `pair_classify` / `coplanar_portion` (issue
!< #33, still used by `fossil_alpha_wrap`) snap and return rounded points, as described below.
!< `pair_cut` needs two exactly non-degenerate triangles (`facet_drop` /= 0); `degenerate_pair_meet` (issue #96, moved
!< here from `fossil_body_distance` by issue #99) decides, exactly, whether a pair with a degenerate triangle meets.
!<
!< `pair_classify` decides how two triangles meet from the exact signs of Shewchuk's `orient3d`
!< (each triangle's vertices against the other's plane), never from floating-point distances:
!<  - PAIR_SEPARATED: no common point (all signs of one triangle strictly on one side, or the two
!<    sections of the planes' line do not overlap);
!<  - PAIR_CROSSING: they share a segment of positive length (proper crossing, a vertex on the
!<    other plane, an edge lying in the other plane);
!<  - PAIR_TOUCHING: they share a single point;
!<  - PAIR_COPLANAR: all six signs are zero (handled by the caller's 2D pass).
!<
!< Snap: a vertex closer to the other plane than `SNAP_REL` (1e-10) times the longest edge of the
!< pair counts as lying on it (sign 0), and a pair whose six vertices all do is coplanar. Inputs
!< meant to be coplanar are rarely exactly so once rotated or scaled in floating point; without
!< the snap, faces coplanar by construction would be classified as near-coplanar crossings with
!< sliver segments. Outside that band the decisions are the exact orient3d signs.
!<
!< Exactness: the decisions are exact up to the snap; the intersection points are floating-point constructions,
!< but **canonical** ones. A point is produced only as
!<  - an input vertex (a zero sign), copied bit for bit;
!<  - an edge of one triangle crossing the plane of the other (`edge_plane_point`), or, when that
!<    edge is coplanar with an edge of the other triangle (exact `orient3d == 0`), the crossing of
!<    the two edges (`edge_edge_point`), or that triangle's vertex when two of its edges qualify.
!< Each formula orders its inputs lexicographically first, so every triangle pair that produces
!< the same geometric point produces the same bits: the facets of both surfaces that share a cut
!< endpoint get identical vertices, which an exact-equality vertex weld then joins.
use fossil_implicit_points, only : implicit_point_t, ip_explicit, ip_lpi, ip_eei, ip_lex_compare, ip_rank
use fossil_predicates,      only : orient3d, orient3d_xyz, orient2d_xy
use penf,                   only : I4P, R8P
use vecfor,                 only : vector_R8P

implicit none
private
public :: PAIR_SEPARATED, PAIR_CROSSING, PAIR_TOUCHING, PAIR_COPLANAR
public :: pair_classify
public :: coplanar_portion
public :: edge_plane_point, edge_edge_point
public :: lex_less
public :: pair_cut, coplanar_cut, facet_drop
public :: degenerate_pair_meet

integer(I4P), parameter :: PAIR_SEPARATED = 0_I4P !< No common point.
integer(I4P), parameter :: PAIR_CROSSING  = 1_I4P !< Common segment of positive length.
integer(I4P), parameter :: PAIR_TOUCHING  = 2_I4P !< Single common point.
integer(I4P), parameter :: PAIR_COPLANAR  = 3_I4P !< Same plane (all six orient3d signs zero, after the snap).
real(R8P),    parameter :: SNAP_REL       = 1.e-10_R8P !< On-plane snap, relative to the pair's longest edge.

contains
   subroutine pair_classify(ta, tb, kind, p, q)
   !< Classify the pair of triangles `ta`, `tb` (3 vertices each) and, for PAIR_CROSSING, return the
   !< common segment (p, q); for PAIR_TOUCHING, p = q = the common point.
   type(vector_R8P), intent(in)  :: ta(3) !< First triangle.
   type(vector_R8P), intent(in)  :: tb(3) !< Second triangle.
   integer(I4P),     intent(out) :: kind  !< PAIR_* code.
   type(vector_R8P), intent(out) :: p, q  !< Common segment (crossing) or point (touching).
   integer(I4P)                  :: sa(3), sb(3)           !< Signs of ta's vertices vs tb's plane, and vice versa.
   type(vector_R8P)              :: xa(2), xb(2)           !< Sections of each triangle with the other's plane.
   integer(I4P)                  :: na, nb                 !< Points in each section.
   type(vector_R8P)              :: d                      !< Direction of the planes' line (ordering only).
   real(R8P)                     :: sxa(2), sxb(2)         !< Parameters along d.
   integer(I4P)                  :: ia_lo, ia_hi, ib_lo, ib_hi, i
   real(R8P)                     :: tol, area2_a, area2_b  !< Snap distance; twice the triangle areas.

   kind = PAIR_SEPARATED
   p = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; q = p
   tol = SNAP_REL * max(longest_edge(ta), longest_edge(tb))
   area2_a = norm(cross(ta(2) - ta(1), ta(3) - ta(1)))
   area2_b = norm(cross(tb(2) - tb(1), tb(3) - tb(1)))
   do i = 1, 3
      sb(i) = snapped_sign(orient3d(ta(1), ta(2), ta(3), tb(i)), area2_a, tol)
      sa(i) = snapped_sign(orient3d(tb(1), tb(2), tb(3), ta(i)), area2_b, tol)
   enddo
   if (all(sb == 0_I4P)) then
      kind = PAIR_COPLANAR
      return
   endif
   if (all(sb == 1_I4P) .or. all(sb == -1_I4P)) return
   if (all(sa == 1_I4P) .or. all(sa == -1_I4P)) return

   call plane_section(t=ta, s=sa, other=tb, x=xa, n=na)
   call plane_section(t=tb, s=sb, other=ta, x=xb, n=nb)
   if (na == 0_I4P .or. nb == 0_I4P) return

   d = cross(ta(2) - ta(1), ta(3) - ta(1)) ! na x nb: only the ordering along the line matters
   d = cross(d, cross(tb(2) - tb(1), tb(3) - tb(1)))
   do i = 1, na ; sxa(i) = dot(xa(i), d) ; enddo
   do i = 1, nb ; sxb(i) = dot(xb(i), d) ; enddo
   ia_lo = 1 ; ia_hi = 1 ; ib_lo = 1 ; ib_hi = 1
   if (na == 2_I4P) then
      if (sxa(2) < sxa(1)) then ; ia_lo = 2 ; else ; ia_hi = 2 ; endif
   endif
   if (nb == 2_I4P) then
      if (sxb(2) < sxb(1)) then ; ib_lo = 2 ; else ; ib_hi = 2 ; endif
   endif
   ! overlap of [a_lo, a_hi] and [b_lo, b_hi] along the line; ties take ta's point
   if (sxa(ia_lo) >= sxb(ib_lo)) then
      p = xa(ia_lo)
      if (sxa(ia_lo) > sxb(ib_hi)) return
   else
      p = xb(ib_lo)
      if (sxb(ib_lo) > sxa(ia_hi)) return
   endif
   if (sxa(ia_hi) <= sxb(ib_hi)) then
      q = xa(ia_hi)
   else
      q = xb(ib_hi)
   endif
   ! vertex merging: a constructed endpoint within the snap of an input vertex is that vertex
   p = snap_to_vertices(p, ta, tb, tol)
   q = snap_to_vertices(q, ta, tb, tol)
   if (same_point(p, q)) then
      kind = PAIR_TOUCHING
      q = p
   else
      kind = PAIR_CROSSING
   endif
   endsubroutine pair_classify

   subroutine plane_section(t, s, other, x, n)
   !< Points of triangle `t` on the plane of `other`, given the exact signs `s` of t's vertices:
   !< vertices with a zero sign (copied), and the crossings of edges with opposite strict signs.
   type(vector_R8P), intent(in)  :: t(3)     !< Triangle.
   integer(I4P),     intent(in)  :: s(3)     !< Signs of t's vertices against other's plane.
   type(vector_R8P), intent(in)  :: other(3) !< Triangle whose plane cuts t.
   type(vector_R8P), intent(out) :: x(2)     !< Section points (at most 2: 3 zeros is coplanar).
   integer(I4P),     intent(out) :: n        !< Number of points.
   integer(I4P)                  :: i, j

   n = 0_I4P
   do i = 1, 3
      if (s(i) == 0_I4P .and. n < 2_I4P) then
         n = n + 1_I4P ; x(n) = t(i)
      endif
   enddo
   do i = 1, 3
      j = mod(i, 3) + 1
      if (s(i) * s(j) < 0_I4P .and. n < 2_I4P) then
         n = n + 1_I4P ; x(n) = crossing_point(u=t(i), v=t(j), other=other)
      endif
   enddo
   endsubroutine plane_section

   function crossing_point(u, v, other) result(x)
   !< The point where edge (u, v) crosses the plane of triangle `other` (u, v strictly on opposite
   !< sides). If the edge is coplanar with an edge of `other` (exact test), the point is on that
   !< edge's line: the canonical edge-edge crossing, or `other`'s vertex if two of its edges qualify.
   type(vector_R8P), intent(in) :: u, v     !< Edge endpoints.
   type(vector_R8P), intent(in) :: other(3) !< Triangle.
   type(vector_R8P)             :: x        !< Crossing point.
   integer(I4P)                 :: k, l, hits, last_k
   real(R8P)                    :: tol, area2 !< Snap distance; twice the reference area.

   hits = 0_I4P ; last_k = 0_I4P
   tol = SNAP_REL * max(longest_edge(other), norm(v - u))
   do k = 1, 3
      l = mod(k, 3) + 1
      ! coplanar lines (k, l) and (u, v): distance of v from the plane (k, l, u) within the snap
      area2 = max(norm(cross(other(l) - other(k), u - other(k))), norm(cross(other(l) - other(k), v - other(k))))
      if (snapped_sign(orient3d(other(k), other(l), u, v), area2, tol) == 0_I4P) then
         hits = hits + 1_I4P
         if (hits == 2_I4P) then
            ! two edges of `other` coplanar with (u, v): the point is their common vertex
            if (mod(last_k, 3) + 1 == k) then
               x = other(k)
            else
               x = other(last_k)
            endif
            return
         endif
         last_k = k
      endif
   enddo
   if (hits == 1_I4P) then
      x = edge_edge_point(u, v, other(last_k), other(mod(last_k, 3) + 1))
   else
      x = edge_plane_point(u, v, other)
   endif
   endfunction crossing_point

   pure function edge_plane_point(u, v, t) result(x)
   !< Canonical crossing of edge (u, v) with the plane of triangle t: endpoints and t's vertices are
   !< ordered lexicographically first, so the result does not depend on their order.
   type(vector_R8P), intent(in) :: u, v  !< Edge endpoints (on opposite sides of the plane).
   type(vector_R8P), intent(in) :: t(3)  !< Triangle.
   type(vector_R8P)             :: x     !< Crossing point.
   type(vector_R8P)             :: a, b, p1, p2, p3, n
   real(R8P)                    :: da, db

   call sorted2(u, v, a, b)
   call sorted3(t(1), t(2), t(3), p1, p2, p3)
   n  = cross(p2 - p1, p3 - p1)
   da = dot(n, a - p1)
   db = dot(n, b - p1)
   x  = a + (da / (da - db)) * (b - a)
   endfunction edge_plane_point

   pure function edge_edge_point(p0, p1, q0, q1) result(x)
   !< Canonical crossing of the coplanar, non-parallel segments (p0, p1) and (q0, q1): each segment's
   !< endpoints and then the two segments are ordered lexicographically, and the 2D solve drops the
   !< axis of the largest component of d1 x d2, all from the four points alone.
   type(vector_R8P), intent(in) :: p0, p1, q0, q1 !< Segments.
   type(vector_R8P)             :: x              !< Crossing point.
   type(vector_R8P)             :: a0, a1, b0, b1, d1, d2, w, r
   real(R8P)                    :: c1(3), c2(3), cr(3), den, tpar
   integer(I4P)                 :: k, i, j

   call sorted2(p0, p1, a0, a1)
   call sorted2(q0, q1, b0, b1)
   if (lex_less(b0, a0) .or. (same_point(a0, b0) .and. lex_less(b1, a1))) then
      d1 = a0 ; a0 = b0 ; b0 = d1
      d1 = a1 ; a1 = b1 ; b1 = d1
   endif
   d1 = a1 - a0 ; d2 = b1 - b0 ; r = b0 - a0
   w  = cross(d1, d2)
   k  = maxloc([abs(w%x), abs(w%y), abs(w%z)], dim=1)
   i  = mod(k, 3) + 1 ; j = mod(k + 1, 3) + 1
   c1 = [d1%x, d1%y, d1%z] ; c2 = [d2%x, d2%y, d2%z] ; cr = [r%x, r%y, r%z]
   den = c1(i) * c2(j) - c1(j) * c2(i)
   if (den == 0._R8P) then
      x = a0 ! parallel: caller's contract violated; return a defined point
      return
   endif
   tpar = (cr(i) * c2(j) - cr(j) * c2(i)) / den
   x = a0 + tpar * d1
   endfunction edge_edge_point

   pure function lex_less(a, b) result(less)
   !< Lexicographic order of points (x, then y, then z).
   type(vector_R8P), intent(in) :: a, b
   logical                      :: less

   if (a%x /= b%x) then
      less = a%x < b%x
   elseif (a%y /= b%y) then
      less = a%y < b%y
   else
      less = a%z < b%z
   endif
   endfunction lex_less

   subroutine coplanar_portion(u, v, t, found, p, q)
   !< Portion of segment (u, v) inside triangle t, both in the same plane (issue #33): the
   !< intersection of a segment with a convex triangle is a segment whose endpoints are among
   !<  - u and v, if inside t or on its boundary;
   !<  - vertices of t lying on (u, v);
   !<  - proper crossings of (u, v) with t's edges (canonical `edge_edge_point`).
   !< The portion is the pair of extreme candidates along (u, v); `found` is false if it is empty or
   !< a single point. The in-plane orientation tests are exact: the points are projected by dropping
   !< the coordinate of the largest component of t's normal (an exact copy of the other two), and
   !< the 2D sign is orient3d of the projected points with (0, 0, 1).
   type(vector_R8P), intent(in)  :: u, v   !< Segment (coplanar with t).
   type(vector_R8P), intent(in)  :: t(3)   !< Triangle.
   logical,          intent(out) :: found  !< True if the portion has positive length.
   type(vector_R8P), intent(out) :: p, q   !< Portion endpoints.
   type(vector_R8P)              :: cand(9), nrm
   integer(I4P)                  :: nc, k, l, m, drop, axis, ilo, ihi
   integer(I4P)                  :: ori, su(3), sv(3), sk, sl
   real(R8P)                     :: par(9), dv(3)
   real(R8P)                     :: tol    !< Vertex-merging distance.

   found = .false.
   p = u ; q = v
   tol  = SNAP_REL * max(longest_edge(t), norm(v - u))
   nrm  = cross(t(2) - t(1), t(3) - t(1))
   drop = maxloc([abs(nrm%x), abs(nrm%y), abs(nrm%z)], dim=1)
   ori  = orient_in_plane(t(1), t(2), t(3), drop, 0._R8P)
   if (ori == 0_I4P) return ! degenerate triangle
   nc = 0_I4P
   do k = 1, 3
      l = mod(k, 3) + 1
      su(k) = ori * orient_in_plane(t(k), t(l), u, drop, tol)
      sv(k) = ori * orient_in_plane(t(k), t(l), v, drop, tol)
   enddo
   if (all(su >= 0_I4P)) call add(u)
   if (all(sv >= 0_I4P)) call add(v)
   dv = [v%x - u%x, v%y - u%y, v%z - u%z]
   axis = maxloc(abs(dv), dim=1)
   do k = 1, 3
      ! vertex of t on the segment
      if (orient_in_plane(u, v, t(k), drop, tol) == 0_I4P .and. between(t(k), u, v, axis)) call add(t(k))
      ! proper crossing with edge (k, l)
      l = mod(k, 3) + 1
      if (su(k) * sv(k) < 0_I4P) then
         sk = orient_in_plane(u, v, t(k), drop, tol)
         sl = orient_in_plane(u, v, t(l), drop, tol)
         if (sk * sl < 0_I4P) call add(edge_edge_point(u, v, t(k), t(l)))
      endif
   enddo
   if (nc < 2_I4P) return
   do m = 1, nc
      par(m) = (cand(m)%x - u%x) * dv(1) + (cand(m)%y - u%y) * dv(2) + (cand(m)%z - u%z) * dv(3)
   enddo
   ilo = minloc(par(1:nc), dim=1)
   ihi = maxloc(par(1:nc), dim=1)
   ! vertex merging (see snap_to_vertices)
   p = snap_to_vertices(cand(ilo), [u, v, u], t, tol)
   q = snap_to_vertices(cand(ihi), [u, v, u], t, tol)
   if (same_point(p, q)) return
   found = .true.
   contains
      subroutine add(x)
      type(vector_R8P), intent(in) :: x
      integer(I4P)                 :: i

      do i = 1, nc
         if (same_point(cand(i), x)) return
      enddo
      nc = nc + 1_I4P
      cand(nc) = x
      endsubroutine add
   endsubroutine coplanar_portion

   ! exact cuts on implicit points (issue #65) ------------------------------------------------------------------------
   subroutine pair_cut(ta, tb, kind, p, q)
   !< `pair_classify` without the snap and with implicit points (issue #65): the signs are the exact orient3d signs,
   !< the section points are constructions (an input vertex, an edge-edge point when the edge is exactly coplanar with
   !< an edge of the other triangle, else a line-plane point), the order along the planes' line is the exact
   !< lexicographic order, and of two coincident candidates the more specific construction is kept (`ip_rank`), so
   !< every pair that produces a point produces the same identity.
   !<
   !< Precondition: both triangles exactly non-degenerate (`facet_drop` /= 0). A degenerate triangle (collinear or
   !< coincident vertices) has no plane: every orient3d against it is 0, so it is classified COPLANAR (as `ta`) or its
   !< vertices all lie "on" the other's plane (as `tb`), and the result is meaningless (issue #96: a segment 0.5 from
   !< a triangle came out PAIR_TOUCHING, and the reverse order PAIR_COPLANAR). `degenerate_pair_meet` decides such pairs
   !< by exact segment / point predicates instead. The boolean arrangement (`collect_pair`) and `find_self_intersections`
   !< pass unsanitised facets here (an open defect found by #96's audit: a boolean with a zero-area sliver as A fails
   !< with CDT_FAILED against a body it does not touch); `sanitize` removes degenerate facets.
   type(vector_R8P),       intent(in)  :: ta(3) !< First triangle.
   type(vector_R8P),       intent(in)  :: tb(3) !< Second triangle.
   integer(I4P),           intent(out) :: kind  !< PAIR_* code.
   type(implicit_point_t), intent(out) :: p, q  !< Common segment (crossing) or point (touching).
   integer(I4P)                        :: sa(3), sb(3), na, nb, i, c
   type(implicit_point_t)              :: xa(2), xb(2), alo, ahi, blo, bhi

   kind = PAIR_SEPARATED
   do i = 1, 3
      sb(i) = isgn(orient3d(ta(1), ta(2), ta(3), tb(i)))
      sa(i) = isgn(orient3d(tb(1), tb(2), tb(3), ta(i)))
   enddo
   if (all(sb == 0_I4P)) then
      kind = PAIR_COPLANAR ; return
   endif
   if (all(sb == 1_I4P) .or. all(sb == -1_I4P)) return
   if (all(sa == 1_I4P) .or. all(sa == -1_I4P)) return
   call section_ip(t=ta, s=sa, other=tb, x=xa, n=na)
   call section_ip(t=tb, s=sb, other=ta, x=xb, n=nb)
   if (na == 0_I4P .or. nb == 0_I4P) return
   call ends(xa, na, alo, ahi)
   call ends(xb, nb, blo, bhi)
   c = ip_lex_compare(alo, blo)
   if (c > 0_I4P) then
      p = alo
   elseif (c < 0_I4P) then
      p = blo
   else
      p = preferred(alo, blo)
   endif
   c = ip_lex_compare(ahi, bhi)
   if (c < 0_I4P) then
      q = ahi
   elseif (c > 0_I4P) then
      q = bhi
   else
      q = preferred(ahi, bhi)
   endif
   c = ip_lex_compare(p, q)
   if (c > 0_I4P) return
   if (c == 0_I4P) then
      kind = PAIR_TOUCHING ; p = preferred(p, q) ; q = p
   else
      kind = PAIR_CROSSING
   endif
   contains
      subroutine ends(x, n, lo, hi)
      type(implicit_point_t), intent(in)  :: x(2)
      integer(I4P),           intent(in)  :: n
      type(implicit_point_t), intent(out) :: lo, hi

      lo = x(1) ; hi = x(1)
      if (n == 2_I4P) then
         if (ip_lex_compare(x(2), x(1)) < 0_I4P) then
            lo = x(2)
         else
            hi = x(2)
         endif
      endif
      endsubroutine ends
   endsubroutine pair_cut

   subroutine section_ip(t, s, other, x, n)
   !< Points of triangle t on the plane of `other` (exact signs s): vertices with a zero sign, crossings of edges with
   !< opposite strict signs.
   type(vector_R8P),       intent(in)  :: t(3), other(3)
   integer(I4P),           intent(in)  :: s(3)
   type(implicit_point_t), intent(out) :: x(2)
   integer(I4P),           intent(out) :: n
   integer(I4P)                        :: i, j

   n = 0_I4P
   do i = 1, 3
      if (s(i) == 0_I4P .and. n < 2_I4P) then
         n = n + 1_I4P ; x(n) = ip_explicit(arr3(t(i)))
      endif
   enddo
   do i = 1, 3
      j = mod(i, 3) + 1
      if (s(i) * s(j) < 0_I4P .and. n < 2_I4P) then
         n = n + 1_I4P ; x(n) = crossing_ip(u=t(i), v=t(j), other=other)
      endif
   enddo
   endsubroutine section_ip

   function crossing_ip(u, v, other) result(x)
   !< Edge (u, v), strictly across the plane of `other`, meets it: on the line of an edge of `other` exactly coplanar
   !< with (u, v), the edge-edge point; on two such lines, their common vertex; else the line-plane point.
   type(vector_R8P), intent(in) :: u, v, other(3)
   type(implicit_point_t)       :: x
   integer(I4P)                 :: k, l, hits, last_k

   hits = 0_I4P ; last_k = 0_I4P
   do k = 1, 3
      l = mod(k, 3) + 1
      if (orient3d(other(k), other(l), u, v) == 0._R8P) then
         hits = hits + 1_I4P
         if (hits == 2_I4P) then
            if (mod(last_k, 3) + 1 == k) then
               x = ip_explicit(arr3(other(k)))
            else
               x = ip_explicit(arr3(other(last_k)))
            endif
            return
         endif
         last_k = k
      endif
   enddo
   if (hits == 1_I4P) then
      x = ip_eei(arr3(u), arr3(v), arr3(other(last_k)), arr3(other(mod(last_k, 3) + 1)))
   else
      x = ip_lpi(arr3(u), arr3(v), arr3(other(1)), arr3(other(2)), arr3(other(3)))
   endif
   endfunction crossing_ip

   subroutine coplanar_cut(u, v, t, found, p, q)
   !< `coplanar_portion` without the snap and with implicit points (issue #65): the portion of segment (u, v) inside
   !< triangle t, all exactly coplanar. Candidates: u and v inside the closed t, vertices of t strictly inside (u, v),
   !< proper crossings with t's edges (edge-edge points); the portion is the pair of extreme candidates in exact
   !< lexicographic order. In-plane signs are exact orient2d of exact coordinate copies in the projection that drops
   !< `facet_drop(t)`.
   type(vector_R8P),       intent(in)  :: u, v   !< Segment.
   type(vector_R8P),       intent(in)  :: t(3)   !< Triangle.
   logical,                intent(out) :: found  !< The portion has positive length.
   type(implicit_point_t), intent(out) :: p, q   !< Portion endpoints.
   type(implicit_point_t)              :: cand(9)
   integer(I4P)                        :: nc, k, l, m, drop, ori, ilo, ihi, su(3), sv(3), ax
   real(R8P)                           :: au(3), av(3), at(3, 3)

   found = .false.
   au = arr3(u) ; av = arr3(v)
   do k = 1, 3 ; at(:, k) = arr3(t(k)) ; enddo
   drop = facet_drop(t)
   if (drop == 0_I4P) return
   ori = o2p(at(:, 1), at(:, 2), at(:, 3), drop)
   nc = 0_I4P
   do k = 1, 3
      l = mod(k, 3) + 1
      su(k) = ori * o2p(at(:, k), at(:, l), au, drop)
      sv(k) = ori * o2p(at(:, k), at(:, l), av, drop)
   enddo
   if (all(su >= 0_I4P)) call add(ip_explicit(au))
   if (all(sv >= 0_I4P)) call add(ip_explicit(av))
   ax = maxloc(abs(av - au), dim=1)
   do k = 1, 3
      l = mod(k, 3) + 1
      if (o2p(au, av, at(:, k), drop) == 0_I4P) then
         if ((at(ax, k) - au(ax)) * (av(ax) - at(ax, k)) > 0._R8P) call add(ip_explicit(at(:, k)))
      endif
      if (su(k) * sv(k) < 0_I4P) then
         if (o2p(au, av, at(:, k), drop) * o2p(au, av, at(:, l), drop) < 0_I4P) &
            call add(ip_eei(au, av, at(:, k), at(:, l)))
      endif
   enddo
   if (nc < 2_I4P) return
   ilo = 1_I4P ; ihi = 1_I4P
   do m = 2, nc
      if (ip_lex_compare(cand(m), cand(ilo)) < 0_I4P) ilo = m
      if (ip_lex_compare(cand(m), cand(ihi)) > 0_I4P) ihi = m
   enddo
   if (ip_lex_compare(cand(ilo), cand(ihi)) == 0_I4P) return
   p = cand(ilo) ; q = cand(ihi) ; found = .true.
   contains
      subroutine add(x)
      type(implicit_point_t), intent(in) :: x
      integer(I4P)                       :: i

      do i = 1, nc
         if (ip_lex_compare(cand(i), x) == 0_I4P) return
      enddo
      nc = nc + 1_I4P ; cand(nc) = x
      endsubroutine add
   endsubroutine coplanar_cut

   pure function facet_drop(t) result(drop)
   !< An axis whose dropping keeps triangle t non-degenerate, exactly: the axes by decreasing |normal component|, the
   !< first with a nonzero exact orient2d of the projection (0 if t is degenerate).
   type(vector_R8P), intent(in) :: t(3) !< Triangle.
   integer(I4P)                 :: drop !< Axis (1..3), 0 if none.
   type(vector_R8P)             :: n
   real(R8P)                    :: a(3, 3), w(3)
   integer(I4P)                 :: k, order(3)

   do k = 1, 3 ; a(:, k) = arr3(t(k)) ; enddo
   n = cross(t(2) - t(1), t(3) - t(1))
   w = abs([n%x, n%y, n%z])
   order(1) = maxloc(w, dim=1) ; w(order(1)) = -1._R8P
   order(2) = maxloc(w, dim=1) ; w(order(2)) = -1._R8P
   order(3) = maxloc(w, dim=1)
   do k = 1, 3
      drop = order(k)
      if (o2p(a(:, 1), a(:, 2), a(:, 3), drop) /= 0_I4P) return
   enddo
   drop = 0_I4P
   endfunction facet_drop

   pure function o2p(a, b, c, drop) result(s)
   !< Exact sign of the orientation of a, b, c in the kept axes (drop + 1, drop + 2) mod 3.
   real(R8P),    intent(in) :: a(3), b(3), c(3)
   integer(I4P), intent(in) :: drop
   integer(I4P)             :: s
   integer(I4P)             :: i, j

   i = mod(drop, 3_I4P) + 1_I4P ; j = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   s = isgn(orient2d_xy(a(i), a(j), b(i), b(j), c(i), c(j)))
   endfunction o2p

   pure function degenerate_pair_meet(a, b, drop_a, drop_b) result(meet)
   !< Exact contact of triangles a and b (vertices as columns) when at least one is exactly degenerate (`facet_drop` = 0:
   !< collinear or coincident vertices; issue #96, shared since issue #99). It complements `pair_cut`, whose
   !< precondition is two non-degenerate triangles (a degenerate one has no plane: every orient3d against it is 0).
   !< A degenerate triangle is the segment between its lexicographic extreme vertices (exactly collinear vertices: the
   !< lexicographic order is their order along the line), a point if they coincide; the contact is then decided by
   !< exact predicates only, no tolerance: segment / point against a regular triangle by orient3d against its plane and
   !< edges (orient2d in the projection dropping its `facet_drop` when in the plane), segment / point against segment /
   !< point by equality, orient3d coplanarity, lexicographic intervals when collinear, orient2d in a projection injective
   !< on their plane otherwise. Closed sets: touching at a single point counts as contact. If neither triangle is
   !< degenerate the result is meaningless (use `pair_cut`).
   real(R8P),    intent(in)           :: a(3, 3), b(3, 3) !< Triangles (at least one exactly degenerate).
   integer(I4P), intent(in), optional :: drop_a           !< `facet_drop` of a, if already known (computed if absent).
   integer(I4P), intent(in), optional :: drop_b           !< `facet_drop` of b, if already known (computed if absent).
   logical                            :: meet             !< They share a point.
   real(R8P)                          :: a0(3), a1(3), b0(3), b1(3)
   integer(I4P)                       :: da, db, k
   type(vector_R8P)                   :: v(3)

   if (present(drop_a)) then
      da = drop_a
   else
      do k = 1, 3 ; v(k) = vector_R8P(a(1, k), a(2, k), a(3, k)) ; enddo
      da = facet_drop(v)
   endif
   if (present(drop_b)) then
      db = drop_b
   else
      do k = 1, 3 ; v(k) = vector_R8P(b(1, k), b(2, k), b(3, k)) ; enddo
      db = facet_drop(v)
   endif
   if (da == 0_I4P .and. db == 0_I4P) then
      call hull_ends(a, a0, a1) ; call hull_ends(b, b0, b1)
      meet = segment_segment_meet(a0, a1, b0, b1)
   elseif (da == 0_I4P) then
      call hull_ends(a, a0, a1)
      meet = segment_triangle_meet(a0, a1, b, db)
   elseif (db == 0_I4P) then
      call hull_ends(b, b0, b1)
      meet = segment_triangle_meet(b0, b1, a, da)
   else
      meet = .false.
   endif
   endfunction degenerate_pair_meet

   pure subroutine hull_ends(t, e0, e1)
   !< The lexicographically smallest and largest vertices of t.
   real(R8P), intent(in)  :: t(3, 3) !< Triangle.
   real(R8P), intent(out) :: e0(3)   !< Smallest vertex.
   real(R8P), intent(out) :: e1(3)   !< Largest vertex.
   integer(I4P)           :: k

   e0 = t(:, 1) ; e1 = t(:, 1)
   do k = 2, 3
      if (lex_lt(t(:, k), e0)) e0 = t(:, k)
      if (lex_lt(e1, t(:, k))) e1 = t(:, k)
   enddo
   endsubroutine hull_ends

   pure function lex_lt(p, q) result(less)
   !< p < q in lexicographic order (x, then y, then z).
   real(R8P), intent(in) :: p(3), q(3) !< Points.
   logical                :: less      !< p < q.

   if (p(1) /= q(1)) then
      less = p(1) < q(1)
   elseif (p(2) /= q(2)) then
      less = p(2) < q(2)
   else
      less = p(3) < q(3)
   endif
   endfunction lex_lt

   pure function segment_triangle_meet(p, q, t, drop) result(meet)
   !< The closed segment (p, q) (a point if p = q) and the regular triangle t (`facet_drop(t)` = drop /= 0) share a
   !< point (exact): an end on the plane inside t; the segment strictly across the plane and its line through the closed
   !< t (orient3d against the three edges, all of one sign or zero); the segment in the plane meeting t (orient2d in the
   !< projection dropping `drop`, injective on the plane).
   real(R8P),    intent(in) :: p(3), q(3), t(3, 3) !< Segment, triangle.
   integer(I4P), intent(in) :: drop                !< `facet_drop(t)`.
   logical                  :: meet                !< They share a point.
   integer(I4P)             :: sp, sq, s(3), k

   sp = isgn(orient3d_xyz(t(:, 1), t(:, 2), t(:, 3), p)) ; sq = isgn(orient3d_xyz(t(:, 1), t(:, 2), t(:, 3), q))
   if (sp * sq > 0_I4P) then
      meet = .false.
   elseif (sp == 0_I4P .and. sq == 0_I4P) then
      meet = in_triangle_2d(p, t, drop) .or. in_triangle_2d(q, t, drop)
      do k = 1, 3
         if (meet) exit
         meet = segments_meet_2d(p, q, t(:, k), t(:, mod(k, 3) + 1), drop)
      enddo
   elseif (sp == 0_I4P) then
      meet = in_triangle_2d(p, t, drop)
   elseif (sq == 0_I4P) then
      meet = in_triangle_2d(q, t, drop)
   else
      do k = 1, 3
         s(k) = isgn(orient3d_xyz(p, q, t(:, k), t(:, mod(k, 3) + 1)))
      enddo
      meet = all(s >= 0_I4P) .or. all(s <= 0_I4P)
   endif
   endfunction segment_triangle_meet

   pure function segment_segment_meet(a0, a1, b0, b1) result(meet)
   !< The closed segments (a0, a1) and (b0, b1), each with its ends in lexicographic order (a point if they coincide),
   !< share a point (exact): points by equality or collinearity; segments not coplanar (orient3d) never; collinear by
   !< their lexicographic intervals; otherwise in a projection injective on their plane (orient2d).
   real(R8P), intent(in) :: a0(3), a1(3), b0(3), b1(3) !< Segments.
   logical               :: meet                       !< They share a point.
   logical               :: c0, c1

   if (all(a0 == a1) .and. all(b0 == b1)) then
      meet = all(a0 == b0)
   elseif (all(a0 == a1)) then
      meet = on_segment_3d(a0, b0, b1)
   elseif (all(b0 == b1)) then
      meet = on_segment_3d(b0, a0, a1)
   elseif (orient3d_xyz(a0, a1, b0, b1) /= 0._R8P) then
      meet = .false.
   else
      c0 = collinear_3d(a0, a1, b0) ; c1 = collinear_3d(a0, a1, b1)
      if (c0 .and. c1) then
         meet = .not. (lex_lt(a1, b0) .or. lex_lt(b1, a0))
      elseif (.not. c0) then
         meet = segments_meet_2d(a0, a1, b0, b1, nonzero_drop(a0, a1, b0))
      else
         meet = segments_meet_2d(a0, a1, b0, b1, nonzero_drop(a0, a1, b1))
      endif
   endif
   endfunction segment_segment_meet

   pure function collinear_3d(p, q, r) result(yes)
   !< p, q, r exactly collinear: the three coordinate projections of their cross product vanish.
   real(R8P), intent(in) :: p(3), q(3), r(3) !< Points.
   logical               :: yes              !< Collinear.

   yes = o2p(p, q, r, 1_I4P) == 0_I4P .and. o2p(p, q, r, 2_I4P) == 0_I4P .and. o2p(p, q, r, 3_I4P) == 0_I4P
   endfunction collinear_3d

   pure function nonzero_drop(p, q, r) result(drop)
   !< The first axis whose dropping keeps p, q, r (not collinear) a triangle (exact orient2d).
   real(R8P), intent(in) :: p(3), q(3), r(3) !< Points.
   integer(I4P)          :: drop             !< Axis.

   do drop = 1_I4P, 2_I4P
      if (o2p(p, q, r, drop) /= 0_I4P) return
   enddo
   drop = 3_I4P
   endfunction nonzero_drop

   pure function on_segment_3d(x, p, q) result(yes)
   !< x on the closed segment (p, q), p /= q (exact): collinear, inside the box of the segment.
   real(R8P), intent(in) :: x(3), p(3), q(3) !< Point, segment.
   logical               :: yes              !< On it.

   yes = collinear_3d(p, q, x)
   if (yes) yes = all(x >= min(p, q)) .and. all(x <= max(p, q))
   endfunction on_segment_3d

   pure function in_triangle_2d(x, t, drop) result(yes)
   !< x (on the plane of t) in the closed t, in the projection dropping `drop` (t not degenerate there).
   real(R8P),    intent(in) :: x(3), t(3, 3) !< Point, triangle.
   integer(I4P), intent(in) :: drop          !< Projection.
   logical                  :: yes           !< Inside or on the boundary.
   integer(I4P)             :: s1, s2, s3

   s1 = o2p(t(:, 1), t(:, 2), x, drop) ; s2 = o2p(t(:, 2), t(:, 3), x, drop) ; s3 = o2p(t(:, 3), t(:, 1), x, drop)
   yes = (s1 >= 0_I4P .and. s2 >= 0_I4P .and. s3 >= 0_I4P) .or. (s1 <= 0_I4P .and. s2 <= 0_I4P .and. s3 <= 0_I4P)
   endfunction in_triangle_2d

   pure function segments_meet_2d(p1, q1, p2, q2, drop) result(yes)
   !< The closed segments (p1, q1) and (p2, q2), coplanar, meet in the projection dropping `drop` (injective on their
   !< plane): a proper crossing, or an end on the other segment (collinear and inside its box).
   real(R8P),    intent(in) :: p1(3), q1(3), p2(3), q2(3) !< Segments.
   integer(I4P), intent(in) :: drop                       !< Projection.
   logical                  :: yes                        !< They meet.
   integer(I4P)             :: d1, d2, d3, d4

   d1 = o2p(p2, q2, p1, drop) ; d2 = o2p(p2, q2, q1, drop) ; d3 = o2p(p1, q1, p2, drop) ; d4 = o2p(p1, q1, q2, drop)
   yes = d1 * d2 < 0_I4P .and. d3 * d4 < 0_I4P
   if (yes) return
   yes = (d1 == 0_I4P .and. in_box_2d(p1, p2, q2, drop)) .or. (d2 == 0_I4P .and. in_box_2d(q1, p2, q2, drop)) .or. &
         (d3 == 0_I4P .and. in_box_2d(p2, p1, q1, drop)) .or. (d4 == 0_I4P .and. in_box_2d(q2, p1, q1, drop))
   endfunction segments_meet_2d

   pure function in_box_2d(x, p, q, drop) result(yes)
   !< x inside the box of (p, q) in the kept coordinates of the projection dropping `drop`.
   real(R8P),    intent(in) :: x(3), p(3), q(3) !< Point, segment.
   integer(I4P), intent(in) :: drop             !< Projection.
   logical                  :: yes              !< Inside.
   integer(I4P)             :: u, v

   u = mod(drop, 3_I4P) + 1_I4P ; v = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   yes = x(u) >= min(p(u), q(u)) .and. x(u) <= max(p(u), q(u)) .and. x(v) >= min(p(v), q(v)) .and. &
         x(v) <= max(p(v), q(v))
   endfunction in_box_2d

   pure function preferred(a, b) result(c)
   !< Of two coincident constructions, the more specific (`ip_rank`), else a.
   type(implicit_point_t), intent(in) :: a, b
   type(implicit_point_t)             :: c

   c = a
   if (ip_rank(b) < ip_rank(a)) c = b
   endfunction preferred

   pure function arr3(a) result(r)
   type(vector_R8P), intent(in) :: a
   real(R8P)                    :: r(3)

   r = [a%x, a%y, a%z]
   endfunction arr3

   pure function isgn(v) result(s)
   real(R8P), intent(in) :: v
   integer(I4P)          :: s

   s = 0_I4P
   if (v > 0._R8P) s = 1_I4P
   if (v < 0._R8P) s = -1_I4P
   endfunction isgn

   ! private helpers
   function orient_in_plane(a, b, c, drop, tol) result(s)
   !< Sign of the 2D orientation of (a, b, c) in the coordinate plane obtained by dropping axis
   !< `drop`: orient3d of the projected points (exact copies of two coordinates, z = 0) and
   !< (0, 0, 1), exact; 0 when c is within `tol` of the line (a, b) (the snap of the module: a
   !< vertex meant to lie on the other mesh's edge is off it by round-off once rotated, and an
   !< exact sign would then call nearly collinear segments crossing and build a far-off point).
   type(vector_R8P), intent(in) :: a, b, c !< Points.
   integer(I4P),     intent(in) :: drop    !< Dropped axis (1, 2, 3).
   real(R8P),        intent(in) :: tol     !< Snap distance (0: exact sign).
   integer(I4P)                 :: s       !< -1, 0, 1.
   type(vector_R8P)             :: pa, pb, pc
   real(R8P)                    :: det, len_ab

   pa = proj(a, drop) ; pb = proj(b, drop) ; pc = proj(c, drop)
   det = orient3d(pa, pb, pc, vector_R8P(0._R8P, 0._R8P, 1._R8P))
   len_ab = sqrt((pb%x - pa%x)**2 + (pb%y - pa%y)**2)
   if (abs(det) <= tol * len_ab) then
      s = 0_I4P
   else
      s = sgn(det)
   endif
   endfunction orient_in_plane

   pure function proj(a, drop) result(r)
   type(vector_R8P), intent(in) :: a
   integer(I4P),     intent(in) :: drop
   type(vector_R8P)             :: r

   select case (drop)
   case (1) ; r%x = a%y ; r%y = a%z
   case (2) ; r%x = a%z ; r%y = a%x
   case default ; r%x = a%x ; r%y = a%y
   endselect
   r%z = 0._R8P
   endfunction proj

   pure function between(w, a, b, axis) result(yes)
   !< w (collinear with a, b) lies on the closed segment (a, b): compared on the axis of the largest
   !< component of b - a.
   type(vector_R8P), intent(in) :: w, a, b
   integer(I4P),     intent(in) :: axis
   logical                      :: yes
   real(R8P)                    :: cw, ca, cb

   select case (axis)
   case (1) ; cw = w%x ; ca = a%x ; cb = b%x
   case (2) ; cw = w%y ; ca = a%y ; cb = b%y
   case default ; cw = w%z ; ca = a%z ; cb = b%z
   endselect
   yes = cw >= min(ca, cb) .and. cw <= max(ca, cb)
   endfunction between

   pure function same_point(a, b) result(same)
   type(vector_R8P), intent(in) :: a, b
   logical                      :: same

   same = a%x == b%x .and. a%y == b%y .and. a%z == b%z
   endfunction same_point

   pure subroutine sorted2(u, v, a, b)
   type(vector_R8P), intent(in)  :: u, v
   type(vector_R8P), intent(out) :: a, b

   if (lex_less(v, u)) then
      a = v ; b = u
   else
      a = u ; b = v
   endif
   endsubroutine sorted2

   pure subroutine sorted3(u, v, w, a, b, c)
   type(vector_R8P), intent(in)  :: u, v, w
   type(vector_R8P), intent(out) :: a, b, c
   type(vector_R8P)              :: t

   a = u ; b = v ; c = w
   if (lex_less(b, a)) then ; t = a ; a = b ; b = t ; endif
   if (lex_less(c, b)) then ; t = b ; b = c ; c = t ; endif
   if (lex_less(b, a)) then ; t = a ; a = b ; b = t ; endif
   endsubroutine sorted3

   pure function snap_to_vertices(x, ta, tb, tol) result(y)
   !< Vertex merging (issue #33): `x` replaced by the nearest vertex of `ta` or `tb` if that vertex
   !< is within `tol`, else unchanged. A point that coincides with an input vertex up to round-off
   !< (a vertex lying on the other mesh's edge, after a rotation) is otherwise built by one pair as
   !< a copy of the vertex and by another as an edge-edge crossing, with different last bits, and
   !< the exact-equality weld then leaves a crack.
   type(vector_R8P), intent(in) :: x            !< Point.
   type(vector_R8P), intent(in) :: ta(3), tb(3) !< Candidate vertices.
   real(R8P),        intent(in) :: tol          !< Merging distance.
   type(vector_R8P)             :: y            !< Merged point.
   real(R8P)                    :: best, d
   integer(I4P)                 :: i

   y = x
   best = tol
   do i = 1, 3
      d = norm(ta(i) - x)
      if (d <= best) then ; best = d ; y = ta(i) ; endif
      d = norm(tb(i) - x)
      if (d <= best) then ; best = d ; y = tb(i) ; endif
   enddo
   endfunction snap_to_vertices

   pure function snapped_sign(det, area2, tol) result(s)
   !< Sign of an orient3d determinant, 0 when the point it measures is within `tol` of the plane:
   !< det / area2 is that point's signed distance (area2 = twice the reference triangle's area).
   real(R8P), intent(in) :: det, area2, tol
   integer(I4P)          :: s

   if (abs(det) <= tol * area2) then
      s = 0_I4P
   else
      s = sgn(det)
   endif
   endfunction snapped_sign

   pure function longest_edge(t) result(l)
   type(vector_R8P), intent(in) :: t(3)
   real(R8P)                    :: l

   l = max(norm(t(2) - t(1)), norm(t(3) - t(2)), norm(t(1) - t(3)))
   endfunction longest_edge

   pure function norm(a) result(n)
   type(vector_R8P), intent(in) :: a
   real(R8P)                    :: n

   n = sqrt(a%x * a%x + a%y * a%y + a%z * a%z)
   endfunction norm

   pure function sgn(v) result(s)
   real(R8P), intent(in) :: v
   integer(I4P)          :: s

   s = 0_I4P
   if (v > 0._R8P) s = 1_I4P
   if (v < 0._R8P) s = -1_I4P
   endfunction sgn

   pure function cross(a, b) result(c)
   type(vector_R8P), intent(in) :: a, b
   type(vector_R8P)             :: c

   c%x = a%y * b%z - a%z * b%y
   c%y = a%z * b%x - a%x * b%z
   c%z = a%x * b%y - a%y * b%x
   endfunction cross

   pure function dot(a, b) result(d)
   type(vector_R8P), intent(in) :: a, b
   real(R8P)                    :: d

   d = a%x * b%x + a%y * b%y + a%z * b%z
   endfunction dot
endmodule fossil_tritri_exact
