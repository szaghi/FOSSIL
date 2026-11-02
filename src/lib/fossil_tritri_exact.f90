!< FOSSIL, triangle-triangle intersection with exact predicates and canonical points (issue #33).
module fossil_tritri_exact
!< FOSSIL, triangle-triangle intersection with exact predicates and canonical points (issue #33).
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
use fossil_predicates, only : orient3d
use penf,              only : I4P, R8P
use vecfor,            only : vector_R8P

implicit none
private
public :: PAIR_SEPARATED, PAIR_CROSSING, PAIR_TOUCHING, PAIR_COPLANAR
public :: pair_classify
public :: coplanar_portion
public :: edge_plane_point, edge_edge_point
public :: lex_less

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

   found = .false.
   p = u ; q = v
   nrm  = cross(t(2) - t(1), t(3) - t(1))
   drop = maxloc([abs(nrm%x), abs(nrm%y), abs(nrm%z)], dim=1)
   ori  = orient_in_plane(t(1), t(2), t(3), drop)
   if (ori == 0_I4P) return ! degenerate triangle
   nc = 0_I4P
   do k = 1, 3
      l = mod(k, 3) + 1
      su(k) = ori * orient_in_plane(t(k), t(l), u, drop)
      sv(k) = ori * orient_in_plane(t(k), t(l), v, drop)
   enddo
   if (all(su >= 0_I4P)) call add(u)
   if (all(sv >= 0_I4P)) call add(v)
   dv = [v%x - u%x, v%y - u%y, v%z - u%z]
   axis = maxloc(abs(dv), dim=1)
   do k = 1, 3
      ! vertex of t on the segment
      if (orient_in_plane(u, v, t(k), drop) == 0_I4P .and. between(t(k), u, v, axis)) call add(t(k))
      ! proper crossing with edge (k, l)
      l = mod(k, 3) + 1
      if (su(k) * sv(k) < 0_I4P) then
         sk = orient_in_plane(u, v, t(k), drop)
         sl = orient_in_plane(u, v, t(l), drop)
         if (sk * sl < 0_I4P) call add(edge_edge_point(u, v, t(k), t(l)))
      endif
   enddo
   if (nc < 2_I4P) return
   do m = 1, nc
      par(m) = (cand(m)%x - u%x) * dv(1) + (cand(m)%y - u%y) * dv(2) + (cand(m)%z - u%z) * dv(3)
   enddo
   ilo = minloc(par(1:nc), dim=1)
   ihi = maxloc(par(1:nc), dim=1)
   if (same_point(cand(ilo), cand(ihi))) return
   p = cand(ilo) ; q = cand(ihi)
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

   ! private helpers
   function orient_in_plane(a, b, c, drop) result(s)
   !< Exact sign of the 2D orientation of (a, b, c) in the coordinate plane obtained by dropping
   !< axis `drop`: orient3d of the projected points (exact copies of two coordinates, z = 0) and
   !< (0, 0, 1).
   type(vector_R8P), intent(in) :: a, b, c !< Points.
   integer(I4P),     intent(in) :: drop    !< Dropped axis (1, 2, 3).
   integer(I4P)                 :: s       !< -1, 0, 1.

   s = sgn(orient3d(proj(a, drop), proj(b, drop), proj(c, drop), vector_R8P(0._R8P, 0._R8P, 1._R8P)))
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
