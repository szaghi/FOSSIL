!< Figure and run of the reference page "Curvature", section "Principal curvatures and directions":
!< - a parametric torus (R = 1, r = 0.45) at 48 x 24, 96 x 48 and 192 x 96 vertices; at the feet of 400 points 0.05
!<   off the surface, `query_many(principal_curvatures=, principal_directions=)` against the exact values: k1 = 1/r
!<   around the tube (d1 along the meridian), k2 = cos(b) / (R + r cos(b)) along it (d2 along the parallel), b the angle
!<   around the tube; the curvature errors and the angle of d1 to the meridian in the foot's tangent plane fall by about
!<   4 per halving of h, the tilt of that plane (the foot normal is the facet's) by 2;
!< - the figure: the 48 x 24 torus with d1 (ember) and d2 (teal) at the feet of the 400 points.
program fig_ref_curvature_principal
use fossil,     only : surface_stl_object, facet_object
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P), RR = 1._R8P, RT = 0.45_R8P
integer(I4P), parameter         :: NP = 400
type(surface_stl_object)        :: torus
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: q(NP), a1(NP), b1(NP), a2(NP), b2(NP)
real(R8P)                       :: d(NP), cp(3, NP), nm(3, NP), k(2, NP), dirs(3, 2, NP), ab(2, NP), r(2), &
                                   e1, e2, ea, et, b, rho, kx, tb(3), tp_(3), c(3)
integer(I4P)                    :: lev, nu, nv, i, m, u

do i = 1, NP ! quasi-random angles (golden ratio sequences): reproducible on any compiler
   ab(:, i) = 2 * PI * [modulo(i * 0.6180339887498949_R8P, 1._R8P), modulo(i * 0.7548776662466927_R8P, 1._R8P)]
   q(i) = vector_R8P((RR + (RT + 0.05_R8P) * cos(ab(2, i))) * cos(ab(1, i)), &
                     (RR + (RT + 0.05_R8P) * cos(ab(2, i))) * sin(ab(1, i)), (RT + 0.05_R8P) * sin(ab(2, i)))
enddo
print '(A)', '  vertices   facets   max|k1 - 1/r|   max|k2 - exact|   d1 vs meridian [deg]: in plane   normal tilt'
do lev = 0, 2
   nu = 48 * 2**lev ; nv = 24 * 2**lev
   call make_torus(nu, nv, f) ; call torus%adopt_facets(facets=f)
   call torus%build_curvature                    ! the per-vertex shape operator (cubic jet over the 2-ring)
   call torus%query_many(points=q, distances=d, closest_points=cp, normals=nm, principal_curvatures=k, &
                         principal_directions=dirs)
   e1 = 0._R8P ; e2 = 0._R8P ; ea = 0._R8P ; et = 0._R8P
   do i = 1, NP
      rho = hypot(cp(1, i), cp(2, i)) ; b = atan2(cp(3, i), rho - RR)
      kx = cos(b) / (RR + RT * cos(b))
      tb = [-sin(b) * cp(1, i) / rho, -sin(b) * cp(2, i) / rho, cos(b)]     ! the meridian direction
      e1 = max(e1, abs(k(1, i) - 1._R8P / RT)) ; e2 = max(e2, abs(k(2, i) - kx))
      tp_ = tb - dot_product(tb, nm(:, i)) * nm(:, i) ; tp_ = tp_ / norm2(tp_) ! the meridian in the foot's tangent plane
      ea = max(ea, acos(min(1._R8P, abs(dot_product(dirs(:, 1, i), tp_)))) * 180._R8P / PI)
      et = max(et, asin(min(1._R8P, abs(dot_product(tb, nm(:, i))))) * 180._R8P / PI)      ! the foot normal's tilt
   enddo
   print '(I10,I9,ES16.3,ES18.3,F21.4,F14.4)', nu * nv, 2 * nu * nv, e1, e2, ea, et
   if (lev == 0) then
      m = 0
      do i = 1, NP
         m = m + 1
         c = cp(:, i) + 0.01_R8P * nm(:, i)
         a1(m) = vec(c - 0.07_R8P * dirs(:, 1, i)) ; b1(m) = vec(c + 0.07_R8P * dirs(:, 1, i))
         a2(m) = vec(c - 0.07_R8P * dirs(:, 2, i)) ; b2(m) = vec(c + 0.07_R8P * dirs(:, 2, i))
      enddo
      call fig_open(u, 'ref-curvature-principal', 'torus 48 x 24: d1 (k1 = 1/r, around the tube) ember, d2 (along it) teal', &
                    -60._R8P, 38._R8P)
      call fig_surface(u, 'body', torus)
      call fig_segments(u, 'patch', a1(1:m), b1(1:m))
      call fig_segments(u, 'probe', a2(1:m), b2(1:m))
      call fig_close(u)
   endif
enddo
print '(A,F7.4,A,F7.4,A,F7.4)', 'k2 exact from ', -1._R8P / (RR - RT), ' (inner equator) to ', 1._R8P / (RR + RT), &
      ' (outer); 1/r = ', 1._R8P / RT

contains
   pure function vec(x) result(v)
   !< Array to vector.
   real(R8P), intent(in) :: x(3)
   type(vector_R8P)      :: v

   v = vector_R8P(x(1), x(2), x(3))
   endfunction vec

   subroutine make_torus(nu, nv, fo)
   !< The parametric torus, nu x nv vertices, two facets per quad (outward normals); seams closed bit for bit.
   integer(I4P),                    intent(in)  :: nu, nv
   type(facet_object), allocatable, intent(out) :: fo(:)
   integer(I4P)                                 :: i, j, n

   allocate(fo(2 * nu * nv))
   n = 0
   do j = 0, nv - 1 ; do i = 0, nu - 1
      n = n + 1
      fo(n)%vertex(1) = tp(i, j, nu, nv) ; fo(n)%vertex(2) = tp(i + 1, j, nu, nv) ; fo(n)%vertex(3) = tp(i + 1, j + 1, nu, nv)
      call fo(n)%compute_metrix
      n = n + 1
      fo(n)%vertex(1) = tp(i, j, nu, nv) ; fo(n)%vertex(2) = tp(i + 1, j + 1, nu, nv) ; fo(n)%vertex(3) = tp(i, j + 1, nu, nv)
      call fo(n)%compute_metrix
   enddo ; enddo
   endsubroutine make_torus

   function tp(i, j, nu, nv) result(v)
   !< Vertex (i, j) of the nu x nv torus, indices modulo the grid.
   integer(I4P), intent(in) :: i, j, nu, nv
   type(vector_R8P)         :: v
   real(R8P)                :: a, b

   a = 2._R8P * PI * mod(i, nu) / nu ; b = 2._R8P * PI * mod(j, nv) / nv
   v = vector_R8P((RR + RT * cos(b)) * cos(a), (RR + RT * cos(b)) * sin(a), RT * sin(b))
   endfunction tp
endprogram fig_ref_curvature_principal
