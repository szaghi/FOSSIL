!run figck-wall-spacing fig_ck_wall_spacing
!fig ck-wall-spacing
!< Figure of the recipe "Extrude an O-grid and check its wall spacing": the NACA0012 section and the O-grid marched
!< from it (every second grid line), seen from above, zoomed on the nose.
program fig_ck_wall_spacing
use fossil,     only : curvilinear_block_2d_t, hyperbolic_grid_2d, geometric_spacing
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: g
type(vector_R8P)             :: p(20000), q(20000), w0(128), w1(128)
real(R8P)                    :: u(129), v(129), xc, yt
integer(I4P)                 :: i, j, n, un

do i = 1, 129
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / 129._R8P))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= 65)
enddo
call hyperbolic_grid_2d(u, v, geometric_spacing(1e-3_R8P, 1.15_R8P, 30_I4P), g)
n = 0
do j = 1, g%nj, 2
   do i = 1, g%ni - 1
      n = n + 1 ; p(n) = vector_R8P(g%u(i, j), g%v(i, j), 0._R8P) ; q(n) = vector_R8P(g%u(i + 1, j), g%v(i + 1, j), 0._R8P)
   enddo
enddo
do i = 1, g%ni, 2
   do j = 1, g%nj - 1
      n = n + 1 ; p(n) = vector_R8P(g%u(i, j), g%v(i, j), 0._R8P) ; q(n) = vector_R8P(g%u(i, j + 1), g%v(i, j + 1), 0._R8P)
   enddo
enddo
do i = 1, 128
   w0(i) = vector_R8P(u(i), v(i), 0._R8P) ; w1(i) = vector_R8P(u(mod(i, 129) + 1), v(mod(i, 129) + 1), 0._R8P)
enddo
call fig_open(un, 'ck-wall-spacing', 'hyperbolic_grid_2d: NACA0012 O-grid', -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.25_R8P, 0._R8P, 0._R8P), 0.45_R8P)
call fig_segments(un, 'ghost', p(1:n), q(1:n))
call fig_segments(un, 'accent', w0, w1)
call fig_close(un)
endprogram fig_ck_wall_spacing
