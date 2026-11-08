!run figck-smooth fig_ck_smooth
!fig ck-smooth-before
!fig ck-smooth-after
!< Figures of the cookbook recipe "Smooth a noisy surface": the marching-cubes unit sphere with its vertices moved
!< radially by up to 3 %, then after 10 Taubin (lambda, mu) pairs, from the same camera.
program fig_ck_smooth
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, SMOOTH_METHOD_TAUBIN
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: f(29, 29, 29), x(29)
integer(I4P)                    :: i, j, k, u

x = [(-1.4_R8P + 0.1_R8P * i, i = 0, 28)]
do k = 1, 29 ; do j = 1, 29 ; do i = 1, 29
   f(i, j, k) = norm2([x(i), x(j), x(k)]) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), x(1), x(1)), bmax=vector_R8P(x(29), x(29), x(29)), &
                        iso=0._R8P, surface=facets)
do i = 1, size(facets) ; do j = 1, 3
   facets(i)%vertex(j) = (1._R8P + 0.03_R8P * sin(1.e3_R8P * (facets(i)%vertex(j)%x + 2 * facets(i)%vertex(j)%y + &
                         3 * facets(i)%vertex(j)%z))) * facets(i)%vertex(j)
enddo ; enddo
call sphere%adopt_facets(facets=facets)
call fig_open(u, 'ck-smooth-before', 'noisy: vertices moved radially by up to 3 %', -50._R8P, 30._R8P)
call fig_surface(u, 'body', sphere)
call fig_close(u)
call sphere%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=10_I4P)
call fig_open(u, 'ck-smooth-after', 'smooth(SMOOTH_METHOD_TAUBIN, iterations=10)', -50._R8P, 30._R8P)
call fig_surface(u, 'body', sphere)
call fig_close(u)
endprogram fig_ck_smooth
