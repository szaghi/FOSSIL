!run figck-principal-curvature fig_ck_principal_curvature
!fig ck-principal-curvature
!< Figure of the cookbook recipe "Principal curvatures and directions at the wall": the ellipsoid of
!< ck_principal_curvature (semi-axes 2, 1.5, 1) with, at the feet of points spread around it, the principal directions:
!< d1 (the larger curvature) in ember, d2 in teal, each a short segment centred on the foot.
program fig_ck_principal_curvature
use fossil,     only : surface_stl_object, facet_object, extract_isosurface
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter         :: NA = 24, NB = 11
type(surface_stl_object)        :: body
type(facet_object), allocatable :: facets(:)
type(vector_R8P)                :: q(NA * NB), a1(NA * NB), b1(NA * NB), a2(NA * NB), b2(NA * NB)
real(R8P)                       :: phi(61, 51, 41), d(NA * NB), cp(3, NA * NB), nm(3, NA * NB), k(2, NA * NB), &
                                   dirs(3, 2, NA * NB), al, be, c(3), e1(3), e2(3)
integer(I4P)                    :: i, j, l, m, u

do l = 1, 41 ; do j = 1, 51 ; do i = 1, 61
   phi(i, j, l) = norm2([(-2.4_R8P + 0.08_R8P * (i - 1)) / 2, (-2._R8P + 0.08_R8P * (j - 1)) / 1.5_R8P, &
                         -1.6_R8P + 0.08_R8P * (l - 1)]) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-2.4_R8P, -2._R8P, -1.6_R8P), &
                        bmax=vector_R8P(2.4_R8P, 2._R8P, 1.6_R8P), iso=0._R8P, surface=facets)
call body%adopt_facets(facets=facets)
call body%build_curvature
m = 0
do j = 1, NB ; do i = 1, NA
   m = m + 1
   al = 2 * acos(-1._R8P) * (i - 0.5_R8P) / NA ; be = acos(-1._R8P) * (j - 0.5_R8P) / NB - acos(-1._R8P) / 2
   q(m) = vector_R8P(2.3_R8P * cos(be) * cos(al), 1.8_R8P * cos(be) * sin(al), 1.3_R8P * sin(be))
enddo ; enddo
call body%query_many(points=q, distances=d, closest_points=cp, normals=nm, principal_curvatures=k, principal_directions=dirs)
do m = 1, NA * NB
   c = cp(:, m) + 0.01_R8P * nm(:, m)          ! lifted off the facets
   e1 = 0.09_R8P * dirs(:, 1, m) ; e2 = 0.09_R8P * dirs(:, 2, m)
   a1(m) = vector_R8P(c(1) - e1(1), c(2) - e1(2), c(3) - e1(3)) ; b1(m) = vector_R8P(c(1) + e1(1), c(2) + e1(2), c(3) + e1(3))
   a2(m) = vector_R8P(c(1) - e2(1), c(2) - e2(2), c(3) - e2(3)) ; b2(m) = vector_R8P(c(1) + e2(1), c(2) + e2(2), c(3) + e2(3))
enddo
print '(A,I0,A,F7.3,A,F7.3,A,F7.3,A,F7.3)', 'feet ', NA * NB, ': k1 from ', minval(k(1, :)), ' to ', maxval(k(1, :)), &
      ', k2 from ', minval(k(2, :)), ' to ', maxval(k(2, :))
call fig_open(u, 'ck-principal-curvature', 'principal directions at the feet: d1 (larger curvature) ember, d2 teal', &
              -55._R8P, 28._R8P)
call fig_surface(u, 'body', body)
call fig_segments(u, 'patch', a1, b1)
call fig_segments(u, 'probe', a2, b2)
call fig_close(u)
endprogram fig_ck_principal_curvature
