!run figck-curvature fig_ck_curvature
!fig ck-curvature
!< Figure of the recipe "Measure the curvature at a point": the marching-cubes unit sphere of the recipe, each facet
!< coloured by the error H - 1 of the mean curvature query_many(curvatures=) returns at its centroid's foot (queried
!< radially from 2 % outside; teal below 1, ember above, saturating at 0.05, about the 99th percentile), the facets
!< where |H - 1| > 0.1 marked by a dot at their centroid (the slivers marching cubes leaves near grid nodes: the jet
!< fit is poor there, |H - 1| up to ~4), and the recipe's query point with its foot.
program fig_ck_curvature
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_segment, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P),   allocatable :: q(:)
real(R8P),          allocatable :: dq(:), hq(:, :)
type(vector_R8P)                :: p(1), foot(1)
real(R8P),          parameter   :: SCALE = 0.05_R8P, OUTLIER = 0.1_R8P
real(R8P)                       :: phi(41, 41, 41), d(1), hk(2, 1), cp(3, 1)
character(80)                   :: caption
integer(I4P)                    :: u, i, j, k, n, status

do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41                                        ! the recipe's sphere
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
p = vector_R8P(0.3_R8P, 0.4_R8P, 2._R8P)                                           ! the recipe's query
call sphere%query_many(points=p, distances=d, closest_points=cp, curvatures=hk)
foot(1)%x = cp(1, 1) ; foot(1)%y = cp(2, 1) ; foot(1)%z = cp(3, 1)
n = sphere%get_facets_number()                                                     ! H at every facet's foot
allocate(q(n), dq(n), hq(2, n))
do i = 1, n
   f => sphere%facet_at(i)
   q(i) = f%centroid * (1.02_R8P / f%centroid%normL2())     ! radially, 2 % outside
enddo
call sphere%query_many(points=q, distances=dq, curvatures=hq)
write(caption, '(A,F6.4,A,F4.2,A,F3.1,A,I0,A)') 'foot H = ', hk(1, 1), '; colour H - 1 to +-', SCALE, &
                                               ', dots |H - 1| > ', OUTLIER, ' (', count(abs(hq(1, :) - 1._R8P) > OUTLIER), ')'
call fig_open(u, 'ck-curvature', trim(caption), -60._R8P, 30._R8P)
do i = 1, n
   f => sphere%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), value=(hq(1, i) - 1._R8P) / SCALE)
enddo
call fig_points(u, 'accent', pack(q / 1.02_R8P, abs(hq(1, :) - 1._R8P) > OUTLIER))
call fig_segment(u, 'foot', p(1), foot(1))
call fig_points(u, 'probe', p)
call fig_points(u, 'foot', foot)
call fig_close(u)
endprogram fig_ck_curvature
