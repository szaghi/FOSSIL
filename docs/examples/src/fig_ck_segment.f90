!run figck-segment fig_ck_segment
!fig ck-segment
!< Figure of the cookbook recipe "Segment a part by its shape diameter": the marching-cubes dumbbell, its facets
!< coloured by the label segment_sdf gives them with two clusters: the thin bar (label 1) in ember, the balls (label 2)
!< in the body colour.
program fig_ck_segment
use fossil,     only : surface_stl_object, facet_object, extract_isosurface
use fossil_fig, only : fig_open, fig_triangle, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)        :: dumbbell
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: g
integer(I4P), allocatable       :: labels(:)
real(R8P)                       :: f(46, 21, 21), x(46), y(21)
integer(I4P)                    :: i, j, k, u
character(80)                   :: caption

x = [(-1.8_R8P + 0.08_R8P * i, i = 0, 45)] ; y = [(-0.8_R8P + 0.08_R8P * i, i = 0, 20)]
do k = 1, 21 ; do j = 1, 21 ; do i = 1, 46
   f(i, j, k) = min(norm2([x(i) + 1._R8P, y(j), y(k)]) - 0.6_R8P, norm2([x(i) - 1._R8P, y(j), y(k)]) - 0.6_R8P, &
                    max(norm2([y(j), y(k)]) - 0.2_R8P, abs(x(i)) - 1._R8P))
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), y(1), y(1)), bmax=vector_R8P(x(46), y(21), y(21)), &
                        iso=0._R8P, surface=facets)
call dumbbell%adopt_facets(facets=facets)
call dumbbell%segment_sdf(facet_labels=labels, num_clusters=2_I4P)
write(caption, '(A,I0,A,I0,A)') 'segment_sdf, 2 clusters: the bar (', count(labels == 1), ' facets) and the balls (', &
                                count(labels == 2), ')'
call fig_open(u, 'ck-segment', trim(caption), -70._R8P, 25._R8P)
do i = 1, dumbbell%get_facets_number()
   g => dumbbell%facet_at(i)
   call fig_triangle(u, merge('patch', 'body ', labels(i) == 1), g%vertex(1), g%vertex(2), g%vertex(3))
enddo
call fig_close(u)
endprogram fig_ck_segment
