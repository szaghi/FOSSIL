program ck_segment
use fossil, only : surface_stl_object, facet_object, extract_isosurface, SDF_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: dumbbell
type(facet_object), allocatable :: facets(:)
integer(I4P), allocatable       :: labels(:)
real(R8P),    allocatable       :: sdf(:)
real(R8P)                       :: f(46, 21, 21), x(46), y(21)
integer(I4P)                    :: i, j, k, status

x = [(-1.8_R8P + 0.08_R8P * i, i = 0, 45)] ; y = [(-0.8_R8P + 0.08_R8P * i, i = 0, 20)]
do k = 1, 21 ; do j = 1, 21 ; do i = 1, 46     ! two balls of radius 0.6 joined by a bar of radius 0.2
   f(i, j, k) = min(norm2([x(i) + 1._R8P, y(j), y(k)]) - 0.6_R8P, norm2([x(i) - 1._R8P, y(j), y(k)]) - 0.6_R8P, &
                    max(norm2([y(j), y(k)]) - 0.2_R8P, abs(x(i)) - 1._R8P))
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(x(1), y(1), y(1)), bmax=vector_R8P(x(46), y(21), y(21)), &
                        iso=0._R8P, surface=facets)
call dumbbell%adopt_facets(facets=facets)
call dumbbell%segment_sdf(facet_labels=labels, sdf=sdf, num_clusters=2_I4P, status=status)
if (status /= SDF_STATUS_OK) error stop 'segment_sdf failed'
print '(A,I0,A,I0)', 'facets ', size(labels), ', unassigned (label 0) ', count(labels == 0)
do k = 1, 2
   print '(A,I0,A,I5,A,F7.4)', 'label ', k, ': ', count(labels == k), ' facets, mean shape diameter ', &
                               sum(sdf, mask=labels == k) / max(1, count(labels == k))
enddo
endprogram ck_segment
