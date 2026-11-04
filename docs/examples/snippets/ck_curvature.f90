program ck_curvature
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(41, 41, 41), d(1), hk(2, 1)
integer(I4P)                    :: i, j, k, status

! a sphere of radius 1: marching cubes on the zero level of |x| - 1, sampled on [-1.5, 1.5]^3
do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
! mean (H) and Gaussian (K) curvature at the foot of a query point
call sphere%query_many(points=[vector_R8P(0.3_R8P, 0.4_R8P, 2._R8P)], distances=d, curvatures=hk)
print '(A,I0)',   'facets ', sphere%get_facets_number()
print '(A,F8.4,A)', 'H ', hk(1, 1), '  (exact 1)'
print '(A,F8.4,A)', 'K ', hk(2, 1), '  (exact 1)'
endprogram ck_curvature
