program ref_decimate
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, DEC_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: fine, coarse
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(41, 41, 41), lower, upper
integer(I4P)                    :: i, j, k, target, status
do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41                     ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call fine%adopt_facets(facets=facets)
call fine%sanitize(fill_holes=.true.)
print '(A,I5,A,F8.5,A)', 'input:          ', fine%get_facets_number(), ' facets, volume ', fine%get_volume(), &
                         ' (sphere: 4.18879)'
coarse = fine
target = fine%get_facets_number()
do k = 1, 3
   target = target / 4
   call coarse%decimate(target_facets=target, status=status)    ! in place: each level from the previous one
   if (status /= DEC_STATUS_OK) error stop 'decimate did not reach the target'
   call coarse%hausdorff(other=fine, lower=lower, upper=upper)    ! certified distance to the input
   print '(A,I5,A,I5,A,F8.5,A,L1,A,F6.4,A,F6.4)', 'target ', target, ': ', coarse%get_facets_number(), &
         ' facets, volume ', coarse%get_volume(), ', closed manifold ', coarse%is_watertight() .and. &
         coarse%is_manifold(), ', Hausdorff ', lower, ' .. ', upper
enddo
endprogram ref_decimate
