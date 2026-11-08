program ref_isotropic_remesh
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, REM_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: cube, sphere, input
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(25, 25, 25), lower, upper
integer(I4P)                    :: i, j, k, status
logical                         :: features
do k = 0, 1
   features = k == 1
   call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call cube%isotropic_remesh(target_length=0.2_R8P, iterations=5_I4P, preserve_features=features, status=status)
   if (status /= REM_STATUS_OK) error stop 'remesh failed'
   print '(A,L1,A,I4,A,F6.4,A,F6.4,A,F7.5,A,F7.5)', 'cube, preserve_features ', features, ': ', &
         cube%get_facets_number(), ' facets, edges ', cube%smallest_edge_len(), ' to ', cube%largest_edge_len(), &
         ', area ', cube%get_area(), ', volume ', cube%get_volume()
enddo
do k = 1, 25 ; do j = 1, 25 ; do i = 1, 25                       ! a unit sphere by marching cubes
   phi(i, j, k) = norm2(-1.5_R8P + 0.125_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call input%adopt_facets(facets=facets)
call input%sanitize(fill_holes=.true.)
print '(A,I4,A,F6.4,A,F6.4,A,F7.5)', 'sphere input:  ', input%get_facets_number(), ' facets, edges ', &
      input%smallest_edge_len(), ' to ', input%largest_edge_len(), ', volume ', input%get_volume()
sphere = input
do k = 1, 3                                                       ! one iteration per call, to watch it
   call sphere%isotropic_remesh(target_length=0.1_R8P, iterations=1_I4P, status=status)
   print '(A,I1,A,I4,A,F6.4,A,F6.4,A,F7.5)', 'iteration ', k, ':   ', sphere%get_facets_number(), &
         ' facets, edges ', sphere%smallest_edge_len(), ' to ', sphere%largest_edge_len(), ', volume ', sphere%get_volume()
enddo
call sphere%hausdorff(other=input, lower=lower, upper=upper)     ! certified distance to the input
print '(A,F6.4,A,F6.4)', 'Hausdorff distance to the input ', lower, ' .. ', upper
endprogram ref_isotropic_remesh
