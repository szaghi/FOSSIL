program ck_markers
use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
real(R8P), allocatable   :: centroids(:,:), areas(:), normals(:,:)
real(R8P)                :: fraction
integer(I4P)             :: status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%surface_markers(dx=0.1_R8P, centroids=centroids, areas=areas, normals=normals, spacing_fraction=fraction, &
                          preserve_features=.true., status=status)        ! for a grid of spacing 0.1
if (status /= STATUS_OK) error stop 'surface_markers failed'
print '(A,I0)', 'markers ', size(areas)
print '(A,F10.6,A,ES10.3)', 'sum of the areas ', sum(areas), ', mean area ', sum(areas) / size(areas)
print '(A,F7.4)', 'spacings within [0.8, 1.2] dx: ', fraction
print '(A,3F8.4,A,3F5.1)', 'marker 1 at ', centroids(:, 1), ', normal ', normals(:, 1)
endprogram ck_markers
