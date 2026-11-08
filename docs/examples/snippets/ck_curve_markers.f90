program ck_curve_markers
use fossil, only : surface_stl_object, curve_object, slice_surface, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: wing
type(curve_object)       :: section
real(R8P), allocatable   :: centroids(:,:), lengths(:), normals(:,:)
real(R8P)                :: fraction
integer(I4P)             :: status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)                    ! the airfoil section
call section%curve_markers(dx=0.05_R8P, centroids=centroids, lengths=lengths, normals=normals, &
                           spacing_fraction=fraction, status=status)  ! for a grid of spacing 0.05
if (status /= STATUS_OK) error stop 'curve_markers failed'
print '(A,I0,A,F9.6,A,F9.6)', 'markers ', size(lengths), ', each of length ', lengths(1), ', total ', sum(lengths)
print '(A,F7.4)', 'chord spacings within [0.8, 1.2] dx: ', fraction
print '(A,2F8.4,A,2F8.4)', 'marker 1 at ', centroids(:, 1), ', normal ', normals(:, 1)
endprogram ck_curve_markers
