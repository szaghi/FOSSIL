!run ck_feature_boxes ck_feature_boxes
program ck_feature_boxes
use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
integer(I4P), parameter  :: N = 13                    ! 13^3 cells of side 0.1 over the unit cube
type(surface_stl_object) :: cube
real(R8P)                :: lo(3, N**3), hi(3, N**3)
logical                  :: touched(N**3)
integer(I4P)             :: i, j, k, c, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%build_features(angle=30._R8P, status=status)  ! the sharp edges (dihedral above 30 degrees), in a tree
if (status /= STATUS_OK) error stop 'build_features failed'
c = 0
do k = 1, N ; do j = 1, N ; do i = 1, N
   c = c + 1 ; lo(:, c) = -0.13_R8P + 0.1_R8P * ([i, j, k] - 1) ; hi(:, c) = lo(:, c) + 0.1_R8P
enddo ; enddo ; enddo
call cube%box_feature_many(lo=lo, hi=hi, touched=touched, status=status)
print '(A,I0,A,I0)', 'cells ', N**3, ', touching a sharp edge ', count(touched)
endprogram ck_feature_boxes
