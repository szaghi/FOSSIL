!run ck_hausdorff ck_hausdorff
program ck_hausdorff
use fossil, only : surface_stl_object, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: wing, turned
real(R8P)                :: lower, upper
integer(I4P)             :: status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
turned = wing                                                 ! a second version: 2 degrees more incidence
call turned%rotate(axis=ez_R8P, angle=-2._R8P * PI / 180._R8P, center=vector_R8P(0._R8P, 0._R8P, 0._R8P))
call turned%analyze
call wing%hausdorff(other=turned, lower=lower, upper=upper, tolerance=1e-6_R8P, status=status)
if (status /= DEV_STATUS_OK) error stop 'hausdorff did not converge'
print '(A,F10.7,A,F10.7,A)', 'Hausdorff distance in [', lower, ',', upper, ']'
print '(A,F10.7)', 'trailing edge moved     ', 0.75_R8P * 2._R8P * sin(PI / 180._R8P)
endprogram ck_hausdorff
