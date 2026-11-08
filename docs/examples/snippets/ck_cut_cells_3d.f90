program ck_cut_cells_3d
use fossil, only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, ez_R8P, vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cube
type(cut_cells_3d_t)     :: cut
real(R8P)                :: x(17)
integer(I4P)             :: i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%rotate(axis=ez_R8P, angle=PI / 6._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))  ! turned off the
call cube%rotate(axis=ex_R8P, angle=PI / 9._R8P, center=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))  ! grid planes
call cube%analyze
x = [(-0.5_R8P + 0.125_R8P * i, i = 0, 16)]                   ! 16^3 cells over [-0.5, 1.5]^3
call cut_cells_3d(cube, x, x, x, cut, status=status)
if (status /= STATUS_OK) error stop 'cut_cells_3d failed'
print '(A,I0,A,I0,A,I0)', 'cells ', cut%nx * cut%ny * cut%nz, ', cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P), &
                          ', solid ', count(cut%alpha == 0._R8P)
print '(A,F16.12)', 'fluid volume         ', sum(cut%alpha) * 0.125_R8P**3
print '(A,F16.12)', 'block - cube volume  ', 8._R8P - cube%get_volume()
print '(A,F16.12)', 'wetted area          ', sum(cut%wet_area)
endprogram ck_cut_cells_3d
