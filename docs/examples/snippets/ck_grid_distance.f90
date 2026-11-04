program ck_grid_distance
use fossil, only : distance_field_object, surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)    :: cube
type(distance_field_object) :: field
real(R8P)                   :: x(41)
integer(I4P)                :: i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
x = [(-0.5_R8P + 0.05_R8P * i, i = 0, 40)] ! 41 nodes on [-0.5, 1.5]
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=cube, status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
print '(A,I0)',   'band cells        ', field%get_band_cells()
print '(A,F8.4)', 'phi at the centre ', field%phi(21, 21, 21)
print '(A,F8.4)', 'phi at a corner   ', field%phi(1, 1, 1)
endprogram ck_grid_distance
