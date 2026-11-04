use fossil, only : surface_stl_object, distance_field_object, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)    :: surface
type(distance_field_object) :: field
real(R8P), allocatable      :: x(:)
real(R8P)                   :: h, d, err_band, err_far
integer(I4P)                :: i, j, k, m, ni, nj, nk, sweeps, fallback, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! uniform core: 31 points, spacing 0.04
h = 0.04_R8P
do m = 1, 6                                                     ! 6 points per side, spacing growing by 1.2
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x, status=status)
