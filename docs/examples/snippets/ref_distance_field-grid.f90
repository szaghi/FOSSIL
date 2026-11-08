use fossil,     only : surface_stl_object, distance_field_object, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)    :: cube
type(distance_field_object) :: field
real(R8P), allocatable      :: x(:), err(:,:,:), q(:,:,:)
real(R8P)                   :: h, e(3), exact, band_err, far_err, foot_err
integer(I4P)                :: i, j, k, m, n, sweeps, fallback, flagged, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! [0, 1]^3
! each axis: 29 nodes with spacing 0.05 on [-0.2, 1.2], then 12 per side growing by 1.15
x = [(-0.2_R8P + 0.05_R8P * i, i = 0, 28)]
h = 0.05_R8P
do m = 1, 12
   h = 1.15_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
n = size(x)
call field%initialize(x=x, y=x, z=x, status=status)
