use fossil,     only : surface_stl_object, distance_field_object, curve_object, slice_surface, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)    :: surface
type(distance_field_object) :: field
type(curve_object)          :: section
real(R8P), allocatable      :: x(:)
real(R8P)                   :: h, uv(2)
type(vector_R8P)            :: a, b
integer(I4P)                :: i, m, ni, nj, nk, status, u, s
character(100)              :: caption

x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! uniform core: 31 points, spacing 0.04
h = 0.04_R8P
do m = 1, 6                                                     ! 6 points per side, spacing growing by 1.2
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x, status=status)
if (status /= STATUS_OK) error stop 'not a valid grid'
call field%get_dims(ni=ni, nj=nj, nk=nk)
print '(A,I0,A,I0,A,I0)', 'grid ', ni, ' x ', nj, ' x ', nk
print '(A,2F9.4)',        'x from, to       ', x(1), x(size(x))
print '(A,2F9.4)',        'spacing core, end ', x(8) - x(7), x(2) - x(1)
