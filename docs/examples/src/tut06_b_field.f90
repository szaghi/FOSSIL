!run tut06 tut06_b_field
!< Tutorial 6, second run: the signed-distance field of the dragon on the stretched grid, in one compute call.
program tut06_b_field
use fossil, only : surface_stl_object, distance_field_object, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface
type(distance_field_object) :: field
real(R8P), allocatable      :: x(:)
real(R8P)                   :: h
integer(I4P)                :: i, m, ni, nj, nk, sweeps, fallback, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! the grid of the first run
h = 0.04_R8P
do m = 1, 6
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x, status=status)
!region compute
call field%compute(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER, sweeps=sweeps, fallback_cells=fallback, &
                   status=status)
if (status /= STATUS_OK) error stop 'distance field failed'
call field%get_dims(ni=ni, nj=nj, nk=nk)
print '(A,I0,A,I0,A,I0)', 'grid         ', ni, ' x ', nj, ' x ', nk
print '(A,F8.4,A,I0,A)',  'band         ', field%get_band_width(), ' half-width, ', field%get_band_cells(), ' cells'
print '(A,I0,A,I0)',      'sweeps       ', sweeps, ', point-query signs ', fallback
print '(A,I0)',           'inside cells ', count(field%phi < 0._R8P)
!endregion compute
endprogram tut06_b_field
