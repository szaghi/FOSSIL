!run tut06 tut06_distance_field
!< Tutorial 6: the signed-distance field of the dragon on a stretched rectilinear grid.
program tut06_distance_field
!region grid
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
!endregion grid
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
!region check
err_band = 0._R8P ; err_far = 0._R8P
do k = 1, nk ; do j = 1, nj ; do i = 1, ni
   d = surface%distance(point=x(i) * ex_R8P + x(j) * ey_R8P + x(k) * ez_R8P, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   if (field%in_band(i, j, k)) then
      err_band = max(err_band, abs(field%phi(i, j, k) - d))
   else
      err_far = max(err_far, abs(field%phi(i, j, k) - d))
   endif
enddo ; enddo ; enddo
print '(A,ES10.3)', 'max error in the band ', err_band
print '(A,ES10.3)', 'max error elsewhere   ', err_far
!endregion check
endprogram tut06_distance_field
