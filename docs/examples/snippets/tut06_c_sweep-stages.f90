use fossil,     only : surface_stl_object, distance_field_object, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter     :: STOPS(3) = [1, 4, 8]
real(R8P),    parameter     :: BIG = 1.e300_R8P             ! above it: no magnitude yet
type(surface_stl_object)    :: surface
type(distance_field_object) :: field, banded, signed, stopped(3)
real(R8P), allocatable      :: x(:), lo(:), hi(:)
real(R8P)                   :: h, xs
integer(I4P)                :: i, m, n, s, sweeps, u
character(96)               :: caption(3)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! the grid of the first run
h = 0.04_R8P
do m = 1, 6
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x)
call field%compute_band(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER)   ! exact, |d| <= band
print '(A,I0,A)', 'compute_band:  ', field%get_band_cells(), ' band cells'
banded = field
call field%compute_signs(surface=surface)                         ! flood fill: a sign for every other node
print '(A,I0,A)', 'compute_signs: ', count(field%phi < 0._R8P), ' nodes inside'
signed = field
call field%sweep(sweeps=sweeps)                                   ! fast sweeping to convergence
print '(A,I0,A)', 'sweep:         ', sweeps, ' directional sweeps to converge'
do s = 1, size(STOPS)                                             ! the signed field swept a few times only
   stopped(s) = signed
   call stopped(s)%sweep(max_sweeps=STOPS(s))
   write(caption(s), '(A,I2,A,I0,A,ES8.1)') 'stopped after ', STOPS(s), ' sweeps: ', &
                                            count(abs(stopped(s)%phi) > BIG), ' nodes unknown, error ', &
                                            maxval(abs(stopped(s)%phi - field%phi), abs(stopped(s)%phi) <= BIG)
   print '(A)', trim(caption(s))
enddo
