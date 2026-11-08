use fossil, only : curve_object, cut_cells_2d_t, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)   :: section
type(cut_cells_2d_t) :: cut
real(R8P)            :: x(49), h
integer(I4P)         :: i, j, k, c(2), p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
h = 0.025_R8P
x = [(-0.6_R8P + h * i, i = 0, 48)]                  ! the same 48 x 48 block as step 1
call section%cut_cells(x=x, y=x, cut=cut)
print '(A,I0,A,I0,A,I0)', 'cells: fluid ', count(cut%alpha == 1._R8P), ', solid ', count(cut%alpha == 0._R8P), &
                          ', cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P)
print '(A,F10.7,A,F10.7)', 'solid area from alpha ', sum(1._R8P - cut%alpha) * h**2, ', area of the curve ', &
                           section%get_area()
