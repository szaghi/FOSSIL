use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)       :: section
type(surface_stl_object) :: surface
type(cut_cells_2d_t)     :: cut2
type(cut_cells_3d_t)     :: cut3
real(R8P)                :: x(49), h
integer(I4P)             :: i, k, c(3), p, status

call section%load_from_file(file_name='dragon-section.dat', status=status)
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]                  ! cell edges: 48 cells per direction, h = 0.025
h = 0.025_R8P
call section%cut_cells(x=x, y=x, cut=cut2)
print '(A,I0,A,I0)', '2-D: cut cells ', count(cut2%alpha > 0._R8P .and. cut2%alpha < 1._R8P), &
                     ', cells with 2+ pieces ', count(cut2%pieces > 1)
print '(A,F10.6,A,F10.6)', '     solid area ', sum(1._R8P - cut2%alpha) * h**2, '  curve area ', section%get_area()
c(1:2) = maxloc(cut2%wet_length)
print '(A,I0,A,I0,A,F8.5,A,4F8.5)', '     cell (', c(1), ',', c(2), ') alpha ', cut2%alpha(c(1), c(2)), &
                                    '  apertures ', cut2%aperture(:, c(1), c(2))
