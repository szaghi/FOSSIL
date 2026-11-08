use fossil, only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(cut_cells_3d_t)     :: cut
type(cut_cells_2d_t)     :: cut2
real(R8P)                :: x(25), h
integer(I4P)             :: i, c(3), p, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)    ! chapter 3
h = 0.05_R8P
x = [(-0.6_R8P + h * i, i = 0, 24)]                  ! 24^3 cells, h = 0.05; z = 0 is the plane x(13)
call cut_cells_3d(surface=surface, x=x, y=x, z=x, cut=cut)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells: fluid ', count(cut%alpha == 1._R8P), ', solid ', count(cut%alpha == 0._R8P), &
                                ', cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P), &
                                ', flagged ', count(cut%pieces < 0)
