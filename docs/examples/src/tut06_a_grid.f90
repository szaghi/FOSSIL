!run tut06-grid tut06_a_grid
!fig tut06-grid
!< Tutorial 6, first run: the stretched rectilinear grid of the field: a uniform core over the domain and six nodes per
!< side growing by 1.2. The figure draws the grid lines of the plane x = 0 and the section of the dragon by it.
program tut06_a_grid
use fossil_fig, only : fig_open, fig_segment, fig_close
!region grid
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
!endregion grid
call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 1_I4P, 0._R8P, section)                 ! the plane x = 0: (u, v) = (y, z)
write(caption, '(A,I0,A,I0,A)') 'grid plane x = 0: ', ni, ' x ', nj, ' nodes, uniform core, stretched sides'
call fig_open(u, 'tut06-grid', trim(caption), 0._R8P, 0._R8P, height=500_I4P)
do i = 1, size(x)
   call fig_segment(u, 'box', x(i) * ey_R8P + x(1) * ez_R8P, x(i) * ey_R8P + x(size(x)) * ez_R8P)
   call fig_segment(u, 'box', x(1) * ey_R8P + x(i) * ez_R8P, x(size(x)) * ey_R8P + x(i) * ez_R8P)
enddo
do s = 1, section%get_vertices_number()
   if (section%get_next(s) == 0_I4P) cycle
   uv = section%get_vertex(s) ; a = 1.e-3_R8P * ex_R8P + uv(1) * ey_R8P + uv(2) * ez_R8P   ! over the lines
   uv = section%get_vertex(section%get_next(s)) ; b = 1.e-3_R8P * ex_R8P + uv(1) * ey_R8P + uv(2) * ez_R8P
   call fig_segment(u, 'hole', a, b)
enddo
call fig_close(u)
endprogram tut06_a_grid
