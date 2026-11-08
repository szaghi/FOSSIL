!< Figure of the cookbook recipe "Measure how much of a surface surrounds a point": the winding number w of the open
!< box (the unit cube without its top face) on the plane x = 0.5, seen side-on, coloured by 2 w - 1 (teal: w near 1,
!< ember: w near 0); the section of the box, open at the top, and the recipe's four points.
program fig_ck_winding_number
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter  :: n = 60
type(surface_stl_object) :: box
real(R8P)                :: y(0:n), z(0:n), w, h, zq(4) = [0.5_R8P, 0.9_R8P, 1.1_R8P, 2._R8P]
integer(I4P)             :: i, j, u

call box%load_from_file(file_name='cube-01.stl', guess_format=.true., &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), clip_max=vector_R8P(2._R8P, 2._R8P, 0.9_R8P))
h = 2.5_R8P / n
y = [(-0.7314_R8P + h * i, i = 0, n)] ; z = [(-0.2814_R8P + h * i, i = 0, n)]   ! no cell centre on a wall
call fig_open(u, 'ck-winding-number', 'winding number of the open box on the plane x = 0.5: 1 teal, 0 ember', &
              0._R8P, 0._R8P, height=520_I4P)
do j = 0, n - 1 ; do i = 0, n - 1
   w = box%winding_number(point=vector_R8P(0.5_R8P, 0.5_R8P * (y(i) + y(i + 1)), 0.5_R8P * (z(j) + z(j + 1))))
   call fig_quad(u, 'field', q(y(i), z(j)), q(y(i + 1), z(j)), q(y(i + 1), z(j + 1)), q(y(i), z(j + 1)), &
                 value=1._R8P - 2._R8P * w)
enddo ; enddo
call fig_segments(u, 'accent', [q(0._R8P, 1._R8P), q(0._R8P, 0._R8P), q(1._R8P, 0._R8P)], &
                               [q(0._R8P, 0._R8P), q(1._R8P, 0._R8P), q(1._R8P, 1._R8P)])
call fig_points(u, 'probe', [(q(0.5_R8P, zq(i)), i = 1, 4)])
do i = 1, 4
   w = box%winding_number(point=vector_R8P(0.5_R8P, 0.5_R8P, zq(i)))
   call fig_label(u, 'accent', q(0.56_R8P, zq(i)), 'w '//trim(fmt(w)))
enddo
call fig_close(u)
contains
   function q(yy, zz) result(p)
   !< A point of the plane x = 0.5, lifted towards the viewer (+x).
   real(R8P), intent(in) :: yy, zz !< In-plane coordinates.
   type(vector_R8P)      :: p      !< Point.

   p = vector_R8P(0.6_R8P, yy, zz)
   endfunction q

   function fmt(v) result(s)
   !< A value with 4 decimals.
   real(R8P), intent(in) :: v !< Value.
   character(8)          :: s !< Text.

   write(s, '(F6.4)') v
   endfunction fmt
endprogram fig_ck_winding_number
