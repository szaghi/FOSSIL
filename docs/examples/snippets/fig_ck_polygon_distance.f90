!< Figure of the cookbook recipe "Compute the signed distance to a polygon": the unit square of ck_polygon_distance
!< (adopt_loops, counter-clockwise), its signed distance (curve%distance at the cell centres, colour tanh(d / 0.8):
!< teal inside, ember outside) and the recipe's three query points joined to their closest points (compute_distance),
!< each labelled with its distance.
program fig_ck_polygon_distance
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_quad, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter :: nc = 40                                ! cells per side
real(R8P),    parameter :: x0 = -0.4_R8P, x1 = 2.4_R8P, lift = 1.e-3_R8P
type(curve_object)      :: square
type(vector_R8P)        :: a(4), b(4), o
real(R8P)               :: q(2, 3), foot(2, 3), d(3), h, c(2)
integer(I4P)            :: i, j, u
character(80)           :: text

call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
q = reshape([0.5_R8P, 0.25_R8P, 2._R8P, 2._R8P, 1.5_R8P, 0.5_R8P], [2, 3])  ! the recipe's three points
do i = 1, 3
   call square%compute_distance(q(:, i), distance=d(i), closest_point=foot(:, i))
enddo

call fig_open(u, 'ck-polygon-distance', 'signed distance to the square, teal inside; three queries, closest points', &
              -90._R8P, 90._R8P, height=520_I4P)
h = (x1 - x0) / nc
do j = 1, nc ; do i = 1, nc
   c = [x0 + (i - 0.5_R8P) * h, x0 + (j - 0.5_R8P) * h]
   call fig_quad(u, 'field', p(c(1) - 0.5_R8P * h, c(2) - 0.5_R8P * h), p(c(1) + 0.5_R8P * h, c(2) - 0.5_R8P * h), &
                 p(c(1) + 0.5_R8P * h, c(2) + 0.5_R8P * h), p(c(1) - 0.5_R8P * h, c(2) + 0.5_R8P * h), &
                 value=tanh(square%distance(c) / 0.8_R8P))
enddo ; enddo
a = [p(0._R8P, 0._R8P, lift), p(1._R8P, 0._R8P, lift), p(1._R8P, 1._R8P, lift), p(0._R8P, 1._R8P, lift)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a, b)                               ! the square
do i = 1, 3
   call fig_segment(u, 'accent', p(q(1, i), q(2, i), lift), p(foot(1, i), foot(2, i), lift))
   write(text, '(F7.4)') d(i)
   call fig_label(u, 'accent', p(q(1, i) + 0.06_R8P, q(2, i) + 0.04_R8P, 2._R8P * lift), trim(adjustl(text)))
enddo
call fig_points(u, 'probe', [(p(q(1, i), q(2, i), 2._R8P * lift), i = 1, 3)])
call fig_points(u, 'foot', [(p(foot(1, i), foot(2, i), 2._R8P * lift), i = 1, 3)])
o = p(x0 + 0.08_R8P, x0 + 0.08_R8P, lift)                          ! axes glyph at a corner
call fig_segment(u, 'accent', o, o + 0.25_R8P * ex_R8P) ; call fig_label(u, 'accent', o + 0.25_R8P * ex_R8P, 'u')
call fig_segment(u, 'accent', o, o + 0.25_R8P * ey_R8P) ; call fig_label(u, 'accent', o + 0.25_R8P * ey_R8P, 'v')
call fig_close(u)
contains
   function p(uu, vv, dz) result(r)
   !< A point of the plane of the curve (z = 0), lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: uu, vv !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: r      !< Point.

   r = uu * ex_R8P + vv * ey_R8P
   if (present(dz)) r = r + dz * ez_R8P
   endfunction p
endprogram fig_ck_polygon_distance
