!< Figure of the cookbook recipe "Moments and quadrature of cut cells": the disc of ck_cut_moments on 16 x 16 cells,
!< seen from above and zoomed on cell (4, 8): the cells coloured by their fluid fraction alpha, the disc's polygon, and
!< the quadrature points of the pieces of the 3 x 3 cells around it (signed triangles from the cell centre, 7 points
!< each, so a point can fall in the solid with a negative weight; a fluid cell takes the 3 x 3 Gauss rule).
program fig_ck_cut_moments
use fossil,     only : curve_object, cut_cells_2d_t
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curve_object)            :: disc
type(cut_cells_2d_t)          :: cut
real(R8P)                     :: u(48), v(48), x(17)
type(vector_R8P), allocatable :: a(:), b(:)
integer(I4P)                  :: i, j, p, q, un

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(0.0625_R8P * i, i = 0, 16)]
call disc%cut_cells(x=x, y=x, cut=cut, moments_order=2_I4P, quadrature=.true.)
call fig_open(un, 'ck-cut-moments', 'quadrature points of the pieces: weight > 0 teal, < 0 ember', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.5_R8P * (x(4) + x(5)), 0.5_R8P * (x(8) + x(9)), 0._R8P), 0.11_R8P)
do j = 1, 16 ; do i = 1, 16
   call fig_quad(un, 'frac', pt(x(i), x(j)), pt(x(i + 1), x(j)), pt(x(i + 1), x(j + 1)), pt(x(i), x(j + 1)), &
                 value=max(cut%alpha(i, j), 0._R8P))
enddo ; enddo
a = [(pt(x(i), 0._R8P, 1.e-3_R8P), i = 1, 17), (pt(0._R8P, x(i), 1.e-3_R8P), i = 1, 17)]
b = [(pt(x(i), 1._R8P, 1.e-3_R8P), i = 1, 17), (pt(1._R8P, x(i), 1.e-3_R8P), i = 1, 17)]
call fig_segments(un, 'box', a, b)
a = [(pt(u(i), v(i), 2.e-3_R8P), i = 1, 48)] ; b = [(pt(u(mod(i, 48) + 1), v(mod(i, 48) + 1), 2.e-3_R8P), i = 1, 48)]
call fig_segments(un, 'accent', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do p = 1, cut%np
   if (any(abs(cut%piece_cell(:, p) - [4, 8]) > 1)) cycle
   do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
      if (cut%quad_weight(q) > 0._R8P) then
         a = [a, pt(cut%quad_point(1, q), cut%quad_point(2, q), 3.e-3_R8P)]
      else
         b = [b, pt(cut%quad_point(1, q), cut%quad_point(2, q), 3.e-3_R8P)]
      endif
   enddo
enddo
call fig_points(un, 'probe', a)
call fig_points(un, 'outside', b)
call fig_close(un)
contains
   function pt(xx, yy, dz) result(r)
   !< A point of the plane, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: xx, yy !< Coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: r      !< Point.

   r = vector_R8P(xx, yy, 0._R8P)
   if (present(dz)) r%z = dz
   endfunction pt
endprogram fig_ck_cut_moments
