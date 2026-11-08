!< Figure of the reference page "Cut cells: solver support": the quadrature rules of the 2-D cut cells of a five-pointed
!< star (the solid, filled) on 4 x 4 cells. Volume points of the cut pieces: teal with a positive weight, ember with a
!< negative one (signed triangles from the cell centre, where a piece is not star-shaped about it); wall points gold.
program fig_ref_cut_solver_quadrature
use fossil,     only : curve_object, cut_cells_2d_t
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_triangle, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: PI = acos(-1._R8P)
type(curve_object)            :: star
type(cut_cells_2d_t)          :: cut
real(R8P)                     :: u(10), v(10), x(5), r, err
logical,          allocatable :: cutp(:), keep(:)
type(vector_R8P), allocatable :: pts(:)
integer(I4P)                  :: i, j, k, p, w, n

do k = 1, 10 ! a five-pointed star about (0.5, 0.5), counter-clockwise
   r = merge(0.42_R8P, 0.17_R8P, mod(k, 2) == 1)
   u(k) = 0.5_R8P + r * cos(PI / 2 + 2 * PI * (k - 1) / 10) ; v(k) = 0.5_R8P + r * sin(PI / 2 + 2 * PI * (k - 1) / 10)
enddo
call star%adopt_loops(u=u, v=v, loop_sizes=[10_I4P])
x = [(0.25_R8P * i, i = 0, 4)]
call star%cut_cells(x=x, y=x, cut=cut, quadrature=.true.)
allocate(cutp(cut%np))
cutp = [(cut%piece_area(p) < 0.0625_R8P - 1e-12_R8P .and. cut%piece_area(p) > 0._R8P, p = 1, cut%np)]
n = size(cut%quad_weight)
print '(A,I0,A,I0,A,I0)', 'pieces ', cut%np, ', cut ', count(cutp), ', volume points ', n
print '(A,I0,A,I0)', 'negative volume weights ', count(cut%quad_weight < 0._R8P), ', wall points ', &
                     size(cut%wall_quad_weight)
err = 0._R8P
do p = 1, cut%np
   err = max(err, abs(sum(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1)) - cut%piece_area(p)))
enddo
print '(A,ES9.2)', 'largest |sum of the weights - piece area| ', err

call fig_open(w, 'ref-cut-solver-quadrature', 'quadrature of the cut pieces: weights > 0 teal, < 0 ember; wall gold', &
              -90._R8P, 90._R8P, height=460_I4P)
do k = 1, 10
   call fig_triangle(w, 'body', q(0.5_R8P, 0.5_R8P), q(u(k), v(k)), q(u(mod(k, 10) + 1), v(mod(k, 10) + 1)))
enddo
call fig_segments(w, 'box', [([(q(x(i), x(j)), i = 1, 4)], j = 1, 5)], [([(q(x(i + 1), x(j)), i = 1, 4)], j = 1, 5)])
call fig_segments(w, 'box', [([(q(x(i), x(j)), j = 1, 4)], i = 1, 5)], [([(q(x(i), x(j + 1)), j = 1, 4)], i = 1, 5)])
allocate(keep(n), pts(n)) ; keep = .false.
do p = 1, cut%np
   if (cutp(p)) keep(cut%quad_start(p):cut%quad_start(p + 1) - 1) = .true.
enddo
pts = [(q(cut%quad_point(1, k), cut%quad_point(2, k), 1e-3_R8P), k = 1, n)] ! lifted over the star
call fig_points(w, 'probe', pack(pts, keep .and. cut%quad_weight > 0._R8P))
call fig_points(w, 'hole', pack(pts, keep .and. cut%quad_weight < 0._R8P))
call fig_points(w, 'foot', [(q(cut%wall_quad_point(1, k), cut%wall_quad_point(2, k), 1e-3_R8P), &
                             k = 1, size(cut%wall_quad_weight))])
call fig_close(w)
contains
   function q(a, b, dz) result(pt)
   !< A point of the plane, lifted by dz toward the viewer.
   real(R8P), intent(in)           :: a, b !< Coordinates.
   real(R8P), intent(in), optional :: dz   !< Lift.
   type(vector_R8P)                :: pt   !< Point.

   pt = vector_R8P(a, b, 0._R8P)
   if (present(dz)) pt%z = dz
   endfunction q
endprogram fig_ref_cut_solver_quadrature
