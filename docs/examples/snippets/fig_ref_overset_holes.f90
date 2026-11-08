!< Figure of the reference page "Overset assembly kernels": a lens 0.008 thick at most (chord 1) in a 241 x 41
!< background block, its holes by three kernels: the exact flags (hole_flags_2d, teal), the holes PEGASUS 5's line of
!< sight misses on a 20 x 20 hole map (ember), and the vertices where X-rays along v through 41 rays differ from the
!< exact flags (gold). The vertical scale is stretched 6 times; the rear half is shown.
program fig_ref_overset_holes
use fossil,     only : curvilinear_block_2d_t, curve_object, block_logical_t, hole_flags_2d, hole_map_2d, xray_flags_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P), STRETCH = 6._R8P
type(curvilinear_block_2d_t)       :: bg(1)
type(curve_object)                 :: lens
type(block_logical_t), allocatable :: exact(:), parity(:), pegasus(:), xray(:)
real(R8P)                          :: u(40), v(40), x(241, 41), y(241, 41)
integer(I4P)                       :: i, j, unresolved, w
logical                            :: rear(241 * 41)

do i = 1, 20 ! upper side from the trailing edge to the leading edge, then the lower side back
   u(i) = 1._R8P - (i - 1) / 19._R8P ; v(i) = 0.004_R8P * sin(PI * u(i)) + 0.0001_R8P
enddo
do i = 21, 40
   u(i) = (i - 21) / 19._R8P ; v(i) = -0.004_R8P * sin(PI * u(i)) - 0.0001_R8P
enddo
call lens%adopt_loops(u=u, v=v, loop_sizes=[40_I4P])
do j = 1, 41
   do i = 1, 241
      x(i, j) = -0.1_R8P + 1.2_R8P * (i - 1) / 240 ; y(i, j) = -0.0503_R8P + 0.1006_R8P * (j - 1) / 40
   enddo
enddo
call bg(1)%initialize(x, y)
call hole_flags_2d(bg, lens, exact)
call hole_map_2d(bg, lens, parity, [20_I4P, 20_I4P])
call hole_map_2d(bg, lens, pegasus, [20_I4P, 20_I4P], pegasus=.true., unresolved=unresolved)
call xray_flags_2d(bg, lens, xray, 41_I4P)
print '(A,I0)', 'exact holes                       ', count(exact(1)%v)
print '(A,I0)', 'hole map, crossing parity: differ ', count(parity(1)%v .neqv. exact(1)%v)
print '(A,I0,A,I0)', 'hole map, PEGASUS 5: holes missed ', count(exact(1)%v .and. .not. pegasus(1)%v), &
                     ', unresolved ', unresolved
print '(A,I0)', 'X-rays (41 rays): differ          ', count(xray(1)%v .neqv. exact(1)%v)

call fig_open(w, 'ref-overset-holes', 'lens: holes (teal), missed by PEGASUS 5 (ember), X-ray errors (gold)', &
              -90._R8P, 90._R8P, height=260_I4P)
rear = reshape(x, [241 * 41]) >= 0.63_R8P .and. abs(reshape(y, [241 * 41])) < 0.012_R8P ! the rear half, near the lens
call fig_segments(w, 'accent', [(p(u(i), v(i)), i = 1, 7)], [(p(u(i + 1), v(i + 1)), i = 1, 7)])
call fig_segments(w, 'accent', [(p(u(i), v(i)), i = 33, 39)], [(p(u(i + 1), v(i + 1)), i = 33, 39)])
call fig_segments(w, 'accent', [p(u(40), v(40))], [p(u(1), v(1))])
call fig_points(w, 'box', pack(grid(), rear .and. .not. exact(1)%v))
call fig_points(w, 'probe', pack(grid(), rear .and. exact(1)%v))
call fig_points(w, 'hole', pack(grid(), rear .and. exact(1)%v .and. .not. pegasus(1)%v))
call fig_points(w, 'foot', pack(grid(), rear .and. (xray(1)%v .neqv. exact(1)%v)))
call fig_close(w)
contains
   function p(a, b) result(q)
   !< A point of the plane, its v stretched.
   real(R8P), intent(in) :: a, b !< Coordinates.
   type(vector_R8P)      :: q    !< Point.

   q = vector_R8P(a, STRETCH * b, 0._R8P)
   endfunction p

   function grid() result(q)
   !< The background vertices, linear index (i fastest).
   type(vector_R8P) :: q(241 * 41) !< Points.
   integer(I4P)     :: k, l

   q = [((p(x(k, l), y(k, l)), k = 1, 241), l = 1, 41)]
   endfunction grid
endprogram fig_ref_overset_holes
