!< Figure of the reference page "Body-fitted grid support", section The wall-spacing report: an O-grid of 64 x 12 cells
!< around the unit circle, first height 0.04 (1 + 0.5 cos theta), growth 1.15, its first layers sheared on the upper half
!< so they leave the wall at an angle. `block_report_2d` with u_tau = 1 and nu = 0.04 (y1+ = 1 for a height 0.04).
!< Top of the circle, zoomed: the first-layer cells coloured by the y1+ of their wall vertex (teal below the target 1,
!< ember above), and for each wall vertex its first grid-line step (gold) and the distance of the same off-wall point
!< to the wall (teal), which the report uses for y1+.
program fig_ref_body_fitted_report
use fossil,     only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                       block_report_t, block_report_2d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter      :: PI = acos(-1._R8P)
integer(I4P), parameter      :: NI = 65, NJ = 13
type(curvilinear_block_2d_t) :: grid(1)
type(topo_record_t)          :: rec
type(curve_object)           :: wall
type(block_report_t)         :: report
real(R8P)                    :: u(NI, NJ), v(NI, NJ), t, r, h, a
type(vector_R8P)             :: s1(NI - 1), s2(NI - 1), f(NI - 1)
integer(I4P)                 :: i, j, k, w, m

do i = 1, NI
   t = 2 * PI * (i - 1) / (NI - 1)
   h = 0.04_R8P * (1._R8P + 0.5_R8P * cos(t)) ; a = 1.5_R8P * max(0._R8P, sin(t))**2
   r = 1._R8P
   do j = 1, NJ
      u(i, j) = r * cos(t + a * min(r - 1._R8P, 0.15_R8P)) ; v(i, j) = r * sin(t + a * min(r - 1._R8P, 0.15_R8P))
      r = r + h * 1.15_R8P**(j - 1)
   enddo
enddo
u(NI, :) = u(1, :) ; v(NI, :) = v(1, :)
call grid(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)
call wall_cutter_2d(grid, wall)
call block_report_2d(grid(1), wall, 0.04_R8P, report, u_tau=1._R8P, target_y1p=1._R8P)
print '(A,I0,A,3F7.4,A,F6.3)', 'wall vertices ', report%nw, ', y1+ min mean max ', report%y1p_min, report%y1p_mean, &
                               report%y1p_max, ', fraction above 1 ', report%above
k = maxloc(report%y1_line / report%y1_wall, 1)
print '(A,I0,A,F8.5,A,F8.5,A,F7.4)', 'largest ratio at wall vertex ', k, ': y1 along the line ', report%y1_line(k), &
                                     ', to the wall ', report%y1_wall(k), ', y1+ ', report%y1p(k)
k = 1
print '(A,I0,A,F8.5,A,F8.5,A,F7.4)', 'unsheared wall vertex ', k, ':  y1 along the line ', report%y1_line(k), &
                                     ', to the wall ', report%y1_wall(k), ', y1+ ', report%y1p(k)
print '(A,F6.3,A,F7.4,A,I0)', 'growth max ', report%growth_max, ', orthogonality min ', report%orth_min, &
                              ', cells not valid ', report%not_valid

call fig_open(w, 'ref-body-fitted-report', 'first layer by y1+ (teal < 1 < ember); gold: grid-line step; '// &
              'teal: wall distance', -90._R8P, 90._R8P, width=760_I4P, height=460_I4P)
call fig_zoom(w, vector_R8P(0.06_R8P, 1.03_R8P, 0._R8P), 0.17_R8P)
do m = 1, report%nw
   i = report%point(1, m)
   if (i >= NI) cycle
   call fig_quad(w, 'field', p(i, 1), p(i + 1, 1), p(i + 1, 2), p(i, 2), max(-1._R8P, min(1._R8P, 2 * (report%y1p(m) - 1))))
enddo
call fig_segments(w, 'box', [((p(i, j), i = 1, NI - 1), j = 1, NJ)], [((p(i + 1, j), i = 1, NI - 1), j = 1, NJ)])
call fig_segments(w, 'box', [((p(i, j), i = 1, NI), j = 2, NJ - 1)], [((p(i, j + 1), i = 1, NI), j = 2, NJ - 1)])
do i = 1, NI - 1   ! the first step along the grid line, and the radial foot of the same off-wall point
   s1(i) = p(i, 1) + vector_R8P(0._R8P, 0._R8P, 0.002_R8P) ; s2(i) = p(i, 2) + vector_R8P(0._R8P, 0._R8P, 0.002_R8P)
   r = sqrt(u(i, 2)**2 + v(i, 2)**2) ; f(i) = vector_R8P(u(i, 2) / r, v(i, 2) / r, 0.002_R8P)
enddo
call fig_segments(w, 'accent', s1, s2)
call fig_segments(w, 'probe', f, s2)
call fig_close(w)
contains
   function p(i, j) result(q)
   !< Vertex (i, j), lifted a little.
   integer(I4P), intent(in) :: i, j !< Indices.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(u(i, j), v(i, j), 0.001_R8P)
   endfunction p
endprogram fig_ref_body_fitted_report
