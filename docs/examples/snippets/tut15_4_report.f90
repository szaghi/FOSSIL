!< Tutorial 15, step 4: the wall-spacing report of the O-grid of step 3 (rebuilt here), with the u_tau the first height
!< was designed with and with the flat-plate estimate at each wall point; the report saved as text. Figure: y1+ along
!< the upper side, chord on the horizontal axis, for the two u_tau.
program tut15_4_report
use fossil,     only : surface_stl_object, curve_object, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, &
                       slice_surface, hyperbolic_grid_2d, geometric_spacing, wall_cutter_2d, block_report_t, &
                       block_report_2d, save_block_report, utau_flat_plate
use fossil_fig, only : fig_open, fig_segments, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)      :: wing
type(curve_object)            :: section, wall
type(curvilinear_block_2d_t)  :: grid(1)
type(topo_record_t)           :: rec
type(block_report_t)          :: report, estimate
real(R8P),        allocatable :: u(:), v(:), orth(:,:)
type(vector_R8P), allocatable :: given(:), flat(:)
real(R8P)                     :: nu, u_inf, u_tau, p(2), q(2)
integer(I4P),     allocatable :: upper(:)
integer(I4P)                  :: i, j, m, first, last, unit, un
character(200)                :: line

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)   ! the O-grid of step 3
call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)
call section%get_loop_range(1_I4P, first, last)
allocate(u(2 * (last - first + 1)), v(2 * (last - first + 1)))
do i = first, last
   p = section%get_vertex(i) ; q = section%get_vertex(merge(first, i + 1_I4P, i == last))
   do j = 0, 1
      u(2 * (i - first) + j + 1) = p(1) + (q(1) - p(1)) * j / 2._R8P
      v(2 * (i - first) + j + 1) = p(2) + (q(2) - p(2)) * j / 2._R8P
   enddo
enddo
nu = 1.e-6_R8P ; u_inf = 1._R8P ; u_tau = utau_flat_plate(u_inf, nu, 0.5_R8P)
call hyperbolic_grid_2d(u, v, geometric_spacing(nu / u_tau, 1.2_R8P, 40_I4P), grid(1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)                                            ! side j = 1 is the wall
call wall_cutter_2d(grid, wall)                                         ! the wall as a curve
call block_report_2d(grid(1), wall, nu, report, u_tau=u_tau, target_y1p=1._R8P)
print '(A,I0,A,3F8.4,A,F6.3)', 'u_tau given, ', report%nw, ' wall points: y1+ min mean max', report%y1p_min, &
                               report%y1p_mean, report%y1p_max, ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
print '(A,ES10.3,A,ES10.3)', '  y1_line there ', report%y1_line(m), ', y1_wall ', report%y1_wall(m)
print '(A,F7.4,A,F7.4,A,F7.4,A,I0)', '  growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                                     ', orthogonality min ', report%orth_min, ', not valid ', report%not_valid
call grid(1)%quality(orthogonality=orth)                                ! the cells behind the minimum (chapter 13)
print '(A,I0,A,I0,A,2I3)', '  cells with orthogonality < 0.5: ', count(orth < 0.5_R8P), ' of ', size(orth), &
                           ', the worst (i, j) =', minloc(orth)
call block_report_2d(grid(1), wall, nu, estimate, u_inf=u_inf, x_le=-0.25_R8P, target_y1p=1._R8P)
print '(A,L1,A,3F8.4,A,F6.3)', 'flat-plate u_tau (estimated ', estimate%utau_estimated, '): y1+ min mean max', &
                               estimate%y1p_min, estimate%y1p_mean, estimate%y1p_max, ', above target ', estimate%above
m = maxloc(estimate%y1p, dim=1)
print '(A,F7.4,A,2F9.5,A)', '  max ', estimate%y1p(m), ' at (x, y) = (', grid(1)%u(estimate%point(1, m), 1), &
                            grid(1)%v(estimate%point(1, m), 1), ')'
call save_block_report(report, 'tut15-report.txt')
open(newunit=unit, file='tut15-report.txt', action='read')
do i = 1, 7
   read(unit, '(A)') line
   if (i /= 5) print '(A)', trim(line)                                  ! line 5, the cell quality, is long
enddo
close(unit)

upper = pack([(i, i = 1, report%nw)], grid(1)%v(report%point(1, :), 1) > 0._R8P) ! the upper side, by chord
do i = 2, size(upper)
   j = i
   do while (j > 1)
      if (x_of(upper(j - 1)) <= x_of(upper(j))) exit
      upper([j - 1, j]) = upper([j, j - 1]) ; j = j - 1
   enddo
enddo
given = [(chart(x_of(upper(i)), report%y1p(upper(i))), i = 1, size(upper))]
flat = [(chart(x_of(upper(i)), estimate%y1p(upper(i))), i = 1, size(upper))]
call fig_open(un, 'tut15-y1p', 'y1+ along the upper side: u_tau given (gold), flat-plate estimate (ember)', &
              -90._R8P, 90._R8P, height=380_I4P)
call fig_segment(un, 'box', chart(-0.25_R8P, 0._R8P), chart(0.75_R8P, 0._R8P))
call fig_segment(un, 'box', chart(-0.25_R8P, 0._R8P), chart(-0.25_R8P, 2.2_R8P))
call fig_segment(un, 'box', chart(-0.25_R8P, 1._R8P), chart(0.75_R8P, 1._R8P))
call fig_label(un, 'box', chart(-0.30_R8P, 0._R8P), '0')
call fig_label(un, 'box', chart(-0.30_R8P, 1._R8P), '1')
call fig_label(un, 'box', chart(-0.30_R8P, 2._R8P), '2')
call fig_label(un, 'box', chart(-0.27_R8P, -0.18_R8P), 'leading edge')
call fig_label(un, 'box', chart(0.62_R8P, -0.18_R8P), 'trailing edge')
call fig_label(un, 'box', chart(0.18_R8P, -0.18_R8P), 'x')
call fig_segments(un, 'accent', given(1:size(given) - 1), given(2:))
call fig_segments(un, 'patch', flat(1:size(flat) - 1), flat(2:))
call fig_points(un, 'accent', given)
call fig_points(un, 'patch', flat)
call fig_close(un)
contains
   real(R8P) function x_of(w)
   !< Chordwise position of wall point w.
   integer(I4P), intent(in) :: w !< Wall point.

   x_of = grid(1)%u(report%point(1, w), 1)
   endfunction x_of

   function chart(x, y1p) result(c)
   !< The chart point of (x, y1+): x across, y1+ up (0.25 per unit).
   real(R8P), intent(in) :: x, y1p !< Chord position, y1+.
   type(vector_R8P)      :: c      !< Point.

   c = vector_R8P(x, 0.25_R8P * y1p, 0._R8P)
   endfunction chart
endprogram tut15_4_report
