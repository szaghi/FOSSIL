use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_box, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(curve_object)      :: section
type(rigid_pose_2d_t)   :: q0, q1
type(swept_cells_2d_t)  :: sw
real(R8P)               :: x(49)
integer(I4P)            :: i

call section%load_from_file(file_name='dragon-section.dat')            ! chapter 7: the section z = 0
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]                           ! 48 x 48 cells of side 0.025
q0 = pose_2d(angle=0._R8P)
q1 = pose_2d(angle=0._R8P, translation=[0.02_R8P, 0._R8P])            ! a pure translation, 0.8 of a cell
call swept_cells_2d(section, q0, q1, x, x, sw)
print '(A,L1,A,I0,A,I0,A,I0)', 'translation: exact ', sw%exact, ', cells swept ', count(sw%swept /= 0._R8P), &
                               ', cleared ', count(sw%cleared), ', covered ', count(sw%covered)
print '(A,ES9.2)', '  sum of the swept areas ', sum(sw%swept)
print '(A,I0,A,ES9.2,A,2I3)', '  cells with residual > 1e-14: ', count(abs(sw%residual) > 1.e-14_R8P), &
                              ', largest ', maxval(abs(sw%residual)), ' in cell', maxloc(abs(sw%residual))
