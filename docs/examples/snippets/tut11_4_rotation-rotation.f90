use fossil,     only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, &
                       pose_from_axis_angle, pose_2d, slice_surface, swept_cells_2d_t, swept_cells_3d_t, &
                       swept_cells_2d, swept_cells_3d
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(swept_cells_2d_t)   :: sw
real(R8P)                :: x(49)
integer(I4P)             :: i, n, status
integer(I4P), allocatable :: pairs(:,:)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 3_I4P, -0.2_R8P, section, merge_collinear=.true., status=status)   ! z = -0.2
call section%self_intersections(pairs, n)
print '(A,I0,A,I0)', 'section z = -0.2: ', section%get_vertices_number(), ' segments, crossing pairs ', n
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
call swept_cells_2d(section, pose_2d(angle=0._R8P), pose_2d(angle=pi / 18._R8P), x, x, sw)   ! 10 degrees
print '(A,L1,A,I0,A,I0)', 'rotation: exact ', sw%exact, ', cells swept ', count(sw%swept /= 0._R8P), &
                          ', largest evaluations per face ', sw%evaluations
print '(A,ES9.2,A,ES9.2)', '  sum of the swept areas ', sum(sw%swept), ', largest |residual| ', maxval(abs(sw%residual))
