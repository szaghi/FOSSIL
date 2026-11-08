use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter  :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, moved, box
type(rigid_pose_t)       :: p1
real(R8P)                :: gap, gap1, ps(3), po(3)
integer(I4P)             :: fa, fb
logical                  :: hit

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! a unit cube beside the dragon
call box%translate(x=0.7_R8P, y=-0.5_R8P, z=-0.5_R8P)
call box%analyze
call surface%distance_to(box, distance=gap, point_self=ps, point_other=po, facet_self=fa, facet_other=fb, &
                         intersecting=hit, pose_self=p1)
print '(A,F8.5,A,I0,A,I0,A,L1)', 'gap to the cube ', gap, ' (facets ', fa, ', ', fb, '), touching ', hit
print '(A,3F8.4,A,3F8.4)', 'closest pair', ps, '  -', po
