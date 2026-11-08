use fossil,     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, inertia_to_world, STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_label, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter  :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, wrap, rest
type(rigid_pose_t)       :: p1
real(R8P)                :: mass, volume, centre(3), inertia(3,3), inertia1(3,3), moments(3), axes(3,3)
integer(I4P)             :: status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%mass_properties(mass=mass, status=status)
print '(A,L1)', 'the placed dragon is a closed body: ', status == STATUS_OK     ! chapter 2: it is not
call surface%alpha_wrap(alpha=0.06_R8P, offset=0.006_R8P, wrapped=wrap, status=status)  ! a closed surrogate
print '(A,I0,A,L1)', 'alpha wrap: ', wrap%get_facets_number(), ' facets, watertight ', wrap%is_watertight()
call wrap%mass_properties(mass=mass, volume=volume, center=centre, inertia=inertia, principal_moments=moments, &
                          principal_axes=axes, density=1000._R8P, status=status)     ! a body of density 1000
if (status /= STATUS_OK) error stop 'the wrap is not closed'
print '(A,ES11.4,A,ES11.4)', 'volume ', volume, ', mass ', mass
print '(A,3F8.4)', 'centre of mass   ', centre
print '(A,3ES11.4)', 'principal moments', moments
