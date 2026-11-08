program ref_collision
use fossil,              only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, motion_bound_3d, first_impact
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
implicit none
real(R8P), parameter        :: PI = acos(-1._R8P), TOL = 1.e-6_R8P
character(14), parameter    :: STATUS_NAME(4) = ['contact       ', 'separate      ', 'initial       ', 'not converged ']
type(surface_stl_object)    :: a, b, centred, cubes(4)
type(rigid_pose_t)          :: a0, a1, rest, p0(4), p1(4)
type(facet_object), pointer :: f
real(R8P), allocatable      :: x(:), y(:), z(:)
real(R8P)                   :: toi, d
integer(I4P)                :: i, k, status, steps, pairs, pair(2), tested

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! [0, 1]^3: the frame origin at a corner
allocate(x(0), y(0), z(0))
do i = 1, a%get_facets_number()
   f => a%facet_at(i)
   x = [x, (f%vertex(k)%x, k = 1, 3)] ; y = [y, (f%vertex(k)%y, k = 1, 3)] ; z = [z, (f%vertex(k)%z, k = 1, 3)]
enddo
a0 = pose_identity()
a1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 2, translation=[1._R8P, 0._R8P, 0._R8P])
print '(A,F8.5)', 'mu, frame origin at a corner:   ', motion_bound_3d(a0, a1, x, y, z)
print '(A,F8.5)', 'mu, frame origin at the centre: ', motion_bound_3d(a0, a1, x - 0.5_R8P, y - 0.5_R8P, z - 0.5_R8P)

call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! b rests on [0, 1]^3
rest = pose_identity()
a0 = pose_identity() ; a0%t = [-2._R8P, 0._R8P, 1._R8P + 4 * TOL]          ! a slides over b's top face, 4 tol above it
a1 = a0 ; a1%t(1) = 2._R8P
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, tolerance=TOL, distance=d, iterations=steps)
print '(3A,F9.6,A,I0,A,ES10.3)', 'plain:       ', STATUS_NAME(status), ' toi ', toi, ', steps ', steps, &
                                 ', distance ', d
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, tolerance=TOL, distance=d, iterations=steps, &
                      directional=.true., pairs=pairs)
print '(3A,F9.6,A,I0,A,ES10.3,A,I0)', 'directional: ', STATUS_NAME(status), ' toi ', toi, ', steps ', steps, &
                                      ', distance ', d, ', pairs ', pairs

call centred%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call centred%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P) ; call centred%analyze ! frame origin at the centre
cubes = centred
do i = 1, 4 ! four cubes on a row, 3 apart; 1 and 2 close in on each other, 3 turns in place, 4 drifts away
   p0(i) = pose_identity() ; p0(i)%t = [3._R8P * (i - 1), 0._R8P, 0._R8P] ; p1(i) = p0(i)
enddo
p1(1)%t(1) = 1.5_R8P ; p1(2)%t(1) = 1.6_R8P
p1(3) = pose_from_axis_angle(axis=[1._R8P, 1._R8P, 0._R8P], angle=PI, translation=p0(3)%t)
p1(4)%t(1) = 12._R8P
call first_impact(cubes, p0, p1, toi, pair, status, tolerance=TOL, tested=tested)
print '(3A,2I2,A,F9.6,A,I0,A)', 'first impact: ', STATUS_NAME(status), ' pair', pair, ', toi ', toi, &
                                ', ', tested, ' of 6 pairs tested'
endprogram ref_collision
