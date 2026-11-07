program ck_ccd
use fossil, only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, CCD_CONTACT
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: a, b, plate
type(rigid_pose_t)       :: a0, a1, rest
real(R8P)                :: toi, d, pa(3), pb(3), d0, d1
integer(I4P)             :: status, steps

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call b%translate(x=3._R8P) ; call b%analyze                     ! b rests at x = 3..4
rest = pose_identity()
! over the step, a turns 90 degrees about z (about its own origin) and slides 4 along x
a0 = pose_identity()
a1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI / 2, translation=[4._R8P, 0._R8P, 0._R8P])
call a%time_of_impact(b, a0, a1, rest, rest, toi, status, distance=d, point_self=pa, point_other=pb, &
                      iterations=steps)
print '(A,L1,A,F9.6,A,I0,A)', 'contact ', status == CCD_CONTACT, ' at s = ', toi, ' (', steps, ' steps)'
print '(A,ES10.3)', 'distance there ', d
print '(A,3F8.4)', 'on a ', pa
print '(A,3F8.4)', 'on b ', pb
! a plate 1e-3 thick crossed within one step: apart at both ends, the contact is still found
call plate%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call plate%resize(x=4._R8P, y=4._R8P, z=1.e-3_R8P) ; call plate%translate(x=-2._R8P, y=-2._R8P) ; call plate%analyze
a0%t = [0._R8P, 0._R8P, 2._R8P] ; a1 = a0 ; a1%t = [0._R8P, 0._R8P, -3._R8P]
call a%distance_to(plate, d0, pose_self=a0) ; call a%distance_to(plate, d1, pose_self=a1)
call a%time_of_impact(plate, a0, a1, rest, rest, toi, status)
print '(A,2F7.3,A,L1,A,F9.6)', 'plate: distance at s = 0, 1 ', d0, d1, '; contact ', status == CCD_CONTACT, ' at s = ', toi
endprogram ck_ccd
