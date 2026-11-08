!run ck_first_impact ck_first_impact
program ck_first_impact
use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, &
                   pose_2d, first_impact, first_impact_2d, CCD_CONTACT
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: cubes(4)
type(curve_object)       :: squares(4)
type(rigid_pose_t)       :: p0(4), p1(4)
type(rigid_pose_2d_t)    :: q0(4), q1(4)
real(R8P)                :: toi
integer(I4P)             :: i, pair(2), status, tested

do i = 1, 4                                                   ! four unit cubes on a row, 1 apart: x = 0, 2, 4, 6
   call cubes(i)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   p0(i) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [2._R8P * (i - 1), 0._R8P, 0._R8P]) ; p1(i) = p0(i)
enddo
p1(2)%t = [2._R8P, 2._R8P, 0._R8P]                            ! cube 2 slides up out of the row
p1(4) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 4, [4.5_R8P, 0._R8P, 0._R8P])  ! cube 4 turns towards 3
call first_impact(cubes, p0, p1, toi, pair, status, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', '3-D: contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
do i = 1, 4                                                   ! the same in 2-D: unit squares, the same motions
   call squares(i)%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
   q0(i) = pose_2d(0._R8P, [2._R8P * (i - 1), 0._R8P]) ; q1(i) = q0(i)
enddo
q1(2) = pose_2d(0._R8P, [2._R8P, 2._R8P]) ; q1(4) = pose_2d(PI / 4, [4.5_R8P, 0._R8P])
call first_impact_2d(squares, q0, q1, toi, pair, status, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', '2-D: contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
endprogram ck_first_impact
