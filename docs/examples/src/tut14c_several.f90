!run tut14c tut14c_several
!< Tutorial 14, run 3: the dragon and three pebbles over one step: the earliest contact (first_impact), the same toi
!< by the loop over all pairs, and the bodies near each other at that instant (proximity_pairs_of).
program tut14c_several
use fossil, only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, &
                   first_impact, proximity_pairs_of, CCD_CONTACT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)  :: dragon, pebble, bodies(4), at_toi(4)
type(rigid_pose_t)        :: poses0(4), poses1(4)
type(vector_R8P)          :: bmin, bmax
real(R8P)                 :: scale, toi, toi_all, s
real(R8P),    allocatable :: gaps(:)
integer(I4P), allocatable :: near(:, :)
integer(I4P)              :: i, j, pair(2), tested, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! the dragon and the pebble of run 1
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-0.5_R8P * (bmin + bmax))
call dragon%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call dragon%resize(x=scale, y=scale, z=scale)
call dragon%analyze
call pebble%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call pebble%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
call pebble%resize(x=0.04_R8P, y=0.04_R8P, z=0.04_R8P)
call pebble%analyze
!region several
bodies(1) = dragon ; poses0(1) = pose_identity() ; poses1(1) = poses0(1)       ! the dragon at rest
bodies(2:4) = pebble                                                             ! three pebbles, three throws
poses0(2) = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=0._R8P, translation=[-0.6_R8P, -0.3_R8P, -0.2_R8P])
poses1(2) = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=PI / 2, translation=[0.6_R8P, -0.3_R8P, -0.2_R8P])
poses0(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=0._R8P, translation=[0._R8P, 0.3_R8P, 0.6_R8P])
poses1(3) = pose_from_axis_angle(axis=[1._R8P, 0._R8P, 0._R8P], angle=PI, translation=[0._R8P, 0.3_R8P, -0.6_R8P])
poses0(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0._R8P, translation=[-0.15_R8P, 0.3_R8P, 0.45_R8P])
poses1(4) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=PI, translation=[1.05_R8P, 0.3_R8P, 0.45_R8P])
call first_impact(bodies, poses0, poses1, toi, pair, status, tolerance=1.e-6_R8P, tested=tested)
print '(A,L1,A,F9.6,A,2I2,A,I0,A)', 'first contact ', status == CCD_CONTACT, ' at s = ', toi, ', bodies', pair, &
                                    ' (', tested, ' of 6 pairs tested)'
!endregion several
!region loop
toi_all = 1._R8P                                                                 ! the same by every pair
do i = 1, 3
   do j = i + 1, 4
      call bodies(i)%time_of_impact(bodies(j), poses0(i), poses1(i), poses0(j), poses1(j), s, status, &
                                     tolerance=1.e-6_R8P)
      if (status == CCD_CONTACT) then
         toi_all = min(toi_all, s)
         print '(A,I0,A,I0,A,F9.6)', '  pair ', i, '-', j, ': contact at s = ', s
      else
         print '(A,I0,A,I0,A)', '  pair ', i, '-', j, ': separate'
      endif
   enddo
enddo
print '(A,L1)', 'equal to the loop over all pairs, bit for bit: ', toi_all == toi
!endregion loop
!region near
do i = 1, 4                                                                      ! the bodies at the first contact
   at_toi(i) = bodies(i)
   call at_toi(i)%move(pose_interpolate(poses0(i), poses1(i), toi))
enddo
call proximity_pairs_of(at_toi, cutoff=0.3_R8P, pairs=near, distances=gaps)
do i = 1, size(gaps)
   print '(A,I0,A,I0,A,ES10.3)', 'within 0.3 at s = toi: bodies ', near(1, i), ' and ', near(2, i), ', gap ', gaps(i)
enddo
!endregion near
endprogram tut14c_several
