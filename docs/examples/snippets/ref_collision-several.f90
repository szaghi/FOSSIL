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
