do i = 1, 64
   cu(i) = 0.5_R8P * cos(2 * PI * (i - 1) / 64) ; cv(i) = 0.5_R8P * sin(2 * PI * (i - 1) / 64)
enddo
call circle%adopt_loops(u=cu, v=cv, loop_sizes=[64_I4P])                       ! a 64-gon of radius 0.5
call square%adopt_loops(u=[-0.5_R8P, 0.5_R8P, 0.5_R8P, -0.5_R8P], v=[-0.5_R8P, -0.5_R8P, 0.5_R8P, 0.5_R8P], &
                        loop_sizes=[4_I4P])                                    ! a unit square about the origin
call circle%distance_to(square, distance=d, point_self=qa, point_other=qb, &
                        pose_other=pose_2d(angle=PI / 4, translation=[1.5_R8P, 0._R8P]))   ! a corner towards the circle
print '(A,F8.5,A,F8.5,A,2F8.5,A,2F8.5)', '2-D: ', d, ' (exact ', 1._R8P - sqrt(0.5_R8P), ') from', qa, ' to', qb
