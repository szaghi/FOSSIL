n = 0
do k = 0, 2 * NF                                             ! the same period in steps half as long
   poses(1) = pose_2d(angle=PI / 6._R8P * sin(PI * mod(k, 2 * NF) / NF), &
                      translation=[0._R8P, 0.3_R8P * sin(PI * k / NF + PI / 2._R8P)])
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   if (k > 0) call exposed_points(old, masks, exposed, violations)
   if (k > 0) n = n + violations
   old = masks
enddo
print '(A,I0,A,I0)', 'steps per period ', 2 * NF, ': holes to field ', n
