call section%load_from_file(file_name='dragon-section.dat')         ! chapter 7: the section z = 0
call section%enter_device                                           ! the segment BVH and its caches
h = 1.2_R8P / nc
b = 0
do j = 1, nc ; do i = 1, nc
   b = b + 1
   lo2(:, b) = -0.6_R8P + [i - 1, j - 1] * h
   hi2(:, b) = lo2(:, b) + h
enddo ; enddo
call section%box_classify_many_device(lo=lo2, hi=hi2, codes=codes2(:, 1), radius=0.02_R8P)
call section%box_classify_many(lo=lo2, hi=hi2, codes=codes2(:, 2), radius=0.02_R8P)
print '(A,I0,A,I0)', 'curve boxes: device vs host mismatches ', count(codes2(:, 1) /= codes2(:, 2)), ' of ', nc**2
b = 0
do j = 1, np ; do i = 1, np
   b = b + 1
   pt(:, b) = -0.6_R8P + 1.2_R8P * [i - 1, j - 1] / (np - 1._R8P)
enddo ; enddo
call section%distance_many_device(points=pt, distances=d(:, 1))     ! signed, negative inside
call section%distance_many(points=pt, distances=d(:, 2))
print '(A,I0,A,I0,A,I0)', 'curve distances: device vs host mismatches ', count(d(:, 1) /= d(:, 2)), ' of ', np**2, &
      ', inside ', count(d(:, 1) < 0._R8P)
call section%exit_device
