call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
h = 0.25_R8P
g = [(-0.1_R8P + h * i, i = 0, 6)]                                     ! 6^3 cells, the block covers the body
call cut_cells_3d(surface=cube, x=g, y=g, z=g, cut=cut3)
print '(A,I0,A,I0,A,I0)', 'cube: cut cells ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1._R8P), &
                          ', multi-piece ', count(cut3%pieces > 1), ', flagged ', count(cut3%pieces < 0)
print '(A,ES9.2)', '  covering sum: solid volume - body volume ', sum(1._R8P - cut3%alpha) * h**3 - cube%get_volume()
res = 0._R8P
do k = 1, 6 ; do j = 1, 6 ; do i = 1, 6                                ! divergence theorem, cell by cell
   res = max(res, maxval(abs([cut3%face_x(i + 1, j, k) - cut3%face_x(i, j, k), &
                              cut3%face_y(i, j + 1, k) - cut3%face_y(i, j, k), &
                              cut3%face_z(i, j, k + 1) - cut3%face_z(i, j, k)] * h**2 - cut3%wet_normal(:, i, j, k))))
enddo ; enddo ; enddo
print '(A,ES9.2)', '  largest divergence residual ', res
print '(A,I0,A,ES10.3)', '  small cells (0 < alpha < 1e-2) ', count(cut3%alpha > 0._R8P .and. cut3%alpha < 1.e-2_R8P), &
                         ', smallest alpha ', minval(cut3%alpha, cut3%alpha > 0._R8P)
