call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%move(pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=0.3_R8P, translation=[0.21_R8P, 0.13_R8P, 0.07_R8P]))
g = [(-0.15_R8P + 0.2_R8P * p, p = 0, 8)]
call cut_cells_3d(cube, g, g, g, cut, status=status, links=links, edges=.true., moments_order=2_I4P, quadrature=.true.)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'status ', status, ', pieces ', cut%np, ', multi-piece cells ', &
                                    count(cut%pieces > 1), ', links ', size(links%piece), ', unmatched ', links%unmatched
print '(A,I0,A,I0)', 'cut edges (x, y, z) ', count(cut%edge_x > 0._R8P .and. cut%edge_x < 1._R8P) + &
                     count(cut%edge_y > 0._R8P .and. cut%edge_y < 1._R8P) + &
                     count(cut%edge_z > 0._R8P .and. cut%edge_z < 1._R8P), ', nodes inside ', count(cut%node_in)
q = 0._R8P
do p = 1, cut%np                                       ! the zeroth moment is the volume
   q = max(q, abs(cut%piece_moments(1, p) - cut%piece_volume(p)))
enddo
print '(A,I0,A,ES9.2,A,I0)', 'moments per piece ', size(cut%piece_moments, 1), ', |m0 - volume| ', q, &
                             ', pieces with quad_status 1 ', count(cut%quad_status /= 0)
