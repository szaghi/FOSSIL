type(surface_stl_object) :: cube
type(swept_cells_3d_t)   :: s3
type(rigid_pose_t)       :: p1
real(R8P)                :: g(11), c(3)

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)    ! [0, 1]^3
g = [(-0.25_R8P + 0.15_R8P * i, i = 0, 10)]                              ! 10^3 cells of side 0.15
c = 0.5_R8P                                                               ! the centre stays: t = c - R c
p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=pi / 12._R8P)
p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=pi / 12._R8P, translation=c - matmul(p1%r, c))
call swept_cells_3d(cube, pose_identity(), p1, g, g, g, s3)
print '(A,L1,A,I0,A,I0)', 'cube, 15 degrees: exact ', s3%exact, ', cells swept ', count(s3%swept /= 0._R8P), &
                          ', cleared ', count(s3%cleared)
print '(A,ES9.2,A,ES9.2)', '  sum of the swept volumes ', sum(s3%swept), ', largest |residual| ', &
                           maxval(abs(s3%residual))
