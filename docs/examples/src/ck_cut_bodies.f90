!run ck_cut_bodies ck_cut_bodies
program ck_cut_bodies
use fossil, only : surface_stl_object, cut_cells_3d_t, piece_bodies_t, rigid_pose_t, pose_from_axis_angle, &
                   merge_posed_bodies, cut_cells_3d_bodies, curve_object, cut_cells_2d_t, rigid_pose_2d_t, pose_2d, &
                   cut_cells_2d_bodies, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: cube(2), merged
type(rigid_pose_t)        :: pose(2)
type(cut_cells_3d_t)      :: cut
type(piece_bodies_t)      :: pb
type(curve_object)        :: disc
type(cut_cells_2d_t)      :: cut2
integer(I4P), allocatable :: facet_body(:), facet_sub_body(:)
real(R8P)                 :: x(17), u(32), v(32)
integer(I4P)              :: i, p, status, nboth
logical                   :: overlap

call cube(1)%load_from_file(file_name='cube-01.stl', guess_format=.true.)
cube(2) = cube(1)
pose(1) = pose_from_axis_angle(axis=[1._R8P, 1._R8P, 0._R8P], angle=0.3_R8P, translation=[-1.3_R8P, 0._R8P, 0._R8P])
pose(2) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.5_R8P, translation=[0.4_R8P, 0._R8P, 0._R8P])
call merge_posed_bodies(cube, pose, merged, facet_body, facet_sub_body, status=status, overlap=overlap)
print '(A,I0,A,L1)', 'merged facets ', merged%get_facets_number(), ', overlap ', overlap
x = [(-2._R8P + 0.25_R8P * i, i = 0, 16)]                     ! 16^3 cells over [-2, 2]^3
call cut_cells_3d_bodies(cube, pose, x, x, x, cut, status=status, piece_bodies=pb)
if (status /= STATUS_OK) error stop 'cut_cells_3d_bodies failed'
nboth = 0
do p = 1, pb%np
   if (pb%start(p + 1) - pb%start(p) == 2) nboth = nboth + 1  ! a piece bounded by both cubes
enddo
print '(A,I0,A,I0,A,I0)', 'pieces ', cut%np, ', wall terms ', pb%nw, ', pieces bounded by both cubes ', nboth
print '(A,F16.12)', 'solid volume         ', sum(1._R8P - cut%alpha) * 0.25_R8P**3
print '(A,F16.12)', 'V1 + V2              ', 2._R8P * cube(1)%get_volume()
pose(2) = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=0.5_R8P, translation=[-0.4_R8P, 0._R8P, 0._R8P])
call cut_cells_3d_bodies(cube, pose, x, x, x, cut, status=status)
print '(A,I0)', 'overlapping cubes: status ', status
do i = 1, 32
   u(i) = 0.2_R8P * cos(0.19634954084936207_R8P * (i - 1)) ; v(i) = 0.2_R8P * sin(0.19634954084936207_R8P * (i - 1))
enddo
call disc%adopt_loops(u=u, v=v, loop_sizes=[32_I4P])           ! 2-D: two discs, each in its own pose
call cut_cells_2d_bodies([disc, disc], [pose_2d(0._R8P, [-0.25_R8P, 0._R8P]), pose_2d(0.1_R8P, [0.25_R8P, 0.02_R8P])], &
                         x(5:13), x(5:13), cut2, status=status, piece_bodies=pb)
print '(A,I0,A,I0)', '2-D: status ', status, ', pieces ', cut2%np
endprogram ck_cut_bodies
