program ck_swept_cells_3d
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, swept_cells_3d, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: box
type(swept_cells_3d_t)   :: sw
real(R8P)                :: x(13)
integer(I4P)             :: i, status

call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call box%translate(delta=vector_R8P(-0.5_R8P, -0.5_R8P, -0.5_R8P))
call box%resize(factor=vector_R8P(0.3_R8P, 0.1_R8P, 0.2_R8P)) ; call box%analyze   ! 0.3 x 0.1 x 0.2 about its origin
x = [(0.083_R8P * i, i = 0, 12)]                              ! 12^3 cells on [0, 0.996]^3
call swept_cells_3d(box, pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [0.4_R8P, 0.45_R8P, 0.47_R8P]), &
                    pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 9, [0.55_R8P, 0.5_R8P, 0.52_R8P]), &
                    x, x, x, sw, status=status)               ! turned 20 degrees about z and moved, in one step
if (status /= STATUS_OK) error stop 'swept_cells_3d failed'
print '(A,L1,A,ES9.2)', 'exact ', sw%exact, ', sum of the swept volumes ', sum(sw%swept)
print '(A,F9.6,A,F9.6)', 'fluid gained ', sum(sw%swept, mask=sw%swept > 0._R8P), ', lost ', &
                         -sum(sw%swept, mask=sw%swept < 0._R8P)
print '(A,I0,A,I0)', 'cells cleared ', count(sw%cleared), ', covered ', count(sw%covered)
print '(A,ES9.2)', 'max |volume balance residual| ', maxval(abs(sw%residual))
endprogram ck_swept_cells_3d
