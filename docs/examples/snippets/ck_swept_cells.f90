program ck_swept_cells
use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(curve_object)     :: plate
type(swept_cells_2d_t) :: sw
real(R8P)              :: x(21)
integer(I4P)           :: i, status

call plate%adopt_loops(u=[-0.15_R8P, 0.15_R8P, 0.15_R8P, -0.15_R8P], v=[-0.05_R8P, -0.05_R8P, 0.05_R8P, 0.05_R8P], &
                       loop_sizes=[4_I4P])                               ! a 0.3 x 0.1 plate about its body origin
x = [(0.05_R8P * i, i = 0, 20)]                                          ! 20 x 20 cells on [0, 1]
call swept_cells_2d(plate, pose_2d(0._R8P, [0.4_R8P, 0.45_R8P]), pose_2d(acos(-1._R8P) / 6, [0.55_R8P, 0.5_R8P]), &
                    x, x, sw, status=status)                             ! turned 30 degrees and moved in one step
if (status /= STATUS_OK) error stop 'swept_cells_2d failed'
print '(A,L1,A,ES9.2)', 'exact ', sw%exact, ', sum of the swept areas ', sum(sw%swept)
print '(A,F8.5,A,F8.5)', 'fluid gained ', sum(sw%swept, mask=sw%swept > 0._R8P), ', lost ', &
                         -sum(sw%swept, mask=sw%swept < 0._R8P)
print '(A,I0,A,I0)', 'cells cleared ', count(sw%cleared), ', covered ', count(sw%covered)
print '(A,ES9.2)', 'max |volume balance residual| ', maxval(abs(sw%residual))
endprogram ck_swept_cells
