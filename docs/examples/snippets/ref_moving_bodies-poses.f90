use fossil, only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_compose, pose_inverse, &
                   pose_interpolate, to_world, to_body, swept_cells_3d_t, swept_cells_3d
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(rigid_pose_t)       :: p0, p1, q, half
type(surface_stl_object) :: cube, moved, rebuilt
type(swept_cells_3d_t)   :: sw
type(vector_R8P)         :: points(27)
real(R8P)                :: d_posed(27), d_moved(27), d_rebuilt(27), x(7), b(3)
integer(I4P)             :: i, j, k, l

p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[1._R8P, 1._R8P, 1._R8P], angle=PI / 6._R8P, translation=[0.1_R8P, -0.05_R8P, 0.02_R8P])
q = pose_compose(pose_inverse(p1), p1)                                          ! p1^-1 o p1 = identity
print '(A,ES9.2)', 'p1^-1 o p1 - identity: ', max(maxval(abs(q%r - p0%r)), maxval(abs(q%t)))
b = to_body(p1, to_world(p1, [0.3_R8P, 0.2_R8P, 0.1_R8P]))
print '(A,ES9.2)', 'to_body(to_world(x)) - x: ', maxval(abs(b - [0.3_R8P, 0.2_R8P, 0.1_R8P]))
half = pose_interpolate(p0, p1, 0.5_R8P)                                       ! the screw motion at mid-step
print '(A,ES9.2,A,ES9.2)', 'R(1/2)^2 - R1: ', maxval(abs(matmul(half%r, half%r) - p1%r)), &  ! angle linear in s
                           ', t(1/2) - t1/2: ', maxval(abs(half%t - 0.5_R8P * p1%t))           ! translation linear in s
q = pose_interpolate(p0, p1, 1._R8P)
print '(A,L1)', 'interpolate(p0, p1, 1) = p1 bit for bit: ', all(q%r == p1%r) .and. all(q%t == p1%t)
