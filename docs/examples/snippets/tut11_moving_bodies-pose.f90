use fossil, only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, &
                   pose_2d, swept_cells_2d_t, swept_cells_2d, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                        ! cell centres per direction
real(R8P),    parameter       :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object)      :: surface, moved
type(curve_object)            :: section
type(rigid_pose_t)            :: p0, p1
type(rigid_pose_2d_t)         :: q0, q1
type(swept_cells_2d_t)        :: sw
type(vector_R8P), allocatable :: centres(:)
real(R8P),        allocatable :: d_posed(:), d_moved(:)
logical,          allocatable :: fresh(:), dead(:)
real(R8P)                     :: h, x(49)
integer(I4P)                  :: i, j, k, c, checked

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / n
allocate(centres(n**3), d_posed(n**3), d_moved(n**3), fresh(n**3), dead(n**3))
c = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   c = c + 1
   centres(c) = (-0.6_R8P + (i - 0.5_R8P) * h) * ex_R8P + (-0.6_R8P + (j - 0.5_R8P) * h) * ey_R8P + &
                (-0.6_R8P + (k - 0.5_R8P) * h) * ez_R8P
enddo ; enddo ; enddo
p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 1._R8P, 0._R8P], angle=pi / 18._R8P, translation=[0.02_R8P, 0._R8P, 0._R8P])
