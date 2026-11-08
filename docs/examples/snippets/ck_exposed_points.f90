program ck_exposed_points
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, hole_flags_2d, &
                   block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, exposed_points, &
                   rigid_pose_2d_t, pose_2d, SELECT_EXPLICIT
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P), H = 0.15_R8P
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: wall
type(overset_mask_t),  allocatable :: m0(:), m1(:), m2(:)
type(block_logical_t), allocatable :: exposed(:)
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                       :: i, j, violations

do j = 1, 13 ; do i = 1, 33                                   ! an O-grid around the unit cylinder (body frame)
   r = 1._R8P + (j - 1) / 12._R8P
   u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
enddo ; enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, 41 ; do i = 1, 41                                   ! a Cartesian background, h = 0.15
   x(i, j) = -3._R8P + H * (i - 1) ; y(i, j) = -3._R8P + H * (j - 1)
enddo ; enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter_2d(blocks(1:1), wall)
call assemble(0._R8P, m0)                                     ! step n: the cylinder at the origin
call assemble(0.6_R8P * H, m1)                                ! step n + 1: moved by 0.6 h along x
call exposed_points(m0, m1, exposed, violations)
print '(A,I0,A,I0)', 'step of 0.6 h: background points exposed ', count(exposed(2)%v), &
                     ', holes become field points ', violations
call assemble(3._R8P * H, m2)                                 ! a step of 3 h instead
call exposed_points(m0, m2, exposed, violations)
print '(A,I0,A,I0)', 'step of 3 h:   background points exposed ', count(exposed(2)%v), &
                     ', holes become field points ', violations
contains
   subroutine assemble(dx, masks)
   !< Holes, donors and classification with the O-grid posed at (dx, 0).
   real(R8P),                         intent(in)  :: dx
   type(overset_mask_t), allocatable, intent(out) :: masks(:)
   type(block_logical_t), allocatable             :: holes(:)
   type(donor_set_t),     allocatable             :: donors(:)
   type(rigid_pose_2d_t)                          :: poses(2)

   poses(1) = pose_2d(0._R8P, [dx, 0._R8P]) ; poses(2) = pose_2d(0._R8P)
   call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 2_I4P, holes=holes)
   endsubroutine assemble
endprogram ck_exposed_points
