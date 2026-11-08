use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY, rigid_pose_2d_t, &
                       pose_2d, to_world, curve_object, wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, &
                       donor_candidates_2d, overset_mask_t, classify_2d, exposed_points, SELECT_EXPLICIT
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33, NF = 12   ! NF steps per period
real(R8P),    parameter            :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(rigid_pose_2d_t)              :: poses(2)
type(block_logical_t), allocatable :: holes(:), exposed(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:), old(:)
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ), angle
integer(I4P)                       :: i, j, k, violations, n

call make_ogrid(u, v)                                        ! the O-grid of step 1, in its body frame
call blocks(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(blocks(1:1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, MJ
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(1)%build_search ; call blocks(2)%build_search    ! once: the O-grid moves by its pose only
call wall_cutter_2d(blocks(1:1), cutter)                     ! once, in the body frame
print '(A)', 'step  pitch  holes  receptors (bg + O-grid)  orphans  exposed  holes to field'
do k = 0, NF
   angle = PI / 6._R8P * sin(2 * PI * mod(k, NF) / NF)     ! pitch 30 degrees, heave 0.3
   poses(1) = pose_2d(angle=angle, translation=[0._R8P, 0.3_R8P * sin(2 * PI * k / NF + PI / 2._R8P)])
   poses(2) = pose_2d(angle=0._R8P)                          ! the background stays
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.], poses=poses, cutter_pose=poses(1))
   call donor_candidates_2d(blocks, donors, poses=poses)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   if (k == 0) old = masks
   call exposed_points(old, masks, exposed, violations)      ! uncovered since the previous step
   print '(I4,F7.1,I7,I8,A,I3,I15,I9,I12)', k, angle * 180._R8P / PI, count(masks(2)%iota == 0), &
         count(masks(2)%iota < 0), ' +', count(masks(1)%iota < 0), count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0), &
         count(exposed(2)%v), violations
   old = masks
enddo
