!run figck-exposed-points fig_ck_exposed_points
!fig ck-exposed-points
!< Figure of the cookbook recipe "Find the points a moving body uncovers": the background of ck_exposed_points around
!< the cylinder's trailing side, the wall at step n (faint) and at step n + 1, moved by 0.6 h (stone); background holes
!< at n + 1 in teal and the exposed points (holes or receptors at n, field points at n + 1) in ember.
program fig_ck_exposed_points
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, hole_flags_2d, &
                   block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, exposed_points, &
                   rigid_pose_2d_t, pose_2d, SELECT_EXPLICIT
use fossil_fig, only : fig_open, fig_zoom, fig_points, fig_segments, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P), H = 0.15_R8P
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: wall
type(overset_mask_t),  allocatable :: m0(:), m1(:), m2(:)
type(block_logical_t), allocatable :: exposed(:)
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                       :: i, j, violations, n, un
type(vector_R8P), allocatable      :: pts(:)
type(vector_R8P)                   :: a0(32), b0(32), a1(32), b1(32)
character(90)                      :: caption

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
allocate(pts(41 * 41)) ; n = 0
do j = 1, 41 ; do i = 1, 41
   n = n + 1 ; pts(n) = vector_R8P(x(i, j), y(i, j), 0._R8P)
enddo ; enddo
a0 = [(vector_R8P(u(i, 1), v(i, 1), 0._R8P), i = 1, 32)] ; b0 = [a0(2:32), a0(1)]
a1 = a0 + vector_R8P(0.6_R8P * H, 0._R8P, 0._R8P) ; b1 = b0 + vector_R8P(0.6_R8P * H, 0._R8P, 0._R8P)
write(caption, '(A,I0,A)') 'a step of 0.6 h: ', count(exposed(2)%v), ' background points exposed (ember), holes in teal'
call fig_open(un, 'ck-exposed-points', trim(caption), -90._R8P, 90._R8P, height=420_I4P)
call fig_zoom(un, vector_R8P(-0.2_R8P, 0._R8P, 0._R8P), 1.45_R8P)
call fig_segments(un, 'ghost', a0, b0)
call fig_segments(un, 'body', a1, b1)
call fig_points(un, 'ghost', pack(pts, m1(2)%iota /= 0 .and. .not. exposed(2)%v))
call fig_points(un, 'inside', pack(pts, m1(2)%iota == 0))
call fig_points(un, 'hole', pack(pts, exposed(2)%v))
call fig_close(un)
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
endprogram fig_ck_exposed_points
