program ck_overset
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, hole_flags_2d, &
                   block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, fringe_levels_2d, &
                   connectivity_t, build_connectivity, SELECT_EXPLICIT, FRINGE_INNER_1
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: blocks(2)
type(topo_record_t)          :: rec
type(curve_object)           :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)         :: conn
real(R8P)                    :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                 :: i, j

do j = 1, 13                                                  ! an O-grid around the unit cylinder
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)                                ! the wall
rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)                                ! the boundary that receives
do j = 1, 41                                                  ! a Cartesian background
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter_2d(blocks(1:1), wall)                        ! the cutter: the O-grid's wall
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
call donor_candidates_2d(blocks, donors)                      ! every vertex located in the other block
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
call fringe_levels_2d(blocks, masks)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
print '(A,I0,A,I0)', 'background: holes ', count(masks(2)%iota == 0), ', receptors ', count(masks(2)%iota < 0)
print '(A,I0,A,I0)', 'O-grid: receptors ', count(masks(1)%iota < 0), ', orphans ', count(masks(1)%orphan > 0) + &
                     count(masks(2)%orphan > 0)
print '(A,I0,A,I0)', 'INNER 1 fringe ', count(masks(2)%fringe == FRINGE_INNER_1), ', receptor list ', conn%n
print '(A,I0,A,3I3,A,2F8.4)', 'first receptor: block ', conn%receptor_block(1), ', donor cell', conn%corner(:, 1), &
                              ', xi', conn%xi(1:2, 1)
endprogram ck_overset
