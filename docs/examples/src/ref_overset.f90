!run ref_overset ref_overset
program ref_overset
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, hole_flags_2d, &
                   block_logical_t, block_real_t, overlap_matrix_2d, donor_set_t, donor_candidates_2d, overset_mask_t, &
                   classify_2d, SELECT_EXPLICIT, SELECT_IMPLICIT, fringe_levels_2d, overlap_width_2d, connectivity_t, &
                   build_connectivity, save_connectivity, FRINGE_INNER_1, FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: wall
type(block_logical_t), allocatable :: holes(:), holes_offset(:)
type(block_real_t),    allocatable :: dist(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:), implicit(:)
type(connectivity_t)               :: conn
integer(I4P),          allocatable :: pair(:,:)
real(R8P),             allocatable :: width(:,:)
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41), r
integer(I4P)                       :: i, j, iterations, circular
character(100)                     :: line

!region blocks
do j = 1, 13 ! an O-grid around the unit cylinder, out to r = 2
   r = 1._R8P + (j - 1) / 12._R8P
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, 41 ! a Cartesian background, spacing 0.15
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
!endregion blocks

!region holes
call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.], distances=dist)
call hole_flags_2d(blocks, wall, holes_offset, offset=0.3_R8P, cuts=[.false., .true.])
print '(A,I0,A,I0,A,F7.4)', 'holes ', count(holes(2)%v), ', with offset 0.3 ', count(holes_offset(2)%v), &
                            ', nearest field point to the wall ', minval(dist(2)%v, mask=.not. holes(2)%v)
!endregion holes

!region donors
call overlap_matrix_2d(blocks, pair)
print '(A,4I2)', 'pair matrix ', pair
call donor_candidates_2d(blocks, donors, pair=pair, wall=wall)
print '(A,I0,A,I0)', 'candidates: O-grid vertices in the background ', size(donors(1)%block), &
                     ', background vertices in the O-grid ', size(donors(2)%block)
print '(A,2F8.4)', 'donor / receptor volume, O-grid vertices: min max ', minval(donors(1)%ratio), maxval(donors(1)%ratio)
!endregion donors

!region classify
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes_offset, iterations=iterations, &
                 circular=circular)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'explicit: receptors ', count(masks(1)%iota < 0), ' + ', count(masks(2)%iota < 0), &
                                    ', orphans ', count(masks(1)%orphan > 0) + count(masks(2)%orphan > 0), &
                                    ', iterations ', iterations, ', circular ', circular
call classify_2d(blocks, donors, implicit, SELECT_IMPLICIT, 3_I4P, iterations=iterations)
print '(A,I0,A,I0,A,I0,A,I0)', 'implicit: receptors ', count(implicit(1)%iota < 0), ' + ', &
                               count(implicit(2)%iota < 0), ', holes ', count(implicit(2)%iota == 0), &
                               ', iterations ', iterations
call fringe_levels_2d(blocks, masks, levels=2_I4P)
print '(A,4I5)', 'fringe INNER 1, INNER 2 (background), OUTER 1, OUTER 2 (O-grid) ', &
                 count(masks(2)%fringe == FRINGE_INNER_1), count(masks(2)%fringe == FRINGE_INNER_2), &
                 count(masks(1)%fringe == FRINGE_OUTER_1), count(masks(1)%fringe == FRINGE_OUTER_2)
call overlap_width_2d(blocks, masks, width)
print '(A,F8.4)', 'overlap width ', width(1, 2)
!endregion classify

!region output
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
call save_connectivity(conn, 'cylinder.conn')
open(newunit=i, file='cylinder.conn', action='read')
do j = 1, 5
   read(i, '(A)') line ; print '(A)', line(1:min(len_trim(line), 100))
enddo
close(i)
!endregion output
endprogram ref_overset
