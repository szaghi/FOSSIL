do j = 1, 33                                                 ! a Cartesian background, spacing 0.15
   do i = 1, 41
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(2)%build_search
block
type(curve_object)                 :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)               :: conn
integer(I4P)                       :: iterations, circular
character(200)                     :: line

call wall_cutter_2d(blocks(1:1), wall)                       ! the cutter: the O-grid's wall edges
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.]) ! cut the background only
print '(A,I0,A,I0)', 'background holes ', count(holes(2)%v), ' of ', size(holes(2)%v)
call donor_candidates_2d(blocks, donors)                     ! every vertex located in the other block
print '(A,I0,A,I0)', 'donor candidates: O-grid ', size(donors(1)%block), ', background ', size(donors(2)%block)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes, iterations=iterations, &
                 circular=circular)
call fringe_levels_2d(blocks, masks)
do r = 1, 2
   print '(A,A,A,I0,A,I0,A,I0,A,I0)', '  ', blocks(r)%name, ': field ', count(masks(r)%iota == 1), ', holes ', &
         count(masks(r)%iota == 0), ', receptors ', count(masks(r)%iota < 0), ', orphans ', count(masks(r)%orphan > 0)
enddo
print '(A,I0,A,I0)', 'fixed point in ', iterations, ' iterations, circular pairs ', circular
print '(A,4(1X,I0))', 'fringe INNER 1, INNER 2, OUTER 1, OUTER 2:', &
      sum([(count(masks(r)%fringe == FRINGE_INNER_1), r = 1, 2)]), sum([(count(masks(r)%fringe == FRINGE_INNER_2), r = 1, 2)]), &
      sum([(count(masks(r)%fringe == FRINGE_OUTER_1), r = 1, 2)]), sum([(count(masks(r)%fringe == FRINGE_OUTER_2), r = 1, 2)])
call overlap_width_2d(blocks, masks, d)
print '(A,F7.4)', 'overlap width ', d(1, 2)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
call save_connectivity(conn, 'tut13-connectivity.txt')
open(newunit=un, file='tut13-connectivity.txt', action='read')
do r = 1, 5                                                  ! the header and the first receptor
   read(un, '(A)') line
   print '(A)', trim(line)
enddo
close(un)
