use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY, curve_object, &
                       wall_cutter_2d, hole_flags_2d, block_logical_t, donor_set_t, donor_candidates_2d, &
                       overset_mask_t, classify_2d, SELECT_EXPLICIT, ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ)
integer(I4P)                       :: i, j, k

call make_ogrid(u, v)                                        ! the O-grid of step 1
do j = 1, MJ                                                 ! the background of step 3
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
call blocks(2)%build_search
do k = 10, 16, 2                                             ! the first k - 1 cell layers of the O-grid only
   call blocks(1)%destroy
   call blocks(1)%initialize(u(:, 1:k), v(:, 1:k), name='thin')
   call detect_topology_2d(blocks(1:1))
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
   call blocks(1)%add_record(rec)
   rec%face = 4 ; rec%lo = [1, k, 1] ; rec%hi = [NI, k, 1] ; rec%name = 'overlap'
   call blocks(1)%add_record(rec)
   call blocks(1)%build_search
   call wall_cutter_2d(blocks(1:1), cutter)
   call hole_flags_2d(blocks, cutter, holes, cuts=[.false., .true.])
   call donor_candidates_2d(blocks, donors)
   call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
   print '(A,I2,A,I2,A,I2,A,I2,A,I2,A)', 'O-grid of ', k - 1, ' layers: receptors ', count(masks(1)%iota < 0), &
         ', orphans ', count(masks(1)%orphan > 0), ' (stencil on a hole ', &
         count(masks(1)%orphan == ORPHAN_STENCIL_HOLE), ', receptor in the stencil ', &
         count(masks(1)%orphan == ORPHAN_RECEPTOR_DONOR), ')'
enddo
