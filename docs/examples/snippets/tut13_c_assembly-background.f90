use fossil,     only : curvilinear_block_2d_t, topo_record_t, detect_topology_2d, TOPO_BOUNDARY
use fossil,     only : curve_object, wall_cutter_2d, hole_flags_2d, hole_map_2d, block_logical_t, donor_set_t, &
                       donor_candidates_2d, overset_mask_t, classify_2d, fringe_levels_2d, overlap_width_2d, &
                       connectivity_t, build_connectivity, save_connectivity, SELECT_EXPLICIT, FRINGE_INNER_1, &
                       FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
use penf,       only : I4P, R8P

implicit none
integer(I4P), parameter            :: NI = 49, NJ = 17, MI = 41, MJ = 33
type(curvilinear_block_2d_t)       :: blocks(2)              ! 1: the O-grid, 2: the background
type(topo_record_t)                :: rec
type(curve_object)                 :: cutter
type(block_logical_t), allocatable :: holes(:), holes_map(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)               :: conn
real(R8P)                          :: u(NI, NJ), v(NI, NJ), x(MI, MJ), y(MI, MJ)
real(R8P),             allocatable :: d(:,:)
integer(I4P)                       :: i, j, b, iterations, circular, un
character(200)                     :: line

call make_ogrid(u, v)                                        ! the O-grid of step 1, with its records
call blocks(1)%initialize(u, v, name='ellipse')
call detect_topology_2d(blocks(1:1))
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [NI, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, NJ, 1] ; rec%hi = [NI, NJ, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, MJ                                                 ! a Cartesian background, spacing 0.15
   do i = 1, MI
      x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -2.4_R8P + 0.15_R8P * (j - 1)
   enddo
enddo
call blocks(2)%initialize(x, y, name='background')
do b = 1, 2
   call blocks(b)%build_search
enddo
