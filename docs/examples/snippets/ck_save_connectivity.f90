program ck_save_connectivity
use fossil, only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, hole_flags_2d, &
                   block_logical_t, donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, connectivity_t, &
                   build_connectivity, save_connectivity, SELECT_EXPLICIT
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t)       :: blocks(2)
type(topo_record_t)                :: rec
type(curve_object)                 :: wall
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
type(connectivity_t)               :: conn
real(R8P)                          :: u(33, 13), v(33, 13), x(41, 41), y(41, 41)
character(200)                     :: line
integer(I4P)                       :: i, j, un, status

do j = 1, 13 ; do i = 1, 33                                    ! the cylinder O-grid of "Assemble an overset grid"
   u(i, j) = (1._R8P + (j - 1) / 12._R8P) * cos(2 * PI * mod(i - 1, 32) / 32)
   v(i, j) = (1._R8P + (j - 1) / 12._R8P) * sin(2 * PI * mod(i - 1, 32) / 32)
enddo ; enddo
call blocks(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)
rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)
do j = 1, 41 ; do i = 1, 41                                    ! and its Cartesian background
   x(i, j) = -3._R8P + 0.15_R8P * (i - 1) ; y(i, j) = -3._R8P + 0.15_R8P * (j - 1)
enddo ; enddo
call blocks(2)%initialize(x, y)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter_2d(blocks(1:1), wall)
call hole_flags_2d(blocks, wall, holes, cuts=[.false., .true.])
call donor_candidates_2d(blocks, donors)
call classify_2d(blocks, donors, masks, SELECT_EXPLICIT, 3_I4P, holes=holes)
call build_connectivity(masks, donors, 3_I4P, 2_I4P, conn)
call save_connectivity(conn, 'ck_overset.txt', status)        ! plain text, one receptor per line
open(newunit=un, file='ck_overset.txt', action='read')
do i = 1, 5
   read(un, '(A)') line ; print '(A)', trim(line)
enddo
close(un)
endprogram ck_save_connectivity
