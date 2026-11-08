program ck_overset_3d
use fossil, only : surface_stl_object, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, block_logical_t, &
                   wall_cutter, hole_flags, donor_set_t, donor_candidates, overset_mask_t, classify, fringe_levels, &
                   SELECT_EXPLICIT, FRINGE_INNER_1, FRINGE_OUTER_1
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_t)          :: blocks(2)
type(topo_record_t)                :: rec
type(surface_stl_object)           :: cutter
type(block_logical_t), allocatable :: holes(:)
type(donor_set_t),     allocatable :: donors(:)
type(overset_mask_t),  allocatable :: masks(:)
real(R8P)                          :: x(25, 13, 6), y(25, 13, 6), z(25, 13, 6), r, th, ph, g(17)
real(R8P)                          :: bx(17, 17, 17), by(17, 17, 17), bz(17, 17, 17)
integer(I4P)                       :: i, j, k, b

do k = 1, 6 ; do j = 1, 13 ; do i = 1, 25                   ! block 1: a shell 1 <= r <= 2 around a sphere
   r = 1._R8P + 0.2_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call blocks(1)%initialize(x, y, z)
rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [25, 13, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)                              ! k = 1: the wall
rec%face = 6 ; rec%lo = [1, 1, 6] ; rec%hi = [25, 13, 6] ; rec%name = 'overlap'
call blocks(1)%add_record(rec)                              ! k = 6: receives from the background
g = [(-2.4_R8P + 0.3_R8P * (i - 1), i = 1, 17)]             ! block 2: a Cartesian background, h = 0.3
do k = 1, 17 ; do j = 1, 17
   bx(:, j, k) = g ; by(:, j, k) = g(j) ; bz(:, j, k) = g(k)
enddo ; enddo
call blocks(2)%initialize(bx, by, bz)
call blocks(1)%build_search ; call blocks(2)%build_search
call wall_cutter(blocks(1:1), cutter)
call hole_flags(blocks, cutter, holes, cuts=[.false., .true.])
call donor_candidates(blocks, donors)                       ! every vertex located in the other block
call classify(blocks, donors, masks, SELECT_EXPLICIT, 2_I4P, holes=holes)   ! stencils of width 2
call fringe_levels(blocks, masks)
do b = 1, 2
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'block ', b, ': field ', count(masks(b)%iota == 1), ', holes ', &
         count(masks(b)%iota == 0), ', receptors ', count(masks(b)%iota < 0), ', orphans ', count(masks(b)%orphan > 0)
enddo
print '(A,I0,A,I0)', 'background INNER 1 fringe ', count(masks(2)%fringe == FRINGE_INNER_1), &
                     ', shell OUTER 1 fringe ', count(masks(1)%fringe == FRINGE_OUTER_1)
endprogram ck_overset_3d
