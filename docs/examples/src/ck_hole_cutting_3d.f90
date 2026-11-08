!run ck_hole_cutting_3d ck_hole_cutting_3d
program ck_hole_cutting_3d
use fossil, only : surface_stl_object, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, block_logical_t, &
                   wall_cutter, hole_flags, hole_map, xray_flags
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter               :: PI = acos(-1._R8P)
type(curvilinear_block_t)          :: blocks(2)
type(topo_record_t)                :: rec
type(surface_stl_object)           :: cutter
type(block_logical_t), allocatable :: exact(:), map(:), xray(:)
real(R8P)                          :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph, g(17)
real(R8P)                          :: bx(17, 17, 17), by(17, 17, 17), bz(17, 17, 17)
integer(I4P)                       :: i, j, k, unresolved

do k = 1, 3 ; do j = 1, 13 ; do i = 1, 25                   ! block 1: a shell 1 <= r <= 1.5 around a sphere
   r = 1._R8P + 0.25_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call blocks(1)%initialize(x, y, z)
rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [25, 13, 1] ; rec%name = 'wall'
call blocks(1)%add_record(rec)                              ! the wall: face 5 (k = 1), r = 1
g = [(-1.53_R8P + 0.19125_R8P * (i - 1), i = 1, 17)]        ! block 2: a Cartesian background, 17^3 vertices
do k = 1, 17 ; do j = 1, 17
   bx(:, j, k) = g ; by(:, j, k) = g(j) ; bz(:, j, k) = g(k)
enddo ; enddo
call blocks(2)%initialize(bx, by, bz)
call wall_cutter(blocks(1:1), cutter)                        ! the wall faces as a closed surface
call hole_flags(blocks, cutter, exact, cuts=[.false., .true.])            ! signed distance < 0, exact
call hole_map(blocks, cutter, map, [12_I4P, 12_I4P, 12_I4P], cuts=[.false., .true.], unresolved=unresolved)
call xray_flags(blocks, cutter, xray, [13_I4P, 13_I4P], cuts=[.false., .true.])  ! 13 x 13 rays along z
print '(A,I0,A,L1)', 'cutter facets ', cutter%get_facets_number(), ', watertight ', cutter%is_watertight()
print '(A,I0)', 'background holes by signed distance ', count(exact(2)%v)
print '(A,I0,A,I0,A,I0)', 'by hole map ', count(map(2)%v), ', differ ', count(map(2)%v .neqv. exact(2)%v), &
                          ', unresolved ', unresolved
print '(A,I0,A,I0)', 'by X-rays   ', count(xray(2)%v), ', differ ', count(xray(2)%v .neqv. exact(2)%v)
endprogram ck_hole_cutting_3d
