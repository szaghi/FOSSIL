!run ck_hole_map ck_hole_map
program ck_hole_map
use fossil, only : surface_stl_object, curve_object, slice_surface, curvilinear_block_2d_t, block_logical_t, &
                   hole_flags_2d, hole_map_2d, xray_flags_2d
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)           :: wing
type(curve_object)                 :: cutter
type(curvilinear_block_2d_t)       :: background(1)
type(block_logical_t), allocatable :: exact(:), map(:), xray(:)
real(R8P)                          :: x(97, 25), y(97, 25)
integer(I4P)                       :: i, j, unresolved

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, cutter)              ! the cutter: the airfoil section
do j = 1, 25 ; do i = 1, 97                                    ! a Cartesian background, h = 0.0125
   x(i, j) = -0.35_R8P + 0.0125_R8P * (i - 1) ; y(i, j) = -0.15_R8P + 0.0125_R8P * (j - 1)
enddo ; enddo
call background(1)%initialize(x, y)
call hole_flags_2d(background, cutter, exact)                 ! signed distance < 0, every vertex
call hole_map_2d(background, cutter, map, [16_I4P, 16_I4P], unresolved=unresolved)   ! 16 x 16 map
call xray_flags_2d(background, cutter, xray, 33_I4P)           ! 33 X-rays along y
print '(A,I0)', 'holes by signed distance ', count(exact(1)%v)
print '(A,I0,A,I0,A,I0)', 'holes by hole map        ', count(map(1)%v), ', differ ', count(map(1)%v .neqv. exact(1)%v), &
                          ', unresolved ', unresolved
print '(A,I0,A,I0)', 'holes by X-rays          ', count(xray(1)%v), ', differ ', count(xray(1)%v .neqv. exact(1)%v)
endprogram ck_hole_map
