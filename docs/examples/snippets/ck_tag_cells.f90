program ck_tag_cells
use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
real(R8P)                :: lo(3, 1000), hi(3, 1000), h
integer(I4P)             :: codes(1000), i, j, k, n

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
h = 0.17_R8P ! 10^3 cells on [-0.35, 1.35]^3
n = 0
do k = 0, 9 ; do j = 0, 9 ; do i = 0, 9
   n = n + 1
   lo(:, n) = -0.35_R8P + h * [i, j, k] ; hi(:, n) = lo(:, n) + h
enddo ; enddo ; enddo
call cube%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.1_R8P)
print '(A,I4)', 'cut     ', count(codes == BOX_CUT)
print '(A,I4)', 'near    ', count(codes == BOX_NEAR)
print '(A,I4)', 'inside  ', count(codes == BOX_INSIDE)
print '(A,I4)', 'outside ', count(codes == BOX_OUTSIDE)
endprogram ck_tag_cells
