use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nb = 8                      ! blocks per direction
integer(I4P), parameter  :: nc = 32                     ! cells per direction of the tagged block
type(surface_stl_object) :: surface
real(R8P)                :: lo(3, nb**3), hi(3, nb**3), x(nc + 1), hb
integer(I4P)             :: codes(nb**3), tags(nc, nc, nc), i, j, k, b

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
hb = 1.2_R8P / nb                                       ! the domain [-0.6, 0.6]^3 cut into 8^3 blocks
b = 0
do k = 1, nb ; do j = 1, nb ; do i = 1, nb
   b = b + 1
   lo(:, b) = -0.6_R8P + [i - 1, j - 1, k - 1] * hb
   hi(:, b) = lo(:, b) + hb
enddo ; enddo ; enddo
call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.05_R8P)
print '(A,I0,A,I0,A,I0,A,I0)', 'blocks: cut ', count(codes == BOX_CUT), ', near ', count(codes == BOX_NEAR), &
                               ', inside ', count(codes == BOX_INSIDE), ', outside ', count(codes == BOX_OUTSIDE)
