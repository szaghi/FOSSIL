use fossil,     only : surface_stl_object, curve_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter   :: nb = 16, nc = 32, np = 64     ! 3-D boxes, 2-D boxes, 2-D points per direction
type(surface_stl_object)  :: surface
type(curve_object)        :: section
real(R8P)                 :: lo(3, nb**3), hi(3, nb**3), lo2(2, nc**2), hi2(2, nc**2), pt(2, np**2), d(np**2, 2), h
integer(I4P)              :: codes(nb**3, 2), codes2(nc**2, 2), i, j, k, b

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%enter_device
h = 1.2_R8P / nb
b = 0
do k = 1, nb ; do j = 1, nb ; do i = 1, nb
   b = b + 1
   lo(:, b) = -0.6_R8P + [i - 1, j - 1, k - 1] * h
   hi(:, b) = lo(:, b) + h
enddo ; enddo ; enddo
call surface%box_classify_many_device(lo=lo, hi=hi, codes=codes(:, 1), radius=0.05_R8P)  ! winding-number sign
call surface%box_classify_many(lo=lo, hi=hi, codes=codes(:, 2), radius=0.05_R8P)
print '(A,I0,A,I0)', 'surface boxes: device vs host mismatches ', count(codes(:, 1) /= codes(:, 2)), ' of ', nb**3
print '(A,I0,A,I0,A,I0,A,I0)', '  cut ', count(codes(:, 1) == BOX_CUT), ', near ', count(codes(:, 1) == BOX_NEAR), &
      ', inside ', count(codes(:, 1) == BOX_INSIDE), ', outside ', count(codes(:, 1) == BOX_OUTSIDE)
call surface%exit_device
