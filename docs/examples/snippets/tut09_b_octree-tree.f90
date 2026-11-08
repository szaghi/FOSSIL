use fossil, only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter   :: LEVELS = 5                     ! finest side 1.2 / 2^5 = 0.0375
type(surface_stl_object)  :: surface
real(R8P),    allocatable :: lo(:,:), hi(:,:), next_lo(:,:), next_hi(:,:)
integer(I4P), allocatable :: codes(:)
integer(I4P)              :: level, n, c, i, j, leaves, leaves_in, leaves_out
real(R8P)                 :: side

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
lo = reshape([-0.6_R8P, -0.6_R8P, -0.6_R8P], [3, 1]) ; hi = -lo          ! the root box, the domain [-0.6, 0.6]^3
leaves = 0 ; leaves_in = 0 ; leaves_out = 0
print '(A)', 'level   side     boxes    cut   near inside outside'
do level = 0, LEVELS
   n = size(lo, 2) ; side = hi(1, 1) - lo(1, 1)
   allocate(codes(n))
   call surface%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.5_R8P * side * sqrt(3._R8P))  ! child diagonal
   print '(I5,F8.4,I10,4I7)', level, side, n, count(codes == BOX_CUT), count(codes == BOX_NEAR), &
                              count(codes == BOX_INSIDE), count(codes == BOX_OUTSIDE)
   leaves_in = leaves_in + count(codes == BOX_INSIDE) ; leaves_out = leaves_out + count(codes == BOX_OUTSIDE)
   if (level == LEVELS) exit
   allocate(next_lo(3, 8 * n), next_hi(3, 8 * n))
   c = 0
   do i = 1, n
      if (codes(i) /= BOX_CUT .and. codes(i) /= BOX_NEAR) cycle            ! inside or outside: a leaf
      do j = 0, 7                                                          ! the eight children
         c = c + 1
         next_lo(:, c) = lo(:, i) + 0.5_R8P * side * [mod(j, 2), mod(j / 2, 2), j / 4]
         next_hi(:, c) = next_lo(:, c) + 0.5_R8P * side
      enddo
   enddo
   leaves = leaves + n - c / 8
   lo = next_lo(:, 1:c) ; hi = next_hi(:, 1:c)
   deallocate(codes, next_lo, next_hi)
enddo
leaves = leaves + n
print '(A,I0,A,I0,A,I0,A,I0)', 'leaves ', leaves, ' (inside ', leaves_in, ', outside ', leaves_out, &
                               '); uniform grid of the finest side: ', 8**LEVELS
