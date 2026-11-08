use fossil,     only : surface_stl_object, curve_object, slice_surface, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE, &
                       STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: LEVELS = 6
type(surface_stl_object)  :: wing
type(curve_object)        :: section
real(R8P),    allocatable :: lo(:,:), hi(:,:), next_lo(:,:), next_hi(:,:), leaf_lo(:,:), leaf_hi(:,:)
real(R8P),    allocatable :: x(:), y(:)
integer(I4P), allocatable :: codes(:), leaf_code(:), tags(:,:), ftags(:,:)
real(R8P)                 :: side
integer(I4P)              :: level, n, nl, i, j, c, status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)                  ! the airfoil in the plane z = 0.5
lo = reshape([-0.5_R8P, -0.75_R8P], [2, 1]) ; hi = reshape([1._R8P, 0.75_R8P], [2, 1])  ! the root box; chord -0.25..0.75
allocate(leaf_lo(2, 4000), leaf_hi(2, 4000), leaf_code(4000))
nl = 0
do level = 0, LEVELS
   n = size(lo, 2)
   side = hi(1, 1) - lo(1, 1)
   allocate(codes(n))
   ! refine where |phi| may be below the children's diagonal (Matsushita & Aoki): cut or near within it
   call section%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.5_R8P * side * sqrt(2._R8P), status=status)
   if (status /= STATUS_OK) error stop 'box_classify_many failed'
   print '(A,I0,A,F8.5,A,I5,A,I4,A,I4,A,I4,A,I4)', 'level ', level, ', side ', side, ': boxes', n, ', cut', &
         count(codes == BOX_CUT), ', near', count(codes == BOX_NEAR), ', inside', count(codes == BOX_INSIDE), &
         ', outside', count(codes == BOX_OUTSIDE)
   allocate(next_lo(2, 4 * n), next_hi(2, 4 * n))
   c = 0
   do i = 1, n
      if (level < LEVELS .and. (codes(i) == BOX_CUT .or. codes(i) == BOX_NEAR)) then
         do j = 0, 3                                                    ! the four children
            c = c + 1
            next_lo(:, c) = lo(:, i) + 0.5_R8P * side * [mod(j, 2), j / 2]
            next_hi(:, c) = next_lo(:, c) + 0.5_R8P * side
         enddo
      else
         nl = nl + 1                                                    ! a leaf
         leaf_lo(:, nl) = lo(:, i) ; leaf_hi(:, nl) = hi(:, i) ; leaf_code(nl) = codes(i)
      endif
   enddo
   lo = next_lo(:, 1:c) ; hi = next_hi(:, 1:c)
   deallocate(codes, next_lo, next_hi)
enddo
print '(A,I0,A,I0,A,I0)', 'leaves ', nl, ', finest cut ', count(leaf_code(1:nl) == BOX_CUT), &
      ', uniform grid of the finest size ', 4**LEVELS
