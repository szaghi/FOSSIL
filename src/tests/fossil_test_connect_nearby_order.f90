!< FOSSIL test: `connect_nearby_vertices` is order-independent (audit #14 S8).
!<
!< Loads the same fan-of-triangles geometry twice — natural facet order and reversed —
!< runs `sanitize` on each, then verifies that the full vertex sets are identical
!< (as sorted multisets) to within a tight numerical tolerance.
!<
!< This test was promoted from `fossil_xfail_connect_nearby_order` once the
!< union-find rewrite of `connect_nearby_vertices` removed the order-dependence.

program fossil_test_connect_nearby_order

use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none

type(surface_stl_object)  :: surface_natural   !< Natural facet order.
type(surface_stl_object)  :: surface_reversed  !< Reversed facet order.
real(R8P), parameter      :: TOL = 1.0e-10_R8P !< Position agreement tolerance.
integer(I4P)              :: n_facets, nv_total
real(R8P), allocatable    :: vn(:,:), vr(:,:)  !< Vertex arrays (3, nv_total).
integer(I4P)              :: f, v, k
real(R8P)                 :: max_diff
logical                   :: all_passed
type(facet_object), pointer :: fp

call surface_natural %load_from_file(file_name='src/tests/s8_fan.stl',          is_ascii=.true.)
call surface_reversed%load_from_file(file_name='src/tests/s8_fan_reversed.stl', is_ascii=.true.)
call surface_natural %sanitize
call surface_reversed%sanitize

n_facets = surface_natural%get_facets_number()
nv_total = 3 * n_facets
allocate(vn(3, nv_total), vr(3, nv_total))

! Collect all vertex positions.
do f = 1, n_facets
   fp => surface_natural%facet_at(f)
   do v = 1, 3
      k = (f - 1) * 3 + v
      vn(1, k) = fp%vertex(v)%x
      vn(2, k) = fp%vertex(v)%y
      vn(3, k) = fp%vertex(v)%z
   enddo
   fp => surface_reversed%facet_at(f)
   do v = 1, 3
      k = (f - 1) * 3 + v
      vr(1, k) = fp%vertex(v)%x
      vr(2, k) = fp%vertex(v)%y
      vr(3, k) = fp%vertex(v)%z
   enddo
enddo

! Sort both vertex arrays by (x, y, z) lexicographic order.
call lex_sort(vn, nv_total)
call lex_sort(vr, nv_total)

! Compare sorted vertex sets element-wise.
max_diff = 0._R8P
do k = 1, nv_total
   max_diff = max(max_diff, &
      sqrt((vn(1,k)-vr(1,k))**2 + (vn(2,k)-vr(2,k))**2 + (vn(3,k)-vr(3,k))**2))
enddo

print '(A,ES12.4)', 'max vertex position diff between orderings: ', max_diff
all_passed = (max_diff <= TOL)
print '(A,L1)', 'Are all tests passed? ', all_passed
if (.not. all_passed) error stop 1
contains
   subroutine lex_sort(v, n)
   !< In-place insertion sort of columns of v(3,n) by (x,y,z) lexicographic order.
   real(R8P),    intent(inout) :: v(3, n)
   integer(I4P), intent(in)    :: n
   real(R8P)                   :: tmp(3)
   integer(I4P)                :: i, j

   do i = 2, n
      tmp = v(:, i)
      j = i - 1
      do while (j >= 1 .and. lex_gt(v(:, j), tmp))
         v(:, j + 1) = v(:, j)
         j = j - 1
      enddo
      v(:, j + 1) = tmp
   enddo
   endsubroutine lex_sort

   pure logical function lex_gt(a, b)
   !< True if (ax,ay,az) > (bx,by,bz) lexicographically.
   real(R8P), intent(in) :: a(3), b(3)
   if (a(1) /= b(1)) then
      lex_gt = a(1) > b(1)
   elseif (a(2) /= b(2)) then
      lex_gt = a(2) > b(2)
   else
      lex_gt = a(3) > b(3)
   endif
   endfunction lex_gt
endprogram fossil_test_connect_nearby_order
