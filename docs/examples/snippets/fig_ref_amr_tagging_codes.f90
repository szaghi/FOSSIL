!< Figure of the reference page "AMR tagging", section "What it does": one layer of boxes through the dragon (the part
!< above it cut away), each classified by one `box_classify_many` batch (cut, near within half a box side, inside,
!< outside), and one box asked alone with `box_intersects_surface`.
program fig_ref_amr_tagging_codes
use fossil,     only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig, only : fig_open, fig_triangle, fig_quad, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: N = 40
type(surface_stl_object)  :: dragon
type(facet_object), pointer :: f
type(vector_R8P)          :: bmin, bmax
real(R8P),    allocatable :: lo(:,:), hi(:,:)
integer(I4P), allocatable :: codes(:)
real(R8P)                 :: h(2), z0, z1, val
integer(I4P)              :: i, j, m, status, u
character(120)            :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
h = [(bmax%x - bmin%x), (bmax%y - bmin%y)] * 1.1_R8P / N
z0 = 0.5_R8P * (bmin%z + bmax%z) - 0.5_R8P * h(1) ; z1 = z0 + h(1)   ! one layer of N x N boxes across the middle
allocate(lo(3, N * N), hi(3, N * N), codes(N * N))
m = 0
do j = 1, N
   do i = 1, N
      m = m + 1
      lo(:, m) = [bmin%x - 0.05_R8P * (bmax%x - bmin%x) + (i - 1) * h(1), &
                  bmin%y - 0.05_R8P * (bmax%y - bmin%y) + (j - 1) * h(2), z0]
      hi(:, m) = [lo(1, m) + h(1), lo(2, m) + h(2), z1]
   enddo
enddo
call dragon%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.5_R8P * h(1), status=status)
if (status /= STATUS_OK) error stop 'box_classify_many failed'
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'boxes ', N * N, ': cut ', count(codes == BOX_CUT), ', near ', &
      count(codes == BOX_NEAR), ', inside ', count(codes == BOX_INSIDE), ', outside ', count(codes == BOX_OUTSIDE)
m = findloc(codes, BOX_CUT, 1)
print '(A,I0,A,L1)', 'box ', m, ' alone: box_intersects_surface = ', dragon%box_intersects_surface(lo(:, m), hi(:, m))
m = findloc(codes, BOX_NEAR, 1)
print '(A,I0,A,L1)', 'box ', m, ' alone: box_intersects_surface = ', dragon%box_intersects_surface(lo(:, m), hi(:, m))

write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'one layer of boxes: ', count(codes == BOX_CUT), ' cut (cream), ', &
      count(codes == BOX_NEAR), ' near (ember), ', count(codes == BOX_INSIDE), ' inside (teal), ', &
      count(codes == BOX_OUTSIDE), ' outside'
call fig_open(u, 'ref-amr-tagging-codes', trim(caption), -60._R8P, 50._R8P, width=820_I4P, height=520_I4P)
do i = 1, dragon%get_facets_number()                    ! the dragon below the layer: the part above is cut away
   f => dragon%facet_at(i)
   if (max(f%vertex(1)%z, f%vertex(2)%z, f%vertex(3)%z) <= z0) &
      call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
do m = 1, N * N
   select case (codes(m))
   case (BOX_CUT)    ; val = 0._R8P
   case (BOX_NEAR)   ; val = 0.6_R8P
   case (BOX_INSIDE) ; val = -1._R8P
   case default      ; cycle
   endselect
   call fig_quad(u, 'field', vector_R8P(lo(1, m), lo(2, m), z0), vector_R8P(hi(1, m), lo(2, m), z0), &
                 vector_R8P(hi(1, m), hi(2, m), z0), vector_R8P(lo(1, m), hi(2, m), z0), val)
enddo
call fig_box(u, 'box', vector_R8P(lo(1, 1), lo(2, 1), z0), vector_R8P(hi(1, N * N), hi(2, N * N), z1))
call fig_close(u)
endprogram fig_ref_amr_tagging_codes
