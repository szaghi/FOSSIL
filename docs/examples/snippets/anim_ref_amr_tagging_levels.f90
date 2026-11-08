!< Animation of the reference page "AMR tagging": a quadtree over the section of the dragon (the plane z through the
!< middle of its box), refined level by level. Each frame classifies the boxes of one level with one
!< `box_classify_many` batch (radius = the children's diagonal): cut (cream) and near (ember) boxes are split into the
!< next level; inside (teal) and outside (outlined) boxes become leaves and stay.
program anim_ref_amr_tagging_levels
use fossil,     only : surface_stl_object, curve_object, slice_surface, BOX_CUT, BOX_NEAR, BOX_INSIDE, STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_quad, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: LEVELS = 6, MAXLEAF = 20000
type(surface_stl_object)  :: dragon
type(curve_object)        :: section
type(vector_R8P)          :: bmin, bmax
type(vector_R8P), allocatable :: sa(:), sb(:)
real(R8P),    allocatable :: lo(:,:), hi(:,:), next_lo(:,:), next_hi(:,:), leaf_lo(:,:), leaf_hi(:,:)
integer(I4P), allocatable :: codes(:), leaf_code(:)
real(R8P)                 :: side, c(2), r(2)
integer(I4P)              :: level, n, nl, i, j, m, status, u
character(120)            :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call slice_surface(dragon, 3_I4P, 0.5_R8P * (bmin%z + bmax%z), section)
print '(A,I0,A,I0,A)', 'section: ', section%get_loops_number(), ' loops, ', section%get_vertices_number(), ' vertices'
c = [0.5_R8P * (bmin%x + bmax%x), 0.5_R8P * (bmin%y + bmax%y)]
side = 1.05_R8P * max(bmax%x - bmin%x, bmax%y - bmin%y)
lo = reshape(c - 0.5_R8P * side, [2, 1]) ; hi = reshape(c + 0.5_R8P * side, [2, 1])
allocate(leaf_lo(2, MAXLEAF), leaf_hi(2, MAXLEAF), leaf_code(MAXLEAF)) ; nl = 0
allocate(sa(section%get_vertices_number()), sb(section%get_vertices_number())) ; m = 0
do i = 1, section%get_vertices_number()
   if (section%get_next(i) == 0_I4P) cycle
   r = section%get_vertex(i) ; m = m + 1 ; sa(m) = p(r(1), r(2), 3.e-3_R8P)
   r = section%get_vertex(section%get_next(i)) ; sb(m) = p(r(1), r(2), 3.e-3_R8P)
enddo

call fig_open(u, 'anim-ref-amr-tagging-levels', '', -90._R8P, 90._R8P, width=820_I4P, height=520_I4P)
call fig_zoom(u, p(c(1), c(2)), 0.3_R8P * side)
call fig_anim(u, 1.2_R8P)
call fig_segments(u, 'accent', sa(1:m), sb(1:m))
do level = 0, LEVELS
   n = size(lo, 2) ; side = hi(1, 1) - lo(1, 1)
   allocate(codes(n))
   call section%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.5_R8P * side * sqrt(2._R8P), status=status)
   if (status /= STATUS_OK) error stop 'box_classify_many failed'
   write(caption, '(A,I0,A,I0,A,I0,A,I0,A,I0)') 'level ', level, ': ', n, ' boxes, cut ', count(codes == BOX_CUT), &
         ', near ', count(codes == BOX_NEAR), ', inside ', count(codes == BOX_INSIDE)
   if (level < LEVELS) then
      write(caption, '(A,A,I0)') trim(caption), '; leaves so far ', nl + count(codes /= BOX_CUT .and. codes /= BOX_NEAR)
   else
      write(caption, '(A,A,I0)') trim(caption), '; all leaves ', nl + n
   endif
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption), seconds=merge(2.8_R8P, 1.2_R8P, level == LEVELS))
   call draw_leaves
   do i = 1, n
      call draw_box(lo(:, i), hi(:, i), codes(i))
   enddo
   allocate(next_lo(2, 4 * n), next_hi(2, 4 * n)) ; j = 0
   do i = 1, n
      if (level < LEVELS .and. (codes(i) == BOX_CUT .or. codes(i) == BOX_NEAR)) then
         do m = 0, 3
            j = j + 1
            next_lo(:, j) = lo(:, i) + 0.5_R8P * side * [mod(m, 2), m / 2]
            next_hi(:, j) = next_lo(:, j) + 0.5_R8P * side
         enddo
      else
         nl = nl + 1 ; leaf_lo(:, nl) = lo(:, i) ; leaf_hi(:, nl) = hi(:, i) ; leaf_code(nl) = codes(i)
      endif
   enddo
   lo = next_lo(:, 1:j) ; hi = next_hi(:, 1:j)
   deallocate(codes, next_lo, next_hi)
enddo
print '(A,I0,A,I0)', 'leaves ', nl, ', uniform grid of the finest size ', 4**LEVELS
call fig_close(u)
contains
   subroutine draw_leaves
   !< The leaves of the previous levels: inside filled teal, outside outlined.
   integer(I4P) :: k

   do k = 1, nl
      if (leaf_code(k) == BOX_INSIDE) then
         call fig_quad(u, 'field', p(leaf_lo(1, k), leaf_lo(2, k)), p(leaf_hi(1, k), leaf_lo(2, k)), &
                       p(leaf_hi(1, k), leaf_hi(2, k)), p(leaf_lo(1, k), leaf_hi(2, k)), -1._R8P)
      endif
   enddo
   call outline(leaf_lo(:, 1:nl), leaf_hi(:, 1:nl))
   endsubroutine draw_leaves

   subroutine draw_box(a, b, code)
   !< A box of the current level: cut cream, near ember, inside teal, outside outlined only.
   real(R8P),    intent(in) :: a(2), b(2)
   integer(I4P), intent(in) :: code

   select case (code)
   case (BOX_CUT)    ; call fig_quad(u, 'field', p(a(1), a(2)), p(b(1), a(2)), p(b(1), b(2)), p(a(1), b(2)), 0._R8P)
   case (BOX_NEAR)   ; call fig_quad(u, 'field', p(a(1), a(2)), p(b(1), a(2)), p(b(1), b(2)), p(a(1), b(2)), 0.6_R8P)
   case (BOX_INSIDE) ; call fig_quad(u, 'field', p(a(1), a(2)), p(b(1), a(2)), p(b(1), b(2)), p(a(1), b(2)), -1._R8P)
   endselect
   call outline(reshape(a, [2, 1]), reshape(b, [2, 1]))
   endsubroutine draw_box

   subroutine outline(a, b)
   !< The sides of the boxes [a(:, k), b(:, k)], just above the plane.
   real(R8P), intent(in)         :: a(:,:), b(:,:)
   type(vector_R8P), allocatable :: e0(:), e1(:)
   integer(I4P)                  :: k

   if (size(a, 2) == 0) return
   allocate(e0(4 * size(a, 2)), e1(4 * size(a, 2)))
   do k = 1, size(a, 2)
      e0(4 * k - 3) = p(a(1, k), a(2, k), 1.e-3_R8P) ; e1(4 * k - 3) = p(b(1, k), a(2, k), 1.e-3_R8P)
      e0(4 * k - 2) = p(b(1, k), a(2, k), 1.e-3_R8P) ; e1(4 * k - 2) = p(b(1, k), b(2, k), 1.e-3_R8P)
      e0(4 * k - 1) = p(b(1, k), b(2, k), 1.e-3_R8P) ; e1(4 * k - 1) = p(a(1, k), b(2, k), 1.e-3_R8P)
      e0(4 * k)     = p(a(1, k), b(2, k), 1.e-3_R8P) ; e1(4 * k)     = p(a(1, k), a(2, k), 1.e-3_R8P)
   enddo
   call fig_segments(u, 'box', e0, e1)
   endsubroutine outline

   function p(px, py, dz) result(q)
   !< The point (px, py, dz) of the figure, dz = 0 by default.
   real(R8P), intent(in)           :: px, py
   real(R8P), intent(in), optional :: dz
   type(vector_R8P)                :: q

   q = vector_R8P(px, py, 0._R8P)
   if (present(dz)) q%z = dz
   endfunction p
endprogram anim_ref_amr_tagging_levels
