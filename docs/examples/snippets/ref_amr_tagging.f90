!< Reference page "AMR tagging": a quadtree refined around the NACA 0012 section by box classification (refine a box
!< that is cut, or near within its children's diagonal), then signed tag bands and corner tags on a stretched block.
program ref_amr_tagging
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segment, fig_close
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
! a block of 48 x 24 cells around the airfoil, stretched away from y = 0 by 1.1
x = [(-0.35_R8P + 1.2_R8P * i / 48._R8P, i = 0, 48)]
allocate(y(25)) ; y(13) = 0._R8P
do j = 1, 12
   y(13 + j) = y(12 + j) + 0.006_R8P * 1.1_R8P**(j - 1)
   y(13 - j) = -y(13 + j)
enddo
allocate(tags(48, 24), ftags(48, 24))
call section%tag_block(x=x, y=y, band=3, tags=tags, feature_band=2, feature_tags=ftags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
do i = -4, 4
   print '(A,I3,A,I4)', 'tag', i, ': cells', count(tags == i)
enddo
print '(A,I0,A,I0)', 'cells holding a corner ', count(ftags == 0), ', within 2 cells of one ', count(ftags <= 2)
call draw_tree
call draw_tags
contains
   subroutine draw_tree
   !< The leaves of the quadtree coloured by their code, and the airfoil.
   integer(I4P)   :: u, k
   real(R8P)      :: value
   character(100) :: caption

   write(caption, '(A,I0,A,I0,A)') 'quadtree of ', LEVELS, ' levels, ', nl, &
                                   ' leaves: cut cream, near ember, inside teal'
   call fig_open(u, 'ref-amr-tagging-quadtree', trim(caption), -90._R8P, 90._R8P, height=500_I4P)
   do k = 1, nl
      select case (leaf_code(k))
      case (BOX_CUT)    ; value = 0._R8P
      case (BOX_NEAR)   ; value = 0.6_R8P
      case (BOX_INSIDE) ; value = -1._R8P
      case default      ; cycle
      endselect
      call fig_quad(u, 'field', p(leaf_lo(1, k), leaf_lo(2, k)), p(leaf_hi(1, k), leaf_lo(2, k)), &
                    p(leaf_hi(1, k), leaf_hi(2, k)), p(leaf_lo(1, k), leaf_hi(2, k)), value)
   enddo
   do k = 1, nl
      if (leaf_code(k) == BOX_OUTSIDE) call outline(u, 'box', leaf_lo(:, k), leaf_hi(:, k))
   enddo
   call draw_section(u)
   call fig_close(u)
   endsubroutine draw_tree

   subroutine draw_tags
   !< The stretched block coloured by the signed tags, the corner cells outlined, and the airfoil.
   integer(I4P)   :: u, i, j
   character(100) :: caption

   write(caption, '(A)') 'tag_block, band 3: cut cream, inside teal, outside ember; corner band outlined'
   call fig_open(u, 'ref-amr-tagging-tags', trim(caption), -90._R8P, 90._R8P, height=420_I4P)
   do j = 1, 24 ; do i = 1, 48
      call fig_quad(u, 'field', p(x(i), y(j)), p(x(i + 1), y(j)), p(x(i + 1), y(j + 1)), p(x(i), y(j + 1)), &
                    real(tags(i, j), R8P) / 4._R8P)
   enddo ; enddo
   do j = 1, 24 ; do i = 1, 48
      if (ftags(i, j) <= 2) call outline(u, 'probe', [x(i), y(j)], [x(i + 1), y(j + 1)])
   enddo ; enddo
   call draw_section(u)
   call fig_close(u)
   endsubroutine draw_tags

   subroutine outline(u, style, a, b)
   !< The four sides of the box [a, b], just above the plane.
   integer(I4P), intent(in) :: u
   character(*), intent(in) :: style
   real(R8P),    intent(in) :: a(2), b(2)

   call fig_segment(u, style, p(a(1), a(2), 1.e-3_R8P), p(b(1), a(2), 1.e-3_R8P))
   call fig_segment(u, style, p(b(1), a(2), 1.e-3_R8P), p(b(1), b(2), 1.e-3_R8P))
   call fig_segment(u, style, p(b(1), b(2), 1.e-3_R8P), p(a(1), b(2), 1.e-3_R8P))
   call fig_segment(u, style, p(a(1), b(2), 1.e-3_R8P), p(a(1), a(2), 1.e-3_R8P))
   endsubroutine outline

   subroutine draw_section(u)
   !< The airfoil, just above the plane.
   integer(I4P), intent(in) :: u
   real(R8P)                :: a(2), b(2)
   integer(I4P)             :: s

   do s = 1, section%get_vertices_number()
      if (section%get_next(s) == 0_I4P) cycle
      a = section%get_vertex(s) ; b = section%get_vertex(section%get_next(s))
      call fig_segment(u, 'accent', p(a(1), a(2), 2.e-3_R8P), p(b(1), b(2), 2.e-3_R8P))
   enddo
   endsubroutine draw_section

   function p(a, b, c) result(q)
   !< The point (a, b, c) of the figure, c = 0 by default.
   real(R8P), intent(in)           :: a, b
   real(R8P), intent(in), optional :: c
   type(vector_R8P)                :: q

   q = vector_R8P(a, b, 0._R8P)
   if (present(c)) q%z = c
   endfunction p
endprogram ref_amr_tagging
