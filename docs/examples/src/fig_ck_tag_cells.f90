!run figck-tag-cells fig_ck_tag_cells
!fig ck-tag-cells
!< Figure of the cookbook recipe "Tag the cells near the wall": the codes box_classify_many gives the recipe's 10^3
!< cells of edge 0.17 on [-0.35, 1.35]^3 around cube-01 (radius 0.1), the layer of cells z in [0.33, 0.50] seen from
!< above. Colours: cut ember, near pale ember, inside teal, outside cream, each kind labelled on
!< one cell; the section of the cube (the wall) is drawn over.
program fig_ck_tag_cells
use fossil,     only : surface_stl_object, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use fossil_fig, only : fig_open, fig_quad, fig_segment, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: kl = 4                                ! the layer drawn (z in [0.33, 0.50])
type(surface_stl_object) :: cube
type(vector_R8P)         :: a(4), b(4), o
real(R8P)                :: lo(3, 1000), hi(3, 1000), h, zl, g, v
integer(I4P)             :: codes(1000), i, j, k, n, u
character(80)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the recipe's cells and call
h = 0.17_R8P
n = 0
do k = 0, 9 ; do j = 0, 9 ; do i = 0, 9
   n = n + 1
   lo(:, n) = -0.35_R8P + h * [i, j, k] ; hi(:, n) = lo(:, n) + h
enddo ; enddo ; enddo
call cube%box_classify_many(lo=lo, hi=hi, codes=codes, radius=0.1_R8P)

zl = -0.35_R8P + h * (kl + 1)                                     ! top of the layer
write(caption, '(A,F4.2,A,F4.2,A)') 'box_classify_many(radius=0.1), the cells z in [', zl - h, ', ', zl, ']'
call fig_open(u, 'ck-tag-cells', trim(caption), -90._R8P, 90._R8P, height=520_I4P)
g = 0.006_R8P                                                     ! a gap between cells, so that each one reads
do j = 0, 9 ; do i = 0, 9
   n = 1 + i + 10 * j + 100 * kl
   select case (codes(n))
   case (BOX_CUT)     ; v = 1._R8P
   case (BOX_NEAR)    ; v = 0.4_R8P
   case (BOX_INSIDE)  ; v = -0.8_R8P
   case (BOX_OUTSIDE) ; v = 0._R8P
   endselect
   call fig_quad(u, 'field', p(lo(1, n) + g, lo(2, n) + g), p(hi(1, n) - g, lo(2, n) + g), &
                 p(hi(1, n) - g, hi(2, n) - g), p(lo(1, n) + g, hi(2, n) - g), value=v)
enddo ; enddo
a = [p(0._R8P, 0._R8P), p(1._R8P, 0._R8P), p(1._R8P, 1._R8P), p(0._R8P, 1._R8P)]
b = [a(2), a(3), a(4), a(1)]
call fig_segments(u, 'accent', a, b)                              ! the section of the cube
do i = 8, 9                                                       ! legend: the code of one cell of each kind
   n = 1 + i + 10 * i + 100 * kl ; call fig_label(u, 'accent', p(lo(1, n) + 0.01_R8P, lo(2, n) + 0.07_R8P), name(codes(n)))
enddo
n = 1 + 2 + 10 * 5 + 100 * kl ; call fig_label(u, 'accent', p(lo(1, n) + 0.01_R8P, lo(2, n) + 0.07_R8P), name(codes(n)))
n = 1 + 4 + 10 * 4 + 100 * kl ; call fig_label(u, 'accent', p(lo(1, n) + 0.01_R8P, lo(2, n) + 0.07_R8P), name(codes(n)))
o = p(-0.32_R8P, -0.32_R8P)                                       ! axes glyph at a corner
call fig_segment(u, 'accent', o, o + 0.2_R8P * ex_R8P) ; call fig_label(u, 'accent', o + 0.2_R8P * ex_R8P, 'x')
call fig_segment(u, 'accent', o, o + 0.2_R8P * ey_R8P) ; call fig_label(u, 'accent', o + 0.2_R8P * ey_R8P, 'y')
call fig_close(u)
contains
   function name(code) result(text)
   !< Name of a box code.
   integer(I4P), intent(in) :: code !< Code.
   character(:), allocatable :: text !< Name.

   select case (code)
   case (BOX_CUT)     ; text = 'cut'
   case (BOX_NEAR)    ; text = 'near'
   case (BOX_INSIDE)  ; text = 'inside'
   case default       ; text = 'outside'
   endselect
   endfunction name

   function p(xx, yy) result(q)
   !< A point of the top plane of the layer.
   real(R8P), intent(in) :: xx, yy !< In-plane coordinates.
   type(vector_R8P)      :: q      !< Point.

   q = xx * ex_R8P + yy * ey_R8P + zl * ez_R8P
   endfunction p
endprogram fig_ck_tag_cells
