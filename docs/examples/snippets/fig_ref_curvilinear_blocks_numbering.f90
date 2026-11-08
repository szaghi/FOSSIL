!< Figure of the reference page "Curvilinear structured blocks", section Data model: a block of 6 x 4 x 3 vertices (a
!< box bent by a sine along i), its grid lines, the index directions from vertex (1, 1, 1), the
!< six faces numbered 1 to 6 (-i, +i, -j, +j, -k, +k) and cell (4, 1, 2) filled, with its eight vertices i..i+1,
!< j..j+1, k..k+1.
program fig_ref_curvilinear_blocks_numbering
use fossil,     only : curvilinear_block_t
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter   :: PI = acos(-1._R8P)
integer(I4P), parameter   :: NI = 6, NJ = 4, NK = 3, CI = 4, CJ = 1, CK = 2
character(7), parameter   :: FACE_NAME(6) = ['1 (-i) ', '2 (+i) ', '3 (-j) ', '4 (+j) ', '5 (-k) ', '6 (+k) ']
type(curvilinear_block_t) :: b
real(R8P)                 :: x(NI, NJ, NK), y(NI, NJ, NK), z(NI, NJ, NK), t, r
real(R8P), allocatable    :: vol(:,:,:)
type(vector_R8P)          :: g, fc
integer(I4P)              :: i, j, k, f, w, lo(3, 6), hi(3, 6)

do k = 1, NK
   do j = 1, NJ
      do i = 1, NI
         t = PI * (i - 1) / (NI - 1)
         x(i, j, k) = 0.8_R8P * (i - 1) ; y(i, j, k) = 0.8_R8P * (j - 1) + 0.25_R8P * sin(t)
         z(i, j, k) = 0.7_R8P * (k - 1) + 0.35_R8P * sin(t) * (1._R8P + 0.3_R8P * (j - 1))
      enddo
   enddo
enddo
call b%initialize(x, y, z)
call b%geometry(volume=vol)
print '(A,I0,A,I0,A,I0,A,I0,A)', 'vertices ', NI, ' x ', NJ, ' x ', NK, ', cells ', b%cells_number(), &
                                 ' (5 x 3 x 2)'
print '(A,F9.6,A,F9.6,A)', 'volume of cell (4, 1, 2) ', vol(CI, CJ, CK), ', block ', sum(vol), ' (right-handed)'

call fig_open(w, 'ref-curvilinear-blocks-numbering', '6 x 4 x 3 vertices; faces 1-6; cell (4,1,2) = vertices 4..5, '// &
              '1..2, 2..3', -62._R8P, 24._R8P, width=760_I4P, height=520_I4P)
call fig_zoom(w, vector_R8P(2._R8P, 1.2_R8P, 0.9_R8P), 2.45_R8P)
call fig_segments(w, 'box', [(((p(i, j, k), i = 1, NI - 1), j = 1, NJ), k = 1, NK)], &
                            [(((p(i + 1, j, k), i = 1, NI - 1), j = 1, NJ), k = 1, NK)])
call fig_segments(w, 'box', [(((p(i, j, k), i = 1, NI), j = 1, NJ - 1), k = 1, NK)], &
                            [(((p(i, j + 1, k), i = 1, NI), j = 1, NJ - 1), k = 1, NK)])
call fig_segments(w, 'box', [(((p(i, j, k), i = 1, NI), j = 1, NJ), k = 1, NK - 1)], &
                            [(((p(i, j, k + 1), i = 1, NI), j = 1, NJ), k = 1, NK - 1)])
i = CI ; j = CJ ; k = CK   ! the cell: its six faces
call fig_quad(w, 'patch', p(i, j, k), p(i, j + 1, k), p(i, j + 1, k + 1), p(i, j, k + 1))
call fig_quad(w, 'patch', p(i + 1, j, k), p(i + 1, j + 1, k), p(i + 1, j + 1, k + 1), p(i + 1, j, k + 1))
call fig_quad(w, 'patch', p(i, j, k), p(i + 1, j, k), p(i + 1, j, k + 1), p(i, j, k + 1))
call fig_quad(w, 'patch', p(i, j + 1, k), p(i + 1, j + 1, k), p(i + 1, j + 1, k + 1), p(i, j + 1, k + 1))
call fig_quad(w, 'patch', p(i, j, k), p(i + 1, j, k), p(i + 1, j + 1, k), p(i, j + 1, k))
call fig_quad(w, 'patch', p(i, j, k + 1), p(i + 1, j, k + 1), p(i + 1, j + 1, k + 1), p(i, j + 1, k + 1))
call fig_points(w, 'foot', [(((p(i, j, k), i = CI, CI + 1), j = CJ, CJ + 1), k = CK, CK + 1)])
call fig_label(w, 'accent', p(CI, CJ, CK + 1) + vector_R8P(-0.3_R8P, -0.3_R8P, 0.25_R8P), 'cell (4,1,2)')
! index directions from vertex (1, 1, 1), along the first grid lines
call fig_segments(w, 'probe', [(p(i, 1, 1), i = 1, NI - 1)], [(p(i + 1, 1, 1), i = 1, NI - 1)])
call fig_segments(w, 'probe', [(p(1, j, 1), j = 1, NJ - 1)], [(p(1, j + 1, 1), j = 1, NJ - 1)])
call fig_segments(w, 'probe', [(p(1, 1, k), k = 1, NK - 1)], [(p(1, 1, k + 1), k = 1, NK - 1)])
call fig_points(w, 'probe', [p(1, 1, 1)])
call fig_label(w, 'probe', p(1, 1, 1) + vector_R8P(-0.2_R8P, -0.2_R8P, -0.2_R8P), 'vertex (1,1,1)')
call fig_label(w, 'probe', p(NI, 1, 1) + vector_R8P(0.2_R8P, -0.1_R8P, -0.1_R8P), 'i')
call fig_label(w, 'probe', p(1, NJ, 1) + vector_R8P(-0.1_R8P, 0.2_R8P, -0.1_R8P), 'j')
call fig_label(w, 'probe', p(1, 1, NK) + vector_R8P(-0.1_R8P, -0.1_R8P, 0.2_R8P), 'k')
lo(:, 1) = [1, 1, 1] ; hi(:, 1) = [1, NJ, NK] ; lo(:, 2) = [NI, 1, 1] ; hi(:, 2) = [NI, NJ, NK]
lo(:, 3) = [1, 1, 1] ; hi(:, 3) = [NI, 1, NK] ; lo(:, 4) = [1, NJ, 1] ; hi(:, 4) = [NI, NJ, NK]
lo(:, 5) = [1, 1, 1] ; hi(:, 5) = [NI, NJ, 1] ; lo(:, 6) = [1, 1, NK] ; hi(:, 6) = [NI, NJ, NK]
! the six faces: labels at their vertex means, pushed outwards (along i: tangentially, j: radially, k: vertically)
do f = 1, 6
   fc = vector_R8P(sum(x(lo(1, f):hi(1, f), lo(2, f):hi(2, f), lo(3, f):hi(3, f))), &
                   sum(y(lo(1, f):hi(1, f), lo(2, f):hi(2, f), lo(3, f):hi(3, f))), &
                   sum(z(lo(1, f):hi(1, f), lo(2, f):hi(2, f), lo(3, f):hi(3, f)))) / &
        real(product(hi(:, f) - lo(:, f) + 1), R8P)
   g = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   select case (f)
   case (1, 2)
      g%x = merge(-1._R8P, 1._R8P, f == 1)
   case (3, 4)
      g%y = merge(-1._R8P, 1._R8P, f == 3) ; if (f == 3) g%x = -1.2_R8P
   case default
      g%z = merge(-1._R8P, 1._R8P, f == 5)
   endselect
   call fig_label(w, 'accent', fc + 0.5_R8P * g, 'face '//trim(FACE_NAME(f)))
enddo
call fig_close(w)
contains
   function p(i, j, k) result(q)
   !< Vertex (i, j, k).
   integer(I4P), intent(in) :: i, j, k !< Indices.
   type(vector_R8P)         :: q       !< Point.

   q = vector_R8P(x(i, j, k), y(i, j, k), z(i, j, k))
   endfunction p
endprogram fig_ref_curvilinear_blocks_numbering
