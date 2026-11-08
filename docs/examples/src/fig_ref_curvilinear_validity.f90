!run figref-curvilinear-validity fig_ref_curvilinear_validity
!fig ref-curvilinear-validity
!< Figure of the reference page "Curvilinear structured blocks": the exact validity of a 2-D block of 20 x 12 cells,
!< a wavy grid with one vertex pushed across its neighbour (a fold) and one edge collapsed (two coincident vertices).
!< Valid cells are coloured by their minimum Jacobian over the mean, from cream (1) to sepia (0); invalid cells are
!< ember, singular ones teal.
program fig_ref_curvilinear_validity
use fossil,     only : curvilinear_block_2d_t, CELL_VALID, CELL_SINGULAR, CELL_INVALID
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: u(21, 13), v(21, 13)
integer(I4P), allocatable    :: st(:,:)
real(R8P),    allocatable    :: sk(:,:), orth(:,:), asp(:,:), mj(:,:)
integer(I4P)                 :: i, j, w

do j = 1, 13
   do i = 1, 21
      u(i, j) = (i - 1) + 0.25_R8P * sin(0.5_R8P * (j - 1)) ; v(i, j) = (j - 1) + 0.25_R8P * sin(0.4_R8P * (i - 1))
   enddo
enddo
u(15, 7) = u(15, 7) + 1.4_R8P                     ! pushed across the next vertex line: a fold
u(5, 4) = u(6, 4) ; v(5, 4) = v(6, 4)             ! an edge collapsed: two coincident vertices
call b%initialize(u, v)
call b%validity(st)
call b%quality(sk, orth, asp, mj)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells ', size(st), ': valid ', count(st == CELL_VALID), ', singular ', &
                               count(st == CELL_SINGULAR), ', invalid ', count(st == CELL_INVALID)
print '(A,F7.4,A,F7.4)', 'min J / mean J of the valid cells: lowest ', minval(mj, st == CELL_VALID), ', highest ', &
                         maxval(mj, st == CELL_VALID)
call fig_open(w, 'ref-curvilinear-validity', 'exact validity: min J / mean (cream 1, sepia 0), invalid ember, singular teal', &
              -90._R8P, 90._R8P, height=360_I4P)
do j = 1, 12
   do i = 1, 20
      select case (st(i, j))
      case (CELL_VALID)
         call fig_quad(w, 'frac', p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1), 1._R8P - max(0._R8P, mj(i, j)))
      case (CELL_SINGULAR)
         call fig_quad(w, 'probe', p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1))
      case default
         call fig_quad(w, 'patch', p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1))
      endselect
   enddo
enddo
call fig_segments(w, 'cell', [((p(i, j), i = 1, 20), j = 1, 13)], [((p(i + 1, j), i = 1, 20), j = 1, 13)])
call fig_segments(w, 'cell', [((p(i, j), i = 1, 21), j = 1, 12)], [((p(i, j + 1), i = 1, 21), j = 1, 12)])
call fig_close(w)
contains
   function p(i, j) result(q)
   !< Vertex (i, j), lifted a little toward the viewer for the lines.
   integer, intent(in) :: i, j !< Indices.
   type(vector_R8P)    :: q    !< Point.

   q = vector_R8P(u(i, j), v(i, j), 0._R8P)
   endfunction p
endprogram fig_ref_curvilinear_validity
