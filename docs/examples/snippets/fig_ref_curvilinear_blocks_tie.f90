!< Figure of the reference page "Curvilinear structured blocks", section Principles: the low-side tie rule on a sheared
!< 2-D block of 5 x 3 cells (u = i + j / 2, v = j: straight edges, dyadic coordinates, so every point below lies
!< exactly on an edge or a vertex). The midpoints of the interior edges and the interior vertices are located by
!< `locate_many`; a segment joins each point to the centre of the cell that owns it, the cell of lower index.
program fig_ref_curvilinear_blocks_tie
use fossil,     only : curvilinear_block_2d_t, LOC_INSIDE
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter      :: NI = 6, NJ = 4
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: u(NI, NJ), v(NI, NJ), pu(64), pv(64)
type(vector_R8P)             :: a(64), c(64)
integer(I4P), allocatable    :: cell(:,:), st(:)
real(R8P),    allocatable    :: xi(:,:)
integer(I4P)                 :: i, j, n, low, w
character(120)               :: caption

do j = 1, NJ
   do i = 1, NI
      u(i, j) = (i - 1) + 0.5_R8P * (j - 1) ; v(i, j) = j - 1
   enddo
enddo
call b%initialize(u, v)
call b%build_search
n = 0
do j = 1, NJ - 1   ! interior vertical edges (between cells i - 1 and i): midpoints
   do i = 2, NI - 1
      n = n + 1 ; pu(n) = 0.5_R8P * (u(i, j) + u(i, j + 1)) ; pv(n) = 0.5_R8P * (v(i, j) + v(i, j + 1))
   enddo
enddo
do j = 2, NJ - 1   ! interior horizontal edges
   do i = 1, NI - 1
      n = n + 1 ; pu(n) = 0.5_R8P * (u(i, j) + u(i + 1, j)) ; pv(n) = 0.5_R8P * (v(i, j) + v(i + 1, j))
   enddo
enddo
do j = 2, NJ - 1   ! interior vertices, shared by four cells
   do i = 2, NI - 1
      n = n + 1 ; pu(n) = u(i, j) ; pv(n) = v(i, j)
   enddo
enddo
call b%locate_many(pu(1:n), pv(1:n), cell, xi, st)
low = 0
do i = 1, n   ! the owner is the lowest cell among those touching the point: xi = 1 along every shared direction
   if (st(i) == LOC_INSIDE .and. all(xi(:, i) > 0._R8P)) low = low + 1
   a(i) = vector_R8P(pu(i), pv(i), 0.002_R8P)
   c(i) = centre(cell(1, i), cell(2, i))
enddo
print '(A,I0,A,I0,A,I0)', 'points on shared edges or vertices ', n, ', inside ', count(st == LOC_INSIDE), &
                          ', owned by the low side (xi > 0) ', low
write(caption, '(A,I0,A)') 'low-side rule: ', n, ' points on shared edges and vertices, each owned by the lower cell'
call fig_open(w, 'ref-curvilinear-blocks-tie', trim(caption), -90._R8P, 90._R8P, width=720_I4P, height=360_I4P)
do j = 1, NJ - 1
   do i = 1, NI - 1
      call fig_label(w, 'box', centre(i, j) + vector_R8P(-0.15_R8P, -0.3_R8P, 0.003_R8P), cid(i, j))
   enddo
enddo
call fig_segments(w, 'body', [((p(i, j), i = 1, NI - 1), j = 1, NJ)], [((p(i + 1, j), i = 1, NI - 1), j = 1, NJ)])
call fig_segments(w, 'body', [((p(i, j), i = 1, NI), j = 1, NJ - 1)], [((p(i, j + 1), i = 1, NI), j = 1, NJ - 1)])
call fig_segments(w, 'accent', a(1:n), c(1:n))
call fig_points(w, 'probe', a(1:n))
call fig_points(w, 'accent', c(1:n))
call fig_close(w)
contains
   function p(i, j) result(q)
   !< Vertex (i, j).
   integer(I4P), intent(in) :: i, j !< Indices.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(u(i, j), v(i, j), 0._R8P)
   endfunction p

   function centre(i, j) result(q)
   !< Vertex mean of cell (i, j), lifted.
   integer(I4P), intent(in) :: i, j !< Cell.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(0.25_R8P * (u(i, j) + u(i + 1, j) + u(i, j + 1) + u(i + 1, j + 1)), &
                  0.25_R8P * (v(i, j) + v(i + 1, j) + v(i, j + 1) + v(i + 1, j + 1)), 0.002_R8P)
   endfunction centre

   function cid(i, j) result(s)
   !< Label "(i,j)".
   integer(I4P), intent(in) :: i, j !< Cell.
   character(5)             :: s    !< Text.

   write(s, '(A,I1,A,I1,A)') '(', i, ',', j, ')'
   endfunction cid
endprogram fig_ref_curvilinear_blocks_tie
