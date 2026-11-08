!< Figure of the reference page "Curvilinear structured blocks", section Point location: a wavy 2-D block of 14 x 9 cells
!< with one vertex pushed across its neighbour (two folded cells, outlined in ember), and a lattice of points over its
!< right part and beyond, located by `locate_many`: LOC_INSIDE in teal, LOC_OUTSIDE in ember, LOC_WARPED (covered
!< only by a folded cell) in gold. For one point in the fold the candidate cells of the robust path are drawn (the
!< cells whose box contains it, tested in index order) and the cell that contains it is filled.
program fig_ref_curvilinear_blocks_locate
use fossil,     only : curvilinear_block_2d_t, CELL_VALID, LOC_INSIDE, LOC_OUTSIDE, LOC_WARPED
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_box, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter      :: NI = 15, NJ = 10, MU = 37, MV = 31
type(curvilinear_block_2d_t) :: b
real(R8P)                    :: u(NI, NJ), v(NI, NJ), pu(MU * MV), pv(MU * MV), q(2)
integer(I4P), allocatable    :: cell(:,:), st(:), valid(:,:), c1(:,:), s1(:)
real(R8P),    allocatable    :: xi(:,:), x1(:,:)
type(vector_R8P), allocatable :: pts(:)
integer(I4P)                 :: i, j, n, w, ncand
character(120)               :: caption

do j = 1, NJ
   do i = 1, NI
      u(i, j) = (i - 1) + 0.3_R8P * sin(0.6_R8P * (j - 1)) ; v(i, j) = (j - 1) + 0.3_R8P * sin(0.5_R8P * (i - 1))
   enddo
enddo
u(10, 5) = u(10, 5) + 1.35_R8P            ! pushed across the next vertex line: two cells fold
call b%initialize(u, v)
call b%validity(valid)
call b%build_search
n = 0
do j = 1, MV
   do i = 1, MU
      n = n + 1 ; pu(n) = 5.05_R8P + 0.27_R8P * (i - 1) + 0.09_R8P * mod(j, 2) ; pv(n) = 0.05_R8P + 0.27_R8P * (j - 1)
   enddo
enddo
call b%locate_many(pu, pv, cell, xi, st)
print '(A,I0,A,I0)', 'cells ', size(valid), ', not valid ', count(valid /= CELL_VALID)
print '(A,I0,A,I0,A,I0,A,I0)', 'points ', n, ': inside ', count(st == LOC_INSIDE), ', outside ', &
                               count(st == LOC_OUTSIDE), ', warped ', count(st == LOC_WARPED)
q = [10._R8P, 4.4_R8P] ; ncand = candidates(q(1), q(2))   ! the shown point, in the fold, inside a valid cell
call b%locate_many([q(1)], [q(2)], c1, x1, s1)
print '(A,2F6.2,A,I0,A,2I3,A,I0)', 'point ', q, ': candidate boxes ', ncand, ', cell ', c1(:, 1), ', status ', s1(1)

write(caption, '(A,I0,A,I0,A,I0,A)') 'located: ', count(st == LOC_INSIDE), ' inside (teal), ', count(st == LOC_OUTSIDE), &
                                     ' outside (ember), ', count(st == LOC_WARPED), ' warped (gold)'
call fig_open(w, 'ref-curvilinear-blocks-locate', trim(caption), -90._R8P, 90._R8P, width=760_I4P, height=520_I4P)
call fig_zoom(w, vector_R8P(10.6_R8P, 4.2_R8P, 0._R8P), 3.4_R8P)
call fig_segments(w, 'box', [((p(i, j), i = 1, NI - 1), j = 1, NJ)], [((p(i + 1, j), i = 1, NI - 1), j = 1, NJ)])
call fig_segments(w, 'box', [((p(i, j), i = 1, NI), j = 1, NJ - 1)], [((p(i, j + 1), i = 1, NI), j = 1, NJ - 1)])
do j = 1, NJ - 1
   do i = 1, NI - 1
      if (valid(i, j) /= CELL_VALID) call fig_segments(w, 'hole', [p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1)], &
                                                       [p(i + 1, j), p(i + 1, j + 1), p(i, j + 1), p(i, j)])
      if (box_has(i, j, q(1), q(2))) then
         call fig_box(w, 'accent', vector_R8P(minval(u(i:i + 1, j:j + 1)), minval(v(i:i + 1, j:j + 1)), 0.001_R8P), &
                      vector_R8P(maxval(u(i:i + 1, j:j + 1)), maxval(v(i:i + 1, j:j + 1)), 0.001_R8P))
      endif
   enddo
enddo
call b%locate_many([q(1)], [q(2)], c1, x1, s1)
i = c1(1, 1) ; j = c1(2, 1)
call fig_quad(w, 'ghost', p(i, j), p(i + 1, j), p(i + 1, j + 1), p(i, j + 1))
allocate(pts(n)) ; pts = [(vector_R8P(pu(i), pv(i), 0.002_R8P), i = 1, n)]
call fig_points(w, 'inside', pack(pts, st == LOC_INSIDE))
call fig_points(w, 'outside', pack(pts, st == LOC_OUTSIDE))
call fig_points(w, 'accent', pack(pts, st == LOC_WARPED))
call fig_points(w, 'foot', [vector_R8P(q(1), q(2), 0.003_R8P)])
write(caption, '(I0,A,I0,A,I0,A)') ncand, ' candidate boxes; only folded cell (', c1(1, 1), ',', c1(2, 1), &
                                   ') holds it: warped'
call fig_label(w, 'accent', vector_R8P(q(1) - 3.2_R8P, q(2) + 1.1_R8P, 0.003_R8P), trim(caption))
call fig_close(w)
contains
   function p(i, j) result(r)
   !< Vertex (i, j).
   integer(I4P), intent(in) :: i, j !< Indices.
   type(vector_R8P)         :: r    !< Point.

   r = vector_R8P(u(i, j), v(i, j), 0._R8P)
   endfunction p

   logical function box_has(i, j, x, y)
   !< The box of cell (i, j) contains (x, y).
   integer(I4P), intent(in) :: i, j !< Cell.
   real(R8P),    intent(in) :: x, y !< Point.

   box_has = x >= minval(u(i:i + 1, j:j + 1)) .and. x <= maxval(u(i:i + 1, j:j + 1)) .and. &
             y >= minval(v(i:i + 1, j:j + 1)) .and. y <= maxval(v(i:i + 1, j:j + 1))
   endfunction box_has

   integer(I4P) function candidates(x, y)
   !< Number of cells whose box contains (x, y).
   real(R8P), intent(in) :: x, y !< Point.
   integer(I4P)          :: ii, jj

   candidates = 0
   do jj = 1, NJ - 1
      do ii = 1, NI - 1
         if (box_has(ii, jj, x, y)) candidates = candidates + 1
      enddo
   enddo
   endfunction candidates
endprogram fig_ref_curvilinear_blocks_locate
