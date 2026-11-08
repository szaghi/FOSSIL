!< Figure of the reference page "Cut cells: solver support", section Principles: what the graph is built on, on 10 x 6
!< cells of 0.1. A square aligned with the grid lines (its right side on the line x = 0.4, its top side on y = 0.4)
!< leaves pieces of zero area where its wall is the high edge of a cell; a thin tilted plate splits cells into two
!< pieces; an open polyline leaves its cells unclosed, so they have no pieces (untrusted). Every piece is a dot at its
!< centroid: teal for a cell with one piece, ember for a cell with two; a zero-area piece is its open edge (gold, the
!< wall on the cell's high edge, its centroid is undefined); the untrusted cells
!< are outlined in ember (2-D: `closed` false, no pieces; 3-D: `pieces = -1`).
program fig_ref_cut_solver_pieces
use fossil,     only : curve_object, cut_cells_2d_t
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P), parameter          :: H = 0.1_R8P, LIFT = 1.e-3_R8P, A = 0.5_R8P
type(curve_object)            :: body
type(cut_cells_2d_t)          :: cut
real(R8P)                     :: x(11), y(7), c(2), t(2), n(2), pu(4), pv(4)
type(vector_R8P), allocatable :: s(:), e(:), cen(:)
logical,          allocatable :: split(:), zero(:)
integer(I4P)                  :: i, j, k, m, p, w

x = [(H * i, i = 0, 10)] ; y = [(H * j, j = 0, 6)]
c = [0.66_R8P, 0.27_R8P] ; t = [cos(A), sin(A)] ; n = [-sin(A), cos(A)]
do k = 1, 4                                                   ! the plate corners, counter-clockwise
   pu(k) = c(1) + merge(0.25_R8P, -0.25_R8P, k == 2 .or. k == 3) * t(1) + merge(0.012_R8P, -0.012_R8P, k >= 3) * n(1)
   pv(k) = c(2) + merge(0.25_R8P, -0.25_R8P, k == 2 .or. k == 3) * t(2) + merge(0.012_R8P, -0.012_R8P, k >= 3) * n(2)
enddo
call body%adopt_loops(u=[0.13_R8P, 0.4_R8P, 0.4_R8P, 0.13_R8P, pu, 0.62_R8P, 0.78_R8P, 0.93_R8P], &
                      v=[0.13_R8P, 0.13_R8P, 0.4_R8P, 0.4_R8P, pv, 0.53_R8P, 0.47_R8P, 0.56_R8P], &
                      loop_sizes=[4_I4P, 4_I4P, 3_I4P], closed=[.true., .true., .false.])
call body%cut_cells(x=x, y=y, cut=cut)
allocate(split(cut%np), zero(cut%np), cen(cut%np))
do p = 1, cut%np
   split(p) = cut%pieces(cut%piece_cell(1, p), cut%piece_cell(2, p)) > 1
   zero(p)  = cut%piece_area(p) <= 0._R8P
   cen(p)   = q(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 2 * LIFT)
enddo
print '(A,I0,A,I0,A,I0)', 'cells ', size(cut%pieces), ', pieces ', cut%np, ', regular (whole) cells ', &
                          count(cut%alpha >= 1._R8P)
print '(A,I0,A,I0,A,I0)', 'cells with two pieces ', count(cut%pieces == 2), ', zero-area pieces ', count(zero), &
                          ', unclosed cells (no pieces) ', count(.not. cut%closed)

call fig_open(w, 'ref-cut-solver-pieces', 'one piece (teal), two pieces (ember), zero-area piece '// &
              '(gold edge), unclosed cell (ember)', -90._R8P, 90._R8P, width=720_I4P, height=460_I4P)
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), i = 1, 10)], j = 1, 7)], &
                  [([(q(x(i + 1), y(j), 0._R8P), i = 1, 10)], j = 1, 7)])
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), j = 1, 6)], i = 1, 11)], &
                  [([(q(x(i), y(j + 1), 0._R8P), j = 1, 6)], i = 1, 11)])
call fig_quad(w, 'body', q(0.13_R8P, 0.13_R8P, 0._R8P), q(0.4_R8P, 0.13_R8P, 0._R8P), q(0.4_R8P, 0.4_R8P, 0._R8P), &
              q(0.13_R8P, 0.4_R8P, 0._R8P))
call fig_quad(w, 'body', q(pu(1), pv(1), 0._R8P), q(pu(2), pv(2), 0._R8P), q(pu(3), pv(3), 0._R8P), &
              q(pu(4), pv(4), 0._R8P))
call fig_segments(w, 'body', [q(0.62_R8P, 0.53_R8P, LIFT), q(0.78_R8P, 0.47_R8P, LIFT)], &
                  [q(0.78_R8P, 0.47_R8P, LIFT), q(0.93_R8P, 0.56_R8P, LIFT)])
allocate(s(4 * size(cut%pieces)), e(4 * size(cut%pieces))) ; m = 0
do j = 1, 6                                                   ! untrusted cells
   do i = 1, 10
      if (cut%closed(i, j)) cycle
      s(m + 1) = q(x(i), y(j), LIFT)         ; e(m + 1) = q(x(i + 1), y(j), LIFT)
      s(m + 2) = q(x(i + 1), y(j), LIFT)     ; e(m + 2) = q(x(i + 1), y(j + 1), LIFT)
      s(m + 3) = q(x(i + 1), y(j + 1), LIFT) ; e(m + 3) = q(x(i), y(j + 1), LIFT)
      s(m + 4) = q(x(i), y(j + 1), LIFT)     ; e(m + 4) = q(x(i), y(j), LIFT)
      m = m + 4
   enddo
enddo
call fig_segments(w, 'hole', s(1:m), e(1:m))
call fig_points(w, 'probe', pack(cen, .not. split .and. .not. zero))
call fig_points(w, 'hole', pack(cen, split .and. .not. zero))
m = 0                                                         ! zero-area pieces: their open edges
do p = 1, cut%np
   if (.not. zero(p)) cycle
   i = cut%piece_cell(1, p) ; j = cut%piece_cell(2, p)
   if (cut%piece_aperture(2, p) > 0._R8P) then
      m = m + 1 ; s(m) = q(x(i + 1), y(j), 3 * LIFT) ; e(m) = q(x(i + 1), y(j + 1), 3 * LIFT)
   endif
   if (cut%piece_aperture(4, p) > 0._R8P) then
      m = m + 1 ; s(m) = q(x(i), y(j + 1), 3 * LIFT) ; e(m) = q(x(i + 1), y(j + 1), 3 * LIFT)
   endif
enddo
call fig_segments(w, 'accent', s(1:m), e(1:m))
call fig_points(w, 'foot', [(0.5_R8P * (s(k) + e(k)), k = 1, m)])
call fig_label(w, 'accent', q(0.36_R8P, 0.45_R8P, 3 * LIFT), 'walls on grid lines')
call fig_label(w, 'accent', q(0.78_R8P, 0.63_R8P, 3 * LIFT), 'open curve')
call fig_label(w, 'accent', q(0.88_R8P, 0.25_R8P, 3 * LIFT), 'thin plate')
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_cut_solver_pieces
