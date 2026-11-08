!run figref-cut-solver-coarsen fig_ref_cut_solver_coarsen
!fig ref-cut-solver-coarsen
!< Figure of the reference page "Cut cells: solver support", section Coarsening: a thin plate (0.02 thick, tilted by
!< 0.45 rad) across 12 x 8 fine cells of 1/12, coarsened at ratio 2 (cut_coarsen_2d). The fine grid is faint, the
!< coarse grid gold; every fine piece is a small teal dot, the fine links inside one coarse cell are teal lines, and
!< every coarse piece, a connected component of those links, is an ember dot at its centroid: two per coarse cell
!< where the plate separates the fine pieces.
program fig_ref_cut_solver_coarsen
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, cut_coarsen_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NX = 12, NY = 8, R = 2
real(R8P),    parameter       :: H = 1._R8P / 12._R8P, A = 0.45_R8P, LIFT = 1.e-3_R8P
type(curve_object)            :: plate
type(cut_cells_2d_t)          :: fine, coarse
type(cut_links_t)             :: links, coarse_links
real(R8P)                     :: x(NX + 1), y(NY + 1), c(2), t(2), n(2), pu(4), pv(4)
integer(I4P),     allocatable :: f2c(:)
type(vector_R8P), allocatable :: s(:), e(:)
integer(I4P)                  :: i, j, k, m, w, p, q2, status

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
c = [0.5_R8P, 0.33_R8P] ; t = [cos(A), sin(A)] ; n = [-sin(A), cos(A)]
do k = 1, 4                                                   ! the plate corners, counter-clockwise
   pu(k) = c(1) + merge(0.5_R8P, -0.5_R8P, k == 2 .or. k == 3) * t(1) + merge(0.01_R8P, -0.01_R8P, k >= 3) * n(1)
   pv(k) = c(2) + merge(0.5_R8P, -0.5_R8P, k == 2 .or. k == 3) * t(2) + merge(0.01_R8P, -0.01_R8P, k >= 3) * n(2)
enddo
call plate%adopt_loops(u=pu, v=pv, loop_sizes=[4_I4P])
call plate%cut_cells(x=x, y=y, cut=fine)
call cut_links_2d(fine, links)
call cut_coarsen_2d(fine, links, x, y, R, coarse, coarse_links, f2c, status)
print '(A,I0,A,I0,A,I0,A,I0)', 'fine: pieces ', fine%np, ', cells with two pieces ', count(fine%pieces == 2), &
                               ', links ', size(links%piece), ', status ', status
print '(A,I0,A,I0,A,I0)', 'coarse (ratio 2): pieces ', coarse%np, ', cells with two pieces ', count(coarse%pieces == 2), &
                          ', links ', size(coarse_links%piece)
print '(A,ES9.2)', 'coarse fluid area - fine fluid area ', sum(coarse%piece_area) - sum(fine%piece_area)

call fig_open(w, 'ref-cut-solver-coarsen', 'ratio 2: fine pieces (teal) joined by the links inside a coarse cell; '// &
              'coarse pieces (ember)', -90._R8P, 90._R8P, width=800_I4P, height=560_I4P)
call fig_segments(w, 'box', [([(pt(x(i), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)], &
                  [([(pt(x(i + 1), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)])
call fig_segments(w, 'box', [([(pt(x(i), y(j), 0._R8P), j = 1, NY)], i = 1, NX + 1)], &
                  [([(pt(x(i), y(j + 1), 0._R8P), j = 1, NY)], i = 1, NX + 1)])
call fig_segments(w, 'accent', [([(pt(x(i), y(j), LIFT), i = 1, NX, R)], j = 1, NY + 1, R)], &
                  [([(pt(x(i + R), y(j), LIFT), i = 1, NX, R)], j = 1, NY + 1, R)])
call fig_segments(w, 'accent', [([(pt(x(i), y(j), LIFT), j = 1, NY, R)], i = 1, NX + 1, R)], &
                  [([(pt(x(i), y(j + R), LIFT), j = 1, NY, R)], i = 1, NX + 1, R)])
call fig_quad(w, 'body', pt(pu(1), pv(1), 0._R8P), pt(pu(2), pv(2), 0._R8P), pt(pu(3), pv(3), 0._R8P), &
              pt(pu(4), pv(4), 0._R8P))
allocate(s(size(links%piece)), e(size(links%piece))) ; m = 0
do p = 1, fine%np
   do k = links%start(p), links%start(p + 1) - 1
      q2 = links%piece(k)
      if (q2 < p) cycle
      if ((fine%piece_cell(1, p) - 1) / R /= (fine%piece_cell(1, q2) - 1) / R .or. &
          (fine%piece_cell(2, p) - 1) / R /= (fine%piece_cell(2, q2) - 1) / R) cycle
      m = m + 1
      s(m) = pt(fine%piece_centroid(1, p), fine%piece_centroid(2, p), 2 * LIFT)
      e(m) = pt(fine%piece_centroid(1, q2), fine%piece_centroid(2, q2), 2 * LIFT)
   enddo
enddo
call fig_segments(w, 'probe', s(1:m), e(1:m))
call fig_points(w, 'probe', [(pt(fine%piece_centroid(1, p), fine%piece_centroid(2, p), 2 * LIFT), p = 1, fine%np)])
call fig_points(w, 'hole', [(pt(coarse%piece_centroid(1, p), coarse%piece_centroid(2, p), 3 * LIFT), p = 1, coarse%np)])
call fig_close(w)
contains
   function pt(a_, b_, dz) result(rr)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: rr         !< Point.

   rr = vector_R8P(a_, b_, dz)
   endfunction pt
endprogram fig_ref_cut_solver_coarsen
