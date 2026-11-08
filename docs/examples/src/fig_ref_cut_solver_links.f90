!run figref-cut-solver-links fig_ref_cut_solver_links
!fig ref-cut-solver-links
!< Figure of the reference page "Cut cells: solver support": the piece graph of a thin plate (0.03 thick, tilted) across
!< 8 x 6 cells of 0.125. Every piece is a dot at its centroid (ember: a cell with two pieces, one on each side of the
!< plate); every link runs from a piece's centroid to the middle of the fluid edge interval it shares with the other
!< piece, then to that piece's centroid (gold). No link crosses the plate.
program fig_ref_cut_solver_links
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: A = 0.35_R8P ! plate angle [rad]
type(curve_object)            :: plate
type(cut_cells_2d_t)          :: cut
type(cut_links_t)             :: links
real(R8P)                     :: x(9), y(7), c(2), t(2), n(2)
type(vector_R8P), allocatable :: s(:), e(:)
logical,          allocatable :: split(:)
type(vector_R8P)              :: mid
integer(I4P)                  :: i, j, k, p, m, w

c = [0.5_R8P, 0.37_R8P] ; t = [cos(A), sin(A)] ; n = [-sin(A), cos(A)]
call plate%adopt_loops(u=[c(1) - 0.45_R8P * t(1) - 0.015_R8P * n(1), c(1) + 0.45_R8P * t(1) - 0.015_R8P * n(1), &
                          c(1) + 0.45_R8P * t(1) + 0.015_R8P * n(1), c(1) - 0.45_R8P * t(1) + 0.015_R8P * n(1)], &
                       v=[c(2) - 0.45_R8P * t(2) - 0.015_R8P * n(2), c(2) + 0.45_R8P * t(2) - 0.015_R8P * n(2), &
                          c(2) + 0.45_R8P * t(2) + 0.015_R8P * n(2), c(2) - 0.45_R8P * t(2) + 0.015_R8P * n(2)], &
                       loop_sizes=[4_I4P])
x = [(0.125_R8P * i, i = 0, 8)] ; y = [(0.125_R8P * i, i = 0, 6)]
call plate%cut_cells(x=x, y=y, cut=cut)
call cut_links_2d(cut, links)
allocate(split(cut%np))
split = [(cut%pieces(cut%piece_cell(1, p), cut%piece_cell(2, p)) > 1, p = 1, cut%np)]
print '(A,I0,A,I0,A,I0,A,I0)', 'cells ', size(cut%pieces), ', pieces ', cut%np, ', cells with two pieces ', &
                               count(cut%pieces == 2), ', links (both directions) ', size(links%piece)
call fig_open(w, 'ref-cut-solver-links', 'piece graph: pieces (dots, ember: two per cell) and links (gold)', &
              -90._R8P, 90._R8P, height=380_I4P)
call fig_segments(w, 'box', [([(q(x(i), y(j)), i = 1, 8)], j = 1, 7)], [([(q(x(i + 1), y(j)), i = 1, 8)], j = 1, 7)])
call fig_segments(w, 'box', [([(q(x(i), y(j)), j = 1, 6)], i = 1, 9)], [([(q(x(i), y(j + 1)), j = 1, 6)], i = 1, 9)])
call fig_quad(w, 'body', q(plate_u(1), plate_v(1)), q(plate_u(2), plate_v(2)), q(plate_u(3), plate_v(3)), &
              q(plate_u(4), plate_v(4)))
allocate(s(2 * size(links%piece)), e(2 * size(links%piece))) ; m = 0
do p = 1, cut%np ! each link drawn from a centroid to the middle of the shared fluid interval to the other centroid
   do k = links%start(p), links%start(p + 1) - 1
      if (links%piece(k) < p) cycle
      mid = shared(p, links%piece(k))
      m = m + 1 ; s(m) = q(cut%piece_centroid(1, p), cut%piece_centroid(2, p)) ; e(m) = mid
      m = m + 1 ; s(m) = mid ; e(m) = q(cut%piece_centroid(1, links%piece(k)), cut%piece_centroid(2, links%piece(k)))
   enddo
enddo
call fig_segments(w, 'accent', s(1:m), e(1:m))
call fig_points(w, 'probe', pack([(q(cut%piece_centroid(1, p), cut%piece_centroid(2, p)), p = 1, cut%np)], .not. split))
call fig_points(w, 'hole', pack([(q(cut%piece_centroid(1, p), cut%piece_centroid(2, p)), p = 1, cut%np)], split))
call fig_close(w)
contains
   function q(a, b) result(pt)
   !< A point of the plane.
   real(R8P), intent(in) :: a, b !< Coordinates.
   type(vector_R8P)      :: pt   !< Point.

   pt = vector_R8P(a, b, 0._R8P)
   endfunction q

   function shared(a, b) result(pt)
   !< The middle of the fluid edge interval that pieces a and b share.
   integer(I4P), intent(in) :: a, b !< Pieces.
   type(vector_R8P)         :: pt   !< Point.
   integer(I4P)             :: ia, ib, id, ii, jj

   pt = q(0._R8P, 0._R8P)
   do ia = cut%piece_iv_start(a), cut%piece_iv_start(a + 1) - 1
      do ib = cut%piece_iv_start(b), cut%piece_iv_start(b + 1) - 1
         if (cut%piece_iv_edge(ia) /= cut%piece_iv_edge(ib) .or. cut%piece_iv_lo(ia) /= cut%piece_iv_lo(ib)) cycle
         id = cut%piece_iv_edge(ia)
         if (id <= 9 * 6) then ! a vertical edge (ii, jj)
            ii = mod(id - 1, 9) + 1 ; jj = (id - 1) / 9 + 1
            pt = q(x(ii), 0.5_R8P * (cut%piece_iv_lo(ia) + cut%piece_iv_hi(ia)))
         else                  ! a horizontal one
            id = id - 9 * 6 ; ii = mod(id - 1, 8) + 1 ; jj = (id - 1) / 8 + 1
            pt = q(0.5_R8P * (cut%piece_iv_lo(ia) + cut%piece_iv_hi(ia)), y(jj))
         endif
         return
      enddo
   enddo
   endfunction shared

   real(R8P) function plate_u(k)
   !< u of the plate corner k.
   integer, intent(in) :: k !< Corner.

   plate_u = c(1) + merge(0.45_R8P, -0.45_R8P, k == 2 .or. k == 3) * t(1) + merge(0.015_R8P, -0.015_R8P, k >= 3) * n(1)
   endfunction plate_u

   real(R8P) function plate_v(k)
   !< v of the plate corner k.
   integer, intent(in) :: k !< Corner.

   plate_v = c(2) + merge(0.45_R8P, -0.45_R8P, k == 2 .or. k == 3) * t(2) + merge(0.015_R8P, -0.015_R8P, k >= 3) * n(2)
   endfunction plate_v
endprogram fig_ref_cut_solver_links
