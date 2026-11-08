!< Figure of the reference page "Cut cells: solver support", section "Neighbourhoods for flux redistribution": a thin
!< plate (0.02 thick, tilted by 0.3 rad) across 10 x 8 cells of 0.1. The middle cell the plate crosses holds two
!< pieces; the flux-redistribution set of each (frd_neighbourhoods) is drawn as lines from the piece to its members:
!< teal for the piece above the plate, ember for the piece below. Neither set crosses the plate, though both lie in the
!< same 3 x 3 block of cells: a member must be reached by a monotone path of links.
program fig_ref_cut_solver_frd
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, frd_t, frd_neighbourhoods
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NX = 10, NY = 8
real(R8P),    parameter       :: H = 0.1_R8P, A = 0.3_R8P, LIFT = 1.e-3_R8P
type(curve_object)            :: plate
type(cut_cells_2d_t)          :: cut
type(cut_links_t)             :: links
type(frd_t)                   :: frd
real(R8P)                     :: x(NX + 1), y(NY + 1), c(2), t(2), n(2), pu(4), pv(4)
type(vector_R8P), allocatable :: s(:), e(:)
integer(I4P)                  :: i, j, k, m, w, ic, jc, p(2)
character(8)                  :: side(2) = ['above   ', 'below   ']
character(5)                  :: style(2) = ['probe', 'hole ']

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
c = [0.5_R8P, 0.4_R8P] ; t = [cos(A), sin(A)] ; n = [-sin(A), cos(A)]
do k = 1, 4                                                   ! the plate corners, counter-clockwise
   pu(k) = c(1) + merge(0.42_R8P, -0.42_R8P, k == 2 .or. k == 3) * t(1) + merge(0.01_R8P, -0.01_R8P, k >= 3) * n(1)
   pv(k) = c(2) + merge(0.42_R8P, -0.42_R8P, k == 2 .or. k == 3) * t(2) + merge(0.01_R8P, -0.01_R8P, k >= 3) * n(2)
enddo
call plate%adopt_loops(u=pu, v=pv, loop_sizes=[4_I4P])
call plate%cut_cells(x=x, y=y, cut=cut)
call cut_links_2d(cut, links)
call frd_neighbourhoods(links, cut%piece_cell, cut%piece_area, x, y, frd)
ic = 6 ; jc = 5                                               ! the cell of (0.55, 0.42)
if (cut%pieces(ic, jc) /= 2) error stop 'the chosen cell does not hold two pieces'
p(1) = cut%first_piece(ic, jc) ; p(2) = p(1) + 1
if (cut%piece_centroid(2, p(1)) < cut%piece_centroid(2, p(2))) p = p([2, 1])   ! p(1) above the plate
print '(A,2I3,A,I0,A,I0)', 'cell', ic, jc, ': pieces ', cut%pieces(ic, jc), ', cells with two pieces ', &
                           count(cut%pieces == 2)
do k = 1, 2
   print '(A,A,A,F6.3,A,I0,A,F6.3,A,ES9.2)', '  piece ', trim(side(k)), ': fraction ', cut%piece_area(p(k)) / H**2, &
         ', FRD members ', frd%start(p(k) + 1) - frd%start(p(k)), ', own weight ', frd%weight(frd%start(p(k))), &
         ', weights sum - 1 ', sum(frd%weight(frd%start(p(k)):frd%start(p(k) + 1) - 1)) - 1._R8P
enddo

call fig_open(w, 'ref-cut-solver-frd', 'FRD sets of the two pieces of one cell: above the plate (teal), below (ember)', &
              -90._R8P, 90._R8P, width=640_I4P, height=500_I4P)
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)], &
                  [([(q(x(i + 1), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)])
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), j = 1, NY)], i = 1, NX + 1)], &
                  [([(q(x(i), y(j + 1), 0._R8P), j = 1, NY)], i = 1, NX + 1)])
call fig_segments(w, 'accent', [q(x(ic - 1), y(jc - 1), LIFT), q(x(ic + 2), y(jc - 1), LIFT), &
                                q(x(ic + 2), y(jc + 2), LIFT), q(x(ic - 1), y(jc + 2), LIFT)], &
                               [q(x(ic + 2), y(jc - 1), LIFT), q(x(ic + 2), y(jc + 2), LIFT), &
                                q(x(ic - 1), y(jc + 2), LIFT), q(x(ic - 1), y(jc - 1), LIFT)])
call fig_quad(w, 'body', q(pu(1), pv(1), 0._R8P), q(pu(2), pv(2), 0._R8P), q(pu(3), pv(3), 0._R8P), &
              q(pu(4), pv(4), 0._R8P))
allocate(s(9), e(9))
do k = 1, 2
   m = 0
   do i = frd%start(p(k)) + 1, frd%start(p(k) + 1) - 1
      m = m + 1
      s(m) = q(cut%piece_centroid(1, p(k)), cut%piece_centroid(2, p(k)), 2 * LIFT)
      e(m) = q(cut%piece_centroid(1, frd%member(i)), cut%piece_centroid(2, frd%member(i)), 2 * LIFT)
   enddo
   call fig_segments(w, trim(style(k)), s(1:m), e(1:m))
   call fig_points(w, trim(style(k)), [s(1), e(1:m)])
enddo
call fig_label(w, 'accent', q(x(ic + 2) + 0.01_R8P, y(jc + 2) + 0.02_R8P, 3 * LIFT), '3 x 3 block')
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_cut_solver_frd
