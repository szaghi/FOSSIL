!run figref-cut-solver-merge fig_ref_cut_solver_merge
!fig ref-cut-solver-merge
!< Figure of the reference page "Cut cells: solver support", section "Merge targets": a solid disc (radius 0.27, a
!< 96-gon) on 10 x 10 cells of 0.1, the cells shaded by their fluid fraction. Every small piece (teal) has a gold line
!< to its single target (merge_targets, tau = 1/2), the linked piece that is not small and best aligned with the wall
!< normal. One diagonal piece shows its normal-directed targets in ember, each labelled with its weight.
program fig_ref_cut_solver_merge
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, merge_t, merge_targets, MERGE_NONE, &
                       MERGE_OK, MERGE_CHAIN
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: N = 10, NS = 96
real(R8P),    parameter       :: H = 0.1_R8P, PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)            :: disc
type(cut_cells_2d_t)          :: cut
type(cut_links_t)             :: links
type(merge_t)                 :: mt
real(R8P)                     :: x(N + 1), u(NS), v(NS)
type(vector_R8P), allocatable :: s(:), e(:)
integer(I4P)                  :: i, j, k, m, w, p, pm
character(16)                 :: lab

x = [(H * i, i = 0, N)]
u = [(0.52_R8P + 0.27_R8P * cos(2 * PI * (k - 1) / NS), k = 1, NS)]
v = [(0.47_R8P + 0.27_R8P * sin(2 * PI * (k - 1) / NS), k = 1, NS)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[NS])
call disc%cut_cells(x=x, y=x, cut=cut)
call cut_links_2d(cut, links)
call merge_targets(links, cut%piece_cell, cut%piece_area, cut%piece_wet_normal, x, x, mt)
print '(A,I0,A,I0,A,I0,A,I0)', 'pieces ', cut%np, ', served ', count(mt%status /= MERGE_NONE), ', direct ', &
                               count(mt%status == MERGE_OK), ', by a chain ', count(mt%status == MERGE_CHAIN)
pm = 0                                                        ! a served piece with three normal-directed targets
do p = 1, cut%np
   if (mt%status(p) == MERGE_NONE) cycle
   if (mt%mix_start(p + 1) - mt%mix_start(p) == 3 .and. cut%piece_centroid(1, p) > 0.6_R8P .and. &
       cut%piece_centroid(2, p) > 0.6_R8P) then
      pm = p ; exit
   endif
enddo
if (pm == 0) error stop 'no diagonal piece'
print '(A,2I3,A,F6.3,A,3F6.3)', 'piece in cell', cut%piece_cell(:, pm), ', fraction ', cut%piece_area(pm) / H**2, &
                                ', normal-directed weights', mt%mix_weight(mt%mix_start(pm):mt%mix_start(pm + 1) - 1)

call fig_open(w, 'ref-cut-solver-merge', 'merge targets of the small pieces (gold); '// &
              'normal-directed targets and weights of one (ember)', -90._R8P, 90._R8P, width=800_I4P, height=560_I4P)
call fig_segments(w, 'box', [([(q(x(i), x(j), 0._R8P), i = 1, N)], j = 1, N + 1)], &
                  [([(q(x(i + 1), x(j), 0._R8P), i = 1, N)], j = 1, N + 1)])
call fig_segments(w, 'box', [([(q(x(i), x(j), 0._R8P), j = 1, N)], i = 1, N + 1)], &
                  [([(q(x(i), x(j + 1), 0._R8P), j = 1, N)], i = 1, N + 1)])
do j = 1, N
   do i = 1, N
      if (cut%alpha(i, j) <= 0._R8P) cycle
      call fig_quad(w, 'frac', q(x(i), x(j), 0._R8P), q(x(i + 1), x(j), 0._R8P), q(x(i + 1), x(j + 1), 0._R8P), &
                    q(x(i), x(j + 1), 0._R8P), value=cut%alpha(i, j))
   enddo
enddo
call fig_segments(w, 'body', [(q(u(k), v(k), LIFT), k = 1, NS)], [(q(u(mod(k, NS) + 1), v(mod(k, NS) + 1), LIFT), k = 1, NS)])
allocate(s(cut%np), e(cut%np)) ; m = 0
do p = 1, cut%np
   if (mt%status(p) == MERGE_NONE .or. mt%target(p) < 1 .or. p == pm) cycle
   m = m + 1
   s(m) = c(p, 2 * LIFT) ; e(m) = c(mt%target(p), 2 * LIFT)
enddo
call fig_segments(w, 'accent', s(1:m), e(1:m))
call fig_points(w, 'foot', e(1:m))
call fig_points(w, 'probe', s(1:m))
m = 0
do k = mt%mix_start(pm), mt%mix_start(pm + 1) - 1
   m = m + 1 ; s(m) = c(pm, 3 * LIFT) ; e(m) = c(mt%mix_piece(k), 3 * LIFT)
   write(lab, '(F5.3)') mt%mix_weight(k)
   call fig_label(w, 'hole', e(m) + 0.6_R8P * (e(m) - s(m)), trim(lab))
enddo
call fig_segments(w, 'hole', s(1:m), e(1:m))
call fig_points(w, 'hole', [s(1), e(1:m)])
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q

   function c(pp, dz) result(r)
   !< The centroid of piece pp, lifted by dz.
   integer(I4P), intent(in) :: pp !< Piece.
   real(R8P),    intent(in) :: dz !< Lift.
   type(vector_R8P)         :: r  !< Point.

   r = vector_R8P(cut%piece_centroid(1, pp), cut%piece_centroid(2, pp), dz)
   endfunction c
endprogram fig_ref_cut_solver_merge
