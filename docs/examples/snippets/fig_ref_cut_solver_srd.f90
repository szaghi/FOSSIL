!< Figure of the reference page "Cut cells: solver support", section "Neighbourhoods for state redistribution": a wavy
!< wall (solid below y = 0.24 + 0.08 sin(3 pi x)) on 12 x 6 cells of 1/12, the cells shaded by their fluid fraction.
!< Every small piece (area below half the cell, teal) has a short ember stroke along its unit wall normal and gold
!< lines to the other members of its neighbourhood (srd_neighbourhoods, tau = 1/2): they lie along the normal.
program fig_ref_cut_solver_srd
use fossil,     only : curve_object, cut_cells_2d_t, cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, &
                       SRD_COMPLETE, SRD_INCOMPLETE
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NX = 12, NY = 6, NW = 96
real(R8P),    parameter       :: H = 1._R8P / 12._R8P, PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)            :: wall
type(cut_cells_2d_t)          :: cut
type(cut_links_t)             :: links
type(srd_t)                   :: srd
real(R8P)                     :: x(NX + 1), y(NY + 1), u(NW + 3), v(NW + 3), nrm(2), sz
type(vector_R8P), allocatable :: s(:), e(:), ps(:), pn(:)
integer(I4P)                  :: i, j, k, m, p, w, nsmall, mmax

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
u(1) = -0.05_R8P ; v(1) = -0.1_R8P ; u(2) = 1.05_R8P ; v(2) = -0.1_R8P   ! counter-clockwise, the solid on the left:
do k = 0, NW                                                  ! under the block, then the wave from right to left
   u(k + 3) = 1.05_R8P - 1.1_R8P * k / NW ; v(k + 3) = 0.24_R8P + 0.08_R8P * sin(3 * PI * u(k + 3))
enddo
call wall%adopt_loops(u=u, v=v, loop_sizes=[NW + 3_I4P])
call wall%cut_cells(x=x, y=y, cut=cut)
call cut_links_2d(cut, links)
call srd_neighbourhoods(links, cut%piece_cell, cut%piece_area, cut%piece_centroid, cut%piece_wet_normal, x, y, srd)
nsmall = count(srd%status == SRD_COMPLETE .or. srd%status == SRD_INCOMPLETE)
mmax = 0
do p = 1, cut%np
   if (srd%status(p) == SRD_COMPLETE) mmax = max(mmax, srd%start(p + 1) - srd%start(p))
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'pieces ', cut%np, ', small ', nsmall, ', incomplete ', count(srd%status == SRD_INCOMPLETE), &
                                ', largest neighbourhood ', mmax
print '(A,ES9.2)', 'sum of the weighted volumes - fluid area ', sum(srd%vhat) - sum(cut%piece_area)
print '(A,F6.3,A,F6.3)', 'eta_nb of the small pieces from ', minval(srd%eta_nb, srd%status == SRD_COMPLETE), ' to ', &
                         maxval(srd%eta_nb, srd%status == SRD_COMPLETE)

call fig_open(w, 'ref-cut-solver-srd', 'small pieces (teal), their wall normals (ember) and SRD members (gold lines)', &
              -90._R8P, 90._R8P, width=720_I4P, height=440_I4P)
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)], &
                  [([(q(x(i + 1), y(j), 0._R8P), i = 1, NX)], j = 1, NY + 1)])
call fig_segments(w, 'box', [([(q(x(i), y(j), 0._R8P), j = 1, NY)], i = 1, NX + 1)], &
                  [([(q(x(i), y(j + 1), 0._R8P), j = 1, NY)], i = 1, NX + 1)])
do j = 1, NY
   do i = 1, NX
      if (cut%alpha(i, j) <= 0._R8P) cycle
      call fig_quad(w, 'frac', q(x(i), y(j), 0._R8P), q(x(i + 1), y(j), 0._R8P), q(x(i + 1), y(j + 1), 0._R8P), &
                    q(x(i), y(j + 1), 0._R8P), value=cut%alpha(i, j))
   enddo
enddo
call fig_segments(w, 'body', [(q(u(k), v(k), LIFT), k = 3, NW + 2)], [(q(u(k + 1), v(k + 1), LIFT), k = 3, NW + 2)])
allocate(s(16 * cut%np), e(16 * cut%np), ps(cut%np), pn(cut%np)) ; m = 0 ; k = 0
do p = 1, cut%np
   if (srd%status(p) /= SRD_COMPLETE .and. srd%status(p) /= SRD_INCOMPLETE) cycle
   do i = srd%start(p) + 1, srd%start(p + 1) - 1
      m = m + 1
      s(m) = q(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 2 * LIFT)
      e(m) = q(cut%piece_centroid(1, srd%member(i)), cut%piece_centroid(2, srd%member(i)), 2 * LIFT)
   enddo
   k = k + 1 ; ps(k) = q(cut%piece_centroid(1, p), cut%piece_centroid(2, p), 3 * LIFT)
   sz = sqrt(sum(cut%piece_wet_normal(:, p)**2)) ; nrm = cut%piece_wet_normal(:, p) / sz    ! into the fluid
   pn(k) = q(cut%piece_centroid(1, p) + 0.35_R8P * H * nrm(1), cut%piece_centroid(2, p) + 0.35_R8P * H * nrm(2), 3 * LIFT)
enddo
call fig_segments(w, 'accent', s(1:m), e(1:m))
call fig_points(w, 'foot', e(1:m))
call fig_segments(w, 'hole', ps(1:k), pn(1:k))
call fig_points(w, 'probe', ps(1:k))
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_cut_solver_srd
