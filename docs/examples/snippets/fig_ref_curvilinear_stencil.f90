!< Figure of the reference page "Curvilinear structured blocks": the 4 x 4 stencil of a point just above the wake cut
!< of a C-grid (the grid of the topology figure). The point lies in a cell of the upper wake, j = 1; the stencil's
!< first row, j = 0, is beyond face 3 and is mapped through the wake-cut record to the lower wake (ember), so the
!< stencil is the one of the grid unwrapped across the cut. Zoomed on the wake.
program fig_ref_curvilinear_stencil
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, LOC_INSIDE
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t) :: c(1)
real(R8P)                    :: u(61, 7), v(61, 7), w0(2, 61), nrm(2, 61), phi, f
real(R8P),    allocatable    :: xi(:,:), wt(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:), nodes(:,:,:), sst(:)
type(vector_R8P)             :: s(500), e(500)
integer(I4P)                 :: i, j, k, m, w
logical                      :: below(16)

do i = 1, 61 ! the C-grid: 13 vertices per wake side, 36 around the ellipse (chord 1, thickness 0.2), 7 layers
   if (i <= 13) then
      w0(:, i) = [3._R8P - 2._R8P * (i - 1) / 12, 0._R8P] ; nrm(:, i) = [0._R8P, -1._R8P]
   elseif (i >= 49) then
      w0(:, i) = [3._R8P - 2._R8P * (61 - i) / 12, 0._R8P] ; nrm(:, i) = [0._R8P, 1._R8P]
   else
      phi = 2 * PI * (i - 13) / 36
      w0(:, i) = [0.5_R8P + 0.5_R8P * cos(phi), -0.1_R8P * sin(phi)]
      nrm(:, i) = [0.1_R8P * cos(phi), -0.5_R8P * sin(phi)] ; nrm(:, i) = nrm(:, i) / norm2(nrm(:, i))
   endif
enddo
do k = 1, 40 ! smooth the marching directions along i, so that the layers do not cross at the trailing edge
   nrm(:, 2:60) = (nrm(:, 1:59) + 2 * nrm(:, 2:60) + nrm(:, 3:61)) / 4
   do i = 2, 60
      nrm(:, i) = nrm(:, i) / norm2(nrm(:, i))
   enddo
enddo
do i = 1, 61
   do j = 1, 7
      u(i, j) = w0(1, i) + 0.15_R8P * (j - 1) * nrm(1, i) ; v(i, j) = w0(2, i) + 0.15_R8P * (j - 1) * nrm(2, i)
   enddo
enddo
call c(1)%initialize(u, v) ; call detect_topology_2d(c) ; call c(1)%build_search
call c(1)%locate_many([2.05_R8P], [0.04_R8P], cell, xi, st)
print '(A,L1,A,2I3,A,2F8.4)', 'inside ', st(1) == LOC_INSIDE, ', cell', cell(:, 1), ', xi', xi(:, 1)
call c(1)%stencils(cell, xi, 4_I4P, nodes, wt, sst)
do k = 1, 16
   below(k) = v(nodes(1, k, 1), nodes(2, k, 1)) < 0._R8P .or. &
              (v(nodes(1, k, 1), nodes(2, k, 1)) == 0._R8P .and. nodes(1, k, 1) <= 13)
enddo
print '(A,4(2I3,2X))', 'stencil nodes (i, j), first row ', (nodes(:, k, 1), k = 1, 4)
print '(A,I0,A,F6.3)', 'nodes mapped below the cut ', count(below), ', weights sum ', sum(wt(:, 1))
do m = 2, 4, 2 ! interpolate exp(x) cos(3 y) with the widths 2 and 4
   call c(1)%stencils(cell, xi, int(m, I4P), nodes, wt, sst)
   f = 0._R8P
   do k = 1, m * m
      f = f + wt(k, 1) * exp(u(nodes(1, k, 1), nodes(2, k, 1))) * cos(3 * v(nodes(1, k, 1), nodes(2, k, 1)))
   enddo
   print '(A,I0,A,ES10.2)', 'width ', m, ': error of exp(x) cos(3 y) ', abs(f - exp(2.05_R8P) * cos(0.12_R8P))
enddo
call fig_open(w, 'ref-curvilinear-stencil', 'a 4 x 4 stencil across the wake cut: mapped nodes in ember', &
              -90._R8P, 90._R8P, height=380_I4P)
call fig_zoom(w, vector_R8P(2.05_R8P, 0.05_R8P, 0._R8P), 0.42_R8P)
m = 0
do j = 1, 7
   do i = 1, 60
      m = m + 1 ; s(m) = p(i, j) ; e(m) = p(i + 1, j)
   enddo
enddo
call fig_segments(w, 'box', s(1:m), e(1:m))
m = 0
do j = 1, 6
   do i = 1, 61
      if (i > 13 .and. i < 49) cycle
      m = m + 1 ; s(m) = p(i, j) ; e(m) = p(i, j + 1)
   enddo
enddo
call fig_segments(w, 'box', s(1:m), e(1:m))
call fig_segments(w, 'accent', [(p(i, 1_I4P), i = 49, 60)], [(p(i + 1, 1_I4P), i = 49, 60)])
call fig_quad(w, 'ghost', p(cell(1, 1), cell(2, 1)), p(cell(1, 1) + 1, cell(2, 1)), p(cell(1, 1) + 1, cell(2, 1) + 1), &
              p(cell(1, 1), cell(2, 1) + 1))
call fig_points(w, 'probe', [(p(nodes(1, k, 1), nodes(2, k, 1)), k = 1, 16)])
call fig_points(w, 'hole', pack([(p(nodes(1, k, 1), nodes(2, k, 1)), k = 1, 16)], below))
call fig_points(w, 'foot', [vector_R8P(2.05_R8P, 0.04_R8P, 0._R8P)])
call fig_label(w, 'accent', vector_R8P(1.72_R8P, -0.02_R8P, 0._R8P), 'wake cut')
call fig_close(w)
contains
   function p(i, j) result(q)
   !< Vertex (i, j).
   integer(I4P), intent(in) :: i, j !< Indices.
   type(vector_R8P)         :: q    !< Point.

   q = vector_R8P(u(i, j), v(i, j), 0._R8P)
   endfunction p
endprogram fig_ref_curvilinear_stencil
