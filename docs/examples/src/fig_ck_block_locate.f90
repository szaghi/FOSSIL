!run figck-block-locate fig_ck_block_locate
!fig ck-block-locate
!< Figure of the cookbook recipe "Locate a point in a curvilinear block": the O-grid of ck_block_locate, the cell that
!< contains the query point next to the seam, and the 16 vertices of its width-4 stencil, on both sides of the seam.
program fig_ck_block_locate
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d
use fossil_fig, only : fig_open, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter         :: PI = acos(-1._R8P), lift = 1.e-3_R8P
type(curvilinear_block_2d_t) :: grid(1)
type(vector_R8P)             :: a(32 * 9 + 33 * 8), b(32 * 9 + 33 * 8), s(16)
real(R8P)                    :: u(33, 9), v(33, 9), r
real(R8P),    allocatable    :: xi(:,:), w(:,:)
integer(I4P), allocatable    :: cell(:,:), st(:), nodes(:,:,:), sst(:)
integer(I4P)                 :: i, j, k, m, uu

do j = 1, 9
   r = 1._R8P + 0.15_R8P * (j - 1)
   do i = 1, 33
      u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call grid(1)%initialize(u, v) ; call detect_topology_2d(grid) ; call grid(1)%build_search
call grid(1)%locate_many([1.62_R8P], [0.03_R8P], cell, xi, st)
call grid(1)%stencils(cell, xi, 4_I4P, nodes, w, sst)

call fig_open(uu, 'ck-block-locate', 'an O-grid, the cell holding a point next to its seam, and the width-4 stencil', &
              -90._R8P, 90._R8P, height=520_I4P)
i = cell(1, 1) ; j = cell(2, 1)
call fig_quad(uu, 'cell', pt(u(i, j), v(i, j)), pt(u(i + 1, j), v(i + 1, j)), pt(u(i + 1, j + 1), v(i + 1, j + 1)), &
              pt(u(i, j + 1), v(i, j + 1)))
m = 0
do j = 1, 9
   do i = 1, 32
      m = m + 1 ; a(m) = pt(u(i, j), v(i, j), lift) ; b(m) = pt(u(i + 1, j), v(i + 1, j), lift)
   enddo
enddo
do j = 1, 8
   do i = 1, 33
      m = m + 1 ; a(m) = pt(u(i, j), v(i, j), lift) ; b(m) = pt(u(i, j + 1), v(i, j + 1), lift)
   enddo
enddo
call fig_segments(uu, 'body', a, b)
do k = 1, 16
   s(k) = pt(u(nodes(1, k, 1), nodes(2, k, 1)), v(nodes(1, k, 1), nodes(2, k, 1)), 2 * lift)
enddo
call fig_points(uu, 'foot', s)
call fig_points(uu, 'probe', [pt(1.62_R8P, 0.03_R8P, 3 * lift)])
call fig_label(uu, 'accent', pt(1.66_R8P, 0.08_R8P, 3 * lift), 'p')
call fig_label(uu, 'accent', pt(2.25_R8P, -0.05_R8P, 3 * lift), 'seam')
call fig_close(uu)
contains
   function pt(px, py, dz) result(q)
   !< A point of the plane z = 0, lifted by dz (optional) towards the viewer.
   real(R8P), intent(in)           :: px, py !< In-plane coordinates.
   real(R8P), intent(in), optional :: dz     !< Lift.
   type(vector_R8P)                :: q      !< Point.

   q = px * ex_R8P + py * ey_R8P
   if (present(dz)) q = q + dz * ez_R8P
   endfunction pt
endprogram fig_ck_block_locate
