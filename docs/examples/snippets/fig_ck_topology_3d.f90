!< Figure of the cookbook recipe "Find the seam and the poles of a 3-D block": the spherical shell of ck_topology_3d,
!< its inner sphere r = 1 in stone and its outer sphere r = 1.5 translucent; the seam faces 1 and 2 (i = 1 and
!< i = 25, the same vertices) in ember, the two collapsed pole lines (faces 3 and 4) in yellow.
program fig_ck_topology_3d
use fossil_fig, only : fig_open, fig_quad, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P)
real(R8P)            :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph
integer(I4P)         :: i, j, k, u

do k = 1, 3 ; do j = 1, 13 ; do i = 1, 25
   r = 1._R8P + 0.25_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:)
x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call fig_open(u, 'ck-topology-3d', 'detect_topology: the seam i = 1 / i = 25 (ember) and the two poles (yellow)', &
              -30._R8P, 20._R8P)
do j = 1, 12 ; do i = 1, 24
   call fig_quad(u, 'body', p(i, j, 1), p(i + 1, j, 1), p(i + 1, j + 1, 1), p(i, j + 1, 1))
   call fig_quad(u, 'ghost', p(i, j, 3), p(i + 1, j, 3), p(i + 1, j + 1, 3), p(i, j + 1, 3))
enddo ; enddo
do k = 1, 2 ; do j = 1, 12
   call fig_quad(u, 'patch', p(1, j, k), p(1, j + 1, k), p(1, j + 1, k + 1), p(1, j, k + 1))
enddo ; enddo
call fig_segment(u, 'accent', p(1, 1, 1), p(1, 1, 3))
call fig_segment(u, 'accent', p(1, 13, 1), p(1, 13, 3))
call fig_label(u, 'accent', p(1, 1, 3) + vector_R8P(0.05_R8P, 0._R8P, 0.05_R8P), 'pole (face 3)')
call fig_label(u, 'accent', p(1, 13, 3) + vector_R8P(0.05_R8P, 0._R8P, -0.15_R8P), 'pole (face 4)')
call fig_label(u, 'patch', p(1, 5, 3) + vector_R8P(0.08_R8P, 0._R8P, 0._R8P), 'seam (faces 1, 2)')
call fig_close(u)
contains
   function p(i, j, k) result(v)
   integer(I4P), intent(in) :: i, j, k
   type(vector_R8P)         :: v

   v = vector_R8P(x(i, j, k), y(i, j, k), z(i, j, k))
   endfunction p
endprogram fig_ck_topology_3d
