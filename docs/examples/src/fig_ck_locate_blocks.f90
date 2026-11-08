!run figck-locate-blocks fig_ck_locate_blocks
!fig ck-locate-blocks
!< Figure of the cookbook recipe "Find which block holds a point": the two blocks of ck_locate_blocks seen from above,
!< the shell's equator (its inner and outer circles and radial lines, stone) and the box (its outline at z = 0,
!< faint); the four points (drawn at their x, y), each labelled with the block locate_in_blocks returns (0: none).
program fig_ck_locate_blocks
use fossil,     only : curvilinear_block_t, locate_in_blocks
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
type(curvilinear_block_t) :: blocks(2)
integer(I4P), allocatable :: block(:), cell(:,:), status(:)
real(R8P),    allocatable :: xi(:,:)
real(R8P)                 :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph, g(9), px(4), py(4), pz(4)
real(R8P)                 :: bx(9, 9, 9), by(9, 9, 9), bz(9, 9, 9)
type(vector_R8P)          :: a(72), b(72), p(4)
integer(I4P)              :: i, j, k, u, n
character(12)             :: text

do k = 1, 3 ; do j = 1, 13 ; do i = 1, 25
   r = 1._R8P + 0.25_R8P * (k - 1) ; th = 2 * PI * (i - 1) / 24 ; ph = PI * (j - 1) / 12
   x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
enddo ; enddo ; enddo
x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; x(:,[1,13],:) = 0._R8P ; y(:,[1,13],:) = 0._R8P
call blocks(1)%initialize(x, y, z)
g = [(-1.6_R8P + 0.4_R8P * (i - 1), i = 1, 9)]
do k = 1, 9 ; do j = 1, 9
   bx(:, j, k) = g ; by(:, j, k) = g(j) ; bz(:, j, k) = g(k)
enddo ; enddo
call blocks(2)%initialize(bx, by, bz)
call blocks(1)%build_search ; call blocks(2)%build_search
px = [0.3_R8P, 1.2_R8P, 1.5_R8P, 2._R8P] ; py = [0.2_R8P, 0.3_R8P, 1.5_R8P, 2._R8P] ; pz = [0.1_R8P, 0.2_R8P, 0._R8P, 2._R8P]
call locate_in_blocks(blocks, px, py, pz, block=block, cell=cell, xi=xi, status=status)
n = 0
do i = 1, 24                                                ! the equator j = 7: circles r = 1, 1.5 and radial lines
   n = n + 1 ; a(n) = vector_R8P(x(i, 7, 1), y(i, 7, 1), 0._R8P) ; b(n) = vector_R8P(x(i + 1, 7, 1), y(i + 1, 7, 1), 0._R8P)
   n = n + 1 ; a(n) = vector_R8P(x(i, 7, 3), y(i, 7, 3), 0._R8P) ; b(n) = vector_R8P(x(i + 1, 7, 3), y(i + 1, 7, 3), 0._R8P)
   n = n + 1 ; a(n) = vector_R8P(x(i, 7, 1), y(i, 7, 1), 0._R8P) ; b(n) = vector_R8P(x(i, 7, 3), y(i, 7, 3), 0._R8P)
enddo
call fig_open(u, 'ck-locate-blocks', 'locate_in_blocks: the block of each point (shell 1, box 2, none 0)', &
              -90._R8P, 90._R8P, height=440_I4P)
call fig_box(u, 'ghost', vector_R8P(-1.6_R8P, -1.6_R8P, 0._R8P), vector_R8P(1.6_R8P, 1.6_R8P, 0._R8P))
call fig_segments(u, 'body', a, b)
do i = 1, 4
   p(i) = vector_R8P(px(i), py(i), 0._R8P)
   write(text, '(A,I0)') 'block ', block(i)
   call fig_label(u, 'foot', p(i) + vector_R8P(0.06_R8P, -0.05_R8P, 0._R8P), trim(text))
enddo
call fig_points(u, 'probe', p)
call fig_close(u)
endprogram fig_ck_locate_blocks
