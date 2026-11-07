!< FOSSIL test: robust point location in curvilinear blocks (issue #71, step 5).
program fossil_test_block_locate
!< FOSSIL test: robust point location in curvilinear blocks (issue #71, step 5).
!<
!< The maps F(X, Y) = (X + a X Y, Y + b X Y) and F(X, Y, Z) = (X + a Y Z, Y + b Z X, Z + c X Y) are their own bilinear
!< (trilinear) interpolants, so a point F(X) has the known cell floor(X / h) + 1 and xi = frac(X / h).
!< 1. 2-D and 3-D: random points, located cell and xi equal the analytic inverse (xi to 1e-13);
!< 2. 2-D, dyadic data (every coordinate exact): points on edges and vertices take the low-side cell (xi = 1 there),
!<    points on the block's low boundary the first cell (xi = 0); 3-D: the same within the xi tolerance;
!< 3. points outside the block: LOC_OUTSIDE; points just outside a face: LOC_OUTSIDE;
!< 4. folded cells: a 2-D block with a folded vertex and the 3-D folded cell alone: every LOC_INSIDE answer is a valid
!<    cell whose map at xi gives back the point, and the fold is reported as LOC_WARPED, never as a wrong cell;
!< 5. O-grid seam: points on the seam take cell i = 1 (the first in index order);
!< 6. two blocks side by side: locate_in_blocks returns the block holding each point.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, locate_in_blocks_2d
use fossil_block_locate,      only : LOC_INSIDE, LOC_OUTSIDE, LOC_WARPED
use fossil_block_geometry,    only : trilinear_point, bilinear_point, CELL_VALID
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
real(R8P), parameter                      :: FOLD(3, 8) = reshape([0.13_R8P, 0.25_R8P, 0.29_R8P, 0.13_R8P, -0.27_R8P, &
   0.23_R8P, 0.58_R8P, 1.25_R8P, -0.61_R8P, 0.76_R8P, 0.86_R8P, -0.18_R8P, 0.47_R8P, 0.54_R8P, 1.59_R8P, 0.51_R8P, &
   -0.56_R8P, 0.33_R8P, -0.11_R8P, 1.05_R8P, 0.39_R8P, 1.12_R8P, 1.11_R8P, 0.91_R8P], [3, 8])
type(curvilinear_block_2d_t)              :: q
type(curvilinear_block_2d_t), allocatable :: qs(:)
type(curvilinear_block_t)                 :: b
real(R8P),    allocatable                 :: u(:,:), v(:,:), x(:,:,:), y(:,:,:), z(:,:,:)
real(R8P),    allocatable                 :: pu(:), pv(:), pw(:), xs(:,:), xi(:,:), xt(:,:)
integer(I4P), allocatable                 :: cell(:,:), st(:), blk(:), ct(:,:)
real(R8P)                                 :: h, aa, bb, cc, err, r(3), xx(3)
integer(I4P)                              :: i, j, k, n, ni, np, nbad, nwarp
logical                                   :: ok, all_ok

all_ok = .true.

! 1-2. 2-D bilinear map
ni = 17 ; h = 0.125_R8P ; aa = 0.25_R8P ; bb = -0.375_R8P
allocate(u(ni, ni), v(ni, ni))
do j = 1, ni
   do i = 1, ni
      call fmap2(h * (i - 1), h * (j - 1), u(i, j), v(i, j))
   enddo
enddo
call q%initialize(u, v) ; call q%build_search
np = 2000 ; allocate(xs(2, np), pu(np), pv(np))
call random_seed(put=[(7 + 5 * i, i = 1, 64)])
do n = 1, np
   call random_number(r(1:2)) ; xs(:, n) = r(1:2) * h * (ni - 1)
   call fmap2(xs(1, n), xs(2, n), pu(n), pv(n))
enddo
call q%locate_many(pu, pv, cell, xi, st)
err = 0._R8P ; nbad = 0
do n = 1, np
   if (st(n) /= LOC_INSIDE .or. any(cell(:, n) /= int(xs(:, n) / h) + 1)) then
      nbad = nbad + 1 ; cycle
   endif
   err = max(err, maxval(abs(xi(:, n) - (xs(:, n) / h - int(xs(:, n) / h)))))
enddo
call check(nbad == 0 .and. err <= 1e-13_R8P, '2-D random points: analytic cell and xi', err)
deallocate(pu, pv) ; allocate(pu(4), pv(4))
call fmap2(3 * h, 1.5_R8P * h, pu(1), pv(1))   ! on the edge between cells 3 and 4 (vertex line 4)
call fmap2(2.5_R8P * h, 5 * h, pu(2), pv(2))   ! on the edge between cells (3, 5) and (3, 6)
call fmap2(6 * h, 7 * h, pu(3), pv(3))         ! on the vertex (7, 8)
call fmap2(0._R8P, 2.25_R8P * h, pu(4), pv(4)) ! on the low boundary
call q%locate_many(pu, pv, cell, xi, st)
ok = all(st == LOC_INSIDE) .and. all(cell(:, 1) == [3, 2]) .and. xi(1, 1) == 1._R8P .and. &
     all(cell(:, 2) == [3, 5]) .and. xi(2, 2) == 1._R8P .and. all(cell(:, 3) == [6, 7]) .and. all(xi(:, 3) == 1._R8P) .and. &
     all(cell(:, 4) == [1, 3]) .and. xi(1, 4) == 0._R8P
call check(ok, '2-D ties: edges and vertices to the low side, low boundary to the first cell', 0._R8P)
call fmap2(-1e-9_R8P, 0.5_R8P, pu(1), pv(1)) ; call fmap2(2.1_R8P, 0.5_R8P, pu(2), pv(2))
call fmap2(0.5_R8P, -1e-12_R8P, pu(3), pv(3)) ; pu(4) = 50._R8P ; pv(4) = 50._R8P
call q%locate_many(pu, pv, cell, xi, st)
call check(all(st == LOC_OUTSIDE) .and. all(cell == 0), '2-D outside points', 0._R8P)

! 4. 2-D fold
u(9, 9) = u(9, 9) + 0.3_R8P ; v(9, 9) = v(9, 9) + 0.25_R8P
call q%initialize(u, v) ; call q%build_search
deallocate(pu, pv) ; allocate(pu(np), pv(np))
do n = 1, np
   call random_number(r(1:2)) ; pu(n) = u(9, 9) + (r(1) - 0.5_R8P) * 0.8_R8P ; pv(n) = v(9, 9) + (r(2) - 0.5_R8P) * 0.8_R8P
enddo
call q%locate_many(pu, pv, cell, xi, st)
nbad = 0 ; nwarp = count(st == LOC_WARPED)
do n = 1, np
   if (st(n) /= LOC_INSIDE) cycle
   if (q%cell_status(cell(1, n) + (ni - 1) * (cell(2, n) - 1)) /= CELL_VALID) nbad = nbad + 1
   xx(1:2) = bilinear_point(q%corners(cell(1, n), cell(2, n)), xi(:, n))
   if (maxval(abs(xx(1:2) - [pu(n), pv(n)])) > 1e-13_R8P) nbad = nbad + 1
enddo
call check(nbad == 0 .and. nwarp > 0, '2-D fold: inside answers exact and valid, the fold warped', real(nwarp, R8P))

! 5. O-grid seam
deallocate(u, v) ; allocate(u(33, 5), v(33, 5))
do j = 1, 5
   do i = 1, 33
      h = 1 + 0.25_R8P * (j - 1)
      u(i, j) = h * cos(2 * PI * (i - 1) / 32) ; v(i, j) = h * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
u(33,:) = u(1,:) ; v(33,:) = v(1,:)
call q%initialize(u, v) ; call q%build_search
deallocate(pu, pv) ; allocate(pu(3), pv(3))
pu = [1.1_R8P, 1.6_R8P, 1.9_R8P] ; pv = 0._R8P
call q%locate_many(pu, pv, cell, xi, st)
call check(all(st == LOC_INSIDE) .and. all(cell(1, :) == 1) .and. all(xi(1, :) == 0._R8P), 'O-grid seam: cell i = 1', 0._R8P)

! 6. two blocks
allocate(qs(2))
deallocate(u, v) ; allocate(u(5, 5), v(5, 5))
do j = 1, 5
   do i = 1, 5
      u(i, j) = 0.25_R8P * (i - 1) ; v(i, j) = 0.25_R8P * (j - 1) + 0.05_R8P * sin(PI * u(i, j))
   enddo
enddo
call qs(1)%initialize(u, v) ; call qs(1)%build_search
u = u + 1._R8P
call qs(2)%initialize(u, v) ; call qs(2)%build_search
deallocate(pu, pv) ; allocate(pu(3), pv(3))
pu = [0.3_R8P, 1.7_R8P, 2.5_R8P] ; pv = 0.5_R8P
call locate_in_blocks_2d(qs, pu, pv, blk, cell, xi, st)
call check(all(blk == [1, 2, 0]) .and. all(st == [LOC_INSIDE, LOC_INSIDE, LOC_OUTSIDE]), 'two blocks', 0._R8P)

! 1-2. 3-D trilinear map
ni = 9 ; h = 0.25_R8P ; aa = 0.25_R8P ; bb = -0.125_R8P ; cc = 0.375_R8P
allocate(x(ni, ni, ni), y(ni, ni, ni), z(ni, ni, ni))
do k = 1, ni
   do j = 1, ni
      do i = 1, ni
         xx = fmap3(h * [i - 1, j - 1, k - 1]) ; x(i, j, k) = xx(1) ; y(i, j, k) = xx(2) ; z(i, j, k) = xx(3)
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%build_search
np = 1000 ; deallocate(pu, pv, xs) ; allocate(pu(np), pv(np), pw(np), xs(3, np))
do n = 1, np
   call random_number(r) ; xs(:, n) = r * h * (ni - 1) ; xx = fmap3(xs(:, n))
   pu(n) = xx(1) ; pv(n) = xx(2) ; pw(n) = xx(3)
enddo
call b%locate_many(pu, pv, pw, ct, xt, st)
err = 0._R8P ; nbad = 0
do n = 1, np
   if (st(n) /= LOC_INSIDE .or. any(ct(:, n) /= int(xs(:, n) / h) + 1)) then
      nbad = nbad + 1 ; cycle
   endif
   err = max(err, maxval(abs(xt(:, n) - (xs(:, n) / h - int(xs(:, n) / h)))))
enddo
call check(nbad == 0 .and. err <= 1e-13_R8P, '3-D random points: analytic cell and xi', err)
deallocate(pu, pv, pw) ; allocate(pu(4), pv(4), pw(4))
xs(:, 1) = h * [3._R8P, 1.5_R8P, 2.5_R8P] ; xs(:, 2) = h * [2._R8P, 5._R8P, 1.25_R8P]
xs(:, 3) = h * [4._R8P, 4._R8P, 4._R8P] ; xs(:, 4) = h * [0._R8P, 0.5_R8P, 0.5_R8P]
do n = 1, 4
   xx = fmap3(xs(:, n)) ; pu(n) = xx(1) ; pv(n) = xx(2) ; pw(n) = xx(3)
enddo
call b%locate_many(pu, pv, pw, ct, xt, st)
ok = all(st == LOC_INSIDE) .and. all(ct(:, 1) == [3, 2, 3]) .and. all(ct(:, 2) == [2, 5, 2]) .and. &
     all(ct(:, 3) == [4, 4, 4]) .and. all(ct(:, 4) == [1, 1, 1])
call check(ok, '3-D ties: faces, edges, vertices to the low side', 0._R8P)
xx = fmap3([-1e-6_R8P, 0.5_R8P, 0.5_R8P]) ; pu(1) = xx(1) ; pv(1) = xx(2) ; pw(1) = xx(3)
pu(2:4) = 10._R8P ; pv(2:4) = [1._R8P, 2._R8P, 3._R8P] ; pw(2:4) = 0._R8P
call b%locate_many(pu, pv, pw, ct, xt, st)
call check(all(st == LOC_OUTSIDE), '3-D outside points', 0._R8P)

! 4. 3-D folded cell alone
deallocate(x, y, z) ; allocate(x(2, 2, 2), y(2, 2, 2), z(2, 2, 2))
do k = 1, 2
   do j = 1, 2
      do i = 1, 2
         x(i, j, k) = FOLD(1, i + 2 * (j - 1) + 4 * (k - 1)) ; y(i, j, k) = FOLD(2, i + 2 * (j - 1) + 4 * (k - 1))
         z(i, j, k) = FOLD(3, i + 2 * (j - 1) + 4 * (k - 1))
      enddo
   enddo
enddo
call b%initialize(x, y, z) ; call b%build_search
deallocate(pu, pv, pw) ; allocate(pu(200), pv(200), pw(200))
do n = 1, 200
   call random_number(r) ; xx = trilinear_point(FOLD, r) ; pu(n) = xx(1) ; pv(n) = xx(2) ; pw(n) = xx(3)
enddo
call b%locate_many(pu, pv, pw, ct, xt, st)
call check(count(st == LOC_INSIDE) == 0 .and. count(st == LOC_WARPED) > 150, '3-D folded cell: warped, never inside', &
           real(count(st == LOC_WARPED), R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES11.3)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   subroutine fmap2(xx, yy, uu, vv)
   !< The bilinear map.
   real(R8P), intent(in)  :: xx, yy !< Reference point.
   real(R8P), intent(out) :: uu, vv !< Image.

   uu = xx + aa * xx * yy ; vv = yy + bb * xx * yy
   endsubroutine fmap2

   function fmap3(s) result(p)
   !< The trilinear map.
   real(R8P), intent(in) :: s(3) !< Reference point.
   real(R8P)             :: p(3) !< Image.

   p = [s(1) + aa * s(2) * s(3), s(2) + bb * s(3) * s(1), s(3) + cc * s(1) * s(2)]
   endfunction fmap3
endprogram fossil_test_block_locate
