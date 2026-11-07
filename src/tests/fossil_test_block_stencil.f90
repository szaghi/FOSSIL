!< FOSSIL test: Lagrange interpolation stencils on curvilinear blocks (issue #71, step 8).
program fossil_test_block_stencil
!< FOSSIL test: Lagrange interpolation stencils on curvilinear blocks (issue #71, step 8).
!<
!< 1. index space: the weights reproduce s1^a s2^b (s3^c) for every degree below q in each direction (q = 2..5, 1e-12);
!< 2. observed order on a smooth curved 2-D grid (N = 16, 32, 64): q = 2, 3, 4 with the order-matched xi converge at
!<    order >= q - 0.3; q = 4 with the trilinear xi stalls near second order (the reason for the order matching);
!<    3-D: q = 4 order-matched at order >= 3.7 (N = 8, 16);
!< 3. across an O-grid seam: stencils on the closed grid (nodes mapped through the seam record) interpolate what the
!<    same stencils give on the unwrapped grid (1e-13), at points on both sides of the seam;
!< 4. at the physical boundary: shifted one-sided (STENCIL_SHIFTED) or narrowed (STENCIL_REDUCED); weights sum to 1.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, detect_topology_2d
use fossil_block_stencil,     only : STENCIL_CENTRED, STENCIL_SHIFTED, STENCIL_REDUCED
use fossil_block_locate,      only : LOC_INSIDE
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q2(:)
type(curvilinear_block_2d_t)              :: w2
type(curvilinear_block_t)                 :: b3
real(R8P),    allocatable                 :: u(:,:), v(:,:), x(:,:,:), y(:,:,:), z(:,:,:), xi(:,:), wt(:,:), xs(:,:)
real(R8P),    allocatable                 :: pu(:), pv(:), uu(:,:), vv(:,:), xi2(:,:), wt2(:,:)
integer(I4P), allocatable                 :: cell(:,:), st(:), nodes(:,:,:), ss(:), cell2(:,:), nodes2(:,:,:), ss2(:)
real(R8P)                                 :: err, e(3), rate, rr(3), h, f1, f2, ex
integer(I4P)                              :: i, j, k, n, q, a, b, c, l, ni, nj, ig
logical                                   :: ok, all_ok

all_ok = .true.
call random_seed(put=[(3 + 11 * i, i = 1, 64)])

! 1. polynomial exactness in index space (2-D and 3-D, xi as given)
ni = 12
call grid2(ni, 0.08_R8P, u, v) ; allocate(q2(1)) ; call q2(1)%initialize(u, v)
allocate(cell(2, 50), xi(2, 50))
do n = 1, 50
   call random_number(rr(1:2)) ; cell(:, n) = 3 + int(rr(1:2) * 6) ; call random_number(rr(1:2)) ; xi(:, n) = rr(1:2)
enddo
err = 0._R8P
do q = 2, 5
   call q2(1)%stencils(cell, xi, q, nodes, wt, ss, order_matched=.false.)
   do n = 1, 50
      do b = 0, q - 1
         do a = 0, q - 1
            ex = (cell(1, n) + xi(1, n)) ** a * (cell(2, n) + xi(2, n)) ** b
            f1 = sum(wt(:, n) * real(nodes(1, :, n), R8P) ** a * real(nodes(2, :, n), R8P) ** b)
            err = max(err, abs(f1 - ex) / max(1._R8P, abs(ex)))
         enddo
      enddo
   enddo
enddo
call check(err <= 1e-12_R8P, '2-D weights reproduce index polynomials of degree < q (q = 2..5)', err)
deallocate(cell, xi)
call grid3(8, x, y, z) ; call b3%initialize(x, y, z)
allocate(cell(3, 20), xi(3, 20))
do n = 1, 20
   call random_number(rr) ; cell(:, n) = 3 + int(rr * 3) ; call random_number(rr) ; xi(:, n) = rr
enddo
err = 0._R8P
call b3%stencils(cell, xi, 4, nodes, wt, ss, order_matched=.false.)
do n = 1, 20
   do c = 0, 3
      do b = 0, 3
         do a = 0, 3
            ex = (cell(1, n) + xi(1, n)) ** a * (cell(2, n) + xi(2, n)) ** b * (cell(3, n) + xi(3, n)) ** c
            f1 = sum(wt(:, n) * real(nodes(1, :, n), R8P) ** a * real(nodes(2, :, n), R8P) ** b * &
                     real(nodes(3, :, n), R8P) ** c)
            err = max(err, abs(f1 - ex) / max(1._R8P, abs(ex)))
         enddo
      enddo
   enddo
enddo
call check(err <= 1e-12_R8P, '3-D weights reproduce index polynomials of degree < 4 per direction', err)
deallocate(cell, xi)

! 2. observed order, 2-D
do q = 2, 4
   do l = 1, 3
      e(l) = interp_error_2d(16 * 2 ** (l - 1), q, .true.)
   enddo
   rate = log(e(2) / e(3)) / log(2._R8P)
   call check(rate >= q - 0.3_R8P, '2-D order-matched, q = ' // achar(48 + q) // ': observed order', rate)
enddo
do l = 1, 3
   e(l) = interp_error_2d(16 * 2 ** (l - 1), 4_I4P, .false.)
enddo
rate = log(e(2) / e(3)) / log(2._R8P)
call check(rate < 2.6_R8P, '2-D q = 4 with the trilinear xi: order stalls near 2', rate)
e(1) = interp_error_3d(8_I4P) ; e(2) = interp_error_3d(16_I4P)
rate = log(e(1) / e(2)) / log(2._R8P)
call check(rate >= 3.7_R8P, '3-D order-matched, q = 4: observed order', rate)

! 3. across an O-grid seam
ni = 41 ; nj = 9
deallocate(u, v) ; allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   h = 1._R8P + 0.1_R8P * (j - 1)
   do i = 1, ni
      u(i, j) = h * cos(2 * PI * (i - 1) / (ni - 1)) ; v(i, j) = 0.8_R8P * h * sin(2 * PI * (i - 1) / (ni - 1))
   enddo
enddo
u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
call q2(1)%initialize(u, v) ; call detect_topology_2d(q2) ; call q2(1)%build_search
allocate(uu(ni + 8, nj), vv(ni + 8, nj)) ! unwrapped: 4 columns of ghosts on each side, copied
do i = 1, ni + 8
   ig = modulo(i - 5 - 1, ni - 1) + 1
   uu(i, :) = u(ig, :) ; vv(i, :) = v(ig, :)
enddo
call w2%initialize(uu, vv) ; call w2%build_search
allocate(pu(40), pv(40))
do n = 1, 40
   call random_number(rr(1:2))
   h = 1.25_R8P + 0.3_R8P * rr(1) ; f1 = (rr(2) - 0.5_R8P) * 0.6_R8P
   pu(n) = h * cos(f1) ; pv(n) = 0.8_R8P * h * sin(f1)
enddo
call q2(1)%locate_many(pu, pv, cell, xi, st)
call w2%locate_many(pu, pv, cell2, xi2, ss2)
ok = all(st == LOC_INSIDE) .and. all(ss2 == LOC_INSIDE)
err = 0._R8P
do q = 2, 5
   call q2(1)%stencils(cell, xi, q, nodes, wt, ss)
   call w2%stencils(cell2, xi2, q, nodes2, wt2, ss2)
   ok = ok .and. all(ss == STENCIL_CENTRED) .and. all(ss2 == STENCIL_CENTRED)
   do n = 1, 40
      f1 = 0._R8P ; f2 = 0._R8P
      do k = 1, q * q
         f1 = f1 + wt(k, n) * fun(u(nodes(1, k, n), nodes(2, k, n)), v(nodes(1, k, n), nodes(2, k, n)))
         f2 = f2 + wt2(k, n) * fun(uu(nodes2(1, k, n), nodes2(2, k, n)), vv(nodes2(1, k, n), nodes2(2, k, n)))
      enddo
      err = max(err, abs(f1 - f2))
   enddo
enddo
call check(ok .and. err <= 1e-13_R8P, 'across the O-grid seam: closed grid = unwrapped grid (q = 2..5)', err)

! 4. physical boundary
deallocate(pu, pv) ; allocate(pu(1), pv(1))
pu = 1.78_R8P ; pv = 0.05_R8P ! in the outermost ring of cells
call q2(1)%locate_many(pu, pv, cell, xi, st)
call q2(1)%stencils(cell, xi, 4, nodes, wt, ss)
ok = ss(1) == STENCIL_SHIFTED .and. abs(sum(wt(:, 1)) - 1._R8P) <= 1e-14_R8P .and. maxval(nodes(2, :, 1)) == nj
call q2(1)%stencils(cell, xi, 4, nodes, wt, ss, reduce=.true.)
ok = ok .and. ss(1) == STENCIL_REDUCED .and. abs(sum(wt(:, 1)) - 1._R8P) <= 1e-14_R8P
call check(ok, 'outer boundary: shifted, or reduced; weights sum to 1', 0._R8P)

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

   pure real(R8P) function fun(px, py) result(f)
   !< Smooth test function.
   real(R8P), intent(in) :: px, py !< Point.

   f = sin(2 * px + 1) * cos(3 * py - 0.5_R8P)
   endfunction fun

   subroutine grid2(n, amp, u, v)
   !< Smooth curved grid: x = X + amp sin(2 pi Y), y = Y + amp sin(2 pi X), on [0, 1]^2.
   integer(I4P),           intent(in)  :: n      !< Cells per side.
   real(R8P),              intent(in)  :: amp    !< Amplitude.
   real(R8P), allocatable, intent(out) :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)                        :: i, j   !< Counters.

   allocate(u(n + 1, n + 1), v(n + 1, n + 1))
   do j = 1, n + 1
      do i = 1, n + 1
         call map2((i - 1._R8P) / n, (j - 1._R8P) / n, amp, u(i, j), v(i, j))
      enddo
   enddo
   endsubroutine grid2

   pure subroutine map2(xx, yy, amp, uu, vv)
   !< The curved map.
   real(R8P), intent(in)  :: xx, yy, amp !< Reference point, amplitude.
   real(R8P), intent(out) :: uu, vv      !< Image.

   uu = xx + amp * sin(2 * PI * yy) ; vv = yy + amp * sin(2 * PI * xx)
   endsubroutine map2

   real(R8P) function interp_error_2d(n, q, matched) result(err)
   !< Max interpolation error of fun at 300 points of [0.3, 0.7]^2 (reference) on the curved grid with n cells.
   integer(I4P), intent(in) :: n, q   !< Cells per side, width.
   logical,      intent(in) :: matched !< Order-matched xi.
   type(curvilinear_block_2d_t) :: g   !< Grid.
   real(R8P),    allocatable    :: gu(:,:), gv(:,:), qu(:), qv(:), gx(:,:), gw(:,:) !< Data.
   integer(I4P), allocatable    :: gc(:,:), gs(:), gn(:,:,:), gss(:) !< Results.
   integer(I4P)                 :: m, k  !< Counters.
   real(R8P)                    :: r(2), f !< Work.

   call grid2(n, 0.08_R8P, gu, gv) ; call g%initialize(gu, gv) ; call g%build_search
   allocate(qu(300), qv(300))
   call random_seed(put=[(5 + 3 * m, m = 1, 64)])
   do m = 1, 300
      call random_number(r) ; r = 0.3_R8P + 0.4_R8P * r ; call map2(r(1), r(2), 0.08_R8P, qu(m), qv(m))
   enddo
   call g%locate_many(qu, qv, gc, gx, gs)
   call g%stencils(gc, gx, q, gn, gw, gss, order_matched=matched)
   err = 0._R8P
   do m = 1, 300
      f = 0._R8P
      do k = 1, q * q
         f = f + gw(k, m) * fun(gu(gn(1, k, m), gn(2, k, m)), gv(gn(1, k, m), gn(2, k, m)))
      enddo
      err = max(err, abs(f - fun(qu(m), qv(m))))
   enddo
   endfunction interp_error_2d

   subroutine grid3(n, x, y, z)
   !< Smooth curved 3-D grid on [0, 1]^3.
   integer(I4P),           intent(in)  :: n         !< Cells per side.
   real(R8P), allocatable, intent(out) :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)                        :: i, j, k   !< Counters.
   real(R8P)                           :: p(3)      !< Image.

   allocate(x(n + 1, n + 1, n + 1), y(n + 1, n + 1, n + 1), z(n + 1, n + 1, n + 1))
   do k = 1, n + 1
      do j = 1, n + 1
         do i = 1, n + 1
            p = map3([i - 1._R8P, j - 1._R8P, k - 1._R8P] / n)
            x(i, j, k) = p(1) ; y(i, j, k) = p(2) ; z(i, j, k) = p(3)
         enddo
      enddo
   enddo
   endsubroutine grid3

   pure function map3(s) result(p)
   !< The curved 3-D map.
   real(R8P), intent(in) :: s(3) !< Reference point.
   real(R8P)             :: p(3) !< Image.

   p = s + 0.06_R8P * [sin(2 * PI * s(2)) * cos(PI * s(3)), sin(2 * PI * s(3)) * cos(PI * s(1)), &
                       sin(2 * PI * s(1)) * cos(PI * s(2))]
   endfunction map3

   real(R8P) function interp_error_3d(n) result(err)
   !< Max interpolation error (q = 4, order-matched) of g = fun(x, y) cos(z) at 60 points of [0.35, 0.65]^3.
   integer(I4P), intent(in)  :: n   !< Cells per side.
   type(curvilinear_block_t) :: g   !< Grid.
   real(R8P),    allocatable :: gx(:,:,:), gy(:,:,:), gz(:,:,:), qx(:), qy(:), qz(:), xx(:,:), ww(:,:) !< Data.
   integer(I4P), allocatable :: gc(:,:), gs(:), gn(:,:,:), gss(:) !< Results.
   integer(I4P)              :: m, k  !< Counters.
   real(R8P)                 :: r(3), p(3), f !< Work.

   call grid3(n, gx, gy, gz) ; call g%initialize(gx, gy, gz) ; call g%build_search
   allocate(qx(60), qy(60), qz(60))
   call random_seed(put=[(9 + 5 * m, m = 1, 64)])
   do m = 1, 60
      call random_number(r) ; p = map3(0.35_R8P + 0.3_R8P * r) ; qx(m) = p(1) ; qy(m) = p(2) ; qz(m) = p(3)
   enddo
   call g%locate_many(qx, qy, qz, gc, xx, gs)
   call g%stencils(gc, xx, 4_I4P, gn, ww, gss)
   err = 0._R8P
   do m = 1, 60
      f = 0._R8P
      do k = 1, 64
         f = f + ww(k, m) * fun(gx(gn(1, k, m), gn(2, k, m), gn(3, k, m)), gy(gn(1, k, m), gn(2, k, m), gn(3, k, m))) * &
             cos(gz(gn(1, k, m), gn(2, k, m), gn(3, k, m)))
      enddo
      err = max(err, abs(f - fun(qx(m), qy(m)) * cos(qz(m))))
   enddo
   endfunction interp_error_3d
endprogram fossil_test_block_stencil
