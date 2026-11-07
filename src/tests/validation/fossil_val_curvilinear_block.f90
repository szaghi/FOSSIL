!< FOSSIL validation: curvilinear blocks at size (issue #71, step 10).
program fossil_val_curvilinear_block
!< FOSSIL validation: curvilinear blocks at size (issue #71, step 10).
!<
!< 1. A spherical-shell O-grid around the unit sphere (129 x 65 x 33 vertices: azimuth with a seam, polar angle with
!<    both poles collapsed, radius stretched to 4): topology (the seam on faces 1 and 2, the poles singular), validity
!<    statistics, 200000 random points in the shell located by the robust path and by the walk from nearby start
!<    cells (equal bit for bit), timings; interpolation order of the width-2 and width-4 order-matched stencils on
!<    the shell (32 and 64 azimuthal cells).
!< 2. A parabolic C-grid (801 x 201 vertices) around a smooth cusped body: the wake cut, 1000000 points located robustly and by the walk from the
!<    previous point along a sorted sweep (equal), timings; orders of q = 2, 3, 4 on three refinements (400 x 100 to
!<    1600 x 400 cells, the slope over the range).
!< Exit status 0 if every identity holds and every observed order is within 0.3 of q.

use fossil, only : curvilinear_block_t, curvilinear_block_2d_t, detect_topology, detect_topology_2d, TOPO_ABUT, &
                   TOPO_SINGULAR, CELL_VALID, CELL_SINGULAR, LOC_INSIDE
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter         :: PI = acos(-1._R8P)
type(curvilinear_block_t)    :: shell(1)
type(curvilinear_block_2d_t) :: cg(1)
real(R8P),    allocatable    :: x(:,:,:), y(:,:,:), z(:,:,:), u(:,:), v(:,:), px(:), py(:), pz(:), x0(:,:), x1(:,:)
integer(I4P), allocatable    :: st(:,:,:), c0(:,:), c1(:,:), s0(:), s1(:), start(:,:)
real(R8P)                    :: r(3), e(3), rate, t
integer(I4P)                 :: n, nj, nf, q, l
integer(I8P)                 :: t0, t1, t2, cr
logical                      :: ok

ok = .true.
call random_seed(put=[(13 + 17 * n, n = 1, 64)])

! 1. spherical shell
call shell_grid(128, 64, 32, x, y, z)
call shell(1)%initialize(x, y, z, name='shell') ; call detect_topology(shell)
print '(A,I0,A,I0,A,I0)', 'shell: records ', shell(1)%nrec, ', abutting ', count(shell(1)%rec(:)%kind == TOPO_ABUT), &
      ', singular ', count(shell(1)%rec(:)%kind == TOPO_SINGULAR)
ok = ok .and. count(shell(1)%rec(:)%kind == TOPO_ABUT) == 2 .and. count(shell(1)%rec(:)%kind == TOPO_SINGULAR) == 2
call system_clock(t0, cr)
call shell(1)%validity(st)
call system_clock(t1)
print '(A,I0,A,I0,A,I0,A,F8.2,A)', 'shell: cells ', size(st), ', valid ', count(st == CELL_VALID), ', singular ', &
      count(st == CELL_SINGULAR), ' (', real(t1 - t0, R8P) / cr, ' s)'
ok = ok .and. count(st == CELL_VALID) + count(st == CELL_SINGULAR) == size(st)
call shell(1)%build_search
allocate(px(200000), py(200000), pz(200000))
do n = 1, 200000
   call random_number(r) ; t = 1.05_R8P + 2.8_R8P * r(1)
   r(2) = acos(1 - 2 * r(2)) ; r(3) = 2 * PI * r(3)
   px(n) = t * sin(r(2)) * cos(r(3)) ; py(n) = t * sin(r(2)) * sin(r(3)) ; pz(n) = t * cos(r(2))
enddo
call system_clock(t0)
call shell(1)%locate_many(px, py, pz, c0, x0, s0)
call system_clock(t1)
allocate(start(3, size(px)))
do n = 1, size(px)
   start(:, n) = max(1, min(c0(:, n) + [2, -1, 1], [127, 63, 31]))
enddo
call shell(1)%locate_many(px, py, pz, c1, x1, s1, walk=.true., start=start, jumps=nj, fallbacks=nf)
call system_clock(t2)
print '(A,I0,A,F8.0,A,F8.0,A,I0,A,I0)', 'shell: inside ', count(s0 == LOC_INSIDE), ', robust ', &
      size(px) / (real(t1 - t0, R8P) / cr), ' pts/s, walk ', size(px) / (real(t2 - t1, R8P) / cr), ' pts/s, jumps ', &
      nj, ', fallbacks ', nf
ok = ok .and. all(c0 == c1) .and. all(x0 == x1) .and. all(s0 == s1) .and. count(s0 == LOC_INSIDE) == size(px)
do q = 2, 4, 2
   e(1) = shell_error(32, q) ; e(2) = shell_error(64, q)
   rate = log(e(1) / e(2)) / log(2._R8P)
   print '(A,I0,A,2ES10.2,A,F6.2)', 'shell: q = ', q, ', errors ', e(1:2), ', order ', rate
   ok = ok .and. rate >= q - 0.3_R8P
enddo
deallocate(px, py, pz, start)

! 2. C-grid
call c_grid(800, 200, u, v)
call cg(1)%initialize(u, v) ; call detect_topology_2d(cg) ; call cg(1)%build_search
print '(A,I0)', 'C-grid: records ', cg(1)%nrec
ok = ok .and. cg(1)%nrec == 2
allocate(px(1000000), py(1000000))
do n = 1, 1000000 ! a sorted sweep: rows of points
   call random_number(r(1:2))
   px(n) = -1.5_R8P + 5.5_R8P * real(mod(n - 1, 1000), R8P) / 1000 + 0.004_R8P * r(1)
   py(n) = -4.5_R8P + 9._R8P * real((n - 1) / 1000, R8P) / 1000 + 0.004_R8P * r(2)
enddo
call system_clock(t0)
call cg(1)%locate_many(px, py, c0, x0, s0)
call system_clock(t1)
call cg(1)%locate_many(px, py, c1, x1, s1, walk=.true., jumps=nj, fallbacks=nf)
call system_clock(t2)
print '(A,I0,A,F10.0,A,F10.0,A,I0,A,I0)', 'C-grid: inside ', count(s0 == LOC_INSIDE), ', robust ', &
      size(px) / (real(t1 - t0, R8P) / cr), ' pts/s, walk ', size(px) / (real(t2 - t1, R8P) / cr), ' pts/s, jumps ', &
      nj, ', fallbacks ', nf
ok = ok .and. all(c0 == c1) .and. all(x0 == x1) .and. all(s0 == s1)
do q = 2, 4
   do l = 1, 3
      e(l) = cgrid_error(100 * 2 ** (l - 1), q)
   enddo
   rate = log(e(1) / e(3)) / log(4._R8P) ! slope over the whole range (max norms over random points are noisy)
   print '(A,I0,A,3ES10.2,A,F6.2)', 'C-grid: q = ', q, ', errors ', e, ', order ', rate
   ok = ok .and. rate >= q - 0.3_R8P
enddo

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine shell_grid(na, np, nr, x, y, z)
   !< Spherical shell: azimuth (na cells, seam copied), polar angle (np cells, poles collapsed), radius 1..4 (nr).
   integer(I4P),           intent(in)  :: na, np, nr !< Cells.
   real(R8P), allocatable, intent(out) :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)                        :: i, j, k  !< Counters.
   real(R8P)                           :: rr, th, ph !< Spherical coordinates.

   allocate(x(na + 1, np + 1, nr + 1), y(na + 1, np + 1, nr + 1), z(na + 1, np + 1, nr + 1))
   do k = 1, nr + 1
      rr = 1._R8P + 3._R8P * ((k - 1._R8P) / nr) ** 1.5_R8P
      do j = 1, np + 1
         ph = PI * (j - 1) / np
         do i = 1, na + 1
            th = 2 * PI * (i - 1) / na
            x(i, j, k) = rr * sin(ph) * cos(th) ; y(i, j, k) = rr * sin(ph) * sin(th) ; z(i, j, k) = rr * cos(ph)
         enddo
         if (j == 1 .or. j == np + 1) then
            x(:, j, k) = 0._R8P ; y(:, j, k) = 0._R8P ; z(:, j, k) = merge(rr, -rr, j == 1)
         endif
      enddo
   enddo
   x(na + 1,:,:) = x(1,:,:) ; y(na + 1,:,:) = y(1,:,:) ; z(na + 1,:,:) = z(1,:,:)
   endsubroutine shell_grid

   real(R8P) function shell_error(n, q) result(err)
   !< Max interpolation error of exp(x/3) cos(y/2) z at 300 points of the shell (r in [1.6, 2.4], away from the poles).
   integer(I4P), intent(in)  :: n, q !< Refinement (cells n, n / 2, n / 4), width.
   type(curvilinear_block_t) :: g(1) !< Grid.
   real(R8P),    allocatable :: gx(:,:,:), gy(:,:,:), gz(:,:,:), qx(:), qy(:), qz(:), xi(:,:), w(:,:) !< Data.
   integer(I4P), allocatable :: gc(:,:), gs(:), nd(:,:,:), ns(:) !< Results.
   integer(I4P)              :: m, k !< Counters.
   real(R8P)                 :: rr(3), rad, f !< Work.

   call shell_grid(n, n / 2, n / 4, gx, gy, gz)
   call g(1)%initialize(gx, gy, gz) ; call detect_topology(g) ; call g(1)%build_search
   allocate(qx(300), qy(300), qz(300))
   call random_seed(put=[(29 + 3 * m, m = 1, 64)])
   do m = 1, 300
      call random_number(rr) ; rad = 1.6_R8P + 0.8_R8P * rr(1) ; rr(2) = 0.5_R8P + (PI - 1) * rr(2)
      rr(3) = 2 * PI * rr(3)
      qx(m) = rad * sin(rr(2)) * cos(rr(3)) ; qy(m) = rad * sin(rr(2)) * sin(rr(3)) ; qz(m) = rad * cos(rr(2))
   enddo
   call g(1)%locate_many(qx, qy, qz, gc, xi, gs)
   call g(1)%stencils(gc, xi, q, nd, w, ns)
   err = 0._R8P
   do m = 1, 300
      f = 0._R8P
      do k = 1, q ** 3
         f = f + w(k, m) * fn3(gx(nd(1, k, m), nd(2, k, m), nd(3, k, m)), gy(nd(1, k, m), nd(2, k, m), nd(3, k, m)), &
                               gz(nd(1, k, m), nd(2, k, m), nd(3, k, m)))
      enddo
      err = max(err, abs(f - fn3(qx(m), qy(m), qz(m))))
   enddo
   endfunction shell_error

   pure real(R8P) function fn3(a, b, c) result(f)
   !< Smooth test function.
   real(R8P), intent(in) :: a, b, c !< Point.

   f = exp(a / 3) * cos(b / 2) * c
   endfunction fn3

   subroutine c_grid(ni, nj, u, v)
   !< Parabolic C-grid z = zeta^2, zeta = s + i (g(s) + t), s in [-2, 2] (ni cells), t in [0, 1.5] (nj cells);
   !< g(s) = 0.3 e^(1 - 1 / (1 - s^2)) on |s| < 1, 0 on the wake (smooth: the interpolation order needs a smooth map).
   integer(I4P),           intent(in)  :: ni, nj !< Cells.
   real(R8P), allocatable, intent(out) :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)                        :: i, j   !< Counters.
   real(R8P)                           :: s, h   !< Parameters.

   allocate(u(ni + 1, nj + 1), v(ni + 1, nj + 1))
   do i = 1, ni / 2 + 1
      s = 2._R8P * (2._R8P * (i - 1) / ni - 1._R8P)
      do j = 1, nj + 1
         h = bump(s) + 1.5_R8P * (j - 1) / nj
         u(i, j) = s ** 2 - h ** 2 ; v(i, j) = 2 * s * h
         u(ni + 2 - i, j) = u(i, j) ; v(ni + 2 - i, j) = -v(i, j)
      enddo
   enddo
   endsubroutine c_grid

   pure real(R8P) function bump(s) result(g)
   !< The body thickness: 0.3 e^(1 - 1 / (1 - s^2)) on |s| < 1, else 0 (C-infinity).
   real(R8P), intent(in) :: s !< Parameter.

   g = 0._R8P
   if (abs(s) < 1._R8P) g = 0.3_R8P * exp(1._R8P - 1._R8P / (1._R8P - s * s))
   endfunction bump

   real(R8P) function cgrid_error(n, q) result(err)
   !< Max interpolation error of sin(x) exp(y / 4) at 400 points of the C-grid (s in [-1.6, 1.6], t in [0.4, 1.1]).
   integer(I4P), intent(in)     :: n, q !< Refinement (4 n x n cells), width.
   type(curvilinear_block_2d_t) :: g(1) !< Grid.
   real(R8P),    allocatable    :: gu(:,:), gv(:,:), qu(:), qv(:), xi(:,:), w(:,:) !< Data.
   integer(I4P), allocatable    :: gc(:,:), gs(:), nd(:,:,:), ns(:) !< Results.
   integer(I4P)                 :: m, k !< Counters.
   real(R8P)                    :: rr(2), s, h, f !< Work.

   call c_grid(4 * n, n, gu, gv)
   call g(1)%initialize(gu, gv) ; call detect_topology_2d(g) ; call g(1)%build_search
   allocate(qu(400), qv(400))
   call random_seed(put=[(31 + 7 * m, m = 1, 64)])
   do m = 1, 400
      call random_number(rr) ; s = -1.6_R8P + 3.2_R8P * rr(1)
      h = bump(s) + 0.4_R8P + 0.7_R8P * rr(2)
      qu(m) = s ** 2 - h ** 2 ; qv(m) = 2 * s * h
   enddo
   call g(1)%locate_many(qu, qv, gc, xi, gs)
   call g(1)%stencils(gc, xi, q, nd, w, ns)
   err = 0._R8P
   do m = 1, 400
      f = 0._R8P
      do k = 1, q * q
         f = f + w(k, m) * sin(gu(nd(1, k, m), nd(2, k, m))) * exp(gv(nd(1, k, m), nd(2, k, m)) / 4)
      enddo
      err = max(err, abs(f - sin(qu(m)) * exp(qv(m) / 4)))
   enddo
   endfunction cgrid_error
endprogram fossil_val_curvilinear_block
