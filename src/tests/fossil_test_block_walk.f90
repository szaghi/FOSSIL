!< FOSSIL test: the walking point location in curvilinear blocks (issue #71, step 6).
program fossil_test_block_walk
!< FOSSIL test: the walking point location in curvilinear blocks (issue #71, step 6).
!<
!< The walk (stencil jumping from the previous point's cell, or from given start cells) must return exactly what the
!< robust path returns: the same cell, the same xi bit for bit and the same status, on
!< 1. a 2-D O-grid around a cylinder (seam record found by detect_topology_2d),
!< 2. a 2-D C-grid (the parabolic grid z = zeta^2 around a cusped body) with its wake cut,
!< 3. a 3-D polar grid with a singular axis and a seam,
!< 4. a 3-D wavy O-grid,
!< with random points in the bounding box (some outside), points on the seams and cuts, and start cells a few cells
!< away from the answer. The counts of jumps and fallbacks are printed.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, detect_topology, detect_topology_2d
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
real(R8P),    allocatable                 :: u(:,:), v(:,:), x(:,:,:), y(:,:,:), z(:,:,:)
integer(I4P)                              :: i, j, k, ni, nj, nk
real(R8P)                                 :: r, t, p(2)
logical                                   :: all_ok

all_ok = .true.
call random_seed(put=[(11 + 7 * i, i = 1, 64)])
allocate(q(1), b(1))

! 1. O-grid
ni = 49 ; nj = 13
allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   r = 1._R8P + 1.5_R8P * ((j - 1._R8P) / (nj - 1)) ** 1.3_R8P
   do i = 1, ni
      t = 2 * PI * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
   enddo
enddo
u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
call q(1)%initialize(u, v) ; call detect_topology_2d(q) ; call q(1)%build_search
call compare_2d(q(1), 'O-grid', [2.5_R8P, 2.5_R8P], [1.3_R8P, 0._R8P])
deallocate(u, v)

! 2. C-grid: z = zeta^2 (injective on Im zeta > 0) of zeta = s + i (g(s) + (j - 1) d), g = 0 on the wake |s| >= 1
ni = 81 ; nj = 11
allocate(u(ni, nj), v(ni, nj))
do i = 1, (ni + 1) / 2
   p(1) = 2._R8P * (2._R8P * (i - 1) / (ni - 1) - 1._R8P) ! s, from -2 to 0
   do j = 1, nj
      r = 0.3_R8P * max(1._R8P - p(1) ** 2, 0._R8P) ** 2 + 0.1_R8P * (j - 1)
      u(i, j) = p(1) ** 2 - r ** 2 ; v(i, j) = 2 * p(1) * r
      u(ni + 1 - i, j) = u(i, j) ; v(ni + 1 - i, j) = -v(i, j) ! s -> -s
   enddo
enddo
call q(1)%initialize(u, v) ; call detect_topology_2d(q) ; call q(1)%build_search
call check(q(1)%nrec == 2, 'C-grid: wake cut found')
call compare_2d(q(1), 'C-grid', [4.5_R8P, 5._R8P], [2.2_R8P, 0._R8P])
deallocate(u, v)

! 3. polar grid
ni = 25 ; nj = 8 ; nk = 6
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         r = 0.2_R8P * (j - 1) * (1 + 0.1_R8P * k)
         x(i, j, k) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; y(i, j, k) = r * sin(2 * PI * (i - 1) / (ni - 1))
         z(i, j, k) = 0.3_R8P * (k - 1)
      enddo
   enddo
enddo
x(ni,:,:) = x(1,:,:) ; y(ni,:,:) = y(1,:,:)
call b(1)%initialize(x, y, z) ; call detect_topology(b) ; call b(1)%build_search
call compare_3d(b(1), 'polar grid', [2._R8P, 2._R8P, 1.5_R8P], [0.7_R8P, 0._R8P, 0.65_R8P])

! 4. wavy O-grid
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         r = (1._R8P + 0.2_R8P * (j - 1)) * (1._R8P + 0.05_R8P * sin(3 * 2 * PI * (i - 1) / (ni - 1) + 0.4_R8P * k))
         x(i, j, k) = r * cos(2 * PI * (i - 1) / (ni - 1)) ; y(i, j, k) = r * sin(2 * PI * (i - 1) / (ni - 1))
         z(i, j, k) = 0.2_R8P * (k - 1) + 0.03_R8P * cos(2 * PI * (i - 1) / (ni - 1)) * j
      enddo
   enddo
enddo
x(ni,:,:) = x(1,:,:) ; y(ni,:,:) = y(1,:,:) ; z(ni,:,:) = z(1,:,:)
call b(1)%initialize(x, y, z) ; call detect_topology(b) ; call b(1)%build_search
call compare_3d(b(1), 'wavy O-grid', [2.6_R8P, 2.6_R8P, 1.2_R8P], [1.5_R8P, 0._R8P, 0.5_R8P])

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine compare_2d(blk, label, half, seam)
   !< Robust against walked location on random points of the box [-half, half] and on a seam point.
   type(curvilinear_block_2d_t), intent(in) :: blk   !< Block.
   character(*),                 intent(in) :: label !< Label.
   real(R8P),                    intent(in) :: half(2), seam(2) !< Box half sizes, a point on a seam.
   integer(I4P), parameter                  :: NP = 3000 !< Points.
   real(R8P)                                :: pu(NP), pv(NP), rr(2) !< Points.
   integer(I4P), allocatable                :: c0(:,:), c1(:,:), c2(:,:), s0(:), s1(:), s2(:), st(:,:) !< Results.
   real(R8P),    allocatable                :: x0(:,:), x1(:,:), x2(:,:) !< Results.
   integer(I4P)                             :: n, nj1, nf1, nj2, nf2 !< Counters.
   logical                                  :: ok    !< Agreement.

   do n = 1, NP
      call random_number(rr) ; pu(n) = (2 * rr(1) - 1) * half(1) ; pv(n) = (2 * rr(2) - 1) * half(2)
   enddo
   pu(1:10) = seam(1) + 0.01_R8P * [(n, n = 0, 9)] ; pv(1:10) = seam(2)
   call blk%locate_many(pu, pv, c0, x0, s0)
   call blk%locate_many(pu, pv, c1, x1, s1, walk=.true., jumps=nj1, fallbacks=nf1)
   allocate(st(2, NP))
   do n = 1, NP
      st(:, n) = max(1, min(c0(:, n) + [3, -2], [blk%ni - 1, blk%nj - 1]))
   enddo
   call blk%locate_many(pu, pv, c2, x2, s2, walk=.true., start=st, jumps=nj2, fallbacks=nf2)
   ok = all(c0 == c1) .and. all(x0 == x1) .and. all(s0 == s1) .and. all(c0 == c2) .and. all(x0 == x2) .and. all(s0 == s2)
   print '(2A,5(A,I0),A)', label, ':', ' points not inside ', count(s0 /= 1), ', jumps ', nj1, ', fallbacks ', nf1, &
         ' (from start cells: jumps ', nj2, ', fallbacks ', nf2, ')'
   call check(ok, label // ': walk = robust path, bit for bit')
   endsubroutine compare_2d

   subroutine compare_3d(blk, label, half, seam)
   !< Robust against walked location on random points of the box [-half, half] (z in [0, 2 half(3)]) and a seam point.
   type(curvilinear_block_t), intent(in) :: blk   !< Block.
   character(*),              intent(in) :: label !< Label.
   real(R8P),                 intent(in) :: half(3), seam(3) !< Box half sizes, a point on a seam.
   integer(I4P), parameter               :: NP = 1500 !< Points.
   real(R8P)                             :: pu(NP), pv(NP), pw(NP), rr(3) !< Points.
   integer(I4P), allocatable             :: c0(:,:), c1(:,:), c2(:,:), s0(:), s1(:), s2(:), st(:,:) !< Results.
   real(R8P),    allocatable             :: x0(:,:), x1(:,:), x2(:,:) !< Results.
   integer(I4P)                          :: n, nj1, nf1, nj2, nf2 !< Counters.
   logical                               :: ok    !< Agreement.

   do n = 1, NP
      call random_number(rr) ; pu(n) = (2 * rr(1) - 1) * half(1) ; pv(n) = (2 * rr(2) - 1) * half(2)
      pw(n) = 2 * rr(3) * half(3)
   enddo
   pu(1:10) = seam(1) + 0.01_R8P * [(n, n = 0, 9)] ; pv(1:10) = seam(2) ; pw(1:10) = seam(3)
   call blk%locate_many(pu, pv, pw, c0, x0, s0)
   call blk%locate_many(pu, pv, pw, c1, x1, s1, walk=.true., jumps=nj1, fallbacks=nf1)
   allocate(st(3, NP))
   do n = 1, NP
      st(:, n) = max(1, min(c0(:, n) + [2, -2, 1], [blk%ni - 1, blk%nj - 1, blk%nk - 1]))
   enddo
   call blk%locate_many(pu, pv, pw, c2, x2, s2, walk=.true., start=st, jumps=nj2, fallbacks=nf2)
   ok = all(c0 == c1) .and. all(x0 == x1) .and. all(s0 == s1) .and. all(c0 == c2) .and. all(x0 == x2) .and. all(s0 == s2)
   print '(2A,5(A,I0),A)', label, ':', ' points not inside ', count(s0 /= 1), ', jumps ', nj1, ', fallbacks ', nf1, &
         ' (from start cells: jumps ', nj2, ', fallbacks ', nf2, ')'
   call check(ok, label // ': walk = robust path, bit for bit')
   endsubroutine compare_3d

   subroutine check(cond, label)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.

   print '(A,L1,2A)', '  ', cond, '  ', label
   all_ok = all_ok .and. cond
   endsubroutine check
endprogram fossil_test_block_walk
