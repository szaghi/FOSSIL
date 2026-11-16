!< FOSSIL test: ghost points of a finite-difference wall closure on a 2-D block (#109 WP1 step 03).
program fossil_test_ghost_points_2d
!< FOSSIL test: ghost points of a finite-difference wall closure on a 2-D block (#109 WP1 step 03).
!<
!< 1. Hand case, the low-side rule: the unit square on lines of step 1/4 through its sides (u, v = -1/2 .. 3/2). The
!<    perturbed points put u = 0 inside and u = 1 outside, so the solid points are {0, 1/4, 1/2, 3/4}^2 (16): with reach 1
!<    the 12 on the boundary of that block are ghosts, with reach 2 all 16; masks and steps of chosen points.
!< 2. The point flags equal the 2-D cut cells' `node_in` bit for bit on the same lines (square, circle, NACA 0012).
!< 3. Brute force: a circle and a NACA 0012 section on generic lines (uniform and stretched, mixed reach), the flags from
!<    `is_point_inside` (no point on the wall) and the ghosts from a direct loop over all points, against
!<    `ghost_points_2d`: count, indexes, coordinates, masks, steps.
!< 4. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path bit for bit on every
!<    case; an empty curve has no ghost on either path.
!< 5. Bad input: one point on an axis, lines not increasing, reach 0, size(reach) /= 2: STATUS_INVALID_INPUT, n = 0.
use fossil, only : curve_object, cut_cells_2d_t, ghost_points_t, ghost_points_2d, STATUS_OK, STATUS_INVALID_INPUT, &
                   GHOST_XM, GHOST_XP, GHOST_YM, GHOST_YP
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P)
type(curve_object)    :: square, circle, naca, empty
type(ghost_points_t)  :: gh
type(cut_cells_2d_t)  :: cut
real(R8P)             :: u(64), v(64), nu(79), nv(79), xc
real(R8P), allocatable :: x(:), y(:)
integer(I4P)          :: i, st, g
logical               :: ok

ok = .true.
! 1. the unit square on lines through its sides
call square%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call ghost_points_2d(square, x, x, [1_I4P, 1_I4P], gh, status=st)
call check('square reach 1: 12 ghosts', st == STATUS_OK .and. gh%n == 12_I4P .and. gh%dims == 2_I4P .and. gh%nk == 1_I4P)
call check('square reach 1: every ghost in {0, 1/4, 1/2, 3/4}^2', all(gh%points >= 0._R8P) .and. all(gh%points <= 0.75_R8P))
g = find(gh, [3_I4P, 3_I4P])                           ! (0, 0): on two sides, inside by the rule
call check('square (0,0): ghost', g > 0)
if (g > 0) call check('square (0,0): mask and steps', gh%mask(g) == ior(GHOST_XM, GHOST_YM) .and. &
                       all(gh%steps(:, g) == [1_I4P, 0_I4P, 1_I4P, 0_I4P]))
g = find(gh, [6_I4P, 6_I4P])                           ! (3/4, 3/4): u = 1 and v = 1 are outside by the rule
call check('square (3/4,3/4): ghost', g > 0)
if (g > 0) call check('square (3/4,3/4): mask and steps', gh%mask(g) == ior(GHOST_XP, GHOST_YP) .and. &
                       all(gh%steps(:, g) == [0_I4P, 1_I4P, 0_I4P, 1_I4P]))
call check('square (1,1/2) is not a ghost (fluid)', find(gh, [7_I4P, 5_I4P]) == 0)
call check('square: ghosts in lexicographic order', ordered(gh))
call twin('square reach 1', square, x, x, [1_I4P, 1_I4P], gh)
call nodes('square', square, x, x)
call ghost_points_2d(square, x, x, [2_I4P, 2_I4P], gh, status=st)
call check('square reach 2: 16 ghosts', st == STATUS_OK .and. gh%n == 16_I4P)
g = find(gh, [4_I4P, 5_I4P])                           ! (1/4, 1/2): fluid 2 steps to -u, 3 to +u, 3 to -v, 2 to +v
if (g > 0) call check('square (1/4,1/2) reach 2: steps', all(gh%steps(:, g) == [2_I4P, 0_I4P, 0_I4P, 2_I4P]))
call check('square (1/4,1/2) reach 2: ghost', g > 0)
call twin('square reach 2', square, x, x, [2_I4P, 2_I4P], gh)
call ghost_points_2d(square, x, x, [1_I4P, 4_I4P], gh, status=st)
call check('square reach (1, 4): 16 ghosts', st == STATUS_OK .and. gh%n == 16_I4P)
call twin('square reach (1, 4)', square, x, x, [1_I4P, 4_I4P], gh)
! 2-3. a circle and a NACA 0012 section, against a brute force
u = [(0.03_R8P + 0.61_R8P * cos(2 * PI * i / 64), i = 0, 63)] ; v = [(-0.02_R8P + 0.61_R8P * sin(2 * PI * i / 64), i = 0, 63)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
x = [(-0.81_R8P + 0.0531_R8P * i, i = 0, 31)] ; y = [(-0.79_R8P + 0.0547_R8P * i, i = 0, 29)]
call brute('circle reach 1', circle, x, y, [1_I4P, 1_I4P])
call brute('circle reach (3, 2)', circle, x, y, [3_I4P, 2_I4P])
call nodes('circle', circle, x, y)
do i = 0, 40                                           ! NACA 0012, counter-clockwise: upper side from the trailing edge
   xc = 0.5_R8P * (1._R8P + cos(PI * i / 40)) ; nu(i + 1) = xc ; nv(i + 1) = thickness(xc)
enddo
do i = 1, 38                                           ! lower side from the nose
   xc = 0.5_R8P * (1._R8P - cos(PI * i / 39)) ; nu(41 + i) = xc ; nv(41 + i) = -thickness(xc)
enddo
call naca%adopt_loops(u=nu, v=nv, loop_sizes=[79_I4P])
x = [(-0.103_R8P + 0.0377_R8P * i, i = 0, 31)] ; y = [(-0.0917_R8P + 0.0119_R8P * i, i = 0, 15)]
call brute('naca reach 1', naca, x, y, [1_I4P, 1_I4P])
call brute('naca reach (2, 3)', naca, x, y, [2_I4P, 3_I4P])
x = [(-0.103_R8P + 1.21_R8P * (real(i, R8P) / 31._R8P)**1.6_R8P, i = 0, 31)]  ! stretched, fine at the nose
call brute('naca stretched, reach 2', naca, x, y, [2_I4P, 2_I4P])
call nodes('naca', naca, x, y)
! 4. an empty curve
call ghost_points_2d(empty, x, y, [1_I4P, 1_I4P], gh, status=st)
call check('empty curve: no ghost', st == STATUS_OK .and. gh%n == 0_I4P .and. allocated(gh%ijk))
call twin('empty curve', empty, x, y, [1_I4P, 1_I4P], gh)
! 5. bad input
call ghost_points_2d(square, [0.5_R8P], y, [1_I4P, 1_I4P], gh, status=st)
call check('one u point: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P .and. .not. allocated(gh%ijk))
call ghost_points_2d(square, x, [0._R8P, 0._R8P, 1._R8P], [1_I4P, 1_I4P], gh, status=st, device=.true.)
call check('repeated v: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
call ghost_points_2d(square, x, y, [0_I4P, 1_I4P], gh, status=st)
call check('reach 0: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
call ghost_points_2d(square, x, y, [1_I4P, 1_I4P, 1_I4P], gh, status=st)
call check('three reaches: refused', st == STATUS_INVALID_INPUT .and. gh%n == 0_I4P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label, pass)
   !< Record one check.
   character(*), intent(in) :: label !< Check.
   logical,      intent(in) :: pass  !< Result.

   if (.not. pass) then
      ok = .false. ; print '(A)', 'FAILED: '//label
   endif
   endsubroutine check

   pure real(R8P) function thickness(t)
   !< NACA 0012 half thickness at t (closed trailing edge).
   real(R8P), intent(in) :: t !< Chord fraction.

   thickness = 0.6_R8P * (0.2969_R8P * sqrt(t) - 0.1260_R8P * t - 0.3516_R8P * t**2 + 0.2843_R8P * t**3 - 0.1036_R8P * t**4)
   endfunction thickness

   integer(I4P) function find(gp, ij) result(g)
   !< Slot of the point ij among the ghosts, 0 if not a ghost.
   type(ghost_points_t), intent(in) :: gp    !< Ghost points.
   integer(I4P),         intent(in) :: ij(2) !< Point.
   integer(I4P)                     :: q     !< Work.

   g = 0_I4P
   do q = 1, gp%n
      if (all(gp%ijk(:, q) == ij)) then
         g = q ; return
      endif
   enddo
   endfunction find

   logical function ordered(gp)
   !< The ghosts in lexicographic order, i fastest.
   type(ghost_points_t), intent(in) :: gp !< Ghost points.
   integer(I4P)                     :: q  !< Work.

   ordered = .true.
   do q = 2, gp%n
      if (gp%ijk(1, q - 1) + gp%ni * gp%ijk(2, q - 1) >= gp%ijk(1, q) + gp%ni * gp%ijk(2, q)) ordered = .false.
   enddo
   endfunction ordered

   integer(I4P) function mismatches(a, b) result(m)
   !< Differences between two ghost sets (coordinates by their bits).
   type(ghost_points_t), intent(in) :: a, b !< Ghost points.

   m = 0_I4P
   if (a%n /= b%n .or. a%dims /= b%dims .or. a%ni /= b%ni .or. a%nj /= b%nj .or. a%nk /= b%nk .or. &
       any(a%reach /= b%reach)) then
      m = 1_I4P ; return
   endif
   m = count(a%ijk /= b%ijk) + count(a%mask /= b%mask) + count(a%steps /= b%steps) + &
       count(transfer(a%points, 1_I8P, size(a%points)) /= transfer(b%points, 1_I8P, size(b%points)))
   endfunction mismatches

   subroutine twin(label, c, x, y, reach, host)
   !< The device path against the host result.
   character(*),         intent(in) :: label      !< Case.
   type(curve_object),   intent(in) :: c          !< Curve.
   real(R8P),            intent(in) :: x(:), y(:) !< Lines.
   integer(I4P),         intent(in) :: reach(2)   !< Reach.
   type(ghost_points_t), intent(in) :: host       !< Host result.
   type(ghost_points_t)             :: dev        !< Device result.
   integer(I4P)                     :: st2, m     !< Work.

   call ghost_points_2d(c, x, y, reach, dev, status=st2, device=.true.)
   m = mismatches(host, dev)
   call check(label//': device twin', st2 == STATUS_OK .and. m == 0_I4P)
   if (m /= 0_I4P) print '(A,I0)', '   mismatches ', m
   endsubroutine twin

   subroutine nodes(label, c, x, y)
   !< The solid points of the ghost search are the cut cells' node states: every ghost is a solid node, and with a
   !< reach covering the block every solid node with a fluid node on its lines is a ghost.
   character(*),       intent(in) :: label      !< Case.
   type(curve_object), intent(in) :: c          !< Curve.
   real(R8P),          intent(in) :: x(:), y(:) !< Lines.
   type(ghost_points_t)           :: gp         !< Ghosts.
   logical, allocatable           :: solid(:,:) !< Solid flags from the ghosts.
   integer(I4P)                   :: q, i, j    !< Work.
   logical                        :: line       !< A fluid node on the lines of (i, j).

   call c%cut_cells(x=x, y=y, cut=cut)
   call ghost_points_2d(c, x, y, [size(x, kind=I4P), size(y, kind=I4P)], gp)
   allocate(solid(size(x), size(y))) ; solid = .false.
   do q = 1, gp%n
      solid(gp%ijk(1, q), gp%ijk(2, q)) = .true.
   enddo
   call check(label//': ghosts are solid nodes', all(cut%node_in .or. .not. solid))
   do j = 1, size(y) ; do i = 1, size(x)
      line = .not. (all(cut%node_in(:, j)) .and. all(cut%node_in(i, :)))
      if (cut%node_in(i, j) .and. line .and. .not. solid(i, j)) then
         call check(label//': solid node with a fluid line is a ghost', .false.) ; return
      endif
   enddo ; enddo
   endsubroutine nodes

   subroutine brute(label, c, x, y, reach)
   !< `ghost_points_2d` against a direct loop over all points with flags from `is_point_inside`; then the device twin.
   character(*),       intent(in) :: label      !< Case.
   type(curve_object), intent(in) :: c          !< Curve.
   real(R8P),          intent(in) :: x(:), y(:) !< Lines.
   integer(I4P),       intent(in) :: reach(2)   !< Reach.
   type(ghost_points_t)           :: gp         !< Module result.
   logical,      allocatable      :: ins(:,:)   !< Brute-force flags.
   integer(I4P)                   :: off(2, 4), q(2), stp(4), msk, n, i, j, d, t, st2, m !< Work.

   allocate(ins(size(x), size(y)))
   do j = 1, size(y) ; do i = 1, size(x)
      ins(i, j) = c%is_point_inside([x(i), y(j)])
   enddo ; enddo
   off = reshape([-1, 0, 1, 0, 0, -1, 0, 1], [2, 4])
   call ghost_points_2d(c, x, y, reach, gp, status=st2)
   call check(label//': status', st2 == STATUS_OK)
   n = 0_I4P ; m = 0_I4P
   do j = 1, size(y) ; do i = 1, size(x)
      if (.not. ins(i, j)) cycle
      msk = 0_I4P ; stp = 0_I4P
      do d = 1, 4
         do t = 1, reach((d + 1) / 2)
            q = [i, j] + t * off(:, d)
            if (any(q < 1) .or. q(1) > size(x) .or. q(2) > size(y)) exit
            if (.not. ins(q(1), q(2))) then
               stp(d) = t ; msk = msk + 2_I4P**(d - 1) ; exit
            endif
         enddo
      enddo
      if (msk == 0_I4P) cycle
      n = n + 1_I4P
      if (n > gp%n) cycle
      if (any(gp%ijk(:, n) /= [i, j]) .or. gp%mask(n) /= msk .or. any(gp%steps(:, n) /= stp) .or. &
          any(gp%points(:, n) /= [x(i), y(j)])) m = m + 1_I4P
   enddo ; enddo
   call check(label//': count = brute force', n == gp%n .and. n > 0_I4P)
   call check(label//': ghosts = brute force', m == 0_I4P)
   if (n /= gp%n .or. m /= 0_I4P) print '(A,3(I0,1X))', '   brute, module, mismatches ', n, gp%n, m
   call twin(label, c, x, y, reach, gp)
   endsubroutine brute
endprogram fossil_test_ghost_points_2d
