!< FOSSIL test: the load samples of a curve (#109 WP2 step 01).
program fossil_test_load_samples_2d
!< FOSSIL test: the load samples of a curve (#109 WP2 step 01).
!<
!< 1. Weights: on a disc with a hole (outer 64-gon counter-clockwise, hole 32-gon clockwise) the weights of every segment
!<    sum to its length, all of them to the perimeter, to round-off, for 1, 2, 3 Gauss points; normals outward (away from
!<    the centre on the outer loop, towards it on the hole), points on the segment, bodies, ids.
!< 2. Exactness: on one segment the rule of q points integrates every monomial u^a v^b of degree 2q - 1 to round-off and
!<    misses one of degree 2q.
!< 3. Visibility: two unit squares 0.25 apart; probes 0.1 out all visible, 0.5 out blocked exactly on the two facing
!<    sides, 1.5 in all blocked; an open polyline: the missing last segment is degenerate (zero weight and normal, probes
!<    not visible).
!< 4. Cells against a brute-force search; twin device = host bit for bit, posed too; posed = `to_world` of the unposed
!<    samples (to round-off, weights and visibility bit for bit) and = a copy moved by `move`; resident lifecycle;
!<    refusals; the empty curve.
use fossil, only : curve_object, load_samples_t, load_samples_2d, load_samples_reset, load_samples_mismatches, &
                   rigid_pose_2d_t, pose_2d, to_world, rotate, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter    :: PI = acos(-1._R8P), EPS = epsilon(1._R8P)
type(curve_object)      :: disc, squares, seg, open_line, moved, empty
type(load_samples_t)    :: sh, sd, sr, sc, sp
type(rigid_pose_2d_t)   :: pose
real(R8P), allocatable  :: u(:), v(:), x(:), y(:)
real(R8P)               :: e, per, ln, a(2), b(2), r
integer(I4P)            :: i, s, k, q, j, st, nbad, nblk
logical                 :: ok, miss

ok = .true.
u = [[(0.5_R8P + 0.4_R8P * cos(2 * PI * i / 64), i = 0, 63)], [(0.5_R8P + 0.15_R8P * cos(-2 * PI * i / 32), i = 0, 31)]]
v = [[(0.5_R8P + 0.4_R8P * sin(2 * PI * i / 64), i = 0, 63)], [(0.5_R8P + 0.15_R8P * sin(-2 * PI * i / 32), i = 0, 31)]]
call disc%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 32_I4P])
pose = pose_2d(0.3_R8P, [0.05_R8P, -0.02_R8P])
moved = disc ; call moved%move(pose)                    ! before enter_device: a host copy
u = [0._R8P, 1._R8P, 1._R8P, 0._R8P, 1.25_R8P, 2.25_R8P, 2.25_R8P, 1.25_R8P]
v = [0._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P]
call squares%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 4_I4P])
! 1. weights, normals, points
do q = 1, 3
   call load_samples_2d(disc, samples=sh, status=st, order=q)
   call check('disc order '//char(48 + q)//': status, q, n', st == STATUS_OK .and. sh%q == q .and. sh%n == 96 * q .and. &
              sh%ne == 96 .and. sh%dims == 2)
   e = 0._R8P ; per = 0._R8P ; nbad = 0
   do i = 1, 96
      a = disc%get_vertex(i) ; b = disc%get_vertex(disc%get_next(i)) ; ln = norm2(b - a) ; per = per + ln
      e = max(e, abs(sum(sh%weight((i - 1) * q + 1:i * q)) - ln) / ln)
      do s = (i - 1) * q + 1, i * q
         r = dot_product(sh%normal(:, s), sh%points(:, s) - 0.5_R8P)
         if ((i <= 64 .and. r <= 0._R8P) .or. (i > 64 .and. r >= 0._R8P)) nbad = nbad + 1
         if (abs(norm2(sh%normal(:, s)) - 1._R8P) > 4._R8P * EPS) nbad = nbad + 1
         if (abs(dot_product(sh%points(:, s) - a, sh%normal(:, s))) > 4._R8P * EPS) nbad = nbad + 1
         if (sh%element(s) /= i .or. sh%body(s) /= 1 .or. sh%patch(s) /= 1) nbad = nbad + 1
      enddo
   enddo
   call check('disc order '//char(48 + q)//': weights of a segment sum to its length (4 eps)', e <= 4._R8P * EPS)
   call check('disc order '//char(48 + q)//': all weights sum to the perimeter', abs(sum(sh%weight) - per) <= 1.e-14_R8P)
   call check('disc order '//char(48 + q)//': outward unit normals (the hole too), points on the segment, ids', nbad == 0)
enddo
! 2. exactness on one segment (an open polyline of two vertices)
call seg%adopt_loops(u=[0.3_R8P, 1.2_R8P], v=[-0.4_R8P, 0.7_R8P], loop_sizes=[2_I4P], closed=[.false.])
do q = 1, 3
   call load_samples_2d(seg, samples=sh, status=st, order=q)
   e = 0._R8P ; miss = .false.
   do i = 0, 2 * q
      do j = 0, 2 * q - i
         r = 0._R8P
         do s = 1, sh%n
            r = r + sh%weight(s) * sh%points(1, s)**i * sh%points(2, s)**j
         enddo
         if (i + j <= 2 * q - 1) then
            e = max(e, abs(r - seg_integral([0.3_R8P, -0.4_R8P], [1.2_R8P, 0.7_R8P], i, j)))
         elseif (abs(r - seg_integral([0.3_R8P, -0.4_R8P], [1.2_R8P, 0.7_R8P], i, j)) > 1.e-8_R8P) then
            miss = .true.
         endif
      enddo
   enddo
   call check('segment, '//char(48 + q)//' points: exact to degree '//char(48 + 2 * q - 1)//' (1e-14), not to the next', &
              e < 1.e-14_R8P .and. miss .and. all(sh%weight(q + 1:2 * q) == 0._R8P))
enddo
! 3. visibility blocked by a second body
call load_samples_2d(squares, samples=sh, status=st, order=3_I4P, lengths=[0.1_R8P, 0.5_R8P, -1.5_R8P, 0._R8P])
nbad = 0 ; nblk = 0
do s = 1, sh%n
   if (.not. sh%visible(1, s) .or. sh%visible(3, s)) nbad = nbad + 1
   if (.not. sh%visible(4, s) .or. any(sh%probe(:, 4, s) /= sh%points(:, s))) nbad = nbad + 1
   if (.not. sh%visible(2, s)) then
      nblk = nblk + 1
      if (.not. ((sh%element(s) == 2 .and. sh%normal(1, s) == 1._R8P) .or. (sh%element(s) == 8 .and. &
                 sh%normal(1, s) == -1._R8P))) nbad = nbad + 1
   endif
   if (any(sh%probe(:, 2, s) /= sh%points(:, s) + 0.5_R8P * sh%normal(:, s))) nbad = nbad + 1
enddo
call check('squares: 0.1 out visible, 1.5 in blocked, length 0 the point', nbad == 0 .and. all(sh%body(13:24) == 2))
call check('squares: 0.5 out blocked on the facing sides only (2 x 3 points)', nblk == 6)
call open_line%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P], v=[0._R8P, 0._R8P, 1._R8P], loop_sizes=[3_I4P], closed=[.false.])
call load_samples_2d(open_line, samples=sd, status=st, order=2_I4P, lengths=[0.1_R8P])
call check('open polyline: the missing last segment degenerate', st == STATUS_OK .and. sd%n == 6 .and. &
           all(sd%weight(5:6) == 0._R8P) .and. all(sd%normal(:, 5:6) == 0._R8P) .and. .not. any(sd%visible(1, 5:6)) .and. &
           all(sd%visible(1, 1:4)) .and. abs(sum(sd%weight) - 2._R8P) <= 4._R8P * EPS)
! 4. cells, twins, posed, resident
x = [(-0.35_R8P + 0.11_R8P * i + 0.003_R8P * i * i, i = 0, 25)] ; y = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call load_samples_2d(squares, x, y, sh, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.3_R8P])
call check('cells: lines, flags, brute force', st == STATUS_OK .and. sh%with_cells .and. sh%ni == 26 .and. sh%nk == 0 .and. &
           count_bad_cells(sh) == 0)
call squares%enter_device ; call disc%enter_device
call load_samples_2d(squares, x, y, sd, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.3_R8P], device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. load_samples_mismatches(sh, sd) == 0)
call load_samples_2d(disc, x, y, sh, status=st, order=3_I4P, lengths=[0.05_R8P, 0.3_R8P], pose=pose)
call load_samples_2d(disc, x, y, sd, status=st, order=3_I4P, lengths=[0.05_R8P, 0.3_R8P], pose=pose, device=.true.)
call check('twin posed: device = host, every field', st == STATUS_OK .and. load_samples_mismatches(sh, sd) == 0 .and. &
           count(.not. sh%visible) > 0)
call load_samples_2d(disc, samples=sp, status=st, order=3_I4P, lengths=[0.05_R8P, 0.3_R8P])
e = 0._R8P
do s = 1, sh%n
   e = max(e, norm2(sh%points(:, s) - to_world(pose, sp%points(:, s))), norm2(sh%normal(:, s) - rotate(pose, sp%normal(:, s))))
   do k = 1, 2
      e = max(e, norm2(sh%probe(:, k, s) - to_world(pose, sp%probe(:, k, s))))
   enddo
enddo
call check('posed: points, probes = to_world, normals rotated (4e-15)', e < 4.e-15_R8P)
call check('posed: weights and visibility = unposed, bit for bit; cells brute force', all(sh%weight == sp%weight) .and. &
           all(sh%visible .eqv. sp%visible) .and. count_bad_cells(sh) == 0)
call load_samples_2d(moved, samples=sc, status=st, order=3_I4P, lengths=[0.05_R8P, 0.3_R8P])
call check('posed: = a moved copy (4e-15), the same visibility', maxval(abs(sc%points - sh%points)) < 4.e-15_R8P .and. &
           maxval(abs(sc%probe - sh%probe)) < 4.e-15_R8P .and. all(sc%visible .eqv. sh%visible))
call load_samples_2d(squares, x, y, sh, status=st, lengths=[0.1_R8P, 0.5_R8P])
call load_samples_2d(squares, x, y, sr, status=st, lengths=[0.1_R8P, 0.5_R8P], resident=.true.)
call check('resident: flagged', st == STATUS_OK .and. sr%on_device)
call sr%update_host
call check('resident: update_host gives the host bits', load_samples_mismatches(sh, sr) == 0 .and. sr%on_device)
sc = sr
call sc%update_host                                     ! a copy: warns on stderr, drops the flag
call check('resident: a copy is not mapped', .not. sc%on_device .and. load_samples_mismatches(sh, sc) == 0)
call sr%exit_device ; call sr%enter_device
call check('resident: exit_device, enter_device', sr%on_device)
call load_samples_2d(squares, samples=sr, status=st, device=.true.)
call check('resident: a new call releases it first', st == STATUS_OK .and. .not. sr%on_device .and. sr%n == 16)
call load_samples_2d(squares, samples=sr, status=st, resident=.true.)
call load_samples_reset(sr)
call check('resident: load_samples_reset', .not. sr%on_device .and. sr%n == 0 .and. .not. allocated(sr%points))
! refusals and the empty curve
call load_samples_2d(squares, samples=sd, status=st, order=4_I4P)
call check('refused: order 4', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_2d(squares, x, samples=sd, status=st)
call check('refused: lines in part', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_2d(squares, x, [1._R8P, 0._R8P], sd, status=st)
call check('refused: a line not increasing', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call squares%exit_device
call load_samples_2d(squares, samples=sd, status=st, device=.true.)
call check('refused: device path of an unmapped curve (stderr message)', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_2d(empty, x, y, sd, status=st, lengths=[0.1_R8P], device=.true.)
call check('empty curve: no sample, both paths', st == STATUS_OK .and. sd%n == 0 .and. allocated(sd%points))
call disc%exit_device

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (.not. cond) then
      print '(A)', 'FAILED: '//label
      ok = .false.
   else
      print '(A)', 'ok: '//label
   endif
   endsubroutine check

   real(R8P) function seg_integral(p0, p1, ia, ib) result(rr)
   !< Exact integral of u^ia v^ib along the segment p0-p1: L int_0^1 (u0 + s du)^ia (v0 + s dv)^ib ds, the integrand
   !< expanded in powers of s.
   real(R8P),    intent(in) :: p0(2), p1(2) !< Ends.
   integer(I4P), intent(in) :: ia, ib       !< Exponents.
   real(R8P)    :: c(0:8), cn(0:8)          !< Polynomial in s.
   integer(I4P) :: d, jj, kk                !< Counters.

   c = 0._R8P ; c(0) = 1._R8P
   do d = 1, 2
      do jj = 1, merge(ia, ib, d == 1)
         cn = 0._R8P
         do kk = 0, 7
            cn(kk) = cn(kk) + p0(d) * c(kk)
            cn(kk + 1) = cn(kk + 1) + (p1(d) - p0(d)) * c(kk)
         enddo
         c = cn
      enddo
   enddo
   rr = 0._R8P
   do kk = 0, 8
      rr = rr + c(kk) / (kk + 1)
   enddo
   rr = rr * norm2(p1 - p0)
   endfunction seg_integral

   integer(I4P) function count_bad_cells(smp) result(nb)
   !< Probes whose cells differ from a brute-force search on x, y (low-side rule, 0 outside).
   type(load_samples_t), intent(in) :: smp !< Samples.
   integer(I4P) :: s, k                    !< Counters.

   nb = 0
   do s = 1, smp%n
      do k = 1, smp%m
         if (smp%probe_cell(1, k, s) /= brute(x, smp%probe(1, k, s))) nb = nb + 1
         if (smp%probe_cell(2, k, s) /= brute(y, smp%probe(2, k, s))) nb = nb + 1
      enddo
   enddo
   endfunction count_bad_cells

   pure integer(I4P) function brute(e, c) result(idx)
   !< The cell i with e(i) < c <= e(i+1) by a linear scan, 0 if none.
   real(R8P), intent(in) :: e(:) !< Lines.
   real(R8P), intent(in) :: c    !< Coordinate.
   integer(I4P) :: i             !< Counter.

   idx = 0
   do i = 1, size(e) - 1
      if (e(i) < c .and. c <= e(i + 1)) idx = i
   enddo
   endfunction brute
endprogram fossil_test_load_samples_2d
