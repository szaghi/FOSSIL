!< FOSSIL test: the 2-D swept geometry of several rigid bodies, Gauss path (#109 WP3 step 04).
program fossil_test_swept_bodies_2d
!< FOSSIL test: the 2-D swept geometry of several rigid bodies, Gauss path (#109 WP3 step 04).
!<
!< `swept_cells_2d_bodies` (`fossil_swept_cells`) over polygonal circles and squares (nq = 3):
!< 1. one body: equal to `swept_cells_2d(..., exact=.false.)` bit for bit (`swept_cells_2d_mismatches`), a 48-gon
!<    turning while translating;
!< 2. two far-apart bodies on their own motions: every cell and face on the first body's side equals the first body
!<    alone, every one on the second's side the second alone, bit for bit, but the residual to round-off (1e-16): at a
!<    node where a cell holds no wall its balance is that of body 1's field (zero to round-off);
!< 3. the swept area is A1 - A0 of the cut of the union at the end poses (`cut_cells_2d_bodies`), bit for bit, and the
!<    swept areas add up to zero over a block holding both bodies;
!< 4. two squares 0.2 apart on one rigid motion (one body frame, the same poses), with cells holding both walls (a
!<    piece bounded by both): faces, alpha and swept equal those of the merged curve on that motion bit for bit, the
!<    residual (the wall form in those cells) the merged one's (the face form) to round-off;
!< 5. refusals: exact absent or .true., poses not matching, bodies overlapping at the end pose.
use fossil,               only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, &
                                 swept_cells_2d_bodies, cut_cells_2d_t, cut_cells_2d_bodies, merge_posed_bodies, &
                                 piece_bodies_t, STATUS_OK, STATUS_INVALID_INPUT
use fossil_swept_compare, only : swept_cells_2d_mismatches
use penf,                 only : I4P, I8P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P)
type(curve_object)    :: circ, two(2)
type(rigid_pose_2d_t) :: p0(2), p1(2)
logical               :: ok

ok = .true.
call polygon(circ, 0.5_R8P, 48_I4P)
two(1) = circ ; two(2) = circ
call one_body_case
call far_apart_case
call pair_case
call refusal_case
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine one_body_case
   !< Case 1: one body equals swept_cells_2d on the Gauss path.
   type(swept_cells_2d_t) :: sa, sb
   real(R8P), allocatable :: x(:), y(:)
   integer(I4P)           :: i, st1, st2, m

   x = [(-0.83_R8P + 0.15_R8P * i, i = 0, 12)] ; y = [(-0.77_R8P + 0.15_R8P * i, i = 0, 10)]
   p0(1) = pose_2d(0.1_R8P, [0.02_R8P, -0.03_R8P]) ; p1(1) = pose_2d(0.45_R8P, [0.19_R8P, 0.08_R8P])
   call swept_cells_2d(circ, p0(1), p1(1), x, y, sa, nq=3_I4P, exact=.false., status=st1)
   call swept_cells_2d_bodies([circ], p0(1:1), p1(1:1), x, y, sb, nq=3_I4P, exact=.false., status=st2)
   m = abs(st1) + abs(st2) + swept_cells_2d_mismatches(sa, sb)
   print '(A,I0,A,ES9.2)', 'one body vs swept_cells_2d: mismatches ', m, ', max |residual| ', maxval(abs(sb%residual))
   ok = ok .and. m == 0_I4P
   endsubroutine one_body_case

   subroutine far_apart_case
   !< Cases 2 and 3: two far-apart bodies on their own motions; the swept areas.
   type(swept_cells_2d_t) :: sw, s1, s2
   type(cut_cells_2d_t)   :: k0, k1
   real(R8P), allocatable :: x(:), y(:)
   real(R8P)              :: total
   integer(I4P)           :: i, j, st, m, nh, nv

   x = [(-0.83_R8P + 0.15_R8P * i, i = 0, 24)] ; y = [(-0.77_R8P + 0.15_R8P * i, i = 0, 10)]
   nh = 12_I4P ! cells 1..12 (u < 0.97): the first circle's side; 13..24: the second's
   p0(1) = pose_2d(0.1_R8P, [0.02_R8P, -0.03_R8P]) ; p1(1) = pose_2d(0.45_R8P, [0.19_R8P, 0.08_R8P])
   p0(2) = pose_2d(-0.2_R8P, [2.1_R8P, 0.01_R8P]) ; p1(2) = pose_2d(-0.6_R8P, [1.93_R8P, -0.06_R8P])
   call swept_cells_2d_bodies(two, p0, p1, x, y, sw, nq=3_I4P, exact=.false., status=st)
   call swept_cells_2d(circ, p0(1), p1(1), x, y, s1, nq=3_I4P, exact=.false.)
   call swept_cells_2d(circ, p0(2), p1(2), x, y, s2, nq=3_I4P, exact=.false.)
   m = abs(st)
   do j = 1, size(y) - 1
      do i = 1, size(x) - 1
         if (i <= nh) then
            m = m + cell_diff(sw, s1, i, j)
         else
            m = m + cell_diff(sw, s2, i, j)
         endif
      enddo
   enddo
   nv = count(sw%swept /= 0._R8P)
   print '(A,I0,A,I0,A)', 'two far-apart bodies vs each body alone: mismatches ', m, ' (', nv, ' cells with a swept area)'
   ok = ok .and. m == 0_I4P .and. nv > 0_I4P
   call cut_cells_2d_bodies(two, p0, x, y, k0, status=st) ; m = abs(st)
   call cut_cells_2d_bodies(two, p1, x, y, k1, status=st) ; m = m + abs(st)
   do j = 1, size(y) - 1
      do i = 1, size(x) - 1
         if (bits(sw%swept(i, j)) /= bits((k1%alpha(i, j) - k0%alpha(i, j)) * (x(i + 1) - x(i)) * (y(j + 1) - y(j)))) &
            m = m + 1_I4P
      enddo
   enddo
   total = sum(sw%swept)
   print '(A,I0,A,ES9.2)', 'swept = A1 - A0 of the union: mismatches ', m, ', |sum of swept| ', abs(total)
   ok = ok .and. m == 0_I4P .and. abs(total) < 1.e-14_R8P
   endsubroutine far_apart_case

   subroutine pair_case
   !< Case 4: two squares 0.2 apart on one rigid motion vs their merged curve.
   type(curve_object)        :: merged, pair(2)
   type(swept_cells_2d_t)    :: sw, sm
   type(cut_cells_2d_t)      :: kc
   type(piece_bodies_t)      :: pb
   type(rigid_pose_2d_t)     :: q(2)
   integer(I4P), allocatable :: sb(:), ss(:)
   real(R8P),    allocatable :: x(:), y(:)
   real(R8P)                 :: dres
   integer(I4P)              :: i, st, m, nmix, p

   call pair(1)%adopt_loops(u=[0._R8P, 0.5_R8P, 0.5_R8P, 0._R8P], v=[0._R8P, 0._R8P, 0.5_R8P, 0.5_R8P], loop_sizes=[4_I4P])
   call pair(2)%adopt_loops(u=[0.7_R8P, 1.2_R8P, 1.2_R8P, 0.7_R8P], v=[0._R8P, 0._R8P, 0.5_R8P, 0.5_R8P], loop_sizes=[4_I4P])
   x = [(-0.45_R8P + 0.3_R8P * i, i = 0, 6)] ; y = [(-0.31_R8P + 0.3_R8P * i, i = 0, 4)]
   q = pose_2d(0._R8P)
   call merge_posed_bodies(pair, q, merged, sb, ss, status=st) ! the pair as one body
   p0 = pose_2d(0._R8P) ; p1 = pose_2d(0.12_R8P, [0.05_R8P, 0.04_R8P])
   call swept_cells_2d(merged, p0(1), p1(1), x, y, sm, nq=3_I4P, exact=.false.)
   call swept_cells_2d_bodies(pair, p0, p1, x, y, sw, nq=3_I4P, exact=.false., status=st)
   m = abs(st)
   call cut_cells_2d_bodies(pair, p0, x, y, kc, status=st, piece_bodies=pb)
   nmix = 0_I4P
   do p = 1, pb%np
      if (pb%start(p + 1) - pb%start(p) >= 2_I4P) nmix = nmix + 1_I4P
   enddo
   m = m + nbits(sw%alpha0, sm%alpha0) + nbits(sw%alpha1, sm%alpha1) + nbits(sw%swept, sm%swept) + &
       nbits(sw%face_u, sm%face_u) + nbits(sw%face_v, sm%face_v)
   dres = maxval(abs(sw%residual - sm%residual))
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2)', 'two squares 0.2 apart on one motion vs the merged curve: mismatches ', m, &
         ', pieces bounded by both ', nmix, ', max |residual difference| ', dres, ', max |residual| ', maxval(abs(sm%residual))
   ok = ok .and. m == 0_I4P .and. nmix > 0_I4P .and. dres < 1.e-16_R8P
   endsubroutine pair_case

   subroutine refusal_case
   !< Case 5: refusals.
   type(swept_cells_2d_t) :: sw
   real(R8P), allocatable :: x(:), y(:)
   integer(I4P)           :: i, st, m

   x = [(-0.83_R8P + 0.15_R8P * i, i = 0, 16)] ; y = [(-0.77_R8P + 0.15_R8P * i, i = 0, 10)]
   p0(1) = pose_2d(0._R8P) ; p1(1) = p0(1)
   p0(2) = pose_2d(0._R8P, [1.2_R8P, 0._R8P]) ; p1(2) = pose_2d(0._R8P, [0.8_R8P, 0._R8P])
   m = 0_I4P
   call swept_cells_2d_bodies(two, p0, p0, x, y, sw, nq=3_I4P, status=st) ! exact absent
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_2d_bodies(two, p0, p0, x, y, sw, nq=3_I4P, exact=.true., status=st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_2d_bodies(two, p0(1:1), p0, x, y, sw, nq=3_I4P, exact=.false., status=st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_2d_bodies(two, p0, p1, x, y, sw, nq=3_I4P, exact=.false., status=st) ! overlap at the end
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   print '(A,I0)', 'refusals (exact absent, exact, poses, overlap): mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine refusal_case

   function cell_diff(a, b, i, j) result(n)
   !< Differing values of cell (i, j) of a and b: the cell arrays and its four faces; the residual to round-off.
   type(swept_cells_2d_t), intent(in) :: a, b !< Swept geometries.
   integer(I4P),           intent(in) :: i, j !< Cell.
   integer(I4P)                       :: n    !< Mismatches.

   n = 0_I4P
   if (bits(a%alpha0(i, j)) /= bits(b%alpha0(i, j))) n = n + 1_I4P
   if (bits(a%alpha1(i, j)) /= bits(b%alpha1(i, j))) n = n + 1_I4P
   if (bits(a%swept(i, j)) /= bits(b%swept(i, j))) n = n + 1_I4P
   if (abs(a%residual(i, j) - b%residual(i, j)) > 1.e-16_R8P) n = n + 1_I4P ! round-off: see the case notes
   if (any(bits(a%wall(:, i, j)) /= bits(b%wall(:, i, j)))) n = n + 1_I4P
   if (bits(a%face_u(i, j)) /= bits(b%face_u(i, j))) n = n + 1_I4P
   if (bits(a%face_u(i + 1, j)) /= bits(b%face_u(i + 1, j))) n = n + 1_I4P
   if (bits(a%face_v(i, j)) /= bits(b%face_v(i, j))) n = n + 1_I4P
   if (bits(a%face_v(i, j + 1)) /= bits(b%face_v(i, j + 1))) n = n + 1_I4P
   if (a%cleared(i, j) .neqv. b%cleared(i, j)) n = n + 1_I4P
   if (a%covered(i, j) .neqv. b%covered(i, j)) n = n + 1_I4P
   endfunction cell_diff

   function nbits(a, b) result(n)
   !< Differing values of two real arrays of one shape (bit patterns).
   real(R8P), intent(in) :: a(:,:), b(:,:) !< Arrays.
   integer(I4P)          :: n              !< Mismatches.

   n = count(transfer(a, 1_I8P, size(a)) /= transfer(b, 1_I8P, size(b)))
   endfunction nbits

   elemental function bits(v) result(b)
   !< Bit pattern of a real.
   real(R8P), intent(in) :: v !< Value.
   integer(I8P)          :: b !< Bits.

   b = transfer(v, 1_I8P)
   endfunction bits

   subroutine polygon(c, r, n)
   !< Regular n-gon of radius r about the origin, counter-clockwise.
   type(curve_object), intent(inout) :: c !< Curve.
   real(R8P),          intent(in)    :: r !< Radius.
   integer(I4P),       intent(in)    :: n !< Vertices.
   real(R8P)                         :: u(n), v(n) !< Vertices.
   integer(I4P)                      :: k !< Work.

   do k = 1, n
      u(k) = r * cos(2._R8P * PI * (k - 1) / n) ; v(k) = r * sin(2._R8P * PI * (k - 1) / n)
   enddo
   call c%adopt_loops(u=u, v=v, loop_sizes=[n])
   endsubroutine polygon
endprogram fossil_test_swept_bodies_2d
