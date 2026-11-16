!< FOSSIL test: the 3-D swept geometry of several rigid bodies, Gauss path (#109 WP3 step 04).
program fossil_test_swept_bodies
!< FOSSIL test: the 3-D swept geometry of several rigid bodies, Gauss path (#109 WP3 step 04).
!<
!< `swept_cells_3d_bodies` (`fossil_swept_cells`) over unit cubes (`bench_bodies`, nq = 2):
!< 1. one body: equal to `swept_cells_3d(..., exact=.false.)` bit for bit (`swept_cells_3d_mismatches`), turning about
!<    a tilted axis while translating;
!< 2. two far-apart bodies on their own motions: every cell and face on the first body's side equals the first body
!<    alone, every one on the second's side the second alone, bit for bit, but the residual to round-off (1e-16): at a
!<    node where a cell holds no wall its balance is that of body 1's field (zero to round-off), where the second body
!<    alone uses its own;
!< 3. the swept volume is V1 - V0 of the cut of the union at the end poses (`cut_cells_3d_bodies`), bit for bit, and the
!<    swept volumes add up to zero over a block holding both bodies;
!< 4. two cubes 0.2 apart, on one rigid motion (one body frame, the same poses), with cells holding both walls (a piece
!<    bounded by both): the faces,
!<    alpha and swept equal those of the merged surface on that motion bit for bit, the residual (the wall form in
!<    those cells) equals the merged one's (the face form) to round-off;
!< 5. refusals: exact absent or .true., poses not matching, cubes overlapping at the end pose.
use fossil,                   only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, &
                                     swept_cells_3d_t, swept_cells_3d, swept_cells_3d_bodies, cut_cells_3d_t, &
                                     cut_cells_3d_bodies, merge_posed_bodies, piece_bodies_t, STATUS_OK, STATUS_INVALID_INPUT
use fossil_swept_3d_compare,  only : swept_cells_3d_mismatches
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P

implicit none

type(body_t), allocatable :: bb(:)
type(surface_stl_object)  :: cube, two(2)
type(rigid_pose_t)        :: p0(2), p1(2)
logical                   :: ok

ok = .true.
call bench_bodies(bb)
call to_surface(bb(1)%soup, cube) ! unit cube [0, 1]^3
two(1) = cube ; two(2) = cube
call one_body_case
call far_apart_case
call pair_case
call refusal_case
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine one_body_case
   !< Case 1: one body equals swept_cells_3d on the Gauss path.
   type(swept_cells_3d_t)   :: sa, sb
   type(surface_stl_object) :: one(1) ! an array, not a constructor (f08/0011 finalisation)
   real(R8P)                :: e(5)
   integer(I4P)             :: i, st1, st2, m

   e = [(-0.23_R8P + 0.45_R8P * i, i = 0, 4)]
   p0(1) = pose_identity()
   p1(1) = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P, [0.17_R8P, -0.11_R8P, 0.09_R8P])
   call swept_cells_3d(cube, p0(1), p1(1), e, e, e, sa, nq=2_I4P, exact=.false., status=st1)
   one(1) = cube
   call swept_cells_3d_bodies(one, p0(1:1), p1(1:1), e, e, e, sb, nq=2_I4P, exact=.false., status=st2)
   m = abs(st1) + abs(st2) + swept_cells_3d_mismatches(sa, sb)
   print '(A,I0,A,ES9.2)', 'one body vs swept_cells_3d: mismatches ', m, ', max |residual| ', maxval(abs(sb%residual))
   ok = ok .and. m == 0_I4P
   endsubroutine one_body_case

   subroutine far_apart_case
   !< Cases 2 and 3: two far-apart bodies on their own motions; the swept volumes.
   type(swept_cells_3d_t) :: sw, s1, s2
   type(cut_cells_3d_t)   :: k0, k1
   real(R8P), allocatable :: x(:), e(:)
   real(R8P)              :: total
   integer(I4P)           :: i, j, k, st, m, nh, nv

   x = [(-0.55_R8P + 0.45_R8P * i, i = 0, 12)] ; e = [(-0.55_R8P + 0.45_R8P * i, i = 0, 5)]
   nh = 6_I4P ! cells 1..6 (x < 2.15): the first cube's side; 7..12: the second's
   p0(1) = pose_identity()
   p1(1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.25_R8P, [0.13_R8P, 0.07_R8P, -0.05_R8P])
   p0(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0.1_R8P, [3.0_R8P, 0.02_R8P, 0.01_R8P])
   p1(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], -0.15_R8P, [2.9_R8P, 0.05_R8P, 0.04_R8P])
   call swept_cells_3d_bodies(two, p0, p1, x, e, e, sw, nq=2_I4P, exact=.false., status=st)
   call swept_cells_3d(cube, p0(1), p1(1), x, e, e, s1, nq=2_I4P, exact=.false.)
   call swept_cells_3d(cube, p0(2), p1(2), x, e, e, s2, nq=2_I4P, exact=.false.)
   m = abs(st)
   do k = 1, size(e) - 1
      do j = 1, size(e) - 1
         do i = 1, size(x) - 1
            if (i <= nh) then
               m = m + cell_diff(sw, s1, i, j, k)
            else
               m = m + cell_diff(sw, s2, i, j, k)
            endif
         enddo
      enddo
   enddo
   nv = count(sw%swept /= 0._R8P)
   print '(A,I0,A,I0,A)', 'two far-apart bodies vs each body alone: mismatches ', m, ' (', nv, ' cells with a swept volume)'
   ok = ok .and. m == 0_I4P .and. nv > 0_I4P
   ! the swept volume is V1 - V0 of the union's cuts, bit for bit; it sums to zero (the solids keep their volume)
   call cut_cells_3d_bodies(two, p0, x, e, e, k0, status=st) ; m = abs(st)
   call cut_cells_3d_bodies(two, p1, x, e, e, k1, status=st) ; m = m + abs(st)
   do k = 1, size(e) - 1
      do j = 1, size(e) - 1
         do i = 1, size(x) - 1
            if (bits(sw%swept(i, j, k)) /= bits((k1%alpha(i, j, k) - k0%alpha(i, j, k)) * (x(i + 1) - x(i)) * &
                                                (e(j + 1) - e(j)) * (e(k + 1) - e(k)))) m = m + 1_I4P
         enddo
      enddo
   enddo
   total = sum(sw%swept)
   print '(A,I0,A,ES9.2)', 'swept = V1 - V0 of the union: mismatches ', m, ', |sum of swept| ', abs(total)
   ok = ok .and. m == 0_I4P .and. abs(total) < 1.e-14_R8P
   endsubroutine far_apart_case

   subroutine pair_case
   !< Case 4: two cubes 0.2 apart on one rigid motion vs their merged surface. The second cube's body frame is the first's
   !< (its vertices shifted by 1.2 along x), so both take the same poses and the screw motion is one rigid motion.
   type(surface_stl_object)  :: merged, pair(2)
   type(swept_cells_3d_t)    :: sw, sm
   type(cut_cells_3d_t)      :: kc
   type(piece_bodies_t)      :: pb
   type(rigid_pose_t)        :: q(2)
   integer(I4P), allocatable :: fb(:), fs(:)
   real(R8P),    allocatable :: x(:), e(:)
   real(R8P)                 :: dres
   integer(I4P)              :: i, st, m, nmix, p

   x = [(-0.45_R8P + 0.45_R8P * i, i = 0, 7)] ; e = [(-0.32_R8P + 0.45_R8P * i, i = 0, 4)]
   pair(1) = cube ; pair(2) = cube
   call pair(2)%move(pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [1.2_R8P, 0._R8P, 0._R8P]))
   q = pose_identity()
   call merge_posed_bodies(pair, q, merged, fb, fs, status=st) ! the pair as one body
   p0 = pose_identity()
   p1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.08_R8P, [0.06_R8P, 0.03_R8P, 0.02_R8P])
   call swept_cells_3d(merged, p0(1), p1(1), x, e, e, sm, nq=2_I4P, exact=.false.)
   call swept_cells_3d_bodies(pair, p0, p1, x, e, e, sw, nq=2_I4P, exact=.false., status=st)
   m = abs(st)
   ! pieces bounded by both bodies (their cells take the wall form)
   call cut_cells_3d_bodies(pair, p0, x, e, e, kc, status=st, piece_bodies=pb)
   nmix = 0_I4P
   do p = 1, pb%np
      if (pb%start(p + 1) - pb%start(p) >= 2_I4P) nmix = nmix + 1_I4P
   enddo
   m = m + nbits(sw%alpha0, sm%alpha0) + nbits(sw%alpha1, sm%alpha1) + nbits(sw%swept, sm%swept) + &
       nbits(sw%face_x, sm%face_x) + nbits(sw%face_y, sm%face_y) + nbits(sw%face_z, sm%face_z)
   dres = maxval(abs(sw%residual - sm%residual))
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2)', 'two cubes 0.2 apart on one motion vs the merged surface: mismatches ', m, &
         ', pieces bounded by both ', nmix, ', max |residual difference| ', dres, ', max |residual| ', maxval(abs(sm%residual))
   ok = ok .and. m == 0_I4P .and. nmix > 0_I4P .and. dres < 1.e-15_R8P
   endsubroutine pair_case

   subroutine refusal_case
   !< Case 5: refusals.
   type(swept_cells_3d_t) :: sw
   real(R8P)              :: e(5)
   integer(I4P)           :: i, st, m

   e = [(-0.23_R8P + 0.45_R8P * i, i = 0, 4)]
   p0(1) = pose_identity() ; p1(1) = p0(1)
   p0(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [1.5_R8P, 0._R8P, 0._R8P])
   p1(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0.5_R8P, 0._R8P, 0._R8P])
   m = 0_I4P
   call swept_cells_3d_bodies(two, p0, p0, e, e, e, sw, nq=2_I4P, status=st) ! exact absent
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_3d_bodies(two, p0, p0, e, e, e, sw, nq=2_I4P, exact=.true., status=st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_3d_bodies(two, p0(1:1), p0, e, e, e, sw, nq=2_I4P, exact=.false., status=st)
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   call swept_cells_3d_bodies(two, p0, p1, e, e, e, sw, nq=2_I4P, exact=.false., status=st) ! overlap at the end
   if (st /= STATUS_INVALID_INPUT) m = m + 1_I4P
   print '(A,I0)', 'refusals (exact absent, exact, poses, overlap): mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine refusal_case

   function cell_diff(a, b, i, j, k) result(n)
   !< Differing values of cell (i, j, k) of a and b: the cell arrays, its low faces (and its high ones on the last row).
   type(swept_cells_3d_t), intent(in) :: a, b    !< Swept geometries.
   integer(I4P),           intent(in) :: i, j, k !< Cell.
   integer(I4P)                       :: n       !< Mismatches.

   n = 0_I4P
   if (bits(a%alpha0(i, j, k)) /= bits(b%alpha0(i, j, k))) n = n + 1_I4P
   if (bits(a%alpha1(i, j, k)) /= bits(b%alpha1(i, j, k))) n = n + 1_I4P
   if (bits(a%swept(i, j, k)) /= bits(b%swept(i, j, k))) n = n + 1_I4P
   if (abs(a%residual(i, j, k) - b%residual(i, j, k)) > 1.e-16_R8P) n = n + 1_I4P ! round-off: see the case notes
   if (any(bits3(a%wall(:, i, j, k)) /= bits3(b%wall(:, i, j, k)))) n = n + 1_I4P
   if (bits(a%face_x(i, j, k)) /= bits(b%face_x(i, j, k))) n = n + 1_I4P
   if (bits(a%face_y(i, j, k)) /= bits(b%face_y(i, j, k))) n = n + 1_I4P
   if (bits(a%face_z(i, j, k)) /= bits(b%face_z(i, j, k))) n = n + 1_I4P
   if (bits(a%face_x(i + 1, j, k)) /= bits(b%face_x(i + 1, j, k))) n = n + 1_I4P
   if (bits(a%face_y(i, j + 1, k)) /= bits(b%face_y(i, j + 1, k))) n = n + 1_I4P
   if (bits(a%face_z(i, j, k + 1)) /= bits(b%face_z(i, j, k + 1))) n = n + 1_I4P
   if (a%cleared(i, j, k) .neqv. b%cleared(i, j, k)) n = n + 1_I4P
   if (a%covered(i, j, k) .neqv. b%covered(i, j, k)) n = n + 1_I4P
   endfunction cell_diff

   function nbits(a, b) result(n)
   !< Differing values of two real arrays of one shape (bit patterns).
   real(R8P), intent(in) :: a(..), b(..) !< Arrays.
   integer(I4P)          :: n            !< Mismatches.

   n = 0_I4P
   select rank (a)
   rank (3)
      select rank (b)
      rank (3)
         n = count(transfer(a, 1_I8P, size(a)) /= transfer(b, 1_I8P, size(b)))
      endselect
   endselect
   endfunction nbits

   elemental function bits(v) result(b)
   !< Bit pattern of a real.
   real(R8P), intent(in) :: v !< Value.
   integer(I8P)          :: b !< Bits.

   b = transfer(v, 1_I8P)
   endfunction bits

   function bits3(v) result(b)
   !< Bit patterns of a 3-vector.
   real(R8P), intent(in) :: v(3) !< Values.
   integer(I8P)          :: b(3) !< Bits.

   b = transfer(v, 1_I8P, 3)
   endfunction bits3
endprogram fossil_test_swept_bodies
