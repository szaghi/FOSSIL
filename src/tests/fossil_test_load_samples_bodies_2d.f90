!< FOSSIL test: the load samples of several curves (#109 WP3 step 02).
program fossil_test_load_samples_bodies_2d
!< FOSSIL test: the load samples of several curves (#109 WP3 step 02).
!<
!< 1. One curve: `load_samples_2d_bodies` equals `load_samples_2d` bit for bit in every single-body field (probes, lines,
!<    order 3; posed; the device path), body = 1, sub_body = the loop body, no second body; an empty second curve is
!<    ignored.
!< 2. Two circles far apart (posed): the samples of curve b are the slots of its single-body set, bit for bit, in curve
!<    order.
!< 3. Clearance: |x - c_other| - R to the chord sag (two of them: the point and the other polygon), the second body the
!<    other circle; the least clearance is the gap to the sag.
!< 4. Visibility: two squares 1/8 apart: on the facing sides a probe shorter than the gap is visible, a longer one is not
!<    (the single-body set sees both), clearance 1/8, second foot across, to round-off.
!< 5. Buoyancy: a circle and a square, posed: per curve F_b = -g A_b and M_b = A_b (x_c,b - c) x (-g) per unit depth to
!<    round-off (`integrate_loads` per_body = the curves).
!< 6. Twin: posed circles 0.06 apart with probes and lines: device = host in every field, bit for bit; resident.
!< 7. Refusals: size(poses), no curve, unmapped curves on the device path.
use fossil,          only : curve_object, load_samples_t, load_samples_2d, load_samples_2d_bodies, load_samples_reset, &
                            load_samples_mismatches, integrate_loads, rigid_pose_2d_t, pose_2d, to_world, STATUS_OK, &
                            STATUS_INVALID_INPUT
use iso_fortran_env, only : int64
use penf,            only : I4P, R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P), R = 0.5_R8P, G8 = 0.125_R8P
integer(I4P), parameter :: NV = 128_I4P
type(curve_object)      :: circle, empty, sqa, sqb
type(curve_object)      :: one(1), two(2), squares(2), mix(2)
type(load_samples_t)    :: sh, sb, sd, sa, sc, sr
type(rigid_pose_2d_t)   :: pa, pb
real(R8P)               :: u(NV), v(NV), e, emin, sag, c(2), gv(2), ar(2), xc(2, 2), fv(2), mv, ftot(2), mtot, sc1
real(R8P), allocatable  :: x(:), y(:), p(:), pf(:,:), pm(:)
integer(I4P)            :: s, st, nbad, nf, b
logical                 :: ok

ok = .true.
u = [(R * cos(2 * PI * s / NV), s = 0, NV - 1)] ; v = [(R * sin(2 * PI * s / NV), s = 0, NV - 1)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[NV])
sag = R * (1._R8P - cos(PI / NV))
x = [(-1.23_R8P + 0.05_R8P * s, s = 0, 49)] ; y = [(-0.71_R8P + 0.05_R8P * s, s = 0, 29)]
! 1. one curve
one(1) = circle
call load_samples_2d(circle, x, y, sh, status=st, order=3_I4P, lengths=[0.05_R8P, 0.2_R8P])
call load_samples_2d_bodies(one, x, y, sb, status=st, order=3_I4P, lengths=[0.05_R8P, 0.2_R8P])
call check('one curve: status, flags', st == STATUS_OK .and. sb%n == sh%n .and. sh%n == 3_I4P * NV .and. &
           sb%with_bodies .and. sb%bodies == 1_I4P)
call check('one curve: = load_samples_2d bit for bit in every single-body field', load_samples_mismatches(sh, strip(sb)) == 0)
call check('one curve: body 1, sub_body the loop body, no second body', all(sb%body == 1_I4P) .and. &
           all(sb%sub_body == sh%body) .and. all(sb%second_body == 0_I4P) .and. all(sb%clearance == huge(1._R8P)))
pa = pose_2d(0.3_R8P, [0.04_R8P, -0.03_R8P])
call load_samples_2d(circle, x, y, sh, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pa)
call load_samples_2d_bodies(one, x, y, sb, status=st, poses=[pa], lengths=[-0.02_R8P, 0.07_R8P])
call check('one curve posed: = load_samples_2d bit for bit', st == STATUS_OK .and. sb%n > 0_I4P .and. &
           load_samples_mismatches(sh, strip(sb)) == 0_I4P)
call circle%enter_device ; call one(1)%enter_device ! map the array elements themselves (a copy is not mapped)
call load_samples_2d(circle, x, y, sd, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pa, device=.true.)
call load_samples_2d_bodies(one, x, y, sc, status=st, poses=[pa], lengths=[-0.02_R8P, 0.07_R8P], device=.true.)
call check('one curve posed, device path: = load_samples_2d(device) and the host bodies set, bit for bit', &
           st == STATUS_OK .and. load_samples_mismatches(sd, strip(sc)) == 0_I4P .and. &
           load_samples_mismatches(sb, sc) == 0_I4P)
call circle%exit_device ; call one(1)%exit_device
two(1) = circle ; two(2) = empty
call load_samples_2d_bodies(two, x, y, sc, status=st, poses=[pa, pa], lengths=[-0.02_R8P, 0.07_R8P])
call check('an empty second curve is ignored', st == STATUS_OK .and. sc%bodies == 2_I4P .and. &
           load_samples_mismatches(sh, strip(sc)) == 0_I4P .and. all(sc%second_body == 0_I4P))
! 2. two circles far apart (gap 0.3), probes shorter than the gap
pa = pose_2d(0.2_R8P, [-0.65_R8P, 0.01_R8P]) ; pb = pose_2d(-0.4_R8P, [0.65_R8P, -0.02_R8P])
two(1) = circle ; two(2) = circle
call load_samples_2d_bodies(two, x, y, sb, status=st, poses=[pa, pb], lengths=[0.05_R8P, 0.12_R8P])
call load_samples_2d(circle, x, y, sa, lengths=[0.05_R8P, 0.12_R8P], pose=pa)
call load_samples_2d(circle, x, y, sc, lengths=[0.05_R8P, 0.12_R8P], pose=pb)
call check('far apart: the two sets in curve order', st == STATUS_OK .and. sb%n == sa%n + sc%n .and. &
           all(sb%body(1:sa%n) == 1_I4P) .and. all(sb%body(sa%n + 1:) == 2_I4P))
call check('far apart: every curve''s slots = its single-body set, bit for bit', same_slice(sb, 0_I4P, sa) .and. &
           same_slice(sb, sa%n, sc))
! 3. clearance
e = 0._R8P ; nbad = 0_I4P ; emin = huge(1._R8P)
do s = 1_I4P, sb%n
   c = pb%t ; if (sb%body(s) == 2_I4P) c = pa%t
   e = max(e, abs(sb%clearance(s) - (norm2(sb%points(:, s) - c) - R)))
   emin = min(emin, sb%clearance(s))
   if (sb%second_body(s) /= 3_I4P - sb%body(s)) nbad = nbad + 1_I4P
   if (abs(norm2(sb%second_foot(:, s) - sb%points(:, s)) - sb%clearance(s)) > 1.e-13_R8P) nbad = nbad + 1_I4P
enddo
call check('clearance: |x - c_other| - R (chord sag), second body and foot', e <= 2.0001_R8P * sag .and. nbad == 0_I4P)
call check('clearance: the least one is the gap 0.3 (chord sag)', abs(emin - (norm2(pb%t - pa%t) - 2 * R)) <= 2.0001 * sag)
! 4. visibility across a narrow gap: two squares 1/8 apart
call sqa%adopt_loops(u=[-0.6_R8P, -0.5_R8P * G8, -0.5_R8P * G8, -0.6_R8P], v=[-0.4_R8P, -0.4_R8P, 0.4_R8P, 0.4_R8P], &
                     loop_sizes=[4_I4P])
call sqb%adopt_loops(u=[0.5_R8P * G8, 0.6_R8P, 0.6_R8P, 0.5_R8P * G8], v=[-0.4_R8P, -0.4_R8P, 0.4_R8P, 0.4_R8P], &
                     loop_sizes=[4_I4P])
squares(1) = sqa ; squares(2) = sqb
call load_samples_2d_bodies(squares, samples=sb, status=st, lengths=[0.1_R8P, 0.2_R8P])
call load_samples_2d(sqa, samples=sa, lengths=[0.1_R8P, 0.2_R8P])
nf = 0_I4P ; nbad = 0_I4P
do s = 1_I4P, sb%n
   if (abs(sb%points(1, s) + sb%normal(1, s) * 0.5_R8P * G8) > 1.e-14_R8P .or. abs(sb%normal(1, s)) < 0.5_R8P) then
      if (.not. all(sb%visible(:, s))) nbad = nbad + 1_I4P ! no other side is shadowed
      cycle
   endif
   nf = nf + 1_I4P ! a sample on a facing side
   if (.not. sb%visible(1, s) .or. sb%visible(2, s)) nbad = nbad + 1_I4P
   if (abs(sb%clearance(s) - G8) > 1.e-14_R8P .or. sb%second_body(s) /= 3_I4P - sb%body(s)) nbad = nbad + 1_I4P
   if (abs(sb%second_foot(1, s) + sb%points(1, s)) > 1.e-14_R8P .or. abs(sb%second_foot(2, s) - sb%points(2, s)) > &
       1.e-14_R8P) nbad = nbad + 1_I4P
   if (sb%body(s) == 1_I4P) then
      if (.not. all(sa%visible(:, s))) nbad = nbad + 1_I4P ! alone, square a does not block it
   endif
enddo
call check('narrow gap: facing sides, the long probe enters the other square, clearance 1/8, second foot across', &
           st == STATUS_OK .and. nf == 4_I4P .and. nbad == 0_I4P)
! 5. buoyancy per curve: a circle and a square, posed
pa = pose_2d(0.4_R8P, [-0.8_R8P, 0.1_R8P]) ; pb = pose_2d(-0.2_R8P, [1.1_R8P, -0.3_R8P])
mix(1) = circle ; mix(2) = sqa
gv = [0.3_R8P, -9.81_R8P] ; c = [0.1_R8P, 0.2_R8P]
ar(1) = 0.5_R8P * NV * R * R * sin(2 * PI / NV) ; ar(2) = (0.6_R8P - 0.5_R8P * G8) * 0.8_R8P
xc(:, 1) = to_world(pa, [0._R8P, 0._R8P]) ; xc(:, 2) = to_world(pb, [-0.5_R8P * (0.6_R8P + 0.5_R8P * G8), 0._R8P])
call load_samples_2d_bodies(mix, samples=sb, status=st, poses=[pa, pb])
allocate(p(sb%n))
do s = 1_I4P, sb%n
   p(s) = 2._R8P + dot_product(gv, sb%points(:, s))
enddo
call integrate_loads(sb, ftot, mtot, status=st, pressure=p, about=c, per_body_force=pf, per_body_moment=pm)
e = 0._R8P
do b = 1_I4P, 2_I4P
   fv = -gv * ar(b)
   mv = ar(b) * ((xc(1, b) - c(1)) * (-gv(2)) - (xc(2, b) - c(2)) * (-gv(1)))
   sc1 = norm2(gv) * ar(b) * 3._R8P
   e = max(e, maxval(abs(pf(:, b) - fv)) / sc1, abs(pm(b) - mv) / sc1)
enddo
call check('buoyancy: per curve F_b = -g A_b, M_b = A_b (x_c,b - c) x (-g) to round-off', st == STATUS_OK .and. &
           size(pf, 2) == 2 .and. e < 1.e-13_R8P)
call check('buoyancy: the totals are the sums of the curves', &
           maxval(abs(ftot - (pf(:, 1) + pf(:, 2)))) < 1.e-13_R8P * norm2(ftot) .and. &
           abs(mtot - (pm(1) + pm(2))) < 1.e-12_R8P * (abs(mtot) + norm2(ftot)))
! 6. twin: posed circles 0.06 apart with probes and lines
pa = pose_2d(0.4_R8P, [-0.53_R8P, 0.013_R8P]) ; pb = pose_2d(-0.2_R8P, [0.53_R8P, -0.011_R8P])
two(1) = circle ; two(2) = circle ; call two(1)%enter_device ; call two(2)%enter_device
call load_samples_2d_bodies(two, x, y, sh, status=st, poses=[pa, pb], lengths=[0.05_R8P, 0.1_R8P])
call load_samples_2d_bodies(two, x, y, sd, status=st, poses=[pa, pb], lengths=[0.05_R8P, 0.1_R8P], device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. sh%n > 0_I4P .and. .not. all(sh%visible) .and. &
           load_samples_mismatches(sh, sd) == 0_I4P)
call load_samples_2d_bodies(two, x, y, sr, status=st, poses=[pa, pb], lengths=[0.05_R8P, 0.1_R8P], resident=.true.)
call sr%update_host
call check('twin: resident, update_host = host', st == STATUS_OK .and. sr%on_device .and. &
           load_samples_mismatches(sh, sr) == 0_I4P)
call load_samples_reset(sr)
call two(1)%exit_device ; call two(2)%exit_device
! 7. refusals
call load_samples_2d_bodies(two, samples=sb, status=st, poses=[pa])
call check('refused: size(poses) /= size(curves)', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)
call load_samples_2d_bodies(two(1:0), samples=sb, status=st)
call check('refused: no curve', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)
call load_samples_2d_bodies(two, samples=sb, status=st, device=.true.)
call check('refused: device path of unmapped curves (stderr message)', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)

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

   function strip(c) result(s)
   !< A bodies set without its several-bodies fields (a single-body set's shapes and flags).
   type(load_samples_t), intent(in) :: c !< Bodies set.
   type(load_samples_t)             :: s !< Stripped copy.

   s = c ; s%on_device = .false. ; s%device_key = 0
   s%with_bodies = .false. ; s%bodies = 0_I4P
   deallocate(s%sub_body, s%clearance, s%second_body, s%second_foot)
   allocate(s%sub_body(0), s%clearance(0), s%second_body(0), s%second_foot(c%dims, 0))
   endfunction strip

   logical function same_slice(c, off, a) result(same)
   !< The slots off+1 ... off+a%n of the bodies set c equal the single-body set a, every field by bits.
   type(load_samples_t), intent(in) :: c, a !< Bodies set, single-body set.
   integer(I4P),         intent(in) :: off  !< Offset.
   integer(I4P) :: n                        !< Samples.

   n = a%n
   same = eqr(reshape(c%points(:, off + 1:off + n), [2 * n]), reshape(a%points, [2 * n])) .and. &
          eqr(c%weight(off + 1:off + n), a%weight) .and. &
          eqr(reshape(c%normal(:, off + 1:off + n), [2 * n]), reshape(a%normal, [2 * n])) .and. &
          all(c%element(off + 1:off + n) == a%element) .and. all(c%patch(off + 1:off + n) == a%patch) .and. &
          all(c%sub_body(off + 1:off + n) == a%body) .and. &
          eqr(reshape(c%probe(:, :, off + 1:off + n), [size(a%probe)]), reshape(a%probe, [size(a%probe)])) .and. &
          all(c%visible(:, off + 1:off + n) .eqv. a%visible) .and. all(c%probe_cell(:, :, off + 1:off + n) == a%probe_cell)
   endfunction same_slice

   logical function eqr(a, b)
   !< Equal bits.
   real(R8P), intent(in) :: a(:), b(:) !< Reals.

   eqr = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
   endfunction eqr
endprogram fossil_test_load_samples_bodies_2d
