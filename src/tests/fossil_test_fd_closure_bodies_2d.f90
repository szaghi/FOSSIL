!< FOSSIL test: the closure table of several bodies on a 2-D block (#109 WP3 step 01).
program fossil_test_fd_closure_bodies_2d
!< FOSSIL test: the closure table of several bodies on a 2-D block (#109 WP3 step 01).
!<
!< 1. One curve: `fd_closure_2d_bodies` equals `fd_closure_2d` bit for bit in every single-body field (curvature, probes,
!<    velocities, ILW; posed; the device path), body = 1, sub_body = the facet_body, no second body; an empty second curve
!<    is ignored.
!< 2. Two circles far apart (posed, own velocities): the ghosts are the union of each body's, every ghost equals the
!<    single-body table of its body bit for bit; clearance |foot - c_other| - R to the chord sag.
!< 3. Narrow gap: two squares 1/8 apart: the facing feet have clearance 1/8 to round-off, the second foot across, shared;
!<    the probe longer than the gap is not visible, the shorter one is; shared_distance below the gap shares none.
!< 4. Overlap: two circles 0.4 apart (R = 1/2): overlapping ghosts inside both, |d| the smaller of the two.
!< 5. Twin: posed circles 0.06 apart with velocities, probes and ILW: device = host in every field, bit for bit;
!<    resident; the velocity V_b + omega_b (-(p_v - c_v), p_u - c_u) of the owner, bit for bit.
!< 6. Refusals: size(poses), size(angular_velocities), no curve, unmapped curves on the device path.
use fossil, only : curve_object, ib_closure_t, fd_closure_2d, fd_closure_2d_bodies, ib_closure_reset, ib_closure_mismatches, &
                   rigid_pose_2d_t, pose_2d, STATUS_OK, STATUS_INVALID_INPUT
use iso_fortran_env, only : int64
use penf,   only : I4P, R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P), R = 0.5_R8P, G8 = 0.125_R8P
integer(I4P), parameter :: NV = 128_I4P
type(curve_object)      :: circle, empty, sqa, sqb
type(curve_object)      :: one(1), two(2), squares(2)
type(ib_closure_t)      :: ch, cb, cd, ca, cc, cr
type(rigid_pose_2d_t)   :: pa, pb, pp(2)
real(R8P)               :: u(NV), v(NV), e, sag, c(2), w2(2), vu, vv, vl(2, 2), cl(2, 2)
real(R8P), allocatable  :: x(:), y(:), xb(:), qd(:)
integer(I4P)            :: g, k, st, nbad, nf, ns, b
logical                 :: ok

ok = .true.
u = [(R * cos(2 * PI * g / NV), g = 0, NV - 1)] ; v = [(R * sin(2 * PI * g / NV), g = 0, NV - 1)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[NV]) ; call circle%build_curvature(status=st)
call check('circle: curvature built', st == STATUS_OK)
sag = R * (1._R8P - cos(PI / NV))
x = [(-1.23_R8P + 0.05_R8P * g, g = 0, 49)] ; y = [(-0.71_R8P + 0.05_R8P * g, g = 0, 29)]
vl = reshape([0.1_R8P, -0.2_R8P, -0.3_R8P, 0.05_R8P], [2, 2]) ; w2 = [2._R8P, -1.5_R8P]
cl = reshape([0.01_R8P, -0.02_R8P, 0.2_R8P, 0.1_R8P], [2, 2])
! 1. one curve
one(1) = circle
call fd_closure_2d(circle, x, y, [2_I4P, 3_I4P], ch, status=st, lengths=[0.05_R8P, 0.1_R8P], linear_velocity=vl(:, 1), &
                   angular_velocity=w2(1), centre=cl(:, 1), ilw=.true.)
call fd_closure_2d_bodies(one, x, y, [2_I4P, 3_I4P], cb, status=st, lengths=[0.05_R8P, 0.1_R8P], &
                          linear_velocities=vl(:, 1:1), angular_velocities=w2(1:1), centres=cl(:, 1:1), ilw=.true.)
call check('one curve: status, flags', st == STATUS_OK .and. cb%n == ch%n .and. ch%n > 50_I4P .and. cb%with_bodies .and. &
           cb%bodies == 1_I4P .and. cb%with_curvature .and. cb%with_velocity .and. cb%with_ilw)
call check('one curve: = fd_closure_2d bit for bit in every single-body field', ib_closure_mismatches(ch, strip(cb)) == 0_I4P)
call check('one curve: body 1, sub_body the facet_body, no second body', all(cb%body == 1_I4P) .and. &
           all(cb%sub_body == ch%body) .and. all(cb%second_body == 0_I4P) .and. .not. any(cb%shared) .and. &
           .not. any(cb%overlap) .and. all((cb%ilw_body == 1_I4P) .eqv. (cb%ilw_element > 0_I4P)))
pa = pose_2d(0.3_R8P, [0.04_R8P, -0.03_R8P])
call fd_closure_2d(circle, x, y, [3_I4P, 1_I4P], ch, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pa, &
                   linear_velocity=vl(:, 1), angular_velocity=w2(1), ilw=.true.)
call fd_closure_2d_bodies(one, x, y, [3_I4P, 1_I4P], cb, status=st, lengths=[-0.02_R8P, 0.07_R8P], poses=[pa], &
                          linear_velocities=vl(:, 1:1), angular_velocities=w2(1:1), ilw=.true.)
call check('one curve posed: = fd_closure_2d bit for bit', st == STATUS_OK .and. cb%n > 0_I4P .and. &
           ib_closure_mismatches(ch, strip(cb)) == 0_I4P)
call circle%enter_device ; call one(1)%enter_device ! map the array elements themselves (a copy is not mapped)
call fd_closure_2d(circle, x, y, [3_I4P, 1_I4P], cd, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pa, &
                   linear_velocity=vl(:, 1), angular_velocity=w2(1), ilw=.true., device=.true.)
call fd_closure_2d_bodies(one, x, y, [3_I4P, 1_I4P], cc, status=st, lengths=[-0.02_R8P, 0.07_R8P], poses=[pa], &
                          linear_velocities=vl(:, 1:1), angular_velocities=w2(1:1), ilw=.true., device=.true.)
call check('one curve posed, device path: = fd_closure_2d(device) and the host bodies table, bit for bit', &
           st == STATUS_OK .and. ib_closure_mismatches(cd, strip(cc)) == 0_I4P .and. ib_closure_mismatches(cb, cc) == 0_I4P)
call circle%exit_device ; call one(1)%exit_device
two(1) = circle ; two(2) = empty
call fd_closure_2d_bodies(two, x, y, [3_I4P, 1_I4P], cc, status=st, lengths=[-0.02_R8P, 0.07_R8P], poses=[pa, pa], &
                          linear_velocities=vl, angular_velocities=w2, ilw=.true.)
call check('an empty second curve is ignored', st == STATUS_OK .and. cc%bodies == 2_I4P .and. &
           ib_closure_mismatches(ch, strip(cc)) == 0_I4P .and. all(cc%second_body == 0_I4P))
! 2. two circles far apart (gap 0.3 > reach h)
pa = pose_2d(0.2_R8P, [-0.65_R8P, 0.01_R8P]) ; pb = pose_2d(-0.4_R8P, [0.65_R8P, -0.02_R8P])
two(1) = circle ; two(2) = circle
call fd_closure_2d_bodies(two, x, y, [2_I4P, 2_I4P], cb, status=st, lengths=[0.05_R8P, 0.12_R8P], poses=[pa, pb], &
                          linear_velocities=vl, angular_velocities=w2, centres=cl, ilw=.true.)
call fd_closure_2d(circle, x, y, [2_I4P, 2_I4P], ca, lengths=[0.05_R8P, 0.12_R8P], pose=pa, linear_velocity=vl(:, 1), &
                   angular_velocity=w2(1), centre=cl(:, 1), ilw=.true.)
call fd_closure_2d(circle, x, y, [2_I4P, 2_I4P], cc, lengths=[0.05_R8P, 0.12_R8P], pose=pb, linear_velocity=vl(:, 2), &
                   angular_velocity=w2(2), centre=cl(:, 2), ilw=.true.)
call check('far apart: the ghosts are the union of each body''s', st == STATUS_OK .and. cb%n == ca%n + cc%n .and. &
           count(cb%body == 1_I4P) == ca%n .and. count(cb%body == 2_I4P) == cc%n .and. cb%overlap_points == 0_I4P)
nbad = 0_I4P
do g = 1_I4P, cb%n
   if (cb%body(g) == 1_I4P) then
      if (.not. same_ghost(cb, g, ca, 1_I4P)) nbad = nbad + 1_I4P
   else
      if (.not. same_ghost(cb, g, cc, 2_I4P)) nbad = nbad + 1_I4P
   endif
enddo
call check('far apart: every ghost = the single-body table of its body, bit for bit', nbad == 0_I4P)
e = 0._R8P ; nbad = 0_I4P
do g = 1_I4P, cb%n
   c = pb%t ; if (cb%body(g) == 2_I4P) c = pa%t
   e = max(e, abs(cb%clearance(g) - (norm2(cb%foot(:, g) - c) - R)))
   if (cb%second_body(g) /= 3_I4P - cb%body(g)) nbad = nbad + 1_I4P
enddo
call check('far apart: clearance |foot - c_other| - R (chord sag), second body', e <= 2.0001_R8P * sag .and. &
           nbad == 0_I4P .and. .not. any(cb%shared) .and. abs(cb%shared_distance - 0.1_R8P) < 1.e-12_R8P)
! 3. narrow gap: two squares 1/8 apart
call sqa%adopt_loops(u=[-0.6_R8P, -0.5_R8P * G8, -0.5_R8P * G8, -0.6_R8P], v=[-0.4_R8P, -0.4_R8P, 0.4_R8P, 0.4_R8P], &
                     loop_sizes=[4_I4P])
call sqb%adopt_loops(u=[0.5_R8P * G8, 0.6_R8P, 0.6_R8P, 0.5_R8P * G8], v=[-0.4_R8P, -0.4_R8P, 0.4_R8P, 0.4_R8P], &
                     loop_sizes=[4_I4P])
squares(1) = sqa ; squares(2) = sqb
xb = [(-0.85_R8P + 0.1_R8P * g, g = 0, 17)]
call fd_closure_2d_bodies(squares, xb, y, [2_I4P, 2_I4P], cb, status=st, lengths=[0.1_R8P, 0.2_R8P])
call fd_closure_2d(sqa, xb, y, [2_I4P, 2_I4P], ca, lengths=[0.1_R8P, 0.2_R8P])
nf = 0_I4P ; nbad = 0_I4P ; ns = 0_I4P
do g = 1_I4P, cb%n
   if (cb%shared(g) .neqv. (cb%clearance(g) < cb%shared_distance)) nbad = nbad + 1_I4P
   if (cb%shared(g)) ns = ns + 1_I4P
   if (cb%body(g) /= 1_I4P .or. abs(cb%foot(1, g) + 0.5_R8P * G8) > 1.e-14_R8P .or. abs(cb%foot(2, g)) > 0.35_R8P) cycle
   nf = nf + 1_I4P
   if (abs(cb%clearance(g) - G8) > 1.e-14_R8P .or. .not. cb%shared(g) .or. cb%second_body(g) /= 2_I4P) nbad = nbad + 1_I4P
   if (abs(cb%second_foot(1, g) - 0.5_R8P * G8) > 1.e-14_R8P .or. abs(cb%second_foot(2, g) - cb%foot(2, g)) > 1.e-14_R8P) &
      nbad = nbad + 1_I4P
   if (.not. cb%visible(1, g) .or. cb%visible(2, g)) nbad = nbad + 1_I4P
   k = ghost_of(ca, cb%ijk(:, g))
   if (k < 1_I4P) then
      nbad = nbad + 1_I4P
   elseif (.not. all(ca%visible(:, k))) then
      nbad = nbad + 1_I4P
   endif
enddo
call check('narrow gap: facing feet, clearance 1/8, second foot, shared, the long probe enters the other square', &
           st == STATUS_OK .and. nf > 8_I4P .and. nbad == 0_I4P .and. ns > nf .and. ns < cb%n)
call fd_closure_2d_bodies(squares, xb, y, [2_I4P, 2_I4P], cb, status=st, shared_distance=0.1_R8P)
call check('narrow gap: shared_distance 0.1 < 1/8 shares none', st == STATUS_OK .and. .not. any(cb%shared))
! 4. overlap: two circles 0.4 apart
pp(1) = pose_2d(0._R8P, [-0.2_R8P, 0._R8P]) ; pp(2) = pose_2d(0._R8P, [0.2_R8P, 0._R8P])
call fd_closure_2d_bodies(two, x, y, [2_I4P, 2_I4P], cb, status=st, poses=pp)
nbad = 0_I4P
allocate(qd(cb%n))
do g = 1_I4P, cb%n
   if (cb%overlap(g)) then
      if (norm2(cb%points(:, g) - pp(1)%t) > R + sag .or. norm2(cb%points(:, g) - pp(2)%t) > R + sag) nbad = nbad + 1_I4P
   else
      if (norm2(cb%points(:, g) - pp(1)%t) < R - sag .and. norm2(cb%points(:, g) - pp(2)%t) < R - sag) nbad = nbad + 1_I4P
   endif
enddo
do b = 1_I4P, 2_I4P
   call circle%query_many(points=cb%points, distances=qd, pose=pp(b))
   do g = 1_I4P, cb%n
      if (abs(cb%distance(g)) > abs(qd(g))) nbad = nbad + 1_I4P
   enddo
enddo
call check('overlap: points inside both counted, overlapping ghosts inside both, |d| the smaller', st == STATUS_OK .and. &
           cb%overlap_points > count(cb%overlap) .and. count(cb%overlap) > 4 .and. nbad == 0_I4P)
! 5. twin: posed circles 0.06 apart with velocities, probes and ILW
pa = pose_2d(0.4_R8P, [-0.53_R8P, 0.013_R8P]) ; pb = pose_2d(-0.2_R8P, [0.53_R8P, -0.011_R8P])
two(1) = circle ; two(2) = circle ; call two(1)%enter_device ; call two(2)%enter_device
call fd_closure_2d_bodies(two, x, y, [2_I4P, 3_I4P], ch, status=st, lengths=[0.05_R8P, 0.1_R8P], poses=[pa, pb], &
                          linear_velocities=vl, angular_velocities=w2, ilw=.true.)
call fd_closure_2d_bodies(two, x, y, [2_I4P, 3_I4P], cd, status=st, lengths=[0.05_R8P, 0.1_R8P], poses=[pa, pb], &
                          linear_velocities=vl, angular_velocities=w2, ilw=.true., device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. ch%n > 0_I4P .and. any(ch%shared) .and. &
           .not. all(ch%visible) .and. any(ch%ilw_body == 2_I4P) .and. ib_closure_mismatches(ch, cd) == 0_I4P)
call fd_closure_2d_bodies(two, x, y, [2_I4P, 3_I4P], cr, status=st, lengths=[0.05_R8P, 0.1_R8P], poses=[pa, pb], &
                          linear_velocities=vl, angular_velocities=w2, ilw=.true., resident=.true.)
call cr%update_host
call check('twin: resident, update_host = host', st == STATUS_OK .and. cr%on_device .and. &
           ib_closure_mismatches(ch, cr) == 0_I4P)
call ib_closure_reset(cr)
nbad = 0_I4P
do g = 1_I4P, ch%n
   b = ch%body(g) ; c = pa%t ; if (b == 2_I4P) c = pb%t
   vu = vl(1, b) - w2(b) * (ch%foot(2, g) - c(2)) ; vv = vl(2, b) + w2(b) * (ch%foot(1, g) - c(1))
   if (vu /= ch%velocity(1, g) .or. vv /= ch%velocity(2, g)) nbad = nbad + 1_I4P
enddo
call check('twin: velocity V_b + omega_b (-(p_v - c_v), p_u - c_u) of the owner, bit for bit', nbad == 0_I4P)
call two(1)%exit_device ; call two(2)%exit_device
! 6. refusals
call fd_closure_2d_bodies(two, x, y, [2_I4P, 2_I4P], cb, status=st, poses=[pa])
call check('refused: size(poses) /= size(curves)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_2d_bodies(two, x, y, [2_I4P, 2_I4P], cb, status=st, angular_velocities=w2(1:1))
call check('refused: size(angular_velocities) /= size(curves)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_2d_bodies(two(1:0), x, y, [2_I4P, 2_I4P], cb, status=st)
call check('refused: no curve', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_2d_bodies(two, x, y, [2_I4P, 2_I4P], cb, status=st, device=.true.)
call check('refused: device path of unmapped curves (stderr message)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)

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
   !< A bodies table without its several-bodies fields (a single-body table's shapes and flags).
   type(ib_closure_t), intent(in) :: c !< Bodies table.
   type(ib_closure_t)             :: s !< Stripped copy.

   s = c ; s%on_device = .false. ; s%device_key = 0
   s%with_bodies = .false. ; s%bodies = 0_I4P ; s%overlap_points = 0_I4P ; s%shared_distance = 0._R8P
   deallocate(s%sub_body, s%clearance, s%second_body, s%second_foot, s%shared, s%overlap, s%ilw_body)
   allocate(s%sub_body(0), s%clearance(0), s%second_body(0), s%second_foot(c%dims, 0), s%shared(0), s%overlap(0), &
            s%ilw_body(2 * c%dims, 0))
   endfunction strip

   integer(I4P) function ghost_of(c, ijk) result(k)
   !< The ghost of table c at indexes ijk (0: none).
   type(ib_closure_t), intent(in) :: c      !< Table.
   integer(I4P),       intent(in) :: ijk(:) !< Indexes.
   integer(I4P) :: g                        !< Counter.

   k = 0_I4P
   do g = 1_I4P, c%n
      if (all(c%ijk(:, g) == ijk)) then
         k = g ; return
      endif
   enddo
   endfunction ghost_of

   logical function same_ghost(cb, g, cs, b) result(same)
   !< Ghost g of the bodies table cb equals the ghost at the same indexes of body b's single-body table cs, every field
   !< by bits (body = b, sub_body = cs's body; ilw_body = b where a crossing).
   type(ib_closure_t), intent(in) :: cb, cs !< Bodies and single-body tables.
   integer(I4P),       intent(in) :: g, b   !< Ghost, body.
   integer(I4P) :: s                        !< Single-body ghost.

   same = .false.
   s = ghost_of(cs, cb%ijk(:, g))
   if (s < 1_I4P) return
   same = eqr(cb%points(:, g), cs%points(:, s)) .and. cb%mask(g) == cs%mask(s) .and. all(cb%steps(:, g) == cs%steps(:, s))
   same = same .and. eqr([cb%distance(g)], [cs%distance(s)]) .and. eqr(cb%foot(:, g), cs%foot(:, s)) .and. &
          eqr(cb%normal(:, g), cs%normal(:, s)) .and. cb%element(g) == cs%element(s) .and. &
          cb%feature(g) == cs%feature(s) .and. cb%patch(g) == cs%patch(s) .and. cb%body(g) == b .and. &
          cb%sub_body(g) == cs%body(s)
   same = same .and. eqr(reshape(cb%probe(:, :, g), [size(cb%probe(:, :, g))]), &
                         reshape(cs%probe(:, :, s), [size(cs%probe(:, :, s))])) .and. &
          all(cb%visible(:, g) .eqv. cs%visible(:, s)) .and. all(cb%probe_cell(:, :, g) == cs%probe_cell(:, :, s))
   same = same .and. eqr(cb%kappa(:, g), cs%kappa(:, s)) .and. eqr(cb%direction(:, 1, g), cs%direction(:, 1, s)) .and. &
          eqr(cb%velocity(:, g), cs%velocity(:, s))
   same = same .and. eqr(cb%ilw_t(:, g), cs%ilw_t(:, s)) .and. all(cb%ilw_element(:, g) == cs%ilw_element(:, s)) .and. &
          eqr(reshape(cb%ilw_normal(:, :, g), [8]), reshape(cs%ilw_normal(:, :, s), [8])) .and. &
          eqr(cb%ilw_kappa(1, :, g), cs%ilw_kappa(1, :, s)) .and. &
          eqr(reshape(cb%ilw_direction(:, 1, :, g), [8]), reshape(cs%ilw_direction(:, 1, :, s), [8])) .and. &
          all((cb%ilw_body(:, g) == b) .eqv. (cs%ilw_element(:, s) > 0_I4P))
   endfunction same_ghost

   logical function eqr(a, b)
   !< Equal bits.
   real(R8P), intent(in) :: a(:), b(:) !< Reals.

   eqr = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
   endfunction eqr
endprogram fossil_test_fd_closure_bodies_2d
