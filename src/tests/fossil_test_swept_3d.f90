!< FOSSIL test: 3-D swept geometry between two poses (issue #41 step 5c).
program fossil_test_swept_3d
!< FOSSIL test: 3-D swept geometry between two poses (issue #41 step 5c).
!<
!< - The unit cube translating across a 6^3 block, exact path (default): the time-averaged face fluid
!<   areas equal an independent analytic integration (breakpoints from the cube's faces, three-point
!<   Gauss on every piece; the solid area of a face is piecewise quadratic with jumps) to round-off;
!<   the uniform-flow residual is round-off; the projected walls sum to 0 over the block.
!< - Exact path on bodies at a rotated pose, translating: the cube and the L-prism (non-convex)
!<   rotated about (1, 1, 1), the 32-face cone: uniform-flow residual at round-off.
!< - Negative control, `exact=.false.` (Gauss-Legendre, no breakpoints): the cube faces approach the
!<   analytic ones only as nq = 2, 4, 8, 16 grows, at about first order (jumps).
!< - The cube rotating by 0.3 rad about the (1, 1, 1) axis through its centre while translating:
!<   the residual falls with nq = 2, 4, 8.
!< - The poses at the Gauss nodes given explicitly (the screw poses) give the same face areas bit
!<   for bit; the 3-D `pose_interpolate` returns the end poses exactly and the mid pose at half the
!<   relative angle.
!< - Refusals: decreasing edges, nq = 0.

use fossil,                   only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, &
                                     swept_cells_3d_t, swept_cells_3d, gauss_legendre_01, STATUS_OK, STATUS_INVALID_INPUT
use fossil_rigid_pose,        only : pose_log_3d
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable       :: b(:)
type(surface_stl_object)        :: s
type(swept_cells_3d_t)          :: sw, sw2
type(surface_stl_object)        :: sb
real(R8P)                       :: e2(7), rex(3), err0, res0
type(rigid_pose_t)              :: p0, p1, pm, pe0, pe1
type(rigid_pose_t), allocatable :: pq(:)
real(R8P), allocatable          :: sq(:), wq(:)
real(R8P)                       :: e(7), err(4), res(4), d(3), ax(3), ang, ang_m, w(3), c(3)
integer(I4P)                    :: i, k, st, st2
logical                         :: ok(7)

call bench_bodies(b)
call to_surface(b(1)%soup, s) ! unit cube [0, 1]^3
e = [(-0.23_R8P + 0.3_R8P * i, i = 0, 6)]
d = [0.31_R8P, 0.17_R8P, 0.23_R8P]
p0 = pose_identity() ; p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, d)
call swept_cells_3d(s, p0, p1, e, e, e, sw, status=st)
err0 = face_error(sw)
res0 = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
w = [sum(sw%wall(1, :, :, :)), sum(sw%wall(2, :, :, :)), sum(sw%wall(3, :, :, :))]
ok(1) = st == STATUS_OK .and. sw%exact .and. err0 < 1.e-14_R8P .and. res0 < 1.e-14_R8P .and. maxval(abs(w)) < 1.e-13_R8P
print '(A,ES9.1,A,ES9.1,A,3ES9.1,A,I0,A,L1)', 'cube translating, exact: max face error vs analytic ', err0, &
      '  max residual ', res0, '  sum of walls ', w, '  evaluations/plane ', sw%evaluations, '  ok ', ok(1)
do k = 1, 4
   call swept_cells_3d(s, p0, p1, e, e, e, sw, nq=2_I4P**k, exact=.false., status=st)
   err(k) = face_error(sw)
   res(k) = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
enddo
ok(6) = st == STATUS_OK .and. .not. sw%exact .and. err(2) < err(1) .and. err(3) < err(2) .and. err(4) < err(3) .and. &
        res(4) < res(1) .and. err(4) > 1.e-4_R8P
print '(A,4ES9.1,A,4ES9.1,A,L1)', 'negative control (Gauss): max face error nq=2,4,8,16 ', err, '  max residual ', res, &
      '  ok ', ok(6)
! exact path on rotated / non-convex bodies, translating
e2 = [(-0.6_R8P + 0.45_R8P * i, i = 0, 6)]
c = [0.5_R8P, 0.5_R8P, 0.5_R8P]
p0 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P)
p0%t = c - matmul(p0%r, c) ; p1 = p0 ; p1%t = p0%t + d
call swept_cells_3d(s, p0, p1, e, e, e, sw, status=st)
rex(1) = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
ok(7) = st == STATUS_OK .and. sw%exact
call to_surface(b(8)%soup, sb) ! L-prism [0, 2]^2 x [0, 1]
c = [1._R8P, 1._R8P, 0.5_R8P]
p0 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P)
p0%t = c - matmul(p0%r, c) ; p1 = p0 ; p1%t = p0%t + 1.7_R8P * d
call swept_cells_3d(sb, p0, p1, e2, e2, e2, sw, status=st)
rex(2) = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
ok(7) = ok(7) .and. st == STATUS_OK .and. sw%exact
call to_surface(b(7)%soup, sb) ! cone, 32 faces
p0 = pose_identity() ; p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0.13_R8P, -0.29_R8P, 0.41_R8P])
call swept_cells_3d(sb, p0, p1, e2, e2, e2, sw, status=st)
rex(3) = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
ok(7) = ok(7) .and. st == STATUS_OK .and. sw%exact .and. maxval(rex) < 1.e-13_R8P
print '(A,3ES9.1,A,L1)', 'exact path, max residual: rotated cube, rotated L-prism, cone ', rex, '  ok ', ok(7)
p0 = pose_identity()
! rotation about (1, 1, 1) through the centre, plus the translation
c = [0.5_R8P, 0.5_R8P, 0.5_R8P]
p1 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P)
p1%t = c - matmul(p1%r, c) + d
do k = 1, 3
   call swept_cells_3d(s, p0, p1, e, e, e, sw, nq=2_I4P**k, status=st)
   res(k) = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
enddo
ok(2) = st == STATUS_OK .and. res(2) < res(1) .and. res(3) < res(2)
print '(A,3ES9.1,A,ES9.1,A,L1)', 'cube rotating and translating: max residual nq=2,4,8 ', res(1:3), '  estimate (nq=8) ', &
      sw%error_estimate, '  ok ', ok(2)
! explicit poses, interpolation
call gauss_legendre_01(3_I4P, sq, wq)
allocate(pq(3))
do i = 1, 3
   pq(i) = pose_interpolate(p0, p1, sq(i))
enddo
call swept_cells_3d(s, p0, p1, e, e, e, sw, nq=3_I4P, status=st)
call swept_cells_3d(s, p0, p1, e, e, e, sw2, nq=3_I4P, poses=pq, status=st2)
ok(3) = st == STATUS_OK .and. st2 == STATUS_OK .and. all(sw%face_x == sw2%face_x) .and. all(sw%face_y == sw2%face_y) .and. &
        all(sw%face_z == sw2%face_z) .and. all(sw2%residual == 0._R8P)
pe0 = pose_interpolate(p0, p1, 0._R8P) ; pe1 = pose_interpolate(p0, p1, 1._R8P) ; pm = pose_interpolate(p0, p1, 0.5_R8P)
call pose_log_3d(p0, p1, ax, ang) ; call pose_log_3d(p0, pm, ax, ang_m)
ok(4) = all(pe0%r == p0%r) .and. all(pe0%t == p0%t) .and. all(pe1%r == p1%r) .and. all(pe1%t == p1%t) .and. &
        abs(ang_m - 0.5_R8P * ang) < 1.e-14_R8P .and. abs(ang - 0.3_R8P) < 1.e-14_R8P
print '(A,L1,A,L1)', 'explicit node poses = screw poses bit for bit: ', ok(3), '  interpolation ends exact, mid angle: ', ok(4)
call swept_cells_3d(s, p0, p1, e(7:1:-1), e, e, sw, status=st)
call swept_cells_3d(s, p0, p1, e, e, e, sw, nq=0_I4P, status=st2)
ok(5) = st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function face_error(sw_) result(emax)
   !< Largest |face area - analytic| over every face of the block (translation by d of the unit cube).
   type(swept_cells_3d_t), intent(in) :: sw_  !< Swept geometry.
   real(R8P)                          :: emax !< Largest error.
   integer                            :: i1, i2, i3 !< Face indexes.

   emax = 0._R8P
   do i3 = 1, 6
      do i2 = 1, 6
         do i1 = 1, 7
            emax = max(emax, abs(sw_%face_x(i1, i2, i3) - face_ref(1, e(i1), e(i2), e(i2 + 1), e(i3), e(i3 + 1))))
            emax = max(emax, abs(sw_%face_y(i2, i1, i3) - face_ref(2, e(i1), e(i3), e(i3 + 1), e(i2), e(i2 + 1))))
            emax = max(emax, abs(sw_%face_z(i2, i3, i1) - face_ref(3, e(i1), e(i2), e(i2 + 1), e(i3), e(i3 + 1))))
         enddo
      enddo
   enddo
   endfunction face_error

   function face_ref(ax_, a, lo1, hi1, lo2, hi2) result(avg)
   !< Analytic time average of the fluid area of the face normal to ax_ at a, spans [lo1, hi1] x
   !< [lo2, hi2] along the cyclic in-plane axes (ax_ + 1, ax_ + 2) for x faces, and the arguments as
   !< the caller orders them otherwise (the area is symmetric in the two spans).
   integer,   intent(in) :: ax_                    !< Face normal axis.
   real(R8P), intent(in) :: a, lo1, hi1, lo2, hi2  !< Plane, spans.
   real(R8P)             :: avg                    !< Average.
   real(R8P)             :: b(10), cand(10), sa, sb, h !< Breakpoints.
   real(R8P)             :: dn, d1, d2             !< Motion normal and along the spans.
   integer               :: n, q, p, a1, a2        !< Work.
   real(R8P), parameter  :: GX(3) = [-sqrt(0.6_R8P), 0._R8P, sqrt(0.6_R8P)], GW(3) = [5._R8P, 8._R8P, 5._R8P] / 9._R8P

   select case (ax_)
   case (1) ; a1 = 2 ; a2 = 3
   case (2) ; a1 = 3 ; a2 = 1
   case default ; a1 = 1 ; a2 = 2
   endselect
   dn = d(ax_) ; d1 = d(a1) ; d2 = d(a2)
   cand = [a / dn, (a - 1._R8P) / dn, lo1 / d1, (lo1 - 1._R8P) / d1, hi1 / d1, (hi1 - 1._R8P) / d1, &
           lo2 / d2, (lo2 - 1._R8P) / d2, hi2 / d2, (hi2 - 1._R8P) / d2]
   n = 0
   do q = 1, 10
      if (cand(q) > 0._R8P .and. cand(q) < 1._R8P) then
         n = n + 1 ; b(n) = cand(q)
      endif
   enddo
   do q = 2, n
      do p = q, 2, -1
         if (b(p) < b(p - 1)) b(p - 1:p) = [b(p), b(p - 1)]
      enddo
   enddo
   avg = 0._R8P
   do q = 0, n
      sa = 0._R8P ; if (q > 0) sa = b(q)
      sb = 1._R8P ; if (q < n) sb = b(q + 1)
      h = sb - sa
      do p = 1, 3
         avg = avg + 0.5_R8P * h * GW(p) * farea([a, lo1, hi1, lo2, hi2, dn, d1, d2], 0.5_R8P * (sa + sb) + 0.5_R8P * h * GX(p))
      enddo
   enddo
   endfunction face_ref

   pure function farea(g, sv) result(ar)
   !< Fluid area of the face at instant sv: the unit cube at s d covers the face plane when
   !< 0 + sv dn < a < 1 + sv dn, over the product of its overlaps with the two spans.
   real(R8P), intent(in) :: g(8), sv !< (a, lo1, hi1, lo2, hi2, dn, d1, d2), instant.
   real(R8P)             :: ar       !< Fluid area.

   ar = (g(3) - g(2)) * (g(5) - g(4))
   if (g(1) > sv * g(6) .and. g(1) < 1._R8P + sv * g(6)) &
      ar = ar - max(0._R8P, min(g(3), 1._R8P + sv * g(7)) - max(g(2), sv * g(7))) * &
                max(0._R8P, min(g(5), 1._R8P + sv * g(8)) - max(g(4), sv * g(8)))
   endfunction farea
endprogram fossil_test_swept_3d
