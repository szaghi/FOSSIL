!< FOSSIL test: conservative advancement on analytic distance functions (issue #73, step 3).
program fossil_test_ccd_advance
!< FOSSIL test: conservative advancement on analytic distance functions (issue #73, step 3).
!<
!< d(s) = max(0, a |s - c| + b (s - c)^2 + e) with its Lipschitz bound mu; every result is certified by sampling d on
!< [0, toi] (10^5 points: d >= tol / 2 there, so no earlier contact):
!< 1. a linear approach at the bound's speed (a = 2, mu = 2): one step, toi <= 0.5, d(toi) <= tol;
!< 2. the same with mu = 4 (slack): toi <= 0.5, steps within ceiling(2 mu / tol), about log2(1 / tol);
!< 3. a parabola touching at s = 0.6 (grazing contact): CCD_CONTACT for tol = 1e-4; with tol = 1e-6 the default cap of
!<    1000 steps is hit (steps ~ 1 / sqrt(tol)): CCD_NOT_CONVERGED, toi a certified lower bound below 0.6;
!< 4. a near miss (d_min = 0.01 > tol = 1e-3): CCD_SEPARATE, toi = 1;
!< 5. tunnelling: a fast V (a = 10) far apart at s = 0 and s = 1 (the end-pose check sees no contact): CCD_CONTACT at
!<    0.3 in one step;
!< 6. d(0) = 0: CCD_INITIAL_CONTACT, toi = 0, no step; mu = 0: CCD_SEPARATE with one call; tolerance 0 or a cap of 0:
!<    CCD_INVALID_INPUT, no call;
!< 7. negative control: a bound below the true speed (mu = 1 for a = 2) misses the contact of case 1.

use fossil,                    only : ccd_advance, CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT, CCD_NOT_CONVERGED, &
                                      CCD_INVALID_INPUT
use fossil_ccd_test_distances, only : analytic_distance_t
use penf,                      only : I4P, R8P

implicit none

type(analytic_distance_t) :: f
real(R8P)                 :: toi, d
integer(I4P)              :: st, it
logical                   :: all_ok

all_ok = .true.

! 1. linear approach at the bound's speed
f = analytic_distance_t(a=2._R8P, c=0.5_R8P)
call ccd_advance(f, 2._R8P, 1e-6_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. it == 1 .and. toi <= 0.5_R8P .and. d <= 1e-6_R8P .and. certified(f, toi, 1e-6_R8P), &
           'linear approach, mu = speed: contact in one step, toi <= 0.5, d <= tol, certified (toi)', toi)

! 2. slack bound
f = analytic_distance_t(a=2._R8P, c=0.5_R8P)
call ccd_advance(f, 4._R8P, 1e-6_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. toi <= 0.5_R8P .and. d <= 1e-6_R8P .and. it <= ceiling(2 * 4._R8P / 1e-6_R8P) .and. &
           it <= 30 .and. certified(f, toi, 1e-6_R8P), 'linear approach, mu = 2 speed: contact, certified (steps)', &
           real(it, R8P))

! 3. grazing contact
f = analytic_distance_t(b=4._R8P, c=0.6_R8P)
call ccd_advance(f, 4.8_R8P, 1e-4_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. toi <= 0.6_R8P .and. d <= 1e-4_R8P .and. certified(f, toi, 1e-4_R8P), &
           'grazing contact, tol 1e-4: contact, certified (steps)', real(it, R8P))
f = analytic_distance_t(b=4._R8P, c=0.6_R8P)
call ccd_advance(f, 4.8_R8P, 1e-6_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_NOT_CONVERGED .and. it == 1000 .and. toi < 0.6_R8P .and. d > 1e-6_R8P .and. &
           certified(f, toi, 1e-6_R8P), 'grazing contact, tol 1e-6: cap hit, toi a certified lower bound (toi)', toi)

! 4. near miss
f = analytic_distance_t(b=4._R8P, c=0.6_R8P, e=0.01_R8P)
call ccd_advance(f, 4.8_R8P, 1e-3_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_SEPARATE .and. toi == 1._R8P .and. certified(f, toi, 1e-3_R8P), &
           'near miss d_min = 0.01: separate, toi = 1 (steps)', real(it, R8P))

! 5. tunnelling
f = analytic_distance_t(a=10._R8P, c=0.3_R8P)
call check(f%exact(0._R8P) > 1._R8P .and. f%exact(1._R8P) > 1._R8P, 'tunnelling: the end poses are far apart (d(0))', &
           f%exact(0._R8P))
call ccd_advance(f, 10._R8P, 1e-6_R8P, toi, st, distance=d, iterations=it)
call check(st == CCD_CONTACT .and. it == 1 .and. toi <= 0.3_R8P .and. certified(f, toi, 1e-6_R8P), &
           'tunnelling: contact found in one step (toi)', toi)

! 6. degenerate inputs
f = analytic_distance_t(a=1._R8P)
call ccd_advance(f, 1._R8P, 1e-6_R8P, toi, st, iterations=it)
call check(st == CCD_INITIAL_CONTACT .and. toi == 0._R8P .and. it == 0, 'd(0) = 0: initial contact, toi = 0', toi)
f = analytic_distance_t(e=1._R8P)
call ccd_advance(f, 0._R8P, 1e-6_R8P, toi, st, iterations=it)
call check(st == CCD_SEPARATE .and. toi == 1._R8P .and. f%calls == 1, 'mu = 0: separate with one call (calls)', &
           real(f%calls, R8P))
f = analytic_distance_t(e=1._R8P)
call ccd_advance(f, 1._R8P, 0._R8P, toi, st)
call check(st == CCD_INVALID_INPUT .and. f%calls == 0, 'tolerance 0: invalid input, no call', real(st, R8P))
call ccd_advance(f, 1._R8P, 1e-6_R8P, toi, st, max_iterations=0_I4P)
call check(st == CCD_INVALID_INPUT .and. f%calls == 0, 'cap 0: invalid input, no call', real(st, R8P))

! 7. negative control
f = analytic_distance_t(a=2._R8P, c=0.5_R8P)
call ccd_advance(f, 1._R8P, 1e-6_R8P, toi, st)
call check(.not. (st == CCD_CONTACT .and. toi <= 0.5_R8P), 'negative control: a bound below the speed misses (status)', &
           real(st, R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES23.15)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   logical function certified(g, t, tol) result(ok)
   !< d >= tol / 2 at 10^5 points of [0, t) (no earlier contact).
   type(analytic_distance_t), intent(in) :: g   !< Function.
   real(R8P),                 intent(in) :: t   !< toi.
   real(R8P),                 intent(in) :: tol !< Tolerance.
   integer(I4P)                          :: j   !< Counter.

   ok = .true.
   do j = 0, 99999
      ok = ok .and. g%exact(t * j / 100000._R8P) >= 0.5_R8P * tol * (1._R8P - 1e-9_R8P)
   enddo
   endfunction certified
endprogram fossil_test_ccd_advance
