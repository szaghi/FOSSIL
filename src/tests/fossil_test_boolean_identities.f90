!< FOSSIL test: boolean identities on the box cases (issue #63, roadmap #58 WP5, step 1).
!<
!< The harness `fossil_boolean_harness` (volume, second moments, area, deviation from the inputs, membership, validity)
!< on the box cases validated so far against analytic volumes (issue #33):
!< 1. 3 of the 6 coplanar cases of fossil_test_boolean_coplanar (offset xy, identical, box on face plane), rotated
!<    (axis (1, 2, 3), angle 0.7);
!< 2. the second seeded configuration of fossil_test_boolean_stress (same generator, same seed; rotated).
!< Regression subset (the 2 s budget); every case, axis-aligned too, is src/tests/validation/fossil_val_boolean_identities.f90.
!< Every pair must pass (nothing new found on these). Area identity only without coplanar face planes; validity only
!< without contact; self-intersections only as rounding folds (depth <= 1e-12 of the diagonal, issue #65).
!< 3. Negative control: the union of the generic case with one facet dropped fails the harness (volume, moments,
!<    deviation or validity), so the harness detects a wrong result.
program fossil_test_boolean_identities

use fossil,                 only : surface_stl_object
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

integer(I4P), parameter :: N_CASES = 2
real(R8P),    parameter :: APART_S(3, 4) = reshape([0.5_R8P, 1._R8P, 1._R8P,   1._R8P, 1._R8P, 1._R8P, &
                                                    1._R8P, 1.5_R8P, 0.5_R8P,  0.5_R8P, 0.5_R8P, 0.5_R8P], [3, 4])
real(R8P),    parameter :: APART_T(3, 4) = reshape([1.25_R8P, 0._R8P, 0._R8P,  0._R8P, -1.5_R8P, 0._R8P, &
                                                    0._R8P, 0._R8P, 1.5_R8P,   1.25_R8P, 1.25_R8P, 1.25_R8P], [3, 4])
type(boolean_report_t)  :: rep
real(R8P)               :: s(3), t(3), axis(3), angle, ov(3)
integer(I8P)            :: seed
integer(I4P)            :: c, k, npass, nall
logical                 :: rotated, contact, coplanar, ok, all_ok

all_ok = .true.

! ---- 1. coplanar cases
npass = 0 ; nall = 0
call pair('offset xy 0.5     ', [0.5_R8P, 0.5_R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], .true., .false., .true.)
call pair('identical         ', [0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], .true., .false., .true.)
call pair('box on face plane ', [0.5_R8P, 0.25_R8P, 0._R8P], [1._R8P, 0.5_R8P, 0.5_R8P], .true., .false., .true.)
ok = npass == nall
print '(A,I0,A,I0,A,L1)', '1. coplanar cases: ', npass, ' of ', nall, ' pass  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. stress configurations
npass = 0 ; nall = 0
seed = 20261102_I8P
do c = 1, N_CASES
   if (c <= N_CASES) then
      do k = 1, 3
         s(k) = 0.5_R8P * real(1 + int(3 * rnd()), R8P)
         t(k) = -0.5_R8P + 0.25_R8P * real(int(7 * rnd()), R8P)
      enddo
   else
      s = APART_S(:, c - N_CASES) ; t = APART_T(:, c - N_CASES)
   endif
   rotated = mod(c, 2) == 0
   axis  = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
   angle = 6.283185307179586_R8P * rnd()
   ov = min(1._R8P, t + s) - max(0._R8P, t)
   contact = any(ov == 0._R8P) .and. .not. any(ov < 0._R8P)
   coplanar = any(t == 0._R8P) .or. any(t == 1._R8P) .or. any(t + s == 0._R8P) .or. any(t + s == 1._R8P)
   if (c == 2) call pair('stress', t, s, rotated, contact, coplanar, axis, angle)
enddo
ok = npass == nall
print '(A,I0,A,I0,A,L1)', '2. stress configurations: ', npass, ' of ', nall, ' pass  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. negative control
block
   type(surface_stl_object) :: a, b
   call boxes(a, b, [0.5_R8P, 0.5_R8P, 0.5_R8P], [1._R8P, 1._R8P, 1._R8P], .false.)
   call check_boolean(a, b, rep, 100_I4P, 7_I8P, .true., .false., corrupt=1_I4P)
   call print_report('   corrupted union  ', rep)
   ok = .not. boolean_pass(rep)
endblock
print '(A,L1)', '3. a union with a facet dropped fails the harness: ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine pair(label, tt, ss, rot, touch, copl, ax, an)
   !< One pair: unit cube and the box ss at tt, optionally rotated; harness, report, count.
   character(*),        intent(in) :: label
   real(R8P),           intent(in) :: tt(3), ss(3)
   logical,             intent(in) :: rot, touch, copl
   real(R8P), optional, intent(in) :: ax(3), an
   type(surface_stl_object)        :: a, b

   call boxes(a, b, tt, ss, rot, ax, an)
   call check_boolean(a, b, rep, 40_I4P, 20261106_I8P, .not. copl, touch)
   nall = nall + 1
   if (boolean_pass(rep)) npass = npass + 1
   call print_report('   '//merge('rot ', '    ', rot)//label, rep)
   endsubroutine pair

   subroutine boxes(a, b, tt, ss, rot, ax, an)
   !< A = unit cube (cube.stl), B = the box of sizes ss at tt; both rotated by (ax, an) (default (1, 2, 3), 0.7).
   type(surface_stl_object), intent(inout) :: a, b
   real(R8P),                intent(in)    :: tt(3), ss(3)
   logical,                  intent(in)    :: rot
   real(R8P), optional,      intent(in)    :: ax(3), an
   real(R8P)                               :: x(3), g

   x = [1._R8P, 2._R8P, 3._R8P] ; g = 0.7_R8P
   if (present(ax)) x = ax
   if (present(an)) g = an
   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%resize(x=ss(1), y=ss(2), z=ss(3))
   call b%translate(x=tt(1), y=tt(2), z=tt(3))
   if (rot) then
      call a%rotate(axis=vector_R8P(x(1), x(2), x(3)), angle=g)
      call b%rotate(axis=vector_R8P(x(1), x(2), x(3)), angle=g)
   endif
   call a%analyze ; call b%analyze
   endsubroutine boxes

   function rnd() result(r)
   !< Park-Miller in (0, 1), the generator of fossil_test_boolean_stress.
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_boolean_identities
