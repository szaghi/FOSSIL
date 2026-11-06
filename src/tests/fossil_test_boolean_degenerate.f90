!< FOSSIL test: boolean identities on near-degenerate configurations (issue #63, roadmap #58 WP5, step 3).
!<
!< One case per family of the validation sweep (src/tests/validation/fossil_val_boolean_degenerate.f90), each checked
!< by the harness (fossil_boolean_harness) on the four ops; every case must pass (a non-OK status is allowed):
!< - the unit cube against itself rotated by 1e-4 and by 1e-8 rad about a generic axis: the first cracked (the
!<   retriangulation lost a facet edge with a cut 3.5e-5 inside it), the second lost 2 % of the volume (a folded
!<   Bowyer-Watson cavity), both with status OK.
!< The third defect of the step (shared sub-triangles misdetected beside crossings) needs a curved body: its case, a
!< 24^3 sphere against itself translated by 1e-6, is gated in the validation sweep (2.2 s for a 12^3 sphere in the
!< debug coverage build, over the budget). The snap, nested and contact families pass before and after the fix: they
!< are in the validation sweep and in fossil_test_boolean_coplanar / fossil_test_boolean_stress.
program fossil_test_boolean_degenerate

use fossil,                 only : surface_stl_object
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

type(surface_stl_object) :: a, b
type(vector_R8P)         :: axis
integer(I4P)             :: npass, nall

npass = 0 ; nall = 0
axis = vector_R8P(0.3_R8P, -0.5_R8P, 0.8_R8P) ; axis = axis / axis%normL2()
call cube(a)
call cube(b) ; call b%rotate(axis=axis, angle=1e-4_R8P) ; call b%analyze
call run('cube rotated 1e-4', a, b, .true., .false.)
call cube(b) ; call b%rotate(axis=axis, angle=1e-8_R8P) ; call b%analyze
call run('cube rotated 1e-8', a, b, .true., .false.)
print '(A,I0,A,I0)', 'cases passing: ', npass, ' of ', nall
print '(A,L1)', 'Are all tests passed? ', npass == nall

contains
   subroutine cube(s)
   !< The unit cube centred at the origin.
   type(surface_stl_object), intent(inout) :: s     !< Cube.
   type(surface_stl_object)                :: fresh !< Reset.

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
   call s%analyze
   endsubroutine cube

   subroutine run(lab, a, b, area_identity, contact)
   !< Check one case, print it, count it.
   character(*),             intent(in) :: lab           !< Case name.
   type(surface_stl_object), intent(in) :: a, b          !< Inputs.
   logical,                  intent(in) :: area_identity !< Check the area identity.
   logical,                  intent(in) :: contact       !< The inputs touch.
   type(boolean_report_t)               :: rep           !< Report.
   character(30)                        :: lab30         !< Padded label.

   call check_boolean(a, b, rep, 40_I4P, 20261106_I8P, area_identity, contact)
   lab30 = lab
   call print_report('   '//lab30, rep)
   nall = nall + 1 ; if (boolean_pass(rep)) npass = npass + 1
   endsubroutine run
endprogram fossil_test_boolean_degenerate
