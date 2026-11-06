!< FOSSIL test: boolean identities on curved bodies (issue #63, roadmap #58 WP5, step 2).
!<
!< Marching-cubes bodies (fossil_boolean_bodies: sphere, torus, capped cylinder, 10^3 samples) and the unit cube, in
!< pairs, B moved by a fixed-seed proper rigid motion (Park-Miller, random axis and angle, translation components within +-0.2)
!< relative to A. The harness checks every identity of the four ops (fossil_boolean_harness); every pair must pass.
!< The full sweep (finer bodies, more motions, the alpha-wrapped bunny and dragon) is
!< src/tests/validation/fossil_val_boolean_curved.f90.
program fossil_test_boolean_curved

use fossil,                 only : surface_stl_object
use fossil_boolean_bodies,  only : BODY_SPHERE, BODY_TORUS, BODY_CYLINDER, mc_body, move_body
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

integer(I4P), parameter  :: N = 10
integer(I4P), parameter  :: KA(2) = [BODY_SPHERE, BODY_CYLINDER]
integer(I4P), parameter  :: KB(2) = [BODY_TORUS, 0_I4P]
character(8), parameter  :: NAMES(0:3) = ['cube    ', 'sphere  ', 'torus   ', 'cylinder']
type(surface_stl_object) :: a, b
type(boolean_report_t)   :: rep
real(R8P)                :: axis(3), angle, t(3)
integer(I8P)             :: seed
integer(I4P)             :: c, npass, nall
character(40)            :: label

seed = 20261106_I8P ; npass = 0 ; nall = 0
do c = 1, size(KA)
   call body(KA(c), a) ; call body(KB(c), b)
   axis = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] ; angle = 6.283185307179586_R8P * rnd()
   t = 0.4_R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
   call move_body(b, axis, angle, t)
   call check_boolean(a, b, rep, 50_I4P, seed, .true., .false.)
   write(label, '(3A)') trim(NAMES(KA(c))), ' x ', trim(NAMES(KB(c)))
   call print_report('   '//label(1:24), rep)
   nall = nall + 1 ; if (boolean_pass(rep)) npass = npass + 1
enddo
print '(A,I0,A,I0)', 'pairs passing: ', npass, ' of ', nall
print '(A,L1)', 'Are all tests passed? ', npass == nall

contains
   subroutine body(kind, s)
   !< A marching-cubes body, or the unit cube centred at the origin (kind 0).
   integer(I4P),             intent(in)    :: kind
   type(surface_stl_object), intent(inout) :: s

   if (kind == 0) then
      call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
      call s%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P)
      call s%analyze
   else
      call mc_body(kind, N, s)
   endif
   endsubroutine body

   function rnd() result(r)
   !< Park-Miller in (0, 1).
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_boolean_curved
