!< FOSSIL test-side module: analytic distance functions for the conservative advancement tests (issue #73).
module fossil_ccd_test_distances
!< FOSSIL test-side module: analytic distance functions for the conservative advancement tests (issue #73).
!<
!< d(s) = max(0, |a (s - c)| + b (s - c)^2 + e) with the Lipschitz bound of each case given by the test; every call is
!< counted.
use fossil, only : ccd_evaluator_t
use penf,   only : I4P, R8P

implicit none
private
public :: analytic_distance_t

type, extends(ccd_evaluator_t) :: analytic_distance_t
   !< d(s) = max(0, a |s - c| + b (s - c)^2 + e).
   real(R8P)    :: a = 0._R8P   !< Linear (V-shaped) coefficient.
   real(R8P)    :: b = 0._R8P   !< Quadratic coefficient.
   real(R8P)    :: c = 0._R8P   !< Centre.
   real(R8P)    :: e = 0._R8P   !< Offset.
   integer(I4P) :: calls = 0_I4P !< Calls.
   contains
      procedure, pass(self) :: distance !< d(s).
      procedure, pass(self) :: exact    !< d(s), not counted.
endtype analytic_distance_t

contains
   subroutine distance(self, s, d)
   !< d(s), counted.
   class(analytic_distance_t), intent(inout) :: self !< Function.
   real(R8P),                  intent(in)    :: s    !< Parameter.
   real(R8P),                  intent(out)   :: d    !< Distance.

   self%calls = self%calls + 1_I4P
   d = self%exact(s)
   endsubroutine distance

   pure function exact(self, s) result(d)
   !< d(s).
   class(analytic_distance_t), intent(in) :: self !< Function.
   real(R8P),                  intent(in) :: s    !< Parameter.
   real(R8P)                              :: d    !< Distance.

   d = max(0._R8P, self%a * abs(s - self%c) + self%b * (s - self%c)**2 + self%e)
   endfunction exact
endmodule fossil_ccd_test_distances
