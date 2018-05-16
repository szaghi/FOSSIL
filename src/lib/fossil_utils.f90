!< FOSSIL, utils library.

module fossil_utils
!< FOSSIL, utils library.

use penf, only : I4P, R4P, R8P, ZeroR4P

implicit none
private
public :: EPS
public :: FRLEN
public :: PI

real(R8P), protected    :: EPS=real(ZeroR4P, R8P)     !< Small EPSILON to avoid rund off errors.
integer(I4P), parameter :: FRLEN=80                   !< Maximum length of fossil STL label string.
real(R8P),    parameter :: PI = 4._R8P * atan(1._R8P) !< Pi greek.
endmodule fossil_utils
