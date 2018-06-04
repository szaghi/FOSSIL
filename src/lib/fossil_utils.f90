!< FOSSIL, utils library.

module fossil_utils
!< FOSSIL, utils library.

use penf, only : I4P, R4P, R8P, ZeroR4P
use vecfor, only : vector_R8P

implicit none
private
public :: EPS
public :: FRLEN
public :: PI
public :: is_inside_bb

real(R8P), protected    :: EPS=real(ZeroR4P, R8P)     !< Small EPSILON to avoid rund off errors.
integer(I4P), parameter :: FRLEN=80                   !< Maximum length of fossil STL label string.
real(R8P),    parameter :: PI = 4._R8P * atan(1._R8P) !< Pi greek.

contains
   pure function is_inside_bb(bmin, bmax, point)
   !< Return the true if point is inside AABB.
   type(vector_R8P), intent(in) :: bmin, bmax   !< Bounding box extents.
   type(vector_R8P), intent(in) :: point        !< Point reference.
   logical                      :: is_inside_bb !< Check result.

   is_inside_bb = ((point%x >= bmin%x.and.point%x <= bmax%x).and.&
                   (point%y >= bmin%y.and.point%y <= bmax%y).and.&
                   (point%z >= bmin%z.and.point%z <= bmax%z))
   endfunction is_inside_bb
endmodule fossil_utils
