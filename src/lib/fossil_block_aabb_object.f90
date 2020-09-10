!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition for cartesian block.

module fossil_block_aabb_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) class definition for cartesian block.

use penf, only : FR8P, I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, normL2_R8P, vector_R8P

implicit none
private
public :: aabb_object

type :: aabb_object
   !< AABB refinment class definition for cartesian block.
   type(vector_R8P) :: bmin !< Minimum point of AABB.
   type(vector_R8P) :: bmax !< Maximum point of AABB.
   contains
      ! public methods
      procedure, pass(self) :: compute_octants !< Compute AABB octants.
      procedure, pass(self) :: is_inside       !< Return the true if point is inside ABB.
      procedure, pass(self) :: vertex          !< Return AABB vertices.
      ! operators
      generic :: assignment(=) => aabb_assign_aabb      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_assign_aabb !< Operator `=`.
endtype aabb_object

contains
   ! public methods
   pure subroutine compute_octants(self, octant)
   !< Return AABB octants.
   class(aabb_object), intent(in)  :: self      !< AABB.
   type(aabb_object),  intent(out) :: octant(8) !< AABB octants.
   type(vector_R8P)                :: vertex(8) !< AABB vertices.
   integer(I4P)                    :: o         !< Counter.

   vertex = self%vertex()
   octant(1)%bmin = self%bmin      ; octant(1)%bmax = 0.5_R8P * (self%bmin + self%bmax)
   octant(8)%bmin = octant(1)%bmax ; octant(8)%bmax = self%bmax
   do o=2, 7 ! loop over remaining octants
      octant(o)%bmin = 0.5_R8P * (self%bmin + vertex(o)) ; octant(o)%bmax = 0.5_R8P * (vertex(o) + self%bmax)
   enddo
   endsubroutine compute_octants

   pure function is_inside(self, point)
   !< Return the true if point is inside ABB.
   class(aabb_object), intent(in) :: self      !< AABB.
   type(vector_R8P),   intent(in) :: point     !< Point reference.
   logical                        :: is_inside !< Check result.

   is_inside = ((point%x >= self%bmin%x.and.point%x <= self%bmax%x).and.&
                (point%y >= self%bmin%y.and.point%y <= self%bmax%y).and.&
                (point%z >= self%bmin%z.and.point%z <= self%bmax%z))
   endfunction is_inside

   pure function vertex(self)
   !< Return AABB vertices.
   class(aabb_object), intent(in) :: self      !< AABB.
   type(vector_R8P)               :: vertex(8) !< AABB vertices.

   vertex(1) = self%bmin
   vertex(2) = self%bmax%x * ex_R8P + self%bmin%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(3) = self%bmin%x * ex_R8P + self%bmax%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(4) = self%bmax%x * ex_R8P + self%bmax%y * ey_R8P + self%bmin%z * ez_R8P
   vertex(5) = self%bmin%x * ex_R8P + self%bmin%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(6) = self%bmax%x * ex_R8P + self%bmin%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(7) = self%bmin%x * ex_R8P + self%bmax%y * ey_R8P + self%bmax%z * ez_R8P
   vertex(8) = self%bmax
   endfunction vertex

   ! operators
   ! =
   pure subroutine aabb_assign_aabb(lhs, rhs)
   !< Operator `=`.
   class(aabb_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_object),  intent(in)    :: rhs !< Right hand side.

   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
   endsubroutine aabb_assign_aabb
endmodule fossil_block_aabb_object
