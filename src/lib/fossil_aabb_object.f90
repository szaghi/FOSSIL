!< FOSSIL,  Axis-Aligned Bounding Box (AABB) class definition.

module fossil_aabb_object
!< FOSSIL,  Axis-Aligned Bounding Box (AABB) class definition.

use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none
private
public :: aabb_object

type :: aabb_object
   !< FOSSIL  Axis-Aligned Bounding Box (AABB) class.
   type(vector_R8P)               :: bmin    !< Minimum point of AABB.
   type(vector_R8P)               :: bmax    !< Maximum point of AABB.
   type(aabb_object), allocatable :: aabb(:) !< Octree of children AABB.
   contains
      ! public methods
      procedure, pass(self) :: closest_point        !< Return closest point on AABB from point reference.
      procedure, pass(self) :: destroy              !< Destroy file.
      procedure, pass(self) :: distance             !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: do_segment_intersect !< Given segment, return true if the segment intersect AABB.
      procedure, pass(self) :: initialize           !< Initialize file.
      ! operators
      generic :: assignment(=) => aabb_assign_aabb       !< Overload `=`.
      procedure, pass(lhs),  private :: aabb_assign_aabb !< Operator `=`.
endtype aabb_object

contains
   ! public methods
   pure function closest_point(self, point) result(closest)
   !< Return closest point on (or in) AABB from point reference.
   class(aabb_object), intent(in) :: self    !< AABB tree.
   type(vector_R8P),   intent(in) :: point   !< Point reference.
   type(vector_R8P)               :: closest !< Closest point on (on in) aabb to point.

   closest = point
   closest%x = max(closest%x, self%bmin%x) ; closest%x = min(closest%x, self%bmax%x)
   closest%y = max(closest%y, self%bmin%y) ; closest%y = min(closest%y, self%bmax%y)
   closest%z = max(closest%z, self%bmin%z) ; closest%z = min(closest%z, self%bmax%z)
   endfunction closest_point

   elemental subroutine destroy(self)
   !< Destroy file.
   class(aabb_object), intent(inout) :: self  !< AABB tree.
   type(aabb_object)                 :: fresh !< Fresh instance of file STL.

   self = fresh
   endsubroutine destroy

   pure function distance(self, point)
   !< Return the (square) distance from point to AABB.
   class(aabb_object), intent(in) :: self     !< AABB tree.
   type(vector_R8P),   intent(in) :: point    !< Point reference.
   real(R8P)                      :: distance !< Distance from point to aabb.

   distance = 0._R8P
   if (point%x < self%bmin%x) distance = distance + (self%bmin%x - point%x) * (self%bmin%x - point%x)
   if (point%y < self%bmin%y) distance = distance + (self%bmin%y - point%y) * (self%bmin%y - point%y)
   if (point%z < self%bmin%z) distance = distance + (self%bmin%z - point%z) * (self%bmin%z - point%z)
   if (point%x > self%bmax%x) distance = distance + (point%x - self%bmax%x) * (point%x - self%bmax%x)
   if (point%y > self%bmax%y) distance = distance + (point%y - self%bmax%y) * (point%y - self%bmax%y)
   if (point%z > self%bmax%z) distance = distance + (point%z - self%bmax%z) * (point%z - self%bmax%z)
   endfunction distance

   pure function do_segment_intersect(self, point_1, point_2) result(do_intersect)
   !< Given segment point_1->point_2, return true if the segment intersect AABB.
   class(aabb_object), intent(in) :: self               !< AABB tree.
   type(vector_R8P),   intent(in) :: point_1, point_2   !< Segment.
   logical                        :: do_intersect       !< Test result.
   type(vector_R8P)               :: aabb_center        !< AABB center.
   type(vector_R8P)               :: aabb_halflength    !< AABB halflength.
   type(vector_R8P)               :: segment_midpoint   !< Segment midpoint.
   type(vector_R8P)               :: segment_halflength !< Segment halflength.
   real(R8P)                      :: adx, ady, adz      !< Separating axis.
   real(R8P), parameter           :: eps=1e-6_R8P       !< Small epsilon to control round off errors.

   do_intersect = .false.
   aabb_center = 0.5_R8P * (self%bmin + self%bmax)
   aabb_halflength = self%bmax - aabb_center
   segment_midpoint = 0.5_R8P * (point_1 + point_2)
   segment_halflength = point_2 - segment_midpoint
   segment_midpoint = segment_midpoint - aabb_center
   adx = abs(segment_halflength%x)
   if (abs(segment_midpoint%x) > aabb_halflength%x + adx) return
   ady = abs(segment_halflength%y)
   if (abs(segment_midpoint%y) > aabb_halflength%y + ady) return
   adz = abs(segment_halflength%z)
   if (abs(segment_midpoint%z) > aabb_halflength%z + adz) return
   adx = adx + eps
   ady = ady + eps
   adz = adz + eps
   if (abs(segment_midpoint%y * segment_halflength%z - segment_midpoint%z * segment_halflength%y) > &
      aabb_halflength%y * adz + aabb_halflength%z * ady) return
   if (abs(segment_midpoint%z * segment_halflength%x - segment_midpoint%x * segment_halflength%z) > &
      aabb_halflength%x * adz + aabb_halflength%z * adx) return
   if (abs(segment_midpoint%x * segment_halflength%y - segment_midpoint%y * segment_halflength%x) > &
      aabb_halflength%x * ady + aabb_halflength%y * adx) return
   ! no separating axis found
   do_intersect = .true.
   endfunction do_segment_intersect

   pure subroutine initialize(self, facet)
   !< Initialize file.
   class(aabb_object), intent(inout) :: self     !< AABB tree.
   type(facet_object), intent(in)    :: facet(:) !< Facets list.

   call self%destroy
   self%bmin%x = minval(facet(:)%bb(1)%x)
   self%bmin%y = minval(facet(:)%bb(1)%y)
   self%bmin%z = minval(facet(:)%bb(1)%z)
   self%bmax%x = maxval(facet(:)%bb(2)%x)
   self%bmax%y = maxval(facet(:)%bb(2)%y)
   self%bmax%z = maxval(facet(:)%bb(2)%z)
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine aabb_assign_aabb(lhs, rhs)
   !< Operator `=`.
   class(aabb_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_object),  intent(in)    :: rhs !< Right hand side.

   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
   if (allocated(lhs%aabb)) deallocate(lhs%aabb)
   if (allocated(rhs%aabb)) allocate(lhs%aabb(1:size(rhs%aabb, dim=1)), source=rhs%aabb)
   endsubroutine aabb_assign_aabb
endmodule fossil_aabb_object
