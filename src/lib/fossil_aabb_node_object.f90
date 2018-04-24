!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.

module fossil_aabb_node_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.
!<
!< This is just a *container* for AABB tree's nodes.

use fossil_aabb_object, only : aabb_object
use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none
private
public :: aabb_node_object

type :: aabb_node_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree-node class.
   type(aabb_object), allocatable :: aabb !< AABB data.
   contains
      ! public methods
      procedure, pass(self) :: closest_point    !< Return closest point on AABB from point reference.
      procedure, pass(self) :: destroy          !< Destroy AABB.
      procedure, pass(self) :: distance         !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: do_ray_intersect !< Return true if AABB is intersected by ray.
      procedure, pass(self) :: initialize       !< Initialize AABB.
      ! operators
      generic :: assignment(=) => aabb_node_assign_aabb_node      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_node_assign_aabb_node !< Operator `=`.
endtype aabb_node_object

contains
   ! public methods
   pure function closest_point(self, point) result(closest)
   !< Return closest point on (or in) AABB from point reference.
   class(aabb_node_object), intent(in) :: self    !< AABB box.
   type(vector_R8P),        intent(in) :: point   !< Point reference.
   type(vector_R8P)                    :: closest !< Closest point on (on in) aabb to point.

   closest = MaxR8P
   if (allocated(self%aabb)) closest = self%aabb%closest_point(point=point)
   endfunction closest_point

   elemental subroutine destroy(self)
   !< Destroy AABB.
   class(aabb_node_object), intent(inout) :: self  !< AABB.
   type(aabb_node_object)                 :: fresh !< Fresh instance of AABB.

   self = fresh
   endsubroutine destroy

   pure function distance(self, point)
   !< Return the (square) distance from point to AABB.
   class(aabb_node_object), intent(in) :: self     !< AABB.
   type(vector_R8P),        intent(in) :: point    !< Point reference.
   real(R8P)                           :: distance !< Distance from point to AABB.

   distance = 0._R8P
   if (allocated(self%aabb)) distance = self%aabb%distance(point=point)
   endfunction distance

   pure function do_ray_intersect(self, ray_origin, ray_direction) result(do_intersect)
   !< Return true if AABB is intersected by ray from origin and oriented as ray direction vector.
   class(aabb_node_object), intent(in) :: self          !< AABB.
   type(vector_R8P),        intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction !< Ray direction.
   logical                             :: do_intersect  !< Test result.

   do_intersect = .false.
   if (allocated(self%aabb)) do_intersect = self%aabb%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)
   endfunction do_ray_intersect

   pure subroutine initialize(self, facet, bmin, bmax)
   !< Initialize AABB.
   class(aabb_node_object), intent(inout)        :: self     !< AABB.
   type(facet_object),      intent(in), optional :: facet(:) !< Facets list.
   type(vector_R8P),        intent(in), optional :: bmin     !< Minimum point of AABB.
   type(vector_R8P),        intent(in), optional :: bmax     !< Maximum point of AABB.

   call self%destroy
   if (present(facet).or.(present(bmin).and.present(bmin))) then
      allocate(self%aabb)
      call self%aabb%initialize(facet=facet, bmin=bmin, bmax=bmax)
   endif
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine aabb_node_assign_aabb_node(lhs, rhs)
   !< Operator `=`.
   class(aabb_node_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_node_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(lhs%aabb)) then
      call lhs%aabb%destroy
      deallocate(lhs%aabb)
   endif
   if (allocated(rhs%aabb)) then
      allocate(lhs%aabb)
      lhs%aabb = rhs%aabb
   endif
   endsubroutine aabb_node_assign_aabb_node
endmodule fossil_aabb_node_object
