!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.

module fossil_aabb_node_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree-node class definition.
!<
!< This is just a *container* for AABB tree's nodes.

use fossil_aabb_object, only : aabb_object
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : vector_R8P

implicit none
private
public :: aabb_node_object

type :: aabb_node_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree-node class.
   private
   type(aabb_object), allocatable :: aabb !< AABB data.
   contains
      ! public methods
      procedure, pass(self) :: add_facets                  !< Add facets to AABB.
      procedure, pass(self) :: bmin                        !< Return AABB bmin.
      procedure, pass(self) :: bmax                        !< Return AABB bmax.
      procedure, pass(self) :: closest_point               !< Return closest point on AABB from point reference.
      procedure, pass(self) :: compute_octants             !< Compute AABB octants.
      procedure, pass(self) :: compute_vertices_nearby     !< Compute vertices nearby.
      procedure, pass(self) :: destroy                     !< Destroy AABB.
      procedure, pass(self) :: distance                    !< Return the (square) distance from point to AABB.
      procedure, pass(self) :: distance_from_facets        !< Return the (square) distance from point to AABB's facets.
      procedure, pass(self) :: do_ray_intersect            !< Return true if AABB is intersected by ray.
      procedure, pass(self) :: get_aabb_facets             !< Get AABB facets list.
      procedure, pass(self) :: has_facets                  !< Return true if AABB has facets.
      procedure, pass(self) :: initialize                  !< Initialize AABB.
      procedure, pass(self) :: is_allocated                !< Return true is node is allocated.
      procedure, pass(self) :: ray_intersections_number    !< Return ray intersections number.
      procedure, pass(self) :: save_geometry_tecplot_ascii !< Save AABB geometry into Tecplot ascii file.
      procedure, pass(self) :: save_facets_into_file_stl   !< Save facets into file STL.
      procedure, pass(self) :: translate                   !< Translate AABB by delta.
      procedure, pass(self) :: union                       !< Make AABB the union of other AABBs.
      procedure, pass(self) :: update_extents              !< Update AABB bounding box extents.
      ! operators
      generic :: assignment(=) => aabb_node_assign_aabb_node      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_node_assign_aabb_node !< Operator `=`.
endtype aabb_node_object

contains
   ! public methods
   pure subroutine add_facets(self, facet_id, facet, is_exclusive)
   !< Add facets to AABB.
   !<
   !< @note Facets added to AABB are removed to facets list that is also returned.
   class(aabb_node_object), intent(inout)        :: self         !< AABB.
   type(list_id_object),    intent(inout)        :: facet_id     !< List of facets IDs.
   type(facet_object),      intent(in)           :: facet(:)     !< Facets list.
   logical,                 intent(in), optional :: is_exclusive !< Sentinel to enable/disable exclusive addition.

   if (allocated(self%aabb)) call self%aabb%add_facets(facet_id=facet_id, facet=facet, is_exclusive=is_exclusive)
   endsubroutine add_facets

   pure function bmin(self)
   !< Return AABB bmin.
   class(aabb_node_object), intent(in) :: self !< AABB box.
   type(vector_R8P)                    :: bmin !< AABB bmin.

   if (allocated(self%aabb)) bmin = self%aabb%bmin
   endfunction bmin

   pure function bmax(self)
   !< Return AABB bmax.
   class(aabb_node_object), intent(in) :: self !< AABB box.
   type(vector_R8P)                    :: bmax !< AABB bmax.

   if (allocated(self%aabb)) bmax = self%aabb%bmax
   endfunction bmax

   pure function closest_point(self, point) result(closest)
   !< Return closest point on (or in) AABB from point reference.
   class(aabb_node_object), intent(in) :: self    !< AABB box.
   type(vector_R8P),        intent(in) :: point   !< Point reference.
   type(vector_R8P)                    :: closest !< Closest point on (on in) aabb to point.

   closest = MaxR8P
   if (allocated(self%aabb)) closest = self%aabb%closest_point(point=point)
   endfunction closest_point

   pure subroutine compute_octants(self, octant)
   !< Return AABB octants.
   class(aabb_node_object), intent(in)  :: self      !< AABB.
   type(aabb_object),       intent(out) :: octant(8) !< AABB octants.
   type(vector_R8P)                     :: vertex(8) !< AABB vertices.
   integer(I4P)                         :: o         !< Counter.

   call self%aabb%compute_octants(octant=octant)
   endsubroutine compute_octants

   pure subroutine compute_vertices_nearby(self, facet, tolerance_to_be_identical, tolerance_to_be_nearby)
   !< Compute vertices nearby.
   class(aabb_node_object), intent(in)    :: self                      !< AABB.
   type(facet_object),      intent(inout) :: facet(:)                  !< Facets list.
   real(R8P),               intent(in)    :: tolerance_to_be_identical !< Tolerance to identify identical vertices.
   real(R8P),               intent(in)    :: tolerance_to_be_nearby    !< Tolerance to identify nearby vertices.

   if (allocated(self%aabb)) call self%aabb%compute_vertices_nearby(facet=facet,                                         &
                                                                    tolerance_to_be_identical=tolerance_to_be_identical, &
                                                                    tolerance_to_be_nearby=tolerance_to_be_nearby)
   endsubroutine compute_vertices_nearby

   elemental subroutine destroy(self)
   !< Destroy AABB.
   class(aabb_node_object), intent(inout) :: self  !< AABB.
   type(aabb_node_object)                 :: fresh !< Fresh instance of AABB.

   if (allocated(self%aabb)) then
      call self%aabb%destroy
      deallocate(self%aabb)
   endif
   endsubroutine destroy

   pure function distance(self, point)
   !< Return the (square) distance from point to AABB.
   class(aabb_node_object), intent(in) :: self     !< AABB.
   type(vector_R8P),        intent(in) :: point    !< Point reference.
   real(R8P)                           :: distance !< Distance from point to AABB.

   distance = MaxR8P
   if (allocated(self%aabb)) distance = self%aabb%distance(point=point)
   endfunction distance

   pure function distance_from_facets(self, facet, point) result(distance)
   !< Return the (square) distance from point to AABB's facets.
   class(aabb_node_object), intent(in) :: self      !< AABB.
   type(facet_object),      intent(in) :: facet(:)  !< Facets list.
   type(vector_R8P),        intent(in) :: point     !< Point reference.
   real(R8P)                           :: distance  !< Distance from point to AABB's facets.

   distance = MaxR8P
   if (allocated(self%aabb)) distance = self%aabb%distance_from_facets(facet=facet, point=point)
   endfunction distance_from_facets

   pure function do_ray_intersect(self, ray_origin, ray_direction) result(do_intersect)
   !< Return true if AABB is intersected by ray from origin and oriented as ray direction vector.
   class(aabb_node_object), intent(in) :: self          !< AABB.
   type(vector_R8P),        intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction !< Ray direction.
   logical                             :: do_intersect  !< Test result.

   do_intersect = .false.
   if (allocated(self%aabb)) do_intersect = self%aabb%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)
   endfunction do_ray_intersect

   pure subroutine get_aabb_facets(self, facet, aabb_facet)
   !< Get AABB facets list.
   class(aabb_node_object), intent(in)               :: self          !< AABB.
   type(facet_object),      intent(in)               :: facet(:)      !< Whole facets list.
   type(facet_object),      intent(out), allocatable :: aabb_facet(:) !< AABB facets list.

   if (allocated(self%aabb)) call self%aabb%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
   endsubroutine get_aabb_facets

   pure function has_facets(self)
   !< Return true if AABB has facets.
   class(aabb_node_object), intent(in) :: self       !< AABB box.
   logical                             :: has_facets !< Check result.

   has_facets = allocated(self%aabb)
   if (has_facets) has_facets = self%aabb%has_facets()
   endfunction has_facets

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

   pure function is_allocated(self)
   !< Return true if node is allocated.
   class(aabb_node_object), intent(in) :: self         !< AABB box.
   logical                             :: is_allocated !< Check result.

   is_allocated = allocated(self%aabb)
   endfunction is_allocated

   pure function ray_intersections_number(self, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number.
   class(aabb_node_object), intent(in) :: self                 !< AABB.
   type(facet_object),      intent(in) :: facet(:)             !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                        :: intersections_number !< Intersection number.

   intersections_number = 0
   if (allocated(self%aabb)) &
      intersections_number = self%aabb%ray_intersections_number(facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
   endfunction ray_intersections_number

   subroutine  save_geometry_tecplot_ascii(self, file_unit, aabb_name)
   !< Save AABB geometry into Tecplot ascii file.
   class(aabb_node_object), intent(in)           :: self       !< AABB.
   integer(I4P),            intent(in)           :: file_unit  !< File unit.
   character(*),            intent(in), optional :: aabb_name  !< Name of AABB.

   if (allocated(self%aabb)) call self%aabb%save_geometry_tecplot_ascii(file_unit=file_unit, aabb_name=aabb_name)
   endsubroutine  save_geometry_tecplot_ascii

   subroutine save_facets_into_file_stl(self, facet, file_name, is_ascii)
   !< Save facets into file STL.
   class(aabb_node_object), intent(in) :: self      !< AABB.
   type(facet_object),      intent(in) :: facet(:)  !< Facets list.
   character(*),            intent(in) :: file_name !< File name.
   logical,                 intent(in) :: is_ascii  !< Sentinel for file format.

   if (allocated(self%aabb)) call self%aabb%save_facets_into_file_stl(facet=facet, file_name=file_name, is_ascii=is_ascii)
   endsubroutine save_facets_into_file_stl

   elemental subroutine translate(self, delta)
   !< Translate AABB by delta.
   class(aabb_node_object), intent(inout) :: self  !< AABB.
   type(vector_R8P),        intent(in)    :: delta !< Delta of translation.

   if (allocated(self%aabb)) call self%aabb%translate(delta=delta)
   endsubroutine translate

   pure subroutine union(self, node, id)
   !< Make AABB the union of other AABBs.
   class(aabb_node_object), intent(inout) :: self    !< AABB.
   type(aabb_node_object),  intent(in)    :: node(:) !< Nodes list.
   integer(I4P),            intent(in)    :: id(:)   !< Nodes ID list.
   integer(I4P)                           :: i       !< Counter.

   call self%destroy
   allocate(self%aabb)
   do i=1, size(id, dim=1)
      if (allocated(node(id(i))%aabb)) call self%aabb%union(other=node(id(i))%aabb)
   enddo
   endsubroutine union

   pure subroutine update_extents(self, facet)
   !< Update AABB bounding box extents.
   class(aabb_node_object), intent(inout) :: self     !< AABB.
   type(facet_object),      intent(in)    :: facet(:) !< Facets list.

   if (allocated(self%aabb)) call self%aabb%update_extents(facet=facet)
   endsubroutine update_extents

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
