!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

module fossil_aabb_tree_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

use fossil_aabb_object, only : aabb_object
use fossil_aabb_node_object, only : aabb_node_object
use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P, R8P, MaxR8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none
private
public :: aabb_tree_object

integer(I4P), parameter :: TREE_RATIO=8 !< Tree refinement ratio, it is assumed to be an **octree**.

type :: aabb_tree_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree class.
   !<
   !< @note The tree is assumed to be an **octree**.
   integer(I4P)                        :: refinement_levels=0 !< Total number of refinement levels used.
   integer(I4P)                        :: nodes_number=0      !< Total number of tree nodes.
   type(aabb_node_object), allocatable :: node(:)             !< AABB tree nodes [0:nodes_number-1].

   contains
      ! public methods
      procedure, pass(self) :: destroy                       !< Destroy AABB tree.
      procedure, pass(self) :: initialize                    !< Initialize AABB tree.
      ! operators
      generic :: assignment(=) => aabb_tree_assign_aabb_tree      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_tree_assign_aabb_tree !< Operator `=`.
endtype aabb_tree_object

contains
   ! public methods
   elemental subroutine destroy(self)
   !< Destroy AABB tree.
   class(aabb_tree_object), intent(inout) :: self  !< AABB tree.
   type(aabb_tree_object)                 :: fresh !< Fresh instance of AABB tree.

   self = fresh
   endsubroutine destroy

   ! pure subroutine initialize(self, refinement_levels, facet, bmin, bmax)
   subroutine initialize(self, refinement_levels, facet, bmin, bmax)
   !< Initialize AABB tree.
   class(aabb_tree_object), intent(inout)        :: self              !< AABB box.
   integer(I4P),            intent(in)           :: refinement_levels !< Total number of refinement levels used.
   type(facet_object),      intent(in), optional :: facet(:)          !< Facets list.
   type(vector_R8P),        intent(in), optional :: bmin              !< Minimum point of AABB.
   type(vector_R8P),        intent(in), optional :: bmax              !< Maximum point of AABB.
   integer(I4P)                                  :: level             !< Counter.
   integer(I4P)                                  :: b, bb, bbb, bbbb  !< Counter.
   type(aabb_node_object)                        :: parent            !< Parent node.
   type(aabb_object)                             :: octant(8)         !< AABB octants.

   call self%destroy
   self%refinement_levels = refinement_levels
   self%nodes_number = nodes_number(refinement_levels=self%refinement_levels)
   allocate(self%node(0:self%nodes_number-1))
   associate(node=>self%node)
      ! inizialize all tree nodes with only the bounding box
      call node(0)%initialize(facet=facet, bmin=bmin, bmax=bmax)
      do level=1, self%refinement_levels                                                             ! loop over refinement levels
         b = first_node(level=level)                                                                 ! first node at level
         do bb=1, nodes_number_at_level(level=level), TREE_RATIO                                     ! loop over nodes at level
            bbb = b + bb - 1                                                                         ! node numeration in tree
            parent = node(parent_node(node=bbb))                                                     ! parent of the current node
            if (allocated(parent%aabb)) then                                                         ! create children nodes
               call parent%aabb%compute_octant(octant)                                               ! compute parent AABB octants
               do bbbb=0, TREE_RATIO-1                                                               ! loop over children
                  call node(bbb+bbbb)%initialize(bmin=octant(bbbb+1)%bmin, bmax=octant(bbbb+1)%bmax) ! initialize node
               enddo
            endif
         enddo
      enddo
      ! fill all tree nodes with facets
   endassociate
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine aabb_tree_assign_aabb_tree(lhs, rhs)
   !< Operator `=`.
   class(aabb_tree_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_tree_object),  intent(in)    :: rhs !< Right hand side.
   integer                                :: b   !< Counter.

   if (allocated(lhs%node)) then
      do b=1, lhs%nodes_number
         call lhs%node%destroy
      enddo
      deallocate(lhs%node)
   endif
   lhs%refinement_levels = rhs%refinement_levels
   lhs%nodes_number = rhs%nodes_number
   if (allocated(rhs%node)) then
      allocate(lhs%node(0:lhs%nodes_number-1))
      do b=0, lhs%nodes_number-1
         lhs%node(b) = rhs%node(b)
      enddo
   endif
   endsubroutine aabb_tree_assign_aabb_tree

   ! non TBP
   pure function first_child_node(node)
   !< Return first child tree node.
   integer(I4P), intent(in) :: node             !< Node queried.
   integer(I4P)             :: first_child_node !< First child tree node.

   first_child_node = node * TREE_RATIO + 1
   endfunction first_child_node

   pure function first_node(level)
   !< Return first tree node at a given level.
   integer(I4P), intent(in) :: level      !< Refinement level queried.
   integer(I4P)             :: first_node !< Number of tree nodes at given level.

   first_node = nodes_number(refinement_levels=level-1)
   endfunction first_node

   pure function nodes_number(refinement_levels)
   !< Return total number of tree nodes given the total number refinement levels used.
   integer(I4P), intent(in) :: refinement_levels !< Total number of refinement levels used.
   integer(I4P)             :: nodes_number      !< Total number of tree nodes.
   integer                  :: level             !< Counter.

   nodes_number = 0
   do level=0, refinement_levels
      nodes_number = nodes_number + nodes_number_at_level(level=level)
   enddo
   endfunction nodes_number

   pure function nodes_number_at_level(level) result(nodes_number)
   !< Return number of tree nodes at a given level.
   integer(I4P), intent(in) :: level        !< Refinement level queried.
   integer(I4P)             :: nodes_number !< Number of tree nodes at given level.

   nodes_number = TREE_RATIO ** (level)
   endfunction nodes_number_at_level

   pure function parent_node(node)
   !< Return parent tree node.
   integer(I4P), intent(in) :: node        !< Node queried.
   integer(I4P)             :: parent_node !< Parent tree node.

   parent_node = (node - 1) / TREE_RATIO
   endfunction parent_node
endmodule fossil_aabb_tree_object
