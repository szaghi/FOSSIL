!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

module fossil_aabb_tree_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

use fossil_aabb_object, only : aabb_object
use fossil_aabb_node_object, only : aabb_node_object
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use penf, only : I4P, R8P, MaxR8P, str
use vecfor, only : vector_R8P

implicit none
private
public :: aabb_tree_object

integer(I4P), parameter :: TREE_RATIO=8 !< Tree refinement ratio, it is assumed to be an **octree**.

type :: ofsm
   !< Octree Finite State Machine class for efficient searching of neiighbors.
   integer(I4P) :: octant=0    !< Octant ID, [1,2,3,4,5,6,7,8].
   integer(I4P) :: direction=0 !< Direction, 0=Halt, 1=X+, 2=X-, 3=Y+, 4=Y-, 5=Z+, 6=Z-.
endtype ofsm

type(ofsm) :: octree_fsm(6,8) = reshape(                                                                   &
                                        [ofsm(2,0), ofsm(2,2), ofsm(3,0), ofsm(3,4), ofsm(5,0), ofsm(5,6), &
                                         ofsm(1,1), ofsm(1,0), ofsm(4,0), ofsm(4,4), ofsm(6,0), ofsm(6,6), &
                                         ofsm(4,0), ofsm(4,2), ofsm(1,3), ofsm(1,0), ofsm(7,0), ofsm(7,6), &
                                         ofsm(3,1), ofsm(3,0), ofsm(2,3), ofsm(2,0), ofsm(8,0), ofsm(8,6), &
                                         ofsm(6,0), ofsm(6,2), ofsm(7,0), ofsm(7,4), ofsm(1,5), ofsm(1,0), &
                                         ofsm(5,1), ofsm(5,0), ofsm(8,0), ofsm(8,4), ofsm(2,5), ofsm(2,0), &
                                         ofsm(8,0), ofsm(8,2), ofsm(5,3), ofsm(5,0), ofsm(3,5), ofsm(3,0), &
                                         ofsm(7,1), ofsm(7,0), ofsm(6,3), ofsm(6,0), ofsm(4,5), ofsm(4,0)], [6,8])

type :: aabb_tree_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree class.
   !<
   !< @note The tree is assumed to be an **octree**.
   !< The octree uses a breath-leaf counting, with the following convention:
   !<
   !<```
   !>                  +----+----+
   !>                 /|   /|   /|
   !>                / |  7 |  8 |
   !>               /  +----*----+
   !>              /  /|/  /|/  /|
   !>             /  / |  5 |  6 |
   !>            /  / /+----+----+
   !>           /  / // / // /  /
   !>          +----+----+/ /  /
   !>         /| / /| / /| /  /
   !>        / |/ //|/ //|/  /
   !>       /  +----*----+  /
   !>      /  /|// /|// /| /
   !>     /  / |/ / |/ / |/
   !>    /  / /+----+----+
   !>   /  / // / // /  /
   !>  +----+----+/ /  /       y(j)  z(k)
   !>  | /  | /  | /  /        ^    ^
   !>  |/ 3/|/ 4/|/  /         |   /
   !>  +----*----+  /          |  /
   !>  | /  | /  | /           | /
   !>  |/ 1 |/ 2 |/            |/
   !>  +----+----+             +------->x(i)
   !<```
   integer(I4P)                        :: refinement_levels=2    !< Total number of refinement levels used.
   integer(I4P)                        :: nodes_number=0         !< Total number of tree nodes.
   type(aabb_node_object), allocatable :: node(:)                !< AABB tree nodes [0:nodes_number-1].
   logical                             :: is_initialized=.false. !< Sentinel to check is AABB tree is initialized.
   contains
      ! public methods
      procedure, pass(self) :: compute_vertices_nearby     !< Compute vertices nearby.
      procedure, pass(self) :: destroy                     !< Destroy AABB tree.
      procedure, pass(self) :: distance                    !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: distance_tree               !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: distribute_facets           !< Distribute facets into AABB nodes.
      procedure, pass(self) :: distribute_facets_tree      !< Distribute facets into AABB nodes.
      procedure, pass(self) :: has_children                !< Return true if node has at least one child allocated.
      procedure, pass(self) :: initialize                  !< Initialize AABB tree.
      procedure, pass(self) :: loop_node                   !< Loop over all nodes.
      procedure, pass(self) :: ray_intersections_number    !< Return ray intersections number.
      procedure, pass(self) :: save_geometry_tecplot_ascii !< Save AABB tree boxes geometry into Tecplot ascii file.
      procedure, pass(self) :: translate                   !< Translate AABB tree by delta.
      ! operators
      generic :: assignment(=) => aabb_tree_assign_aabb_tree      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_tree_assign_aabb_tree !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: distance_node                 !< Return the (minimum) distance from point to node AABB tree.
      procedure, pass(self), private :: ray_intersections_number_node !< Return ray intersections number into a node of AABB tree.
endtype aabb_tree_object

contains
   ! public methods
   pure subroutine compute_vertices_nearby(self, facet, tolerance_to_be_identical, tolerance_to_be_nearby)
   !< Compute vertices nearby.
   class(aabb_tree_object), intent(in)    :: self                      !< AABB tree.
   type(facet_object),      intent(inout) :: facet(:)                  !< Facets list.
   real(R8P),               intent(in)    :: tolerance_to_be_identical !< Tolerance to identify identical vertices.
   real(R8P),               intent(in)    :: tolerance_to_be_nearby    !< Tolerance to identify nearby vertices.
   integer(I4P)                           :: level                     !< Counter.
   integer(I4P)                           :: b, bb, bbb                !< Counter.

   if (self%nodes_number > 0) then
      level=self%refinement_levels                                      ! check only max refinement level
      b = first_node(level=level)                                       ! first node at level
      do bb=1, nodes_number_at_level(level=level)                       ! loop over nodes at level
         bbb = b + bb - 1                                               ! node numeration in tree
         call self%node(bbb)%compute_vertices_nearby(facet=facet,                                         &
                                                     tolerance_to_be_identical=tolerance_to_be_identical, &
                                                     tolerance_to_be_nearby=tolerance_to_be_nearby)
      enddo
   endif
   endsubroutine compute_vertices_nearby

   elemental subroutine destroy(self)
   !< Destroy AABB tree.
   class(aabb_tree_object), intent(inout) :: self  !< AABB tree.
   type(aabb_tree_object)                 :: fresh !< Fresh instance of AABB tree.

   self = fresh
   endsubroutine destroy

   pure function distance(self, facet, point)
   !< Compute the (minimum) distance from a point to the triangulated surface.
   class(aabb_tree_object), intent(in) :: self            !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)        !< Facets list.
   type(vector_R8P),        intent(in) :: point           !< Point coordinates.
   real(R8P)                           :: distance        !< Minimum distance from point to the triangulated surface.
   real(R8P), allocatable              :: distance_(:)    !< Minimum distance, temporary buffer.
   integer(I4P), allocatable           :: aabb_closest(:) !< Index of closest AABB.
   integer(I4P)                        :: level           !< Counter.
   integer(I4P)                        :: b, bb, bbb      !< Counter.

   associate(node=>self%node)
      allocate(distance_(0:self%refinement_levels))
      allocate(aabb_closest(0:self%refinement_levels))
      distance_ = MaxR8P
      aabb_closest = -1
      do level=0, self%refinement_levels                  ! loop over refinement levels
         b = first_node(level=level)                      ! first node at finest level
         do bb=1, nodes_number_at_level(level=level)      ! loop over nodes at level
            bbb = b + bb - 1                              ! node numeration in tree
            if (node(bbb)%is_allocated()) then
               distance = node(bbb)%distance(point=point) ! node distance
               if (distance <= distance_(level)) then
                  distance_(level) = distance             ! update minimum distance
                  aabb_closest(level) = bbb               ! store closest node
               endif
            endif
         enddo
      enddo
      distance = MaxR8P
      do level=0, self%refinement_levels
         if (aabb_closest(level) >= 0) then
            distance = min(distance, node(aabb_closest(level))%distance_from_facets(facet=facet, point=point))
         endif
      enddo
   endassociate
   endfunction distance

   function distance_tree(self, facet, point) result(distance)
   !< Compute the (minimum) distance from a point to the triangulated surface.
   class(aabb_tree_object), intent(in) :: self         !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)     !< Facets list.
   type(vector_R8P),        intent(in) :: point        !< Point coordinates.
   real(R8P)                           :: distance     !< Minimum distance from point to the triangulated surface.

   distance = self%distance_node(n=0, facet=facet, point=point)
   endfunction distance_tree

   pure subroutine distribute_facets(self, facet, is_exclusive, do_update_extents)
   !< Distribute facets into AABB nodes.
   class(aabb_tree_object), intent(inout)        :: self               !< AABB tree.
   type(facet_object),      intent(in)           :: facet(:)           !< Facets list.
   logical,                 intent(in), optional :: is_exclusive       !< Sentinel to enable/disable exclusive addition.
   logical,                 intent(in), optional :: do_update_extents  !< Sentinel to enable/disable AABB extents update.
   logical                                       :: do_update_extents_ !< Sentinel to enable/disable AABB extents update, local var.
   type(list_id_object)                          :: facet_id           !< List of facets IDs.
   integer(I4P)                                  :: level              !< Counter.
   integer(I4P)                                  :: b, bb, bbb         !< Counter.

   do_update_extents_ = .true. ; if (present(do_update_extents)) do_update_extents_ = do_update_extents
   associate(node=>self%node)
      call facet_id%initialize(id=facet%id)

      ! add facets to nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (facet_id%ids_number > 0) then                       ! check if facets list still has facets
               call node(bbb)%add_facets(facet_id=facet_id, &
                                         facet=facet,       &
                                         is_exclusive=is_exclusive) ! add facets to node and prune added facets from list
            endif
         enddo
      enddo

      ! destroy void nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (.not.node(bbb)%has_facets()) call node(bbb)%destroy ! destroy void node
         enddo
      enddo

      ! update AABB extents
      if (do_update_extents_) then
         do level=self%refinement_levels, 0, -1           ! loop over refinement levels
            b = first_node(level=level)                   ! first node at level
            do bb=1, nodes_number_at_level(level=level)   ! loop over nodes at level
               bbb = b + bb - 1                           ! node numeration in tree
               call node(bbb)%update_extents(facet=facet) ! update extents
            enddo
         enddo
      endif
   endassociate
   endsubroutine distribute_facets

   pure subroutine distribute_facets_tree(self, facet)
   !< Distribute facets into AABB nodes.
   class(aabb_tree_object), intent(inout) :: self       !< AABB tree.
   type(facet_object),      intent(in)    :: facet(:)   !< Facets list.
   type(list_id_object)                   :: facet_id   !< List of facets IDs.
   integer(I4P)                           :: level      !< Counter.
   integer(I4P)                           :: b, bb, bbb !< Counter.
   integer(I4P)                           :: parent     !< Parent node index.

   associate(node=>self%node)
      ! add facets to nodes
      call facet_id%initialize(id=facet%id)                             ! initialize facets IDs list
      if (facet_id%ids_number > 0) then                                 ! check if facets list still has facets
         call node(0)%add_facets(facet_id=facet_id, &
                                 facet=facet,       &
                                 is_exclusive=.false.)                  ! add facets to root node
         do level=1, self%refinement_levels                             ! loop over refinement levels
            b = first_node(level=level)                                 ! first node at level
            do bb=1, nodes_number_at_level(level=level), TREE_RATIO     ! loop over nodes at level
               parent = parent_node(node=b + bb - 1)                    ! parent of the current node
               if (node(parent)%is_allocated()) then                    ! check if parent exist
                  facet_id = node(parent)%facet_id()                    ! store parent facets IDs list
                  if (facet_id%ids_number > 0) then                     ! check if facets list still has facets
                     do bbb=b + bb - 1, b + bb -1 + TREE_RATIO - 1
                        call node(bbb)%add_facets(facet_id=facet_id, &
                                                  facet=facet,       &
                                                  is_exclusive=.true.) ! add facets to node
                     enddo
                  endif
               endif
            enddo
         enddo
      endif

      ! destroy void nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (.not.node(bbb)%has_facets()) call node(bbb)%destroy ! destroy void node
         enddo
      enddo

      ! update AABB extents
      ! do level=self%refinement_levels, 0, -1           ! loop over refinement levels
      !    b = first_node(level=level)                   ! first node at level
      !    do bb=1, nodes_number_at_level(level=level)   ! loop over nodes at level
      !       bbb = b + bb - 1                           ! node numeration in tree
      !       call node(bbb)%update_extents(facet=facet) ! update extents
      !    enddo
      ! enddo
   endassociate
   endsubroutine distribute_facets_tree

   pure function has_children(self, node)
   !< Return true if node has at least one child allocated.
   class(aabb_tree_object), intent(in) :: self         !< AABB tree.
   integer(I4P),            intent(in) :: node         !< Node queried.
   logical                             :: has_children !< Check result.
   integer                             :: n, nn        !< Counter.

   has_children = .false.
   n = first_child_node(node=node)
   if (n<=self%nodes_number - TREE_RATIO + 1) then
      do nn=n, n + TREE_RATIO - 1
         if (self%node(nn)%is_allocated()) then
            has_children = .true.
            return
         endif
      enddo
   endif
   endfunction has_children

   pure subroutine initialize(self, refinement_levels, facet, largest_edge_len, bmin, bmax, do_facets_distribute, is_exclusive, &
                              do_update_extents)
   !< Initialize AABB tree.
   class(aabb_tree_object), intent(inout)        :: self                  !< AABB tree.
   integer(I4P),            intent(in), optional :: refinement_levels     !< AABB refinement levels.
   type(facet_object),      intent(in), optional :: facet(:)              !< Facets list.
   real(R8P),               intent(in), optional :: largest_edge_len      !< Largest edge lenght.
   type(vector_R8P),        intent(in), optional :: bmin                  !< Minimum point of AABB.
   type(vector_R8P),        intent(in), optional :: bmax                  !< Maximum point of AABB.
   logical,                 intent(in), optional :: do_facets_distribute  !< Sentinel to enable/disable facets distribution.
   logical,                 intent(in), optional :: is_exclusive          !< Sentinel to enable/disable exclusive addition.
   logical,                 intent(in), optional :: do_update_extents     !< Sentinel to enable/disable AABB extents update.
   integer(I4P)                                  :: refinement_levels_    !< AABB refinement levels, local variable.
   logical                                       :: do_facets_distribute_ !< Sentinel to enable/dis. facets distribution, local var.
   integer(I4P)                                  :: level                 !< Counter.
   integer(I4P)                                  :: b, bb, bbb, bbbb      !< Counter.
   integer(I4P)                                  :: parent                !< Parent node index.
   type(aabb_object)                             :: octant(8)             !< AABB octants.

   refinement_levels_ = self%refinement_levels
   call self%destroy
   self%refinement_levels = refinement_levels_ ; if (present(refinement_levels)) self%refinement_levels = refinement_levels
   do_facets_distribute_ = .true. ; if (present(do_facets_distribute)) do_facets_distribute_ = do_facets_distribute

   if (self%refinement_levels >= 0) then
      self%nodes_number = nodes_number(refinement_levels=self%refinement_levels)
      allocate(self%node(0:self%nodes_number-1))
      call self%node(0)%initialize(facet=facet, bmin=bmin, bmax=bmax)
      levels_loop: do level=1, self%refinement_levels                             ! loop over refinement levels
         b = first_node(level=level)                                              ! first node at level
         do bb=1, nodes_number_at_level(level=level), TREE_RATIO                  ! loop over nodes at level
            bbb = b + bb - 1                                                      ! node numeration in tree
            parent = parent_node(node=bbb)                                        ! parent of the current node
            if (self%node(parent)%is_allocated()) then                            ! create children nodes
               call self%node(parent)%compute_octants(octant=octant)              ! compute parent AABB octants
               if (present(largest_edge_len)) then
                  if (largest_edge_len > octant(1)%median()) then                 ! check if refinement has sense
                     ! a further refinement does not have sense
                     self%refinement_levels = level - 1                           ! set rifinement to the previous one
                     exit levels_loop                                             ! exi loop
                  endif
               endif
               do bbbb=0, TREE_RATIO-1                                            ! loop over children
                  call self%node(bbb+bbbb)%initialize(bmin=octant(bbbb+1)%bmin, &
                                                      bmax=octant(bbbb+1)%bmax)   ! initialize node
               enddo
            endif
         enddo
      enddo levels_loop
      if (present(facet).and.(do_facets_distribute_)) call self%distribute_facets_tree(facet=facet)
      self%is_initialized = .true.
   endif
   endsubroutine initialize

   function loop_node(self, facet, aabb_facet, b, l) result(again)
   !< Loop over all nodes.
   !<
   !< @note Impure function: return data of each allocated node exploiting saved local counter.
   class(aabb_tree_object),         intent(in)            :: self          !< AABB tree.
   type(facet_object),              intent(in),  optional :: facet(:)      !< Whole facets list.
   integer(I4P),                    intent(out), optional :: b             !< Current AABB ID.
   integer(I4P),                    intent(out), optional :: l             !< Current AABB level.
   type(facet_object), allocatable, intent(out), optional :: aabb_facet(:) !< AABB facets list.
   logical                                                :: again         !< Flag continuing the loop.
   integer(I4P), save                                     :: bb = -1       !< AABB ID counter.
   integer(I4P)                                           :: bbb           !< Counter.

   again = .false.
   if (allocated(self%node)) then
      if (bb==-1) then
         ! get first allocated node
         do bbb=0, self%nodes_number - 1
            if (self%node(bbb)%is_allocated()) then
               again = .true.
               if (present(facet).and.present(aabb_facet)) call self%node(bbb)%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
               exit
            endif
         enddo
         bb = bbb
      elseif (bb<self%nodes_number - 1) then
         do bbb=bb+1, self%nodes_number - 1
            if (self%node(bbb)%is_allocated()) then
               again = .true.
               if (present(facet).and.present(aabb_facet)) call self%node(bbb)%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
               exit
            endif
         enddo
         bb = bbb
      else
         bb = -1
         again = .false.
      endif
   endif
   if (present(b)) b = bb
   if (present(l)) l = level(b)
   endfunction loop_node

   function ray_intersections_number(self, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number.
   class(aabb_tree_object), intent(in) :: self                 !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)             !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                        :: intersections_number !< Intersection number.
   ! integer(I4P)                        :: level                !< Counter.
   ! integer(I4P)                        :: b, bb, bbb           !< Counter.

   intersections_number = self%ray_intersections_number_node(n=0, facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
   ! intersections_number = 0
   ! associate(node=>self%node)
   !    do level=0, self%refinement_levels                  ! loop over refinement levels
   !       b = first_node(level=level)                      ! first node at finest level
   !       do bb=1, nodes_number_at_level(level=level)      ! loop over nodes at level
   !          bbb = b + bb - 1                              ! node numeration in tree
   !          if (node(bbb)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) then
   !             intersections_number = intersections_number + &
   !                                    node(bbb)%ray_intersections_number(facet=facet, &
   !                                                                       ray_origin=ray_origin, ray_direction=ray_direction)
   !          endif
   !       enddo
   !    enddo
   ! endassociate
   endfunction ray_intersections_number

   subroutine save_geometry_tecplot_ascii(self, file_name)
   !< Save AABB tree boxes geometry into Tecplot ascii file.
   class(aabb_tree_object), intent(in) :: self       !< AABB tree.
   character(*),            intent(in) :: file_name  !< File name.
   integer(I4P)                        :: file_unit  !< File unit.
   integer(I4P)                        :: level      !< Counter.
   integer(I4P)                        :: b, bb, bbb !< Counter.

   associate(node=>self%node)
      if (self%is_initialized) then
         open(newunit=file_unit, file=trim(adjustl(file_name)))
         write(file_unit, '(A)') 'VARIABLES=x y z'
         do level=0, self%refinement_levels
            b = first_node(level=level)
            do bb=1, nodes_number_at_level(level=level)
               bbb = b + bb - 1
               call node(bbb)%save_geometry_tecplot_ascii(file_unit=file_unit, aabb_name='aabb-l_'//trim(str(level, .true.))//&
                                                                                             '-b_'//trim(str(bbb, .true.)))
            enddo
         enddo
         close(file_unit)
      endif
   endassociate
   endsubroutine save_geometry_tecplot_ascii

   elemental subroutine translate(self, delta)
   !< Translate AABB tree by delta.
   class(aabb_tree_object), intent(inout) :: self  !< AABB.
   type(vector_R8P),        intent(in)    :: delta !< Delta of translation.
   integer(I4P)                           :: n     !< Counter.

   if (self%nodes_number > 0) then
      do n=0, self%nodes_number - 1
         call self%node(n)%translate(delta=delta)
      enddo
   endif
   endsubroutine translate

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
   lhs%is_initialized = rhs%is_initialized
   endsubroutine aabb_tree_assign_aabb_tree

   ! private methods
   recursive function distance_node(self, n, facet, point) result(distance)
   !< Return the (minimum) distance from a point to a node of AABB tree.
   class(aabb_tree_object), intent(in) :: self         !< AABB tree.
   integer(I4P),            intent(in) :: n            !< Current AABB node.
   type(facet_object),      intent(in) :: facet(:)     !< Facets list.
   type(vector_R8P),        intent(in) :: point        !< Point coordinates.
   real(R8P)                           :: distance     !< Minimum distance from point to the triangulated surface.
   real(R8P)                           :: distance_    !< Minimum distance, temporary buffer.
   integer(I4P)                        :: aabb_closest !< Closest AABB children node.
   integer(I4P)                        :: fcn          !< First AABB child node.
   integer(I4P)                        :: i            !< Counter.

   associate(node=>self%node)
      if (self%has_children(node=n)) then
         distance_ = MaxR8P                             ! initialize distance of current level-nodes
         aabb_closest = -1                              ! initialize closest node index
         fcn = first_child_node(node=n)                 ! first child node
         do i=fcn, fcn + TREE_RATIO - 1                 ! loop over all children nodes
            if (node(i)%is_allocated()) then            ! check if node is allocated
               distance = node(i)%distance(point=point) ! node distance
               if (distance <= distance_) then          ! check for the new minimum
                  distance_ = distance                  ! update minimum distance
                  aabb_closest = i                      ! store closest node
               endif
            endif
         enddo
         distance = self%distance_node(n=aabb_closest, facet=facet, point=point) ! return distance from closest AABB child node
      else
         distance = node(n)%distance_from_facets(facet=facet, point=point)       ! no children: return distance from current node
      endif
   endassociate
   endfunction distance_node

   recursive function ray_intersections_number_node(self, n, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number into a node of AABB tree.
   class(aabb_tree_object), intent(in) :: self                 !< AABB tree.
   integer(I4P),            intent(in) :: n                    !< Current AABB node.
   type(facet_object),      intent(in) :: facet(:)             !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                        :: intersections_number !< Intersection number.
   integer(I4P)                        :: fcn                  !< First AABB child node.
   integer(I4P)                        :: i                    !< Counter.

   intersections_number = 0
   associate(node=>self%node)
      if (node(n)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) then      ! check if ray intersect AABB
         if (self%has_children(node=n)) then                                                      ! check if AABB has children
            fcn = first_child_node(node=n)                                                        ! first child node
            do i=fcn, fcn + TREE_RATIO - 1                                                        ! loop over all children nodes
               intersections_number = intersections_number +                                    & ! sum chidren intersections
                                      self%ray_intersections_number_node(n=i,                   &
                                                                         facet=facet,           &
                                                                         ray_origin=ray_origin, &
                                                                         ray_direction=ray_direction)
            enddo
         else
            ! there are not children, return intersection of current AABB leaf
            intersections_number = node(n)%ray_intersections_number(facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
         endif
      endif
   endassociate
   endfunction ray_intersections_number_node

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
   integer(I4P)             :: first_node !< First tree node at given level.

   first_node = nodes_number(refinement_levels=level-1)
   endfunction first_node

   pure function last_node(level)
   !< Return last tree node at a given level.
   integer(I4P), intent(in) :: level     !< Refinement level queried.
   integer(I4P)             :: last_node !< Last tree node at given level.

   last_node = first_node(level) + nodes_number_at_level(level) - 1
   endfunction last_node

   pure function level(node)
   !< Return level given a node id.
   integer(I4P), intent(in) :: node  !< Node queried.
   integer(I4P)             :: level !< Level of given node.
   integer(I4P)             :: n     !< Counter.

   level = 0
   n = node
   do while (n /= 0)
      n = (n - 1) / TREE_RATIO
      level = level + 1
   enddo
   endfunction level

   pure function local_id(node)
   !< Return local ID of node.
   integer(I4P), intent(in) :: node     !< Node queried.
   integer(I4P)             :: local_id !< Local ID.

   local_id = node - ((node - 1) / TREE_RATIO) * TREE_RATIO
   endfunction local_id

   pure function location_code(node)
   !< Return location code of node.
   integer(I4P), intent(in)  :: node             !< Node queried.
   integer(I4P), allocatable :: location_code(:) !< Location code.
   integer(I4P)              :: node_level       !< Node level.
   integer(I4P)              :: parent           !< Parent node.
   integer(I4P)              :: l                !< Counter.

   node_level = level(node=node)
   if (node_level > 0) then
      allocate(location_code(1:node_level))
      parent = node
      do l=node_level, 1, -1
         location_code(l) = local_id(parent)
         parent = parent_node(node=node)
      enddo
   else
      allocate(location_code(1))
      location_code(1) = 0
   endif
   endfunction location_code

   pure function global_id(location_code) result(node)
   !< Return node ID given a location code.
   integer(I4P), intent(in)  :: location_code(:) !< Location code.
   integer(I4P)              :: node             !< Node global ID.
   integer(I4P)              :: l                !< Counter.

   node = 0
   if (location_code(1) == 0) then
      return
   else
      node = location_code(1)
      do l=2, size(location_code, dim=1)
         node = first_child_node(node)
         node = node + location_code(l) - 1
      enddo
   endif
   endfunction global_id

   pure subroutine next_location_code(location_code, direction, next_code, next_direction)
   !< Return the node next along a given direction.
   integer(I4P),              intent(in)  :: location_code(:) !< Location code queried.
   integer(I4P),              intent(in)  :: direction        !< Direction, 0=Halt, 1=X+, 2=X-, 3=Y+, 4=Y-, 5=Z+, 6=Z-.
   integer(I4P), allocatable, intent(out) :: next_code(:)     !< Next location code of node along given direction.
   integer(I4P),              intent(out) :: next_direction   !< Next direction.
   integer(I4P)                           :: direction_       !< Direction, local variable.
   integer(I4P)                           :: o                !< Counter.

   direction_ = direction
   next_code = location_code
   do o=size(location_code, dim=1), 1, -1
      if (direction_ > 0.and.location_code(o) > 0) then
         next_code(o) = octree_fsm(direction_, location_code(o))%octant
         direction_ = octree_fsm(direction_, location_code(o))%direction
      else
         exit
      endif
   enddo
   next_direction = direction_
   endsubroutine next_location_code

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

   pure function parent_at_level(node, parent_level) result(parent)
   !< Return parent tree node at a given level.
   integer(I4P), intent(in) :: node         !< Node.
   integer(I4P), intent(in) :: parent_level !< Level.
   integer(I4P)             :: parent       !< Parent.
   integer(I4P)             :: n            !< Counter.

   parent = node
   do n=1, level(node) - parent_level
     parent = (parent - 1) / TREE_RATIO
   enddo
   endfunction parent_at_level

   pure function parent_node(node)
   !< Return parent tree node.
   integer(I4P), intent(in) :: node        !< Node queried.
   integer(I4P)             :: parent_node !< Parent tree node.

   parent_node = (node - 1) / TREE_RATIO
   endfunction parent_node

   pure function siblings(node) result(sbs)
   !< Return siblings of a given node.
   integer(I4P), intent(in) :: node     !< Node queried.
   integer(I4P)             :: sbs(1:7) !< Nodes sibling.
   integer(I4P)             :: myid     !< Local node ID into siblings list.
   integer(I4P)             :: i, s     !< Counter.

   myid = local_id(node=node)
   s = 0
   do i=1, TREE_RATIO
      if (i /= myid) then
         s = s + 1
         sbs(s) = ((node - 1) / TREE_RATIO) * TREE_RATIO + i
      endif
   enddo
   endfunction siblings

   pure function str_location_code(code)
   !< Return string of location code of node.
   integer(I4P), intent(in)      :: code(:)           !< Node location code.
   character(len=:), allocatable :: str_location_code !< String of Location code.
   integer(I4P)                  :: c                 !< Counter.

   str_location_code = ''
   do c=1, size(code, dim=1)
      str_location_code = str_location_code//trim(str(no_sign=.true., n=code(c)))
   enddo
   endfunction str_location_code
endmodule fossil_aabb_tree_object
