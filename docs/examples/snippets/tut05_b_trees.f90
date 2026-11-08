!< Tutorial 5, second run: the same batch on the two kinds of AABB tree. The figures draw the boxes of the BVH four
!< levels below the root and of the octree two levels below (its leaves at the default depth), built over the facets
!< of the placed dragon by fossil_aabb_tree_object as analyze builds them (the surface's own tree is private): the BVH
!< splits the facets, the octree splits the space.
program tut05_b_trees
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_aabb_node_object, only : aabb_node_object, bvh_bmin_x_oac, bvh_bmin_y_oac, bvh_bmin_z_oac, &
                                    bvh_bmax_x_oac, bvh_bmax_y_oac, bvh_bmax_z_oac
use fossil_facet_object,     only : facet_object
use fossil_fig,              only : fig_open, fig_surface, fig_box, fig_close
use fossil,                  only : surface_stl_object, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16                      ! points per direction
type(surface_stl_object)      :: bvh, octree
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_bvh(:), d_oct(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p

allocate(points(n**3), d_bvh(n**3), d_oct(n**3))
h = 1.2_R8P / (n - 1)                                        ! the cube [-0.6, 0.6]^3
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call bvh%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)       ! default AABB_TREE_SAH_BVH
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call bvh%distance_many(points=points, distances=d_bvh, is_signed=.true., is_square_root=.true.)
call octree%distance_many(points=points, distances=d_oct, is_signed=.true., is_square_root=.true.)
print '(A,I0)',     'points                ', size(points)
print '(A,ES10.3)', 'max ||bvh| - |oct||  ', maxval(abs(abs(d_bvh) - abs(d_oct)))
print '(A,I0)',     'bvh, oct signs differ ', count(d_bvh * d_oct < 0._R8P)
call draw('tut05-bvh', AABB_TREE_SAH_BVH, 'SAH BVH', 4_I4P)
call draw('tut05-octree', AABB_TREE_OCTREE, 'octree', 2_I4P)
contains
   subroutine draw(name, kind, title, depth)
   !< The boxes of a tree of the given kind `depth` levels below the root, over the dragon.
   character(*), intent(in)        :: name, title !< Figure name, tree name.
   integer(I4P), intent(in)        :: kind        !< Tree kind.
   integer(I4P), intent(in)        :: depth       !< Levels below the root.
   type(aabb_tree_object)          :: tree
   type(facet_object), allocatable :: facet(:)
   type(aabb_node_object), pointer :: node
   integer(I4P), allocatable       :: level(:), next(:)
   integer(I4P)                    :: child(8), nchild, f, c, l, u
   real(R8P)                       :: longest
   character(96)                   :: caption

   allocate(facet(bvh%get_facets_number()))
   longest = 0._R8P
   do f = 1, size(facet)
      facet(f) = bvh%facet_at(f)
      longest = max(longest, edge(facet(f)%vertex(1), facet(f)%vertex(2)), edge(facet(f)%vertex(2), facet(f)%vertex(3)), &
                    edge(facet(f)%vertex(3), facet(f)%vertex(1)))
   enddo
   call tree%set_tree_kind(kind)
   call tree%initialize(facet=facet, largest_edge_len=longest)
   level = [0_I4P]                                             ! the root, then its children, depth times
   do l = 1, depth
      allocate(next(0))
      do f = 1, size(level)
         call tree%enumerate_children(level(f), child, nchild)
         next = [next, child(1:nchild)]
      enddo
      call move_alloc(next, level)
   enddo
   write(caption, '(A,A,I0,A,I0,A,I0,A)') title, ': ', size(level), ' boxes ', depth, ' levels below the root, ', &
                                          tree%get_nodes_number(), ' nodes in all'
   print '(A)', trim(caption)
   call fig_open(u, name, trim(caption), -60._R8P, 22._R8P)
   call fig_surface(u, 'body', bvh)
   do c = 1, size(level)
      node => tree%node_at(level(c))
      call fig_box(u, 'accent', bvh_bmin_x_oac(node) * ex_R8P + bvh_bmin_y_oac(node) * ey_R8P + bvh_bmin_z_oac(node) * ez_R8P, &
                   bvh_bmax_x_oac(node) * ex_R8P + bvh_bmax_y_oac(node) * ey_R8P + bvh_bmax_z_oac(node) * ez_R8P)
   enddo
   call fig_close(u)
   endsubroutine draw

   pure function edge(a, b) result(l)
   !< Length of the edge a-b.
   type(vector_R8P), intent(in) :: a, b !< Ends.
   real(R8P)                    :: l    !< Length.

   l = norm2([b%x - a%x, b%y - a%y, b%z - a%z])
   endfunction edge
endprogram tut05_b_trees
