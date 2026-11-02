!< FOSSIL, assert the SAH BVH depth cap keeps the traversal stack in bounds.
!<
!< The distance traversal (CPU and device) uses a fixed explicit stack of
!< MAX_TRAVERSAL_STACK = 128 entries with an unchecked push; a binary BVH of depth
!< D needs D + 1 entries. This mesh of 1001 self-similar triangles (triangle k at
!< x = 2^k, size 0.25 * 2^k, k in [-500, 500]) makes the bucketed SAH peel a few
!< triangles per level: uncapped it builds a tree of depth 208, and
!< `surface%distance` overflowed the stack and crashed. With the cap the tree
!< must stop at depth MAX_BVH_DEPTH = 127, and the BVH distance must stay
!< bit-exact against a brute-force scan over every facet.
!<
!< The extreme dynamic range is unavoidable: the 16-bucket SAH peels about 4-5
!< binary octaves of scale per level, so depth > 127 needs ~600 octaves. At the
!< ends of the range the facet metrics leave the normal floating-point range
!< (`det = a*c - b^2` scales as size^4: it underflows below k ~ -268 and
!< overflows above k ~ 258), where the facet kernel's result depends on FMA
!< contraction and is not meaningful. The exactness check therefore queries
!< only scales 2^j with |j| <= 250, whose nearby facets have normal metrics;
!< the depth check covers the whole mesh.

program fossil_test_bvh_depth_cap

use fossil,                  only : surface_stl_object, facet_object, AABB_TREE_SAH_BVH
use fossil_aabb_node_object, only : aabb_node_object
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

integer(I4P), parameter         :: K_MIN = -500_I4P    !< Smallest exponent of the triangle positions.
integer(I4P), parameter         :: K_MAX =  500_I4P    !< Largest exponent of the triangle positions.
integer(I4P), parameter         :: EXPECTED_DEPTH = 127_I4P !< MAX_BVH_DEPTH = MAX_TRAVERSAL_STACK - 1.
type(surface_stl_object)        :: surface              !< Pathological surface.
type(facet_object), allocatable :: facets(:)            !< Facets handed to the surface.
type(facet_object), pointer     :: f                    !< Facet iterator.
type(vector_R8P)                :: p                    !< Query point.
real(R8P)                       :: x, sz                !< Triangle position and size.
real(R8P)                       :: d_tree, d_brute, d   !< Squared distances.
integer(I4P)                    :: k, i, j, depth       !< Counters, tree depth.
integer(I4P)                    :: n_queries, n_mismatch !< Query statistics.
logical                         :: depth_ok, exact_ok   !< Checks.

allocate(facets(K_MAX - K_MIN + 1))
do k = K_MIN, K_MAX
   x  = 2._R8P ** k
   sz = 0.25_R8P * x
   i  = k - K_MIN + 1
   facets(i)%vertex(1) = vector_R8P(x,      0._R8P, 0._R8P)
   facets(i)%vertex(2) = vector_R8P(x + sz, 0._R8P, 0._R8P)
   facets(i)%vertex(3) = vector_R8P(x,      sz,     0._R8P)
enddo
call surface%adopt_facets(facets=facets, aabb_tree_kind=AABB_TREE_SAH_BVH)

depth = tree_depth(node_index=0_I4P, level=0_I4P)
depth_ok = (depth == EXPECTED_DEPTH)

! Query points at scales 2^j, |j| <= 250 (facet metrics in the normal range).
n_queries  = 0_I4P
n_mismatch = 0_I4P
do j = -250_I4P, 250_I4P, 10_I4P
   p = vector_R8P(1.3_R8P * 2._R8P ** j, 0.1_R8P * 2._R8P ** j, 0.2_R8P * 2._R8P ** j)
   d_tree  = surface%distance(point=p)
   d_brute = huge(1._R8P)
   do i = 1_I4P, surface%get_facets_number()
      f => surface%facet_at(i)
      call f%compute_distance(point=p, distance=d)
      d_brute = min(d_brute, d)
   enddo
   n_queries = n_queries + 1_I4P
   if (d_tree /= d_brute) n_mismatch = n_mismatch + 1_I4P
enddo
exact_ok = (n_mismatch == 0_I4P)

print '(A,I0,A,I0,A,I0,A,I0,A)', 'facets=', surface%get_facets_number(), ' nodes=', surface%aabb%get_nodes_number(), &
                                 ' depth=', depth, ' (expected ', EXPECTED_DEPTH, ')'
print '(A,I0,A,I0)', 'BVH vs brute-force d^2 mismatches: ', n_mismatch, ' of ', n_queries
print '(A,L1)', 'Are all tests passed? ', (depth_ok .and. exact_ok)
if (.not. (depth_ok .and. exact_ok)) error stop 1

contains
   recursive function tree_depth(node_index, level) result(depth_max)
   !< Return the depth of the deepest leaf below `node_index`, which sits at `level`.
   integer(I4P), intent(in)        :: node_index !< Node index.
   integer(I4P), intent(in)        :: level      !< Depth of the node.
   integer(I4P)                    :: depth_max  !< Deepest leaf depth.
   type(aabb_node_object), pointer :: node       !< Node.
   integer(I4P)                    :: left, right !< Children.

   node  => surface%aabb%node_at(node_index)
   left  = node%get_left_child()
   right = node%get_right_child()
   if (left == 0_I4P .and. right == 0_I4P) then
      depth_max = level
   else
      depth_max = max(tree_depth(left, level + 1_I4P), tree_depth(right, level + 1_I4P))
   endif
   endfunction tree_depth
endprogram fossil_test_bvh_depth_cap
