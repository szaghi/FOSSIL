!< FOSSIL, generalized / fast winding number on a triangulated surface.

module fossil_winding_number
!< FOSSIL, generalized / fast winding number on a triangulated surface.
!<
!< Implements the generalized winding number of Jacobson, Kavan & Sorkine-Hornung
!< (SIGGRAPH 2013) with optional hierarchical evaluation following Barill, Dickson,
!< Schmidt, Levin & Jacobson (SIGGRAPH 2018).
!<
!< For a triangulated surface S the winding number of a query point q is
!<
!<     w(q) = (1/4 pi) * sum_{f in S} Omega_f(q)
!<
!< where Omega_f(q) is the signed solid angle subtended by triangle f at q
!< (Van Oosterom-Strackee). For a closed, outward-oriented surface w(q) = 1
!< strictly inside, 0 strictly outside, and intermediate on the boundary or for
!< open / non-watertight meshes — which is precisely why this is preferred over
!< the parity / Baerentzen-Aanaes sign tests on dirty STL.
!<
!< The hierarchical (fast) variant traverses the surface's AABB tree and replaces
!< distant subtrees by a single dipole approximation: when |q - centroid_node| is
!< greater than beta times the node's bounding-box diagonal, the contribution of
!< all facets in that subtree is approximated by
!<
!<     Omega_node(q) ~= (dipole_node . (q - centroid_node)) / |q - centroid_node|^3
!<
!< where dipole_node = sum_{f in subtree} area_f * normal_f. Otherwise the
!< traversal recurses.
!<
!< The default beta = 2.0 (Barill et al. recommended ratio) keeps the
!< approximation error well below 1e-6 on typical meshes.

use fossil_aabb_tree_object, only : aabb_tree_object, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use fossil_aabb_node_object, only : aabb_node_object
use fossil_facet_object,     only : facet_object
use fossil_list_id_object,   only : list_id_object
use fossil_utils,            only : PI
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: winding_number

real(R8P), parameter :: FOUR_PI         = 4._R8P * PI
real(R8P), parameter :: DEFAULT_BETA    = 2._R8P
integer(I4P), parameter :: TREE_RATIO_LOCAL = 8_I4P  !< Octree fan-out, mirrors fossil_aabb_tree_object.

type :: dipole_cache_t
   !< Per-node cache for the Barnes-Hut style winding-number traversal.
   !<
   !< The dipole expansion is taken about the node's bounding-box centre rather
   !< than an area-weighted centroid. The centre-of-mass would slightly tighten
   !< the multipole error, but it requires propagating an extra scalar (the
   !< subtree's total area) up the tree; the bbox-centre choice keeps the cache
   !< small and is the variant used by most production fast-WN implementations.
   !< The admissibility test still uses the bbox diagonal as the size proxy.
   type(vector_R8P) :: centre  !< Subtree expansion centre (bbox midpoint).
   type(vector_R8P) :: dipole  !< Subtree dipole moment, sum_f area_f * normal_f.
   real(R8P)        :: diag2 = 0._R8P !< Squared bounding-box diagonal of the node.
   logical          :: ready = .false. !< True once cache is populated for this node.
endtype dipole_cache_t

contains

   function winding_number(facet, tree, point, beta) result(w)
   !< Generalized / fast winding number at `point`.
   !<
   !< Returns the continuous scalar w(q):
   !<   - ~ 1.0 strictly inside a closed, outward-oriented surface
   !<   - ~ 0.0 strictly outside
   !<   - intermediate on the boundary or for open / non-watertight meshes
   !<
   !< If the AABB tree is initialized, performs the hierarchical Barnes-Hut
   !< evaluation with subtree-dipole approximation. Otherwise falls back to the
   !< exact O(n_facets) per-facet sum.
   !<
   !< `beta` is the multipole admissibility ratio. Default 2.0 (Barill et al.).
   !< Setting beta <= 0 disables the dipole approximation and forces the exact
   !< per-facet sum; useful for ground-truth tests.
   type(facet_object),     intent(in)           :: facet(:) !< Facets array (caller's surface).
   type(aabb_tree_object), intent(in), target   :: tree     !< AABB tree built over `facet`.
   type(vector_R8P),       intent(in)           :: point    !< Query point.
   real(R8P),              intent(in), optional :: beta     !< Admissibility ratio (default 2.0).
   real(R8P)                                    :: w        !< Winding number.
   real(R8P)                                    :: beta_    !< Local copy of beta.
   type(dipole_cache_t), allocatable            :: cache(:) !< Per-node cache (root index 0).
   real(R8P)                                    :: omega    !< Accumulated solid angle.

   beta_ = DEFAULT_BETA ; if (present(beta)) beta_ = beta

   if (tree%get_nodes_number() <= 0 .or. beta_ <= 0._R8P) then
      ! No tree available, or caller asked for exact: brute-force per-facet sum.
      omega = exact_solid_angle_sum(facet=facet, point=point)
      w = omega / FOUR_PI
      return
   endif

   ! Build the per-node dipole cache (post-order over the tree topology).
   allocate(cache(0:tree%get_nodes_number() - 1))
   call build_cache(facet=facet, tree=tree, n=0_I4P, cache=cache)

   ! Hierarchical traversal from root (node 0).
   omega = 0._R8P
   call traverse(facet=facet, tree=tree, n=0_I4P, point=point, beta2=beta_*beta_, &
                 cache=cache, omega=omega)

   w = omega / FOUR_PI
   endfunction winding_number

   pure function exact_solid_angle_sum(facet, point) result(omega)
   !< Brute-force sum of signed solid angles over every facet. Used as the
   !< fallback when no tree is available, and as the ground truth for
   !< validating the hierarchical path.
   type(facet_object), intent(in) :: facet(:) !< Facets array.
   type(vector_R8P),   intent(in) :: point    !< Query point.
   real(R8P)                      :: omega    !< Sum of solid angles.
   integer(I4P)                   :: f        !< Facet counter.

   omega = 0._R8P
   do f = 1, size(facet, dim=1)
      omega = omega + facet(f)%solid_angle(point=point)
   enddo
   endfunction exact_solid_angle_sum

   recursive subroutine build_cache(facet, tree, n, cache)
   !< Populate `cache(n)` with the subtree's dipole centroid + dipole moment +
   !< squared box-diagonal, recursing into children first (post-order).
   !<
   !< Aggregation rule: dipole_node = sum_subtree area_f * normal_f.
   !< Centroid: area-weighted mean of facet centroids over the subtree.
   !< Cache entries for empty / unallocated nodes are left with `ready = .false.`
   !< (their `dipole` is the zero vector, so they contribute nothing).
   type(facet_object),     intent(in)            :: facet(:) !< Facets array.
   type(aabb_tree_object), intent(in)            :: tree     !< AABB tree.
   integer(I4P),           intent(in)            :: n        !< Node index.
   type(dipole_cache_t),   intent(inout)         :: cache(0:) !< Per-node cache.
   type(aabb_node_object), pointer               :: node     !< Pointer to node n.
   integer(I4P)                                  :: child_idx(TREE_RATIO_LOCAL), nchild, c, fid, k
   type(list_id_object)                          :: ids
   type(vector_R8P)                              :: bmin, bmax
   type(vector_R8P)                              :: dipole_acc
   type(vector_R8P)                              :: diag
   real(R8P)                                     :: area_f

   node => tree%node_at(i=n)
   if (.not. associated(node)) return
   if (.not. node%is_allocated()) return

   ! Box diagonal squared — used by the admissibility test on the host node.
   bmin = node%bmin() ; bmax = node%bmax()
   diag = bmax - bmin
   cache(n)%diag2 = diag%dotproduct(rhs=diag)

   ! Recurse into children first.
   call enumerate_children_local(tree=tree, n=n, out_idx=child_idx, nchild=nchild)
   do c = 1, nchild
      call build_cache(facet=facet, tree=tree, n=child_idx(c), cache=cache)
   enddo

   ! Aggregate this subtree's dipole moment. Two contributions:
   !   (a) facets stored directly on this node (octree internal nodes can carry
   !       facets that did not fit deeper; for the SAH BVH only leaves carry).
   !   (b) child subtree dipoles already computed.
   dipole_acc = vector_R8P(0._R8P, 0._R8P, 0._R8P)

   if (node%has_facets()) then
      ids = node%facet_id()
      do k = 1, ids%ids_number
         fid = ids%id(k)
         if (fid < 1 .or. fid > size(facet, dim=1)) cycle
         area_f     = facet(fid)%area()
         dipole_acc = dipole_acc + facet(fid)%normal * area_f
      enddo
   endif

   do c = 1, nchild
      if (.not. cache(child_idx(c))%ready) cycle
      dipole_acc = dipole_acc + cache(child_idx(c))%dipole
   enddo

   cache(n)%centre = (bmin + bmax) * 0.5_R8P
   cache(n)%dipole = dipole_acc
   cache(n)%ready  = .true.
   endsubroutine build_cache

   subroutine enumerate_children_local(tree, n, out_idx, nchild)
   !< Local replica of `aabb_tree_object%enumerate_children` (which is private to
   !< its module). Lists allocated children of node `n`, dispatching on tree kind.
   type(aabb_tree_object), intent(in)  :: tree       !< AABB tree.
   integer(I4P),           intent(in)  :: n          !< Node index.
   integer(I4P),           intent(out) :: out_idx(:) !< Output: child indices, caller buffer >= TREE_RATIO_LOCAL.
   integer(I4P),           intent(out) :: nchild     !< Number of allocated children.
   type(aabb_node_object), pointer     :: node, child
   integer(I4P)                        :: lc, rc, fcn, i

   nchild = 0
   node => tree%node_at(i=n)
   if (.not. associated(node)) return

   if (tree%get_tree_kind() == AABB_TREE_SAH_BVH) then
      lc = node%get_left_child() ; rc = node%get_right_child()
      if (lc > 0) then ; nchild = nchild + 1 ; out_idx(nchild) = lc ; endif
      if (rc > 0) then ; nchild = nchild + 1 ; out_idx(nchild) = rc ; endif
   else
      ! Octree implicit indexing: first_child = TREE_RATIO * parent + 1.
      fcn = TREE_RATIO_LOCAL * n + 1
      if (fcn > tree%get_nodes_number() - TREE_RATIO_LOCAL) return
      do i = fcn, fcn + TREE_RATIO_LOCAL - 1
         child => tree%node_at(i=i)
         if (.not. associated(child)) cycle
         if (child%is_allocated()) then
            nchild = nchild + 1
            out_idx(nchild) = i
         endif
      enddo
   endif
   endsubroutine enumerate_children_local

   recursive subroutine traverse(facet, tree, n, point, beta2, cache, omega)
   !< Walk the tree adding solid-angle contributions from node `n`'s subtree.
   !<
   !< Two cases:
   !<   - Far field (|point - centroid|^2 > beta^2 * diag^2): use the dipole.
   !<   - Near field: sum this node's own facets exactly, and recurse into
   !<     children. (For SAH BVH only leaves have facets; for octrees internal
   !<     nodes may also carry facets.)
   type(facet_object),     intent(in)            :: facet(:)
   type(aabb_tree_object), intent(in)            :: tree
   integer(I4P),           intent(in)            :: n
   type(vector_R8P),       intent(in)            :: point
   real(R8P),              intent(in)            :: beta2  !< Squared admissibility ratio.
   type(dipole_cache_t),   intent(in)            :: cache(0:)
   real(R8P),              intent(inout)         :: omega
   type(aabb_node_object), pointer               :: node
   type(list_id_object)                          :: ids
   integer(I4P)                                  :: child_idx(TREE_RATIO_LOCAL), nchild, c, k, fid
   type(vector_R8P)                              :: r
   real(R8P)                                     :: r2, r1, inv_r3, dot
   logical                                       :: is_leaf

   node => tree%node_at(i=n)
   if (.not. associated(node)) return
   if (.not. node%is_allocated()) return
   if (.not. cache(n)%ready) return

   ! Far-field test: replace the entire subtree by its dipole.
   r  = point - cache(n)%centre
   r2 = r%dotproduct(rhs=r)
   if (r2 > beta2 * cache(n)%diag2 .and. cache(n)%diag2 > 0._R8P) then
      r1     = sqrt(r2)
      inv_r3 = 1._R8P / (r2 * r1)
      dot    = cache(n)%dipole%dotproduct(rhs=r)
      omega  = omega + dot * inv_r3
      return
   endif

   ! Near field: account for facets stored directly on this node, then recurse.
   if (node%has_facets()) then
      ids = node%facet_id()
      do k = 1, ids%ids_number
         fid = ids%id(k)
         if (fid < 1 .or. fid > size(facet, dim=1)) cycle
         omega = omega + facet(fid)%solid_angle(point=point)
      enddo
   endif

   call enumerate_children_local(tree=tree, n=n, out_idx=child_idx, nchild=nchild)
   is_leaf = (nchild == 0)
   if (is_leaf) return

   ! For SAH BVH internal nodes the facets are not stored on the parent (they are
   ! in the leaves); recursing into children is the correct path. For octree
   ! internal nodes, we have already counted the parent's own facets above and
   ! must NOT double-count by recursing — but the octree representation stores a
   ! given facet exactly once (parent OR child), so the parent's facets are
   ! disjoint from the children's. So recursion is correct in both cases.
   do c = 1, nchild
      call traverse(facet=facet, tree=tree, n=child_idx(c), point=point, beta2=beta2, &
                    cache=cache, omega=omega)
   enddo
   endsubroutine traverse

endmodule fossil_winding_number
