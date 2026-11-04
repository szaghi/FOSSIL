!< FOSSIL, test the persistent fast-winding-number cache `surface%wn_payload` (issue #29, step 2).
!<
!< For the SAH BVH and the octree, on cube-02 and dragon (bunny and dragon-fine dropped for the
!< per-push time budget, issue #45):
!<   1. built by analyze: allocated, indexed like aabb%node, every present node ready;
!<   2. each node: dipole = sum of area * normal over its own facets if it is a leaf, else the
!<      sum of its children's dipoles (recomputed independently here); centre = bbox centre;
!<      diag2 = squared bbox diagonal;
!<   3. root dipole = sum of area * normal over all facets, to 1e-10 of the total area (absolute:
!<      the sum is ~0 on closed surfaces; the bound covers summation-order rounding);
!<   4. a later analyze rebuilds the same tree and values as load; an octree has exactly
!<      1 + 8 + ... + 8**depth nodes (issue #46); destroy deallocates the cache.
!< Multi-counting (the octree stores every facet once per level; a cache adding a node's own facets
!< to its children's dipoles) is caught by check 2 at the non-root internal octree nodes, whose own
!< facets form an open patch (dipole O(10) on dragon): hence the asserted dragon octree depth >= 2.
!< Check 3 alone catches it only on an open mesh (~6% on bunny; ~1e-17 on the nearly closed dragon).

program fossil_test_wn_payload

use fossil,                  only : surface_stl_object, facet_object, AABB_TREE_SAH_BVH, AABB_TREE_OCTREE
use fossil_aabb_node_object, only : aabb_node_object
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

logical :: ok(4)

ok(1) = check('src/tests/cube-02.stl', AABB_TREE_SAH_BVH, 'cube-02 BVH   ', 0_I4P)
ok(2) = check('src/tests/dragon.stl',  AABB_TREE_SAH_BVH, 'dragon  BVH   ', 0_I4P)
ok(3) = check('src/tests/cube-02.stl', AABB_TREE_OCTREE,  'cube-02 octree', 0_I4P)
ok(4) = check('src/tests/dragon.stl',  AABB_TREE_OCTREE,  'dragon  octree', 2_I4P)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(file_name, tree_kind, label, min_depth) result(passed)
   character(*), intent(in)    :: file_name, label
   integer(I4P), intent(in)    :: tree_kind
   integer(I4P), intent(in)    :: min_depth !< Minimum tree depth (root = 0) required for the case to bite.
   logical                     :: passed
   type(surface_stl_object)    :: s
   type(facet_object), pointer :: f
   type(aabb_node_object), pointer :: nd
   type(vector_R8P)            :: total, dip, bmin, bmax
   real(R8P)                   :: area_tot, d2, err_root
   real(R8P), allocatable      :: first(:,:)
   integer(I4P)                :: i, n, k, fid, nn, nr, n_bad, depth
   logical                     :: rebuilt_ok, destroyed_ok

   call s%load_from_file(file_name=file_name, guess_format=.true., aabb_tree_kind=tree_kind)
   nn = s%aabb%get_nodes_number()
   passed = allocated(s%wn_payload)
   if (.not. passed) then
      print '(A,A)', label, ': wn_payload not allocated'
      return
   endif
   passed = (lbound(s%wn_payload, 1) == 0 .and. ubound(s%wn_payload, 1) == nn - 1)
   depth = tree_depth(s, 0_I4P)
   passed = passed .and. (depth >= min_depth)

   n_bad = 0
   do n = nn - 1, 0, -1
      nd => s%aabb%node_at(n)
      if (.not. nd%is_allocated()) cycle
      if (.not. s%wn_payload(n)%ready) n_bad = n_bad + 1
      bmin = nd%bmin() ; bmax = nd%bmax()
      d2 = (bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2
      if (s%wn_payload(n)%diag2 /= d2) n_bad = n_bad + 1
      if (s%wn_payload(n)%cx /= (bmin%x + bmax%x) * 0.5_R8P) n_bad = n_bad + 1
      dip = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      if (nd%has_facets() .and. n_children(s, n) == 0) then
         do k = 1, nd%facet_id_count()
            fid = nd%facet_id_at(k)
            f => s%facet_at(fid)
            dip = dip + f%normal * f%area()
         enddo
      endif
      call add_children(s, n, dip)
      if (s%wn_payload(n)%dx /= dip%x .or. s%wn_payload(n)%dy /= dip%y .or. s%wn_payload(n)%dz /= dip%z) &
         n_bad = n_bad + 1
   enddo

   total = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; area_tot = 0._R8P
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      total = total + f%normal * f%area()
      area_tot = area_tot + f%area()
   enddo
   err_root = max(abs(s%wn_payload(0)%dx - total%x), abs(s%wn_payload(0)%dy - total%y), &
                  abs(s%wn_payload(0)%dz - total%z))

   ! The rebuild check compares the tree built by load with the one a later analyze builds: the same nodes and
   ! values (issue #46: the octree kept the node array of the requested depth when the largest-edge rule lowered it,
   ! 585 nodes against 73 on dragon). An octree has exactly 1 + 8 + ... + 8**depth nodes.
   allocate(first(3, 0:nn - 1))
   first(1, :) = s%wn_payload(:)%dx ; first(2, :) = s%wn_payload(:)%dy ; first(3, :) = s%wn_payload(:)%dz
   if (tree_kind == AABB_TREE_OCTREE) passed = passed .and. (nn == (8_I4P**(depth + 1) - 1_I4P) / 7_I4P)
   call s%analyze(aabb_tree_kind=tree_kind)
   nr = s%aabb%get_nodes_number()
   rebuilt_ok = nr == nn
   if (rebuilt_ok) rebuilt_ok = all(first(1, :) == s%wn_payload(:)%dx) .and. &
                                all(first(2, :) == s%wn_payload(:)%dy) .and. all(first(3, :) == s%wn_payload(:)%dz)
   call s%destroy
   destroyed_ok = .not. allocated(s%wn_payload)

   passed = passed .and. (n_bad == 0) .and. (err_root <= 1.e-10_R8P * area_tot) .and. rebuilt_ok .and. destroyed_ok
   print '(A,A,I0,A,I0,A,I0,A,ES10.3,A,L1,A,L1,A,L1)', label, ': nodes=', nn, ' depth=', depth, ' bad=', n_bad, &
         ' root err/area=', err_root / area_tot, ' rebuilt=', rebuilt_ok, ' destroyed=', destroyed_ok, ' ok=', passed
   endfunction check

   function n_children(s, n) result(nchild)
   !< Number of present children of node n.
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: n
   integer(I4P)                         :: nchild, child_idx(8)

   call s%aabb%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
   endfunction n_children

   recursive function tree_depth(s, n) result(depth)
   !< Depth of the subtree rooted at node n (a leaf has depth 0).
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: n
   integer(I4P)                         :: depth, child_idx(8), nchild, c

   call s%aabb%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
   depth = 0
   do c = 1, nchild
      depth = max(depth, 1 + tree_depth(s, child_idx(c)))
   enddo
   endfunction tree_depth

   subroutine add_children(s, n, dip)
   !< Add the cached dipoles of node n's children, in enumerate_children order.
   type(surface_stl_object), intent(in)    :: s
   integer(I4P),             intent(in)    :: n
   type(vector_R8P),         intent(inout) :: dip
   integer(I4P)                            :: child_idx(8), nchild, c

   call s%aabb%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
   do c = 1, nchild
      if (.not. s%wn_payload(child_idx(c))%ready) cycle
      dip%x = dip%x + s%wn_payload(child_idx(c))%dx
      dip%y = dip%y + s%wn_payload(child_idx(c))%dy
      dip%z = dip%z + s%wn_payload(child_idx(c))%dz
   enddo
   endsubroutine add_children
endprogram fossil_test_wn_payload
