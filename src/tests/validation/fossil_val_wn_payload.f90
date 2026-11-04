!< VALIDATION version (issue #46): the test as it was before issue #45, on the large meshes (bunny, dragon-fine),
!< run by scripts/run_validation.sh at -O2. The fast regression variant is src/tests/fossil_test_wn_payload.f90.
!<
!< FOSSIL, test the persistent fast-winding-number cache `surface%wn_payload` (issue #29, step 2).
!<
!< For the SAH BVH and the octree, on cube-02, bunny and dragon-fine:
!<   1. built by analyze: allocated, indexed like aabb%node, every present node ready;
!<   2. each node: dipole = sum of area * normal over its own facets if it is a leaf, else the
!<      sum of its children's dipoles (recomputed independently here); centre = bbox centre;
!<      diag2 = squared bbox diagonal;
!<   3. root dipole = sum of area * normal over all facets, to 1e-10 of the total area (absolute:
!<      the sum is ~0 on closed surfaces; the bound covers summation-order rounding). This is the
!<      check that catches multi-counting: the octree stores every facet once per level, so a
!<      cache adding a node's own facets and its children's dipoles is off by ~6% on bunny;
!<   4. a second analyze rebuilds the same values; destroy deallocates the cache.

program fossil_val_wn_payload

use fossil,                  only : surface_stl_object, facet_object, AABB_TREE_SAH_BVH, AABB_TREE_OCTREE
use fossil_aabb_node_object, only : aabb_node_object
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

logical :: ok(6)

ok(1) = check('src/tests/cube-02.stl',     AABB_TREE_SAH_BVH, 'cube-02     BVH   ')
ok(2) = check('src/tests/bunny.stl',       AABB_TREE_SAH_BVH, 'bunny       BVH   ')
ok(3) = check('src/tests/dragon-fine.stl', AABB_TREE_SAH_BVH, 'dragon-fine BVH   ')
ok(4) = check('src/tests/cube-02.stl',     AABB_TREE_OCTREE,  'cube-02     octree')
ok(5) = check('src/tests/bunny.stl',       AABB_TREE_OCTREE,  'bunny       octree')
ok(6) = check('src/tests/dragon-fine.stl', AABB_TREE_OCTREE,  'dragon-fine octree')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(file_name, tree_kind, label) result(passed)
   character(*), intent(in)    :: file_name, label
   integer(I4P), intent(in)    :: tree_kind
   logical                     :: passed
   type(surface_stl_object)    :: s
   type(facet_object), pointer :: f
   type(aabb_node_object), pointer :: nd
   type(vector_R8P)            :: total, dip, bmin, bmax
   real(R8P)                   :: area_tot, d2, err_root
   real(R8P), allocatable      :: first(:,:)
   integer(I4P)                :: i, n, k, fid, nn, n_bad
   logical                     :: rebuilt_ok, destroyed_ok

   call s%load_from_file(file_name=file_name, guess_format=.true., aabb_tree_kind=tree_kind)
   nn = s%aabb%get_nodes_number()
   passed = allocated(s%wn_payload)
   if (.not. passed) then
      print '(A,A)', label, ': wn_payload not allocated'
      return
   endif
   passed = (lbound(s%wn_payload, 1) == 0 .and. ubound(s%wn_payload, 1) == nn - 1)

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

   allocate(first(3, 0:nn - 1))
   first(1, :) = s%wn_payload(:)%dx ; first(2, :) = s%wn_payload(:)%dy ; first(3, :) = s%wn_payload(:)%dz
   call s%analyze(aabb_tree_kind=tree_kind)
   rebuilt_ok = all(first(1, :) == s%wn_payload(:)%dx) .and. all(first(2, :) == s%wn_payload(:)%dy) .and. &
                all(first(3, :) == s%wn_payload(:)%dz)
   call s%destroy
   destroyed_ok = .not. allocated(s%wn_payload)

   passed = passed .and. (n_bad == 0) .and. (err_root <= 1.e-10_R8P * area_tot) .and. rebuilt_ok .and. destroyed_ok
   print '(A,A,I0,A,I0,A,ES10.3,A,L1,A,L1,A,L1)', label, ': nodes=', nn, ' bad=', n_bad, ' root err/area=', &
         err_root / area_tot, ' rebuilt=', rebuilt_ok, ' destroyed=', destroyed_ok, ' ok=', passed
   endfunction check

   function n_children(s, n) result(nchild)
   !< Number of present children of node n.
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: n
   integer(I4P)                         :: nchild, child_idx(8)

   call s%aabb%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
   endfunction n_children

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
endprogram fossil_val_wn_payload
