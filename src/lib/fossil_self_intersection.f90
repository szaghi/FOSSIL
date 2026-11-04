!< FOSSIL, self-intersection detection on a triangulated surface (issue #18 §1.2).

module fossil_self_intersection
!< FOSSIL, self-intersection detection on a triangulated surface.
!<
!< Detects pairs of facets in the same surface whose triangles cross
!< transversally (Möller tri-tri intersection). Reports the facet ID pair plus
!< the intersection segment endpoints. Facets sharing an edge meet only along
!< it (unless coplanar and folded over each other) and are skipped. Facets
!< sharing one vertex are not skipped (issue #46): they may still cross away
!< from it, and are classified exactly; touching only at the shared vertex is
!< not reported.
!<
!< Algorithm:
!<
!<   1. Broad phase: tree-vs-tree traversal over the surface's AABB tree against
!<      itself. At each node-pair we prune if the boxes don't overlap; otherwise
!<      we recurse into the larger node's children. At leaf-pairs we enumerate
!<      candidate facet-pairs (i, j) with i < j whose facet bboxes overlap.
!<      Tree-kind-agnostic — uses the public `enumerate_children` dispatcher.
!<
!<   2. Adjacency filter: count the shared vertices. Vertex identity is
!<      determined by the pool's `vertex_id` (when the pool is in use) or by
!<      coordinate equality within EPS otherwise. Two or three shared: skip.
!<
!<   3. Narrow phase. No shared vertex: `facet%intersect_facet` (Möller 1997),
!<      the segment endpoints of a transverse intersection. One shared vertex:
!<      the exact `pair_classify` (orient3d); a crossing (a common segment of
!<      positive length) is reported, a touch (the vertex alone) is not.
!<      Coplanar pairs sharing a vertex (unit normals parallel to 1e-12 rad)
!<      are still skipped: an overlap of two coplanar facets around a common
!<      vertex is not detected, nor a crossing at an angle below 1e-12 rad.
!<
!< Resolution (retriangulation along intersection segments and interior pruning
!< via the generalized winding number) is left to the §1.1 booleans PR, which
!< already needs the same machinery — see issue #18 for the rationale.

use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,     only : facet_object
use fossil_list_id_object,   only : list_id_object
use fossil_predicates,       only : orient3d
use fossil_tritri_exact,     only : pair_classify, PAIR_CROSSING
use fossil_utils,            only : EPS
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: intersection_pair_t
public :: find_self_intersections

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P  !< Buffer size for `enumerate_children` (octree fan-out).
real(R8P),    parameter :: COPLANAR_SIN2 = 1.e-24_R8P !< Facets sharing a vertex with sin^2(normal angle) below: coplanar.

type :: intersection_pair_t
   !< Single self-intersection record: an unordered pair of facet IDs (`a < b`)
   !< plus the intersection segment endpoints in world coordinates.
   integer(I4P)     :: a = 0_I4P  !< First facet id (always the smaller of the two).
   integer(I4P)     :: b = 0_I4P  !< Second facet id.
   type(vector_R8P) :: p          !< Segment start.
   type(vector_R8P) :: q          !< Segment end.
endtype intersection_pair_t

contains

   subroutine find_self_intersections(facet, tree, pairs, status)
   !< Find all self-intersecting facet pairs in the surface.
   !<
   !< On return:
   !<   - `pairs` is allocated with exactly the number of intersections found
   !<     (zero-length allocation if none, never unallocated)
   !<   - `status` is set to 0 on success
   !<
   !< Falls back to brute-force O(N^2) facet pairs when the AABB tree is not
   !< initialized; on a built tree the broad phase is O(N log N) on typical
   !< inputs.
   type(facet_object),                            intent(in)               :: facet(:) !< Facets array.
   type(aabb_tree_object),                        intent(in), target       :: tree     !< AABB tree built over `facet`.
   type(intersection_pair_t), allocatable,        intent(out)              :: pairs(:) !< Output list of intersection records.
   integer(I4P),                                  intent(out), optional    :: status   !< 0 on success.
   type(intersection_pair_t), allocatable                                  :: buf(:)
   integer(I4P)                                                            :: nfound

   if (present(status)) status = 0_I4P

   ! Start with a small accumulator and double on demand. Nearly all clean
   ! meshes finish with zero pairs, so the initial size is intentionally small.
   allocate(buf(16))
   nfound = 0_I4P

   if (tree%get_nodes_number() <= 0) then
      call brute_force_pairs(facet=facet, buf=buf, nfound=nfound)
   else
      call traverse_node_pair(facet=facet, tree=tree, na=0_I4P, nb=0_I4P, &
                              buf=buf, nfound=nfound)
   endif

   ! Trim the accumulator down to the exact size.
   allocate(pairs(nfound))
   if (nfound > 0_I4P) pairs(1:nfound) = buf(1:nfound)
   deallocate(buf)
   endsubroutine find_self_intersections

   subroutine brute_force_pairs(facet, buf, nfound)
   !< O(N^2) fallback used when no AABB tree is available. Same adjacency
   !< filter and narrow-phase as the BVH path.
   type(facet_object),                     intent(in)     :: facet(:)
   type(intersection_pair_t), allocatable, intent(inout)  :: buf(:)
   integer(I4P),                           intent(inout)  :: nfound
   integer(I4P)                                           :: i, j, n
   type(vector_R8P)                                       :: p, q
   logical                                                :: hit

   n = size(facet, dim=1)
   do i = 1, n - 1
      do j = i + 1, n
         if (.not. bboxes_overlap(facet(i), facet(j))) cycle
         call test_pair(fa=facet(i), fb=facet(j), hit=hit, p=p, q=q)
         if (hit) call push_pair(buf=buf, nfound=nfound, a=i, b=j, p=p, q=q)
      enddo
   enddo
   endsubroutine brute_force_pairs

   recursive subroutine traverse_node_pair(facet, tree, na, nb, buf, nfound)
   !< Tree-vs-tree traversal. Visits node-pairs `(na, nb)` with `na <= nb` to
   !< avoid double-counting. Prunes when the two AABBs don't overlap; recurses
   !< by descending the *larger* of the two nodes (volume-based heuristic).
   !<
   !< When `na == nb` we additionally enumerate the single-node "self" case to
   !< pick up intra-node facet pairs that don't reach a leaf-leaf descent.
   type(facet_object),                     intent(in)    :: facet(:)
   type(aabb_tree_object),                 intent(in)    :: tree
   integer(I4P),                           intent(in)    :: na, nb
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)
   integer(I4P),                           intent(inout) :: nfound
   type(aabb_node_object), pointer                       :: node_a, node_b
   integer(I4P)                                          :: ca(MAX_CHILDREN), cb(MAX_CHILDREN), nca, ncb
   integer(I4P)                                          :: i, j
   type(vector_R8P)                                      :: amin, amax, bmin, bmax
   real(R8P)                                             :: va, vb

   node_a => tree%node_at(i=na)
   node_b => tree%node_at(i=nb)
   if (.not. associated(node_a) .or. .not. associated(node_b)) return
   if (.not. node_a%is_allocated() .or. .not. node_b%is_allocated()) return

   ! Box-vs-box overlap.
   amin = node_a%bmin() ; amax = node_a%bmax()
   bmin = node_b%bmin() ; bmax = node_b%bmax()
   if (amax%x < bmin%x .or. bmax%x < amin%x .or. &
       amax%y < bmin%y .or. bmax%y < amin%y .or. &
       amax%z < bmin%z .or. bmax%z < amin%z) return

   call tree%enumerate_children(n=na, out_idx=ca, nchild=nca)
   call tree%enumerate_children(n=nb, out_idx=cb, nchild=ncb)

   if (nca == 0 .and. ncb == 0) then
      ! Leaf-leaf: enumerate candidate facet pairs.
      call enumerate_leaf_pair(facet=facet, tree=tree, na=na, nb=nb, &
                               buf=buf, nfound=nfound)
      return
   endif

   ! Descend the larger node. If only one of the two has children, descend that.
   if (nca == 0) then
      ! a is leaf; recurse into b's children paired with a.
      do j = 1, ncb
         call traverse_node_pair(facet=facet, tree=tree, na=na, nb=cb(j), &
                                 buf=buf, nfound=nfound)
      enddo
      return
   endif
   if (ncb == 0) then
      do i = 1, nca
         call traverse_node_pair(facet=facet, tree=tree, na=ca(i), nb=nb, &
                                 buf=buf, nfound=nfound)
      enddo
      return
   endif

   ! Both have children: descend the one with the larger bbox volume to keep
   ! the recursion tree balanced.
   va = (amax%x - amin%x) * (amax%y - amin%y) * (amax%z - amin%z)
   vb = (bmax%x - bmin%x) * (bmax%y - bmin%y) * (bmax%z - bmin%z)
   if (va >= vb) then
      do i = 1, nca
         call traverse_node_pair(facet=facet, tree=tree, na=ca(i), nb=nb, &
                                 buf=buf, nfound=nfound)
      enddo
   else
      do j = 1, ncb
         call traverse_node_pair(facet=facet, tree=tree, na=na, nb=cb(j), &
                                 buf=buf, nfound=nfound)
      enddo
   endif
   endsubroutine traverse_node_pair

   subroutine enumerate_leaf_pair(facet, tree, na, nb, buf, nfound)
   !< Enumerate candidate facet pairs at the leaf-pair (na, nb).
   !<
   !< When `na == nb` (the same leaf), enumerate i < j over its facet ids.
   !< When `na /= nb` (two different leaves), enumerate the cross-product but
   !< canonicalize each pair to (min, max) and reject pairs where i >= j to
   !< avoid double-counting (since the recursion may visit (na,nb) and (nb,na)).
   !< The recursion convention is na <= nb at the root, but during descent both
   !< orderings can arise — the canonicalization here is the safety net.
   type(facet_object),                     intent(in)    :: facet(:)
   type(aabb_tree_object),                 intent(in)    :: tree
   integer(I4P),                           intent(in)    :: na, nb
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)
   integer(I4P),                           intent(inout) :: nfound
   type(aabb_node_object), pointer                       :: node_a, node_b
   type(list_id_object)                                  :: ids_a, ids_b
   integer(I4P)                                          :: i, j, fi, fj, lo, hi
   type(vector_R8P)                                      :: p, q
   logical                                               :: hit

   node_a => tree%node_at(i=na)
   node_b => tree%node_at(i=nb)
   if (.not. node_a%has_facets() .or. .not. node_b%has_facets()) return

   ids_a = node_a%facet_id()
   ids_b = node_b%facet_id()

   do i = 1, ids_a%ids_number
      fi = ids_a%id(i)
      if (fi < 1 .or. fi > size(facet, dim=1)) cycle
      do j = 1, ids_b%ids_number
         fj = ids_b%id(j)
         if (fj < 1 .or. fj > size(facet, dim=1)) cycle
         if (fi == fj) cycle
         lo = min(fi, fj) ; hi = max(fi, fj)
         ! Canonical ordering — only test once per pair regardless of which
         ! leaf-pair direction we arrived from.
         if (na == nb .and. fi >= fj) cycle  ! intra-leaf: i < j only
         if (na /= nb .and. fi /= lo) cycle  ! inter-leaf: only the (lo, hi) direction
         if (.not. bboxes_overlap(facet(lo), facet(hi))) cycle
         call test_pair(fa=facet(lo), fb=facet(hi), hit=hit, p=p, q=q)
         if (hit) call push_pair(buf=buf, nfound=nfound, a=lo, b=hi, p=p, q=q)
      enddo
   enddo
   endsubroutine enumerate_leaf_pair

   pure function bboxes_overlap(fa, fb) result(yes)
   !< Cheap AABB overlap test for two facet bounding boxes (already cached on
   !< the facet via `compute_metrix`).
   type(facet_object), intent(in) :: fa, fb
   logical                        :: yes

   yes = .not. (fa%bb(2)%x < fb%bb(1)%x .or. fb%bb(2)%x < fa%bb(1)%x .or. &
                fa%bb(2)%y < fb%bb(1)%y .or. fb%bb(2)%y < fa%bb(1)%y .or. &
                fa%bb(2)%z < fb%bb(1)%z .or. fb%bb(2)%z < fa%bb(1)%z)
   endfunction bboxes_overlap

   subroutine test_pair(fa, fb, hit, p, q)
   !< Narrow phase of one candidate pair (see the module header, steps 2 and 3).
   type(facet_object), intent(in)  :: fa, fb !< The pair.
   logical,            intent(out) :: hit    !< True if the facets cross.
   type(vector_R8P),   intent(out) :: p, q   !< Intersection segment, if hit.
   integer(I4P)                    :: kind   !< PAIR_* code of the exact classification.

   hit = .false.
   select case (shared_vertices_number(fa, fb))
   case (0_I4P)
      call fa%intersect_facet(other=fb, p=p, q=q, intersects=hit)
   case (1_I4P)
      ! the shared vertex lies on both planes: if the other two vertices of either facet are strictly on one side
      ! of the other's plane, the facets meet at that vertex alone. This exact filter skips orient3d on the shared
      ! vertex (always zero, the slow exact path), which every pair of a fan around a vertex would take.
      ! coplanar (to 1e-12 rad: unit normals parallel, a vertex in common) are skipped, see the module header; the
      ! test is cheap, while orient3d on coplanar points takes its slow exact path (fans of flat faces)
      if ((fa%normal%y * fb%normal%z - fa%normal%z * fb%normal%y)**2 + &
          (fa%normal%z * fb%normal%x - fa%normal%x * fb%normal%z)**2 + &
          (fa%normal%x * fb%normal%y - fa%normal%y * fb%normal%x)**2 <= COPLANAR_SIN2) return
      if (off_plane_one_side(t=fa%vertex, s=fb%vertex)) return
      if (off_plane_one_side(t=fb%vertex, s=fa%vertex)) return
      call pair_classify(ta=fa%vertex, tb=fb%vertex, kind=kind, p=p, q=q)
      hit = kind == PAIR_CROSSING
   case default
      p = fa%vertex(1) ; q = p ! sharing an edge (or a duplicate): they meet along it
   endselect
   endsubroutine test_pair

   function off_plane_one_side(t, s) result(yes)
   !< True if the vertices of `s` not shared with `t` are all strictly on one side of the plane of `t` (exact).
   type(vector_R8P), intent(in) :: t(3) !< Triangle giving the plane.
   type(vector_R8P), intent(in) :: s(3) !< Triangle sharing one vertex with `t`.
   logical                      :: yes
   integer(I4P)                 :: i, npos, nneg, nfree
   real(R8P)                    :: o

   npos = 0_I4P ; nneg = 0_I4P ; nfree = 0_I4P
   do i = 1, 3
      if (vertex_coincident(s(i), t(1)) .or. vertex_coincident(s(i), t(2)) .or. vertex_coincident(s(i), t(3))) cycle
      nfree = nfree + 1_I4P
      o = orient3d(t(1), t(2), t(3), s(i))
      if (o > 0._R8P) npos = npos + 1_I4P
      if (o < 0._R8P) nneg = nneg + 1_I4P
   enddo
   yes = nfree > 0_I4P .and. (npos == nfree .or. nneg == nfree)
   endfunction off_plane_one_side

   pure function shared_vertices_number(fa, fb) result(n)
   !< Number of vertices of `fa` shared with `fb` (0 to 3).
   !<
   !< Fast path: when the vertex pool is in use, both facets carry non-zero
   !< `vertex_id` triples — equality of any pair of ids is exact and cheap.
   !< Slow path: when the pool is not in use (or any id is unset), fall back
   !< to coordinate comparison within EPS.
   type(facet_object), intent(in) :: fa, fb
   integer(I4P)                   :: n
   integer(I4P)                   :: i, j
   logical                        :: pool_a, pool_b

   n = 0_I4P
   pool_a = all(fa%vertex_id /= 0_I4P)
   pool_b = all(fb%vertex_id /= 0_I4P)
   do i = 1, 3
      do j = 1, 3
         if (pool_a .and. pool_b) then
            if (fa%vertex_id(i) == fb%vertex_id(j)) then
               n = n + 1_I4P ; exit
            endif
         elseif (vertex_coincident(fa%vertex(i), fb%vertex(j))) then
            n = n + 1_I4P ; exit
         endif
      enddo
   enddo
   endfunction shared_vertices_number

   pure function vertex_coincident(va, vb) result(yes)
   !< Coordinate-equality test for two vertices, within EPS.
   type(vector_R8P), intent(in) :: va, vb
   logical                      :: yes

   yes = abs(va%x - vb%x) <= EPS .and. &
         abs(va%y - vb%y) <= EPS .and. &
         abs(va%z - vb%z) <= EPS
   endfunction vertex_coincident

   subroutine push_pair(buf, nfound, a, b, p, q)
   !< Append a record to `buf`, doubling its capacity when full.
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)
   integer(I4P),                           intent(inout) :: nfound
   integer(I4P),                           intent(in)    :: a, b
   type(vector_R8P),                       intent(in)    :: p, q
   type(intersection_pair_t), allocatable                :: tmp(:)
   integer(I4P)                                          :: cap

   cap = size(buf, dim=1)
   if (nfound == cap) then
      allocate(tmp(2 * cap))
      tmp(1:cap) = buf(1:cap)
      call move_alloc(from=tmp, to=buf)
   endif
   nfound       = nfound + 1
   buf(nfound)%a = a
   buf(nfound)%b = b
   buf(nfound)%p = p
   buf(nfound)%q = q
   endsubroutine push_pair

endmodule fossil_self_intersection
