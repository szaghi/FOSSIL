!< FOSSIL, self-intersection detection on a triangulated surface (issue #18 §1.2).

module fossil_self_intersection
!< FOSSIL, self-intersection detection on a triangulated surface.
!<
!< Detects pairs of facets in the same surface whose triangles cross, or
!< overlap in a common plane, decided exactly. Reports the facet ID pair, the
!< kind (SELFINT_*) and the intersection segment endpoints. Facets sharing an edge meet only along
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
!<      exact coordinate equality otherwise. Two or three shared: skip.
!<
!<   3. Narrow phase. No shared vertex (issue #69): `pair_cut` (exact orient3d
!<      signs, implicit section points; it replaced the floating Möller test,
!<      which skipped every coplanar overlap, and vertices EPS apart no longer
!<      count as one): a crossing (common segment of positive
!<      length) is reported, a touch (one common point) only with
!<      `include_touching`; an exactly coplanar pair is decided in the
!<      projection (exact orient2d): overlapping interiors are reported
!<      (SELFINT_COPLANAR), boundaries meeting only with `include_touching`.
!<      One shared vertex:
!<      an exact wedge test on orient3d signs of input vertices (issue #55; it
!<      replaced pair_classify, whose on-plane snap reported crossings on
!<      near-coplanar pairs touching at the vertex only); a crossing (a common
!<      segment of positive length, along an edge included) is reported, a
!<      touch (the vertex alone) is not.
!<      Exactly coplanar pairs sharing a vertex (issue #69; skipped to 1e-12
!<      rad before) are decided in the projection: overlapping interiors are
!<      reported (SELFINT_COPLANAR), a common edge segment from the vertex is a
!<      touch (include_touching); the vertex alone is adjacency.
!<
!< Resolution (retriangulation along intersection segments and interior pruning
!< via the generalized winding number) is left to the §1.1 booleans PR, which
!< already needs the same machinery — see issue #18 for the rationale.

use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object,     only : facet_object
use fossil_list_id_object,   only : list_id_object
use fossil_implicit_points,  only : implicit_point_t, ip_approx
use fossil_predicates,       only : orient2d_xy, orient3d
use fossil_tritri_exact,     only : pair_cut, facet_drop, PAIR_CROSSING, PAIR_TOUCHING, PAIR_COPLANAR
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: intersection_pair_t
public :: find_self_intersections

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P  !< Buffer size for `enumerate_children` (octree fan-out).
integer(I4P), parameter, public :: SELFINT_CROSSING = 1_I4P !< The facets share a segment of positive length.
integer(I4P), parameter, public :: SELFINT_TOUCHING = 2_I4P !< They meet at one point, or along a coplanar boundary.
integer(I4P), parameter, public :: SELFINT_COPLANAR = 3_I4P !< Coplanar, interiors overlapping (positive area).
integer(I4P), parameter, public :: SELFINT_EDGE_FOLD = 4_I4P !< Sharing an edge, coplanar, folded over it.
integer(I4P), parameter, public :: SELFINT_SHEETS = 5_I4P    !< Adjacent around a non-manifold edge with the same direction:
                                                             !< the sheets through the edge cross there.
integer(I4P), parameter :: COPLANAR_APART = 0_I4P, COPLANAR_TOUCH = 1_I4P, COPLANAR_OVERLAP = 2_I4P

type :: intersection_pair_t
   !< Single self-intersection record: an unordered pair of facet IDs (`a < b`)
   !< plus the intersection segment endpoints in world coordinates.
   integer(I4P)     :: a = 0_I4P  !< First facet id (always the smaller of the two).
   integer(I4P)     :: b = 0_I4P  !< Second facet id.
   integer(I4P)     :: kind = 1_I4P !< SELFINT_CROSSING, SELFINT_TOUCHING (include_touching), SELFINT_COPLANAR,
                                    !< SELFINT_EDGE_FOLD, SELFINT_SHEETS.
   type(vector_R8P) :: p          !< Segment start.
   type(vector_R8P) :: q          !< Segment end.
endtype intersection_pair_t

contains

   subroutine find_self_intersections(facet, tree, pairs, status, include_touching)
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
   logical,                                       intent(in),  optional    :: include_touching !< Also report pairs
                                                                                    !< meeting at one point (default no).
   type(intersection_pair_t), allocatable                                  :: buf(:)
   integer(I4P)                                                            :: nfound
   logical                                                                 :: touching

   if (present(status)) status = 0_I4P
   touching = .false. ; if (present(include_touching)) touching = include_touching

   ! Start with a small accumulator and double on demand. Nearly all clean
   ! meshes finish with zero pairs, so the initial size is intentionally small.
   allocate(buf(16))
   nfound = 0_I4P

   if (tree%get_nodes_number() <= 0) then
      call brute_force_pairs(facet=facet, touching=touching, buf=buf, nfound=nfound)
   else
      call traverse_node_pair(facet=facet, tree=tree, na=0_I4P, nb=0_I4P, touching=touching, &
                              buf=buf, nfound=nfound)
   endif

   ! Sheets crossing along non-manifold edges: a pass over the edge fans (issue #69).
   call edge_fans(facet=facet, buf=buf, nfound=nfound)

   ! Trim the accumulator down to the exact size.
   allocate(pairs(nfound))
   if (nfound > 0_I4P) pairs(1:nfound) = buf(1:nfound)
   deallocate(buf)
   endsubroutine find_self_intersections

   subroutine brute_force_pairs(facet, touching, buf, nfound)
   !< O(N^2) fallback used when no AABB tree is available. Same adjacency
   !< filter and narrow-phase as the BVH path.
   type(facet_object),                     intent(in)     :: facet(:)
   logical,                                intent(in)     :: touching
   type(intersection_pair_t), allocatable, intent(inout)  :: buf(:)
   integer(I4P),                           intent(inout)  :: nfound
   integer(I4P)                                           :: i, j, n, kind
   type(vector_R8P)                                       :: p, q
   logical                                                :: hit

   n = size(facet, dim=1)
   do i = 1, n - 1
      do j = i + 1, n
         if (.not. bboxes_overlap(facet(i), facet(j))) cycle
         call test_pair(fa=facet(i), fb=facet(j), touching=touching, hit=hit, kind=kind, p=p, q=q)
         if (hit) call push_pair(buf=buf, nfound=nfound, a=i, b=j, kind=kind, p=p, q=q)
      enddo
   enddo
   endsubroutine brute_force_pairs

   recursive subroutine traverse_node_pair(facet, tree, na, nb, touching, buf, nfound)
   !< Tree-vs-tree traversal. Visits node-pairs `(na, nb)` with `na <= nb` to
   !< avoid double-counting. Prunes when the two AABBs don't overlap; recurses
   !< by descending the *larger* of the two nodes (volume-based heuristic).
   !<
   !< When `na == nb` we additionally enumerate the single-node "self" case to
   !< pick up intra-node facet pairs that don't reach a leaf-leaf descent.
   type(facet_object),                     intent(in)    :: facet(:)
   type(aabb_tree_object),                 intent(in)    :: tree
   integer(I4P),                           intent(in)    :: na, nb
   logical,                                intent(in)    :: touching
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
      call enumerate_leaf_pair(facet=facet, tree=tree, na=na, nb=nb, touching=touching, &
                               buf=buf, nfound=nfound)
      return
   endif

   ! Descend the larger node. If only one of the two has children, descend that.
   if (nca == 0) then
      ! a is leaf; recurse into b's children paired with a.
      do j = 1, ncb
         call traverse_node_pair(facet=facet, tree=tree, na=na, nb=cb(j), &
                                 touching=touching, buf=buf, nfound=nfound)
      enddo
      return
   endif
   if (ncb == 0) then
      do i = 1, nca
         call traverse_node_pair(facet=facet, tree=tree, na=ca(i), nb=nb, &
                                 touching=touching, buf=buf, nfound=nfound)
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
                                 touching=touching, buf=buf, nfound=nfound)
      enddo
   else
      do j = 1, ncb
         call traverse_node_pair(facet=facet, tree=tree, na=na, nb=cb(j), &
                                 touching=touching, buf=buf, nfound=nfound)
      enddo
   endif
   endsubroutine traverse_node_pair

   subroutine enumerate_leaf_pair(facet, tree, na, nb, touching, buf, nfound)
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
   logical,                                intent(in)    :: touching
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)
   integer(I4P),                           intent(inout) :: nfound
   type(aabb_node_object), pointer                       :: node_a, node_b
   type(list_id_object)                                  :: ids_a, ids_b
   integer(I4P)                                          :: i, j, fi, fj, lo, hi, kind
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
         call test_pair(fa=facet(lo), fb=facet(hi), touching=touching, hit=hit, kind=kind, p=p, q=q)
         if (hit) call push_pair(buf=buf, nfound=nfound, a=lo, b=hi, kind=kind, p=p, q=q)
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

   subroutine test_pair(fa, fb, touching, hit, kind, p, q)
   !< Narrow phase of one candidate pair (see the module header, steps 2 and 3).
   type(facet_object), intent(in)  :: fa, fb   !< The pair.
   logical,            intent(in)  :: touching !< Report pairs meeting at one point.
   logical,            intent(out) :: hit      !< True if the pair is reported.
   integer(I4P),       intent(out) :: kind     !< SELFINT_* code, if hit.
   type(vector_R8P),   intent(out) :: p, q     !< Intersection segment (see intersection_pair_t), if hit.
   logical                         :: coplanar

   hit = .false. ; kind = SELFINT_CROSSING
   p = fa%vertex(1) ; q = p
   select case (shared_vertices_number(fa, fb))
   case (0_I4P)
      call disjoint_pair(ta=fa%vertex, tb=fb%vertex, touching=touching, hit=hit, kind=kind, p=p, q=q)
   case (1_I4P)
      ! the shared vertex lies on both planes: if the other two vertices of either facet are strictly on one side
      ! of the other's plane, the facets meet at that vertex alone. This exact filter skips orient3d on the shared
      ! vertex (always zero, the slow exact path), which every pair of a fan around a vertex would take. A pair on
      ! one axis-aligned plane is coplanar without orient3d (flat fans: orient3d's zeros take its slow exact path).
      coplanar = axis_coplanar(fa%vertex, fb%vertex)
      if (.not. coplanar) then
         if (off_plane_one_side(t=fa%vertex, s=fb%vertex)) return
         if (off_plane_one_side(t=fb%vertex, s=fa%vertex)) return
         call shared_vertex_crossing(ta=fa%vertex, tb=fb%vertex, hit=hit, p=p, q=q, coplanar=coplanar)
      endif
      if (coplanar) then
         ! exactly coplanar (issue #69; skipped to 1e-12 rad before): a defect iff the interiors overlap, a touch if
         ! they meet beyond the shared vertex along a common edge segment (with include_touching)
         p = (fa%vertex(1) + fa%vertex(2) + fa%vertex(3)) / 3._R8P
         q = (fb%vertex(1) + fb%vertex(2) + fb%vertex(3)) / 3._R8P
         if (.not. touching .and. .not. boxes_interiors_meet(fa%vertex, fb%vertex)) return
         select case (coplanar_relation(fa%vertex, fb%vertex))
         case (COPLANAR_OVERLAP)
            hit = .true. ; kind = SELFINT_COPLANAR
         case (COPLANAR_TOUCH)
            if (touching .and. edge_overlap_at_vertex(fa%vertex, fb%vertex)) then
               hit = .true. ; kind = SELFINT_TOUCHING
            endif
         endselect
      endif
   case (2_I4P)
      ! sharing an edge: they meet along it, a defect only if coplanar and folded over it (issue #69); sheets crossing
      ! at a non-manifold edge are found by `edge_fans`
      call edge_fold(ta=fa%vertex, tb=fb%vertex, hit=hit, p=p, q=q)
      if (hit) kind = SELFINT_EDGE_FOLD
   case default
      p = fa%vertex(1) ; q = p ! a duplicate facet (three shared vertices): left to remove_duplicate_facets
   endselect
   endsubroutine test_pair

   subroutine edge_fold(ta, tb, hit, p, q)
   !< Two triangles sharing an edge (e0, e1) are folded over it when they are exactly coplanar (orient3d of the free
   !< vertex of tb on the plane of ta is zero) and their free vertices lie on the same side of the edge line in that
   !< plane (exact orient2d in the projection that keeps ta non-degenerate): their interiors overlap.
   type(vector_R8P), intent(in)  :: ta(3), tb(3) !< The triangles.
   logical,          intent(out) :: hit          !< Folded.
   type(vector_R8P), intent(out) :: p, q         !< The shared edge, if folded.
   integer(I4P)                  :: fa, fb, i, j, drop, u, v, n
   logical                       :: sa(3), sb(3)
   real(R8P)                     :: e0(3), e1(3), wa(3), wb(3)

   hit = .false. ; p = ta(1) ; q = p
   sa = .false. ; sb = .false. ; n = 0_I4P
   do i = 1, 3
      do j = 1, 3
         if (vertex_coincident(ta(i), tb(j))) then
            sa(i) = .true. ; sb(j) = .true. ; n = n + 1_I4P
         endif
      enddo
   enddo
   if (n /= 2_I4P) return
   fa = findloc(sa, .false., dim=1) ; fb = findloc(sb, .false., dim=1)
   if (.not. axis_coplanar(ta, tb)) then
      if (orient3d(ta(1), ta(2), ta(3), tb(fb)) /= 0._R8P) return
   endif
   drop = facet_drop(ta)
   if (drop == 0_I4P) return
   u = mod(drop, 3_I4P) + 1_I4P ; v = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   i = mod(fa, 3_I4P) + 1_I4P ; j = mod(fa + 1_I4P, 3_I4P) + 1_I4P
   e0 = [ta(i)%x, ta(i)%y, ta(i)%z] ; e1 = [ta(j)%x, ta(j)%y, ta(j)%z]
   wa = [ta(fa)%x, ta(fa)%y, ta(fa)%z] ; wb = [tb(fb)%x, tb(fb)%y, tb(fb)%z]
   i = o2([e0(u), e0(v)], [e1(u), e1(v)], [wa(u), wa(v)])
   j = o2([e0(u), e0(v)], [e1(u), e1(v)], [wb(u), wb(v)])
   hit = i /= 0_I4P .and. i == j
   if (hit) then
      p%x = e0(1) ; p%y = e0(2) ; p%z = e0(3) ; q%x = e1(1) ; q%y = e1(2) ; q%z = e1(3)
   endif
   endsubroutine edge_fold

   subroutine edge_fans(facet, buf, nfound)
   !< Sheets crossing along a non-manifold edge (issue #69). Around an edge (p0, p1) the incident facets are ordered by
   !< the angle of their half-planes, exactly (orient3d on input vertices), and each traverses the edge p0 -> p1
   !< (+1) or p1 -> p0 (-1). Embedded oriented sheets through the edge alternate around it (the winding number of
   !< the wedges stays in {w, w + 1}); two cyclically adjacent facets with the same direction bound a wedge wound
   !< twice: the sheets cross there, and the pair is reported (SELFINT_SHEETS, p-q the edge). Only fans of 4 or more
   !< facets with as many +1 as -1 are examined (an unbalanced fan is an orientation or boundary defect, not a
   !< crossing), and a fan with two facets in the same half-plane is skipped (coplanar folds, reported pairwise).
   !< Vertex identity: exact coordinates.
   type(facet_object),                     intent(in)    :: facet(:) !< Facets.
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)   !< Accumulator.
   integer(I4P),                           intent(inout) :: nfound   !< Records in buf.
   real(R8P),    allocatable                             :: xyz(:,:)
   integer(I4P), allocatable                             :: vid(:), ord(:), ekey(:,:), eord(:)
   integer(I4P)                                          :: n, i, k, a, b, first, last

   n = size(facet, dim=1)
   if (n < 4_I4P) return
   allocate(xyz(3, 3 * n), vid(3 * n))
   do i = 1, n
      do k = 1, 3
         xyz(:, 3 * (i - 1) + k) = [facet(i)%vertex(k)%x, facet(i)%vertex(k)%y, facet(i)%vertex(k)%z]
      enddo
   enddo
   ! vertex ids: equal coordinates, equal id
   call sort_columns_real(xyz, ord)
   k = 0_I4P
   do i = 1, 3 * n
      if (i == 1) then
         k = 1_I4P
      elseif (any(xyz(:, ord(i)) /= xyz(:, ord(i - 1)))) then
         k = k + 1_I4P
      endif
      vid(ord(i)) = k
   enddo
   ! edges (min id, max id, facet, direction +1 if the facet goes min -> max)
   allocate(ekey(4, 3 * n))
   do i = 1, n
      do k = 1, 3
         a = vid(3 * (i - 1) + k) ; b = vid(3 * (i - 1) + mod(k, 3_I4P) + 1_I4P)
         ekey(:, 3 * (i - 1) + k) = [min(a, b), max(a, b), i, merge(1_I4P, -1_I4P, a < b)]
      enddo
   enddo
   call sort_columns_int(ekey(1:3, :), eord)
   first = 1_I4P
   do while (first <= 3 * n)
      last = first
      do while (last < 3 * n)
         if (any(ekey(1:2, eord(last + 1)) /= ekey(1:2, eord(first)))) exit
         last = last + 1_I4P
      enddo
      if (last - first + 1 >= 4_I4P .and. ekey(1, eord(first)) /= ekey(2, eord(first))) &
         call one_fan(xyz, vid, ekey, eord(first:last), buf, nfound)
      first = last + 1_I4P
   enddo
   endsubroutine edge_fans

   subroutine one_fan(xyz, vid, ekey, slots, buf, nfound)
   !< Order the facets of one edge fan around the edge and report adjacent pairs with the same direction (edge_fans).
   real(R8P),                              intent(in)    :: xyz(:,:)  !< Vertex coordinates, 3 slots per facet.
   integer(I4P),                           intent(in)    :: vid(:)    !< Vertex ids of the slots.
   integer(I4P),                           intent(in)    :: ekey(:,:) !< Edge records (min id, max id, facet, dir).
   integer(I4P),                           intent(in)    :: slots(:)  !< The edge records of this fan.
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)    !< Accumulator.
   integer(I4P),                           intent(inout) :: nfound    !< Records in buf.
   integer(I4P)                                          :: m, r, t, u, c, tmp
   integer(I4P)                                          :: half(size(slots)), fan(size(slots)), dirs(size(slots))
   real(R8P)                                             :: p0(3), p1(3), w(3, size(slots)), aux(3), wt(3)
   type(vector_R8P)                                      :: vp, vq
   logical                                               :: tie
   real(R8P)                                             :: sg

   m = size(slots)
   if (sum(ekey(4, slots)) /= 0_I4P) return
   do t = 1, m
      fan(t) = ekey(3, slots(t)) ; dirs(t) = ekey(4, slots(t))
      do u = 1, 3
         c = 3 * (fan(t) - 1) + u
         if (vid(c) == ekey(1, slots(t))) then
            p0 = xyz(:, c)
         elseif (vid(c) == ekey(2, slots(t))) then
            p1 = xyz(:, c)
         else
            w(:, t) = xyz(:, c)
         endif
      enddo
   enddo
   ! halves about the reference half-plane (facet 1): 1 (0, pi), 2 (pi, 2 pi), 0 on it, 3 opposite (pi)
   aux = p0 + cross3(p1 - p0, w(:, 1) - p0)
   if (o3(p0, p1, w(:, 1), aux) == 0._R8P) return
   if (o3(p0, p1, aux, w(:, 1)) == 0._R8P) return
   half(1) = 0_I4P
   do t = 2, m
      sg = o3(p0, p1, w(:, 1), w(:, t))
      if (sg > 0._R8P) then
         half(t) = 1_I4P
      elseif (sg < 0._R8P) then
         half(t) = 2_I4P
      elseif (o3(p0, p1, aux, w(:, t)) == o3(p0, p1, aux, w(:, 1))) then
         return ! the same half-plane as facet 1: a tie
      else
         half(t) = 3_I4P
      endif
   enddo
   ! insertion sort around the edge; any tie (two facets in one half-plane) skips the fan
   do t = 2, m
      u = t
      do while (u > 1)
         if (.not. fan_before(half(u), half(u - 1), p0, p1, w(:, u), w(:, u - 1), tie)) then
            if (tie) return
            exit
         endif
         tmp = fan(u) ; fan(u) = fan(u - 1) ; fan(u - 1) = tmp
         tmp = dirs(u) ; dirs(u) = dirs(u - 1) ; dirs(u - 1) = tmp
         tmp = half(u) ; half(u) = half(u - 1) ; half(u - 1) = tmp
         wt = w(:, u) ; w(:, u) = w(:, u - 1) ; w(:, u - 1) = wt
         u = u - 1
      enddo
   enddo
   do t = 2, m
      if (.not. fan_before(half(t - 1), half(t), p0, p1, w(:, t - 1), w(:, t), tie)) return
   enddo
   vp%x = p0(1) ; vp%y = p0(2) ; vp%z = p0(3) ; vq%x = p1(1) ; vq%y = p1(2) ; vq%z = p1(3)
   do t = 1, m
      r = mod(t, m) + 1_I4P
      if (dirs(t) == dirs(r)) call push_pair(buf=buf, nfound=nfound, a=min(fan(t), fan(r)), b=max(fan(t), fan(r)), &
                                             kind=SELFINT_SHEETS, p=vp, q=vq)
   enddo
   endsubroutine one_fan

   function fan_before(h1, h2, p0, p1, w1, w2, tie) result(yes)
   !< Free vertex w1 (half h1) strictly before w2 (half h2) around the edge p0-p1; tie: the same half-plane.
   integer(I4P), intent(in)  :: h1, h2           !< Halves (0 reference, 1 (0, pi), 3 pi, 2 (pi, 2 pi)).
   real(R8P),    intent(in)  :: p0(3), p1(3)     !< Edge.
   real(R8P),    intent(in)  :: w1(3), w2(3)     !< Free vertices.
   logical,      intent(out) :: tie              !< Same half-plane.
   logical                   :: yes
   integer(I4P), parameter   :: RANKS(0:3) = [0_I4P, 1_I4P, 3_I4P, 2_I4P]
   real(R8P)                 :: o

   tie = .false.
   if (RANKS(h1) /= RANKS(h2)) then
      yes = RANKS(h1) < RANKS(h2) ; return
   endif
   if (h1 == 0_I4P .or. h1 == 3_I4P) then
      tie = .true. ; yes = .false. ; return
   endif
   o = o3(p0, p1, w1, w2)
   tie = o == 0._R8P
   yes = o > 0._R8P
   endfunction fan_before

   function o3(a1, a2, a3, a4) result(sg)
   !< Exact sign of orient3d as -1, 0, 1 (real), on plain coordinates.
   real(R8P), intent(in) :: a1(3), a2(3), a3(3), a4(3)
   real(R8P)             :: sg
   type(vector_R8P)      :: v1, v2, v3, v4
   real(R8P)             :: d

   v1%x = a1(1) ; v1%y = a1(2) ; v1%z = a1(3) ; v2%x = a2(1) ; v2%y = a2(2) ; v2%z = a2(3)
   v3%x = a3(1) ; v3%y = a3(2) ; v3%z = a3(3) ; v4%x = a4(1) ; v4%y = a4(2) ; v4%z = a4(3)
   d = orient3d(v1, v2, v3, v4)
   sg = 0._R8P
   if (d > 0._R8P) sg = 1._R8P
   if (d < 0._R8P) sg = -1._R8P
   endfunction o3

   pure function cross3(x, y) result(z)
   !< Cross product.
   real(R8P), intent(in) :: x(3), y(3)
   real(R8P)             :: z(3)

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross3

   subroutine sort_columns_real(a, ord)
   !< Stable merge sort of the columns of a, lexicographic: ord(1) is the smallest column.
   real(R8P),                 intent(in)  :: a(:,:) !< Columns to sort.
   integer(I4P), allocatable, intent(out) :: ord(:) !< Permutation.
   integer(I4P), allocatable              :: tmp(:)
   integer(I4P)                           :: n, width, lo, mid, hi, i, j, k

   n = size(a, dim=2)
   allocate(ord(n), tmp(n))
   ord = [(i, i = 1, n)]
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (col_less_real(a(:, ord(j)), a(:, ord(i)))) then
               tmp(k) = ord(j) ; j = j + 1_I4P
            else
               tmp(k) = ord(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         if (i <= mid) tmp(k:hi) = ord(i:mid)
         if (j <= hi) tmp(k:hi) = ord(j:hi)
         ord(lo:hi) = tmp(lo:hi)
         lo = lo + 2_I4P * width
      enddo
      width = 2_I4P * width
   enddo
   endsubroutine sort_columns_real

   subroutine sort_columns_int(a, ord)
   !< Stable merge sort of the columns of a, lexicographic: ord(1) is the smallest column.
   integer(I4P),              intent(in)  :: a(:,:) !< Columns to sort.
   integer(I4P), allocatable, intent(out) :: ord(:) !< Permutation.
   integer(I4P), allocatable              :: tmp(:)
   integer(I4P)                           :: n, width, lo, mid, hi, i, j, k

   n = size(a, dim=2)
   allocate(ord(n), tmp(n))
   ord = [(i, i = 1, n)]
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (col_less_int(a(:, ord(j)), a(:, ord(i)))) then
               tmp(k) = ord(j) ; j = j + 1_I4P
            else
               tmp(k) = ord(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         if (i <= mid) tmp(k:hi) = ord(i:mid)
         if (j <= hi) tmp(k:hi) = ord(j:hi)
         ord(lo:hi) = tmp(lo:hi)
         lo = lo + 2_I4P * width
      enddo
      width = 2_I4P * width
   enddo
   endsubroutine sort_columns_int

   pure function col_less_real(x, y) result(less)
   !< x < y lexicographically.
   real(R8P), intent(in) :: x(:), y(:)
   logical               :: less
   integer(I4P)          :: c

   less = .false.
   do c = 1, size(x)
      if (x(c) < y(c)) then
         less = .true. ; return
      elseif (x(c) > y(c)) then
         return
      endif
   enddo
   endfunction col_less_real

   pure function col_less_int(x, y) result(less)
   !< x < y lexicographically.
   integer(I4P), intent(in) :: x(:), y(:)
   logical                  :: less
   integer(I4P)             :: c

   less = .false.
   do c = 1, size(x)
      if (x(c) < y(c)) then
         less = .true. ; return
      elseif (x(c) > y(c)) then
         return
      endif
   enddo
   endfunction col_less_int

   subroutine disjoint_pair(ta, tb, touching, hit, kind, p, q)
   !< Exact narrow phase of two triangles with no vertex in common (issue #69): `pair_cut` (exact orient3d signs,
   !< implicit section points) for the general position, and for an exactly coplanar pair the relation of the two
   !< triangles in the projection that drops `facet_drop` (exact orient2d): interiors overlapping is a defect
   !< (SELFINT_COPLANAR), boundaries meeting without area is a touch. p, q: the rounded common segment (crossing), the
   !< rounded common point twice (touching), the two centroids (coplanar overlap).
   type(vector_R8P), intent(in)  :: ta(3), tb(3) !< The triangles.
   logical,          intent(in)  :: touching     !< Report pairs meeting at one point (or along a boundary, coplanar).
   logical,          intent(out) :: hit          !< True if reported.
   integer(I4P),     intent(out) :: kind         !< SELFINT_* code, if hit.
   type(vector_R8P), intent(out) :: p, q         !< See above.
   type(implicit_point_t)        :: ip, iq
   integer(I4P)                  :: pk, rel

   hit = .false. ; kind = SELFINT_CROSSING
   p = ta(1) ; q = p
   if (axis_coplanar(ta, tb)) then
      pk = PAIR_COPLANAR ! on one axis-aligned plane: exact, without pair_cut's orient3d zeros
   else
      call pair_cut(ta, tb, pk, ip, iq)
   endif
   select case (pk)
   case (PAIR_CROSSING)
      hit = .true.
      p = to_vector(ip_approx(ip)) ; q = to_vector(ip_approx(iq))
   case (PAIR_TOUCHING)
      if (.not. touching) return
      hit = .true. ; kind = SELFINT_TOUCHING
      p = to_vector(ip_approx(ip)) ; q = p
   case (PAIR_COPLANAR)
      if (.not. touching .and. .not. boxes_interiors_meet(ta, tb)) return ! cannot overlap with positive area
      rel = coplanar_relation(ta, tb)
      if (rel == COPLANAR_OVERLAP) then
         hit = .true. ; kind = SELFINT_COPLANAR
      elseif (rel == COPLANAR_TOUCH .and. touching) then
         hit = .true. ; kind = SELFINT_TOUCHING
      endif
      p = (ta(1) + ta(2) + ta(3)) / 3._R8P ; q = (tb(1) + tb(2) + tb(3)) / 3._R8P
   endselect
   endsubroutine disjoint_pair

   function coplanar_relation(ta, tb) result(rel)
   !< Relation of two exactly coplanar triangles (exact orient2d in the projection dropping `facet_drop(ta)`; the
   !< projection keeps incidences of coplanar points). Separating-axis test with the six edge lines: the interiors are
   !< disjoint iff one triangle lies in the closed outer half-plane of an edge of the other (COPLANAR_OVERLAP
   !< otherwise); the closed triangles are disjoint iff it lies in the open outer half-plane (COPLANAR_TOUCH if only
   !< the closed test separates them, COPLANAR_APART if the open one does). A degenerate triangle (zero area) is apart.
   type(vector_R8P), intent(in) :: ta(3), tb(3) !< The triangles.
   integer(I4P)                 :: rel          !< COPLANAR_APART, COPLANAR_TOUCH, COPLANAR_OVERLAP.
   real(R8P)                    :: a(2, 3), b(2, 3)
   integer(I4P)                 :: drop, sa, sb
   logical                      :: closed_sep, open_sep

   rel = COPLANAR_APART
   drop = facet_drop(ta)
   if (drop == 0_I4P) return
   call project(ta, a) ; call project(tb, b)
   sa = o2(a(:, 1), a(:, 2), a(:, 3)) ; sb = o2(b(:, 1), b(:, 2), b(:, 3))
   if (sa == 0_I4P .or. sb == 0_I4P) return
   closed_sep = .false. ; open_sep = .false.
   call separate(a, sa, b)
   call separate(b, sb, a)
   if (open_sep) then
      rel = COPLANAR_APART
   elseif (closed_sep) then
      rel = COPLANAR_TOUCH
   else
      rel = COPLANAR_OVERLAP
   endif
   contains
      subroutine project(t, x)
      !< Kept coordinates (drop + 1, drop + 2) mod 3 of t.
      type(vector_R8P), intent(in)  :: t(3)
      real(R8P),        intent(out) :: x(2, 3)
      real(R8P)                     :: c(3)
      integer(I4P)                  :: i

      do i = 1, 3
         c = [t(i)%x, t(i)%y, t(i)%z]
         x(1, i) = c(mod(drop, 3_I4P) + 1_I4P) ; x(2, i) = c(mod(drop + 1_I4P, 3_I4P) + 1_I4P)
      enddo
      endsubroutine project

      subroutine separate(t, st, o)
      !< Does an edge line of t (orientation st) separate o (closed and open sense)?
      real(R8P),    intent(in) :: t(2, 3), o(2, 3)
      integer(I4P), intent(in) :: st
      integer(I4P)             :: e, k, sg(3)

      do e = 1, 3
         do k = 1, 3
            sg(k) = st * o2(t(:, e), t(:, mod(e, 3_I4P) + 1_I4P), o(:, k)) ! > 0: inner side of the edge
         enddo
         if (all(sg <= 0_I4P)) closed_sep = .true.
         if (all(sg < 0_I4P)) open_sep = .true.
      enddo
      endsubroutine separate
   endfunction coplanar_relation

   pure function boxes_interiors_meet(ta, tb) result(yes)
   !< The interiors of the bounding boxes of two coplanar triangles meet in their plane (exact comparisons). Coplanar
   !< triangles overlap with positive area only if they do. An axis separates when the closed intervals of the two
   !< triangles at most touch; an axis on which both intervals are a single point (the plane normal to it) does not.
   type(vector_R8P), intent(in) :: ta(3), tb(3) !< The triangles.
   logical                      :: yes

   yes = .not. (apart(min(ta(1)%x, ta(2)%x, ta(3)%x), max(ta(1)%x, ta(2)%x, ta(3)%x),   &
                      min(tb(1)%x, tb(2)%x, tb(3)%x), max(tb(1)%x, tb(2)%x, tb(3)%x)) .or. &
                apart(min(ta(1)%y, ta(2)%y, ta(3)%y), max(ta(1)%y, ta(2)%y, ta(3)%y),   &
                      min(tb(1)%y, tb(2)%y, tb(3)%y), max(tb(1)%y, tb(2)%y, tb(3)%y)) .or. &
                apart(min(ta(1)%z, ta(2)%z, ta(3)%z), max(ta(1)%z, ta(2)%z, ta(3)%z),   &
                      min(tb(1)%z, tb(2)%z, tb(3)%z), max(tb(1)%z, tb(2)%z, tb(3)%z)))
   contains
      pure logical function apart(alo, ahi, blo, bhi)
      !< Intervals [alo, ahi], [blo, bhi] at most touching, unless both are one point.
      real(R8P), intent(in) :: alo, ahi, blo, bhi

      apart = .false.
      if (alo == ahi .and. blo == bhi) return
      apart = ahi <= blo .or. bhi <= alo
      endfunction apart
   endfunction boxes_interiors_meet

   pure function axis_coplanar(ta, tb) result(yes)
   !< The six vertices share one coordinate: both triangles lie on one axis-aligned plane (exact).
   type(vector_R8P), intent(in) :: ta(3), tb(3) !< The triangles.
   logical                      :: yes

   yes = (ta(2)%x == ta(1)%x .and. ta(3)%x == ta(1)%x .and. tb(1)%x == ta(1)%x .and. tb(2)%x == ta(1)%x .and. &
          tb(3)%x == ta(1)%x) .or. &
         (ta(2)%y == ta(1)%y .and. ta(3)%y == ta(1)%y .and. tb(1)%y == ta(1)%y .and. tb(2)%y == ta(1)%y .and. &
          tb(3)%y == ta(1)%y) .or. &
         (ta(2)%z == ta(1)%z .and. ta(3)%z == ta(1)%z .and. tb(1)%z == ta(1)%z .and. tb(2)%z == ta(1)%z .and. &
          tb(3)%z == ta(1)%z)
   endfunction axis_coplanar

   function o2(a, b, c) result(sg)
   !< Exact sign (-1, 0, 1) of the 2-D orientation of a, b, c (0 at once if two of them coincide).
   real(R8P), intent(in) :: a(2), b(2), c(2)
   integer(I4P)          :: sg
   real(R8P)             :: d

   sg = 0_I4P
   if (all(c == a) .or. all(c == b) .or. all(a == b)) return
   d = orient2d_xy(a(1), a(2), b(1), b(2), c(1), c(2))
   sg = 0_I4P
   if (d > 0._R8P) sg = 1_I4P
   if (d < 0._R8P) sg = -1_I4P
   endfunction o2

   pure function to_vector(x) result(v)
   !< A vector_R8P from three reals.
   real(R8P), intent(in) :: x(3)
   type(vector_R8P)      :: v

   v%x = x(1) ; v%y = x(2) ; v%z = x(3)
   endfunction to_vector

   subroutine shared_vertex_crossing(ta, tb, hit, p, q, coplanar)
   !< Exact crossing test of two triangles sharing one vertex s, neither with its two free vertices strictly on one
   !< side of the other's plane (issue #55; pair_classify's on-plane snap invented crossings on near-coplanar pairs).
   !< ta meets the plane of tb in a segment from s to a point x on the planes' common line: x is a free vertex of ta
   !< on that plane, or lies on the segment from a free vertex q off it to the other one, a. They cross (share a
   !< segment of positive length, along an edge of tb included) iff the ray s->x lies in the closed wedge of tb at s.
   !< The plane (s, b, q) cuts the plane of tb along the line s-b; x = q + t (a - q) with t > 0, so orient3d(s, b, q, x)
   !< has the sign of orient3d(s, b, q, a): every decision is the sign of orient3d on input vertices. ta in the plane
   !< of tb (both free vertices on it) is coplanar: skipped, as in the module header. p, q (s and the nearer section
   !< end, floating point) only report the segment.
   type(vector_R8P), intent(in)  :: ta(3), tb(3) !< The triangles.
   logical,          intent(out) :: hit          !< True if they cross.
   type(vector_R8P), intent(out) :: p, q         !< Common segment, if hit.
   logical,          intent(out) :: coplanar     !< The triangles are exactly coplanar (not decided here).
   type(vector_R8P)              :: s, qa, aa, b1, b2, x, y
   real(R8P)                     :: o(2), side1, side2, ref1, ref2
   integer(I4P)                  :: is, js, i, m, fa(2), fb(2)

   hit = .false. ; coplanar = .false.
   p = ta(1) ; q = p
   is = 0_I4P ; js = 0_I4P
   do i = 1, 3
      do m = 1, 3
         if (vertex_coincident(ta(i), tb(m))) then ; is = i ; js = m ; endif
      enddo
   enddo
   if (is == 0_I4P) return
   s = ta(is)
   fa = pack([1_I4P, 2_I4P, 3_I4P], [1_I4P, 2_I4P, 3_I4P] /= is)
   fb = pack([1_I4P, 2_I4P, 3_I4P], [1_I4P, 2_I4P, 3_I4P] /= js)
   b1 = tb(fb(1)) ; b2 = tb(fb(2))
   o = [orient3d(tb(1), tb(2), tb(3), ta(fa(1))), orient3d(tb(1), tb(2), tb(3), ta(fa(2)))]
   if (o(1) == 0._R8P .and. o(2) == 0._R8P) then ! coplanar: decided by the caller
      coplanar = .true. ; return
   endif
   if (o(1) /= 0._R8P) then ; qa = ta(fa(1)) ; aa = ta(fa(2)) ; else ; qa = ta(fa(2)) ; aa = ta(fa(1)) ; endif
   ! sides of x (through a) and references (the other edge's free vertex) about the planes (s, b1, qa), (s, b2, qa)
   side1 = sign3(orient3d(s, b1, qa, aa)) ; ref1 = sign3(orient3d(s, b1, qa, b2))
   side2 = sign3(orient3d(s, b2, qa, aa)) ; ref2 = sign3(orient3d(s, b2, qa, b1))
   hit = (side1 == ref1 .or. side1 == 0._R8P) .and. (side2 == ref2 .or. side2 == 0._R8P) .and. &
         .not. (side1 == 0._R8P .and. side2 == 0._R8P)
   if (.not. hit) return
   ! report: s and the nearer of the two section ends (floating point, reporting only)
   x = section_end(ta, fa, tb)
   y = section_end(tb, fb, ta)
   if (norm2([y%x - s%x, y%y - s%y, y%z - s%z]) < norm2([x%x - s%x, x%y - s%y, x%z - s%z])) x = y
   p = s ; q = x
   contains
      function section_end(t, f, o) result(e)
      !< The end of the segment where t (free vertices f) meets the plane of o: a free vertex on it, or the crossing
      !< of the edge between the free vertices.
      type(vector_R8P), intent(in) :: t(3), o(3)
      integer(I4P),     intent(in) :: f(2)
      type(vector_R8P)             :: e
      real(R8P)                    :: d1, d2

      d1 = orient3d(o(1), o(2), o(3), t(f(1))) ; d2 = orient3d(o(1), o(2), o(3), t(f(2)))
      if (d1 == 0._R8P) then
         e = t(f(1))
      elseif (d2 == 0._R8P .or. d1 * d2 > 0._R8P) then
         e = t(f(2))
      else
         e = t(f(1)) + (d1 / (d1 - d2)) * (t(f(2)) - t(f(1)))
      endif
      endfunction section_end

      pure function sign3(v) result(r)
      !< Sign of v as -1, 0, 1.
      real(R8P), intent(in) :: v
      real(R8P)             :: r

      r = 0._R8P
      if (v > 0._R8P) r = 1._R8P
      if (v < 0._R8P) r = -1._R8P
      endfunction sign3
   endsubroutine shared_vertex_crossing

   function edge_overlap_at_vertex(ta, tb) result(yes)
   !< Two exactly coplanar triangles sharing one vertex s, interiors disjoint: do they share more than s? Two wedges at
   !< s with disjoint interiors share another point only along a common edge ray: an edge of each, collinear through s
   !< (exact orient2d zero in the projection) and pointing the same way (exact coordinate comparisons).
   type(vector_R8P), intent(in) :: ta(3), tb(3) !< The triangles.
   logical                      :: yes
   real(R8P)                    :: a(3, 3), b(3, 3), s(3), ea(3), eb(3)
   integer(I4P)                 :: is, js, i, j, drop, u, v

   yes = .false.
   is = 0_I4P ; js = 0_I4P
   do i = 1, 3
      do j = 1, 3
         if (vertex_coincident(ta(i), tb(j))) then ; is = i ; js = j ; endif
      enddo
   enddo
   if (is == 0_I4P) return
   drop = facet_drop(ta)
   if (drop == 0_I4P) return
   u = mod(drop, 3_I4P) + 1_I4P ; v = mod(drop + 1_I4P, 3_I4P) + 1_I4P
   do i = 1, 3
      a(:, i) = [ta(i)%x, ta(i)%y, ta(i)%z] ; b(:, i) = [tb(i)%x, tb(i)%y, tb(i)%z]
   enddo
   s = a(:, is)
   do i = 1, 3
      if (i == is) cycle
      ea = a(:, i)
      do j = 1, 3
         if (j == js) cycle
         eb = b(:, j)
         if (o2([s(u), s(v)], [ea(u), ea(v)], [eb(u), eb(v)]) /= 0_I4P) cycle
         ! collinear through s: same direction iff, on the axis where they differ from s, both lie on the same side
         if (same_side(s(u), ea(u), eb(u)) .and. same_side(s(v), ea(v), eb(v))) then
            yes = .true. ; return
         endif
      enddo
   enddo
   contains
      pure function same_side(c, x, y) result(ok)
      !< x and y on the same side of c (both equal to c counts as the same side).
      real(R8P), intent(in) :: c, x, y
      logical               :: ok

      ok = (x > c .and. y > c) .or. (x < c .and. y < c) .or. (x == c .and. y == c)
      endfunction same_side
   endfunction edge_overlap_at_vertex

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
   !< to exact coordinate equality.
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
   !< Exact coordinate equality of two vertices (issue #69: within EPS before, so vertices 1e-8 apart counted as one
   !< and their facets were skipped).
   type(vector_R8P), intent(in) :: va, vb
   logical                      :: yes

   yes = va%x == vb%x .and. va%y == vb%y .and. va%z == vb%z
   endfunction vertex_coincident

   subroutine push_pair(buf, nfound, a, b, kind, p, q)
   !< Append a record to `buf`, doubling its capacity when full.
   type(intersection_pair_t), allocatable, intent(inout) :: buf(:)
   integer(I4P),                           intent(inout) :: nfound
   integer(I4P),                           intent(in)    :: a, b, kind
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
   buf(nfound)%kind = kind
   buf(nfound)%p = p
   buf(nfound)%q = q
   endsubroutine push_pair

endmodule fossil_self_intersection
