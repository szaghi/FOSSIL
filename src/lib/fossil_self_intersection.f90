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
integer(I4P), parameter :: COPLANAR_APART = 0_I4P, COPLANAR_TOUCH = 1_I4P, COPLANAR_OVERLAP = 2_I4P
real(R8P),    parameter :: COPLANAR_SIN2 = 1.e-24_R8P !< Facets sharing a vertex with sin^2(normal angle) below: coplanar.

type :: intersection_pair_t
   !< Single self-intersection record: an unordered pair of facet IDs (`a < b`)
   !< plus the intersection segment endpoints in world coordinates.
   integer(I4P)     :: a = 0_I4P  !< First facet id (always the smaller of the two).
   integer(I4P)     :: b = 0_I4P  !< Second facet id.
   integer(I4P)     :: kind = 1_I4P !< SELFINT_CROSSING, SELFINT_TOUCHING (include_touching), SELFINT_COPLANAR.
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

   hit = .false. ; kind = SELFINT_CROSSING
   p = fa%vertex(1) ; q = p
   select case (shared_vertices_number(fa, fb))
   case (0_I4P)
      call disjoint_pair(ta=fa%vertex, tb=fb%vertex, touching=touching, hit=hit, kind=kind, p=p, q=q)
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
      call shared_vertex_crossing(ta=fa%vertex, tb=fb%vertex, hit=hit, p=p, q=q)
   case default
      p = fa%vertex(1) ; q = p ! sharing an edge (or a duplicate): they meet along it
   endselect
   endsubroutine test_pair

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
   call pair_cut(ta, tb, pk, ip, iq)
   select case (pk)
   case (PAIR_CROSSING)
      hit = .true.
      p = to_vector(ip_approx(ip)) ; q = to_vector(ip_approx(iq))
   case (PAIR_TOUCHING)
      if (.not. touching) return
      hit = .true. ; kind = SELFINT_TOUCHING
      p = to_vector(ip_approx(ip)) ; q = p
   case (PAIR_COPLANAR)
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

   function o2(a, b, c) result(sg)
   !< Exact sign (-1, 0, 1) of the 2-D orientation of a, b, c.
   real(R8P), intent(in) :: a(2), b(2), c(2)
   integer(I4P)          :: sg
   real(R8P)             :: d

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

   subroutine shared_vertex_crossing(ta, tb, hit, p, q)
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
   type(vector_R8P)              :: s, qa, aa, b1, b2, x, y
   real(R8P)                     :: o(2), side1, side2, ref1, ref2
   integer(I4P)                  :: is, js, i, m, fa(2), fb(2)

   hit = .false.
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
   if (o(1) == 0._R8P .and. o(2) == 0._R8P) return ! coplanar
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
