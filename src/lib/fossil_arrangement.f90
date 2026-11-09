!< FOSSIL, mesh-arrangement scaffolding for booleans (issue #18 §1.1 stage 1).

module fossil_arrangement
!< FOSSIL, mesh-arrangement scaffolding for booleans (issue #18 §1.1 stage 1).
!<
!< The arrangement of two surfaces A and B for the booleans (Zhou, Grinspun, Zorin, Jacobson, *Mesh Arrangements for
!< Solid Geometry*, SIGGRAPH 2016), exact since issue #65 (Cherchi, Livesu, Scateni, Attene 2020, *Fast and robust mesh
!< arrangements using floating-point arithmetic*): no tolerance anywhere.
!<
!< 1. `arrangement_initialize`: A's and B's facets, owners, empty cut lists.
!< 2. `arrangement_collect_intersections`: every pair of facets whose boxes overlap (tree-vs-tree) is classified
!<    exactly (`pair_cut`, `coplanar_cut` of `fossil_tritri_exact`); a crossing adds a segment between two implicit
!<    points (`fossil_implicit_points`) to both facets; an exactly coplanar pair records the partner and cuts each
!<    facet along the portions of the other's edges inside it. The registry then merges the points: one index per
!<    geometric point, exactly (identity, then exact coordinates), the most specific construction kept, rounded once;
!<    distinct points get distinct rounded coordinates (`separate_rounded`). Two cuts of one facet made by two
!<    transversal facets that cross inside it (exactly, `ip_orient2d`) meet at the point of the three planes
!<    (IPK_TPI, issue #91): it joins the three facets as an isolated point, and the triangulator splits both cuts there.
!<    Points inside a facet edge go to the facet across it (T-junctions, exact collinearity).
!< 3. `arrangement_retriangulate`: each cut facet is triangulated by `cdt_facet` (exact orient2d on implicit points),
!<    with its cuts as constraints; sub-triangles keep the parent's orientation and carry their registry vertices and
!<    which of their edges lie on a cut (for the selection by patches in `fossil_boolean`).
!<
!< Exactly degenerate facets (zero area, `facet_drop` = 0; issues #98, #99): `pair_cut` needs two non-degenerate
!< triangles, so a pair with a degenerate facet is decided by `degenerate_pair_meet` (exact) instead. No contact: the
!< pair is skipped (no cut, no partner) and the facet passes through the arrangement uncut, as in the input. Contact
!< (touching or crossing): the arrangement is flagged (`degenerate_contact`) and not retriangulated
!< (ARR_STATUS_DEGENERATE): a zero-area facet cannot be cut, and its neighbours' cuts along it would leave a
!< T-junction. The boolean reports it as BOOL_STATUS_INVALID_INPUT (`sanitize` the inputs first).
!<
!< Classification, selection and stitching are in `fossil_boolean`.

use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_cdt_implicit,     only : cdt_facet, CDTI_OK
use fossil_facet_object,     only : facet_object
use fossil_implicit_points,  only : implicit_point_t, ip_explicit, ip_approx, ip_orient2d, ip_lex_compare, &
                                    ip_identity_compare, ip_rank, ip_compare, ip_tpi
use fossil_list_id_object,   only : list_id_object
use fossil_tritri_exact,     only : pair_cut, coplanar_cut, facet_drop, degenerate_pair_meet, PAIR_CROSSING, &
                                    PAIR_COPLANAR
use fossil_utils,            only : EPS
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: cut_list_t
public :: arrangement_t
public :: arrangement_initialize
public :: arrangement_collect_intersections
public :: arrangement_collect_pairs
public :: arrangement_retriangulate
public :: coplanar_partner
public :: coplanar_cover
public :: cut_sources
public :: project_to_plane
public :: lift_from_plane
public :: ARR_STATUS_OK, ARR_STATUS_CDT_FAILED, ARR_STATUS_DEGENERATE

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P  !< Buffer size for `enumerate_children` (octree fan-out).

! Status codes for arrangement_retriangulate.
integer(I4P), parameter :: ARR_STATUS_OK         = 0_I4P  !< Retriangulation completed successfully.
integer(I4P), parameter :: ARR_STATUS_CDT_FAILED = 1_I4P  !< CDT could not recover one or more constraint segments.
integer(I4P), parameter :: ARR_STATUS_DEGENERATE = 2_I4P  !< A zero-area facet meets a facet it is paired with (#99).

type :: cut_list_t
   !< Per-facet cuts (issue #65): segments between implicit points of the arrangement's registry (`arrangement_t%pts`),
   !< by index; after `arrangement_collect_intersections` the indices are those of the merged registry (one index per
   !< geometric point, exactly).
   integer(I4P)              :: n_segments = 0_I4P  !< Number of intersection segments stored.
   integer(I4P), allocatable :: seg(:, :)           !< Segment endpoints (2, capacity), registry indices.
   integer(I4P), allocatable :: src(:)              !< Facet of the other surface that made each segment.
   integer(I4P)              :: n_steiner = 0_I4P   !< Number of isolated points (issue #33).
   integer(I4P), allocatable :: steiner(:)          !< Points on the facet's edges coming from the neighbour across them.
   integer(I4P)              :: n_coplanar = 0_I4P  !< Number of coplanar partners (issue #63).
   integer(I4P), allocatable :: coplanar(:)         !< Global indices of the facets of the other surface exactly coplanar
                                                    !< with this one (PAIR_COPLANAR).
endtype cut_list_t

type :: arrangement_t
   !< In-progress mesh arrangement of two surfaces A and B.
   !<
   !< Layout: `facet(1:nA)` are A's facets in their original order (owner = 1),
   !< `facet(nA+1:nA+nB)` are B's facets (owner = 2). The `cut(:)` array is
   !< parallel to `facet(:)`. The `tree_a` / `tree_b` references hold the AABB
   !< trees built over A and B respectively — used by `arrangement_collect_intersections`
   !< for the cross-mesh broad phase. They are **not owned** by this module:
   !< we keep references for traversal but expect the caller (the boolean
   !< driver in surface_stl_object) to manage their lifetime.
   integer(I4P)                       :: n_a = 0_I4P    !< Number of facets from surface A.
   integer(I4P)                       :: n_b = 0_I4P    !< Number of facets from surface B.
   type(facet_object), allocatable    :: facet(:)       !< Concatenated facets [A's then B's].
   integer(I4P),       allocatable    :: owner(:)       !< 1 = facet originated in A, 2 = in B.
   type(cut_list_t),   allocatable    :: cut(:)         !< Per-facet intersection-segment list.
   type(aabb_tree_object), pointer    :: tree_a => null() !< Reference to A's AABB tree (not owned).
   type(aabb_tree_object), pointer    :: tree_b => null() !< Reference to B's AABB tree (not owned).
   integer(I4P)                       :: n_pts = 0_I4P  !< Points in the registry.
   type(implicit_point_t), allocatable :: pts(:)        !< Registry of the cut points (implicit, issue #65).
   real(R8P),          allocatable    :: xyz(:, :)      !< Rounded coordinates of the registry points (3, n_pts).
   integer(I4P),       allocatable    :: corner(:, :)   !< Registry indices of the vertices of every facet (3, n_a + n_b).
   integer(I4P),       allocatable    :: drop(:)        !< `facet_drop` of every facet, 0 if exactly degenerate (#99).
   logical                            :: degenerate_contact = .false. !< A degenerate facet meets its pair (#99).
endtype arrangement_t

contains

   subroutine arrangement_initialize(arr, facet_a, tree_a, facet_b, tree_b)
   !< Populate `arr` from the two input surfaces' facet arrays + AABB trees.
   !<
   !< Copies both facet arrays into `arr%facet` (so subsequent retriangulation
   !< edits don't mutate the inputs), tags each facet's owner, allocates empty
   !< per-facet cut lists. Stores references to the input trees for the
   !< broad-phase step.
   !<
   !< @note Facets are deep-copied (via Fortran's intrinsic derived-type
   !<       assignment, which deep-copies the allocatable components).
   !<       For 100k-facet meshes this is ~10 MB of copy; acceptable for the
   !<       MVP. A future zero-copy variant would take pointers into the
   !<       caller's facet arrays.
   type(arrangement_t),    intent(out)        :: arr      !< Arrangement to populate.
   type(facet_object),     intent(in), target :: facet_a(:) !< Facets of surface A.
   type(aabb_tree_object), intent(in), target :: tree_a   !< AABB tree of surface A.
   type(facet_object),     intent(in), target :: facet_b(:) !< Facets of surface B.
   type(aabb_tree_object), intent(in), target :: tree_b   !< AABB tree of surface B.
   integer(I4P)                               :: n_total, i

   arr%n_a = size(facet_a, dim=1)
   arr%n_b = size(facet_b, dim=1)
   n_total = arr%n_a + arr%n_b

   allocate(arr%facet(n_total))
   allocate(arr%owner(n_total))
   allocate(arr%cut(n_total))
   allocate(arr%drop(n_total))

   do i = 1, arr%n_a
      arr%facet(i) = facet_a(i)
      arr%owner(i) = 1_I4P
   enddo
   do i = 1, arr%n_b
      arr%facet(arr%n_a + i) = facet_b(i)
      arr%owner(arr%n_a + i) = 2_I4P
   enddo
   do i = 1, n_total
      arr%drop(i) = facet_drop(arr%facet(i)%vertex)
   enddo

   arr%tree_a => tree_a
   arr%tree_b => tree_b
   endsubroutine arrangement_initialize

   subroutine arrangement_collect_intersections(arr)
   !< For every (A-facet, B-facet) pair whose AABBs overlap, run the Möller
   !< tri-tri test and append the intersection segment to **both** facets'
   !< cut lists.
   !<
   !< Cross-mesh broad phase: tree-vs-tree traversal between A's tree and B's
   !< tree. Same algorithmic shape as §1.2's self-intersection broad phase but
   !< with the trees coming from different surfaces, so:
   !<   - the "i < j" canonicalisation does not apply (A-facets and B-facets
   !<     are in disjoint id spaces and every pair is asymmetric)
   !<   - the adjacency filter does not apply (vertices of A and B are
   !<     in different coordinate frames; even if they coincide spatially they
   !<     are different topological vertices and intersections are real)
   !<
   !< Both endpoints of each intersection segment are pushed; if a future stage
   !< needs to deduplicate near-coincident segments shared by adjacent facets,
   !< it has all the data it needs in the cut lists.
   type(arrangement_t), intent(inout) :: arr

   if (arr%n_a == 0 .or. arr%n_b == 0) return
   if (.not. associated(arr%tree_a) .or. .not. associated(arr%tree_b)) return
   if (arr%tree_a%get_nodes_number() <= 0 .or. arr%tree_b%get_nodes_number() <= 0) then
      ! Brute-force fallback when either tree is empty.
      call collect_brute_force(arr=arr)
      call merge_points(arr=arr)
      call add_triple_points(arr=arr)
      call propagate_edge_points(arr=arr)
      return
   endif

   call traverse_cross_pair(arr=arr, na=0_I4P, nb=0_I4P)
   call merge_points(arr=arr)
   call add_triple_points(arr=arr)
   call propagate_edge_points(arr=arr)
   endsubroutine arrangement_collect_intersections

   subroutine arrangement_collect_pairs(arr, pa, pb)
   !< The arrangement of one surface with itself (issue #78): `arr` holds the surface as A (B empty) and the pairs
   !< (pa(k), pb(k)) are its facet pairs that cross or overlap in a common plane (`find_self_intersections`). Each
   !< pair is cut exactly as a pair of A and B (`collect_pair`), then the registry is merged and the T-junctions closed
   !< as in `arrangement_collect_intersections`. Facets meeting only along a shared edge or at a shared vertex add no
   !< cut.
   type(arrangement_t), intent(inout) :: arr    !< Arrangement (A = the surface, no B).
   integer(I4P),        intent(in)    :: pa(:)  !< First facet of each pair.
   integer(I4P),        intent(in)    :: pb(:)  !< Second facet of each pair.
   integer(I4P)                       :: k

   do k = 1, size(pa)
      call collect_pair(arr=arr, gi=pa(k), gj=pb(k))
   enddo
   call merge_points(arr=arr)
   call add_triple_points(arr=arr)
   call propagate_edge_points(arr=arr)
   endsubroutine arrangement_collect_pairs

   subroutine collect_brute_force(arr)
   !< O(nA * nB) fallback used only when trees are not initialized.
   type(arrangement_t), intent(inout) :: arr
   integer(I4P)                       :: i, j, gi, gj

   do i = 1, arr%n_a
      gi = i
      do j = 1, arr%n_b
         gj = arr%n_a + j
         if (.not. bboxes_overlap(arr%facet(gi), arr%facet(gj))) cycle
         call collect_pair(arr=arr, gi=gi, gj=gj)
      enddo
   enddo
   endsubroutine collect_brute_force

   recursive subroutine traverse_cross_pair(arr, na, nb)
   !< Cross-mesh tree-vs-tree traversal. `na` is a node index in `arr%tree_a`,
   !< `nb` is a node index in `arr%tree_b`. Prunes when AABBs disjoint, descends
   !< the larger box otherwise; at leaf-leaf, enumerates candidate facet pairs.
   type(arrangement_t), intent(inout) :: arr
   integer(I4P),        intent(in)    :: na, nb
   type(aabb_node_object), pointer    :: node_a, node_b
   integer(I4P)                       :: ca(MAX_CHILDREN), cb(MAX_CHILDREN), nca, ncb
   integer(I4P)                       :: i, j
   type(vector_R8P)                   :: amin, amax, bmin, bmax
   real(R8P)                          :: va, vb

   node_a => arr%tree_a%node_at(i=na)
   node_b => arr%tree_b%node_at(i=nb)
   if (.not. associated(node_a) .or. .not. associated(node_b)) return
   if (.not. node_a%is_allocated() .or. .not. node_b%is_allocated()) return

   amin = node_a%bmin() ; amax = node_a%bmax()
   bmin = node_b%bmin() ; bmax = node_b%bmax()
   if (amax%x < bmin%x .or. bmax%x < amin%x .or. &
       amax%y < bmin%y .or. bmax%y < amin%y .or. &
       amax%z < bmin%z .or. bmax%z < amin%z) return

   call arr%tree_a%enumerate_children(n=na, out_idx=ca, nchild=nca)
   call arr%tree_b%enumerate_children(n=nb, out_idx=cb, nchild=ncb)

   if (nca == 0 .and. ncb == 0) then
      call enumerate_leaf_cross_pair(arr=arr, na=na, nb=nb)
      return
   endif

   if (nca == 0) then
      do j = 1, ncb
         call traverse_cross_pair(arr=arr, na=na, nb=cb(j))
      enddo
      return
   endif
   if (ncb == 0) then
      do i = 1, nca
         call traverse_cross_pair(arr=arr, na=ca(i), nb=nb)
      enddo
      return
   endif

   ! Descend the larger box to keep the recursion balanced.
   va = (amax%x - amin%x) * (amax%y - amin%y) * (amax%z - amin%z)
   vb = (bmax%x - bmin%x) * (bmax%y - bmin%y) * (bmax%z - bmin%z)
   if (va >= vb) then
      do i = 1, nca
         call traverse_cross_pair(arr=arr, na=ca(i), nb=nb)
      enddo
   else
      do j = 1, ncb
         call traverse_cross_pair(arr=arr, na=na, nb=cb(j))
      enddo
   endif
   endsubroutine traverse_cross_pair

   subroutine enumerate_leaf_cross_pair(arr, na, nb)
   !< Enumerate facet pairs at the leaf-pair (na in tree_a, nb in tree_b).
   !< Each pair's facets are appended into the arrangement as global indices:
   !< A's facet i is at global index i; B's facet j is at global index n_a + j.
   type(arrangement_t), intent(inout) :: arr
   integer(I4P),        intent(in)    :: na, nb
   type(aabb_node_object), pointer    :: node_a, node_b
   type(list_id_object)               :: ids_a, ids_b
   integer(I4P)                       :: i, j, fa, fb, gi, gj

   node_a => arr%tree_a%node_at(i=na)
   node_b => arr%tree_b%node_at(i=nb)
   if (.not. node_a%has_facets() .or. .not. node_b%has_facets()) return

   ids_a = node_a%facet_id()
   ids_b = node_b%facet_id()

   do i = 1, ids_a%ids_number
      fa = ids_a%id(i)
      if (fa < 1 .or. fa > arr%n_a) cycle
      gi = fa
      do j = 1, ids_b%ids_number
         fb = ids_b%id(j)
         if (fb < 1 .or. fb > arr%n_b) cycle
         gj = arr%n_a + fb
         if (.not. bboxes_overlap(arr%facet(gi), arr%facet(gj))) cycle
         call collect_pair(arr=arr, gi=gi, gj=gj)
      enddo
   enddo
   endsubroutine enumerate_leaf_cross_pair

   subroutine register_point(arr, x, id)
   !< Append an implicit point to the registry (no lookup: `merge_points` merges afterwards).
   type(arrangement_t),    intent(inout) :: arr !< Arrangement.
   type(implicit_point_t), intent(in)    :: x   !< Point.
   integer(I4P),           intent(out)   :: id  !< Its raw index.
   type(implicit_point_t), allocatable   :: tmp(:)

   if (.not. allocated(arr%pts)) allocate(arr%pts(64))
   if (arr%n_pts == size(arr%pts)) then
      allocate(tmp(2 * size(arr%pts))) ; tmp(1:arr%n_pts) = arr%pts(1:arr%n_pts) ; call move_alloc(tmp, arr%pts)
   endif
   arr%n_pts = arr%n_pts + 1_I4P
   arr%pts(arr%n_pts) = x
   id = arr%n_pts
   endsubroutine register_point

   subroutine merge_points(arr)
   !< One registry index per geometric point, exactly (issue #65). The vertices of every facet are registered too (the
   !< selection connects sub-triangles by registry indices). Points of one identity (`ip_same`) merge first (sort by identity), then points of
   !< different identities and the same exact coordinates (sort by exact lexicographic order, `ip_lex_compare`):
   !< for clean inputs the classification already gives one point one identity, this catches the rest. A merged
   !< point keeps its most specific construction (`ip_rank`: an input vertex first), the same for every facet, and
   !< its coordinates are rounded once (`ip_approx`).
   type(arrangement_t), intent(inout)  :: arr
   integer(I4P)                        :: n_total, gi, k

   n_total = arr%n_a + arr%n_b
   allocate(arr%corner(3, n_total)) ; arr%corner = 0_I4P
   do gi = 1, n_total
      do k = 1, 3
         call register_point(arr, ip_explicit(xyz_of(arr%facet(gi)%vertex(k))), arr%corner(k, gi))
      enddo
   enddo
   call merge_registry(arr)
   endsubroutine merge_points

   subroutine merge_registry(arr)
   !< The merge of `merge_points` on the registry as it is (raw points appended after an earlier merge included): one
   !< index per geometric point, the segments, corners and isolated points remapped, coordinates rounded once.
   type(arrangement_t), intent(inout)  :: arr
   integer(I4P), allocatable           :: idx(:), grp(:), ulist(:), gidx(:), fid(:), best(:), map(:)
   type(implicit_point_t), allocatable :: merged(:)
   integer(I4P)                        :: n_total, gi, k, n, u, nf, j

   n_total = arr%n_a + arr%n_b
   n = arr%n_pts
   if (allocated(arr%xyz)) deallocate(arr%xyz)
   allocate(arr%xyz(3, 0))
   if (n == 0_I4P) return
   allocate(idx(n), grp(n), ulist(n), map(n))
   idx = [(k, k = 1, n)]
   call sort_points(arr%pts, idx, 1_I4P)
   u = 0_I4P
   do k = 1, n
      if (k == 1) then
         u = 1_I4P ; ulist(1) = idx(1)
      elseif (ip_identity_compare(arr%pts(idx(k - 1)), arr%pts(idx(k))) /= 0_I4P) then
         u = u + 1_I4P ; ulist(u) = idx(k)
      endif
      grp(idx(k)) = u
   enddo
   allocate(gidx(u), fid(u), best(u))
   gidx = ulist(1:u)
   call sort_points(arr%pts, gidx, 2_I4P)
   nf = 0_I4P
   do k = 1, u
      if (k == 1) then
         nf = 1_I4P ; best(1) = gidx(1)
      elseif (ip_lex_compare(arr%pts(gidx(k - 1)), arr%pts(gidx(k))) /= 0_I4P) then
         nf = nf + 1_I4P ; best(nf) = gidx(k)
      elseif (ip_rank(arr%pts(gidx(k))) < ip_rank(arr%pts(best(nf)))) then
         best(nf) = gidx(k)
      endif
      fid(grp(gidx(k))) = nf
   enddo
   do k = 1, n
      map(k) = fid(grp(k))
   enddo
   allocate(merged(nf))
   do k = 1, nf
      merged(k) = arr%pts(best(k))
   enddo
   call move_alloc(merged, arr%pts)
   arr%n_pts = nf
   deallocate(arr%xyz) ; allocate(arr%xyz(3, nf))
   do k = 1, nf
      arr%xyz(:, k) = ip_approx(arr%pts(k))
   enddo
   call separate_rounded(arr)
   do gi = 1, n_total
      do j = 1, arr%cut(gi)%n_segments
         arr%cut(gi)%seg(:, j) = map(arr%cut(gi)%seg(:, j))
      enddo
      do k = 1, 3
         if (arr%corner(k, gi) > 0_I4P) arr%corner(k, gi) = map(arr%corner(k, gi))
      enddo
      do j = 1, arr%cut(gi)%n_steiner
         arr%cut(gi)%steiner(j) = map(arr%cut(gi)%steiner(j))
      enddo
   enddo
   endsubroutine merge_registry

   subroutine add_triple_points(arr)
   !< Cuts of one facet that cross (issue #91). A facet gi cut by two facets g, h transversal to it (not coplanar
   !< partners) has two cuts on the lines gi ^ g and gi ^ h; when they cross at a point inside both (exact: the ends of
   !< each strictly on either side of the other in gi's projection, `ip_orient2d`) that point is the meeting point of
   !< the three planes, which no pair of facets constructs: an IPK_TPI. It joins gi, g and h as an isolated point
   !< (the triangulator splits every cut through it) and the registry is merged again, so the three facets share it,
   !< and a TPI equal to a point of the registry (an LPI where the crossing is on a facet edge) takes that point. No
   !< crossing: nothing changes.
   type(arrangement_t), intent(inout) :: arr
   integer(I4P)                       :: gi, k, l, drop, a, b, c, d, g, h, id, n_added

   n_added = 0_I4P
   do gi = 1, arr%n_a + arr%n_b
      if (arr%cut(gi)%n_segments < 2_I4P) cycle
      drop = facet_drop(arr%facet(gi)%vertex)
      if (drop == 0_I4P) cycle
      do k = 1, arr%cut(gi)%n_segments - 1
         g = arr%cut(gi)%src(k)
         if (is_partner(arr, gi, g)) cycle
         a = arr%cut(gi)%seg(1, k) ; b = arr%cut(gi)%seg(2, k)
         if (a == b) cycle
         do l = k + 1, arr%cut(gi)%n_segments
            h = arr%cut(gi)%src(l)
            if (h == g) cycle
            if (is_partner(arr, gi, h)) cycle
            c = arr%cut(gi)%seg(1, l) ; d = arr%cut(gi)%seg(2, l)
            if (c == d .or. c == a .or. c == b .or. d == a .or. d == b) cycle
            if (ip_orient2d(arr%pts(a), arr%pts(b), arr%pts(c), drop) * &
                ip_orient2d(arr%pts(a), arr%pts(b), arr%pts(d), drop) >= 0_I4P) cycle
            if (ip_orient2d(arr%pts(c), arr%pts(d), arr%pts(a), drop) * &
                ip_orient2d(arr%pts(c), arr%pts(d), arr%pts(b), drop) >= 0_I4P) cycle
            call register_point(arr, ip_tpi(plane_of(arr, gi), plane_of(arr, g), plane_of(arr, h)), id)
            call append_steiner(cut=arr%cut(gi), id=id)
            call append_steiner(cut=arr%cut(g), id=id)
            call append_steiner(cut=arr%cut(h), id=id)
            n_added = n_added + 1_I4P
         enddo
      enddo
   enddo
   if (n_added > 0_I4P) call merge_registry(arr)
   endsubroutine add_triple_points

   pure function is_partner(arr, f, o) result(yes)
   !< Facet o is a coplanar partner of facet f (its cuts in f lie on o's edges, not on a plane through f).
   type(arrangement_t), intent(in) :: arr  !< Arrangement.
   integer(I4P),        intent(in) :: f, o !< Facets.
   logical                         :: yes  !< Partner.

   yes = .false.
   if (arr%cut(f)%n_coplanar > 0_I4P) yes = any(arr%cut(f)%coplanar(1:arr%cut(f)%n_coplanar) == o)
   endfunction is_partner

   pure function plane_of(arr, f) result(w)
   !< The three vertices of facet f (columns).
   type(arrangement_t), intent(in) :: arr     !< Arrangement.
   integer(I4P),        intent(in) :: f       !< Facet.
   real(R8P)                       :: w(3, 3) !< Vertices.
   integer(I4P)                    :: v

   do v = 1, 3
      w(:, v) = xyz_of(arr%facet(f)%vertex(v))
   enddo
   endfunction plane_of

   subroutine separate_rounded(arr)
   !< Distinct registry points get distinct rounded coordinates (issue #65). Two exact points closer than an ulp can
   !< round to the same doubles (the cut points of faces coplanar by construction but off by round-off); the welded
   !< result would join their edges (non-manifold, or a crack). Of the points sharing one rounded triple, the most
   !< specific construction keeps it (an input vertex never moves) and each other one moves its x coordinate by as
   !< many ulps as its rank in the group, until every triple is unique: the topology of the output is the exact one,
   !< the geometry moves by a few ulps (the order of the rounding itself).
   type(arrangement_t), intent(inout) :: arr
   integer(I4P), allocatable          :: ord(:)
   integer(I4P)                       :: n, i, j, k, m, pass, t
   logical                            :: clash

   n = arr%n_pts
   if (n < 2_I4P) return
   allocate(ord(n))
   do pass = 1, 64
      ord = [(k, k = 1, n)]
      call sort_xyz(arr%xyz, ord)
      clash = .false.
      i = 1_I4P
      do while (i <= n)
         j = i
         do while (j < n)
            if (any(arr%xyz(:, ord(j + 1)) /= arr%xyz(:, ord(i)))) exit
            j = j + 1_I4P
         enddo
         if (j > i) then
            clash = .true.
            do k = i + 1, j ! most specific first (stable: then by index)
               t = ord(k) ; m = k - 1_I4P
               do while (m >= i)
                  if (ip_rank(arr%pts(ord(m))) <= ip_rank(arr%pts(t))) exit
                  ord(m + 1) = ord(m) ; m = m - 1_I4P
               enddo
               ord(m + 1) = t
            enddo
            do k = i + 1, j
               do m = 1, k - i
                  arr%xyz(1, ord(k)) = nearest(arr%xyz(1, ord(k)), 1._R8P)
               enddo
            enddo
         endif
         i = j + 1_I4P
      enddo
      if (.not. clash) return
   enddo
   endsubroutine separate_rounded

   subroutine sort_xyz(xyz, key)
   !< Bottom-up merge sort of the indices `key` by the lexicographic order of the triples xyz(:, key).
   real(R8P),    intent(in)    :: xyz(:, :)
   integer(I4P), intent(inout) :: key(:)
   integer(I4P), allocatable   :: work(:)
   integer(I4P)                :: n, width, lo, mid, hi, i, j, k

   n = size(key)
   allocate(work(n))
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (less(key(j), key(i))) then
               work(k) = key(j) ; j = j + 1_I4P
            else
               work(k) = key(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i <= mid)
            work(k) = key(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j <= hi)
            work(k) = key(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         lo = lo + 2_I4P * width
      enddo
      key = work
      width = 2_I4P * width
   enddo
   contains
      pure function less(a, b) result(l)
      integer(I4P), intent(in) :: a, b
      logical                  :: l
      integer(I4P)             :: c

      l = .false.
      do c = 1, 3
         if (xyz(c, a) /= xyz(c, b)) then
            l = xyz(c, a) < xyz(c, b) ; return
         endif
      enddo
      endfunction less
   endsubroutine sort_xyz

   subroutine sort_points(pts, key, mode)
   !< Stable bottom-up merge sort of the indices `key` into `pts`: mode 1 by identity, mode 2 by exact coordinates.
   type(implicit_point_t), intent(in)    :: pts(:)
   integer(I4P),           intent(inout) :: key(:)
   integer(I4P),           intent(in)    :: mode
   integer(I4P), allocatable             :: work(:)
   integer(I4P)                          :: n, width, lo, mid, hi, i, j, k

   n = size(key)
   allocate(work(n))
   width = 1_I4P
   do while (width < n)
      lo = 1_I4P
      do while (lo <= n)
         mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
         i = lo ; j = mid + 1_I4P ; k = lo
         do while (i <= mid .and. j <= hi)
            if (cmp(key(j), key(i)) < 0_I4P) then
               work(k) = key(j) ; j = j + 1_I4P
            else
               work(k) = key(i) ; i = i + 1_I4P
            endif
            k = k + 1_I4P
         enddo
         do while (i <= mid)
            work(k) = key(i) ; i = i + 1_I4P ; k = k + 1_I4P
         enddo
         do while (j <= hi)
            work(k) = key(j) ; j = j + 1_I4P ; k = k + 1_I4P
         enddo
         lo = lo + 2_I4P * width
      enddo
      key = work
      width = 2_I4P * width
   enddo
   contains
      function cmp(a, b) result(c)
      integer(I4P), intent(in) :: a, b
      integer(I4P)             :: c

      if (mode == 1_I4P) then
         c = ip_identity_compare(pts(a), pts(b))
      else
         c = ip_lex_compare(pts(a), pts(b))
      endif
      endfunction cmp
   endsubroutine sort_points

   subroutine propagate_edge_points(arr)
   !< Close the T-junctions (issues #33, #65). A cut point lying inside an edge of its facet splits that edge in this
   !< facet's triangulation; the facet across the edge (same surface, `fcon_edge`) gets it as an isolated point, so
   !< both sides split the edge at the same registry points. "Inside the edge" is exact: the point is in the closed
   !< facet, is not one of its vertices (registry index), and is collinear with the edge (`ip_orient2d` = 0).
   type(arrangement_t), intent(inout) :: arr      !< Arrangement.
   integer(I4P), allocatable          :: ids(:)
   integer(I4P)                       :: gi, k, e, nb, gn, n, drop
   type(implicit_point_t)             :: ce(3)

   do gi = 1, arr%n_a + arr%n_b
      if (arr%cut(gi)%n_segments == 0_I4P) cycle
      call unique_ids(arr%cut(gi), ids, n, with_steiner=.false.)
      drop = facet_drop(arr%facet(gi)%vertex)
      if (drop == 0_I4P) cycle
      do k = 1, 3
         ce(k) = arr%pts(arr%corner(k, gi))
      enddo
      do k = 1, n
         if (any(arr%corner(:, gi) == ids(k))) cycle
         do e = 1, 3
            if (ip_orient2d(ce(e), ce(mod(e, 3) + 1), arr%pts(ids(k)), drop) /= 0_I4P) cycle
            nb = arr%facet(gi)%fcon_edge(e)
            if (nb > 0_I4P) then
               gn = nb ; if (arr%owner(gi) == 2_I4P) gn = arr%n_a + nb
               if (gn >= 1_I4P .and. gn <= arr%n_a + arr%n_b) call append_steiner(cut=arr%cut(gn), id=ids(k))
            endif
            exit
         enddo
      enddo
   enddo
   endsubroutine propagate_edge_points

   subroutine unique_ids(cut, ids, n, with_steiner)
   !< Sorted distinct registry indices of a facet's segment endpoints (and isolated points).
   type(cut_list_t),          intent(in)  :: cut
   integer(I4P), allocatable, intent(out) :: ids(:)
   integer(I4P),              intent(out) :: n
   logical,                   intent(in)  :: with_steiner
   integer(I4P), allocatable              :: all_ids(:)
   integer(I4P)                           :: m, k, i, t

   m = 2_I4P * cut%n_segments
   if (with_steiner) m = m + cut%n_steiner
   allocate(all_ids(m))
   if (cut%n_segments > 0_I4P) all_ids(1:2 * cut%n_segments) = reshape(cut%seg(:, 1:cut%n_segments), [2 * cut%n_segments])
   if (with_steiner .and. cut%n_steiner > 0_I4P) all_ids(2 * cut%n_segments + 1:m) = cut%steiner(1:cut%n_steiner)
   do k = 2, m ! insertion sort (short lists)
      t = all_ids(k) ; i = k - 1_I4P
      do while (i >= 1)
         if (all_ids(i) <= t) exit
         all_ids(i + 1) = all_ids(i) ; i = i - 1_I4P
      enddo
      all_ids(i + 1) = t
   enddo
   allocate(ids(m))
   n = 0_I4P
   do k = 1, m
      if (n > 0) then
         if (ids(n) == all_ids(k)) cycle
      endif
      n = n + 1_I4P ; ids(n) = all_ids(k)
   enddo
   endsubroutine unique_ids

   subroutine append_steiner(cut, id)
   !< Append an isolated point (registry index) to a cut list, no duplicates.
   type(cut_list_t), intent(inout) :: cut !< Cut list.
   integer(I4P),     intent(in)    :: id  !< Registry index.
   integer(I4P), allocatable       :: tmp(:)

   if (cut%n_steiner > 0) then
      if (any(cut%steiner(1:cut%n_steiner) == id)) return
   endif
   if (.not. allocated(cut%steiner)) then
      allocate(cut%steiner(4))
   elseif (cut%n_steiner == size(cut%steiner)) then
      allocate(tmp(2 * size(cut%steiner)))
      tmp(1:cut%n_steiner) = cut%steiner(1:cut%n_steiner)
      call move_alloc(from=tmp, to=cut%steiner)
   endif
   cut%n_steiner = cut%n_steiner + 1_I4P
   cut%steiner(cut%n_steiner) = id
   endsubroutine append_steiner

   subroutine append_partner(cut, g)
   !< Append a coplanar partner to a cut list (no duplicates).
   type(cut_list_t), intent(inout) :: cut !< Cut list.
   integer(I4P),     intent(in)    :: g   !< Global index of the partner facet.
   integer(I4P), allocatable       :: tmp(:)

   if (cut%n_coplanar > 0) then
      if (any(cut%coplanar(1:cut%n_coplanar) == g)) return
   endif
   if (.not. allocated(cut%coplanar)) then
      allocate(cut%coplanar(4))
   elseif (cut%n_coplanar == size(cut%coplanar)) then
      allocate(tmp(2 * size(cut%coplanar)))
      tmp(1:cut%n_coplanar) = cut%coplanar(1:cut%n_coplanar)
      call move_alloc(from=tmp, to=cut%coplanar)
   endif
   cut%n_coplanar = cut%n_coplanar + 1_I4P
   cut%coplanar(cut%n_coplanar) = g
   endsubroutine append_partner

   function coplanar_partner(arr, gi, pid) result(partner)
   !< The coplanar partner of facet gi (a facet of the other surface exactly coplanar with it) containing the
   !< sub-triangle of gi with registry vertices `pid`, or 0 (issues #63, #65): the first of `coplanar_cover`.
   type(arrangement_t), intent(in) :: arr     !< Arrangement.
   integer(I4P),        intent(in) :: gi      !< Parent facet of the sub-triangle.
   integer(I4P),        intent(in) :: pid(3)  !< Registry indices of the sub-triangle's vertices.
   integer(I4P)                    :: partner !< Partner index, 0 if none.
   integer(I4P), allocatable       :: g(:)
   integer(I4P)                    :: ng

   call coplanar_cover(arr, gi, pid, g, ng, first=.true.)
   partner = 0_I4P
   if (ng > 0_I4P) partner = g(1)
   endfunction coplanar_partner

   subroutine coplanar_cover(arr, gi, pid, g, ng, first)
   !< The coplanar partners of facet gi (facets exactly coplanar with it, recorded by `collect_pair`) containing the
   !< sub-triangle of gi with registry vertices `pid`, in the order they were recorded (issues #63, #65, #78). The
   !< arrangement cut gi along the portions of its partners' edges inside it, so a sub-triangle lies inside a partner
   !< or outside it: inside iff its three vertices are in the closed partner (a convex set), by exact `ip_orient2d` in
   !< the projection of gi.
   type(arrangement_t),       intent(in)           :: arr    !< Arrangement.
   integer(I4P),              intent(in)           :: gi     !< Parent facet of the sub-triangle.
   integer(I4P),              intent(in)           :: pid(3) !< Registry indices of the sub-triangle's vertices.
   integer(I4P), allocatable, intent(out)          :: g(:)   !< Partners containing it.
   integer(I4P),              intent(out)          :: ng     !< Their number.
   logical,                   intent(in), optional :: first  !< Stop at the first one (default .false.).
   integer(I4P)                                    :: k, drop, og, j, e
   type(implicit_point_t)                          :: gv(3)
   logical                                         :: inside

   allocate(g(max(1_I4P, arr%cut(gi)%n_coplanar))) ; ng = 0_I4P
   if (arr%cut(gi)%n_coplanar == 0_I4P) return
   if (any(pid <= 0_I4P)) return
   drop = facet_drop(arr%facet(gi)%vertex)
   if (drop == 0_I4P) return
   do k = 1, arr%cut(gi)%n_coplanar
      do j = 1, 3
         gv(j) = ip_explicit(xyz_of(arr%facet(arr%cut(gi)%coplanar(k))%vertex(j)))
      enddo
      og = ip_orient2d(gv(1), gv(2), gv(3), drop)
      if (og == 0_I4P) cycle
      inside = .true.
      do j = 1, 3
         do e = 1, 3
            if (og * ip_orient2d(gv(e), gv(mod(e, 3) + 1), arr%pts(pid(j)), drop) < 0_I4P) then
               inside = .false. ; exit
            endif
         enddo
         if (.not. inside) exit
      enddo
      if (inside) then
         ng = ng + 1_I4P ; g(ng) = arr%cut(gi)%coplanar(k)
         if (present(first)) then
            if (first) return
         endif
      endif
   enddo
   endsubroutine coplanar_cover

   subroutine cut_sources(arr, gi, ip, iq, g, ng)
   !< Facets of the other surface whose cut segments in facet gi contain the edge (ip, iq) (registry indices): the
   !< segment and the edge are exactly collinear in gi's projection and the edge's ends are in the closed segment
   !< (exact coordinate order). One facet when the edge lies inside it, two when it lies on their common edge.
   type(arrangement_t), intent(in)  :: arr    !< Arrangement.
   integer(I4P),        intent(in)  :: gi     !< Facet.
   integer(I4P),        intent(in)  :: ip, iq !< Edge.
   integer(I4P),        intent(out) :: g(4)   !< Source facets (distinct).
   integer(I4P),        intent(out) :: ng     !< Their number (at most 4).
   integer(I4P)                     :: k, a, b, drop, ax, sab

   g = 0_I4P ; ng = 0_I4P
   drop = facet_drop(arr%facet(gi)%vertex)
   if (drop == 0_I4P) return
   do k = 1, arr%cut(gi)%n_segments
      a = arr%cut(gi)%seg(1, k) ; b = arr%cut(gi)%seg(2, k)
      if (ng > 0_I4P) then
         if (any(g(1:ng) == arr%cut(gi)%src(k))) cycle
      endif
      if (.not. (ip == a .or. ip == b .or. ip_orient2d(arr%pts(a), arr%pts(b), arr%pts(ip), drop) == 0_I4P)) cycle
      if (.not. (iq == a .or. iq == b .or. ip_orient2d(arr%pts(a), arr%pts(b), arr%pts(iq), drop) == 0_I4P)) cycle
      ax = mod(drop, 3_I4P) + 1_I4P
      sab = ip_compare(arr%pts(a), arr%pts(b), ax)
      if (sab == 0_I4P) then
         ax = mod(drop + 1_I4P, 3_I4P) + 1_I4P ; sab = ip_compare(arr%pts(a), arr%pts(b), ax)
      endif
      if (.not. within(ip) .or. .not. within(iq)) cycle
      ng = ng + 1_I4P ; g(ng) = arr%cut(gi)%src(k)
      if (ng == 4_I4P) return
   enddo
   contains
      function within(c) result(yes)
      integer(I4P), intent(in) :: c
      logical                  :: yes

      yes = c == a .or. c == b
      if (yes) return
      yes = ip_compare(arr%pts(a), arr%pts(c), ax) == sab .and. ip_compare(arr%pts(c), arr%pts(b), ax) == sab
      endfunction within
   endsubroutine cut_sources

   pure function xyz_of(a) result(r)
   type(vector_R8P), intent(in) :: a
   real(R8P)                    :: r(3)

   r = [a%x, a%y, a%z]
   endfunction xyz_of

   subroutine collect_pair(arr, gi, gj)
   !< Intersect facets gi (of A) and gj (of B) exactly (`pair_cut`, issue #65) and append a crossing segment to both
   !< facets' cut lists, its endpoints registered as implicit points. Touching pairs (one common point) add no cut. An
   !< exactly coplanar pair cuts each facet along the portions of the other facet's edges inside it (`coplanar_cut`).
   !< A pair with an exactly degenerate facet (`arr%drop` = 0, issue #99) is never cut: no contact
   !< (`degenerate_pair_meet`, exact) skips it, contact flags the arrangement (`degenerate_contact`).
   type(arrangement_t), intent(inout) :: arr    !< Arrangement.
   integer(I4P),        intent(in)    :: gi, gj !< Global facet indices (A, B).
   integer(I4P)                       :: kind   !< PAIR_* code.
   type(implicit_point_t)             :: p, q   !< Common segment.
   integer(I4P)                       :: e, ip, iq
   logical                            :: found  !< Coplanar portion found.

   if (arr%drop(gi) == 0_I4P .or. arr%drop(gj) == 0_I4P) then
      if (degenerate_pair_meet(plane_of(arr, gi), plane_of(arr, gj), drop_a=arr%drop(gi), drop_b=arr%drop(gj))) &
         arr%degenerate_contact = .true.
      return
   endif
   call pair_cut(ta=arr%facet(gi)%vertex, tb=arr%facet(gj)%vertex, kind=kind, p=p, q=q)
   if (kind == PAIR_CROSSING) then
      call register_point(arr, p, ip) ; call register_point(arr, q, iq)
      call append_segment(cut=arr%cut(gi), ip=ip, iq=iq, g=gj)
      call append_segment(cut=arr%cut(gj), ip=ip, iq=iq, g=gi)
   elseif (kind == PAIR_COPLANAR) then
      call append_partner(cut=arr%cut(gi), g=gj)
      call append_partner(cut=arr%cut(gj), g=gi)
      do e = 1, 3
         call coplanar_cut(u=arr%facet(gj)%vertex(e), v=arr%facet(gj)%vertex(mod(e, 3) + 1), &
                           t=arr%facet(gi)%vertex, found=found, p=p, q=q)
         if (found) then
            call register_point(arr, p, ip) ; call register_point(arr, q, iq)
            call append_segment(cut=arr%cut(gi), ip=ip, iq=iq, g=gj)
         endif
         call coplanar_cut(u=arr%facet(gi)%vertex(e), v=arr%facet(gi)%vertex(mod(e, 3) + 1), &
                           t=arr%facet(gj)%vertex, found=found, p=p, q=q)
         if (found) then
            call register_point(arr, p, ip) ; call register_point(arr, q, iq)
            call append_segment(cut=arr%cut(gj), ip=ip, iq=iq, g=gi)
         endif
      enddo
   endif
   endsubroutine collect_pair

   pure function bboxes_overlap(fa, fb) result(yes)
   !< Cheap AABB overlap test for two facet bounding boxes.
   type(facet_object), intent(in) :: fa, fb
   logical                        :: yes

   yes = .not. (fa%bb(2)%x < fb%bb(1)%x .or. fb%bb(2)%x < fa%bb(1)%x .or. &
                fa%bb(2)%y < fb%bb(1)%y .or. fb%bb(2)%y < fa%bb(1)%y .or. &
                fa%bb(2)%z < fb%bb(1)%z .or. fb%bb(2)%z < fa%bb(1)%z)
   endfunction bboxes_overlap

   subroutine append_segment(cut, ip, iq, g)
   !< Append a segment (registry indices) made by facet g of the other surface to a cut list, growing the buffers.
   type(cut_list_t), intent(inout) :: cut    !< Cut list to extend.
   integer(I4P),     intent(in)    :: ip, iq !< Segment endpoints.
   integer(I4P),     intent(in)    :: g      !< Source facet (global index).
   integer(I4P), allocatable       :: tmp(:, :), tmp1(:)

   if (.not. allocated(cut%seg)) then
      allocate(cut%seg(2, 4), cut%src(4))
   elseif (cut%n_segments == size(cut%seg, dim=2)) then
      allocate(tmp(2, 2 * cut%n_segments), tmp1(2 * cut%n_segments))
      tmp(:, 1:cut%n_segments) = cut%seg(:, 1:cut%n_segments)
      tmp1(1:cut%n_segments) = cut%src(1:cut%n_segments)
      call move_alloc(from=tmp, to=cut%seg) ; call move_alloc(from=tmp1, to=cut%src)
   endif
   cut%n_segments = cut%n_segments + 1_I4P
   cut%seg(:, cut%n_segments) = [ip, iq]
   cut%src(cut%n_segments) = g
   endsubroutine append_segment

   pure subroutine project_to_plane(facet, p3d, u, v)
   !< Project a 3D point `p3d` (assumed on or near `facet`'s plane) into 2D
   !< coordinates (u, v) using the facet-local orthonormal frame:
   !<
   !<   - u-axis = E12 / |E12|            (along the first edge of the facet)
   !<   - v-axis = normal x u             (in-plane, perpendicular to u)
   !<   - origin = facet%vertex(1)
   !<
   !< This is the projection used by the next-stage triangulator: the facet's
   !< three vertices map to (0, 0), (|E12|, 0), and the projected v3 — fan-
   !< triangulation in 2D is then a routine convex-polygon subdivision.
   !<
   !< @note Out-of-plane components are silently dropped. Caller is responsible
   !<       for ensuring `p3d` lies on the facet's plane within tolerance
   !<       (true by construction for cut-list endpoints, which were computed
   !<       as the plane-plane intersection segment).
   type(facet_object), intent(in)  :: facet  !< Reference facet.
   type(vector_R8P),   intent(in)  :: p3d    !< 3D point to project.
   real(R8P),          intent(out) :: u      !< u coordinate in facet-local frame.
   real(R8P),          intent(out) :: v      !< v coordinate in facet-local frame.
   type(vector_R8P)                :: e12_hat, v_axis, dp
   real(R8P)                       :: len_e12

   len_e12 = facet%E12%normL2()
   if (len_e12 <= EPS) then
      u = 0._R8P ; v = 0._R8P
      return
   endif
   e12_hat = facet%E12 * (1._R8P / len_e12)
   v_axis  = facet%normal%crossproduct(rhs=e12_hat)
   dp      = p3d - facet%vertex(1)
   u = dp%dotproduct(rhs=e12_hat)
   v = dp%dotproduct(rhs=v_axis)
   endsubroutine project_to_plane

   pure subroutine lift_from_plane(facet, u, v, p3d)
   !< Inverse of `project_to_plane`: reconstruct a 3D point on `facet`'s plane
   !< from its 2D facet-local coordinates (u, v).
   type(facet_object), intent(in)  :: facet  !< Reference facet.
   real(R8P),          intent(in)  :: u, v   !< Facet-local 2D coordinates.
   type(vector_R8P),   intent(out) :: p3d    !< Reconstructed 3D point.
   type(vector_R8P)                :: e12_hat, v_axis
   real(R8P)                       :: len_e12

   len_e12 = facet%E12%normL2()
   if (len_e12 <= EPS) then
      p3d = facet%vertex(1)
      return
   endif
   e12_hat = facet%E12 * (1._R8P / len_e12)
   v_axis  = facet%normal%crossproduct(rhs=e12_hat)
   p3d     = facet%vertex(1) + e12_hat * u + v_axis * v
   endsubroutine lift_from_plane

   subroutine arrangement_retriangulate(arr, sub_facet, sub_owner, sub_source, status, sub_pid, sub_cut)
   !< Split every facet with cuts (segments, isolated points or coplanar partners) by the exact constrained
   !< triangulation on implicit points (`fossil_cdt_implicit`, issue #65); the other facets pass through unchanged.
   !<
   !< Output (parallel arrays):
   !<   - `sub_facet(:)`  — every sub-triangle or pass-through facet, metrics computed;
   !<   - `sub_owner(:)`  — 1 if it comes from surface A, 2 from B;
   !<   - `sub_source(:)` — global index of its parent facet in `arr%facet`;
   !<   - `sub_pid(3, :)` — registry indices of its vertices;
   !<   - `sub_cut(3, :)` — edge k (opposite vertex k) lies on a cut (a constraint of the facet's triangulation).
   !<
   !< Per facet: its three vertices are points 1, 2, 3, then the distinct registry points of its cuts; the constraints
   !< are its cut segments (split at the points on them by the triangulator). The projection drops an axis where the
   !< facet is exactly non-degenerate (`facet_drop`). A sub-triangle has the orientation of its parent (the
   !< triangulator keeps that of points 1, 2, 3 exactly); its vertices are the registry points rounded once
   !< (`arr%xyz`), the same bits in every facet that shares them. A triangulation failure (crossing constraints from a
   !< self-intersecting input) returns ARR_STATUS_CDT_FAILED. Exactly degenerate facets (issue #99) are never cut: an
   !< arrangement with `degenerate_contact`, or a degenerate facet with a cut (a cut point of a neighbour inside one of
   !< its edges, propagated as an isolated point), returns ARR_STATUS_DEGENERATE and no sub-triangles; a degenerate
   !< facet without cuts passes through. A sliver closing a T-junction whose only isolated point is its own middle
   !< vertex (the facet across the T-junction is split there, by the cuts of a self-arrangement along the T-junction's
   !< collinear edges) is dropped: the split closes the T-junction, the sliver would be a fin on a non-manifold edge.
   type(arrangement_t),             intent(in)            :: arr
   type(facet_object), allocatable, intent(out)           :: sub_facet(:)
   integer(I4P),       allocatable, intent(out)           :: sub_owner(:)
   integer(I4P),       allocatable, intent(out)           :: sub_source(:)
   integer(I4P),                    intent(out), optional :: status
   integer(I4P),       allocatable, intent(out), optional :: sub_pid(:, :)
   logical,            allocatable, intent(out), optional :: sub_cut(:, :)
   integer(I4P)                                           :: gi, i, k, n_total, n, nl, ntri, st, drop, a, b
   integer(I4P), allocatable                              :: ids(:), gid(:), lsegs(:, :), tris(:, :)
   type(implicit_point_t), allocatable                    :: lp(:)
   type(facet_object)                                     :: f
   type(facet_object), allocatable                        :: buf_facet(:)
   integer(I4P),       allocatable                        :: buf_owner(:), buf_source(:), buf_pid(:, :)
   logical,            allocatable                        :: buf_cut(:, :), tcon(:, :), tmp_c(:, :)
   integer(I4P)                                           :: cap, used

   if (present(status)) status = ARR_STATUS_OK
   if (arr%degenerate_contact) then
      if (present(status)) status = ARR_STATUS_DEGENERATE
      allocate(sub_facet(0), sub_owner(0), sub_source(0))
      if (present(sub_pid)) allocate(sub_pid(3, 0))
      if (present(sub_cut)) allocate(sub_cut(3, 0))
      return
   endif
   n_total = arr%n_a + arr%n_b
   cap = max(64_I4P, n_total)
   allocate(buf_facet(cap), buf_owner(cap), buf_source(cap), buf_pid(3, cap), buf_cut(3, cap))
   used = 0_I4P
   do gi = 1, n_total
      if (arr%drop(gi) == 0_I4P) then ! exactly degenerate (issue #99): never cut
         if (.not. uncut_degenerate(arr, gi)) then
            if (present(status)) status = ARR_STATUS_DEGENERATE
            return
         endif
         if (arr%cut(gi)%n_steiner > 0_I4P) cycle ! the facet across its T-junction split at its middle vertex: dropped
         call grow_cut
         call push(buf_facet, buf_owner, buf_source, used, cap, arr%facet(gi), arr%owner(gi), gi, buf_pid)
         buf_pid(:, used) = arr%corner(:, gi) ; buf_cut(:, used) = .false.
         cycle
      endif
      if (arr%cut(gi)%n_segments == 0_I4P .and. arr%cut(gi)%n_steiner == 0_I4P .and. &
          arr%cut(gi)%n_coplanar == 0_I4P) then
         call grow_cut
         call push(buf_facet, buf_owner, buf_source, used, cap, arr%facet(gi), arr%owner(gi), gi, buf_pid)
         buf_pid(:, used) = arr%corner(:, gi) ; buf_cut(:, used) = .false.
         cycle
      endif
      call unique_ids(arr%cut(gi), ids, n, with_steiner=.true.)
      if (allocated(lp)) deallocate(lp)
      if (allocated(gid)) deallocate(gid)
      allocate(lp(3 + n), gid(3 + n))
      do i = 1, 3
         lp(i) = ip_explicit(xyz_of(arr%facet(gi)%vertex(i)))
         gid(i) = arr%corner(i, gi)
      enddo
      nl = 3_I4P
      do k = 1, n
         if (any(gid(1:3) == ids(k))) cycle
         nl = nl + 1_I4P ; lp(nl) = arr%pts(ids(k)) ; gid(nl) = ids(k)
      enddo
      if (allocated(lsegs)) deallocate(lsegs)
      allocate(lsegs(2, max(1_I4P, arr%cut(gi)%n_segments)))
      k = 0_I4P
      do i = 1, arr%cut(gi)%n_segments
         a = local(arr%cut(gi)%seg(1, i)) ; b = local(arr%cut(gi)%seg(2, i))
         if (a == b .or. a == 0_I4P .or. b == 0_I4P) cycle
         k = k + 1_I4P ; lsegs(:, k) = [a, b]
      enddo
      drop = facet_drop(arr%facet(gi)%vertex)
      call cdt_facet(pts=lp(1:nl), drop=drop, segs=lsegs(:, 1:k), tris=tris, ntri=ntri, status=st, tcon=tcon)
      if (st /= CDTI_OK) then
         if (present(status)) status = ARR_STATUS_CDT_FAILED
         return
      endif
      do k = 1, ntri
         do i = 1, 3
            if (tris(i, k) <= 3_I4P) then
               f%vertex(i) = arr%facet(gi)%vertex(tris(i, k))
            else
               f%vertex(i)%x = arr%xyz(1, gid(tris(i, k)))
               f%vertex(i)%y = arr%xyz(2, gid(tris(i, k)))
               f%vertex(i)%z = arr%xyz(3, gid(tris(i, k)))
            endif
         enddo
         call f%compute_metrix
         f%patch = arr%facet(gi)%patch ! a sub-triangle keeps its parent's patch (issue #66)
         call grow_cut
         call push(buf_facet, buf_owner, buf_source, used, cap, f, arr%owner(gi), gi, buf_pid)
         buf_pid(:, used) = gid(tris(:, k)) ; buf_cut(:, used) = tcon(:, k)
      enddo
   enddo
   allocate(sub_facet(used), sub_owner(used), sub_source(used))
   sub_facet (1:used) = buf_facet (1:used)
   sub_owner (1:used) = buf_owner (1:used)
   sub_source(1:used) = buf_source(1:used)
   if (present(sub_pid)) then
      allocate(sub_pid(3, used)) ; sub_pid = buf_pid(:, 1:used)
   endif
   if (present(sub_cut)) then
      allocate(sub_cut(3, used)) ; sub_cut = buf_cut(:, 1:used)
   endif
   contains
      subroutine grow_cut
      !< Keep `buf_cut` as long as the accumulator (grown by `push` when full).
      if (used == cap) then
         allocate(tmp_c(3, 2 * cap)) ; tmp_c(:, 1:used) = buf_cut(:, 1:used) ; call move_alloc(tmp_c, buf_cut)
      endif
      endsubroutine grow_cut

      pure function local(id) result(l)
      !< Local index of registry point id in the facet's list (0 if absent).
      integer(I4P), intent(in) :: id
      integer(I4P)             :: l

      do l = 1, nl
         if (gid(l) == id) return
      enddo
      l = 0_I4P
      endfunction local
   endsubroutine arrangement_retriangulate

   pure function uncut_degenerate(arr, gi) result(yes)
   !< Facet gi is exactly degenerate and uncut (issue #99): no segments, no partners, and its isolated points are its
   !< own vertices (the middle vertex of a sliver closing a T-junction, propagated back to it when the facet across
   !< the T-junction is split there).
   type(arrangement_t), intent(in) :: arr !< Arrangement.
   integer(I4P),        intent(in) :: gi  !< Facet.
   logical                         :: yes !< Degenerate and uncut.
   integer(I4P)                    :: k

   yes = arr%drop(gi) == 0_I4P .and. arr%cut(gi)%n_segments == 0_I4P .and. arr%cut(gi)%n_coplanar == 0_I4P
   if (.not. yes) return
   do k = 1, arr%cut(gi)%n_steiner
      yes = any(arr%corner(:, gi) == arr%cut(gi)%steiner(k))
      if (.not. yes) return
   enddo
   endfunction uncut_degenerate

   subroutine push(buf_facet, buf_owner, buf_source, used, cap, f, owner, source, buf_pid)
   !< Append (f, owner, source) to the parallel-array accumulator, doubling
   !< capacity on overflow.
   type(facet_object), allocatable, intent(inout) :: buf_facet(:)
   integer(I4P),       allocatable, intent(inout) :: buf_owner(:), buf_source(:)
   integer(I4P),                    intent(inout) :: used, cap
   type(facet_object),              intent(in)    :: f
   integer(I4P),                    intent(in)    :: owner, source
   integer(I4P),       allocatable, intent(inout) :: buf_pid(:, :)
   type(facet_object), allocatable                :: tmp_f(:)
   integer(I4P),       allocatable                :: tmp_i(:), tmp_p(:, :)

   if (used == cap) then
      allocate(tmp_p(3, 2 * cap))
      tmp_p(:, 1:used) = buf_pid(:, 1:used)
      call move_alloc(from=tmp_p, to=buf_pid)
      cap = 2 * cap
      allocate(tmp_f(cap))
      tmp_f(1:used) = buf_facet(1:used)
      call move_alloc(from=tmp_f, to=buf_facet)
      allocate(tmp_i(cap))
      tmp_i(1:used) = buf_owner(1:used)
      call move_alloc(from=tmp_i, to=buf_owner)
      allocate(tmp_i(cap))
      tmp_i(1:used) = buf_source(1:used)
      call move_alloc(from=tmp_i, to=buf_source)
   endif
   used = used + 1
   buf_facet (used) = f
   buf_owner (used) = owner
   buf_source(used) = source
   endsubroutine push

endmodule fossil_arrangement
