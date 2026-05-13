!< FOSSIL, mesh-arrangement scaffolding for booleans (issue #18 §1.1 stage 1).

module fossil_arrangement
!< FOSSIL, mesh-arrangement scaffolding for booleans (issue #18 §1.1 stage 1).
!<
!< This module implements the **input side** of the Mesh Arrangements pipeline
!< of Zhou, Grinspun, Zorin & Jacobson (*Mesh Arrangements for Solid Geometry*,
!< SIGGRAPH 2016): given two input surfaces A and B, compute for every facet
!< the list of intersection segments where it crosses any facet of the other
!< surface. The output of this stage — a per-facet list of "cuts" — is the
!< input to the triangulation + WN-tagging + selection stages that follow.
!<
!< Pipeline overview (full §1.1 plan, for context):
!<
!<   stage 1 (this module):
!<     A, B  -->  arrangement_initialize  -->  arrangement (facets + owners + empty cuts)
!<     arrangement  -->  arrangement_collect_intersections  -->  cuts populated
!<
!<   stage 2 (next PR — triangulation):
!<     for each facet f with non-empty cuts:
!<       project f and its cuts to 2D (project_to_plane / lift_from_plane)
!<       fan-triangulate from a Steiner point, inserting cut endpoints as
!<         constraints (or run a CDT — see staging note in issue #18)
!<     output: a flat list of sub-triangles, each tagged by its source facet
!<
!<   stage 3 (next PR — tagging):
!<     for each sub-triangle t with centroid c:
!<       tag_A(t) = inside(c, A)  via winding_number(c, A) > 0.5
!<       tag_B(t) = inside(c, B)  via winding_number(c, B) > 0.5
!<
!<   stage 4 (next PR — selection):
!<     UNION:       keep t if (tag_A xor tag_B) or (boundary on outside of other)
!<     INTERSECT:   keep t if (tag_A and tag_B) ...
!<     DIFFERENCE:  keep t if (from A and not in B) or (from B and in A, flipped)
!<     SYMDIFF:     keep t if exactly one of (tag_A, tag_B) is true
!<
!<   stage 5 (next PR — stitch):
!<     adopt selected triangles back into a surface_stl_object, run sanitize_normals.
!<
!< This module now provides stages 1 and 2 (initialization, intersection
!< collection, retriangulation). Stages 3-5 (WN tagging, op-specific
!< selection, stitching) live in the boolean driver — see
!< `surface_stl_object%boolean`.

use fossil_aabb_node_object, only : aabb_node_object
use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_dt,               only : triangulation_t, cdt_build, &
                                    DT_STATUS_OK, DT_STATUS_RECOVERY_FAILED, &
                                    DT_STATUS_CROSSING_CONSTRAINTS
use fossil_facet_object,     only : facet_object
use fossil_list_id_object,   only : list_id_object
use fossil_utils,            only : EPS
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: cut_list_t
public :: arrangement_t
public :: arrangement_initialize
public :: arrangement_collect_intersections
public :: arrangement_retriangulate
public :: project_to_plane
public :: lift_from_plane
public :: ARR_STATUS_OK, ARR_STATUS_CDT_FAILED

integer(I4P), parameter :: MAX_CHILDREN = 8_I4P  !< Buffer size for `enumerate_children` (octree fan-out).

! Status codes for arrangement_retriangulate.
integer(I4P), parameter :: ARR_STATUS_OK         = 0_I4P  !< Retriangulation completed successfully.
integer(I4P), parameter :: ARR_STATUS_CDT_FAILED = 1_I4P  !< CDT could not recover one or more constraint segments.

type :: cut_list_t
   !< Per-facet list of intersection-segment endpoints in **3D world** coordinates.
   !<
   !< Each consecutive pair `(point(2*k-1), point(2*k))` is one intersection
   !< segment. We do not store segments as a separate `segment_t` type because
   !< the next-stage triangulator only ever consumes the flat point list (each
   !< endpoint becomes a Steiner point in the facet plane), so the pair
   !< grouping is implicit in the index parity.
   !<
   !< `n_segments` counts segments, not points; total points in `point` is
   !< `2 * n_segments` (with `point` over-allocated to amortize append cost).
   integer(I4P)                  :: n_segments = 0_I4P  !< Number of intersection segments stored.
   type(vector_R8P), allocatable :: point(:)            !< Flat 3D point buffer; pairs of consecutive entries form segments.
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

   do i = 1, arr%n_a
      arr%facet(i) = facet_a(i)
      arr%owner(i) = 1_I4P
   enddo
   do i = 1, arr%n_b
      arr%facet(arr%n_a + i) = facet_b(i)
      arr%owner(arr%n_a + i) = 2_I4P
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
      return
   endif

   call traverse_cross_pair(arr=arr, na=0_I4P, nb=0_I4P)
   endsubroutine arrangement_collect_intersections

   subroutine collect_brute_force(arr)
   !< O(nA * nB) fallback used only when trees are not initialized.
   type(arrangement_t), intent(inout) :: arr
   integer(I4P)                       :: i, j, gi, gj
   type(vector_R8P)                   :: p, q
   logical                            :: hit

   do i = 1, arr%n_a
      gi = i
      do j = 1, arr%n_b
         gj = arr%n_a + j
         if (.not. bboxes_overlap(arr%facet(gi), arr%facet(gj))) cycle
         call arr%facet(gi)%intersect_facet(other=arr%facet(gj), p=p, q=q, intersects=hit)
         if (hit) then
            call append_segment(cut=arr%cut(gi), p=p, q=q)
            call append_segment(cut=arr%cut(gj), p=p, q=q)
         endif
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
   type(vector_R8P)                   :: p, q
   logical                            :: hit

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
         call arr%facet(gi)%intersect_facet(other=arr%facet(gj), p=p, q=q, intersects=hit)
         if (hit) then
            call append_segment(cut=arr%cut(gi), p=p, q=q)
            call append_segment(cut=arr%cut(gj), p=p, q=q)
         endif
      enddo
   enddo
   endsubroutine enumerate_leaf_cross_pair

   pure function bboxes_overlap(fa, fb) result(yes)
   !< Cheap AABB overlap test for two facet bounding boxes.
   type(facet_object), intent(in) :: fa, fb
   logical                        :: yes

   yes = .not. (fa%bb(2)%x < fb%bb(1)%x .or. fb%bb(2)%x < fa%bb(1)%x .or. &
                fa%bb(2)%y < fb%bb(1)%y .or. fb%bb(2)%y < fa%bb(1)%y .or. &
                fa%bb(2)%z < fb%bb(1)%z .or. fb%bb(2)%z < fa%bb(1)%z)
   endfunction bboxes_overlap

   subroutine append_segment(cut, p, q)
   !< Append a 3D intersection segment (p, q) to a cut list, growing the
   !< buffer geometrically (initial capacity 4 segments = 8 points, doubling
   !< on overflow).
   type(cut_list_t), intent(inout) :: cut  !< Cut list to extend.
   type(vector_R8P), intent(in)    :: p, q !< Segment endpoints.
   type(vector_R8P), allocatable   :: tmp(:)
   integer(I4P)                    :: cap, used

   used = 2 * cut%n_segments
   if (.not. allocated(cut%point)) then
      allocate(cut%point(8))  ! capacity 4 segments
   else
      cap = size(cut%point, dim=1)
      if (used + 2 > cap) then
         allocate(tmp(2 * cap))
         tmp(1:used) = cut%point(1:used)
         call move_alloc(from=tmp, to=cut%point)
      endif
   endif
   cut%point(used + 1) = p
   cut%point(used + 2) = q
   cut%n_segments      = cut%n_segments + 1
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

   subroutine arrangement_retriangulate(arr, sub_facet, sub_owner, sub_source, status)
   !< Walk every facet in `arr`; for those with non-empty cuts, retriangulate
   !< them along the cut segments using the constrained Delaunay triangulator
   !< (`fossil_dt%cdt_build`). Facets with no cuts pass through unchanged.
   !<
   !< Output is a flat parallel-array trio:
   !<   - `sub_facet(:)`  — every emitted sub-triangle (or pass-through facet),
   !<                       with `compute_metrix` already called.
   !<   - `sub_owner(:)`  — 1 if the sub-triangle came from surface A, 2 from B
   !<                       (inherited from the parent facet's owner tag).
   !<   - `sub_source(:)` — global index of the parent facet in `arr%facet`,
   !<                       so callers can trace each sub-triangle back to its
   !<                       origin if needed.
   !<
   !< Per-facet algorithm:
   !<   1. Project the 3 facet vertices into the facet's local 2D frame
   !<      (project_to_plane). They become CDT input points 1, 2, 3.
   !<   2. For each cut segment (p, q): project both endpoints to 2D, dedup
   !<      against existing input points (within EPS in 2D), append the
   !<      surviving new points; record (idx_p, idx_q) as a CDT constraint.
   !<   3. Run cdt_build with the 2D points + segment constraints. Status
   !<      DT_STATUS_RECOVERY_FAILED bubbles up as ARR_STATUS_CDT_FAILED so
   !<      the boolean driver can decide how to handle it.
   !<   4. For each output sub-triangle, lift its 3 vertices back to 3D via
   !<      lift_from_plane, build a fresh facet_object, copy the parent's
   !<      normal (the lift preserves the plane, so all sub-facets share the
   !<      parent's normal direction), call compute_metrix, append to output.
   !<
   !< Pass-through (empty cuts) emits the parent facet unchanged. We deep-copy
   !< via intrinsic assignment so downstream selection can mutate without
   !< aliasing back into `arr`.
   type(arrangement_t),             intent(in)            :: arr
   type(facet_object), allocatable, intent(out)           :: sub_facet(:)
   integer(I4P),       allocatable, intent(out)           :: sub_owner(:)
   integer(I4P),       allocatable, intent(out)           :: sub_source(:)
   integer(I4P),                    intent(out), optional :: status
   integer(I4P)                                           :: gi, i, k, n_total, n_pts, n_segs, n_cuts
   integer(I4P)                                           :: cdt_status, sub_v(3), n_sub
   real(R8P), allocatable                                 :: pts2d(:,:)
   integer(I4P), allocatable                              :: segs(:,:)
   type(triangulation_t)                                  :: tri
   type(vector_R8P)                                       :: p3d
   real(R8P)                                              :: u, v
   integer(I4P)                                           :: idx_p, idx_q
   ! Output accumulator, doubled on demand.
   type(facet_object), allocatable                        :: buf_facet(:)
   integer(I4P),       allocatable                        :: buf_owner(:), buf_source(:)
   integer(I4P)                                           :: cap, used

   if (present(status)) status = ARR_STATUS_OK

   ! Initial accumulator capacity — most facets pass through, so we grow
   ! geometrically from the original facet count.
   n_total = arr%n_a + arr%n_b
   cap = max(64_I4P, n_total)
   allocate(buf_facet(cap), buf_owner(cap), buf_source(cap))
   used = 0_I4P

   do gi = 1, n_total
      n_cuts = arr%cut(gi)%n_segments

      if (n_cuts == 0_I4P) then
         ! Pass-through.
         call push(buf_facet, buf_owner, buf_source, used, cap, &
                   arr%facet(gi), arr%owner(gi), gi)
         cycle
      endif

      ! --- Build CDT input: 3 facet vertices + cut endpoints (deduped). ---
      ! Maximum points = 3 + 2 * n_cuts ; segments = n_cuts.
      n_pts = 0_I4P ; n_segs = 0_I4P
      if (allocated(pts2d)) deallocate(pts2d)
      if (allocated(segs))  deallocate(segs)
      allocate(pts2d(2, 3 + 2 * n_cuts))
      allocate(segs(2, n_cuts))

      ! Three facet vertices first (so they keep stable indices 1, 2, 3).
      do i = 1, 3
         call project_to_plane(facet=arr%facet(gi), p3d=arr%facet(gi)%vertex(i), &
                               u=pts2d(1, i), v=pts2d(2, i))
      enddo
      n_pts = 3_I4P

      ! Cut endpoints: project, dedup against existing points (within EPS),
      ! record (start_idx, end_idx) as a constraint segment.
      do k = 1, n_cuts
         call project_to_plane(facet=arr%facet(gi), p3d=arr%cut(gi)%point(2*k - 1), u=u, v=v)
         idx_p = find_or_append_2d(pts2d=pts2d, n_pts=n_pts, x=u, y=v)
         call project_to_plane(facet=arr%facet(gi), p3d=arr%cut(gi)%point(2*k), u=u, v=v)
         idx_q = find_or_append_2d(pts2d=pts2d, n_pts=n_pts, x=u, y=v)
         if (idx_p == idx_q) cycle  ! degenerate zero-length cut after dedup
         n_segs = n_segs + 1
         segs(1, n_segs) = idx_p
         segs(2, n_segs) = idx_q
      enddo

      ! --- Run CDT. ---
      call cdt_build(tri=tri, points=pts2d(:, 1:n_pts), segments=segs(:, 1:n_segs), status=cdt_status)
      if (cdt_status /= DT_STATUS_OK) then
         ! Recovery failed for this facet's cuts. We could fall back to the
         ! original facet (better than nothing) or propagate the error. The
         ! boolean driver treats this as a global failure so the user knows
         ! the result is incorrect rather than silently degraded.
         if (present(status)) status = ARR_STATUS_CDT_FAILED
         return
      endif

      ! --- Emit sub-triangles. ---
      n_sub = tri%num_triangles()
      do k = 1, n_sub
         call tri%triangle_vertices(k=k, v=sub_v)
         call build_subfacet_from_uv(arr=arr, gi=gi, pts2d=pts2d, sub_v=sub_v, &
                                     buf_facet=buf_facet, buf_owner=buf_owner, &
                                     buf_source=buf_source, used=used, cap=cap)
      enddo
   enddo

   ! Trim accumulators to exact size.
   allocate(sub_facet (used))
   allocate(sub_owner (used))
   allocate(sub_source(used))
   sub_facet (1:used) = buf_facet (1:used)
   sub_owner (1:used) = buf_owner (1:used)
   sub_source(1:used) = buf_source(1:used)
   endsubroutine arrangement_retriangulate

   function find_or_append_2d(pts2d, n_pts, x, y) result(idx)
   !< Look up (x, y) in `pts2d(:, 1:n_pts)` within tolerance; if found return
   !< its index, otherwise append and return the new index.
   !<
   !< Tolerance: scale-aware — `max(|coord|, 1) * DEDUP_REL_TOL`. Cannot use
   !< the project-wide `EPS` because it is hardcoded to 0 in `fossil_utils`
   !< (a known pre-existing oddity), which would defeat the whole purpose of
   !< this dedup pass. The relative tolerance `1e-10` accommodates the
   !< inevitable round-off in `project_to_plane` while still distinguishing
   !< genuinely different points whose 2D coordinates differ at the eighth
   !< significant digit or further apart.
   real(R8P),    intent(inout) :: pts2d(:, :)
   integer(I4P), intent(inout) :: n_pts
   real(R8P),    intent(in)    :: x, y
   integer(I4P)                :: idx
   integer(I4P)                :: i
   real(R8P)                   :: tol_x, tol_y
   real(R8P), parameter        :: DEDUP_REL_TOL = 1.0e-10_R8P

   tol_x = max(abs(x), 1._R8P) * DEDUP_REL_TOL
   tol_y = max(abs(y), 1._R8P) * DEDUP_REL_TOL
   do i = 1, n_pts
      if (abs(pts2d(1, i) - x) <= tol_x .and. abs(pts2d(2, i) - y) <= tol_y) then
         idx = i ; return
      endif
   enddo
   n_pts = n_pts + 1
   pts2d(1, n_pts) = x
   pts2d(2, n_pts) = y
   idx = n_pts
   endfunction find_or_append_2d

   subroutine build_subfacet_from_uv(arr, gi, pts2d, sub_v, buf_facet, buf_owner, &
                                     buf_source, used, cap)
   !< Lift one CDT sub-triangle (vertex indices into pts2d) back to 3D as a
   !< facet, inheriting the parent facet's normal, then push to the output
   !< accumulator.
   type(arrangement_t),             intent(in)    :: arr
   integer(I4P),                    intent(in)    :: gi
   real(R8P),                       intent(in)    :: pts2d(:, :)
   integer(I4P),                    intent(in)    :: sub_v(3)
   type(facet_object), allocatable, intent(inout) :: buf_facet(:)
   integer(I4P),       allocatable, intent(inout) :: buf_owner(:), buf_source(:)
   integer(I4P),                    intent(inout) :: used, cap
   type(facet_object)                             :: f
   type(vector_R8P)                               :: p3d
   integer(I4P)                                   :: i

   do i = 1, 3
      call lift_from_plane(facet=arr%facet(gi), &
                           u=pts2d(1, sub_v(i)), v=pts2d(2, sub_v(i)), p3d=p3d)
      f%vertex(i) = p3d
   enddo
   ! Inherit parent's normal direction (the lift preserves the plane, so the
   ! geometric normal will be ±parent's normal; copying it sets the orientation
   ! correctly for sub-triangles emitted in the same CCW sense as the CDT).
   f%normal = arr%facet(gi)%normal
   call f%compute_metrix
   ! `compute_metrix` recomputes the normal from vertex order; if it ended up
   ! anti-parallel to the parent (i.e. the CDT emitted the sub-triangle in
   ! CW order relative to the parent's outward normal), flip to restore.
   if (f%normal%dotproduct(rhs=arr%facet(gi)%normal) < 0._R8P) then
      call f%reverse_normal
      call f%compute_metrix
   endif
   call push(buf_facet, buf_owner, buf_source, used, cap, f, arr%owner(gi), gi)
   endsubroutine build_subfacet_from_uv

   subroutine push(buf_facet, buf_owner, buf_source, used, cap, f, owner, source)
   !< Append (f, owner, source) to the parallel-array accumulator, doubling
   !< capacity on overflow.
   type(facet_object), allocatable, intent(inout) :: buf_facet(:)
   integer(I4P),       allocatable, intent(inout) :: buf_owner(:), buf_source(:)
   integer(I4P),                    intent(inout) :: used, cap
   type(facet_object),              intent(in)    :: f
   integer(I4P),                    intent(in)    :: owner, source
   type(facet_object), allocatable                :: tmp_f(:)
   integer(I4P),       allocatable                :: tmp_i(:)

   if (used == cap) then
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
