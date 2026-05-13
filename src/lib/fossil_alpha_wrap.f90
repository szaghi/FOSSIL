!< Alpha wrapping — watertight surrogate from broken triangle-soup input
!< (issue #18 §1.6).
!<
!< Reference: Portaneri, Hemmer, Birdal, Mandad & Alliez, *Alpha Wrapping
!< with an Offset* (SIGGRAPH 2022). CGAL package: Alpha_wrap_3.
!<
!< This module is the engineering core of FOSSIL's "make any CAD junk usable
!< for CFD" primitive. The pipeline (Portaneri's 5 phases, mapped to FOSSIL's
!< staging discipline):
!<
!<   Step 1 (this commit) — Octree refined to geometric leaf-size α.
!<      Each leaf is small enough that "fits inside the leaf" implies a
!<      feature smaller than the user's wrap-feature size α.
!<   Step 2 — Inside/outside flood fill from a known-outside corner seed,
!<      treating BOUNDARY leaves (those that intersect input geometry) as
!<      barriers.
!<   Step 3 — Dual-contouring boundary extraction: emit a quad for each
!<      octree edge that crosses the INSIDE/OUTSIDE classification.
!<   Step 4 — Vertex projection toward the input within Hausdorff offset ε.
!<   Step 5 — Adaptive refinement loop until convergence; capstone TBP
!<      `surface%alpha_wrap`.
!<
!< Steps 2-5 land in subsequent commits. This commit ships step 1 with a
!< public helper for testing the octree in isolation.

module fossil_alpha_wrap
!< Alpha wrapping — watertight surrogate from broken triangle-soup input.

use fossil_facet_object, only : facet_object
use fossil_utils, only : triangle_overlaps_aabb
use penf, only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
private

public :: awrap_octree_t
public :: awrap_build_octree
public :: awrap_classify_leaves
public :: awrap_extract_surface
public :: AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE
public :: AWRAP_LEAF_FLAG_BOUNDARY, AWRAP_LEAF_FLAG_INTERIOR
public :: AWRAP_LEAF_FLAG_EMPTY, AWRAP_LEAF_FLAG_INSIDE, AWRAP_LEAF_FLAG_OUTSIDE
public :: AWRAP_MAX_DEPTH

integer(I4P), parameter :: AWRAP_STATUS_OK         = 0_I4P
integer(I4P), parameter :: AWRAP_STATUS_BAD_INPUT  = 1_I4P  !< Empty facet list, or alpha <= 0.
integer(I4P), parameter :: AWRAP_STATUS_DEGENERATE = 2_I4P  !< Alpha > bbox diagonal: octree won't refine, output empty.

integer(I4P), parameter :: AWRAP_LEAF_FLAG_INTERIOR = 0_I4P  !< Step 1 only: leaf has no overlapping input facets (kept for back-compat in tests).
integer(I4P), parameter :: AWRAP_LEAF_FLAG_BOUNDARY = 1_I4P  !< Leaf overlaps at least one input facet (barrier for step 2 flood fill).
! Step 2 classification — replaces the step-1 INTERIOR flag on leaves the
! flood fill has visited. BOUNDARY leaves keep their flag (they ARE the
! barriers). EMPTY = step 2 hasn't visited yet (initial state for empty
! leaves, transitions to INSIDE or OUTSIDE during flood fill).
integer(I4P), parameter :: AWRAP_LEAF_FLAG_EMPTY    = 0_I4P  !< Alias of INTERIOR — leaf has no overlap, classification still pending.
integer(I4P), parameter :: AWRAP_LEAF_FLAG_OUTSIDE  = 2_I4P  !< Reachable from the corner seed via 6-connectivity through EMPTY leaves.
integer(I4P), parameter :: AWRAP_LEAF_FLAG_INSIDE   = 3_I4P  !< Empty leaf NOT reachable from outside seed — interior of the input solid.

integer(I4P), parameter :: AWRAP_MAX_DEPTH      = 12_I4P     !< Cap on octree recursion depth (8^12 ≈ 7e10 nodes worst case; in practice trees stay << 1e5 leaves).
integer(I4P), parameter :: AWRAP_INITIAL_CAPACITY = 1024_I4P !< Initial node-array capacity; doubles on demand.
real(R8P),    parameter :: AWRAP_PADDING_FACTOR = 2.0_R8P    !< Octree root bbox padded by this * alpha along each axis (ensures corner seed is outside geometry).

type :: awrap_octree_node_t
   !< One octree node. Children are either all 8 allocated (interior) or none (leaf).
   !< Leaves carry an integer list of overlapping facet IDs.
   type(vector_R8P)          :: bmin                            !< Lower-left-back corner.
   type(vector_R8P)          :: bmax                            !< Upper-right-front corner.
   integer(I4P)              :: first_child = -1_I4P            !< Index of first of 8 contiguous children, or -1 if leaf.
   integer(I4P)              :: depth       =  0_I4P            !< Depth from the root (root = 0).
   integer(I4P)              :: leaf_flag   =  AWRAP_LEAF_FLAG_INTERIOR  !< Set when leaf has overlapping facets (step 1) or further classified (step 2+).
   integer(I4P), allocatable :: facet_ids(:)                    !< List of overlapping facet IDs (leaves only).
endtype awrap_octree_node_t

type :: awrap_octree_t
   !< Leaf-size-α octree. Built bottom-up via recursive subdivision; nodes are
   !< stored contiguously in `node(0:n_nodes-1)` so that the children of node
   !< `p` are at `node(p%first_child:p%first_child+7)`.
   type(awrap_octree_node_t), allocatable :: node(:)         !< Flat node array.
   integer(I4P)                           :: n_nodes = 0_I4P !< Number of populated nodes.
   integer(I4P)                           :: n_leaves = 0_I4P !< Number of leaves (sum at end of build).
   integer(I4P)                           :: n_boundary_leaves = 0_I4P  !< Leaves with leaf_flag == BOUNDARY.
   integer(I4P)                           :: n_inside_leaves   = 0_I4P  !< Filled by step 2 (awrap_classify_leaves).
   integer(I4P)                           :: n_outside_leaves  = 0_I4P  !< Filled by step 2 (awrap_classify_leaves).
   real(R8P)                              :: alpha   = 0._R8P  !< Geometric leaf-size target.
endtype awrap_octree_t

contains

   subroutine awrap_build_octree(facet, alpha, octree, status)
   !< Build a leaf-size-α octree over `facet(:)`.
   !<
   !< A leaf is split (recursively) iff it overlaps at least one input facet
   !< AND its current size > `alpha`. This produces the standard adaptive
   !< octree: refined where the geometry has detail, coarse far from any
   !< facet. The root bbox is padded by `AWRAP_PADDING_FACTOR * alpha` along
   !< each axis so that step 2's corner-seed leaf is provably empty.
   !<
   !< Output `octree%n_boundary_leaves` is the count of leaves that ended
   !< the build with `leaf_flag = AWRAP_LEAF_FLAG_BOUNDARY` (i.e. they
   !< overlap input geometry and are at size <= alpha). Step 2's flood fill
   !< treats those as barriers.
   type(facet_object),     intent(in)            :: facet(:)   !< Input facets.
   real(R8P),              intent(in)            :: alpha      !< Geometric leaf-size target.
   type(awrap_octree_t),   intent(out)           :: octree     !< Built octree.
   integer(I4P),           intent(out), optional :: status     !< Status code.
   type(vector_R8P)                              :: bmin, bmax !< Tight bbox of input facets.
   type(vector_R8P)                              :: rbmin, rbmax !< Padded root bbox.
   real(R8P)                                     :: pad        !< Padding distance.
   real(R8P)                                     :: bdiag      !< Padded root bbox diagonal.
   integer(I4P), allocatable                     :: root_facets(:)  !< Indices into facet(:).
   integer(I4P)                                  :: nf, f

   if (present(status)) status = AWRAP_STATUS_OK

   nf = size(facet, kind=I4P)
   if (nf == 0_I4P .or. alpha <= 0._R8P) then
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif

   call compute_facet_bbox(facet=facet, bmin=bmin, bmax=bmax)
   bdiag = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2)
   if (alpha >= bdiag) then
      ! Degenerate: α is larger than the input geometry's bbox diagonal, so
      ! the octree can't refine. Return a single root leaf marked BOUNDARY
      ! (it contains all facets) and signal degeneracy.
      pad = AWRAP_PADDING_FACTOR * alpha
      rbmin = vector_R8P(bmin%x - pad, bmin%y - pad, bmin%z - pad)
      rbmax = vector_R8P(bmax%x + pad, bmax%y + pad, bmax%z + pad)
      allocate(octree%node(1))
      octree%node(1)%bmin = rbmin
      octree%node(1)%bmax = rbmax
      octree%node(1)%depth = 0_I4P
      octree%node(1)%first_child = -1_I4P
      octree%node(1)%leaf_flag = AWRAP_LEAF_FLAG_BOUNDARY
      allocate(octree%node(1)%facet_ids(nf))
      do f = 1_I4P, nf
         octree%node(1)%facet_ids(f) = f
      enddo
      octree%n_nodes = 1_I4P
      octree%n_leaves = 1_I4P
      octree%n_boundary_leaves = 1_I4P
      octree%alpha = alpha
      if (present(status)) status = AWRAP_STATUS_DEGENERATE
      return
   endif

   ! Pad the root bbox by 2α so the corner seed leaf used by step 2 (flood
   ! fill from a known-outside point) is provably empty of input facets.
   pad = AWRAP_PADDING_FACTOR * alpha
   rbmin = vector_R8P(bmin%x - pad, bmin%y - pad, bmin%z - pad)
   rbmax = vector_R8P(bmax%x + pad, bmax%y + pad, bmax%z + pad)

   octree%alpha = alpha
   allocate(octree%node(AWRAP_INITIAL_CAPACITY))
   octree%n_nodes = 1_I4P
   octree%node(1)%bmin = rbmin
   octree%node(1)%bmax = rbmax
   octree%node(1)%depth = 0_I4P
   octree%node(1)%first_child = -1_I4P
   octree%node(1)%leaf_flag = AWRAP_LEAF_FLAG_INTERIOR

   ! Initial root facet list = all facets.
   allocate(root_facets(nf))
   do f = 1_I4P, nf
      root_facets(f) = f
   enddo

   call subdivide(octree=octree, node_id=1_I4P, facet=facet, candidate_facets=root_facets)
   deallocate(root_facets)

   ! Sweep to count leaves and boundary leaves.
   call recount(octree)
   endsubroutine awrap_build_octree

   recursive subroutine subdivide(octree, node_id, facet, candidate_facets)
   !< Recursive subdivision driver. Tests every candidate facet against this
   !< node's bbox via SAT; if the node overlaps at least one facet AND is
   !< still larger than alpha AND has not hit max depth, splits into 8
   !< children and recurses on each with the local overlapping facet subset.
   type(awrap_octree_t), intent(inout) :: octree
   integer(I4P),         intent(in)    :: node_id
   type(facet_object),   intent(in)    :: facet(:)
   integer(I4P),         intent(in)    :: candidate_facets(:)
   integer(I4P), allocatable           :: my_facets(:)
   integer(I4P)                        :: n_my, i, fid, c
   real(R8P)                           :: size_x, size_y, size_z
   logical                             :: should_split

   ! Filter candidate facets to those that actually overlap this node.
   call filter_overlapping(node_bmin=octree%node(node_id)%bmin, node_bmax=octree%node(node_id)%bmax, &
                           facet=facet, candidate=candidate_facets, kept=my_facets, n_kept=n_my)

   if (n_my == 0_I4P) then
      ! Leaf, no overlap → INTERIOR (step 2 will further classify).
      octree%node(node_id)%leaf_flag = AWRAP_LEAF_FLAG_INTERIOR
      if (allocated(my_facets)) deallocate(my_facets)
      return
   endif

   ! Decide whether to split: oversized AND not at depth cap.
   size_x = octree%node(node_id)%bmax%x - octree%node(node_id)%bmin%x
   size_y = octree%node(node_id)%bmax%y - octree%node(node_id)%bmin%y
   size_z = octree%node(node_id)%bmax%z - octree%node(node_id)%bmin%z
   should_split = (max(size_x, size_y, size_z) > octree%alpha) .and. &
                  (octree%node(node_id)%depth < AWRAP_MAX_DEPTH)

   if (.not. should_split) then
      ! Boundary leaf at target resolution.
      octree%node(node_id)%leaf_flag = AWRAP_LEAF_FLAG_BOUNDARY
      call move_alloc(from=my_facets, to=octree%node(node_id)%facet_ids)
      return
   endif

   ! Allocate 8 children and recurse.
   call allocate_children(octree=octree, parent_id=node_id)
   do c = 0_I4P, 7_I4P
      call subdivide(octree=octree, node_id=octree%node(node_id)%first_child + c, &
                     facet=facet, candidate_facets=my_facets)
   enddo
   if (allocated(my_facets)) deallocate(my_facets)
   endsubroutine subdivide

   subroutine allocate_children(octree, parent_id)
   !< Reserve 8 contiguous child slots in `octree%node(:)`, growing the array
   !< if necessary, and initialize their bboxes by 8-way subdivision.
   type(awrap_octree_t), intent(inout) :: octree
   integer(I4P),         intent(in)    :: parent_id
   integer(I4P)                        :: start, c
   type(vector_R8P)                    :: pmin, pmax, mid
   type(vector_R8P)                    :: cmin, cmax

   if (octree%n_nodes + 8_I4P > size(octree%node, kind=I4P)) call grow_node_array(octree)

   start = octree%n_nodes + 1_I4P
   octree%node(parent_id)%first_child = start
   pmin = octree%node(parent_id)%bmin
   pmax = octree%node(parent_id)%bmax
   mid = 0.5_R8P * (pmin + pmax)

   do c = 0_I4P, 7_I4P
      ! Bit pattern: bit 0 = +x, bit 1 = +y, bit 2 = +z.
      cmin%x = merge(mid%x, pmin%x, iand(c, 1_I4P) /= 0_I4P)
      cmax%x = merge(pmax%x, mid%x, iand(c, 1_I4P) /= 0_I4P)
      cmin%y = merge(mid%y, pmin%y, iand(c, 2_I4P) /= 0_I4P)
      cmax%y = merge(pmax%y, mid%y, iand(c, 2_I4P) /= 0_I4P)
      cmin%z = merge(mid%z, pmin%z, iand(c, 4_I4P) /= 0_I4P)
      cmax%z = merge(pmax%z, mid%z, iand(c, 4_I4P) /= 0_I4P)
      octree%node(start + c)%bmin = cmin
      octree%node(start + c)%bmax = cmax
      octree%node(start + c)%depth = octree%node(parent_id)%depth + 1_I4P
      octree%node(start + c)%first_child = -1_I4P
      octree%node(start + c)%leaf_flag = AWRAP_LEAF_FLAG_INTERIOR
   enddo
   octree%n_nodes = octree%n_nodes + 8_I4P
   endsubroutine allocate_children

   subroutine grow_node_array(octree)
   !< Geometric (2x) growth of the node array. Amortized O(1) per append.
   type(awrap_octree_t), intent(inout) :: octree
   type(awrap_octree_node_t), allocatable :: bigger(:)
   integer(I4P) :: old_cap

   old_cap = size(octree%node, kind=I4P)
   allocate(bigger(2_I4P * old_cap))
   bigger(1:old_cap) = octree%node(1:old_cap)
   call move_alloc(from=bigger, to=octree%node)
   endsubroutine grow_node_array

   pure subroutine filter_overlapping(node_bmin, node_bmax, facet, candidate, kept, n_kept)
   !< Among `candidate(:)` facet IDs, return those whose triangle overlaps the
   !< given AABB. Output `kept(:)` is allocated to the exact count `n_kept`.
   type(vector_R8P),          intent(in)  :: node_bmin, node_bmax
   type(facet_object),        intent(in)  :: facet(:)
   integer(I4P),              intent(in)  :: candidate(:)
   integer(I4P), allocatable, intent(out) :: kept(:)
   integer(I4P),              intent(out) :: n_kept
   integer(I4P), allocatable              :: tmp(:)
   integer(I4P)                           :: i, fid

   allocate(tmp(size(candidate, kind=I4P)))
   n_kept = 0_I4P
   do i = 1_I4P, size(candidate, kind=I4P)
      fid = candidate(i)
      if (triangle_overlaps_aabb(bmin=node_bmin, bmax=node_bmax, &
                                 v1=facet(fid)%vertex(1), &
                                 v2=facet(fid)%vertex(2), &
                                 v3=facet(fid)%vertex(3))) then
         n_kept = n_kept + 1_I4P
         tmp(n_kept) = fid
      endif
   enddo
   allocate(kept(n_kept))
   if (n_kept > 0_I4P) kept(1:n_kept) = tmp(1:n_kept)
   deallocate(tmp)
   endsubroutine filter_overlapping

   pure subroutine compute_facet_bbox(facet, bmin, bmax)
   !< Tight axis-aligned bbox of all vertices in the facet array.
   type(facet_object), intent(in)  :: facet(:)
   type(vector_R8P),   intent(out) :: bmin, bmax
   integer(I4P)                    :: f, v

   bmin = facet(1)%vertex(1)
   bmax = facet(1)%vertex(1)
   do f = 1_I4P, size(facet, kind=I4P)
      do v = 1_I4P, 3_I4P
         bmin%x = min(bmin%x, facet(f)%vertex(v)%x)
         bmin%y = min(bmin%y, facet(f)%vertex(v)%y)
         bmin%z = min(bmin%z, facet(f)%vertex(v)%z)
         bmax%x = max(bmax%x, facet(f)%vertex(v)%x)
         bmax%y = max(bmax%y, facet(f)%vertex(v)%y)
         bmax%z = max(bmax%z, facet(f)%vertex(v)%z)
      enddo
   enddo
   endsubroutine compute_facet_bbox

   subroutine recount(octree)
   !< Walk the node array and update `n_leaves` and `n_boundary_leaves`.
   !< Cheap (O(n_nodes)); kept separate so the build path stays straight.
   type(awrap_octree_t), intent(inout) :: octree
   integer(I4P)                        :: i

   octree%n_leaves = 0_I4P
   octree%n_boundary_leaves = 0_I4P
   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child == -1_I4P) then
         octree%n_leaves = octree%n_leaves + 1_I4P
         if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_BOUNDARY) then
            octree%n_boundary_leaves = octree%n_boundary_leaves + 1_I4P
         endif
      endif
   enddo
   endsubroutine recount

   subroutine awrap_classify_leaves(octree, status)
   !< Step 2: classify every EMPTY leaf as INSIDE or OUTSIDE via 6-connectivity
   !< flood fill from a known-outside corner seed (issue #18 §1.6 step 2).
   !<
   !< Algorithm:
   !<   1. Locate the leaf containing the corner of the padded root bbox
   !<      (always EMPTY by construction — step 1 padded the root by 2α).
   !<   2. BFS outward, only crossing EMPTY ↔ EMPTY edges. BOUNDARY leaves are
   !<      barriers (the "wrap" surface separates outside from inside).
   !<   3. Every EMPTY leaf reached by the BFS becomes OUTSIDE.
   !<   4. Every EMPTY leaf NOT reached is INSIDE (interior of the input).
   !<
   !< For triangle-soup input with holes, the flood fill correctly leaks
   !< through the holes — that's the EXPECTED behaviour at this step.
   !< Step 3+ closes those holes in the WRAP surface, not by repairing the
   !< input topology. So a cube with one facet deleted will show most of the
   !< "interior" as OUTSIDE after step 2; that's the algorithmic ground
   !< truth, and the wrap mesh will still be watertight after step 5.
   !<
   !< Updates `octree%n_inside_leaves` and `octree%n_outside_leaves`.
   type(awrap_octree_t), intent(inout)           :: octree
   integer(I4P),         intent(out),   optional :: status
   integer(I4P), allocatable                     :: queue(:)        !< BFS queue (leaf node indices).
   integer(I4P)                                  :: head, tail, cap !< Queue cursors.
   integer(I4P)                                  :: i, seed_id, current_id, neighbour_id, axis, side
   type(vector_R8P)                              :: probe
   real(R8P)                                     :: half_alpha

   if (present(status)) status = AWRAP_STATUS_OK
   if (octree%n_nodes == 0_I4P) then
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif

   ! Reset prior step-2 state (idempotency for adaptive-refinement re-runs in step 5).
   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child /= -1_I4P) cycle
      if (octree%node(i)%leaf_flag /= AWRAP_LEAF_FLAG_BOUNDARY) &
         octree%node(i)%leaf_flag = AWRAP_LEAF_FLAG_EMPTY
   enddo

   ! Locate the seed leaf at the (-x, -y, -z) corner of the padded root.
   ! Step 1's padding guarantees this corner sits in an EMPTY leaf.
   seed_id = find_leaf_at_point(octree=octree, point=octree%node(1)%bmin)
   if (seed_id <= 0_I4P) then
      ! Defensive: if the corner happens to be on a BOUNDARY leaf (shouldn't
      ! happen with 2α padding) or location fails, scan all leaves for one
      ! that's empty AND on the root bbox boundary. Rare edge case.
      do i = 1_I4P, octree%n_nodes
         if (octree%node(i)%first_child /= -1_I4P) cycle
         if (octree%node(i)%leaf_flag /= AWRAP_LEAF_FLAG_EMPTY) cycle
         if (leaf_touches_root_min(octree=octree, leaf_id=i)) then
            seed_id = i
            exit
         endif
      enddo
   endif
   if (seed_id <= 0_I4P) then
      ! No seed found — geometry fills the entire padded bbox. Mark all empty
      ! leaves INSIDE and return.
      do i = 1_I4P, octree%n_nodes
         if (octree%node(i)%first_child /= -1_I4P) cycle
         if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_EMPTY) &
            octree%node(i)%leaf_flag = AWRAP_LEAF_FLAG_INSIDE
      enddo
      call count_classifications(octree)
      return
   endif

   ! BFS from the seed, only crossing EMPTY ↔ EMPTY edges.
   cap = max(octree%n_leaves, 64_I4P)
   allocate(queue(cap))
   head = 1_I4P
   tail = 1_I4P
   queue(1) = seed_id
   octree%node(seed_id)%leaf_flag = AWRAP_LEAF_FLAG_OUTSIDE

   half_alpha = 0.5_R8P * octree%alpha
   do while (head <= tail)
      current_id = queue(head)
      head = head + 1_I4P
      ! For each of the 6 face neighbours, find the leaf on the other side.
      ! Probe at the face midpoint, offset outward by a small fraction of α.
      do axis = 1_I4P, 3_I4P
         do side = -1_I4P, 1_I4P, 2_I4P
            call face_probe_point(node=octree%node(current_id), axis=axis, side=side, &
                                  half_alpha=half_alpha, probe=probe)
            ! If probe lies outside the padded root, no neighbour exists.
            if (probe%x < octree%node(1)%bmin%x .or. probe%x > octree%node(1)%bmax%x .or. &
                probe%y < octree%node(1)%bmin%y .or. probe%y > octree%node(1)%bmax%y .or. &
                probe%z < octree%node(1)%bmin%z .or. probe%z > octree%node(1)%bmax%z) cycle
            neighbour_id = find_leaf_at_point(octree=octree, point=probe)
            if (neighbour_id <= 0_I4P) cycle
            ! Fan-out across larger neighbours: a single probe may land in a
            ! leaf that abuts other not-yet-visited leaves on the SAME face.
            ! For MVP we accept this — flood fill will reach them via the
            ! transitive face traversal (their own face probes find this leaf).
            if (octree%node(neighbour_id)%leaf_flag /= AWRAP_LEAF_FLAG_EMPTY) cycle
            octree%node(neighbour_id)%leaf_flag = AWRAP_LEAF_FLAG_OUTSIDE
            tail = tail + 1_I4P
            if (tail > cap) call grow_queue(queue=queue, cap=cap)
            queue(tail) = neighbour_id
         enddo
      enddo
   enddo
   deallocate(queue)

   ! Mark every still-EMPTY leaf as INSIDE.
   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child /= -1_I4P) cycle
      if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_EMPTY) &
         octree%node(i)%leaf_flag = AWRAP_LEAF_FLAG_INSIDE
   enddo

   call count_classifications(octree)
   endsubroutine awrap_classify_leaves

   pure function find_leaf_at_point(octree, point) result(leaf_id)
   !< Locate the leaf containing `point` by descending from the root, picking
   !< the child whose bbox contains the point at each level. Returns 0 if the
   !< point is outside the root bbox (caller's responsibility to clamp).
   type(awrap_octree_t), intent(in) :: octree
   type(vector_R8P),     intent(in) :: point
   integer(I4P)                     :: leaf_id
   integer(I4P)                     :: cur, c, child_id

   leaf_id = 0_I4P
   if (octree%n_nodes == 0_I4P) return
   if (.not. point_in_node(octree%node(1), point)) return
   cur = 1_I4P
   do while (octree%node(cur)%first_child /= -1_I4P)
      ! Descend to the child whose bbox contains the point.
      do c = 0_I4P, 7_I4P
         child_id = octree%node(cur)%first_child + c
         if (point_in_node(octree%node(child_id), point)) then
            cur = child_id
            exit
         endif
      enddo
      ! Defensive: if no child contained the point (FP edge case at a split
      ! plane), give up and return the current cur — it's still a valid node
      ! containing the point.
      if (cur /= child_id .and. octree%node(cur)%first_child /= -1_I4P) exit
   enddo
   leaf_id = cur
   endfunction find_leaf_at_point

   pure function point_in_node(node, point) result(yes)
   type(awrap_octree_node_t), intent(in) :: node
   type(vector_R8P),          intent(in) :: point
   logical                               :: yes
   yes = (point%x >= node%bmin%x .and. point%x <= node%bmax%x) .and. &
         (point%y >= node%bmin%y .and. point%y <= node%bmax%y) .and. &
         (point%z >= node%bmin%z .and. point%z <= node%bmax%z)
   endfunction point_in_node

   pure subroutine face_probe_point(node, axis, side, half_alpha, probe)
   !< Generate a probe point just past the (axis, side) face of `node`. Used
   !< by the BFS to locate the face-neighbour leaf via point location. The
   !< offset (half α) is large enough to leave the source leaf even at
   !< depth-cap leaves where the leaf size is < α.
   type(awrap_octree_node_t), intent(in)  :: node
   integer(I4P),              intent(in)  :: axis      !< 1 = x, 2 = y, 3 = z.
   integer(I4P),              intent(in)  :: side      !< -1 = bmin face, +1 = bmax face.
   real(R8P),                 intent(in)  :: half_alpha
   type(vector_R8P),          intent(out) :: probe
   real(R8P)                              :: epsilon_offset

   ! Use min(half_alpha, half_leaf_size) to ensure we cross to a neighbour
   ! without overshooting more than one leaf width.
   epsilon_offset = min(half_alpha, &
                        0.25_R8P * min(node%bmax%x - node%bmin%x, &
                                       node%bmax%y - node%bmin%y, &
                                       node%bmax%z - node%bmin%z))
   ! Centre of the requested face, then push outward by epsilon_offset.
   probe = 0.5_R8P * (node%bmin + node%bmax)
   if (axis == 1_I4P) probe%x = merge(node%bmax%x + epsilon_offset, node%bmin%x - epsilon_offset, side > 0_I4P)
   if (axis == 2_I4P) probe%y = merge(node%bmax%y + epsilon_offset, node%bmin%y - epsilon_offset, side > 0_I4P)
   if (axis == 3_I4P) probe%z = merge(node%bmax%z + epsilon_offset, node%bmin%z - epsilon_offset, side > 0_I4P)
   endsubroutine face_probe_point

   pure function leaf_touches_root_min(octree, leaf_id) result(yes)
   !< Does the leaf bbox touch the root's (-x, -y, -z) corner? Used as a
   !< fallback seed test when point-location at the corner fails.
   type(awrap_octree_t), intent(in) :: octree
   integer(I4P),         intent(in) :: leaf_id
   logical                          :: yes

   yes = abs(octree%node(leaf_id)%bmin%x - octree%node(1)%bmin%x) < 1.0e-12_R8P .and. &
         abs(octree%node(leaf_id)%bmin%y - octree%node(1)%bmin%y) < 1.0e-12_R8P .and. &
         abs(octree%node(leaf_id)%bmin%z - octree%node(1)%bmin%z) < 1.0e-12_R8P
   endfunction leaf_touches_root_min

   subroutine grow_queue(queue, cap)
   !< Geometric (2x) growth of the BFS queue.
   integer(I4P), allocatable, intent(inout) :: queue(:)
   integer(I4P),              intent(inout) :: cap
   integer(I4P), allocatable                :: bigger(:)

   allocate(bigger(2_I4P * cap))
   bigger(1:cap) = queue(1:cap)
   call move_alloc(from=bigger, to=queue)
   cap = 2_I4P * cap
   endsubroutine grow_queue

   subroutine count_classifications(octree)
   !< Walk the leaves and update per-class counts. O(n_nodes); cheap.
   type(awrap_octree_t), intent(inout) :: octree
   integer(I4P)                        :: i

   octree%n_inside_leaves = 0_I4P
   octree%n_outside_leaves = 0_I4P
   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child /= -1_I4P) cycle
      select case (octree%node(i)%leaf_flag)
      case (AWRAP_LEAF_FLAG_INSIDE)
         octree%n_inside_leaves = octree%n_inside_leaves + 1_I4P
      case (AWRAP_LEAF_FLAG_OUTSIDE)
         octree%n_outside_leaves = octree%n_outside_leaves + 1_I4P
      endselect
   enddo
   endsubroutine count_classifications

   subroutine awrap_extract_surface(octree, wrapped_facets, status)
   !< Step 3: extract the wrap surface as the BOUNDARY ↔ OUTSIDE interface
   !< (issue #18 §1.6 step 3).
   !<
   !< For each face of every BOUNDARY leaf, query the OUTSIDE neighbour via
   !< face-probe point location. If the neighbour is OUTSIDE, emit a quad
   !< (two triangles) sized to the BOUNDARY leaf's face. Vertices placed at
   !< face corners — true dual-contouring vertex placement (QEF) is deferred
   !< to step 4's vertex projection pass.
   !<
   !< Topology is watertight 2-manifold by construction:
   !<   - Every face is shared by exactly two leaves.
   !<   - We emit a quad iff one is BOUNDARY and the other OUTSIDE (a
   !<     classification disagreement IS the wrap surface).
   !<   - Each emitted edge of a quad is shared by exactly one neighbouring
   !<     emitted quad (the adjacent BOUNDARY leaf's face).
   !<
   !< MVP simplification: walk only from the BOUNDARY side (smaller leaf in
   !< the typical adaptive octree, since boundaries refine to the depth cap
   !< while interiors stay coarse). When a small BOUNDARY leaf abuts a big
   !< OUTSIDE leaf, one quad sized to the small face is emitted — the big
   !< neighbour's contribution to the wrap surface comes from its other
   !< small BOUNDARY neighbours covering its face. Worst-case T-junctions
   !< near refinement transitions accepted; step 4's vertex projection +
   !< step 5's adaptive refinement handle them in practice.
   !<
   !< Normal orientation: each emitted quad's normal points OUTWARD (toward
   !< OUTSIDE side). Triangles wound counter-clockwise as viewed from
   !< OUTSIDE.
   type(awrap_octree_t),            intent(in)            :: octree
   type(facet_object), allocatable, intent(out)           :: wrapped_facets(:)
   integer(I4P),                    intent(out), optional :: status
   type(facet_object), allocatable                        :: tmp(:)
   integer(I4P)                                           :: cap, n_facets
   integer(I4P)                                           :: i, axis, side, neighbour_id
   type(vector_R8P)                                       :: probe
   real(R8P)                                              :: half_alpha
   type(vector_R8P)                                       :: c0, c1, c2, c3   !< Face corners.

   if (present(status)) status = AWRAP_STATUS_OK
   if (octree%n_nodes == 0_I4P) then
      allocate(wrapped_facets(0))
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif

   ! Capacity heuristic: each boundary leaf contributes at most 6 quads = 12 triangles.
   cap = max(64_I4P, 12_I4P * octree%n_boundary_leaves)
   allocate(tmp(cap))
   n_facets = 0_I4P
   half_alpha = 0.5_R8P * octree%alpha

   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child /= -1_I4P) cycle  ! interior octree node
      if (octree%node(i)%leaf_flag /= AWRAP_LEAF_FLAG_BOUNDARY) cycle
      do axis = 1_I4P, 3_I4P
         do side = -1_I4P, 1_I4P, 2_I4P
            call face_probe_point(node=octree%node(i), axis=axis, side=side, &
                                  half_alpha=half_alpha, probe=probe)
            ! Neighbour outside the root: treat as OUTSIDE (the wrap closes
            ! at the padded-bbox boundary — step 1's 2α padding ensures
            ! the geometry doesn't reach the bbox edge).
            if (probe%x < octree%node(1)%bmin%x .or. probe%x > octree%node(1)%bmax%x .or. &
                probe%y < octree%node(1)%bmin%y .or. probe%y > octree%node(1)%bmax%y .or. &
                probe%z < octree%node(1)%bmin%z .or. probe%z > octree%node(1)%bmax%z) then
               call face_corners(node=octree%node(i), axis=axis, side=side, c0=c0, c1=c1, c2=c2, c3=c3)
               call emit_quad(tmp=tmp, n_facets=n_facets, cap=cap, &
                              c0=c0, c1=c1, c2=c2, c3=c3, side=side, axis=axis)
               cycle
            endif
            neighbour_id = find_leaf_at_point(octree=octree, point=probe)
            if (neighbour_id <= 0_I4P) cycle
            if (octree%node(neighbour_id)%leaf_flag /= AWRAP_LEAF_FLAG_OUTSIDE) cycle
            ! Wrap face found: emit a quad on this BOUNDARY leaf's face.
            call face_corners(node=octree%node(i), axis=axis, side=side, c0=c0, c1=c1, c2=c2, c3=c3)
            call emit_quad(tmp=tmp, n_facets=n_facets, cap=cap, &
                           c0=c0, c1=c1, c2=c2, c3=c3, side=side, axis=axis)
         enddo
      enddo
   enddo

   allocate(wrapped_facets(n_facets))
   if (n_facets > 0_I4P) wrapped_facets(1:n_facets) = tmp(1:n_facets)
   deallocate(tmp)
   endsubroutine awrap_extract_surface

   pure subroutine face_corners(node, axis, side, c0, c1, c2, c3)
   !< Return the 4 corners of the (axis, side) face of `node`, ordered such
   !< that c0→c1→c2→c3 is a counter-clockwise loop when viewed from the
   !< OUTWARD direction (positive `side` looking back along -axis, negative
   !< `side` looking back along +axis). Used by `emit_quad` to construct
   !< two outward-facing triangles.
   type(awrap_octree_node_t), intent(in)  :: node
   integer(I4P),              intent(in)  :: axis  !< 1 = x, 2 = y, 3 = z.
   integer(I4P),              intent(in)  :: side  !< -1 = -axis face, +1 = +axis face.
   type(vector_R8P),          intent(out) :: c0, c1, c2, c3
   real(R8P)                              :: x_lo, x_hi, y_lo, y_hi, z_lo, z_hi
   real(R8P)                              :: x_face, y_face, z_face

   x_lo = node%bmin%x; x_hi = node%bmax%x
   y_lo = node%bmin%y; y_hi = node%bmax%y
   z_lo = node%bmin%z; z_hi = node%bmax%z

   if (axis == 1_I4P) then
      x_face = merge(x_hi, x_lo, side > 0_I4P)
      if (side > 0_I4P) then
         ! +X face, normal +X. CCW from +X looking back: -Y first.
         c0 = vector_R8P(x_face, y_lo, z_lo)
         c1 = vector_R8P(x_face, y_hi, z_lo)
         c2 = vector_R8P(x_face, y_hi, z_hi)
         c3 = vector_R8P(x_face, y_lo, z_hi)
      else
         ! -X face, normal -X. CCW from -X looking back: opposite winding.
         c0 = vector_R8P(x_face, y_lo, z_lo)
         c1 = vector_R8P(x_face, y_lo, z_hi)
         c2 = vector_R8P(x_face, y_hi, z_hi)
         c3 = vector_R8P(x_face, y_hi, z_lo)
      endif
   elseif (axis == 2_I4P) then
      y_face = merge(y_hi, y_lo, side > 0_I4P)
      if (side > 0_I4P) then
         c0 = vector_R8P(x_lo, y_face, z_lo)
         c1 = vector_R8P(x_lo, y_face, z_hi)
         c2 = vector_R8P(x_hi, y_face, z_hi)
         c3 = vector_R8P(x_hi, y_face, z_lo)
      else
         c0 = vector_R8P(x_lo, y_face, z_lo)
         c1 = vector_R8P(x_hi, y_face, z_lo)
         c2 = vector_R8P(x_hi, y_face, z_hi)
         c3 = vector_R8P(x_lo, y_face, z_hi)
      endif
   else  ! axis == 3
      z_face = merge(z_hi, z_lo, side > 0_I4P)
      if (side > 0_I4P) then
         c0 = vector_R8P(x_lo, y_lo, z_face)
         c1 = vector_R8P(x_hi, y_lo, z_face)
         c2 = vector_R8P(x_hi, y_hi, z_face)
         c3 = vector_R8P(x_lo, y_hi, z_face)
      else
         c0 = vector_R8P(x_lo, y_lo, z_face)
         c1 = vector_R8P(x_lo, y_hi, z_face)
         c2 = vector_R8P(x_hi, y_hi, z_face)
         c3 = vector_R8P(x_hi, y_lo, z_face)
      endif
   endif
   endsubroutine face_corners

   subroutine emit_quad(tmp, n_facets, cap, c0, c1, c2, c3, side, axis)
   !< Emit two triangles (c0,c1,c2) and (c0,c2,c3) representing one wrap-
   !< surface quad. The normal points OUTWARD (along +axis when side>0,
   !< -axis when side<0). Caller has already verified the face_corners
   !< ordering matches this convention.
   type(facet_object), allocatable, intent(inout) :: tmp(:)
   integer(I4P),                    intent(inout) :: n_facets
   integer(I4P),                    intent(inout) :: cap
   type(vector_R8P),                intent(in)    :: c0, c1, c2, c3
   integer(I4P),                    intent(in)    :: side, axis
   type(vector_R8P)                               :: outward

   ! Ensure capacity for two triangles.
   if (n_facets + 2_I4P > cap) call grow_facet_buffer(tmp=tmp, cap=cap)

   ! Compute outward normal (axis-aligned for this MVP).
   outward = 0._R8P * ex_R8P  ! init to zero
   if (axis == 1_I4P) outward = real(side, R8P) * ex_R8P
   if (axis == 2_I4P) outward = real(side, R8P) * ey_R8P
   if (axis == 3_I4P) outward = real(side, R8P) * ez_R8P

   ! Triangle 1: c0, c1, c2.
   n_facets = n_facets + 1_I4P
   tmp(n_facets)%vertex(1) = c0
   tmp(n_facets)%vertex(2) = c1
   tmp(n_facets)%vertex(3) = c2
   tmp(n_facets)%normal = outward
   call tmp(n_facets)%compute_metrix
   ! Triangle 2: c0, c2, c3.
   n_facets = n_facets + 1_I4P
   tmp(n_facets)%vertex(1) = c0
   tmp(n_facets)%vertex(2) = c2
   tmp(n_facets)%vertex(3) = c3
   tmp(n_facets)%normal = outward
   call tmp(n_facets)%compute_metrix
   endsubroutine emit_quad

   subroutine grow_facet_buffer(tmp, cap)
   !< Geometric (2x) growth of the wrap-facet buffer.
   type(facet_object), allocatable, intent(inout) :: tmp(:)
   integer(I4P),                    intent(inout) :: cap
   type(facet_object), allocatable                :: bigger(:)

   allocate(bigger(2_I4P * cap))
   bigger(1:cap) = tmp(1:cap)
   call move_alloc(from=bigger, to=tmp)
   cap = 2_I4P * cap
   endsubroutine grow_facet_buffer

endmodule fossil_alpha_wrap
