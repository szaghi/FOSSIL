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
use vecfor, only : vector_R8P

implicit none
private

public :: awrap_octree_t
public :: awrap_build_octree
public :: AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE
public :: AWRAP_LEAF_FLAG_BOUNDARY, AWRAP_LEAF_FLAG_INTERIOR
public :: AWRAP_MAX_DEPTH

integer(I4P), parameter :: AWRAP_STATUS_OK         = 0_I4P
integer(I4P), parameter :: AWRAP_STATUS_BAD_INPUT  = 1_I4P  !< Empty facet list, or alpha <= 0.
integer(I4P), parameter :: AWRAP_STATUS_DEGENERATE = 2_I4P  !< Alpha > bbox diagonal: octree won't refine, output empty.

integer(I4P), parameter :: AWRAP_LEAF_FLAG_INTERIOR = 0_I4P  !< Leaf has no overlapping input facets.
integer(I4P), parameter :: AWRAP_LEAF_FLAG_BOUNDARY = 1_I4P  !< Leaf overlaps at least one input facet.

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

endmodule fossil_alpha_wrap
