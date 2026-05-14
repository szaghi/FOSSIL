!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

module fossil_aabb_tree_object
!< FOSSIL, Axis-Aligned Bounding Box (AABB) tree class definition.
!<
!< @note The tree is assumed to be an **octree**.

use fossil_aabb_object, only : aabb_object
use fossil_aabb_node_object, only : aabb_node_object
use fossil_facet_object, only : facet_object, triangle_point_distance
use fossil_list_id_object, only : list_id_object
use fossil_ray_query, only : ray_hit_t
use penf, only : I4P, R8P, MaxR8P, str
use vecfor, only : vector_R8P

implicit none
private
public :: aabb_tree_object
public :: AABB_USE_INDEX, AABB_USE_BRUTE_FORCE
public :: AABB_AUTO_REFINEMENT
public :: AABB_TREE_OCTREE, AABB_TREE_SAH_BVH

integer(I4P), parameter :: TREE_RATIO=8 !< Tree refinement ratio, it is assumed to be an **octree**.

! Flattened distance-only facet payload for the SAH BVH (issue #19 §B1).
!
! The distance traversal's innermost kernel needs only the triangle metrix
! `[v1, e12, e13, a, b, c, det]` — 13 reals out of a ~400-byte `facet_object`.
! Carrying that as a contiguous, leaf-grouped SoA removes the
! `node -> aabb -> facet_id%id -> facet(fid)` indirection chain (three pointer
! hops + a random-stride gather) and turns each leaf's facet scan into a
! stride-1 sweep. `facet_id` keeps the global facet index so the winning facet
! can be mapped back for the pseudo-normal sign step.
!
! The payload is built once after the BVH is built (`build_distance_payload`),
! by walking leaves in node order so each leaf's facets land contiguously; the
! leaf records its `(payload_first, payload_count)` slice. Only the SAH BVH
! populates it — the octree path is untouched and still scans via `facet_id`.
type :: facet_distance_payload
   !< Packed per-facet triangle metrix for the BVH distance traversal.
   real(R8P)    :: v1(3)   = 0._R8P !< Triangle first vertex.
   real(R8P)    :: e12(3)  = 0._R8P !< Edge 1-2, `v2 - v1`.
   real(R8P)    :: e13(3)  = 0._R8P !< Edge 1-3, `v3 - v1`.
   real(R8P)    :: a       = 0._R8P !< `e12 . e12`.
   real(R8P)    :: b       = 0._R8P !< `e12 . e13`.
   real(R8P)    :: c       = 0._R8P !< `e13 . e13`.
   real(R8P)    :: det     = 0._R8P !< Gram determinant `a*c - b*b`.
   integer(I4P) :: facet_id = 0_I4P !< Global facet index (maps the winner back for sign determination).
endtype facet_distance_payload

! dispatch knob exposed through `use_index`/`set_use_index`: lets callers force the brute-force
! distance/ray-intersection path for benchmarking without lying about `is_initialized`.
logical, parameter :: AABB_USE_INDEX       = .true.  !< Use the AABB octree for queries.
logical, parameter :: AABB_USE_BRUTE_FORCE = .false. !< Force brute-force scan over all facets.

! Tree-kind selector. The default `AABB_TREE_SAH_BVH` is a binary BVH built by
! triangle-partition with the surface-area heuristic — adapts to triangle density
! rather than spatial uniformity. Empirical benchmark (commit "switch default to
! SAH BVH"): on dragon-fine (24k facets) at a 32^3 query grid, the BVH was ~110x
! faster than the octree at the same query, with ~40% faster build time too.
! `AABB_TREE_OCTREE` selects the original 8-way space-partitioning octree, kept
! for benchmarking and as a fallback for any user who needs the legacy behaviour.
integer(I4P), parameter :: AABB_TREE_OCTREE  = 0_I4P
integer(I4P), parameter :: AABB_TREE_SAH_BVH = 1_I4P

! Sentinel for auto-tuned refinement: pass this where an explicit depth would go,
! and `initialize` picks the depth from the facet count via `auto_refinement_levels`.
! Only meaningful for the octree path; the SAH BVH self-tunes depth from cost.
integer(I4P), parameter :: AABB_AUTO_REFINEMENT = -1_I4P
! Heuristic parameters (private; the literature suggests 16-32 facets per leaf is
! a broad sweet spot for triangle-mesh BVHs, and FOSSIL's flat empirical curve
! tolerates any choice in that range to within 5%).
integer(I4P), parameter :: AUTO_TARGET_FACETS_PER_LEAF = 64_I4P  !< Heuristic facet-per-leaf target.
integer(I4P), parameter :: AUTO_MIN_LEVELS             = 1_I4P   !< Floor: at least one split.
integer(I4P), parameter :: AUTO_MAX_LEVELS             = 6_I4P   !< Cap: deeper octrees waste space on 2D-in-3D meshes.

type :: ofsm
   !< Octree Finite State Machine class for efficient searching of neiighbors.
   integer(I4P) :: octant=0    !< Octant ID, [1,2,3,4,5,6,7,8].
   integer(I4P) :: direction=0 !< Direction, 0=Halt, 1=X+, 2=X-, 3=Y+, 4=Y-, 5=Z+, 6=Z-.
endtype ofsm

type(ofsm) :: octree_fsm(6,8) = reshape(                                                                   &
                                        [ofsm(2,0), ofsm(2,2), ofsm(3,0), ofsm(3,4), ofsm(5,0), ofsm(5,6), &
                                         ofsm(1,1), ofsm(1,0), ofsm(4,0), ofsm(4,4), ofsm(6,0), ofsm(6,6), &
                                         ofsm(4,0), ofsm(4,2), ofsm(1,3), ofsm(1,0), ofsm(7,0), ofsm(7,6), &
                                         ofsm(3,1), ofsm(3,0), ofsm(2,3), ofsm(2,0), ofsm(8,0), ofsm(8,6), &
                                         ofsm(6,0), ofsm(6,2), ofsm(7,0), ofsm(7,4), ofsm(1,5), ofsm(1,0), &
                                         ofsm(5,1), ofsm(5,0), ofsm(8,0), ofsm(8,4), ofsm(2,5), ofsm(2,0), &
                                         ofsm(8,0), ofsm(8,2), ofsm(5,3), ofsm(5,0), ofsm(3,5), ofsm(3,0), &
                                         ofsm(7,1), ofsm(7,0), ofsm(6,3), ofsm(6,0), ofsm(4,5), ofsm(4,0)], [6,8])

type :: aabb_tree_object
   !< FOSSIL Axis-Aligned Bounding Box (AABB) tree class.
   !<
   !< @note The tree is assumed to be an **octree**.
   !< The octree uses a breath-leaf counting, with the following convention:
   !<
   !<```
   !>                  +----+----+
   !>                 /|   /|   /|
   !>                / |  7 |  8 |
   !>               /  +----*----+
   !>              /  /|/  /|/  /|
   !>             /  / |  5 |  6 |
   !>            /  / /+----+----+
   !>           /  / // / // /  /
   !>          +----+----+/ /  /
   !>         /| / /| / /| /  /
   !>        / |/ //|/ //|/  /
   !>       /  +----*----+  /
   !>      /  /|// /|// /| /
   !>     /  / |/ / |/ / |/
   !>    /  / /+----+----+
   !>   /  / // / // /  /
   !>  +----+----+/ /  /       y(j)  z(k)
   !>  | /  | /  | /  /        ^    ^
   !>  |/ 3/|/ 4/|/  /         |   /
   !>  +----*----+  /          |  /
   !>  | /  | /  | /           | /
   !>  |/ 1 |/ 2 |/            |/
   !>  +----+----+             +------->x(i)
   !<```
   private
   integer(I4P)                        :: tree_kind=AABB_TREE_SAH_BVH            !< Selects between the SAH BVH (default) and the legacy octree.
   integer(I4P)                        :: refinement_levels=AABB_AUTO_REFINEMENT !< Octree depth (AABB_AUTO_REFINEMENT = auto-tune); unused by SAH BVH.
   integer(I4P)                        :: nodes_number=0         !< Total number of tree nodes.
   type(aabb_node_object), allocatable :: node(:)                !< AABB tree nodes [0:nodes_number-1].
   type(facet_distance_payload), allocatable :: payload(:)       !< Flattened distance-only facet SoA, leaf-grouped (SAH BVH only, issue #19 §B1).
   logical                             :: is_initialized=.false. !< Sentinel to check is AABB tree is initialized.
   logical                             :: use_index=.true.       !< Dispatch knob: .true.=use index tree, .false.=brute-force scan.
   contains
      ! read-only accessors (pure, inlined at -O2, zero data copy)
      procedure, pass(self) :: get_refinement_levels !< Return refinement_levels.
      procedure, pass(self) :: get_tree_kind         !< Return tree_kind (octree vs SAH BVH).
      procedure, pass(self) :: get_nodes_number      !< Return nodes_number.
      procedure, pass(self) :: node_at               !< Return pointer to node(i); null() if out of range.
      procedure, pass(self) :: get_is_initialized    !< Return is_initialized.
      procedure, pass(self) :: get_use_index         !< Return use_index dispatch knob.
      ! mutators
      procedure, pass(self) :: set_refinement_levels !< Set refinement_levels (resets initialization).
      procedure, pass(self) :: set_tree_kind         !< Set tree_kind (resets initialization).
      procedure, pass(self) :: set_use_index         !< Set use_index dispatch knob.
      ! public methods
      procedure, pass(self) :: compute_vertices_nearby     !< Compute vertices nearby.
      procedure, pass(self) :: destroy                     !< Destroy AABB tree.
      procedure, pass(self) :: distance                    !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: distance_tree               !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: distance_tree_with_region   !< Return distance + closest facet id + Voronoi region (for signed distance).
      procedure, pass(self) :: distribute_facets           !< Distribute facets into AABB nodes.
      procedure, pass(self) :: distribute_facets_tree      !< Distribute facets into AABB nodes.
      procedure, pass(self) :: has_children                !< Return true if node has at least one child allocated.
      procedure, pass(self) :: initialize                  !< Initialize AABB tree.
      procedure, pass(self) :: loop_node                   !< Loop over all nodes.
      procedure, pass(self) :: ray_intersections_number    !< Return ray intersections number.
      procedure, pass(self) :: intersect_ray_all_tree      !< Tree-accelerated all-hits ray query (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_first_tree    !< Tree-accelerated first-hit ray query (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_any_tree      !< Tree-accelerated any-hit ray query with early exit (issue #18 §2.5).
      procedure, pass(self) :: save_geometry_tecplot_ascii !< Save AABB tree boxes geometry into Tecplot ascii file.
      procedure, pass(self) :: translate                   !< Translate AABB tree by delta.
      ! operators
      generic :: assignment(=) => aabb_tree_assign_aabb_tree      !< Overload `=`.
      procedure, pass(lhs), private :: aabb_tree_assign_aabb_tree !< Operator `=`.
      ! finaliser (releases node array even when wrapped in arrays-of-trees)
      final :: aabb_tree_finalize
      ! private methods
      procedure, pass(self), private :: build_bvh_sah                 !< Build a SAH BVH over the given facet list.
      procedure, pass(self), private :: build_distance_payload        !< Flatten BVH leaf facets into the contiguous distance SoA (issue #19 §B1).
      procedure, pass(self) :: enumerate_children                    !< List the allocated children of a node (octree or BVH).
      procedure, pass(self), private :: distance_node                 !< Update best squared distance by recursing into a node.
      procedure, pass(self), private :: distance_node_with_region     !< Update (best d^2, best facet id, best region) by recursing into a node.
      procedure, pass(self), private :: scan_payload                  !< Stride-1 unsigned leaf scan over the distance payload (issue #19 §B1).
      procedure, pass(self), private :: scan_payload_with_facet       !< Stride-1 signed leaf scan over the distance payload (issue #19 §B1).
      procedure, pass(self), private :: ray_intersections_number_node !< Return ray intersections number into a node of AABB tree.
      procedure, pass(self), private :: intersect_ray_all_node        !< Recursive driver behind intersect_ray_all_tree.
      procedure, pass(self), private :: intersect_ray_first_node      !< Recursive driver behind intersect_ray_first_tree.
      procedure, pass(self), private :: intersect_ray_any_node        !< Recursive driver behind intersect_ray_any_tree.
endtype aabb_tree_object

contains
   pure function auto_refinement_levels(facets_number) result(levels)
   !< Pick a refinement depth that targets `AUTO_TARGET_FACETS_PER_LEAF` facets per
   !< leaf in a fully-populated octree.
   !<
   !< Formula: `levels = ceil(log8(N / target))`, clamped to [AUTO_MIN_LEVELS, AUTO_MAX_LEVELS].
   !< Computed in integer arithmetic (no logs) by walking 8^k until the product
   !< exceeds N/target.
   !<
   !< Why these defaults: empirical sweeps on cube/naca/dragon/dragon-fine showed
   !< a near-flat timing curve vs depth (within 5%) once the best-first traversal
   !< from step 1 was in place. AUTO_TARGET_FACETS_PER_LEAF=64 places the formula
   !< right around the empirical sweet spot for every size we measured. The cap at
   !< AUTO_MAX_LEVELS=6 prevents pathologically deep trees on huge meshes; the
   !< floor at 1 keeps the smallest trees non-degenerate.
   integer(I4P), intent(in) :: facets_number  !< Mesh facet count.
   integer(I4P)             :: levels         !< Auto-picked refinement depth.
   integer(I4P)             :: cap, k

   cap = max(facets_number, 1_I4P) / AUTO_TARGET_FACETS_PER_LEAF
   levels = 0
   k      = 1
   do while (k < cap)
      k      = k * TREE_RATIO
      levels = levels + 1
      if (levels >= AUTO_MAX_LEVELS) exit
   enddo
   if (levels < AUTO_MIN_LEVELS) levels = AUTO_MIN_LEVELS
   endfunction auto_refinement_levels

   ! accessors (pure, scalar return — inlined by gfortran/ifort at -O2, no data copy)
   pure function get_refinement_levels(self) result(n)
   !< Return refinement_levels.
   class(aabb_tree_object), intent(in) :: self !< AABB tree.
   integer(I4P)                        :: n    !< Refinement levels.

   n = self%refinement_levels
   endfunction get_refinement_levels

   pure function get_tree_kind(self) result(kind)
   !< Return tree_kind (`AABB_TREE_OCTREE` or `AABB_TREE_SAH_BVH`).
   class(aabb_tree_object), intent(in) :: self !< AABB tree.
   integer(I4P)                        :: kind !< Tree-kind selector.

   kind = self%tree_kind
   endfunction get_tree_kind

   pure function get_nodes_number(self) result(n)
   !< Return nodes_number.
   class(aabb_tree_object), intent(in) :: self !< AABB tree.
   integer(I4P)                        :: n    !< Nodes number.

   n = self%nodes_number
   endfunction get_nodes_number

   function node_at(self, i) result(p)
   !< Return a pointer to `node(i)` (range `[0, nodes_number-1]`), or `null()` if `i`
   !< is out of range. Intended as a read-only view: external code should treat the
   !< pointer as const. Used by tests and diagnostic code that walks the tree.
   class(aabb_tree_object), target, intent(in) :: self !< AABB tree.
   integer(I4P),                    intent(in) :: i    !< Node index.
   type(aabb_node_object), pointer             :: p    !< Pointer to node(i), or null().

   p => null()
   if (.not. allocated(self%node)) return
   if (i < 0 .or. i > self%nodes_number - 1) return
   p => self%node(i)
   endfunction node_at

   pure function get_is_initialized(self) result(yes)
   !< Return is_initialized.
   class(aabb_tree_object), intent(in) :: self !< AABB tree.
   logical                             :: yes  !< Initialization status.

   yes = self%is_initialized
   endfunction get_is_initialized

   pure function get_use_index(self) result(yes)
   !< Return use_index dispatch knob.
   class(aabb_tree_object), intent(in) :: self !< AABB tree.
   logical                             :: yes  !< .true. if octree dispatch is enabled.

   yes = self%use_index
   endfunction get_use_index

   ! mutators (the only externally-permitted writes)
   elemental subroutine set_refinement_levels(self, refinement_levels)
   !< Set refinement_levels.
   !<
   !< @note This does not (re)build the tree. The caller must call `initialize` afterwards if a
   !< rebuild is desired. We mark the tree as not-initialized to keep the invariant honest.
   class(aabb_tree_object), intent(inout) :: self              !< AABB tree.
   integer(I4P),            intent(in)    :: refinement_levels !< Refinement levels.

   self%refinement_levels = refinement_levels
   self%is_initialized    = .false.
   endsubroutine set_refinement_levels

   elemental subroutine set_tree_kind(self, tree_kind)
   !< Set the tree_kind selector. Like `set_refinement_levels`, this invalidates any existing
   !< build state; the caller must `initialize` afterwards. Out-of-range values raise `error stop`.
   class(aabb_tree_object), intent(inout) :: self      !< AABB tree.
   integer(I4P),            intent(in)    :: tree_kind !< `AABB_TREE_OCTREE` or `AABB_TREE_SAH_BVH`.

   select case (tree_kind)
   case (AABB_TREE_OCTREE, AABB_TREE_SAH_BVH)
      self%tree_kind = tree_kind
   case default
      error stop 'aabb_tree_object%set_tree_kind: unknown tree_kind (valid: AABB_TREE_OCTREE=0, AABB_TREE_SAH_BVH=1)'
   end select
   self%is_initialized = .false.
   endsubroutine set_tree_kind

   elemental subroutine set_use_index(self, use_index)
   !< Set use_index dispatch knob.
   !<
   !< Use `AABB_USE_INDEX` / `AABB_USE_BRUTE_FORCE` for readability. This does not change the
   !< tree's allocation state — it only affects which dispatch path distance/ray queries take.
   class(aabb_tree_object), intent(inout) :: self      !< AABB tree.
   logical,                 intent(in)    :: use_index !< .true.=octree, .false.=brute-force.

   self%use_index = use_index
   endsubroutine set_use_index

   ! finaliser
   subroutine aabb_tree_finalize(self)
   !< Release node storage and reset state.
   type(aabb_tree_object), intent(inout) :: self !< AABB tree.

   if (allocated(self%node))    deallocate(self%node)
   if (allocated(self%payload)) deallocate(self%payload)
   self%nodes_number   = 0
   self%is_initialized = .false.
   endsubroutine aabb_tree_finalize

   ! public methods
   pure subroutine compute_vertices_nearby(self, facet, tolerance_to_be_nearby)
   !< Compute vertices nearby (loose tolerance only; see facet variant).
   class(aabb_tree_object), intent(in)    :: self                   !< AABB tree.
   type(facet_object),      intent(inout) :: facet(:)               !< Facets list.
   real(R8P),               intent(in)    :: tolerance_to_be_nearby !< Tolerance to identify nearby vertices.
   integer(I4P)                           :: level                  !< Counter.
   integer(I4P)                           :: b, bb, bbb             !< Counter.

   if (self%nodes_number > 0) then
      level=self%refinement_levels                                      ! check only max refinement level
      b = first_node(level=level)                                       ! first node at level
      do bb=1, nodes_number_at_level(level=level)                       ! loop over nodes at level
         bbb = b + bb - 1                                               ! node numeration in tree
         call self%node(bbb)%compute_vertices_nearby(facet=facet, tolerance_to_be_nearby=tolerance_to_be_nearby)
      enddo
   endif
   endsubroutine compute_vertices_nearby

   elemental subroutine destroy(self)
   !< Destroy AABB tree.
   class(aabb_tree_object), intent(inout) :: self  !< AABB tree.
   type(aabb_tree_object)                 :: fresh !< Fresh instance of AABB tree.

   self = fresh
   endsubroutine destroy

   pure function distance(self, facet, point)
   !< Compute the (minimum) distance from a point to the triangulated surface.
   class(aabb_tree_object), intent(in) :: self            !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)        !< Facets list.
   type(vector_R8P),        intent(in) :: point           !< Point coordinates.
   real(R8P)                           :: distance        !< Minimum distance from point to the triangulated surface.
   real(R8P), allocatable              :: distance_(:)    !< Minimum distance, temporary buffer.
   integer(I4P), allocatable           :: aabb_closest(:) !< Index of closest AABB.
   integer(I4P)                        :: level           !< Counter.
   integer(I4P)                        :: b, bb, bbb      !< Counter.

   associate(node=>self%node)
      allocate(distance_(0:self%refinement_levels))
      allocate(aabb_closest(0:self%refinement_levels))
      distance_ = MaxR8P
      aabb_closest = -1
      do level=0, self%refinement_levels                  ! loop over refinement levels
         b = first_node(level=level)                      ! first node at finest level
         do bb=1, nodes_number_at_level(level=level)      ! loop over nodes at level
            bbb = b + bb - 1                              ! node numeration in tree
            if (node(bbb)%is_allocated()) then
               distance = node(bbb)%distance(point=point) ! node distance
               if (distance <= distance_(level)) then
                  distance_(level) = distance             ! update minimum distance
                  aabb_closest(level) = bbb               ! store closest node
               endif
            endif
         enddo
      enddo
      distance = MaxR8P
      do level=0, self%refinement_levels
         if (aabb_closest(level) >= 0) then
            distance = min(distance, node(aabb_closest(level))%distance_from_facets(facet=facet, point=point))
         endif
      enddo
   endassociate
   endfunction distance

   function distance_tree(self, facet, point) result(distance)
   !< Compute the (minimum) distance from a point to the triangulated surface.
   !<
   !< Best-first descent: maintain a running best squared distance and recurse only
   !< into children whose box squared distance is below the current best. All
   !< arithmetic is in squared distance — `aabb%distance` and `facet%compute_distance`
   !< both return d^2, so the comparison is unit-consistent.
   class(aabb_tree_object), intent(in) :: self         !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)     !< Facets list.
   type(vector_R8P),        intent(in) :: point        !< Point coordinates.
   real(R8P)                           :: distance     !< Minimum squared distance.

   distance = MaxR8P
   call self%distance_node(n=0, facet=facet, point=point, best=distance)
   endfunction distance_tree

   subroutine distance_tree_with_region(self, facet, point, distance, closest_facet, closest_region)
   !< Best-first traversal that also reports the closest facet id and Voronoi region.
   !<
   !< Used by the pseudo-normal signed-distance path: caller forms
   !< `sign( (point - closest_point) . pseudo_normal_of(closest_facet, closest_region) )`.
   !< Same pruning logic as `distance_tree`, but each leaf update also records which
   !< facet and which region of that facet were the new best.
   class(aabb_tree_object), intent(in)  :: self            !< AABB tree.
   type(facet_object),      intent(in)  :: facet(:)        !< Facets list.
   type(vector_R8P),        intent(in)  :: point           !< Point coordinates.
   real(R8P),               intent(out) :: distance        !< Minimum squared distance.
   integer(I4P),            intent(out) :: closest_facet   !< Facet id holding the closest point (0 if none found).
   integer(I4P),            intent(out) :: closest_region  !< Voronoi region tag of the closest point.

   distance       = MaxR8P
   closest_facet  = 0_I4P
   closest_region = 0_I4P
   call self%distance_node_with_region(n=0, facet=facet, point=point, &
                                       best=distance, best_facet=closest_facet, best_region=closest_region)
   endsubroutine distance_tree_with_region

   pure subroutine distribute_facets(self, facet, is_exclusive, do_update_extents)
   !< Distribute facets into AABB nodes.
   class(aabb_tree_object), intent(inout)        :: self               !< AABB tree.
   type(facet_object),      intent(in)           :: facet(:)           !< Facets list.
   logical,                 intent(in), optional :: is_exclusive       !< Sentinel to enable/disable exclusive addition.
   logical,                 intent(in), optional :: do_update_extents  !< Sentinel to enable/disable AABB extents update.
   logical                                       :: do_update_extents_ !< Sentinel to enable/disable AABB extents update, local var.
   type(list_id_object)                          :: facet_id           !< List of facets IDs.
   integer(I4P)                                  :: level              !< Counter.
   integer(I4P)                                  :: b, bb, bbb         !< Counter.

   do_update_extents_ = .true. ; if (present(do_update_extents)) do_update_extents_ = do_update_extents
   associate(node=>self%node)
      call facet_id%initialize(id=facet%id)

      ! add facets to nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (facet_id%ids_number > 0) then                       ! check if facets list still has facets
               call node(bbb)%add_facets(facet_id=facet_id, &
                                         facet=facet,       &
                                         is_exclusive=is_exclusive) ! add facets to node and prune added facets from list
            endif
         enddo
      enddo

      ! destroy void nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (.not.node(bbb)%has_facets()) call node(bbb)%destroy ! destroy void node
         enddo
      enddo

      ! update AABB extents
      if (do_update_extents_) then
         do level=self%refinement_levels, 0, -1           ! loop over refinement levels
            b = first_node(level=level)                   ! first node at level
            do bb=1, nodes_number_at_level(level=level)   ! loop over nodes at level
               bbb = b + bb - 1                           ! node numeration in tree
               call node(bbb)%update_extents(facet=facet) ! update extents
            enddo
         enddo
      endif
   endassociate
   endsubroutine distribute_facets

   pure subroutine distribute_facets_tree(self, facet)
   !< Distribute facets into AABB nodes.
   class(aabb_tree_object), intent(inout) :: self       !< AABB tree.
   type(facet_object),      intent(in)    :: facet(:)   !< Facets list.
   type(list_id_object)                   :: facet_id   !< List of facets IDs.
   integer(I4P)                           :: level      !< Counter.
   integer(I4P)                           :: b, bb, bbb !< Counter.
   integer(I4P)                           :: parent     !< Parent node index.

   associate(node=>self%node)
      ! add facets to nodes
      call facet_id%initialize(id=facet%id)                             ! initialize facets IDs list
      if (facet_id%ids_number > 0) then                                 ! check if facets list still has facets
         call node(0)%add_facets(facet_id=facet_id, &
                                 facet=facet,       &
                                 is_exclusive=.false.)                  ! add facets to root node
         do level=1, self%refinement_levels                             ! loop over refinement levels
            b = first_node(level=level)                                 ! first node at level
            do bb=1, nodes_number_at_level(level=level), TREE_RATIO     ! loop over nodes at level
               parent = parent_node(node=b + bb - 1)                    ! parent of the current node
               if (node(parent)%is_allocated()) then                    ! check if parent exist
                  facet_id = node(parent)%facet_id()                    ! store parent facets IDs list
                  if (facet_id%ids_number > 0) then                     ! check if facets list still has facets
                     do bbb=b + bb - 1, b + bb -1 + TREE_RATIO - 1
                        call node(bbb)%add_facets(facet_id=facet_id, &
                                                  facet=facet,       &
                                                  is_exclusive=.true.) ! add facets to node
                     enddo
                  endif
               endif
            enddo
         enddo
      endif

      ! destroy void nodes
      do level=self%refinement_levels, 0, -1                        ! loop over refinement levels
         b = first_node(level=level)                                ! first node at level
         do bb=1, nodes_number_at_level(level=level)                ! loop over nodes at level
            bbb = b + bb - 1                                        ! node numeration in tree
            if (.not.node(bbb)%has_facets()) call node(bbb)%destroy ! destroy void node
         enddo
      enddo

      ! update AABB extents
      ! do level=self%refinement_levels, 0, -1           ! loop over refinement levels
      !    b = first_node(level=level)                   ! first node at level
      !    do bb=1, nodes_number_at_level(level=level)   ! loop over nodes at level
      !       bbb = b + bb - 1                           ! node numeration in tree
      !       call node(bbb)%update_extents(facet=facet) ! update extents
      !    enddo
      ! enddo
   endassociate
   endsubroutine distribute_facets_tree

   pure function has_children(self, node)
   !< Return true if node has at least one child allocated.
   class(aabb_tree_object), intent(in) :: self         !< AABB tree.
   integer(I4P),            intent(in) :: node         !< Node queried.
   logical                             :: has_children !< Check result.
   integer                             :: n, nn        !< Counter.

   has_children = .false.
   n = first_child_node(node=node)
   if (n<=self%nodes_number - TREE_RATIO + 1) then
      do nn=n, n + TREE_RATIO - 1
         if (self%node(nn)%is_allocated()) then
            has_children = .true.
            return
         endif
      enddo
   endif
   endfunction has_children

   pure subroutine initialize(self, refinement_levels, tree_kind, facet, largest_edge_len, bmin, bmax, do_facets_distribute, &
                              is_exclusive, do_update_extents)
   !< Initialize AABB tree.
   class(aabb_tree_object), intent(inout)        :: self                  !< AABB tree.
   integer(I4P),            intent(in), optional :: refinement_levels     !< AABB refinement levels.
   integer(I4P),            intent(in), optional :: tree_kind             !< AABB_TREE_OCTREE (default) or AABB_TREE_SAH_BVH.
   type(facet_object),      intent(in), optional :: facet(:)              !< Facets list.
   real(R8P),               intent(in), optional :: largest_edge_len      !< Largest edge lenght.
   type(vector_R8P),        intent(in), optional :: bmin                  !< Minimum point of AABB.
   type(vector_R8P),        intent(in), optional :: bmax                  !< Maximum point of AABB.
   logical,                 intent(in), optional :: do_facets_distribute  !< Sentinel to enable/disable facets distribution.
   logical,                 intent(in), optional :: is_exclusive          !< Sentinel to enable/disable exclusive addition.
   logical,                 intent(in), optional :: do_update_extents     !< Sentinel to enable/disable AABB extents update.
   integer(I4P)                                  :: refinement_levels_    !< AABB refinement levels, local variable.
   integer(I4P)                                  :: tree_kind_            !< AABB tree kind, local variable.
   logical                                       :: do_facets_distribute_ !< Sentinel to enable/dis. facets distribution, local var.
   integer(I4P)                                  :: level                 !< Counter.
   integer(I4P)                                  :: b, bb, bbb, bbbb      !< Counter.
   integer(I4P)                                  :: parent                !< Parent node index.
   type(aabb_object)                             :: octant(8)             !< AABB octants.

   refinement_levels_ = self%refinement_levels
   tree_kind_         = self%tree_kind
   call self%destroy
   self%refinement_levels = refinement_levels_ ; if (present(refinement_levels)) self%refinement_levels = refinement_levels
   self%tree_kind         = tree_kind_         ; if (present(tree_kind))         self%tree_kind         = tree_kind
   do_facets_distribute_ = .true. ; if (present(do_facets_distribute)) do_facets_distribute_ = do_facets_distribute

   ! Resolve the AABB_AUTO_REFINEMENT sentinel into a concrete depth. Auto-tune
   ! needs the facet count, so it requires `facet` to be present. If not present,
   ! fall back to AUTO_MIN_LEVELS so the tree still builds (an octree without
   ! facets has no useful depth anyway — it is just an empty hierarchy).
   if (self%refinement_levels == AABB_AUTO_REFINEMENT) then
      if (present(facet)) then
         self%refinement_levels = auto_refinement_levels(facets_number=size(facet))
      else
         self%refinement_levels = AUTO_MIN_LEVELS
      endif
   endif

   select case (self%tree_kind)
   case (AABB_TREE_SAH_BVH)
      ! Binary BVH with bucketed surface-area heuristic. Self-tunes its depth;
      ! `refinement_levels` is ignored. Requires `facet` to be present (no
      ! triangles -> nothing to partition); `bmin`/`bmax` and `largest_edge_len`
      ! are not used.
      if (present(facet)) then
         call self%build_bvh_sah(facet=facet)
      else
         self%is_initialized = .true.   ! degenerate: empty tree
      endif

   case default  ! AABB_TREE_OCTREE — the original 8-way space-partitioning path.
      if (self%refinement_levels >= 0) then
         self%nodes_number = nodes_number(refinement_levels=self%refinement_levels)
         allocate(self%node(0:self%nodes_number-1))
         call self%node(0)%initialize(facet=facet, bmin=bmin, bmax=bmax)
         levels_loop: do level=1, self%refinement_levels                             ! loop over refinement levels
            b = first_node(level=level)                                              ! first node at level
            do bb=1, nodes_number_at_level(level=level), TREE_RATIO                  ! loop over nodes at level
               bbb = b + bb - 1                                                      ! node numeration in tree
               parent = parent_node(node=bbb)                                        ! parent of the current node
               if (self%node(parent)%is_allocated()) then                            ! create children nodes
                  call self%node(parent)%compute_octants(octant=octant)              ! compute parent AABB octants
                  if (present(largest_edge_len)) then
                     if (largest_edge_len > octant(1)%median()) then                 ! check if refinement has sense
                        ! a further refinement does not have sense
                        self%refinement_levels = level - 1                           ! set rifinement to the previous one
                        exit levels_loop                                             ! exi loop
                     endif
                  endif
                  do bbbb=0, TREE_RATIO-1                                            ! loop over children
                     call self%node(bbb+bbbb)%initialize(bmin=octant(bbbb+1)%bmin, &
                                                         bmax=octant(bbbb+1)%bmax)   ! initialize node
                  enddo
               endif
            enddo
         enddo levels_loop
         if (present(facet).and.(do_facets_distribute_)) call self%distribute_facets_tree(facet=facet)
         self%is_initialized = .true.
      endif
   end select
   endsubroutine initialize

   pure subroutine build_bvh_sah(self, facet)
   !< Build a binary BVH over the given facet list using bucketed surface-area heuristic.
   !<
   !< Algorithm — top-down, recursive:
   !<   1. Compute the union of triangle bboxes (the "node bbox") and the bbox of
   !<      triangle centroids (the "centroid bbox", which guides binning).
   !<   2. If N <= LEAF_TARGET, write a leaf: store all facet ids in the node's
   !<      facet_id list, mark left_child = right_child = 0, return.
   !<   3. Else, on the longest axis of the centroid bbox, bin centroids into
   !<      `BVH_BUCKETS` buckets. For each of `BVH_BUCKETS - 1` candidate splits,
   !<      compute SAH cost:
   !<          cost(split) = T_TRAV
   !<                      + (SA_L / SA_parent) * N_L
   !<                      + (SA_R / SA_parent) * N_R
   !<      Pick the minimum. If best split cost > N (the leaf cost with T_INT=1),
   !<      write a leaf anyway — recursing would not help.
   !<   4. Partition indices by which side of the split each centroid falls on
   !<      (two-pass: count, then fill).
   !<   5. Allocate two child nodes by incrementing `next_idx` (a shared counter
   !<      threaded through the recursion via a module-private state object), and
   !<      recurse into each.
   !<
   !< Storage: builds into a temporary `nodes_tmp` array sized for the worst case
   !< (4 nodes per leaf-target's worth of facets, well above the 2N/T-1 maximum for
   !< a balanced tree). After the build, `move_alloc`s into `self%node(:)` with the
   !< exact final size.
   class(aabb_tree_object), intent(inout) :: self                 !< AABB tree.
   type(facet_object),      intent(in)    :: facet(:)             !< Facets list.
   type(aabb_node_object), allocatable    :: nodes_tmp(:)         !< Working buffer for the build.
   integer(I4P),           allocatable    :: indices(:)           !< Facet indices owned by each pending build call.
   integer(I4P)                           :: n_facets, n_max, n_used, f

   call self%destroy
   self%tree_kind = AABB_TREE_SAH_BVH

   n_facets = size(facet)
   if (n_facets <= 0) then
      self%is_initialized = .true.
      return
   endif

   ! Worst-case node count: a perfectly-unbalanced tree splitting one triangle at
   ! a time would need 2*N - 1 nodes. SAH never produces that, but pick a safe
   ! upper bound and shrink after the build via move_alloc.
   n_max = max(2 * n_facets, 8_I4P)
   allocate(nodes_tmp(0:n_max - 1))

   ! Seed the recursion with all facet indices [1..n_facets].
   allocate(indices(n_facets))
   do f = 1, n_facets
      indices(f) = f
   enddo

   n_used = 1   ! node 0 is the root, already accounted for
   call build_node(nodes=nodes_tmp, idx=indices, lo=1_I4P, hi=n_facets,                                 &
                   facet=facet, this_node=0_I4P, next_idx=n_used, n_max=n_max)

   ! Compact to exact size and hand off to self.
   allocate(self%node(0:n_used - 1))
   do f = 0, n_used - 1
      self%node(f) = nodes_tmp(f)
   enddo
   self%nodes_number   = n_used
   deallocate(nodes_tmp, indices)

   ! Flatten the leaf facets into the contiguous distance payload (issue #19 §B1).
   call self%build_distance_payload(facet=facet)

   self%is_initialized = .true.
   endsubroutine build_bvh_sah

   pure subroutine build_distance_payload(self, facet)
   !< Flatten every BVH leaf's facets into the tree's contiguous
   !< `facet_distance_payload` array, leaf-grouped (issue #19 §B1).
   !<
   !< Walks the node array in index order; for each leaf (a node with
   !< `payload_count`-eligible facets — i.e. `facet_id%ids_number > 0` and no
   !< children) it appends that leaf's facets to `self%payload(:)` contiguously
   !< and records the leaf's `(payload_first, payload_count)` slice. After this
   !< pass the distance traversal never touches `facet_id` again: a leaf scan is
   !< a stride-1 sweep over `self%payload(first : first+count-1)`.
   !<
   !< The total payload length equals the sum of leaf facet counts. For the SAH
   !< BVH this equals `n_facets` exactly (the partition is a permutation — every
   !< facet lands in exactly one leaf), but the code does not rely on that: it
   !< sizes the array from the actual leaf-count sum so it stays correct if the
   !< builder ever changes.
   class(aabb_tree_object), intent(inout) :: self        !< AABB tree.
   type(facet_object),      intent(in)    :: facet(:)    !< Facets list (read-only source of the triangle metrix).
   integer(I4P)                           :: n, total, cursor, k, fid !< Counters.

   if (allocated(self%payload)) deallocate(self%payload)
   if (self%nodes_number <= 0) return

   ! First pass: total payload length = sum of leaf facet counts.
   total = 0_I4P
   do n = 0, self%nodes_number - 1
      if (self%node(n)%get_left_child() == 0 .and. self%node(n)%get_right_child() == 0) then
         total = total + self%node(n)%facet_id_count()
      endif
   enddo
   if (total <= 0_I4P) return
   allocate(self%payload(total))

   ! Second pass: emit each leaf's facets contiguously, record its slice.
   cursor = 1_I4P
   do n = 0, self%nodes_number - 1
      if (self%node(n)%get_left_child() /= 0 .or. self%node(n)%get_right_child() /= 0) cycle
      associate (count => self%node(n)%facet_id_count())
         if (count <= 0_I4P) cycle
         call self%node(n)%set_payload_slice(first=cursor, count=count)
         do k = 1_I4P, count
            fid = self%node(n)%facet_id_at(k)
            self%payload(cursor)%v1     = [facet(fid)%vertex(1)%x, facet(fid)%vertex(1)%y, facet(fid)%vertex(1)%z]
            self%payload(cursor)%e12    = [facet(fid)%E12%x, facet(fid)%E12%y, facet(fid)%E12%z]
            self%payload(cursor)%e13    = [facet(fid)%E13%x, facet(fid)%E13%y, facet(fid)%E13%z]
            self%payload(cursor)%a      = facet(fid)%a
            self%payload(cursor)%b      = facet(fid)%b
            self%payload(cursor)%c      = facet(fid)%c
            self%payload(cursor)%det    = facet(fid)%det
            self%payload(cursor)%facet_id = fid
            cursor = cursor + 1_I4P
         enddo
      end associate
   enddo
   endsubroutine build_distance_payload

   pure recursive subroutine build_node(nodes, idx, lo, hi, facet, this_node, next_idx, n_max)
   !< One recursive build step. Operates on the slice `idx(lo:hi)` of facet ids:
   !< chooses to write `this_node` as a leaf or to split and recurse.
   !<
   !< `idx` is permuted in place during partitioning — on entry `idx(lo:hi)` lists
   !< the facets to be handled, on exit those slots have been rearranged so that
   !< the left subtree owns `idx(lo:mid)` and the right owns `idx(mid+1:hi)`.
   type(aabb_node_object), intent(inout) :: nodes(0:)        !< Working node buffer.
   integer(I4P),           intent(inout) :: idx(:)           !< Facet-id permutation.
   integer(I4P),           intent(in)    :: lo, hi           !< Inclusive range of `idx` owned by this call.
   type(facet_object),     intent(in)    :: facet(:)         !< Facets list (read-only).
   integer(I4P),           intent(in)    :: this_node        !< Index of the node we are writing into `nodes(:)`.
   integer(I4P),           intent(inout) :: next_idx         !< Next free slot in `nodes(:)` (incremented as children are allocated).
   integer(I4P),           intent(in)    :: n_max            !< Upper bound on nodes_tmp size; safety guard.
   integer(I4P), parameter :: BVH_LEAF_TARGET = 64_I4P       !< Below this, write a leaf without trying to split.
   integer(I4P), parameter :: BVH_BUCKETS     = 16_I4P       !< SAH bucket count per axis.
   real(R8P),    parameter :: BVH_T_TRAV      = 0.125_R8P    !< Traversal cost relative to one triangle test (PBRT-typical).
   type(vector_R8P)   :: node_bmin, node_bmax       !< Union of triangle bboxes for this node.
   type(vector_R8P)   :: centroid_bmin, centroid_bmax !< Bbox of triangle centroids (used for binning).
   real(R8P)          :: extent(3), axis_max, sa_parent, best_cost, leaf_cost, c, split_pos
   real(R8P)          :: bucket_bmin(3, BVH_BUCKETS), bucket_bmax(3, BVH_BUCKETS)
   integer(I4P)       :: bucket_count(BVH_BUCKETS)
   real(R8P)          :: left_bmin(3), left_bmax(3), right_bmin(3), right_bmax(3)
   real(R8P)          :: pre_bmin(3, BVH_BUCKETS), pre_bmax(3, BVH_BUCKETS)
   integer(I4P)       :: pre_count(BVH_BUCKETS)
   integer(I4P)       :: n, b, split_axis, best_axis, best_split, k, left_count
   integer(I4P)       :: i, j, tmp, mid
   real(R8P)          :: bucket_size, axis_bmin, axis_bmax
   integer(I4P), allocatable :: leaf_ids(:)

   n = hi - lo + 1

   ! Step 1: compute node bbox and centroid bbox.
   call compute_node_bboxes(facet, idx, lo, hi, node_bmin, node_bmax, centroid_bmin, centroid_bmax)
   call write_node_bbox(nodes(this_node), node_bmin, node_bmax)

   ! Step 2: small enough -> leaf.
   if (n <= BVH_LEAF_TARGET) then
      call write_leaf(nodes(this_node), idx(lo:hi), n)
      return
   endif

   ! If the centroid bbox is degenerate (all centroids coincide), no split can
   ! improve the cost. Write a leaf even if N is large.
   extent(1) = centroid_bmax%x - centroid_bmin%x
   extent(2) = centroid_bmax%y - centroid_bmin%y
   extent(3) = centroid_bmax%z - centroid_bmin%z
   axis_max = max(extent(1), extent(2), extent(3))
   if (axis_max <= 0._R8P) then
      call write_leaf(nodes(this_node), idx(lo:hi), n)
      return
   endif

   ! Step 3: bucketed SAH. Try each of the three axes; pick the best split overall.
   best_cost  = real(n, R8P) ! leaf cost = N * T_INT (T_INT = 1)
   leaf_cost  = best_cost
   best_axis  = 0
   best_split = 0
   sa_parent  = surface_area(node_bmin, node_bmax)
   if (sa_parent <= 0._R8P) sa_parent = 1._R8P   ! degenerate; avoid division by zero

   do split_axis = 1, 3
      if (extent(split_axis) <= 0._R8P) cycle    ! skip flat axes (binning is undefined)

      axis_bmin = component(centroid_bmin, split_axis)
      axis_bmax = component(centroid_bmax, split_axis)
      bucket_size = (axis_bmax - axis_bmin) / real(BVH_BUCKETS, R8P)

      ! Bin centroids into BVH_BUCKETS buckets.
      bucket_count = 0
      do b = 1, BVH_BUCKETS
         bucket_bmin(:, b) =  huge(0._R8P)
         bucket_bmax(:, b) = -huge(0._R8P)
      enddo
      do k = lo, hi
         c = component(facet(idx(k))%centroid, split_axis) - axis_bmin
         b = min(BVH_BUCKETS, max(1_I4P, int(c / bucket_size, I4P) + 1_I4P))
         bucket_count(b) = bucket_count(b) + 1
         call expand_bbox(bucket_bmin(:, b), bucket_bmax(:, b), facet(idx(k))%bb(1), facet(idx(k))%bb(2))
      enddo

      ! Sweep left -> right prefix: pre_count(b) = sum of counts in buckets 1..b,
      ! pre_bmin/bmax(b) = bbox union of buckets 1..b. Used to evaluate left side
      ! at each of the BVH_BUCKETS - 1 candidate split positions.
      pre_count(1) = bucket_count(1)
      pre_bmin(:, 1) = bucket_bmin(:, 1)
      pre_bmax(:, 1) = bucket_bmax(:, 1)
      do b = 2, BVH_BUCKETS
         pre_count(b)   = pre_count(b - 1) + bucket_count(b)
         pre_bmin(:, b) = min(pre_bmin(:, b - 1), bucket_bmin(:, b))
         pre_bmax(:, b) = max(pre_bmax(:, b - 1), bucket_bmax(:, b))
      enddo

      ! Sweep right -> left, evaluating split between bucket b and bucket b+1.
      right_bmin =  huge(0._R8P)
      right_bmax = -huge(0._R8P)
      do b = BVH_BUCKETS, 2, -1
         right_bmin = min(right_bmin, bucket_bmin(:, b))
         right_bmax = max(right_bmax, bucket_bmax(:, b))
         left_count = pre_count(b - 1)
         left_bmin  = pre_bmin(:, b - 1)
         left_bmax  = pre_bmax(:, b - 1)
         if (left_count == 0 .or. left_count == n) cycle  ! degenerate split
         c = BVH_T_TRAV                                                     &
           + (surface_area_arr(left_bmin,  left_bmax)  / sa_parent) * real(left_count, R8P) &
           + (surface_area_arr(right_bmin, right_bmax) / sa_parent) * real(n - left_count, R8P)
         if (c < best_cost) then
            best_cost  = c
            best_axis  = split_axis
            best_split = b - 1   ! split between bucket (b-1) and bucket b
         endif
      enddo
   enddo

   ! Step 4: if no split beats the leaf cost, give up and write a leaf.
   if (best_axis == 0 .or. best_cost >= leaf_cost) then
      call write_leaf(nodes(this_node), idx(lo:hi), n)
      return
   endif

   ! Step 4b: partition idx(lo:hi) by which side of the split each centroid falls on.
   axis_bmin = component(centroid_bmin, best_axis)
   axis_bmax = component(centroid_bmax, best_axis)
   bucket_size = (axis_bmax - axis_bmin) / real(BVH_BUCKETS, R8P)
   split_pos = axis_bmin + bucket_size * real(best_split, R8P)
   ! Two-pointer Hoare-style partition: move left-belongers below `j`, right-belongers above.
   i = lo
   j = hi
   do
      do while (i <= j)
         if (component(facet(idx(i))%centroid, best_axis) >= split_pos) exit
         i = i + 1
      enddo
      do while (i <= j)
         if (component(facet(idx(j))%centroid, best_axis) <  split_pos) exit
         j = j - 1
      enddo
      if (i >= j) exit
      tmp    = idx(i)
      idx(i) = idx(j)
      idx(j) = tmp
      i = i + 1
      j = j - 1
   enddo
   mid = i - 1
   ! Defensive: if the partition came out degenerate (can happen with ties on the
   ! split plane), fall back to a count-balanced split. Should be vanishingly rare
   ! because we already rejected `left_count == 0 || == n` above, but the bin->plane
   ! conversion can re-introduce a tie.
   if (mid < lo .or. mid >= hi) mid = lo + n / 2 - 1

   ! Step 5: allocate two children and recurse.
   if (next_idx + 1 >= n_max) error stop 'aabb_tree%build_node: node buffer exhausted (raise n_max)'
   call nodes(this_node)%set_children(left_child=next_idx, right_child=next_idx + 1_I4P)
   block
      integer(I4P) :: left_idx, right_idx
      left_idx  = next_idx
      right_idx = next_idx + 1_I4P
      next_idx  = next_idx + 2_I4P
      call build_node(nodes=nodes, idx=idx, lo=lo,    hi=mid, facet=facet, &
                      this_node=left_idx,  next_idx=next_idx, n_max=n_max)
      call build_node(nodes=nodes, idx=idx, lo=mid+1, hi=hi,  facet=facet, &
                      this_node=right_idx, next_idx=next_idx, n_max=n_max)
   end block
   contains
      pure function component(v, axis) result(c)
      !< Return v%x / v%y / v%z by integer axis index (1, 2, 3).
      type(vector_R8P), intent(in) :: v
      integer(I4P),     intent(in) :: axis
      real(R8P)                    :: c
      select case (axis)
      case (1) ; c = v%x
      case (2) ; c = v%y
      case (3) ; c = v%z
      end select
      endfunction component
   endsubroutine build_node

   pure subroutine compute_node_bboxes(facet, idx, lo, hi, node_bmin, node_bmax, centroid_bmin, centroid_bmax)
   !< Compute the union of triangle bboxes (for SAH surface areas) and the
   !< bbox of triangle centroids (for binning) over `idx(lo:hi)`.
   type(facet_object), intent(in)  :: facet(:)
   integer(I4P),       intent(in)  :: idx(:)
   integer(I4P),       intent(in)  :: lo, hi
   type(vector_R8P),   intent(out) :: node_bmin, node_bmax
   type(vector_R8P),   intent(out) :: centroid_bmin, centroid_bmax
   integer(I4P)                    :: k, f

   node_bmin     = vector_R8P( huge(0._R8P),  huge(0._R8P),  huge(0._R8P))
   node_bmax     = vector_R8P(-huge(0._R8P), -huge(0._R8P), -huge(0._R8P))
   centroid_bmin = node_bmin
   centroid_bmax = node_bmax
   do k = lo, hi
      f = idx(k)
      node_bmin%x = min(node_bmin%x, facet(f)%bb(1)%x)
      node_bmin%y = min(node_bmin%y, facet(f)%bb(1)%y)
      node_bmin%z = min(node_bmin%z, facet(f)%bb(1)%z)
      node_bmax%x = max(node_bmax%x, facet(f)%bb(2)%x)
      node_bmax%y = max(node_bmax%y, facet(f)%bb(2)%y)
      node_bmax%z = max(node_bmax%z, facet(f)%bb(2)%z)
      centroid_bmin%x = min(centroid_bmin%x, facet(f)%centroid%x)
      centroid_bmin%y = min(centroid_bmin%y, facet(f)%centroid%y)
      centroid_bmin%z = min(centroid_bmin%z, facet(f)%centroid%z)
      centroid_bmax%x = max(centroid_bmax%x, facet(f)%centroid%x)
      centroid_bmax%y = max(centroid_bmax%y, facet(f)%centroid%y)
      centroid_bmax%z = max(centroid_bmax%z, facet(f)%centroid%z)
   enddo
   endsubroutine compute_node_bboxes

   pure subroutine write_node_bbox(node, bmin, bmax)
   !< Initialize this node's underlying aabb_object with the given bbox.
   !< Reuses aabb_object%initialize so the allocatable component is allocated
   !< consistently with octree-built nodes.
   type(aabb_node_object), intent(inout) :: node
   type(vector_R8P),       intent(in)    :: bmin, bmax

   call node%initialize(bmin=bmin, bmax=bmax)
   endsubroutine write_node_bbox

   pure subroutine write_leaf(node, ids, n)
   !< Write a leaf node: bbox already set by `write_node_bbox`; populate facet_id
   !< list with the given facet ids and zero the child links. Uses the dedicated
   !< `set_facet_ids` method on aabb_node_object — `add_facets` would re-filter
   !< by centroid-inside-bbox, which is unnecessary here (the partition already
   !< did that work) and slow.
   type(aabb_node_object), intent(inout) :: node
   integer(I4P),           intent(in)    :: ids(:)
   integer(I4P),           intent(in)    :: n

   call node%set_facet_ids(ids=ids(1:n))
   call node%set_children(left_child=0_I4P, right_child=0_I4P)
   endsubroutine write_leaf

   pure subroutine expand_bbox(bmin, bmax, addmin, addmax)
   !< Grow [bmin, bmax] (passed as length-3 arrays for tight inlining inside the
   !< SAH bucket loops) to include the bbox [addmin, addmax].
   real(R8P),        intent(inout) :: bmin(3), bmax(3)
   type(vector_R8P), intent(in)    :: addmin, addmax

   bmin(1) = min(bmin(1), addmin%x); bmax(1) = max(bmax(1), addmax%x)
   bmin(2) = min(bmin(2), addmin%y); bmax(2) = max(bmax(2), addmax%y)
   bmin(3) = min(bmin(3), addmin%z); bmax(3) = max(bmax(3), addmax%z)
   endsubroutine expand_bbox

   pure function surface_area(bmin, bmax) result(sa)
   !< Surface area of an axis-aligned bbox: 2 * (dx*dy + dy*dz + dz*dx).
   type(vector_R8P), intent(in) :: bmin, bmax
   real(R8P)                    :: sa
   real(R8P)                    :: dx, dy, dz

   dx = max(bmax%x - bmin%x, 0._R8P)
   dy = max(bmax%y - bmin%y, 0._R8P)
   dz = max(bmax%z - bmin%z, 0._R8P)
   sa = 2._R8P * (dx * dy + dy * dz + dz * dx)
   endfunction surface_area

   pure function surface_area_arr(bmin, bmax) result(sa)
   !< Surface area variant taking length-3 arrays (used inside the SAH bucket loop
   !< where bboxes are accumulated component-wise).
   real(R8P), intent(in) :: bmin(3), bmax(3)
   real(R8P)             :: sa
   real(R8P)             :: dx, dy, dz

   dx = max(bmax(1) - bmin(1), 0._R8P)
   dy = max(bmax(2) - bmin(2), 0._R8P)
   dz = max(bmax(3) - bmin(3), 0._R8P)
   sa = 2._R8P * (dx * dy + dy * dz + dz * dx)
   endfunction surface_area_arr

   function loop_node(self, facet, aabb_facet, b, l) result(again)
   !< Loop over all nodes.
   !<
   !< @note Impure function: return data of each allocated node exploiting saved local counter.
   class(aabb_tree_object),         intent(in)            :: self          !< AABB tree.
   type(facet_object),              intent(in),  optional :: facet(:)      !< Whole facets list.
   integer(I4P),                    intent(out), optional :: b             !< Current AABB ID.
   integer(I4P),                    intent(out), optional :: l             !< Current AABB level.
   type(facet_object), allocatable, intent(out), optional :: aabb_facet(:) !< AABB facets list.
   logical                                                :: again         !< Flag continuing the loop.
   integer(I4P), save                                     :: bb = -1       !< AABB ID counter.
   integer(I4P)                                           :: bbb           !< Counter.

   again = .false.
   if (allocated(self%node)) then
      if (bb==-1) then
         ! get first allocated node
         do bbb=0, self%nodes_number - 1
            if (self%node(bbb)%is_allocated()) then
               again = .true.
               if (present(facet).and.present(aabb_facet)) call self%node(bbb)%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
               exit
            endif
         enddo
         bb = bbb
      elseif (bb<self%nodes_number - 1) then
         do bbb=bb+1, self%nodes_number - 1
            if (self%node(bbb)%is_allocated()) then
               again = .true.
               if (present(facet).and.present(aabb_facet)) call self%node(bbb)%get_aabb_facets(facet=facet, aabb_facet=aabb_facet)
               exit
            endif
         enddo
         bb = bbb
      else
         bb = -1
         again = .false.
      endif
   endif
   if (present(b)) b = bb
   if (present(l)) l = level(b)
   endfunction loop_node

   function ray_intersections_number(self, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number.
   class(aabb_tree_object), intent(in) :: self                 !< AABB tree.
   type(facet_object),      intent(in) :: facet(:)             !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin           !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction        !< Ray direction.
   integer(I4P)                        :: intersections_number !< Intersection number.
   ! integer(I4P)                        :: level                !< Counter.
   ! integer(I4P)                        :: b, bb, bbb           !< Counter.

   intersections_number = self%ray_intersections_number_node(n=0, facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
   ! intersections_number = 0
   ! associate(node=>self%node)
   !    do level=0, self%refinement_levels                  ! loop over refinement levels
   !       b = first_node(level=level)                      ! first node at finest level
   !       do bb=1, nodes_number_at_level(level=level)      ! loop over nodes at level
   !          bbb = b + bb - 1                              ! node numeration in tree
   !          if (node(bbb)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) then
   !             intersections_number = intersections_number + &
   !                                    node(bbb)%ray_intersections_number(facet=facet, &
   !                                                                       ray_origin=ray_origin, ray_direction=ray_direction)
   !          endif
   !       enddo
   !    enddo
   ! endassociate
   endfunction ray_intersections_number

   subroutine save_geometry_tecplot_ascii(self, file_name)
   !< Save AABB tree boxes geometry into Tecplot ascii file.
   class(aabb_tree_object), intent(in) :: self       !< AABB tree.
   character(*),            intent(in) :: file_name  !< File name.
   integer(I4P)                        :: file_unit  !< File unit.
   integer(I4P)                        :: level      !< Counter.
   integer(I4P)                        :: b, bb, bbb !< Counter.

   associate(node=>self%node)
      if (self%is_initialized) then
         open(newunit=file_unit, file=trim(adjustl(file_name)))
         write(file_unit, '(A)') 'VARIABLES=x y z'
         do level=0, self%refinement_levels
            b = first_node(level=level)
            do bb=1, nodes_number_at_level(level=level)
               bbb = b + bb - 1
               call node(bbb)%save_geometry_tecplot_ascii(file_unit=file_unit, aabb_name='aabb-l_'//trim(str(level, .true.))//&
                                                                                             '-b_'//trim(str(bbb, .true.)))
            enddo
         enddo
         close(file_unit)
      endif
   endassociate
   endsubroutine save_geometry_tecplot_ascii

   elemental subroutine translate(self, delta)
   !< Translate AABB tree by delta.
   class(aabb_tree_object), intent(inout) :: self  !< AABB.
   type(vector_R8P),        intent(in)    :: delta !< Delta of translation.
   integer(I4P)                           :: n     !< Counter.

   if (self%nodes_number > 0) then
      do n=0, self%nodes_number - 1
         call self%node(n)%translate(delta=delta)
      enddo
   endif
   endsubroutine translate

   ! operators
   ! =
   pure subroutine aabb_tree_assign_aabb_tree(lhs, rhs)
   !< Operator `=`.
   class(aabb_tree_object), intent(inout) :: lhs !< Left hand side.
   type(aabb_tree_object),  intent(in)    :: rhs !< Right hand side.
   integer                                :: b   !< Counter.

   if (allocated(lhs%node)) then
      do b=1, lhs%nodes_number
         call lhs%node%destroy
      enddo
      deallocate(lhs%node)
   endif
   if (allocated(lhs%payload)) deallocate(lhs%payload)
   lhs%tree_kind         = rhs%tree_kind
   lhs%refinement_levels = rhs%refinement_levels
   lhs%nodes_number      = rhs%nodes_number
   if (allocated(rhs%node)) then
      allocate(lhs%node(0:lhs%nodes_number-1))
      do b=0, lhs%nodes_number-1
         lhs%node(b) = rhs%node(b)
      enddo
   endif
   if (allocated(rhs%payload)) lhs%payload = rhs%payload
   lhs%is_initialized = rhs%is_initialized
   lhs%use_index      = rhs%use_index
   endsubroutine aabb_tree_assign_aabb_tree

   ! private methods
   pure subroutine enumerate_children(self, n, out_idx, nchild)
   !< List the allocated children of node `n` into `out_idx(1:nchild)`, in no
   !< particular order. Dispatches on `self%tree_kind`:
   !<
   !<  - AABB_TREE_OCTREE: children are the up-to-8 consecutive nodes starting at
   !<    `first_child_node(n)`. Empty slots are silently skipped via the
   !<    `is_allocated()` test; the original `has_children` predicate is implicit
   !<    in `nchild == 0` on return.
   !<
   !<  - AABB_TREE_SAH_BVH: children are the up-to-2 explicit indices stored in
   !<    `left_child` / `right_child`. A leaf has both at 0 -> `nchild = 0`.
   !<
   !< Sized for the octree (8 slots) so traversal callers can use a single fixed
   !< buffer regardless of tree kind. The BVH path uses only the first 2.
   class(aabb_tree_object), intent(in)  :: self        !< AABB tree.
   integer(I4P),            intent(in)  :: n           !< Node index.
   integer(I4P),            intent(out) :: out_idx(:)  !< Output: child indices (caller-supplied buffer, >= TREE_RATIO slots).
   integer(I4P),            intent(out) :: nchild      !< Number of allocated children found.
   integer(I4P)                         :: fcn, i, lc, rc

   nchild = 0
   if (self%tree_kind == AABB_TREE_SAH_BVH) then
      lc = self%node(n)%get_left_child()
      rc = self%node(n)%get_right_child()
      if (lc > 0) then ; nchild = nchild + 1 ; out_idx(nchild) = lc ; endif
      if (rc > 0) then ; nchild = nchild + 1 ; out_idx(nchild) = rc ; endif
   else
      ! Octree: implicit indexing.
      fcn = first_child_node(node=n)
      if (fcn > self%nodes_number - TREE_RATIO) return    ! out of range -> leaf
      do i = fcn, fcn + TREE_RATIO - 1
         if (self%node(i)%is_allocated()) then
            nchild = nchild + 1
            out_idx(nchild) = i
         endif
      enddo
   endif
   endsubroutine enumerate_children

   recursive subroutine distance_node(self, n, facet, point, best)
   !< Update `best` squared distance by visiting node `n` and pruning descendants.
   !<
   !< Algorithm (best-first BVH traversal with d^2 pruning):
   !<  1. Test the node's own facets against `best` (every node may carry facets;
   !<     a non-leaf can still hold facets that did not fit deeper levels).
   !<  2. Gather allocated children with their AABB squared distance, sort
   !<     ascending so the most promising subtree is visited first — this makes
   !<     `best` tighten quickly and prune siblings.
   !<  3. Recurse only into children whose box d^2 is < current `best`.
   !<
   !< Correctness depends on the fact that an AABB's squared distance to a point
   !< is a lower bound on the squared distance to any facet inside that AABB.
   !< This is why a child whose box d^2 already exceeds the best can be skipped:
   !< its facets cannot improve the answer.
   class(aabb_tree_object), intent(in)    :: self                       !< AABB tree.
   integer(I4P),            intent(in)    :: n                          !< Current AABB node.
   type(facet_object),      intent(in)    :: facet(:)                   !< Facets list.
   type(vector_R8P),        intent(in)    :: point                      !< Point coordinates.
   real(R8P),               intent(inout) :: best                       !< Running best squared distance.
   real(R8P)                              :: facet_d2                   !< Distance from this node's facets.
   real(R8P)                              :: child_d2(TREE_RATIO)       !< Per-child box d^2 (sized for octree).
   integer(I4P)                           :: child_idx(TREE_RATIO)      !< Per-child node index (sized for octree).
   integer(I4P)                           :: nchild                     !< Number of allocated children.
   integer(I4P)                           :: i, j, swap_i               !< Counter.
   real(R8P)                              :: swap_d                     !< Sort helper.

   associate(node => self%node)
      ! 1. node's own facets — internal nodes can still carry facets (octree).
      !    BVH: stride-1 sweep over the flattened distance payload (issue #19 §B1).
      !    Octree: legacy facet_id-indexed scan into the global facet array.
      if (self%tree_kind == AABB_TREE_SAH_BVH) then
         call self%scan_payload(n=n, point=point, best=best)
      else
         facet_d2 = node(n)%distance_from_facets(facet=facet, point=point)
         if (facet_d2 < best) best = facet_d2
      endif

      ! 2. enumerate allocated children. Tree-kind-agnostic via the helper.
      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0) return
      do i = 1, nchild
         child_d2(i) = node(child_idx(i))%distance(point=point)
      enddo

      ! 3. insertion sort by ascending box d^2 (tiny array, branch-predictable).
      do i = 2, nchild
         swap_d  = child_d2(i)
         swap_i  = child_idx(i)
         j = i - 1
         do while (j >= 1)
            if (child_d2(j) <= swap_d) exit
            child_d2(j + 1)  = child_d2(j)
            child_idx(j + 1) = child_idx(j)
            j = j - 1
         enddo
         child_d2(j + 1)  = swap_d
         child_idx(j + 1) = swap_i
      enddo

      ! 4. recurse with pruning — once a child's box d^2 >= best, all remaining
      !    (sorted) children are at least as far, so we can stop early.
      do i = 1, nchild
         if (child_d2(i) >= best) exit
         call self%distance_node(n=child_idx(i), facet=facet, point=point, best=best)
      enddo
   end associate
   endsubroutine distance_node

   pure subroutine scan_payload(self, n, point, best)
   !< Stride-1 scan of leaf `n`'s slice of the flattened distance payload,
   !< updating `best` squared distance (issue #19 §B1, unsigned path).
   !<
   !< A no-op on internal BVH nodes (`payload_count == 0`). The payload carries
   !< the packed triangle metrix contiguously, so this is a sequential sweep with
   !< no `node -> aabb -> facet_id -> facet(fid)` indirection.
   class(aabb_tree_object), intent(in)    :: self    !< AABB tree.
   integer(I4P),            intent(in)    :: n       !< Node index.
   type(vector_R8P),        intent(in)    :: point   !< Point coordinates.
   real(R8P),               intent(inout) :: best    !< Running best squared distance.
   integer(I4P)                           :: first, count, k, idx !< Slice bounds and counters.
   real(R8P)                              :: d2      !< Candidate squared distance.
   type(vector_R8P)                       :: closest !< Discarded (unsigned path).
   integer(I4P)                           :: region  !< Discarded (unsigned path).

   count = self%node(n)%get_payload_count()
   if (count <= 0_I4P) return
   first = self%node(n)%get_payload_first()
   do k = 0_I4P, count - 1_I4P
      idx = first + k
      associate (p => self%payload(idx))
         call triangle_point_distance(v1=vector_R8P(p%v1(1), p%v1(2), p%v1(3)),       &
                                      e12=vector_R8P(p%e12(1), p%e12(2), p%e12(3)),   &
                                      e13=vector_R8P(p%e13(1), p%e13(2), p%e13(3)),   &
                                      a=p%a, b=p%b, c=p%c, det=p%det, point=point,    &
                                      distance=d2, closest=closest, region=region)
      end associate
      if (d2 < best) best = d2
   enddo
   endsubroutine scan_payload

   recursive subroutine distance_node_with_region(self, n, facet, point, best, best_facet, best_region)
   !< Same best-first descent as `distance_node`, but tracks the facet id and Voronoi region
   !< of the closest point (needed for pseudo-normal sign determination).
   class(aabb_tree_object), intent(in)    :: self                       !< AABB tree.
   integer(I4P),            intent(in)    :: n                          !< Current AABB node.
   type(facet_object),      intent(in)    :: facet(:)                   !< Facets list.
   type(vector_R8P),        intent(in)    :: point                      !< Point coordinates.
   real(R8P),               intent(inout) :: best                       !< Running best squared distance.
   integer(I4P),            intent(inout) :: best_facet                 !< Running best facet id.
   integer(I4P),            intent(inout) :: best_region                !< Running best Voronoi region tag.
   real(R8P)                              :: child_d2(TREE_RATIO)       !< Per-child box d^2 (sized for octree).
   integer(I4P)                           :: child_idx(TREE_RATIO)      !< Per-child node index (sized for octree).
   integer(I4P)                           :: nchild                     !< Number of allocated children.
   integer(I4P)                           :: i, j, swap_i               !< Counter.
   real(R8P)                              :: swap_d                     !< Sort helper.

   associate(node => self%node)
      ! BVH: stride-1 sweep over the flattened distance payload, recording the
      ! winning facet id (issue #19 §B1). The Voronoi region is NOT computed here
      ! — only the winner's region matters and the surface layer recomputes it
      ! once on the final closest facet. Octree: legacy facet_id-indexed scan.
      if (self%tree_kind == AABB_TREE_SAH_BVH) then
         call self%scan_payload_with_facet(n=n, point=point, best=best, best_facet=best_facet)
      else
         call node(n)%update_best_from_facets(facet=facet, point=point, &
                                              best=best, best_facet=best_facet, best_region=best_region)
      endif

      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0) return
      do i = 1, nchild
         child_d2(i) = node(child_idx(i))%distance(point=point)
      enddo

      do i = 2, nchild
         swap_d = child_d2(i)
         swap_i = child_idx(i)
         j = i - 1
         do while (j >= 1)
            if (child_d2(j) <= swap_d) exit
            child_d2(j + 1)  = child_d2(j)
            child_idx(j + 1) = child_idx(j)
            j = j - 1
         enddo
         child_d2(j + 1)  = swap_d
         child_idx(j + 1) = swap_i
      enddo

      do i = 1, nchild
         if (child_d2(i) >= best) exit
         call self%distance_node_with_region(n=child_idx(i), facet=facet, point=point, &
                                             best=best, best_facet=best_facet, best_region=best_region)
      enddo
   end associate
   endsubroutine distance_node_with_region

   pure subroutine scan_payload_with_facet(self, n, point, best, best_facet)
   !< Stride-1 scan of leaf `n`'s payload slice, updating `best` squared distance
   !< and the winning facet id (issue #19 §B1, signed/pseudo-normal path).
   !<
   !< Unlike the octree's `update_best_from_facets`, the Voronoi region is not
   !< tracked per facet here: only the *winning* facet's region is needed, and the
   !< surface layer (`compute_distance`) recomputes it once via
   !< `facet(best_facet)%compute_distance_with_region` on the final closest facet.
   !< Dropping the per-facet region work is part of issue #19 §B5; doing it here
   !< for free since the payload kernel is being written fresh.
   class(aabb_tree_object), intent(in)    :: self       !< AABB tree.
   integer(I4P),            intent(in)    :: n          !< Node index.
   type(vector_R8P),        intent(in)    :: point      !< Point coordinates.
   real(R8P),               intent(inout) :: best       !< Running best squared distance.
   integer(I4P),            intent(inout) :: best_facet !< Running best facet id.
   integer(I4P)                           :: first, count, k, idx !< Slice bounds and counters.
   real(R8P)                              :: d2         !< Candidate squared distance.
   type(vector_R8P)                       :: closest    !< Discarded (region recomputed on the winner).
   integer(I4P)                           :: region     !< Discarded (region recomputed on the winner).

   count = self%node(n)%get_payload_count()
   if (count <= 0_I4P) return
   first = self%node(n)%get_payload_first()
   do k = 0_I4P, count - 1_I4P
      idx = first + k
      associate (p => self%payload(idx))
         call triangle_point_distance(v1=vector_R8P(p%v1(1), p%v1(2), p%v1(3)),       &
                                      e12=vector_R8P(p%e12(1), p%e12(2), p%e12(3)),   &
                                      e13=vector_R8P(p%e13(1), p%e13(2), p%e13(3)),   &
                                      a=p%a, b=p%b, c=p%c, det=p%det, point=point,    &
                                      distance=d2, closest=closest, region=region)
         if (d2 < best) then
            best       = d2
            best_facet = p%facet_id
         endif
      end associate
   enddo
   endsubroutine scan_payload_with_facet

   recursive function ray_intersections_number_node(self, n, facet, ray_origin, ray_direction) result(intersections_number)
   !< Return ray intersections number into a node of the AABB tree.
   !<
   !< Tree-kind-agnostic: child enumeration goes through `enumerate_children`, so
   !< both the octree (up to 8 children) and the BVH (up to 2 children) traverse
   !< through the same code path.
   class(aabb_tree_object), intent(in) :: self                       !< AABB tree.
   integer(I4P),            intent(in) :: n                          !< Current AABB node.
   type(facet_object),      intent(in) :: facet(:)                   !< Facets list.
   type(vector_R8P),        intent(in) :: ray_origin                 !< Ray origin.
   type(vector_R8P),        intent(in) :: ray_direction              !< Ray direction.
   integer(I4P)                        :: intersections_number       !< Intersection number.
   integer(I4P)                        :: child_idx(TREE_RATIO)      !< Child indices buffer (sized for octree).
   integer(I4P)                        :: nchild                     !< Number of allocated children.
   integer(I4P)                        :: i                          !< Counter.

   intersections_number = 0
   associate(node=>self%node)
      if (.not. node(n)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) return
      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0) then
         ! Leaf: return intersection count of this node's own facets.
         intersections_number = node(n)%ray_intersections_number(facet=facet, ray_origin=ray_origin, ray_direction=ray_direction)
      else
         do i = 1, nchild
            intersections_number = intersections_number +                                    &
                                   self%ray_intersections_number_node(n=child_idx(i),        &
                                                                      facet=facet,           &
                                                                      ray_origin=ray_origin, &
                                                                      ray_direction=ray_direction)
         enddo
      endif
   endassociate
   endfunction ray_intersections_number_node

   subroutine intersect_ray_all_tree(self, facet, ray_origin, ray_direction, hits)
   !< Tree-accelerated all-hits ray query (issue #18 §2.5).
   !<
   !< Returns every facet hit by the half-ray `t >= 0`, sorted ascending by `t`.
   !< Equivalent contract to `ray_intersect_all_flat` (sorted, no duplicates), only
   !< faster: prunes whole subtrees whose AABB the ray misses.
   class(aabb_tree_object),      intent(in)  :: self           !< AABB tree.
   type(facet_object),           intent(in)  :: facet(:)       !< Facets list.
   type(vector_R8P),             intent(in)  :: ray_origin     !< Ray origin.
   type(vector_R8P),             intent(in)  :: ray_direction  !< Ray direction.
   type(ray_hit_t), allocatable, intent(out) :: hits(:)        !< Sorted hit records.
   type(ray_hit_t), allocatable              :: tmp(:)         !< Pre-sort buffer.
   integer(I4P)                              :: n_hits         !< Running hit count.
   integer(I4P)                              :: n_overflow     !< Hits dropped because tmp was full (re-runs).
   integer(I4P)                              :: cap            !< Capacity of tmp.

   cap = max(16_I4P, size(facet, kind=I4P) / 8_I4P)  ! guess: ~12.5% of facets are hit
   do
      allocate(tmp(cap))
      n_hits = 0_I4P
      n_overflow = 0_I4P
      call self%intersect_ray_all_node(n=0_I4P, facet=facet, ray_origin=ray_origin, &
                                       ray_direction=ray_direction, tmp=tmp, n_hits=n_hits, &
                                       n_overflow=n_overflow)
      if (n_overflow == 0_I4P) exit
      ! Re-run with a buffer big enough to hold the hits we observed plus the overflow.
      cap = (n_hits + n_overflow) * 2_I4P
      deallocate(tmp)
   enddo
   allocate(hits(n_hits))
   if (n_hits > 0_I4P) hits(1:n_hits) = tmp(1:n_hits)
   call sort_hits_by_t(hits)
   endsubroutine intersect_ray_all_tree

   recursive subroutine intersect_ray_all_node(self, n, facet, ray_origin, ray_direction, &
                                                tmp, n_hits, n_overflow)
   !< Recursive descent under `intersect_ray_all_tree`. At each node: prune via the
   !< slab test; at a leaf, delegate the per-facet sweep to `node%intersect_ray_all_facets`.
   !< Tree-kind-agnostic via `enumerate_children`.
   class(aabb_tree_object),      intent(in)    :: self           !< AABB tree.
   integer(I4P),                 intent(in)    :: n              !< Current AABB node.
   type(facet_object),           intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),             intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),             intent(in)    :: ray_direction  !< Ray direction.
   type(ray_hit_t), allocatable, intent(inout) :: tmp(:)         !< Hit accumulator.
   integer(I4P),                 intent(inout) :: n_hits         !< Running hit count.
   integer(I4P),                 intent(inout) :: n_overflow     !< Capacity-overflow counter.
   integer(I4P)                                :: child_idx(TREE_RATIO)
   integer(I4P)                                :: nchild         !< Number of allocated children.
   integer(I4P)                                :: i              !< Loop counter.
   real(R8P)                                   :: t_near, t_far  !< Slab interval.
   logical                                     :: hits_box       !< Slab-test result.

   associate(node=>self%node)
      if (.not. node(n)%is_allocated()) return
      call node(n)%ray_slab_interval(ray_origin=ray_origin, ray_direction=ray_direction, &
                                     t_near=t_near, t_far=t_far, do_intersect=hits_box)
      if (.not. hits_box) return
      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0_I4P) then
         call node(n)%intersect_ray_all_facets(facet=facet, ray_origin=ray_origin, &
                                               ray_direction=ray_direction, tmp=tmp, &
                                               n_hits=n_hits, n_overflow=n_overflow)
      else
         do i = 1_I4P, nchild
            call self%intersect_ray_all_node(n=child_idx(i), facet=facet, ray_origin=ray_origin, &
                                             ray_direction=ray_direction, tmp=tmp, n_hits=n_hits, &
                                             n_overflow=n_overflow)
         enddo
      endif
   endassociate
   endsubroutine intersect_ray_all_node

   subroutine intersect_ray_first_tree(self, facet, ray_origin, ray_direction, hit, has_hit)
   !< Tree-accelerated closest-hit ray query (issue #18 §2.5).
   !<
   !< Visits children sorted by their slab `t_near` (best-first); prunes any child
   !< whose `t_near >= t_best`. Caller checks `has_hit` before consuming `hit`.
   class(aabb_tree_object), intent(in)  :: self           !< AABB tree.
   type(facet_object),      intent(in)  :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)  :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)  :: ray_direction  !< Ray direction.
   type(ray_hit_t),         intent(out) :: hit            !< Closest hit (valid only if has_hit).
   logical,                 intent(out) :: has_hit        !< True if any facet was hit.
   real(R8P)                            :: t_best         !< Running best t (start at +inf).

   has_hit = .false.
   t_best = MaxR8P
   call self%intersect_ray_first_node(n=0_I4P, facet=facet, ray_origin=ray_origin, &
                                      ray_direction=ray_direction, hit=hit, has_hit=has_hit, t_best=t_best)
   endsubroutine intersect_ray_first_tree

   recursive subroutine intersect_ray_first_node(self, n, facet, ray_origin, ray_direction, &
                                                  hit, has_hit, t_best)
   !< Best-first recursive descent. Children are sorted by their slab `t_near` so
   !< the closest box is visited first; once `t_best` shrinks below a sibling's
   !< `t_near`, that sibling and its subtree are pruned.
   class(aabb_tree_object), intent(in)    :: self           !< AABB tree.
   integer(I4P),            intent(in)    :: n              !< Current AABB node.
   type(facet_object),      intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)    :: ray_direction  !< Ray direction.
   type(ray_hit_t),         intent(inout) :: hit            !< Running closest hit.
   logical,                 intent(inout) :: has_hit        !< True once any hit recorded.
   real(R8P),               intent(inout) :: t_best         !< Running best t (used for pruning).
   integer(I4P)                           :: child_idx(TREE_RATIO)
   real(R8P)                              :: child_tnear(TREE_RATIO)
   integer(I4P)                           :: nchild
   integer(I4P)                           :: i, j
   integer(I4P)                           :: swap_i
   real(R8P)                              :: swap_d
   real(R8P)                              :: t_near, t_far
   logical                                :: hits_box

   associate(node=>self%node)
      if (.not. node(n)%is_allocated()) return
      call node(n)%ray_slab_interval(ray_origin=ray_origin, ray_direction=ray_direction, &
                                     t_near=t_near, t_far=t_far, do_intersect=hits_box)
      if (.not. hits_box) return
      if (t_near >= t_best) return  ! whole subtree is past the running best
      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0_I4P) then
         call node(n)%intersect_ray_first_facets(facet=facet, ray_origin=ray_origin, &
                                                 ray_direction=ray_direction, hit=hit, &
                                                 has_hit=has_hit, t_best=t_best)
      else
         ! Score children by their own slab t_near, then visit ascending.
         do i = 1_I4P, nchild
            if (.not. node(child_idx(i))%is_allocated()) then
               child_tnear(i) = MaxR8P
               cycle
            endif
            call node(child_idx(i))%ray_slab_interval(ray_origin=ray_origin, ray_direction=ray_direction, &
                                                      t_near=child_tnear(i), t_far=t_far, do_intersect=hits_box)
            if (.not. hits_box) child_tnear(i) = MaxR8P
         enddo
         ! Insertion sort (nchild <= 8).
         do i = 2_I4P, nchild
            swap_d = child_tnear(i)
            swap_i = child_idx(i)
            j = i - 1_I4P
            do while (j >= 1_I4P)
               if (child_tnear(j) <= swap_d) exit
               child_tnear(j + 1_I4P) = child_tnear(j)
               child_idx(j + 1_I4P)   = child_idx(j)
               j = j - 1_I4P
            enddo
            child_tnear(j + 1_I4P) = swap_d
            child_idx(j + 1_I4P)   = swap_i
         enddo
         do i = 1_I4P, nchild
            if (child_tnear(i) >= t_best) exit  ! all remaining are pruned
            call self%intersect_ray_first_node(n=child_idx(i), facet=facet, ray_origin=ray_origin, &
                                               ray_direction=ray_direction, hit=hit, has_hit=has_hit, &
                                               t_best=t_best)
         enddo
      endif
   endassociate
   endsubroutine intersect_ray_first_node

   subroutine intersect_ray_any_tree(self, facet, ray_origin, ray_direction, max_t, found)
   !< Tree-accelerated any-hit ray query with early exit (issue #18 §2.5).
   !<
   !< Stops as soon as ANY facet on the surface is hit at `0 <= t <= max_t`.
   !< Use cases: shadow rays (occlusion), point-in-mesh quick-rejection. Caller
   !< passes `max_t = MaxR8P` to mean "any hit anywhere along the forward ray".
   !<
   !< Depth-first traversal — no best-first ordering needed because we stop on
   !< the first hit found, regardless of `t`. Subtrees whose AABB the ray
   !< misses or whose AABB-entry is past `max_t` are pruned.
   class(aabb_tree_object), intent(in)  :: self           !< AABB tree.
   type(facet_object),      intent(in)  :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)  :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)  :: ray_direction  !< Ray direction.
   real(R8P),               intent(in)  :: max_t          !< Upper t bound (parametric).
   logical,                 intent(out) :: found          !< True iff at least one hit in [0, max_t].

   found = .false.
   call self%intersect_ray_any_node(n=0_I4P, facet=facet, ray_origin=ray_origin, &
                                    ray_direction=ray_direction, max_t=max_t, found=found)
   endsubroutine intersect_ray_any_tree

   recursive subroutine intersect_ray_any_node(self, n, facet, ray_origin, ray_direction, max_t, found)
   !< Depth-first early-exit recursion. Returns immediately once `found` is set.
   class(aabb_tree_object), intent(in)    :: self           !< AABB tree.
   integer(I4P),            intent(in)    :: n              !< Current AABB node.
   type(facet_object),      intent(in)    :: facet(:)       !< Facets list.
   type(vector_R8P),        intent(in)    :: ray_origin     !< Ray origin.
   type(vector_R8P),        intent(in)    :: ray_direction  !< Ray direction.
   real(R8P),               intent(in)    :: max_t          !< Upper t bound.
   logical,                 intent(inout) :: found          !< Early-exit flag.
   integer(I4P)                           :: child_idx(TREE_RATIO)
   integer(I4P)                           :: nchild
   integer(I4P)                           :: i
   real(R8P)                              :: t_near, t_far
   logical                                :: hits_box

   if (found) return
   associate(node=>self%node)
      if (.not. node(n)%is_allocated()) return
      call node(n)%ray_slab_interval(ray_origin=ray_origin, ray_direction=ray_direction, &
                                     t_near=t_near, t_far=t_far, do_intersect=hits_box)
      if (.not. hits_box) return
      if (t_near > max_t) return  ! whole subtree begins past max_t
      call self%enumerate_children(n=n, out_idx=child_idx, nchild=nchild)
      if (nchild == 0_I4P) then
         call node(n)%intersect_ray_any_facets(facet=facet, ray_origin=ray_origin, &
                                               ray_direction=ray_direction, max_t=max_t, found=found)
      else
         do i = 1_I4P, nchild
            if (found) return
            call self%intersect_ray_any_node(n=child_idx(i), facet=facet, ray_origin=ray_origin, &
                                             ray_direction=ray_direction, max_t=max_t, found=found)
         enddo
      endif
   endassociate
   endsubroutine intersect_ray_any_node

   pure subroutine sort_hits_by_t(hits)
   !< Insertion sort hits ascending by `t`. Hit counts are typically small (~2..20),
   !< so the O(n^2) constant beats library-call overhead.
   type(ray_hit_t), intent(inout) :: hits(:)
   type(ray_hit_t)                :: key
   integer(I4P)                   :: i, j

   do i = 2_I4P, size(hits, kind=I4P)
      key = hits(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (hits(j)%t <= key%t) exit
         hits(j + 1_I4P) = hits(j)
         j = j - 1_I4P
      enddo
      hits(j + 1_I4P) = key
   enddo
   endsubroutine sort_hits_by_t

   ! non TBP
   pure function first_child_node(node)
   !< Return first child tree node.
   integer(I4P), intent(in) :: node             !< Node queried.
   integer(I4P)             :: first_child_node !< First child tree node.

   first_child_node = node * TREE_RATIO + 1
   endfunction first_child_node

   pure function first_node(level)
   !< Return first tree node at a given level.
   integer(I4P), intent(in) :: level      !< Refinement level queried.
   integer(I4P)             :: first_node !< First tree node at given level.

   first_node = nodes_number(refinement_levels=level-1)
   endfunction first_node

   pure function last_node(level)
   !< Return last tree node at a given level.
   integer(I4P), intent(in) :: level     !< Refinement level queried.
   integer(I4P)             :: last_node !< Last tree node at given level.

   last_node = first_node(level) + nodes_number_at_level(level) - 1
   endfunction last_node

   pure function level(node)
   !< Return level given a node id.
   integer(I4P), intent(in) :: node  !< Node queried.
   integer(I4P)             :: level !< Level of given node.
   integer(I4P)             :: n     !< Counter.

   level = 0
   n = node
   do while (n /= 0)
      n = (n - 1) / TREE_RATIO
      level = level + 1
   enddo
   endfunction level

   pure function local_id(node)
   !< Return local ID of node.
   integer(I4P), intent(in) :: node     !< Node queried.
   integer(I4P)             :: local_id !< Local ID.

   local_id = node - ((node - 1) / TREE_RATIO) * TREE_RATIO
   endfunction local_id

   pure function location_code(node)
   !< Return location code of node.
   integer(I4P), intent(in)  :: node             !< Node queried.
   integer(I4P), allocatable :: location_code(:) !< Location code.
   integer(I4P)              :: node_level       !< Node level.
   integer(I4P)              :: parent           !< Parent node.
   integer(I4P)              :: l                !< Counter.

   node_level = level(node=node)
   if (node_level > 0) then
      allocate(location_code(1:node_level))
      parent = node
      do l=node_level, 1, -1
         location_code(l) = local_id(parent)
         parent = parent_node(node=node)
      enddo
   else
      allocate(location_code(1))
      location_code(1) = 0
   endif
   endfunction location_code

   pure function global_id(location_code) result(node)
   !< Return node ID given a location code.
   integer(I4P), intent(in)  :: location_code(:) !< Location code.
   integer(I4P)              :: node             !< Node global ID.
   integer(I4P)              :: l                !< Counter.

   node = 0
   if (location_code(1) == 0) then
      return
   else
      node = location_code(1)
      do l=2, size(location_code, dim=1)
         node = first_child_node(node)
         node = node + location_code(l) - 1
      enddo
   endif
   endfunction global_id

   pure subroutine next_location_code(location_code, direction, next_code, next_direction)
   !< Return the node next along a given direction.
   integer(I4P),              intent(in)  :: location_code(:) !< Location code queried.
   integer(I4P),              intent(in)  :: direction        !< Direction, 0=Halt, 1=X+, 2=X-, 3=Y+, 4=Y-, 5=Z+, 6=Z-.
   integer(I4P), allocatable, intent(out) :: next_code(:)     !< Next location code of node along given direction.
   integer(I4P),              intent(out) :: next_direction   !< Next direction.
   integer(I4P)                           :: direction_       !< Direction, local variable.
   integer(I4P)                           :: o                !< Counter.

   direction_ = direction
   next_code = location_code
   do o=size(location_code, dim=1), 1, -1
      if (direction_ > 0.and.location_code(o) > 0) then
         next_code(o) = octree_fsm(direction_, location_code(o))%octant
         direction_ = octree_fsm(direction_, location_code(o))%direction
      else
         exit
      endif
   enddo
   next_direction = direction_
   endsubroutine next_location_code

   pure function nodes_number(refinement_levels)
   !< Return total number of tree nodes given the total number refinement levels used.
   integer(I4P), intent(in) :: refinement_levels !< Total number of refinement levels used.
   integer(I4P)             :: nodes_number      !< Total number of tree nodes.
   integer                  :: level             !< Counter.

   nodes_number = 0
   do level=0, refinement_levels
      nodes_number = nodes_number + nodes_number_at_level(level=level)
   enddo
   endfunction nodes_number

   pure function nodes_number_at_level(level) result(nodes_number)
   !< Return number of tree nodes at a given level.
   integer(I4P), intent(in) :: level        !< Refinement level queried.
   integer(I4P)             :: nodes_number !< Number of tree nodes at given level.

   nodes_number = TREE_RATIO ** (level)
   endfunction nodes_number_at_level

   pure function parent_at_level(node, parent_level) result(parent)
   !< Return parent tree node at a given level.
   integer(I4P), intent(in) :: node         !< Node.
   integer(I4P), intent(in) :: parent_level !< Level.
   integer(I4P)             :: parent       !< Parent.
   integer(I4P)             :: n            !< Counter.

   parent = node
   do n=1, level(node) - parent_level
     parent = (parent - 1) / TREE_RATIO
   enddo
   endfunction parent_at_level

   pure function parent_node(node)
   !< Return parent tree node.
   integer(I4P), intent(in) :: node        !< Node queried.
   integer(I4P)             :: parent_node !< Parent tree node.

   parent_node = (node - 1) / TREE_RATIO
   endfunction parent_node

   pure function siblings(node) result(sbs)
   !< Return siblings of a given node.
   integer(I4P), intent(in) :: node     !< Node queried.
   integer(I4P)             :: sbs(1:7) !< Nodes sibling.
   integer(I4P)             :: myid     !< Local node ID into siblings list.
   integer(I4P)             :: i, s     !< Counter.

   myid = local_id(node=node)
   s = 0
   do i=1, TREE_RATIO
      if (i /= myid) then
         s = s + 1
         sbs(s) = ((node - 1) / TREE_RATIO) * TREE_RATIO + i
      endif
   enddo
   endfunction siblings

   pure function str_location_code(code)
   !< Return string of location code of node.
   integer(I4P), intent(in)      :: code(:)           !< Node location code.
   character(len=:), allocatable :: str_location_code !< String of Location code.
   integer(I4P)                  :: c                 !< Counter.

   str_location_code = ''
   do c=1, size(code, dim=1)
      str_location_code = str_location_code//trim(str(no_sign=.true., n=code(c)))
   enddo
   endfunction str_location_code
endmodule fossil_aabb_tree_object
