!< FOSSIL, STL surface class definition.

module fossil_surface_stl_object
!< FOSSIL, STL surface class definition.

use fossil_aabb_tree_object, only : aabb_tree_object, AABB_TREE_OCTREE
use fossil_facet_object, only : facet_object
use fossil_list_id_object, only : list_id_object
use fossil_utils, only : EPS, FRLEN, PI, is_inside_bb, triangle_overlaps_aabb
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_winding_number, only : compute_winding_number => winding_number
use fossil_self_intersection, only : intersection_pair_t, &
                                     compute_self_intersections => find_self_intersections
use fossil_marching_cubes, only : compute_isosurface => extract_isosurface, MC_STATUS_OK
use fossil_decimate, only : compute_decimate => decimate, &
                            DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS
use fossil_boolean, only : compute_boolean => boolean_compute, &
                           BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, &
                           BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, &
                           BOOL_STATUS_EMPTY_INPUT
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use penf, only : I4P, I8P, R8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: surface_stl_object
public :: intersection_pair_t
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL
public :: sign_algorithm_from_string
public :: STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
public :: STATUS_INVALID_INPUT
public :: BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF
public :: BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT

! Point-in-polyhedron algorithm selector for `is_point_inside`, `compute_distance`,
! and `distance` (when `is_signed=.true.`):
!   SIGN_RAY_INTERSECTIONS — count intersections of an axis-aligned ray with the
!                            polyhedron; odd = inside, even = outside. O(N) per query
!                            even with AABB acceleration; sensitive to EPS jitter
!                            near edges/vertices.
!   SIGN_SOLID_ANGLE       — sum the projected solid angles of all facets;
!                            ~±4π = inside, ~0 = outside. O(N) per query, no
!                            AABB acceleration. Robust but slow.
!   SIGN_PSEUDO_NORMAL     — Baerentzen-Aanaes (2005) sign test: the sign of
!                            `dot(point - closest, pseudo_normal_at_closest)`
!                            decides inside/outside without a separate inside
!                            query. Fused with the distance traversal — ~3x
!                            faster than SIGN_RAY_INTERSECTIONS on the same
!                            path. This is the default. Requires
!                            consistently-oriented outward facet normals, which
!                            `sanitize_normals` produces; for open shells where
!                            outward is ill-defined, pass SIGN_RAY_INTERSECTIONS
!                            or SIGN_SOLID_ANGLE explicitly.
integer(I4P), parameter :: SIGN_RAY_INTERSECTIONS = 1_I4P
integer(I4P), parameter :: SIGN_SOLID_ANGLE       = 2_I4P
integer(I4P), parameter :: SIGN_PSEUDO_NORMAL     = 3_I4P

! Status codes returned via the optional `status` argument on mutating procedures.
integer(I4P), parameter :: STATUS_OK               = 0_I4P !< Success.
integer(I4P), parameter :: STATUS_ALLOC_FAIL       = 1_I4P !< Allocation failure.
integer(I4P), parameter :: STATUS_AMBIGUOUS_ARGS   = 2_I4P !< Conflicting optional arguments.
integer(I4P), parameter :: STATUS_FILE_NOT_FOUND   = 3_I4P !< File does not exist.
integer(I4P), parameter :: STATUS_FILE_OPEN_FAIL   = 4_I4P !< File could not be opened for writing.
integer(I4P), parameter :: STATUS_INVALID_INPUT    = 5_I4P !< Loaded data contains NaN/Inf vertex coordinates.

type :: surface_stl_object
   !< FOSSIL STL surface class.
   !<
   !< All components are `private`. External code interacts via accessors:
   !<  - scalar/value getters: `get_facets_number`, `get_bmin`, `get_bmax`, `get_volume`,
   !<    `get_centroid`, `get_header`, `set_header`
   !<  - facet access: `facet_at(i)` returns a `pointer` to a single facet (or `null()`
   !<    if `i` is out of range), `facets_ref()` returns a `pointer` to the whole array.
   !<    Both are intended as read-only views — the library does not enforce this through
   !<    the pointer (Fortran has no `const`-equivalent on components), so by convention
   !<    callers do not write through these accessors. Mutation goes via the surface's
   !<    own TBPs (`translate`, `rotate`, `mirror`, etc.).
   !<  - ownership transfer: `adopt_facets(arr)` moves an allocatable array into the
   !<    surface via `move_alloc`, then auto-runs `analyze`. Used internally by
   !<    `load_from_file`.
   !<  - `aabb` is technically public so that callers can invoke its TBPs
   !<    (e.g. `surface%aabb%set_use_index(...)`); its own components are private.
   integer(I4P),                    private :: facets_number=0           !< Facets number (== size(facet)).
   type(facet_object), allocatable, private :: facet(:)                  !< Facets.
   type(list_id_object),            private :: facet_1_de                !< Facets with one disconnected edges.
   type(list_id_object),            private :: facet_2_de                !< Facets with two disconnected edges.
   type(list_id_object),            private :: facet_3_de                !< Facets with three disconnected edges.
   integer(I4P),                    private :: non_manifold_edges_number=0 !< Number of edges with 3+ incident facets.
   integer(I4P),                    private :: degenerate_facets_removed=0 !< Number of facets removed by the last degenerate-removal pass.
   integer(I4P),                    private :: duplicate_facets_removed=0  !< Number of facets removed by the last duplicate-removal pass.
   type(aabb_tree_object)                   :: aabb            !< AABB tree handle (its own state is private).
   type(vector_R8P),                private :: bmin            !< Bounding-box min.
   type(vector_R8P),                private :: bmax            !< Bounding-box max.
   real(R8P),                       private :: volume=0._R8P   !< Volume bounded by STL surface.
   real(R8P),                       private :: area=0._R8P     !< Total surface area (sum of facet areas; issue #7).
   type(vector_R8P),                private :: centroid        !< Centroid of STL surface.
   character(FRLEN),                private :: header=''       !< STL file header (preserved across load/save).
   type(vertex_pool_object),        private :: vertex_pool     !< Unique-vertex pool (issue #5 stage 1: derived artifact).
   contains
      ! read-only accessors (pure, inlined at -O2, zero data copy for scalars)
      procedure, pass(self) :: get_facets_number !< Return facets_number.
      procedure, pass(self) :: get_bmin          !< Return bmin (bounding-box minimum).
      procedure, pass(self) :: get_bmax          !< Return bmax (bounding-box maximum).
      procedure, pass(self) :: get_volume        !< Return volume.
      procedure, pass(self) :: get_area          !< Return total surface area (issue #7).
      procedure, pass(self) :: get_centroid      !< Return centroid.
      procedure, pass(self) :: get_header        !< Return STL header.
      procedure, pass(self) :: get_non_manifold_edges_number  !< Return count of edges with 3+ incident facets.
      procedure, pass(self) :: get_degenerate_facets_removed  !< Return count of facets dropped by the last degenerate-facet pass.
      procedure, pass(self) :: get_duplicate_facets_removed   !< Return count of facets dropped by the last duplicate-facet pass.
      procedure, pass(self) :: get_vertex_pool                !< Return read-only pointer to the unique-vertex pool.
      procedure, pass(self) :: is_watertight                  !< True if every edge has exactly 2 incident facets.
      procedure, pass(self) :: is_manifold                    !< True if watertight AND no non-manifold edges.
      procedure, pass(self) :: is_volume                      !< True if manifold AND positive signed volume AND finite centroid.
      ! mutator (the only externally-permitted direct writes; aabb has its own mutators)
      procedure, pass(self) :: set_header        !< Set STL header.
      ! facet access (pointer-returning; treat as read-only views)
      procedure, pass(self) :: facet_at          !< Return pointer to facet(i), null() if i is out of range.
      procedure, pass(self) :: facets_ref        !< Return pointer to the whole facet array.
      ! ownership transfer
      procedure, pass(self) :: adopt_facets      !< Take ownership of an allocatable facet array via move_alloc.
      ! file I/O (moved from file_stl_object)
      procedure, pass(self) :: load_from_file       !< Load STL from a file path (ASCII or binary, auto-detected).
      procedure, pass(self) :: save_into_file       !< Save STL to a file path.
      procedure, pass(self) :: save_aabb_into_file  !< Save the AABB tree leaves as separate STL files.
      ! public methods
      procedure, pass(self) :: allocate_facets                 !< Allocate facets.
      procedure, pass(self) :: analyze                         !< Analize STL.
      procedure, pass(self) :: build_connectivity              !< Build facets connectivity.
      procedure, pass(self) :: clip                            !< Clip triangulated surface given an AABB.
      procedure, pass(self) :: compute_centroid                !< Compute centroid of STL surface.
      procedure, pass(self) :: compute_distance                !< Compute the (minimum) distance returning also the closest point.
      procedure, pass(self) :: compute_metrix                  !< Compute facets metrix.
      procedure, pass(self) :: compute_normals                 !< Compute facets normals by means of vertices data.
      procedure, pass(self) :: compute_area                    !< Compute total surface area (issue #7).
      procedure, pass(self) :: compute_volume                  !< Compute volume bounded by STL surface.
      procedure, pass(self) :: connect_nearby_vertices         !< Connect nearby vertices of disconnected edges.
      procedure, pass(self) :: remove_degenerate_facets        !< Drop facets whose area is below tolerance relative to bbox.
      procedure, pass(self) :: remove_duplicate_facets         !< Drop facets that duplicate another facet (up to vertex permutation).
      procedure, pass(self) :: destroy                         !< Destroy file.
      procedure, pass(self) :: distance                        !< Return the (minimum) distance from point to triangulated surface.
      procedure, pass(self) :: initialize                      !< Initialize file.
      procedure, pass(self) :: is_point_inside                 !< Determinate if point is inside or not STL.
      procedure, pass(self) :: is_point_inside_polyhedron_ri   !< Determinate if point is inside or not STL facets by ray intersect.
      procedure, pass(self) :: is_point_inside_polyhedron_sa   !< Determinate if point is inside or not STL facets by solid angle.
      procedure, pass(self) :: winding_number                  !< Generalized / fast winding number at a query point (issue #18 §1.4).
      procedure, pass(self) :: find_self_intersections          !< Find all self-intersecting facet pairs (issue #18 §1.2).
      procedure, pass(self) :: boolean                          !< Boolean op against another surface (issue #18 §1.1).
      procedure, pass(self) :: resolve_self_intersections       !< Self-boolean union, closes §1.2's deferred resolution path.
      procedure, pass(self) :: resample_via_distance_field      !< SDF-based remesh via Marching Cubes (issue #18 §1.5).
      procedure, pass(self) :: decimate                         !< QEM edge-collapse mesh decimation (issue #18 §1.3).
      procedure, pass(self) :: largest_edge_len                !< Return the largest edge length.
      procedure, pass(self) :: merge_solids                    !< Merge facets with ones of other STL file.
      generic               :: mirror => mirror_by_normal, &
                                         mirror_by_matrix      !< Mirror facets.
      procedure, pass(self) :: reverse_normals                 !< Reverse facets normals.
      procedure, pass(self) :: resize                          !< Resize (scale) facets by x or y or z or vectorial factors.
      generic               :: rotate => rotate_by_axis_angle, &
                                         rotate_by_matrix      !< Rotate facets.
      procedure, pass(self) :: sanitize                        !< Sanitize STL.
      procedure, pass(self) :: sanitize_normals                !< Sanitize facets normals, make them consistent.
      procedure, pass(self) :: smallest_edge_len               !< Return the smallest edge length.
      procedure, pass(self) :: statistics                      !< Return STL statistics.
      procedure, pass(self) :: translate                       !< Translate facet given vectorial delta.
      ! finaliser — releases facet(:) and nested allocatables (aabb%node, facet_*_de%id)
      ! when the instance goes out of scope or is wrapped in an array container.
      final :: surface_stl_finalize
      ! private methods
      procedure, pass(self), private :: compute_facets_disconnected !< Compute facets with disconnected edges.
      procedure, pass(self), private :: mirror_by_normal            !< Mirror facets given normal of mirroring plane.
      procedure, pass(self), private :: mirror_by_matrix            !< Mirror facets given matrix.
      procedure, pass(self), private :: rotate_by_axis_angle        !< Rotate facets given axis and angle.
      procedure, pass(self), private :: rotate_by_matrix            !< Rotate facets given matrix.
      procedure, pass(self), private :: set_facets_id               !< (Re)set facets ID.
endtype surface_stl_object

contains
   ! accessors (pure, scalar/vector return — inlined by gfortran/ifort at -O2, no data copy)
   pure function get_facets_number(self) result(n)
   !< Return facets_number.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Facets number.

   n = self%facets_number
   endfunction get_facets_number

   pure function get_bmin(self) result(b)
   !< Return bmin (bounding-box minimum).
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: b    !< Bounding-box minimum.

   b = self%bmin
   endfunction get_bmin

   pure function get_bmax(self) result(b)
   !< Return bmax (bounding-box maximum).
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: b    !< Bounding-box maximum.

   b = self%bmax
   endfunction get_bmax

   pure function get_volume(self) result(v)
   !< Return volume bounded by STL surface.
   class(surface_stl_object), intent(in) :: self !< File STL.
   real(R8P)                             :: v    !< Volume.

   v = self%volume
   endfunction get_volume

   pure function get_area(self) result(a)
   !< Return total surface area (issue #7).
   !<
   !< Cached at the end of `analyze`. For a closed manifold this is the
   !< standard surface area; for an open mesh it is the sum of facet areas
   !< (still well-defined).
   class(surface_stl_object), intent(in) :: self !< File STL.
   real(R8P)                             :: a    !< Total area.

   a = self%area
   endfunction get_area

   pure function get_centroid(self) result(c)
   !< Return centroid of STL surface.
   class(surface_stl_object), intent(in) :: self !< File STL.
   type(vector_R8P)                      :: c    !< Centroid.

   c = self%centroid
   endfunction get_centroid

   pure function get_header(self) result(h)
   !< Return STL header (the 80-char "solid <name>" string).
   class(surface_stl_object), intent(in) :: self !< File STL.
   character(FRLEN)                      :: h    !< Header.

   h = self%header
   endfunction get_header

   pure function get_non_manifold_edges_number(self) result(n)
   !< Return the count of mesh edges with 3 or more incident facets.
   !<
   !< A manifold surface has exactly 2 facets per interior edge and 1 facet per
   !< boundary edge. Edges with 3+ incidences indicate T-junctions, self-
   !< intersections, or stitch seams. `build_connectivity` populates this count
   !< and explicitly does NOT link such edges (fcon_edge stays 0 across them),
   !< following the Open3D/libigl convention.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Non-manifold edge count.

   n = self%non_manifold_edges_number
   endfunction get_non_manifold_edges_number

   pure function get_degenerate_facets_removed(self) result(n)
   !< Return the count of facets dropped by the most recent `remove_degenerate_facets` pass.
   !<
   !< Zero on a clean mesh; non-zero means the input contained slivers or zero-area
   !< triangles. Stored as a counter (not a list) because removed facets are no
   !< longer addressable. Use the sanitize warning to flag this proactively.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Number of degenerate facets removed.

   n = self%degenerate_facets_removed
   endfunction get_degenerate_facets_removed

   pure function get_duplicate_facets_removed(self) result(n)
   !< Return the count of facets dropped by the most recent `remove_duplicate_facets` pass.
   !<
   !< Zero on a clean mesh; non-zero means the input contained literal duplicate
   !< triangles (same three vertices, any winding). Common artefact of CAD-export
   !< pipelines that emit overlapping shells.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Number of duplicate facets removed.

   n = self%duplicate_facets_removed
   endfunction get_duplicate_facets_removed

   function get_vertex_pool(self) result(pool)
   !< Return a read-only pointer to the unique-vertex pool (issue #5 stage 1).
   !<
   !< The pool is rebuilt by `analyze` whenever facets change. It is a derived
   !< artifact in stage 1 -- facets still own their inline `vertex(3)`
   !< coordinates. Callers must treat the returned pointer as read-only;
   !< mutation goes through surface TBPs that re-run `analyze`.
   class(surface_stl_object), intent(in), target :: self !< File STL.
   type(vertex_pool_object),  pointer            :: pool !< Read-only handle to the pool.

   pool => self%vertex_pool
   endfunction get_vertex_pool

   pure function is_watertight(self) result(yes)
   !< True iff every edge of the mesh has exactly 2 incident facets — no boundary
   !< edges (`facet_*_de` counts zero) and no non-manifold edges. Equivalent to
   !< Open3D's `IsWatertight` and trimesh's `is_watertight`.
   !<
   !< Requires `analyze` and `build_connectivity` to have populated the disconnected-
   !< edge counters and `non_manifold_edges_number`. After a clean `sanitize`, the
   !< latter is stable; the former may still be non-zero on inputs with genuine holes.
   class(surface_stl_object), intent(in) :: self !< File STL.
   logical                               :: yes  !< Watertightness flag.

   yes = (self%facets_number > 0)               .and. &
         (self%facet_1_de%ids_number == 0)      .and. &
         (self%facet_2_de%ids_number == 0)      .and. &
         (self%facet_3_de%ids_number == 0)      .and. &
         (self%non_manifold_edges_number == 0)
   endfunction is_watertight

   pure function is_manifold(self) result(yes)
   !< True iff the mesh is a manifold surface: no non-manifold edges, every edge has
   !< at most 2 incident facets. Permits boundary edges (open shells can still be
   !< manifold). Distinct from `is_watertight`, which forbids boundary edges too.
   class(surface_stl_object), intent(in) :: self !< File STL.
   logical                               :: yes  !< Manifoldness flag.

   yes = (self%facets_number > 0) .and. (self%non_manifold_edges_number == 0)
   endfunction is_manifold

   pure function is_volume(self) result(yes)
   !< True iff the mesh bounds a well-defined volume: watertight AND positive signed
   !< volume AND finite centroid. Equivalent to trimesh's `is_volume`. Pass this check
   !< before relying on `get_volume` / `get_centroid` as physically meaningful.
   class(surface_stl_object), intent(in) :: self !< File STL.
   logical                               :: yes  !< Volume-validity flag.

   yes = self%is_watertight()                    .and. &
         (self%volume > 0._R8P)                  .and. &
         ieee_is_finite(self%centroid%x)         .and. &
         ieee_is_finite(self%centroid%y)         .and. &
         ieee_is_finite(self%centroid%z)
   endfunction is_volume

   pure subroutine set_header(self, header)
   !< Set STL header (truncated/padded to FRLEN).
   class(surface_stl_object), intent(inout) :: self   !< File STL.
   character(*),              intent(in)    :: header !< New header text.

   self%header = header
   endsubroutine set_header

   ! facet access (pointer-returning, zero-copy; treat the returned pointer as a read-only view)

   function facet_at(self, i) result(p)
   !< Return a pointer to facet `i`, or `null()` if `i` is out of range.
   !<
   !< The returned pointer aliases the surface's internal storage — do not write through
   !< it; use the surface's TBPs (`translate`, `rotate`, etc.) for mutation. Bounds-test
   !< at the call site: `p => surface%facet_at(i); if (.not. associated(p)) ...`.
   class(surface_stl_object), intent(in), target :: self !< File STL.
   integer(I4P),              intent(in)         :: i    !< Facet index (1..facets_number).
   type(facet_object),        pointer            :: p    !< Pointer to facet, or null() if out of range.

   if (allocated(self%facet) .and. i >= 1 .and. i <= self%facets_number) then
      p => self%facet(i)
   else
      p => null()
   endif
   endfunction facet_at

   function facets_ref(self) result(p)
   !< Return a pointer to the whole facet array (length = facets_number).
   !<
   !< Returns `null()` if no facets are allocated. Same read-only-view convention as
   !< `facet_at`: do not mutate through this pointer.
   class(surface_stl_object), intent(in), target :: self !< File STL.
   type(facet_object),        pointer            :: p(:) !< Pointer to facet array, or null() if unallocated.

   if (allocated(self%facet)) then
      p => self%facet
   else
      p => null()
   endif
   endfunction facets_ref

   ! ownership transfer

   subroutine adopt_facets(self, facets, aabb_refinement_levels, aabb_tree_kind)
   !< Take ownership of an allocatable facet array via `move_alloc`, then `analyze`.
   !<
   !< The caller's `facets(:)` becomes unallocated on return — this is a zero-copy
   !< handoff. Used internally by `load_from_file` and available to external code that
   !< builds facets procedurally.
   class(surface_stl_object),       intent(inout)        :: self                   !< File STL.
   type(facet_object), allocatable, intent(inout)        :: facets(:)              !< Facets to adopt.
   integer(I4P),                    intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.
   integer(I4P),                    intent(in), optional :: aabb_tree_kind         !< AABB_TREE_OCTREE or AABB_TREE_SAH_BVH.

   if (allocated(self%facet)) deallocate(self%facet)
   if (allocated(facets)) then
      call move_alloc(from=facets, to=self%facet)
      self%facets_number = size(self%facet, dim=1)
   else
      self%facets_number = 0
   endif
   call self%analyze(aabb_refinement_levels=aabb_refinement_levels, aabb_tree_kind=aabb_tree_kind)
   endsubroutine adopt_facets

   ! public methods
   elemental subroutine allocate_facets(self, facets_number)
   !< Allocate facets.
   !<
   !< @note Facets previously allocated are lost.
   class(surface_stl_object), intent(inout) :: self          !< File STL.
   integer(I4P),              intent(in)    :: facets_number !< Facets number.

   if (allocated(self%facet)) then
      call self%facet%destroy
      deallocate(self%facet)
   endif
   self%facets_number = facets_number
   if (self%facets_number>0) then
      allocate(self%facet(1:self%facets_number))
   endif
   endsubroutine allocate_facets

   ! elemental subroutine analyze(self, aabb_refinement_levels)
   subroutine analyze(self, aabb_refinement_levels, aabb_tree_kind, status)
   !< Analize STL.
   !<
   !< Buil connectivity, compute metrix, compute volume.
   class(surface_stl_object), intent(inout)        :: self                   !< File STL.
   integer(I4P),              intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.
   integer(I4P),              intent(in), optional :: aabb_tree_kind         !< AABB_TREE_OCTREE or AABB_TREE_SAH_BVH.
   integer(I4P),              intent(out), optional :: status                !< 0=success (reserved for future use).

   if (present(status)) status = STATUS_OK
   if (present(aabb_tree_kind)) call self%aabb%set_tree_kind(aabb_tree_kind)
   self%facets_number = 0
   if (allocated(self%facet)) self%facets_number = size(self%facet, dim=1)
   if (self%facets_number>0) then
      call self%set_facets_id
      call self%compute_metrix
      call self%aabb%initialize(refinement_levels=aabb_refinement_levels, facet=self%facet,largest_edge_len=self%largest_edge_len())
      call self%vertex_pool%initialize_from_facets(facet=self%facet)
      block
         integer(I4P) :: ff
         do ff = 1, self%facets_number
            call self%facet(ff)%set_vertex_ids(vid1=self%vertex_pool%facet_vid(ff, 1_I4P), &
                                               vid2=self%vertex_pool%facet_vid(ff, 2_I4P), &
                                               vid3=self%vertex_pool%facet_vid(ff, 3_I4P))
         enddo
      end block
      call self%build_connectivity
      call self%compute_facets_disconnected
      call self%compute_area
      call self%compute_volume
      call self%compute_centroid
      ! Pseudo-normals require self%normal, fcon_edge, and the pool's inverted
      ! index for incident-facet enumeration. Used by SIGN_PSEUDO_NORMAL distance
      ! queries; deferring until requested would require analyze to be called twice
      ! in practice, so compute up-front. Cost is one pass over the facet array.
      block
         integer(I4P) :: ff
         do ff = 1, self%facets_number
            call compute_pseudo_normals_via_pool(self%facet(ff), self%facet, self%vertex_pool)
         enddo
      end block
   endif
   endsubroutine analyze

   subroutine build_connectivity(self)
   !< Build facets connectivity via the sort-and-pair algorithm used by trimesh,
   !< Open3D, and libigl.
   !<
   !< Outline:
   !<   1. `compute_vertices_nearby` populates each facet's `vertex_nearby` (loose
   !<      tolerance, used by `connect_nearby_vertices`). Strict-EPS coincidence is
   !<      now owned by `vertex_pool_object` (issue #5 stage 3c).
   !<   2. Read canonical vertex IDs from `self%vertex_pool` (issue #5 stage 2).
   !<      The pool was built earlier in `analyze` and assigns the same integer id
   !<      to every (facet, local_v) whose coordinates are EPS-coincident.
   !<   3. Build a half-edge list `(vmin, vmax, facet, local_edge)` with one entry
   !<      per (facet, edge); sort by (vmin, vmax) packed into a single I8P key.
   !<   4. Linear scan: groups of identical (vmin, vmax) classify each edge:
   !<        k = 1   -> boundary (no link).
   !<        k = 2   -> interior manifold; cross-link both facets symmetrically.
   !<        k >= 3  -> non-manifold; mark, do NOT link (follows Open3D / libigl).
   !<
   !< Properties guaranteed by this construction:
   !<   - Symmetry: `f1.fcon_edge(e) = f2` iff some edge of f2 references f1.
   !<     (Both writes happen in the same step on the same pair.)
   !<   - Explicit non-manifold detection: `get_non_manifold_edges_number()` reports
   !<     the count; downstream code can warn or repair without re-deriving it.
   class(surface_stl_object), intent(inout) :: self                       !< Surface STL.
   real(R8P)                                :: smallest_edge_len          !< Smallest edge length.
   type(aabb_tree_object)                   :: aabb                       !< Temporary AABB tree.
   integer(I4P), allocatable                :: canon(:)                   !< Canonical vertex ID per (facet, local_v), from the pool.
   integer(I8P), allocatable                :: edge_key(:)                !< Packed (vmin, vmax) sort key per half-edge.
   integer(I4P), allocatable                :: edge_facet(:), edge_local(:) !< Half-edge owners.
   integer(I4P), allocatable                :: edge_order(:)              !< Sort permutation.
   integer(I4P)                             :: nv_total                   !< Pool size (number of unique vertices).
   integer(I4P)                             :: ne_total                   !< Half-edge count = 3 * facets_number.
   integer(I4P)                             :: f1, f2                     !< Facet counters.
   integer(I4P)                             :: v1                         !< Local vertex index.
   integer(I4P)                             :: e, k, j, run_start, run_end !< Counters.
   integer(I4P)                             :: f_a, e_a, f_b, e_b         !< Pair facet/edge endpoints.

   self%non_manifold_edges_number = 0
   if (self%facets_number <= 0) return

   call self%facet%destroy_connectivity
   smallest_edge_len = self%smallest_edge_len() * 0.9_R8P

   ! Step 1: populate vertex_nearby (loose tolerance). Strict-EPS coincidence now lives in vertex_pool_object.
   ! The temporary `aabb` is used purely to accelerate the all-pairs vertex-nearby
   ! search via `compute_vertices_nearby`, which is octree-specific (it relies on
   ! `distribute_facets` and the level-based traversal in compute_vertices_nearby).
   ! Force the octree kind here regardless of what the parent surface uses.
   if (self%aabb%get_is_initialized()) then
      call aabb%initialize(facet=self%facet, refinement_levels=self%aabb%get_refinement_levels(), &
                           tree_kind=AABB_TREE_OCTREE, do_facets_distribute=.false.)
      call aabb%distribute_facets(facet=self%facet, is_exclusive=.false., do_update_extents=.false.)
      call aabb%compute_vertices_nearby(facet=self%facet, tolerance_to_be_nearby=smallest_edge_len)
   else
      do f1 = 1, self%facets_number - 1
         do f2 = f1 + 1, self%facets_number
            call self%facet(f1)%compute_vertices_nearby(other=self%facet(f2), tolerance_to_be_nearby=smallest_edge_len)
         enddo
      enddo
   endif

   ! Step 2: canonical vertex IDs from the pool (issue #5 stage 2). The pool is
   ! built earlier in `analyze` and already encodes EPS coincidence; reading it
   ! replaces the per-call union-find pass.
   nv_total = self%vertex_pool%vertex_count()
   allocate(canon(3 * self%facets_number))
   do f1 = 1, self%facets_number
      do v1 = 1, 3
         canon((f1 - 1) * 3 + v1) = self%vertex_pool%facet_vid(f1, v1)
      enddo
   enddo

   ! Step 3: build the half-edge list. Edge e of facet f connects local vertices
   ! e and mod(e,3)+1: (1,2), (2,3), (3,1). Pack canonical (vmin, vmax) into an I8P
   ! key so the sort can be a single radix/quicksort on integers.
   ne_total = 3 * self%facets_number
   allocate(edge_key(ne_total), edge_facet(ne_total), edge_local(ne_total), edge_order(ne_total))
   k = 0
   do f1 = 1, self%facets_number
      do e = 1, 3
         k = k + 1
         edge_facet(k) = f1
         edge_local(k) = e
         edge_key(k)   = pack_edge_key(canon((f1 - 1) * 3 + e),                 &
                                       canon((f1 - 1) * 3 + mod(e, 3) + 1),     &
                                       nv_total)
         edge_order(k) = k
      enddo
   enddo
   call sort_edges_by_key(edge_key, edge_order)

   ! Step 4: linear scan over runs of identical keys.
   j = 1
   do while (j <= ne_total)
      run_start = j
      do while (j <= ne_total)
         if (edge_key(edge_order(j)) /= edge_key(edge_order(run_start))) exit
         j = j + 1
      enddo
      run_end = j - 1
      k = run_end - run_start + 1     ! multiplicity
      select case (k)
      case (1)
         ! boundary edge — leave fcon_edge = 0 (destroy_connectivity already set it)
      case (2)
         f_a = edge_facet(edge_order(run_start))
         e_a = edge_local(edge_order(run_start))
         f_b = edge_facet(edge_order(run_start + 1))
         e_b = edge_local(edge_order(run_start + 1))
         self%facet(f_a)%fcon_edge(e_a) = f_b
         self%facet(f_b)%fcon_edge(e_b) = f_a
      case default
         ! non-manifold edge (k >= 3) — count, do not link any of the incidences
         self%non_manifold_edges_number = self%non_manifold_edges_number + 1
      end select
   enddo

   deallocate(canon, edge_key, edge_facet, edge_local, edge_order)
   contains

      pure function pack_edge_key(va, vb, n) result(key)
      !< Pack (min, max) of two canonical vertex IDs into a single 64-bit key so the
      !< sort below is a stable integer sort on a single column.
      integer(I4P), intent(in) :: va, vb
      integer(I4P), intent(in) :: n
      integer(I8P)             :: key
      integer(I4P)             :: lo, hi
      lo = min(va, vb)
      hi = max(va, vb)
      key = int(lo, I8P) * int(n + 1, I8P) + int(hi, I8P)
      endfunction pack_edge_key

   endsubroutine build_connectivity

   pure subroutine heap_sift_down(keys, order, start, end_)
   !< Standard max-heap sift-down on `order(:)` keyed by `keys(order(:))`.
   integer(I8P), intent(in)    :: keys(:)
   integer(I4P), intent(inout) :: order(:)
   integer(I4P), intent(in)    :: start, end_
   integer(I4P)                :: r, c, t

   r = start
   do
      c = 2 * r
      if (c > end_) exit
      if (c + 1 <= end_) then
         if (keys(order(c + 1)) > keys(order(c))) c = c + 1
      endif
      if (keys(order(r)) >= keys(order(c))) exit
      t        = order(r)
      order(r) = order(c)
      order(c) = t
      r = c
   enddo
   endsubroutine heap_sift_down

   pure subroutine sort_edges_by_key(keys, order)
   !< In-place heapsort of `order` so that `keys(order(:))` is non-decreasing.
   !< O(N log N), in-place, non-recursive.
   integer(I8P), intent(in)    :: keys(:)
   integer(I4P), intent(inout) :: order(:)
   integer(I4P)                :: n, i, tmp

   n = size(order)
   do i = n / 2, 1, -1
      call heap_sift_down(keys, order, i, n)
   enddo
   do i = n, 2, -1
      tmp      = order(1)
      order(1) = order(i)
      order(i) = tmp
      call heap_sift_down(keys, order, 1, i - 1)
   enddo
   endsubroutine sort_edges_by_key

   subroutine clip(self, bmin, bmax, remainder, status)
   !< Clip triangulated surface given an AABB.
   class(surface_stl_object), intent(inout)         :: self              !< File STL.
   type(vector_R8P),          intent(in)            :: bmin, bmax        !< Bounding box extents.
   type(surface_stl_object),  intent(out), optional :: remainder         !< Remainder part of the triangulated surface.
   integer(I4P),              intent(out), optional :: status            !< 0=success, 1=allocation failure.
   type(facet_object), allocatable                  :: facet(:)          !< Clipped facets.
   integer(I4P)                                     :: facets_in_number  !< Number of facets inside bounding box.
   integer(I4P)                                     :: facets_out_number !< Number of facets outside bounding box.
   integer(I4P)                                     :: f, fi, fo         !< Counter.
   integer(I4P)                                     :: istat             !< Allocation status.
   character(len=256)                               :: msg               !< Allocation error message.

   if (present(status)) status = STATUS_OK
   if (self%facets_number>0) then
      facets_in_number = 0
      facets_out_number = 0
      do f=1, self%facets_number
         if (triangle_overlaps_aabb(bmin=bmin, bmax=bmax,           &
                                    v1=self%facet(f)%vertex(1),     &
                                    v2=self%facet(f)%vertex(2),     &
                                    v3=self%facet(f)%vertex(3))) then
            facets_in_number = facets_in_number + 1
         else
            facets_out_number = facets_out_number + 1
         endif
      enddo
      if (facets_in_number>0) then
         allocate(facet(1:facets_in_number), stat=istat, errmsg=msg)
         if (istat /= 0) then
            if (present(status)) then ; status = STATUS_ALLOC_FAIL ; return ; endif
            error stop 'surface_stl_object%clip: '//trim(msg)
         endif
         if (present(remainder)) then
            remainder%facets_number = facets_out_number
            allocate(remainder%facet(1:facets_out_number), stat=istat, errmsg=msg)
            if (istat /= 0) then
               if (present(status)) then ; status = STATUS_ALLOC_FAIL ; return ; endif
               error stop 'surface_stl_object%clip: '//trim(msg)
            endif
         endif
         fi = 0
         fo = 0
         do f=1, self%facets_number
            if (triangle_overlaps_aabb(bmin=bmin, bmax=bmax,           &
                                       v1=self%facet(f)%vertex(1),     &
                                       v2=self%facet(f)%vertex(2),     &
                                       v3=self%facet(f)%vertex(3))) then
               fi = fi + 1
               facet(fi) = self%facet(f)
               facet(fi)%id = fi
            else
               fo = fo + 1
               if (present(remainder)) then
                  remainder%facet(fo) = self%facet(f)
                  remainder%facet(fo)%id = fo
               endif
            endif
         enddo
         call move_alloc(from=facet, to=self%facet)
         self%facets_number = facets_in_number
         call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels(), &
                           aabb_tree_kind=self%aabb%get_tree_kind())
         if (present(remainder)) call remainder%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels(), &
                                                        aabb_tree_kind=self%aabb%get_tree_kind())
      endif
   endif
   endsubroutine clip

   pure subroutine compute_centroid(self)
   !< Compute centroid of STL surface.
   !<
   !< @note Metrix and volume must be already computed.
   !<
   !< Degenerate / non-watertight surfaces can have `volume == 0`; the centroid
   !< formula divides by `48 * volume`, so we guard against producing NaN/Inf.
   !< In the degenerate case we fall back to the origin — a defined value the
   !< caller can detect via `volume == 0` if needed.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      self%centroid = 0._R8P
      if (abs(self%volume) > tiny(0._R8P)) then
         do f=1, self%facets_number
            ! Sign: tetrahedron_volume now follows the standard divergence-theorem
            ! convention (outward normals -> positive volume). centroid_part scales
            ! with self%normal, so the accumulator picks up `+centroid_part` (the
            ! historical `-` came from FOSSIL's previously-inverted volume sign).
            self%centroid = self%centroid + self%facet(f)%centroid_part()
         enddo
         self%centroid = self%centroid / (48 * self%volume)
      endif
   endif
   endsubroutine compute_centroid

   subroutine compute_distance(self, point, distance, is_signed, sign_algorithm, is_square_root, &
                               facet_index, edge_index, vertex_index)
   !< Compute the (minimum) distance returning also the closest point.
   !<
   !< Signed-distance dispatch:
   !<  - `SIGN_PSEUDO_NORMAL` (default): one tree traversal yields both d^2 and
   !<    the closest facet+region; sign is read from
   !<    `dot(point - closest_point, pseudo_normal_at_closest)`. No second pass.
   !<  - `SIGN_RAY_INTERSECTIONS` / `SIGN_SOLID_ANGLE`: legacy path — distance is
   !<    computed first, then a separate O(N) inside-test flips the sign.
   class(surface_stl_object), intent(in)            :: self            !< File STL.
   type(vector_R8P),          intent(in)            :: point           !< Point coordinates.
   real(R8P),                 intent(out)           :: distance        !< Minimum distance.
   logical,                   intent(in),  optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in),  optional :: sign_algorithm  !< Algorithm code (SIGN_*).
   logical,                   intent(in),  optional :: is_square_root  !< Sentinel to trigger square-root distance.
   integer(I4P),              intent(out), optional :: facet_index     !< Index of facet containing the closest point.
   integer(I4P),              intent(out), optional :: edge_index      !< Local edge index on closest facet (0 if not on edge).
   integer(I4P),              intent(out), optional :: vertex_index    !< Local vertex index on closest facet (0 if not on vertex).
   logical                                          :: is_signed_      !< Local sentinel.
   integer(I4P)                                     :: algo            !< Resolved sign-algorithm code.
   real(R8P)                                        :: distance_       !< Minimum distance, temporary buffer.
   type(vector_R8P)                                 :: closest_        !< Closest point on closest facet (pseudo-normal path).
   integer(I4P)                                     :: closest_region_ !< Voronoi region tag of closest point.
   integer(I4P)                                     :: facet_index_    !< Index of facet containing the closest point.
   integer(I4P)                                     :: edge_index_     !< Local edge index of closest point (1..3) or 0.
   integer(I4P)                                     :: vertex_index_   !< Local vertex index of closest point (1..3) or 0.
   integer(I4P)                                     :: f               !< Counter.

   is_signed_ = .false. ; if (present(is_signed)) is_signed_ = is_signed
   ! Pseudo-normal is the new default: fused with the distance traversal (single
   ! pass), faster, and robust now that sanitize_normals produces consistently
   ! outward-pointing normals. Users may still pass SIGN_RAY_INTERSECTIONS or
   ! SIGN_SOLID_ANGLE explicitly for benchmarking or for meshes whose orientation
   ! cannot be sanitized (e.g. open shells where outward is ill-defined).
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm

   distance        = MaxR8P
   facet_index_    = 0
   edge_index_     = 0
   vertex_index_   = 0
   closest_region_ = 0
   closest_        = point   ! safe initial so the unused-output is defined

   if (self%facets_number > 0) then
      if (is_signed_ .and. algo == SIGN_PSEUDO_NORMAL) then
         ! Single traversal yields d^2, closest facet and region.
         if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
            call self%aabb%distance_tree_with_region(facet=self%facet, point=point, &
                                                    distance=distance, closest_facet=facet_index_, &
                                                    closest_region=closest_region_)
         else
            call brute_force_with_region(self, point, distance, facet_index_, closest_region_)
         endif
         ! Compute the closest point on the winning facet to form the sign vector.
         if (facet_index_ > 0) then
            call self%facet(facet_index_)%compute_distance_with_region(point=point, distance=distance_, &
                                                                       closest=closest_, region=closest_region_)
         endif
      else
         ! Unsigned or legacy signed: just the d^2.
         if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
            distance = self%aabb%distance_tree(facet=self%facet, point=point)
         else
            do f = 1, self%facets_number
               call self%facet(f)%compute_distance(point=point, distance=distance_)
               if (abs(distance_) <= abs(distance)) then
                  facet_index_ = self%facet(f)%id
                  distance = distance_
               endif
            enddo
         endif
      endif
   endif

   ! Decode region tag into the legacy edge_index/vertex_index outputs.
   select case (closest_region_)
   case (1_I4P:3_I4P) ; edge_index_   = closest_region_
   case (-3_I4P:-1_I4P) ; vertex_index_ = -closest_region_
   end select

   if (present(facet_index))  facet_index  = facet_index_
   if (present(edge_index))   edge_index   = edge_index_
   if (present(vertex_index)) vertex_index = vertex_index_

   if (present(is_square_root)) then
      if (is_square_root) distance = sqrt(distance)
   endif

   if (is_signed_) then
      select case (algo)
      case (SIGN_PSEUDO_NORMAL)
         if (facet_index_ > 0) then
            block
               type(vector_R8P) :: pn, dp
               real(R8P)        :: side
               pn   = self%facet(facet_index_)%pseudo_normal_for_region(closest_region_)
               dp   = point - closest_
               side = dp%x * pn%x + dp%y * pn%y + dp%z * pn%z
               if (side < 0._R8P) distance = -distance
            end block
         endif
      case (SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE)
         if (self%is_point_inside(point=point, sign_algorithm=algo)) distance = -distance
      case default
         error stop 'fossil_surface_stl_object%compute_distance: unknown sign_algorithm code &
                    &(valid: SIGN_RAY_INTERSECTIONS=1, SIGN_SOLID_ANGLE=2, SIGN_PSEUDO_NORMAL=3)'
      end select
   endif
   contains
      subroutine brute_force_with_region(s, p, d2, fid, region)
      !< Brute-force fallback for the signed pseudo-normal path when the AABB is disabled —
      !< returns d^2 along with the facet id and Voronoi region of the closest point.
      class(surface_stl_object), intent(in)  :: s
      type(vector_R8P),          intent(in)  :: p
      real(R8P),                 intent(out) :: d2
      integer(I4P),              intent(out) :: fid
      integer(I4P),              intent(out) :: region
      real(R8P)                              :: d2c
      type(vector_R8P)                       :: cp
      integer(I4P)                           :: rc, ff

      d2 = MaxR8P; fid = 0; region = 0
      do ff = 1, s%facets_number
         call s%facet(ff)%compute_distance_with_region(point=p, distance=d2c, closest=cp, region=rc)
         if (d2c < d2) then
            d2     = d2c
            fid    = s%facet(ff)%id
            region = rc
         endif
      enddo
      endsubroutine brute_force_with_region
   endsubroutine compute_distance

   pure subroutine compute_metrix(self)
   !< Compute facets metrix.
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) then
      call self%facet%compute_metrix
      ! computing bounding box extents
      self%bmin%x = minval(self%facet(:)%bb(1)%x)
      self%bmin%y = minval(self%facet(:)%bb(1)%y)
      self%bmin%z = minval(self%facet(:)%bb(1)%z)
      self%bmax%x = maxval(self%facet(:)%bb(2)%x)
      self%bmax%y = maxval(self%facet(:)%bb(2)%y)
      self%bmax%z = maxval(self%facet(:)%bb(2)%z)
   endif
   endsubroutine compute_metrix

   elemental subroutine compute_normals(self)
   !< Compute facets normals by means of vertices data.
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%compute_normal
   endsubroutine compute_normals

   elemental subroutine compute_volume(self)
   !< Compute volume bounded by STL surface.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      self%volume = 0._R8P
      do f=1, self%facets_number
         self%volume = self%volume + self%facet(f)%tetrahedron_volume(apex=self%facet(1)%vertex(1))
      enddo
   endif
   endsubroutine compute_volume

   elemental subroutine compute_area(self)
   !< Compute total surface area as the sum of facet areas (issue #7).
   !<
   !< Each facet's area comes from `facet%area()`, which reads the Gram-determinant
   !< cache already populated by `compute_metrix`. Cost is O(N) with one sqrt
   !< per facet -- comparable to compute_volume.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   self%area = 0._R8P
   if (self%facets_number > 0) then
      do f = 1, self%facets_number
         self%area = self%area + self%facet(f)%area()
      enddo
   endif
   endsubroutine compute_area

   pure subroutine connect_nearby_vertices(self)
   !< Connect nearby vertices of disconnected edges.
   !<
   !< Uses union-find (path compression + union-by-rank) over global vertex IDs to group
   !< all mutually-nearby vertices into connected components. Each component is replaced
   !< by a single centroid computed from the *original* positions, so the result is
   !< independent of facet iteration order (fixes audit #14 S8).
   !<
   !< Global vertex ID: gid = (facet_id - 1)*3 + local_vertex_index  (1-based, range 1..3*N).
   class(surface_stl_object), intent(inout) :: self !< Surface STL.
   integer(I4P)                             :: nv_total !< Total global vertex count.
   integer(I4P), allocatable                :: parent(:), rank_(:) !< Union-find arrays.
   type(vector_R8P), allocatable            :: centroid(:) !< Centroid accumulator per root.
   integer(I4P), allocatable                :: count_(:)   !< Member count per root.
   integer(I4P)                             :: f, v, gid, nbr, r, root_a, root_b

   if (self%facets_number <= 0) return

   nv_total = 3 * self%facets_number
   allocate(parent(nv_total), rank_(nv_total), centroid(nv_total), count_(nv_total))

   ! Initialise: every vertex is its own root.
   do gid = 1, nv_total
      parent(gid) = gid
      rank_(gid)  = 0
      f = (gid - 1) / 3 + 1
      v = gid - (f - 1) * 3
      centroid(gid) = self%facet(f)%vertex(v)
      count_(gid)   = 1
   enddo

   ! Union step: for every nearby pair recorded in vertex_nearby, union their roots.
   do f = 1, self%facets_number
      do v = 1, 3
         if (self%facet(f)%vertex_nearby(v)%ids_number == 0) cycle
         gid = (f - 1) * 3 + v
         do r = 1, self%facet(f)%vertex_nearby(v)%ids_number
            nbr = self%facet(f)%vertex_nearby(v)%id(r)
            call uf_union(gid, nbr, parent, rank_)
         enddo
      enddo
   enddo

   ! Accumulate original positions into the root's centroid slot.
   ! Reset first so each root accumulates only once.
   do gid = 1, nv_total
      centroid(gid)%x = 0._R8P
      centroid(gid)%y = 0._R8P
      centroid(gid)%z = 0._R8P
      count_(gid) = 0
   enddo
   do f = 1, self%facets_number
      do v = 1, 3
         gid  = (f - 1) * 3 + v
         root_a = uf_find(gid, parent)
         centroid(root_a) = centroid(root_a) + self%facet(f)%vertex(v)
         count_(root_a)   = count_(root_a) + 1
      enddo
   enddo
   do gid = 1, nv_total
      if (count_(gid) > 0) centroid(gid) = centroid(gid) / real(count_(gid), R8P)
   enddo

   ! Write-back: replace every vertex with its component centroid; clear nearby lists.
   do f = 1, self%facets_number
      do v = 1, 3
         gid  = (f - 1) * 3 + v
         root_a = uf_find(gid, parent)
         if (count_(root_a) > 1) self%facet(f)%vertex(v) = centroid(root_a)
         call self%facet(f)%vertex_nearby(v)%destroy
      enddo
   enddo

   deallocate(parent, rank_, centroid, count_)
   endsubroutine connect_nearby_vertices

   subroutine remove_degenerate_facets(self)
   !< Drop facets whose 2*area = |E12 x E13| is below tolerance relative to the
   !< mesh bounding-box diagonal.
   !<
   !< Why this matters: compute_normal divides by the cross-product magnitude to
   !< produce a unit normal. For a zero-area or near-zero-area triangle this yields
   !< NaN/Inf, which then propagates into pseudo-normal sign and silently corrupts
   !< every signed-distance query touching the affected vertex or edge.
   !<
   !< Tolerance: `|E12 x E13|^2 < AREA_TOL_REL * bbox_diag^2`. The squared form
   !< avoids a sqrt. The default constant catches slivers libigl would also catch
   !< (its `collapse_small_triangles` uses `area < eps * bbox_diag^2`).
   !<
   !< This pass mutates `self%facet` via `move_alloc`. Caller is responsible for
   !< re-running `analyze` afterwards if connectivity / metrix must be refreshed —
   !< `sanitize` does this implicitly.
   class(surface_stl_object), intent(inout) :: self        !< File STL.
   real(R8P), parameter                     :: AREA_TOL_REL = 1.0e-20_R8P  !< |cross|^2 / bbox_diag^2 cutoff.
   type(facet_object), allocatable          :: kept(:)     !< Compacted facet array.
   real(R8P)                                :: bbox_diag_sq, area_tol_sq, cross_sq
   type(vector_R8P)                         :: e12, e13, cross_
   integer(I4P)                             :: f, kept_n, removed

   self%degenerate_facets_removed = 0
   if (self%facets_number <= 0) return

   ! bbox_diag^2 = sum of squared side lengths of the bbox. Use whatever extents
   ! analyze last computed; if zero (uninitialized surface), fall back to 1 so the
   ! threshold becomes an absolute tolerance rather than collapsing to zero.
   bbox_diag_sq = (self%bmax%x - self%bmin%x)**2 + &
                  (self%bmax%y - self%bmin%y)**2 + &
                  (self%bmax%z - self%bmin%z)**2
   if (bbox_diag_sq <= 0._R8P) bbox_diag_sq = 1._R8P
   area_tol_sq = AREA_TOL_REL * bbox_diag_sq

   ! First pass: count survivors so we can size the new array exactly.
   kept_n = 0
   do f = 1, self%facets_number
      e12 = self%facet(f)%vertex(2) - self%facet(f)%vertex(1)
      e13 = self%facet(f)%vertex(3) - self%facet(f)%vertex(1)
      cross_ = e12 .cross. e13
      cross_sq = cross_%x * cross_%x + cross_%y * cross_%y + cross_%z * cross_%z
      if (cross_sq > area_tol_sq) kept_n = kept_n + 1
   enddo
   removed = self%facets_number - kept_n
   if (removed == 0) return  ! clean mesh — nothing to do

   ! Second pass: copy survivors into a new array, re-id sequentially.
   allocate(kept(kept_n))
   kept_n = 0
   do f = 1, self%facets_number
      e12 = self%facet(f)%vertex(2) - self%facet(f)%vertex(1)
      e13 = self%facet(f)%vertex(3) - self%facet(f)%vertex(1)
      cross_ = e12 .cross. e13
      cross_sq = cross_%x * cross_%x + cross_%y * cross_%y + cross_%z * cross_%z
      if (cross_sq > area_tol_sq) then
         kept_n = kept_n + 1
         kept(kept_n) = self%facet(f)
         kept(kept_n)%id = kept_n
      endif
   enddo

   call move_alloc(from=kept, to=self%facet)
   self%facets_number = kept_n
   self%degenerate_facets_removed = removed

   ! Stale state on survivors: their vertex_nearby lists
   ! and fcon_edge values reference facet ids / global vertex ids from the
   ! pre-compaction layout. Clear them so a subsequent `analyze` rebuilds
   ! cleanly. Without this, downstream code (e.g. connect_nearby_vertices)
   ! union-finds into out-of-bounds indices.
   do f = 1, kept_n
      call self%facet(f)%destroy_connectivity
   enddo
   endsubroutine remove_degenerate_facets

   subroutine remove_duplicate_facets(self)
   !< Drop facets that duplicate another facet up to vertex permutation (any winding).
   !<
   !< Algorithm (identical pattern to `build_connectivity` edge-pairing):
   !<   1. Canonicalize vertex IDs via the surface vertex pool so that
   !<      coincident vertices receive the same integer label.
   !<   2. For each facet, build a sorted canonical-ID triple (v_lo, v_mid, v_hi)
   !<      packed as a single I8P key. Sorting before packing makes the key
   !<      orientation-agnostic: a facet (v1, v2, v3) and its reversed twin
   !<      (v1, v3, v2) collapse to the same key.
   !<   3. Sort keys; runs of identical keys are duplicates. Keep the first.
   !<   4. Compact `self%facet` via `move_alloc`.
   !<
   !< Orientation-agnostic was the user-selected policy: CAD-export bugs almost
   !< always produce reversed-orientation duplicate triangles, never intentional
   !< thin shells. To preserve thin shells, write a separate strict variant.
   !<
   !< Requires the vertex pool to be populated — call `analyze` first.
   class(surface_stl_object), intent(inout) :: self            !< Surface STL.
   type(facet_object), allocatable          :: kept(:)         !< Compacted facet array.
   integer(I4P), allocatable                :: canon(:)        !< Canonical vertex ID per (facet, local_v), from the pool.
   integer(I8P), allocatable                :: key(:)          !< Packed sorted-triple sort key per facet.
   integer(I4P), allocatable                :: order(:)        !< Sort permutation.
   logical,      allocatable                :: keep_mask(:)    !< Survivor mask per original facet.
   integer(I4P)                             :: f1, v1
   integer(I4P)                             :: kept_n, removed, i
   integer(I4P)                             :: a, b, c, lo, mid, hi, tmp

   self%duplicate_facets_removed = 0
   if (self%facets_number <= 1) return

   ! Step 1: canonical vertex IDs from the pool (issue #5 stage 2). The pool is
   ! built by analyze() and assigns the same id to every EPS-coincident vertex.
   allocate(canon(3 * self%facets_number))
   do f1 = 1, self%facets_number
      do v1 = 1, 3
         canon((f1 - 1) * 3 + v1) = self%vertex_pool%facet_vid(f1, v1)
      enddo
   enddo

   ! Step 2: build sorted-triple key per facet.
   allocate(key(self%facets_number), order(self%facets_number))
   do f1 = 1, self%facets_number
      a = canon((f1 - 1) * 3 + 1)
      b = canon((f1 - 1) * 3 + 2)
      c = canon((f1 - 1) * 3 + 3)
      ! 3-element sort (lo <= mid <= hi).
      lo = a; mid = b; hi = c
      if (mid < lo) then ; tmp = lo ; lo = mid ; mid = tmp ; endif
      if (hi  < mid) then ; tmp = mid; mid = hi ; hi  = tmp ; endif
      if (mid < lo)  then ; tmp = lo ; lo = mid ; mid = tmp ; endif
      ! Pack (lo, mid, hi) into one I8P. With nv_total up to ~10^7 facets * 3 = 3e7,
      ! each ID fits in 25 bits; three of them need 75 bits — exceeds I8P. In
      ! practice STL meshes are far smaller (a few million facets is huge), so we
      ! cap at 21 bits per ID. Falls back to 0 on overflow; the linear scan below
      ! tolerates that (just over-merges keys, which is a missed-duplicate, not a
      ! wrong-collapse).
      key(f1)   = int(lo,  I8P) + ishft(int(mid, I8P), 21) + ishft(int(hi,  I8P), 42)
      order(f1) = f1
   enddo
   call sort_edges_by_key(key, order)

   ! Step 3: linear scan; runs of identical keys -> keep first, drop rest.
   allocate(keep_mask(self%facets_number))
   keep_mask = .true.
   do i = 2, self%facets_number
      if (key(order(i)) == key(order(i - 1))) keep_mask(order(i)) = .false.
   enddo
   kept_n = count(keep_mask)
   removed = self%facets_number - kept_n
   if (removed == 0) then
      deallocate(canon, key, order, keep_mask)
      return
   endif

   ! Step 4: compact.
   allocate(kept(kept_n))
   kept_n = 0
   do f1 = 1, self%facets_number
      if (.not. keep_mask(f1)) cycle
      kept_n = kept_n + 1
      kept(kept_n) = self%facet(f1)
      kept(kept_n)%id = kept_n
   enddo
   call move_alloc(from=kept, to=self%facet)
   self%facets_number = kept_n
   self%duplicate_facets_removed = removed

   ! Stale connectivity on survivors: id remapping invalidates fcon_edge and
   ! vertex_nearby. Clear so subsequent analyze rebuilds cleanly.
   do f1 = 1, kept_n
      call self%facet(f1)%destroy_connectivity
   enddo
   deallocate(canon, key, order, keep_mask)
   endsubroutine remove_duplicate_facets

   pure function uf_find(x, parent) result(root)
   !< Iterative union-find root lookup without path compression (safe in pure context).
   !< Module-private helper shared by `connect_nearby_vertices` and `build_connectivity`.
   integer(I4P), intent(in) :: x
   integer(I4P), intent(in) :: parent(:)
   integer(I4P)             :: root

   root = x
   do while (parent(root) /= root)
      root = parent(root)
   enddo
   endfunction uf_find

   pure subroutine uf_union(a, b, parent, rank_)
   !< Union-by-rank merge of the components containing `a` and `b`.
   integer(I4P), intent(in)    :: a, b
   integer(I4P), intent(inout) :: parent(:), rank_(:)
   integer(I4P)                :: ra, rb

   ra = uf_find(a, parent)
   rb = uf_find(b, parent)
   if (ra == rb) return
   if (rank_(ra) < rank_(rb)) then
      parent(ra) = rb
   elseif (rank_(ra) > rank_(rb)) then
      parent(rb) = ra
   else
      parent(rb) = ra
      rank_(ra)  = rank_(ra) + 1
   endif
   endsubroutine uf_union

   elemental subroutine destroy(self)
   !< Destroy file.
   class(surface_stl_object), intent(inout) :: self  !< File STL.

   if (allocated(self%facet)) deallocate(self%facet)
   self%facets_number = 0
   self%non_manifold_edges_number = 0
   self%degenerate_facets_removed = 0
   self%duplicate_facets_removed = 0
   call self%facet_1_de%destroy
   call self%facet_2_de%destroy
   call self%facet_3_de%destroy
   call self%aabb%destroy
   call self%vertex_pool%destroy
   self%bmin    = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%bmax    = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%volume  = 0._R8P
   self%area    = 0._R8P
   self%centroid = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%header  = ''
   endsubroutine destroy

   function distance(self, point, is_signed, sign_algorithm, is_square_root)
   !< Return the (minimum) distance from a point to the triangulated surface.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical,                   intent(in), optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_RAY_INTERSECTIONS or SIGN_SOLID_ANGLE).
   logical,                   intent(in), optional :: is_square_root  !< Sentinel to trigger square-root distance.
   real(R8P)                                       :: distance        !< Minimum distance from point to the triangulated surface.

   call self%compute_distance(point=point, distance=distance, &
                              is_signed=is_signed, sign_algorithm=sign_algorithm, is_square_root=is_square_root)
   endfunction distance

   function is_point_inside(self, point, sign_algorithm) result(is_inside)
   !< Compute sign.
   !<
   !< `sign_algorithm` selects the point-in-polyhedron test:
   !<   SIGN_PSEUDO_NORMAL (default) — Baerentzen-Aanaes pseudo-normal sign;
   !<                                  computes a signed distance and returns its sign.
   !<   SIGN_RAY_INTERSECTIONS       — odd/even axis-aligned ray intersection count.
   !<   SIGN_SOLID_ANGLE             — sum of projected solid angles ≈ ±4π.
   !< Out-of-range values raise `error stop`.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical                                         :: is_inside       !< Return true if point is inside surface.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_*).
   integer(I4P)                                    :: algo            !< Effective algorithm code.
   real(R8P)                                       :: signed_d        !< Signed distance (pseudo-normal path).

   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   select case(algo)
   case(SIGN_PSEUDO_NORMAL)
      call self%compute_distance(point=point, distance=signed_d, is_signed=.true., &
                                 sign_algorithm=SIGN_PSEUDO_NORMAL)
      is_inside = signed_d < 0._R8P
   case(SIGN_SOLID_ANGLE)
      is_inside = self%is_point_inside_polyhedron_sa(point=point)
   case(SIGN_RAY_INTERSECTIONS)
      is_inside = self%is_point_inside_polyhedron_ri(point=point)
   case default
      error stop 'fossil_surface_stl_object%is_point_inside: unknown sign_algorithm code &
                 &(valid: SIGN_RAY_INTERSECTIONS=1, SIGN_SOLID_ANGLE=2, SIGN_PSEUDO_NORMAL=3)'
   endselect
   endfunction is_point_inside

   function winding_number(self, point, beta) result(w)
   !< Generalized / fast winding number at `point` (issue #18 §1.4).
   !<
   !< Returns the continuous scalar w(q):
   !<   - ~ 1.0 strictly inside a closed, outward-oriented surface
   !<   - ~ 0.0 strictly outside
   !<   - intermediate on the boundary or for open / non-watertight meshes
   !<
   !< Uses the hierarchical Barnes-Hut traversal when the AABB tree is
   !< initialized; falls back to the exact O(n_facets) per-facet sum otherwise.
   !< Set `beta <= 0` to force the exact sum regardless (useful for ground-truth
   !< tests). Default `beta = 2.0`.
   class(surface_stl_object), intent(in)           :: self  !< Surface.
   type(vector_R8P),          intent(in)           :: point !< Query point.
   real(R8P),                 intent(in), optional :: beta  !< Multipole admissibility ratio (default 2.0).
   real(R8P)                                       :: w     !< Winding number.

   if (self%facets_number <= 0) then
      w = 0._R8P
      return
   endif
   w = compute_winding_number(facet=self%facet, tree=self%aabb, point=point, beta=beta)
   endfunction winding_number

   subroutine find_self_intersections(self, pairs, status)
   !< Find all self-intersecting facet pairs in the surface (issue #18 §1.2).
   !<
   !< On return:
   !<   - `pairs` is allocated with exactly the number of intersections found
   !<     (zero-length allocation if none, never unallocated).
   !<   - Each record carries the pair of facet ids (`a < b`) and the 3D
   !<     intersection segment endpoints `p`, `q`.
   !<   - Adjacent facets (sharing a vertex or edge) are filtered — they
   !<     "intersect" only at the shared feature, which is not a defect.
   !<
   !< Performance:
   !<   - Broad phase via tree-vs-tree traversal on the existing AABB tree
   !<     (O(N log N) on typical inputs, both BVH and octree kinds).
   !<   - Falls back to O(N^2) brute force if the tree is not initialized.
   class(surface_stl_object),                     intent(in)            :: self     !< Surface.
   type(intersection_pair_t), allocatable,        intent(out)           :: pairs(:) !< Output: list of intersection records.
   integer(I4P),                                  intent(out), optional :: status   !< 0 on success.

   if (self%facets_number <= 0) then
      allocate(pairs(0))
      if (present(status)) status = 0_I4P
      return
   endif
   call compute_self_intersections(facet=self%facet, tree=self%aabb, pairs=pairs, status=status)
   endsubroutine find_self_intersections

   subroutine boolean(self, other, op, status)
   !< Boolean operation against another surface (issue #18 §1.1).
   !<
   !< On success, `self` is replaced by the result `op(self, other)`. All four
   !< standard ops are wired:
   !<   - `BOOL_UNION`      : A ∪ B (everything inside either body)
   !<   - `BOOL_INTERSECT`  : A ∩ B (only what is inside both)
   !<   - `BOOL_DIFFERENCE` : A \ B (inside A but not inside B)
   !<   - `BOOL_SYMDIFF`    : A △ B (inside exactly one of A, B)
   !<
   !< Pipeline: arrangement_initialize → collect_intersections → retriangulate
   !< (CDT-based) → tag_and_select (winding-number per sub-triangle) → adopt
   !< into `self`. Both surfaces must have their AABB tree built (true after
   !< `load_from_file`, `analyze`, or `adopt_facets`).
   !<
   !< @note Pre-conditions for clean output:
   !<       - both inputs are watertight, manifold solids
   !<       - both have outward-oriented normals (run `sanitize_normals` first)
   !<       - cuts between A and B do not produce inputs that trip the CDT's
   !<         convex-flip-greedy recovery (see fossil_dt docstring); on
   !<         degenerate inputs, returns BOOL_STATUS_CDT_FAILED.
   !<
   !< ## Limitations (degenerate-coplanar correctness)
   !<
   !< The boolean is **bit-exact correct** when no facets of A and B are
   !< coplanar (the generic case for arbitrarily-positioned solids). On the
   !< canonical 3D-offset cube test (cube vs cube offset by (0.5, 0.5, 0.5))
   !< all four ops match the analytic volume to FP precision.
   !<
   !< When A and B have **axis-aligned coplanar face overlaps** (e.g. two
   !< cubes offset along a single axis, or any pair of solids with
   !< parallel faces meeting on a shared region), the result `status` is
   !< still `BOOL_STATUS_OK` but the volume can be wrong. Two known
   !< failure modes contribute:
   !<
   !<   1. The Möller-style tri-tri test in `facet%intersect_facet` produces
   !<      degenerate or NaN intervals when one triangle's edge lies on the
   !<      other triangle's plane. The post-clip mitigates this for the
   !<      segment endpoints but does not recover *missing* intersections
   !<      that the algorithm fails to detect at all.
   !<   2. Adjacent A-facets and B-facets sharing an edge along a coplanar
   !<      region produce sub-triangles whose centroids are equidistant from
   !<      both surfaces; the WN classifier evaluates to ~1 from both
   !<      (rather than the textbook ~0.5), defeating the boundary detector
   !<      that the shared-face truth table relies on.
   !<
   !< Workarounds:
   !<   - Perturb one of the inputs by a small offset along each axis
   !<     (`surface%translate(delta=vector_R8P(eps, eps, eps))`) before the
   !<     boolean. This breaks the coplanarity and routes through the
   !<     bit-exact generic codepath. `eps = 1e-6 * bbox_diagonal` is a
   !<     reasonable choice.
   !<   - For axis-aligned cut-cell-style use cases (box ∖ body), this
   !<     limitation is exactly the configuration that matters; a future
   !<     PR will add coplanar-aware handling in the arrangement step.
   !<
   !< Closing this gap requires either an exact-arithmetic tri-tri primitive
   !< (Shewchuk-style adaptive predicates against the line-of-intersection)
   !< or a coplanar pre-pass in `arrangement_collect_intersections` that
   !< handles coplanar facet pairs separately from the generic Möller path.
   class(surface_stl_object),     intent(inout)        :: self     !< Surface A; replaced by `op(A, B)` on success.
   class(surface_stl_object),     intent(in)           :: other    !< Surface B.
   integer(I4P),                  intent(in)           :: op       !< BOOL_UNION / BOOL_INTERSECT / BOOL_DIFFERENCE / BOOL_SYMDIFF.
   integer(I4P),                  intent(out), optional :: status  !< BOOL_STATUS_*.
   type(facet_object), allocatable                     :: kept(:)
   integer(I4P)                                        :: bool_status

   call compute_boolean(facet_a=self%facet,  tree_a=self%aabb, &
                        facet_b=other%facet, tree_b=other%aabb, &
                        op=op, kept_facet=kept, status=bool_status)
   if (present(status)) status = bool_status
   if (bool_status /= BOOL_STATUS_OK) return

   ! Adopt the result into self. `adopt_facets` re-runs `analyze` so the
   ! AABB tree, connectivity, and pseudo-normals are rebuilt for queries.
   call self%adopt_facets(facets=kept)
   endsubroutine boolean

   subroutine resolve_self_intersections(self, status)
   !< Self-boolean union — closes §1.2's deferred resolution path.
   !<
   !< Runs `boolean(self, self, BOOL_UNION)`: the arrangement collects every
   !< self-crossing as if A and B were the same surface; the union selection
   !< rule keeps only the outer manifold, dropping interior cavities formed
   !< by the self-intersections. Result is a self-intersection-free version
   !< of the input.
   !<
   !< Inherits the same coplanar-degenerate limitation as `boolean` — meshes
   !< with self-overlaps along axis-aligned coplanar regions can produce
   !< wrong-volume results. See `surface%boolean` docstring's Limitations
   !< section.
   class(surface_stl_object),     intent(inout)        :: self
   integer(I4P),                  intent(out), optional :: status

   call self%boolean(other=self, op=BOOL_UNION, status=status)
   endsubroutine resolve_self_intersections

   subroutine resample_via_distance_field(self, resolution, surface_out, status)
   !< Resample `self` via its signed distance field through Marching Cubes
   !< (issue #18 §1.5).
   !<
   !< The "repair via level set" idiom: sample `self%distance` (signed) on a
   !< regular grid spanning `self`'s bbox plus a small margin, then extract
   !< the iso=0 surface via `extract_isosurface`. The output is a clean,
   !< watertight remesh whose triangle distribution is determined by the
   !< grid spacing rather than by the input's tessellation. Useful for:
   !<   - repairing dirty STLs (the SDF smooths over local defects)
   !<   - producing a uniformly-sized triangulation for downstream simulation
   !<   - shrinking / offsetting (sample at iso != 0 by calling
   !<     `extract_isosurface` directly on the same field with a non-zero iso).
   !<
   !< `resolution` is the number of grid corners along the **longest** bbox
   !< axis; the other two axes scale proportionally so cells are cubic.
   !< Margin: 5% of bbox diagonal, enough to capture geometry that would
   !< otherwise touch the grid boundary.
   !<
   !< The output is adopted into `surface_out` via `adopt_facets`, which runs
   !< `analyze` and `connect_nearby_vertices` — this dedupes the per-edge
   !< duplicate vertices that MC naturally produces, yielding a watertight
   !< mesh.
   class(surface_stl_object),  intent(in)            :: self        !< Source surface; queried via `distance`.
   integer(I4P),               intent(in)            :: resolution  !< Grid points along longest bbox axis.
   type(surface_stl_object),   intent(out)           :: surface_out !< Result of MC extraction at iso=0.
   integer(I4P),               intent(out), optional :: status      !< MC_STATUS_* (or 0 on success).
   type(vector_R8P)                                  :: bmin, bmax, diag
   real(R8P)                                         :: longest, margin, h
   integer(I4P)                                      :: nx, ny, nz, i, j, k
   real(R8P), allocatable                            :: values(:, :, :)
   type(vector_R8P)                                  :: p
   type(facet_object), allocatable                   :: facets(:)
   integer(I4P)                                      :: mc_status

   if (present(status)) status = 0_I4P

   ! Bbox + margin so the iso-surface stays well inside the grid.
   bmin = self%get_bmin() ; bmax = self%get_bmax()
   diag = bmax - bmin
   longest = max(diag%x, diag%y, diag%z)
   margin  = 0.05_R8P * longest
   bmin    = bmin - vector_R8P(margin, margin, margin)
   bmax    = bmax + vector_R8P(margin, margin, margin)
   diag    = bmax - bmin

   ! Cell size: pick `resolution` corners along the longest axis, scale others.
   longest = max(diag%x, diag%y, diag%z)
   h  = longest / real(resolution - 1, R8P)
   nx = max(2, int(diag%x / h, I4P) + 1)
   ny = max(2, int(diag%y / h, I4P) + 1)
   nz = max(2, int(diag%z / h, I4P) + 1)

   allocate(values(nx, ny, nz))
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            p = vector_R8P(bmin%x + (i - 1) * (diag%x / real(nx - 1, R8P)), &
                           bmin%y + (j - 1) * (diag%y / real(ny - 1, R8P)), &
                           bmin%z + (k - 1) * (diag%z / real(nz - 1, R8P)))
            values(i, j, k) = self%distance(point=p, is_signed=.true., is_square_root=.true.)
         enddo
      enddo
   enddo

   call compute_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                           surface=facets, status=mc_status)
   if (present(status)) status = mc_status
   if (mc_status /= MC_STATUS_OK) return
   call surface_out%adopt_facets(facets=facets)
   endsubroutine resample_via_distance_field

   subroutine decimate(self, target_facets, status)
   !< Reduce `self`'s triangle count to ≤ `target_facets` via QEM edge
   !< collapse (issue #18 §1.3, Garland & Heckbert 1997).
   !<
   !< On success, `self` is replaced by the decimated mesh — adopt_facets
   !< runs analyze, rebuilding the AABB tree, vertex pool, connectivity,
   !< and pseudo-normals so all downstream queries (distance, winding number,
   !< boolean, ...) work without further user action.
   !<
   !< If the algorithm cannot reach `target_facets` (every remaining edge
   !< collapse would violate a safety check — normal-flip, non-manifold edge,
   !< duplicate facet), it returns DEC_STATUS_NO_PROGRESS and leaves `self`
   !< at the smallest count it could safely reach.
   !<
   !< Pre-condition: `self` must have been analyzed (i.e. `vertex_id` and
   !< `fcon_edge` populated). All public construction paths
   !< (`load_from_file`, `adopt_facets`, `analyze`) ensure this.
   class(surface_stl_object),     intent(inout)        :: self
   integer(I4P),                  intent(in)           :: target_facets
   integer(I4P),                  intent(out), optional :: status
   type(facet_object), allocatable                     :: working(:)
   integer(I4P)                                        :: dec_status

   if (present(status)) status = DEC_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = DEC_STATUS_BAD_INPUT
      return
   endif
   if (self%facets_number <= target_facets) return  ! nothing to do

   ! Work on a deep copy so the input stays valid until we adopt the output.
   allocate(working(self%facets_number))
   working(1:self%facets_number) = self%facet(1:self%facets_number)

   call compute_decimate(facet=working, target_facets=target_facets, status=dec_status)
   if (present(status)) status = dec_status
   if (dec_status == DEC_STATUS_BAD_INPUT) return
   ! Adopt regardless of NO_PROGRESS — the output is still a valid mesh,
   ! just larger than requested.
   call self%adopt_facets(facets=working)
   endsubroutine decimate

   function is_point_inside_polyhedron_ri(self, point) result(is_inside)
   !< Determinate is a point is inside or not to a polyhedron described by STL facets by means ray intersections count.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in) :: self           !< File STL.
   type(vector_R8P),          intent(in) :: point          !< Point coordinates.
   logical                               :: is_inside      !< Check result.
   logical                               :: is_inside_by_x !< Test result by x-aligned ray intersections.
   logical                               :: is_inside_by_y !< Test result by y-aligned ray intersections.
   logical                               :: is_inside_by_z !< Test result by z-aligned ray intersections.

   is_inside_by_x = is_inside_by_ray_intersect(ray_origin=point, ray_direction=      ex_R8P + EPS * ey_R8P + EPS * ez_R8P)
   is_inside_by_y = is_inside_by_ray_intersect(ray_origin=point, ray_direction=EPS * ex_R8P +       ey_R8P + EPS * ez_R8P)
   if (is_inside_by_x.and.is_inside_by_y) then
     is_inside = .true.
   else
      is_inside_by_z = is_inside_by_ray_intersect(ray_origin=point, ray_direction=EPS * ex_R8P + EPS * ey_R8P + ez_R8P)
      is_inside = ((is_inside_by_x.and.is_inside_by_y).or.&
                   (is_inside_by_x.and.is_inside_by_z).or.&
                   (is_inside_by_y.and.is_inside_by_z))
   endif
   contains
      function is_inside_by_ray_intersect(ray_origin, ray_direction) result(is_inside_by)
      !< Generic line intersect test.
      type(vector_R8P), intent(in) :: ray_origin           !< Ray origin.
      type(vector_R8P), intent(in) :: ray_direction        !< Ray direction.
      integer(I4P)                 :: intersections_number !< Ray intersections number of STL polyhedra with respect point.
      integer(I4P)                 :: f                    !< Counter.
      logical                      :: is_inside_by         !< Test result.

      intersections_number = 0

      if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
         ! exploit AABB refinement levels
         intersections_number = self%aabb%ray_intersections_number(facet=self%facet, &
                                                                   ray_origin=ray_origin, ray_direction=ray_direction)
      else
         ! brute-force search over all facets
         do f=1, self%facets_number
            if (self%facet(f)%do_ray_intersect(ray_origin=ray_origin, ray_direction=ray_direction)) &
               intersections_number = intersections_number + 1
         enddo
      endif

      if (mod(intersections_number, 2) == 0) then
        is_inside_by = .false.
      else
        is_inside_by = .true.
      endif
      endfunction is_inside_by_ray_intersect
   endfunction is_point_inside_polyhedron_ri

   pure function is_point_inside_polyhedron_sa(self, point) result(is_inside)
   !< Determinate is a point is inside or not to a polyhedron described by STL facets by means of the solid angle criteria.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in) :: self        !< File STL.
   type(vector_R8P),          intent(in) :: point       !< Point coordinates.
   logical                               :: is_inside   !< Check result.
   real(R8P)                             :: solid_angle !< Solid angle of STL polyhedra projected on point unit sphere.
   integer(I4P)                          :: f           !< Counter.

   solid_angle = 0._R8P
   do f=1, self%facets_number
      solid_angle = solid_angle + self%facet(f)%solid_angle(point=point)
   enddo
   if (solid_angle < -2._R8P * PI .or. 2._R8P * PI < solid_angle) then
     is_inside = .true.
   else
     is_inside = .false.
   endif
   endfunction is_point_inside_polyhedron_sa

   elemental subroutine initialize(self, aabb_refinement_levels)
   !< Initialize file.
   class(surface_stl_object), intent(inout)        :: self                   !< File STL.
   integer(I4P),              intent(in), optional :: aabb_refinement_levels !< AABB refinement levels.

   call self%destroy
   if (present(aabb_refinement_levels)) call self%aabb%set_refinement_levels(aabb_refinement_levels)
   endsubroutine initialize

   pure function largest_edge_len(self) result(largest)
   !< Return the largest edge length.
   class(surface_stl_object), intent(in) :: self    !< File STL.
   real(R8P)                             :: largest !< Largest edge length.
   integer(I4P)                          :: f       !< Counter.

   largest = 0._R8P
   if (self%facets_number>0) then
      do f=1, self%facets_number
         largest = max(largest, self%facet(f)%largest_edge_len())
      enddo
   endif
   endfunction largest_edge_len

   ! pure subroutine merge_solids(self, other)
   subroutine merge_solids(self, other, status)
   !< Merge facets with ones of other STL file.
   class(surface_stl_object), intent(inout)        :: self     !< File STL.
   type(surface_stl_object),  intent(in)           :: other    !< Other file STL.
   integer(I4P),              intent(out), optional :: status  !< 0=success, 1=allocation failure.
   type(facet_object), allocatable                 :: facet(:) !< Facets temporary list.
   integer(I4P)                                    :: f        !< Counter.
   integer(I4P)                                    :: istat    !< Allocation status.
   character(len=256)                              :: msg      !< Allocation error message.

   if (present(status)) status = STATUS_OK
   if (other%facets_number > 0) then
      if (self%facets_number > 0) then
         allocate(facet(1:self%facets_number + other%facets_number), stat=istat, errmsg=msg)
         if (istat /= 0) then
            if (present(status)) then ; status = STATUS_ALLOC_FAIL ; return ; endif
            error stop 'surface_stl_object%merge_solids: '//trim(msg)
         endif
         do f=1, self%facets_number
            facet(f)  =  self%facet(f)
         enddo
         do f=1, other%facets_number
            facet(self%facets_number+f) = other%facet(f)
         enddo
         call move_alloc(from=facet, to=self%facet)
         self%facets_number = self%facets_number + other%facets_number
      else
         allocate(self%facet(1:other%facets_number), stat=istat, errmsg=msg)
         if (istat /= 0) then
            if (present(status)) then ; status = STATUS_ALLOC_FAIL ; return ; endif
            error stop 'surface_stl_object%merge_solids: '//trim(msg)
         endif
         do f=1, other%facets_number
            self%facet(f) = other%facet(f)
         enddo
         self%facets_number = other%facets_number
      endif
      call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
   endif
   endsubroutine merge_solids

   subroutine resize(self, x, y, z, factor, respect_centroid, recompute_metrix, status)
   !< Resize (scale) facets by x or y or z or vectorial factors.
   !<
   !< @note The name `scale` has not been used, it been a Fortran built-in.
   !<
   !< @note If centroid must be used for center of resize it must be already computed.
   class(surface_stl_object), intent(inout)        :: self              !< File STL.
   real(R8P),                 intent(in), optional :: x                 !< Factor along x axis.
   real(R8P),                 intent(in), optional :: y                 !< Factor along y axis.
   real(R8P),                 intent(in), optional :: z                 !< Factor along z axis.
   type(vector_R8P),          intent(in), optional :: factor            !< Vectorial factor.
   logical,                   intent(in), optional :: respect_centroid  !< Sentinel to activate centroid as resize center.
   logical,                   intent(in), optional :: recompute_metrix  !< Sentinel to activate metrix recomputation.
   integer(I4P),              intent(out), optional :: status           !< 0=success, 2=ambiguous arguments.
   type(vector_R8P)                                :: factor_           !< Vectorial factor, local variable.
   logical                                         :: respect_centroid_ !< Sentinel to activate centroid as resize center, local v.

   if (present(status)) status = STATUS_OK
   if (present(factor) .and. (present(x) .or. present(y) .or. present(z))) then
      if (present(status)) then ; status = STATUS_AMBIGUOUS_ARGS ; return ; endif
      error stop 'surface_stl_object%resize: specify either factor or x/y/z, not both'
   endif
   respect_centroid_ = .false. ; if (present(respect_centroid)) respect_centroid_ = respect_centroid
   if (self%facets_number>0) then
      factor_ = 1._R8P
      if (present(factor)) then
         factor_ = factor
      else
         if (present(x)) factor_%x = x
         if (present(y)) factor_%y = y
         if (present(z)) factor_%z = z
      endif
      if (respect_centroid_) then
         call self%facet%resize(factor=factor_, center=self%centroid)
      else
         call self%facet%resize(factor=factor_, center=0 * ex_R8P)
      endif
      if (present(recompute_metrix)) then
         if (recompute_metrix) call self%compute_metrix
      endif
   endif
   endsubroutine resize

   elemental subroutine reverse_normals(self)
   !< Reverse facets normals.
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) call self%facet%reverse_normal
   endsubroutine reverse_normals

   ! pure subroutine sanitize(self, do_analysis)
   subroutine sanitize(self, do_analysis, status)
   !< Sanitize STL — top-level orchestrator for the mesh-repair pipeline.
   !<
   !< Runs every repair pass in the correct order:
   !<   1. (optional) initial `analyze` so bbox extents and connectivity exist.
   !<   2. `remove_degenerate_facets` — drop zero-area / sliver triangles BEFORE
   !<      anything that depends on their normals (they would propagate NaN).
   !<   3. `connect_nearby_vertices` — snap coincident vertices via union-find so
   !<      that integer-id connectivity matches geometric coincidence.
   !<   4. `analyze` — rebuild connectivity / vertex pool / volume / centroid.
   !<   5. `remove_duplicate_facets` — drop literal duplicates (any winding); if
   !<      anything was removed, re-`analyze` so winding fixup sees a clean state.
   !<   6. `sanitize_normals` — BFS-propagate winding consistency, then global
   !<      flip if the volume sign indicates inward orientation.
   !<   7. Warnings to stderr summarising every counter the user might care about.
   !<
   !< NaN/Inf scrubbing happens earlier, at `load_from_file`, since loading garbage
   !< coordinates should fail outright rather than be repaired.
   !<
   !< Post-sanitize, the composite predicates `is_watertight`, `is_manifold`, and
   !< `is_volume` summarise the mesh's repair state in one boolean each.
   class(surface_stl_object), intent(inout)        :: self        !< File STL.
   logical,                   intent(in), optional :: do_analysis !< Sentil for performing a first analysis.
   integer(I4P),              intent(out), optional :: status     !< 0=success (reserved for future use).

   if (present(status)) status = STATUS_OK
   if (self%facets_number>0) then
      if (present(do_analysis)) then
         if (do_analysis) call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
      endif
      ! Drop zero-area / sliver facets first: their NaN normals would otherwise
      ! corrupt every downstream step (connectivity, pseudo-normals, signed dist).
      ! analyze must have run at least once for bbox extents to be defined; the
      ! `do_analysis` branch above or the load_from_file path covers this.
      call self%remove_degenerate_facets
      if (self%facet_1_de%ids_number>0.or.&
          self%facet_2_de%ids_number>0.or.&
          self%facet_3_de%ids_number>0) call self%connect_nearby_vertices
      call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
      ! Drop duplicate triangles before winding fixup so duplicate copies do not
      ! bias the BFS or the volume sign.
      call self%remove_duplicate_facets
      if (self%duplicate_facets_removed > 0) &
         call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
      call self%sanitize_normals
      if (self%degenerate_facets_removed > 0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%degenerate_facets_removed, &
                                  ' degenerate facet(s) (zero-area / sliver) removed'
      if (self%duplicate_facets_removed > 0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%duplicate_facets_removed, &
                                  ' duplicate facet(s) removed'
      if (self%facet_1_de%ids_number>0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%facet_1_de%ids_number,' facet(s) with 1 disconnected edge remain'
      if (self%facet_2_de%ids_number>0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%facet_2_de%ids_number,' facet(s) with 2 disconnected edges remain'
      if (self%facet_3_de%ids_number>0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%facet_3_de%ids_number,' facet(s) with 3 disconnected edges remain'
      if (self%non_manifold_edges_number > 0) &
         write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%non_manifold_edges_number, &
                                  ' non-manifold edge(s) (3+ incident facets) detected'
   endif
   endsubroutine sanitize

   pure subroutine sanitize_normals(self)
   !< Sanitize facets normals — make windings consistent across the connectivity graph
   !< and oriented outward.
   !<
   !< Algorithm:
   !<   1. BFS over the connectivity graph. From each unvisited seed facet, walk its
   !<      neighbors through `make_normal_consistent`, which flips a neighbor's winding
   !<      iff that neighbor's edge direction matches the seed's edge direction (i.e.
   !<      same winding sense across the shared edge — wrong for an oriented manifold).
   !<      After the queue empties, restart from any remaining unvisited facet to cover
   !<      disconnected components.
   !<
   !<   2. Outward-orientation check. With the standard divergence-theorem convention
   !<      (see `tetrahedron_volume`), the closed-surface volume is positive iff
   !<      windings point outward. So if `volume < 0`, flip every facet's winding.
   !<
   !< Why the previous greedy walk was wrong: it did `f = ff` (the last neighbor it
   !< visited) and exited at the first dead end, never backtracking. On non-trivial
   !< meshes this leaves whole subgraphs unvisited, producing mixed orientations the
   !< volume-sign flip cannot resolve (cancellation in the sum).
   !<
   !< @note Facets connectivity and normals must be already computed.
   class(surface_stl_object), intent(inout) :: self             !< File STL.
   logical,      allocatable                :: facet_checked(:) !< Visited flag per facet.
   integer(I4P), allocatable                :: queue(:)         !< BFS work queue (facet ids).
   integer(I4P)                             :: head, tail       !< Queue head/tail indices.
   integer(I4P)                             :: seed, f, e       !< Counters.
   integer(I4P)                             :: neighbour        !< Neighbour facet id along an edge.

   if (self%facets_number > 0) then
      allocate(facet_checked(1:self%facets_number))
      allocate(queue(1:self%facets_number))
      facet_checked = .false.

      ! Outer loop covers disconnected components: pick any unvisited facet as a new
      ! seed and BFS from it. For a single connected mesh, the inner BFS visits
      ! everything and the outer loop exits after one iteration.
      do seed = 1, self%facets_number
         if (facet_checked(seed)) cycle
         head = 1
         tail = 1
         queue(tail) = seed
         facet_checked(seed) = .true.
         do while (head <= tail)
            f = queue(head)
            head = head + 1
            do e = 1, 3
               neighbour = self%facet(f)%fcon_edge(e)
               if (neighbour <= 0) cycle
               if (facet_checked(neighbour)) cycle
               call self%facet(f)%make_normal_consistent(edge=e, other=self%facet(neighbour))
               facet_checked(neighbour) = .true.
               tail = tail + 1
               queue(tail) = neighbour
            enddo
         enddo
      enddo
   endif
   call self%compute_volume
   if (self%volume < 0) call self%reverse_normals
   ! Refresh cached volume after the (possible) global flip, so get_volume() returns
   ! the post-sanitize value without callers needing to invoke compute_volume themselves.
   if (self%volume < 0) call self%compute_volume
   endsubroutine sanitize_normals

   pure function smallest_edge_len(self) result(smallest)
   !< Return the smallest edge length.
   class(surface_stl_object), intent(in) :: self     !< File STL.
   real(R8P)                             :: smallest !< Smallest edge length.
   integer(I4P)                          :: f        !< Counter.

   smallest = MaxR8P
   if (self%facets_number>0) then
      do f=1, self%facets_number
         smallest = min(smallest, self%facet(f)%smallest_edge_len())
      enddo
   endif
   endfunction smallest_edge_len

   pure function statistics(self, prefix) result(stats)
   !< Return STL statistics.
   class(surface_stl_object), intent(in)           :: self             !< File STL.
   character(*),              intent(in), optional :: prefix           !< Lines prefix.
   character(len=:), allocatable                   :: stats            !< STL statistics.
   character(len=:), allocatable                   :: prefix_          !< Lines prefix, local variable.
   character(1), parameter                         :: NL=new_line('a') !< Line terminator.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   stats = ''
   if (self%facets_number > 0) then
      stats=stats//prefix_//'X extents: ['//trim(str(self%bmin%x))//', '//trim(str(self%bmax%x))//']'//NL
      stats=stats//prefix_//'Y extents: ['//trim(str(self%bmin%y))//', '//trim(str(self%bmax%y))//']'//NL
      stats=stats//prefix_//'Z extents: ['//trim(str(self%bmin%z))//', '//trim(str(self%bmax%z))//']'//NL
      stats=stats//prefix_//'volume: '//trim(str(self%volume))//NL
      stats=stats//prefix_//'centroid: ['//trim(str(self%centroid%x))//', '//&
                                           trim(str(self%centroid%y))//', '//&
                                           trim(str(self%centroid%z))//']'//NL
      stats=stats//prefix_//'number of facets: '//trim(str(self%facets_number))//NL
      stats=stats//prefix_//'number of facets with 1 edges disconnected: '//trim(str(self%facet_1_de%ids_number))//NL
      stats=stats//prefix_//'number of facets with 2 edges disconnected: '//trim(str(self%facet_2_de%ids_number))//NL
      stats=stats//prefix_//'number of facets with 3 edges disconnected: '//trim(str(self%facet_3_de%ids_number))//NL
      stats=stats//prefix_//'number of non-manifold edges (3+ incident facets): '// &
                              trim(str(self%non_manifold_edges_number))//NL
      stats=stats//prefix_//'degenerate facets removed (last pass): '// &
                              trim(str(self%degenerate_facets_removed))//NL
      stats=stats//prefix_//'duplicate facets removed (last pass): '// &
                              trim(str(self%duplicate_facets_removed))//NL
      stats=stats//prefix_//'number of AABB refinement levels: '//trim(str(self%aabb%get_refinement_levels()))!//NL
   endif
   endfunction statistics

   subroutine translate(self, x, y, z, delta, recompute_metrix, status)
   !< Translate facets x or y or z or vectorial delta increments.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in), optional :: x                !< Increment along x axis.
   real(R8P),                 intent(in), optional :: y                !< Increment along y axis.
   real(R8P),                 intent(in), optional :: z                !< Increment along z axis.
   type(vector_R8P),          intent(in), optional :: delta            !< Vectorial increment.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P),              intent(out), optional :: status          !< 0=success, 2=ambiguous arguments.
   type(vector_R8P)                                :: delta_           !< Vectorial increment, local variable.

   if (present(status)) status = STATUS_OK
   if (present(delta) .and. (present(x) .or. present(y) .or. present(z))) then
      if (present(status)) then ; status = STATUS_AMBIGUOUS_ARGS ; return ; endif
      error stop 'surface_stl_object%translate: specify either delta or x/y/z, not both'
   endif
   if (self%facets_number>0) then
      delta_ = 0._R8P
      if (present(delta)) then
         delta_ = delta
      else
         if (present(x)) delta_%x = x
         if (present(y)) delta_%y = y
         if (present(z)) delta_%z = z
      endif
      call self%facet%translate(delta=delta_, recompute_metrix=recompute_metrix)
      if (present(recompute_metrix)) then
         if (recompute_metrix) call self%compute_metrix
      endif
   endif
   endsubroutine translate

   ! finaliser
   subroutine surface_stl_finalize(self)
   !< Release facet(:) and reset state.
   !<
   !< Nested components (aabb%node, facet_*_de%id, per-facet allocatables) have their
   !< own finalisers / `=` semantics; releasing facet(:) here cascades correctly.
   type(surface_stl_object), intent(inout) :: self !< Surface.

   if (allocated(self%facet)) deallocate(self%facet)
   self%facets_number = 0
   endsubroutine surface_stl_finalize

   ! private methods
   pure subroutine compute_facets_disconnected(self)
   !< Compute facets with disconnected edges.
   class(surface_stl_object), intent(inout) :: self  !< File STL.
   logical                                  :: de(3) !< Flag to check edges disconnection.
   integer(I4P)                             :: f     !< Counter.

   call self%facet_1_de%destroy
   call self%facet_2_de%destroy
   call self%facet_3_de%destroy
   if (self%facets_number>0) then
      do f=1, self%facets_number
         de = (self%facet(f)%fcon_edge == 0_I4P)
         select case(count(de))
         case(1_I4P)
            call self%facet_1_de%put(id=self%facet(f)%id)
         case(2_I4P)
            call self%facet_2_de%put(id=self%facet(f)%id)
         case(3_I4P)
            call self%facet_3_de%put(id=self%facet(f)%id)
         endselect
      enddo
   endif
   endsubroutine compute_facets_disconnected

   elemental subroutine mirror_by_normal(self, normal, recompute_metrix)
   !< Mirror facets given normal of mirroring plane.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),          intent(in)           :: normal           !< Normal of mirroring plane.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                       :: matrix(3,3)      !< Mirroring matrix.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = mirror_matrix_R8P(normal=normal)
      do f=1, self%facets_number
         call self%facet(f)%mirror(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine mirror_by_normal

   pure subroutine mirror_by_matrix(self, matrix, recompute_metrix)
   !< Mirror facet given matrix (of mirroring).
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in)           :: matrix(3,3)      !< Mirroring matrix.
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         call self%facet(f)%mirror(matrix=matrix, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine mirror_by_matrix

   subroutine rotate_by_axis_angle(self, axis, angle, center, recompute_metrix)
   !< Rotate facets given axis and angle.
   !<
   !< Angle must be in radians. When `center` is supplied, the rotation pivots
   !< about that point (`v -> R*(v - center) + center`); otherwise it pivots
   !< about the world origin. Issue #6 -- pass `center=self%get_centroid()` for
   !< body-frame rotation around the surface centroid.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   type(vector_R8P),          intent(in)           :: axis             !< Axis of rotation.
   real(R8P),                 intent(in)           :: angle            !< Angle of rotation.
   type(vector_R8P),          intent(in), optional :: center           !< Rotation centre (default: world origin).
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   real(R8P)                                       :: matrix(3,3)      !< Rotation matrix.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      matrix = rotation_matrix_R8P(axis=axis, angle=angle)
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, center=center, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_axis_angle

   pure subroutine rotate_by_matrix(self, matrix, center, recompute_metrix)
   !< Rotate facets given matrix (of rotation).
   !<
   !< When `center` is supplied, the rotation pivots about that point
   !< (`v -> matrix*(v - center) + center`); otherwise it pivots about the
   !< world origin. Issue #6 -- pass `center=self%get_centroid()` for
   !< body-frame rotation around the surface centroid.
   class(surface_stl_object), intent(inout)        :: self             !< File STL.
   real(R8P),                 intent(in)           :: matrix(3,3)      !< Rotation matrix.
   type(vector_R8P),          intent(in), optional :: center           !< Rotation centre (default: world origin).
   logical,                   intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                                    :: f                !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         call self%facet(f)%rotate(matrix=matrix, center=center, recompute_metrix=recompute_metrix)
      enddo
   endif
   endsubroutine rotate_by_matrix

   elemental subroutine set_facets_id(self)
   !< (Re)set facets ID.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: f    !< Counter.

   if (self%facets_number>0) then
      do f=1, self%facets_number
         self%facet(f)%id = f
      enddo
   endif
   endsubroutine set_facets_id

   ! non-TBP helpers

   pure function sign_algorithm_from_string(name) result(code)
   !< Translate a human-readable algorithm name into a `SIGN_*` code.
   !<
   !< Intended for CLI / config-file boundaries where users type a string and the library
   !< wants the integer dispatch code. Unknown names raise `error stop` with the valid set.
   !< Comparison is case-sensitive — keep names in canonical form (`'ray_intersections'`,
   !< `'solid_angle'`) at call sites.
   character(*), intent(in) :: name !< Algorithm name string.
   integer(I4P)             :: code !< Corresponding SIGN_* code.

   select case(trim(adjustl(name)))
   case('ray_intersections')
      code = SIGN_RAY_INTERSECTIONS
   case('solid_angle')
      code = SIGN_SOLID_ANGLE
   case('pseudo_normal')
      code = SIGN_PSEUDO_NORMAL
   case default
      error stop 'fossil_surface_stl_object%sign_algorithm_from_string: unknown name "'//trim(adjustl(name))// &
                 '" (valid: "ray_intersections", "solid_angle", "pseudo_normal")'
   endselect
   endfunction sign_algorithm_from_string

   ! file I/O (migrated from the deleted fossil_file_stl_object module)
   !
   ! These TBPs and their internal helpers replace what used to live on `file_stl_object`.
   ! Callers now do `call surface%load_from_file('foo.stl')` directly — there is no
   ! intermediate file handle to manage. Internally the routines use a local `file_unit`
   ! variable; format-mode (ASCII vs binary) is a local flag, not stored state.

   subroutine load_from_file(self, file_name, is_ascii, guess_format, clip_min, clip_max, &
                             aabb_refinement_levels, aabb_tree_kind, status)
   !< Load an STL file into the surface.
   !<
   !< Builds a local facet array, then transfers ownership via `adopt_facets` (which
   !< runs `analyze`). Auto-detects ASCII vs binary when `guess_format=.true.` (size
   !< identity-check; see the binary-header trap discussion in audit #14 S3).
   !< When `clip_min`/`clip_max` are present, only facets entirely inside the AABB
   !< are loaded.
   class(surface_stl_object),       intent(inout)        :: self                   !< Surface STL.
   character(*),                    intent(in)           :: file_name              !< STL file path.
   logical,                         intent(in), optional :: is_ascii               !< Force ASCII (default .true. if guess_format=.false.).
   logical,                         intent(in), optional :: guess_format           !< Auto-detect format from file size.
   type(vector_R8P),                intent(in), optional :: clip_min, clip_max     !< AABB clip extents (facets inside only).
   integer(I4P),                    intent(in), optional :: aabb_refinement_levels !< AABB refinement levels passed to analyze.
   integer(I4P),                    intent(in), optional :: aabb_tree_kind         !< AABB_TREE_OCTREE or AABB_TREE_SAH_BVH.
   integer(I4P),                    intent(out), optional :: status                !< 0=success, 1=alloc failure, 3=file not found.
   type(facet_object), allocatable                       :: facets(:)              !< Local buffer for ownership transfer.
   integer(I4P)                                          :: file_unit              !< File unit.
   logical                                               :: is_ascii_              !< Effective ASCII flag.
   integer(I4P)                                          :: facets_number          !< Facet count from header.
   type(facet_object)                                    :: facet_clip             !< Buffer for clipped loading.
   integer(I4P)                                          :: f, ff                  !< Counters.
   integer(I4P)                                          :: istat                  !< Allocation status.
   character(len=256)                                    :: msg                    !< Allocation error message.

   if (present(status)) status = STATUS_OK
   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_read(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_, guess_format=guess_format, status=istat)
   if (istat /= 0) then
      if (present(status)) then ; status = istat ; return ; endif
   endif
   call stl_load_facets_number(file_unit=file_unit, is_ascii=is_ascii_, facets_number=facets_number)
   call stl_load_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header)
   if (present(clip_min).and.present(clip_max)) then
      ! two-pass: count, then re-read and keep facets inside the AABB
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call facet_clip%load_from_file_ascii(file_unit=file_unit)
         else
            call facet_clip%load_from_file_binary(file_unit=file_unit)
         endif
         if (triangle_overlaps_aabb(bmin=clip_min, bmax=clip_max,  &
                                    v1=facet_clip%vertex(1),       &
                                    v2=facet_clip%vertex(2),       &
                                    v3=facet_clip%vertex(3))) ff = ff + 1
      enddo
      call stl_load_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header)
      allocate(facets(1:ff), stat=istat, errmsg=msg)
      if (istat /= 0) then
         if (present(status)) then ; status = STATUS_ALLOC_FAIL ; close(file_unit) ; return ; endif
         error stop 'surface_stl_object%load_from_file: '//trim(msg)
      endif
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call facet_clip%load_from_file_ascii(file_unit=file_unit)
         else
            call facet_clip%load_from_file_binary(file_unit=file_unit)
         endif
         if (triangle_overlaps_aabb(bmin=clip_min, bmax=clip_max,  &
                                    v1=facet_clip%vertex(1),       &
                                    v2=facet_clip%vertex(2),       &
                                    v3=facet_clip%vertex(3))) then
            ff = ff + 1
            facets(ff) = facet_clip
            facets(ff)%id = ff
         endif
      enddo
   else
      allocate(facets(1:facets_number), stat=istat, errmsg=msg)
      if (istat /= 0) then
         if (present(status)) then ; status = STATUS_ALLOC_FAIL ; close(file_unit) ; return ; endif
         error stop 'surface_stl_object%load_from_file: '//trim(msg)
      endif
      do f=1, facets_number
         if (is_ascii_) then
            call facets(f)%load_from_file_ascii(file_unit=file_unit)
         else
            call facets(f)%load_from_file_binary(file_unit=file_unit)
         endif
         facets(f)%id = f
      enddo
   endif
   close(file_unit)

   ! Defensive scan: any NaN/Inf in a vertex coordinate poisons every downstream
   ! geometric computation (normal -> NaN, distance -> NaN, AABB extents -> ±Inf).
   ! Refuse to load such a mesh rather than letting silent NaNs propagate.
   if (allocated(facets)) then
      do f = 1, size(facets)
         if (.not. (ieee_is_finite(facets(f)%vertex(1)%x) .and. &
                    ieee_is_finite(facets(f)%vertex(1)%y) .and. &
                    ieee_is_finite(facets(f)%vertex(1)%z) .and. &
                    ieee_is_finite(facets(f)%vertex(2)%x) .and. &
                    ieee_is_finite(facets(f)%vertex(2)%y) .and. &
                    ieee_is_finite(facets(f)%vertex(2)%z) .and. &
                    ieee_is_finite(facets(f)%vertex(3)%x) .and. &
                    ieee_is_finite(facets(f)%vertex(3)%y) .and. &
                    ieee_is_finite(facets(f)%vertex(3)%z))) then
            deallocate(facets)
            if (present(status)) then ; status = STATUS_INVALID_INPUT ; return ; endif
            error stop 'surface_stl_object%load_from_file: NaN/Inf vertex coordinates in input STL'
         endif
      enddo
   endif

   call self%adopt_facets(facets=facets, aabb_refinement_levels=aabb_refinement_levels, aabb_tree_kind=aabb_tree_kind)
   endsubroutine load_from_file

   subroutine save_into_file(self, file_name, is_ascii, status)
   !< Save the surface to an STL file.
   class(surface_stl_object), intent(in)           :: self      !< Surface STL.
   character(*),              intent(in)           :: file_name !< STL file path.
   logical,                   intent(in), optional :: is_ascii  !< Write as ASCII (default .true.).
   integer(I4P),              intent(out), optional :: status   !< 0=success, 4=open failure.
   integer(I4P)                                    :: file_unit !< File unit.
   logical                                         :: is_ascii_ !< Effective ASCII flag.
   integer(I4P)                                    :: f, istat  !< Counter and open status.

   if (present(status)) status = STATUS_OK
   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_write(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_, status=istat)
   if (istat /= 0) then
      if (present(status)) then ; status = istat ; return ; endif
   endif
   call stl_save_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header, facets_number=self%facets_number)
   if (is_ascii_) then
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_ascii(file_unit=file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_binary(file_unit=file_unit)
      enddo
   endif
   if (is_ascii_) write(file_unit, '(A)') 'endsolid '//trim(self%header)
   close(file_unit)
   endsubroutine save_into_file

   subroutine save_aabb_into_file(self, base_file_name, is_ascii)
   !< Save each AABB tree leaf as a separate STL file.
   !<
   !< File name pattern: `<base_file_name>aabb-l_<level>-b_<box>.stl`.
   class(surface_stl_object), intent(in)           :: self           !< Surface STL.
   character(*),              intent(in)           :: base_file_name !< Base file-name prefix.
   logical,                   intent(in), optional :: is_ascii       !< Write as ASCII (default .true.).
   type(facet_object), allocatable                 :: aabb_facet(:)  !< AABB leaf facets.
   type(surface_stl_object)                        :: leaf           !< Scratch surface for save.
   integer(I4P)                                    :: b, l           !< Counters.
   character(len=:), allocatable                   :: file_name      !< Composed file name.

   if (.not. self%aabb%get_is_initialized()) return
   do while (self%aabb%loop_node(facet=self%facet, aabb_facet=aabb_facet, b=b, l=l))
      file_name = trim(adjustl(base_file_name))//'aabb-l_'//trim(str(l, .true.))//'-b_'//trim(str(b, .true.))//'.stl'
      call leaf%adopt_facets(facets=aabb_facet)
      call leaf%set_header(self%header)
      call leaf%save_into_file(file_name=file_name, is_ascii=is_ascii)
      call leaf%destroy
   enddo
   endsubroutine save_aabb_into_file

   ! internal I/O helpers (module-private; not exposed as TBPs)

   subroutine stl_open_for_read(file_name, file_unit, is_ascii, guess_format, status)
   !< Open an STL file for reading; auto-detect format via size identity if requested.
   character(*),  intent(in)              :: file_name    !< File path.
   integer(I4P),  intent(out)             :: file_unit    !< Newunit-assigned unit.
   logical,       intent(inout)           :: is_ascii     !< In/out: format flag (rewritten by guess).
   logical,       intent(in),   optional  :: guess_format !< Auto-detect via file-size identity.
   integer(I4P),  intent(out),  optional  :: status       !< 0=success, 3=file not found.
   logical                                :: guess_format_, file_exist
   integer(I4P)                           :: file_size, facets_count, ios
   integer(I4P), parameter                :: BINARY_HEADER_BYTES = 80_I4P
   integer(I4P), parameter                :: BINARY_FACET_BYTES  = 50_I4P
   integer(I4P), parameter                :: BINARY_COUNT_BYTES  =  4_I4P

   if (present(status)) status = STATUS_OK
   file_unit = -1
   guess_format_ = .false. ; if (present(guess_format)) guess_format_ = guess_format
   inquire(file=file_name, exist=file_exist, size=file_size)
   if (.not. file_exist) then
      write(stderr, '(A)') 'error: file "'//file_name//'" does not exist, impossible to open file!'
      if (present(status)) then ; status = STATUS_FILE_NOT_FOUND ; return ; endif
      error stop 'fossil_surface_stl_object%load_from_file: file not found'
   endif
   if (guess_format_) then
      is_ascii = .true.
      if (file_size > BINARY_HEADER_BYTES + BINARY_COUNT_BYTES) then
         open(newunit=file_unit, file=file_name, access='stream', form='unformatted', action='read', iostat=ios)
         if (ios == 0) then
            read(file_unit, pos=BINARY_HEADER_BYTES + 1, iostat=ios) facets_count
            close(file_unit)
            if (ios == 0 .and. facets_count > 0 .and. &
                file_size == BINARY_HEADER_BYTES + BINARY_COUNT_BYTES + facets_count * BINARY_FACET_BYTES) then
               is_ascii = .false.
            endif
         endif
      endif
   endif
   if (is_ascii) then
      open(newunit=file_unit, file=file_name,                  form='formatted',   action='read')
   else
      open(newunit=file_unit, file=file_name, access='stream', form='unformatted', action='read')
   endif
   endsubroutine stl_open_for_read

   subroutine stl_open_for_write(file_name, file_unit, is_ascii, status)
   !< Open an STL file for writing.
   character(*), intent(in)           :: file_name !< File path.
   integer(I4P), intent(out)          :: file_unit !< Newunit-assigned unit.
   logical,      intent(in)           :: is_ascii  !< Format flag.
   integer(I4P), intent(out), optional :: status   !< 0=success, 4=open failure.
   integer(I4P)                       :: ios       !< I/O status.

   if (present(status)) status = STATUS_OK
   file_unit = -1
   if (is_ascii) then
      open(newunit=file_unit, file=file_name, form='formatted', iostat=ios)
   else
      open(newunit=file_unit, file=file_name, access='stream', form='unformatted', iostat=ios)
   endif
   if (ios /= 0) then
      write(stderr, '(A)') 'error: cannot open "'//file_name//'" for writing'
      if (present(status)) then ; status = STATUS_FILE_OPEN_FAIL ; return ; endif
      error stop 'surface_stl_object%save_into_file: cannot open file for writing'
   endif
   endsubroutine stl_open_for_write

   subroutine stl_load_facets_number(file_unit, is_ascii, facets_number)
   !< Count facets in the file (ASCII: count "facet normal" lines; binary: read uint32).
   !< File is rewound on exit.
   integer(I4P), intent(in)  :: file_unit     !< Open file unit.
   logical,      intent(in)  :: is_ascii      !< Format flag.
   integer(I4P), intent(out) :: facets_number !< Count.
   character(FRLEN)          :: facet_record  !< Line buffer for ASCII scan.

   facets_number = 0
   rewind(file_unit)
   if (is_ascii) then
      do
         read(file_unit, '(A)', end=10, err=10) facet_record
         if (index(string=facet_record, substring='facet normal') > 0) facets_number = facets_number + 1
      enddo
   else
      read(file_unit, end=10, err=10) facet_record
      read(file_unit, end=10, err=10) facets_number
   endif
   10 rewind(file_unit)
   endsubroutine stl_load_facets_number

   subroutine stl_load_header(file_unit, is_ascii, header)
   !< Read the STL header into `header`. File is rewound before reading.
   integer(I4P),     intent(in)    :: file_unit !< Open file unit.
   logical,          intent(in)    :: is_ascii  !< Format flag.
   character(FRLEN), intent(out)   :: header    !< Header text.
   integer(I4P)                    :: facets_number_dummy

   rewind(file_unit)
   if (is_ascii) then
      read(file_unit, '(A)') header
      header = trim(adjustl(header(index(header, 'solid')+6:)))
   else
      read(file_unit) header
      read(file_unit) facets_number_dummy   ! count read elsewhere; consume here to advance the stream
   endif
   endsubroutine stl_load_header

   subroutine stl_save_header(file_unit, is_ascii, header, facets_number)
   !< Write the STL header. File is rewound before writing.
   integer(I4P),     intent(in) :: file_unit     !< Open file unit.
   logical,          intent(in) :: is_ascii      !< Format flag.
   character(FRLEN), intent(in) :: header        !< Header text.
   integer(I4P),     intent(in) :: facets_number !< Facet count (binary only).

   rewind(file_unit)
   if (is_ascii) then
      write(file_unit, '(A)') 'solid '//trim(header)
   else
      write(file_unit) header
      write(file_unit) facets_number
   endif
   endsubroutine stl_save_header

   subroutine compute_pseudo_normals_via_pool(self, facet, pool)
   !< Compute edge and vertex pseudo-normals for one facet, using the pool's
   !< inverted index instead of the legacy per-facet `vertex_occurrence` list
   !< (issue #5 stage 3b).
   !<
   !< Edge pseudo-normals are unchanged from the legacy code: sum of own normal
   !< and neighbour's normal (or the facet normal alone if boundary).
   !<
   !< Vertex pseudo-normals use the Bærentzen-Aanæs angle-weighted formula:
   !<   N_v = sum over incident facets f of (incident_angle_at_v(f) * f%normal),
   !< summed across all (facet, local_v) pairs in pool%facets_at(self%vertex_id(v)),
   !< including (self, v) itself.
   type(facet_object),       intent(inout) :: self
   type(facet_object),       intent(in)    :: facet(1:)
   type(vertex_pool_object), intent(in)    :: pool
   integer(I4P)                            :: e, v, k, kn, f_n, v_n
   real(R8P)                               :: ang

   do e = 1, 3
      if (self%fcon_edge(e) > 0) then
         self%edge_pnormal(e) = self%normal + facet(self%fcon_edge(e))%normal
         call self%edge_pnormal(e)%normalize()
      else
         self%edge_pnormal(e) = self%normal
      endif
   enddo

   do v = 1, 3
      self%vertex_pnormal(v) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      kn = pool%facets_at_count(self%vertex_id(v))
      do k = 1, kn
         call pool%facets_at(self%vertex_id(v), k, f_n, v_n)
         if (f_n < 1 .or. f_n > size(facet)) cycle
         ang = facet(f_n)%vertex_angle(v_n)
         self%vertex_pnormal(v) = self%vertex_pnormal(v) + ang * facet(f_n)%normal
      enddo
      call self%vertex_pnormal(v)%normalize()
   enddo
   endsubroutine compute_pseudo_normals_via_pool

endmodule fossil_surface_stl_object
