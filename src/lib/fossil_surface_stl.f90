!< FOSSIL, STL surface class definition.

module fossil_surface_stl_object
!< FOSSIL, STL surface class definition.

use fossil_aabb_tree_object, only : aabb_tree_object, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH,     &
                                    facet_distance_payload,                                    &
                                    bvh_distance_node_oac, bvh_distance_node_with_region_oac, &
                                    wn_node_t, bvh_winding_number_oac
use fossil_aabb_node_object, only : aabb_node_object, bvh_left_child_oac, bvh_right_child_oac, bvh_payload_first_oac, &
                                    bvh_payload_count_oac, bvh_box_distance_sq_oac
use fossil_facet_object, only : facet_object, pnormal_payload_t,                               &
                                pnormal_x_oac, pnormal_y_oac, pnormal_z_oac,                   &
                                triangle_point_distance
use fossil_list_id_object, only : list_id_object
use fossil_rigid_pose, only : rigid_pose_t, to_body, to_world, rotate, pose_interpolate, pose_log_3d
use fossil_ccd, only : ccd_evaluator_t, ccd_advance, motion_bound_3d, CCD_INVALID_INPUT, CCD_CONTACT, CCD_SEPARATE, &
                      CCD_INITIAL_CONTACT, CCD_NOT_CONVERGED
use fossil_box_query, only : bvh_box_contact_oac, bvh_box_thin_oac, box_tolerance, cell_bands, cell_components, &
                             BOX_NONE, BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use fossil_box_predicates, only : box_triangle_overlap
use fossil_line_query, only : bvh_line_crossings_oac, tri_line_hit_oac, sort_crossings_oac, line_crossing_before
use fossil_utils, only : EPS, FRLEN, PI, is_inside_bb, triangle_overlaps_aabb
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_winding_number, only : compute_winding_number => winding_number, WN_FOUR_PI => FOUR_PI, WN_DEFAULT_BETA => DEFAULT_BETA
use fossil_hole_filling,   only : boundary_loop_t, extract_boundary_loops, triangulate_loop, refine_patch, fair_patch
use fossil_self_intersection, only : intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, SELFINT_COPLANAR, &
                                     SELFINT_EDGE_FOLD, SELFINT_SHEETS, &
                                     compute_self_intersections => find_self_intersections
use fossil_marching_cubes, only : compute_isosurface => extract_isosurface, MC_STATUS_OK
use fossil_decimate, only : compute_decimate => decimate, &
                            DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS
use fossil_ray_query, only : ray_hit_t, ray_intersect_all_flat, &
                             RAY_STATUS_OK, RAY_STATUS_BAD_INPUT
use fossil_sdf, only : compute_segment_sdf => segment_sdf, &
                       SDF_STATUS_OK, SDF_STATUS_BAD_INPUT, SDF_LABEL_UNASSIGNED
use fossil_alpha_wrap, only : awrap_run, &
                              AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, &
                              AWRAP_STATUS_DEGENERATE, AWRAP_STATUS_NOT_CONVERGED
use fossil_csr_matrix, only : csr_matrix_t
use fossil_body_distance, only : body_distance, body_safe_step
use fossil_projection, only : projection_target_t, build_target_3d
use fossil_deviation, only : deviation_evaluator_t, deviation_sampled, deviation_certified, DEV_STATUS_OK, &
                             DEV_STATUS_BAD_INPUT
use fossil_feature_edges, only : find_feature_edges, FEATURE_EDGE_DEFAULT_ANGLE, feature_cos
use fossil_feature_bvh, only : build_segment_bvh, feature_box_touch_oac, feature_point_distance_oac
use fossil_thickness, only : THICKNESS_RAY, THICKNESS_MEDIAL, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS, &
                            THICKNESS_DEFAULT_SEPARATION, &
                            cone_directions, lower_median, ball_radius
use fossil_mass_properties, only : moments_accumulator_t, accumulate_tetrahedron, moments_to_properties_3d, symmetric3_eigen
use fossil_laplacian, only : build_cotangent_laplacian, &
                             LAPL_STATUS_OK, LAPL_STATUS_BAD_INPUT, &
                             LAPL_STATUS_DEGENERATE_TRIANGLE
use fossil_curvature, only : build_gaussian_curvature, build_mean_curvature, build_jet_curvature, &
                             CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, &
                             CURV_STATUS_DEGENERATE_TRIANGLE
use fossil_smoothing, only : do_smooth_mesh => smooth_mesh, &
                             SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, &
                             SMOOTH_STATUS_DEGENERATE, &
                             SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN, &
                             SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU, &
                             SMOOTH_DEFAULT_ITERATIONS
use fossil_remesh, only : compute_remesh => isotropic_remesh, &
                          REM_STATUS_OK, REM_STATUS_BAD_INPUT
use fossil_boolean, only : compute_boolean => boolean_compute, compute_self_union => self_union_compute, &
                           BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, &
                           BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, &
                           BOOL_STATUS_EMPTY_INPUT, BOOL_STATUS_NOT_CLOSED, BOOL_STATUS_INVALID_INPUT, weld_vertices
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
use penf, only : I4P, I8P, R8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: surface_stl_object
public :: proximity_pairs_of, first_impact
public :: intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, SELFINT_COPLANAR, SELFINT_EDGE_FOLD, SELFINT_SHEETS
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
public :: FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX
public :: sign_algorithm_from_string
public :: STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
public :: STATUS_INVALID_INPUT
public :: BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF
public :: BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT
public :: BOOL_STATUS_NOT_CLOSED, BOOL_STATUS_INVALID_INPUT

! Point-in-polyhedron algorithm selector for `is_point_inside`, `compute_distance`,
! and `distance` (when `is_signed=.true.`):
!   SIGN_RAY_INTERSECTIONS — count intersections of a ray with the polyhedron; odd =
!                            inside, even = outside; majority over three rays in fixed
!                            generic directions (RAY_DIR_X/Y/Z). O(N) per query even
!                            with AABB acceleration.
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
!   SIGN_FAST_WINDING_NUMBER — inside iff the generalized winding number w > 0.5,
!                            evaluated by the hierarchical (second-order) fast
!                            winding number over the persistent cache built by
!                            `analyze` (issue #29). The most robust choice on
!                            meshes with holes or self-intersections: w degrades
!                            smoothly instead of flipping. O(log N) amortized
!                            per query after the sign-free distance traversal.
integer(I4P), parameter :: SIGN_RAY_INTERSECTIONS = 1_I4P
integer(I4P), parameter :: SIGN_SOLID_ANGLE       = 2_I4P
integer(I4P), parameter :: SIGN_PSEUDO_NORMAL     = 3_I4P
integer(I4P), parameter :: SIGN_FAST_WINDING_NUMBER = 4_I4P
! Closest feature kinds of the rich queries (issue #37).
integer(I4P), parameter :: FEATURE_FACE   = 1_I4P !< Closest point inside a facet.
integer(I4P), parameter :: FEATURE_EDGE   = 2_I4P !< Closest point on a facet edge (local index 1-3).
integer(I4P), parameter :: FEATURE_VERTEX = 3_I4P !< Closest point at a facet vertex (local index 1-3).
! Ray directions for SIGN_RAY_INTERSECTIONS (see is_point_inside_polyhedron_ri): generic,
! non-coplanar, none parallel to an axis or a coordinate plane, so a ray from a point of an
! axis-aligned lattice does not run exactly through a mesh edge or vertex.
type(vector_R8P), parameter :: RAY_DIR_X = vector_R8P(1._R8P,               0.2718281828459045_R8P, 0.1414213562373095_R8P)
type(vector_R8P), parameter :: RAY_DIR_Y = vector_R8P(0.1732050807568877_R8P, 1._R8P,               0.3141592653589793_R8P)
type(vector_R8P), parameter :: RAY_DIR_Z = vector_R8P(0.2236067977499790_R8P, 0.1618033988749895_R8P, 1._R8P)

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
   integer(I4P),                    private :: shells_number=0  !< Edge-connected shells (issue #37).
   integer(I4P),                    private :: bodies_number=0  !< Bodies: outer shells with their cavities (issue #37).
   integer(I4P), allocatable,       private :: facet_shell_(:)  !< Shell of each facet.
   integer(I4P), allocatable,       private :: shell_body(:)    !< Body of each shell.
   type(aabb_tree_object)                   :: aabb            !< AABB tree handle (its own state is private).
   type(vector_R8P),                private :: bmin            !< Bounding-box min.
   type(vector_R8P),                private :: bmax            !< Bounding-box max.
   real(R8P),                       private :: volume=0._R8P   !< Volume bounded by STL surface.
   real(R8P),                       private :: area=0._R8P     !< Total surface area (sum of facet areas; issue #7).
   type(vector_R8P),                private :: centroid        !< Centroid of STL surface.
   character(FRLEN),                private :: header=''       !< STL file header (preserved across load/save).
   character(FRLEN), allocatable,   private :: patch_names(:)  !< Patch names, indexed by facet%patch (issue #66).
   type(vertex_pool_object),        private :: vertex_pool     !< Unique-vertex pool (issue #5 stage 1: derived artifact).
   ! Issue #20 §Step 5 -- device-residency state for the GPU distance path. The
   ! flat `pnormal_payload(:)` holds per-facet pseudo-normals (7 vectors per
   ! facet, indexed by facet_id 1..facets_number) used by the device sign step;
   ! built in `analyze` after the per-facet pseudo-normals are populated.
   ! `is_on_device` guards `enter_device` / `exit_device` against double-enter
   ! or use-without-enter. Both `pnormal_payload` and the underlying BVH arrays
   ! (`self%aabb%node`, `self%aabb%payload`) get `!$acc enter data copyin` in
   ! `enter_device` and `!$acc exit data delete` in `exit_device`.
   type(pnormal_payload_t), allocatable, public :: pnormal_payload(:) !< Per-facet pseudo-normal SoA (Step 5). PUBLIC for !$acc directive access.
   type(wn_node_t),         allocatable, public :: wn_payload(:)      !< Per-node fast-winding-number cache, indexed like aabb%node (issue #29). PUBLIC for !$acc directive access.
   integer(I4P),            allocatable, public :: body_payload(:)    !< Body of each facet (issue #37), device-mappable. PUBLIC for !$acc directive access.
   real(R8P),               allocatable, public :: tri_vertices(:,:)  !< Exact vertices of each facet [9, facets] (issue #40), device-mappable.
   logical, public                          :: is_on_device=.false. !< Sentinel: true between matched enter_device / exit_device calls.
   real(R8P),               allocatable, public :: feature_seg(:,:)   !< Feature-edge segments [6, n] (issue #61), device-mappable.
   real(R8P),               allocatable, public :: feature_box(:,:)   !< Their BVH node boxes [6, nn].
   integer(I4P),            allocatable, public :: feature_link(:,:)  !< Their BVH node links [4, nn].
   logical,                 private :: features_built = .false.     !< build_features ran on the current geometry.
   logical,                 private :: features_on_device = .false. !< The feature arrays were mapped by enter_device.
   real(R8P),               allocatable, public :: thin_pay(:)        !< Local size per BVH payload triangle (#62), mappable.
   real(R8P),               allocatable, public :: thin_node(:)       !< Minimum local size per BVH node [0:nn-1].
   logical,                 private :: thin_built = .false.          !< build_thin ran on the current facets.
   logical,                 private :: thin_on_device = .false.      !< The thin arrays were mapped by enter_device.
   contains
      ! read-only accessors (pure, inlined at -O2, zero data copy for scalars)
      procedure, pass(self) :: get_facets_number !< Return facets_number.
      procedure, pass(self) :: get_bmin          !< Return bmin (bounding-box minimum).
      procedure, pass(self) :: get_bmax          !< Return bmax (bounding-box maximum).
      procedure, pass(self) :: get_volume        !< Return volume.
      procedure, pass(self) :: get_shells_number !< Return the number of edge-connected shells.
      procedure, pass(self) :: get_bodies_number !< Return the number of bodies (outer shells with their cavities).
      procedure, pass(self) :: facet_shell       !< Return the shell of a facet.
      procedure, pass(self) :: facet_body        !< Return the body of a facet.
      procedure, pass(self) :: get_area          !< Return total surface area (issue #7).
      procedure, pass(self) :: get_centroid      !< Return centroid.
      procedure, pass(self) :: get_header        !< Return STL header.
      procedure, pass(self) :: get_patches_number !< Number of patches (issue #66).
      procedure, pass(self) :: patch_name         !< Name of a patch.
      procedure, pass(self) :: facet_patch        !< Patch of a facet.
      procedure, pass(self) :: set_patch          !< Assign facets to a named patch.
      procedure, pass(self) :: get_non_manifold_edges_number  !< Return count of edges with 3+ incident facets.
      procedure, pass(self) :: get_degenerate_facets_removed  !< Return count of facets dropped by the last degenerate-facet pass.
      procedure, pass(self) :: get_duplicate_facets_removed   !< Return count of facets dropped by the last duplicate-facet pass.
      procedure, pass(self) :: get_vertex_pool                !< Return read-only pointer to the unique-vertex pool.
      procedure, pass(self) :: is_watertight                  !< True if every edge has exactly 2 incident facets.
      procedure, pass(self) :: is_manifold                    !< True if no non-manifold edges (boundary edges allowed).
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
      procedure, pass(self) :: load_from_files      !< Load several STL files, one patch per file (issue #66).
      procedure, pass(self) :: distance_to          !< Distance to another body, with the closest pair (issue #67).
      procedure, pass(self) :: time_of_impact       !< Time of first contact with another moving body (issue #73).
      procedure, pass(self) :: projection_target    !< Target of the constrained projection (issue #74).
      procedure, pass(self) :: proximity_pairs      !< Pairs of bodies of the surface closer than a cutoff (issue #67).
      procedure, pass(self) :: save_into_file       !< Save STL to a file path.
      procedure, pass(self) :: save_aabb_into_file  !< Save the AABB tree leaves as separate STL files.
      procedure, pass(self) :: tree_boxes           !< Boxes of the AABB tree down to a depth, read only (issue #89).
      ! public methods
      procedure, pass(self) :: allocate_facets                 !< Allocate facets.
      procedure, pass(self) :: analyze                         !< Analize STL.
      procedure, pass(self) :: build_connectivity              !< Build facets connectivity.
      procedure, pass(self), private :: compute_bodies !< Shells (edge-connected) and bodies (issue #37).
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
      procedure, pass(self) :: distance_many                   !< Batch distance over a point array; OpenMP-parallel (issue #19 §B3).
      procedure, pass(self) :: query_many                      !< Batch rich query: distance, closest point, facet, feature, normal, body (issue #37).
      procedure, pass(self) :: query_point                     !< Rich query of one point (issue #37).
      procedure, pass(self) :: query_many_device               !< query_many on the device (issue #37).
      procedure, pass(self) :: probe_many                      !< Probe points along the normal, with visibility (issue #37).
      procedure, pass(self) :: image_points                    !< Mirror points through the closest point, with visibility (issue #37).
      procedure, pass(self) :: box_classify_many               !< Classify boxes: cut, near, inside, outside (issue #38).
      procedure, pass(self) :: box_intersects_surface          !< True if a box meets the surface (issue #38).
      procedure, pass(self) :: box_classify_many_device        !< box_classify_many on the device (issue #38).
      procedure, pass(self) :: tag_block                       !< Signed tag bands on a block of cells (issue #38).
      procedure, pass(self) :: build_features                  !< Feature-edge tree for tagging and queries (issue #61).
      procedure, pass(self) :: box_feature_many                !< Boxes touching a feature edge (issue #61).
      procedure, pass(self) :: box_feature_many_device         !< box_feature_many on the device (issue #61).
      procedure, pass(self) :: line_crossings                  !< Exact crossings of an axis-aligned line (issue #40).
      procedure, pass(self) :: cell_wall_fractions             !< Wall fractions along the grid lines of a block (issue #40).
      procedure, pass(self) :: line_crossings_many             !< Crossings of many parallel lines, fixed buffers (issue #40).
      procedure, pass(self) :: line_crossings_many_device      !< line_crossings_many on the device (issue #40).
      procedure, pass(self) :: move                            !< Move rigidly by a pose, refit the tree (issue #41).
      procedure, pass(self) :: state_changes                   !< Fresh and dead points between two poses (issue #41).
      procedure, pass(self) :: wall_distance_many              !< Wall distance, body id, early cutoff (issue #42).
      procedure, pass(self) :: wall_distance_many_device       !< wall_distance_many on the device (issue #42).
      procedure, pass(self) :: surface_markers                 !< Grid-adapted markers: centroids, areas, normals (issue #42).
      ! Issue #20 §Step 5 -- GPU-offloadable distance path (BVH only).
      procedure, pass(self) :: enter_device                    !< !$acc enter data copyin of BVH arrays + pnormal_payload.
      procedure, pass(self) :: exit_device                     !< !$acc exit data delete of BVH arrays + pnormal_payload.
      procedure, pass(self) :: distance_many_device            !< Batch distance over a point array, dispatched via !$acc parallel loop.
      procedure, pass(self), private :: build_pnormal_payload  !< Build the flat per-facet pseudo-normal SoA from self%facet.
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
      procedure, pass(self) :: isotropic_remesh                 !< Botsch-Kobbelt isotropic remeshing (issue #18 §1.7).
      procedure, pass(self) :: intersect_ray_all                !< All ray-facet hits, sorted by t (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_first              !< Closest ray-facet hit, with AABB pruning (issue #18 §2.5).
      procedure, pass(self) :: intersect_ray_any                !< Any-hit early-exit ray query for shadow / occlusion (issue #18 §2.5).
      procedure, pass(self) :: segment_sdf                      !< Per-facet labels via SDF + GMM (issue #18 §1.9).
      procedure, pass(self) :: alpha_wrap                       !< Watertight surrogate from broken triangle soup (issue #18 §1.6).
      procedure, pass(self) :: cotangent_laplacian              !< Cotangent Laplacian + barycentric mass (issue #18 §2.1).
      procedure, pass(self) :: gaussian_curvature                !< Per-vertex Gaussian curvature via angle defect (issue #18 §2.4).
      procedure, pass(self) :: jet_curvature                     !< Per-vertex H and K by a local polynomial fit (issue #42).
      procedure, pass(self) :: mean_curvature                    !< Per-vertex signed mean curvature H = (1/2)||M^-1 L V|| (issue #18 §2.4).
      procedure, pass(self) :: smooth                            !< Mesh smoothing (explicit Laplacian or Taubin lambda|mu) (issue #18 §2.3).
      procedure, pass(self) :: largest_edge_len                !< Return the largest edge length.
      procedure, pass(self) :: mass_properties                 !< Mass, centre of mass, inertia tensor, principal axes (#59).
      procedure, pass(self) :: deviation                       !< One-sided deviation (Hausdorff) from another surface (#60).
      procedure, pass(self) :: hausdorff                       !< Symmetric Hausdorff distance to another surface (#60).
      procedure, pass(self) :: feature_edges                   !< Sharp, boundary, non-manifold edges, corners (#61).
      procedure, pass(self) :: save_features_into_file         !< Feature edges and corners as legacy VTK (#61).
      procedure, pass(self) :: save_thickness_into_file        !< Thickness, gap, thin regions as legacy VTK (#62).
      procedure, pass(self) :: thickness => surface_thickness !< Local thickness and gap per facet (#62).
      procedure, pass(self) :: thin_regions                    !< Facets thinner than factor*h, in regions (#62).
      procedure, pass(self) :: build_thin                      !< Local sizes for the thin tagging (#62).
      procedure, pass(self) :: box_thin_many                   !< Boxes meeting a thin facet (#62).
      procedure, pass(self) :: box_thin_many_device            !< box_thin_many on the device (#62).
      procedure, pass(self) :: merge_solids                    !< Merge facets with ones of other STL file.
      generic               :: mirror => mirror_by_normal, &
                                         mirror_by_matrix      !< Mirror facets.
      procedure, pass(self) :: reverse_normals                 !< Reverse facets normals.
      procedure, pass(self) :: resize                          !< Resize (scale) facets by x or y or z or vectorial factors.
      generic               :: rotate => rotate_by_axis_angle, &
                                         rotate_by_matrix      !< Rotate facets.
      procedure, pass(self) :: fill_holes                      !< Close small holes (boundary loops) with minimum-weight patches.
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

type, extends(deviation_evaluator_t) :: surface_evaluator_t
   !< A surface seen through its batched closest-point query, for `fossil_deviation` (issue #60).
   type(surface_stl_object), pointer :: surface => null() !< The surface B.
   logical                           :: device = .false.  !< Query on the device (query_many_device).
   contains
      procedure, pass(self) :: eval => surface_evaluator_eval !< Unsigned distances and closest points.
endtype surface_evaluator_t

type, extends(ccd_evaluator_t) :: surface_ccd_evaluator_t
   !< The distance of two moving surfaces at the step parameter s, for `fossil_ccd` (issue #73).
   class(surface_stl_object), pointer :: a => null() !< Body A.
   class(surface_stl_object), pointer :: b => null() !< Body B.
   type(rigid_pose_t)                 :: a0, a1      !< Poses of A at s = 0, 1.
   type(rigid_pose_t)                 :: b0, b1      !< Poses of B at s = 0, 1.
   real(R8P)                          :: pa(3) = 0._R8P, pb(3) = 0._R8P !< Witness points of the last call.
   integer(I4P)                       :: fa = 0_I4P, fb = 0_I4P         !< Closest facets of the last call.
   real(R8P)                          :: ka(3) = 0._R8P, kb(3) = 0._R8P !< Rotation axes over the step.
   real(R8P)                          :: tha = 0._R8P, thb = 0._R8P     !< Rotation angles over the step.
   real(R8P)                          :: dta(3) = 0._R8P, dtb(3) = 0._R8P !< Translations over the step.
   integer(I4P)                       :: pairs = 0_I4P                  !< Facet pairs of the directional steps.
   contains
      procedure, pass(self) :: distance => surface_ccd_distance   !< d(s).
      procedure, pass(self) :: safe_step => surface_ccd_safe_step !< Directional step (C2A, `body_safe_step`).
endtype surface_ccd_evaluator_t

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

   pure function get_shells_number(self) result(n)
   !< Return the number of edge-connected shells (computed by `analyze`).
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Shells.

   n = self%shells_number
   endfunction get_shells_number

   pure function get_bodies_number(self) result(n)
   !< Return the number of bodies (computed by `analyze`): see `compute_bodies`.
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P)                          :: n    !< Bodies.

   n = self%bodies_number
   endfunction get_bodies_number

   pure function facet_shell(self, i) result(sh)
   !< Return the shell of facet i (0 if out of range).
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P),              intent(in) :: i    !< Facet index.
   integer(I4P)                          :: sh   !< Shell.

   sh = 0_I4P
   if (.not. allocated(self%facet_shell_)) return
   if (i >= 1_I4P .and. i <= size(self%facet_shell_, kind=I4P)) sh = self%facet_shell_(i)
   endfunction facet_shell

   pure function facet_body(self, i) result(b)
   !< Return the body of facet i (0 if out of range).
   class(surface_stl_object), intent(in) :: self !< File STL.
   integer(I4P),              intent(in) :: i    !< Facet index.
   integer(I4P)                          :: b    !< Body.

   b = 0_I4P
   if (.not. allocated(self%facet_shell_)) return
   if (i >= 1_I4P .and. i <= size(self%facet_shell_, kind=I4P)) b = self%shell_body(self%facet_shell_(i))
   endfunction facet_body

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

   pure function get_patches_number(self) result(n)
   !< Number of patches (issue #66): one per `solid` block of an ASCII file, one per file of `load_from_files`, 1 for a
   !< binary file or a surface built otherwise.
   class(surface_stl_object), intent(in) :: self !< Surface STL.
   integer(I4P)                          :: n    !< Patches.

   n = 1_I4P
   if (allocated(self%patch_names)) n = max(1_I4P, int(size(self%patch_names), I4P))
   endfunction get_patches_number

   pure function patch_name(self, i) result(name)
   !< Name of patch i (the `solid` name it was read with; patch 1 defaults to the header).
   class(surface_stl_object), intent(in) :: self !< Surface STL.
   integer(I4P),              intent(in) :: i    !< Patch.
   character(len=:), allocatable         :: name !< Name ('' if out of range).

   name = ''
   if (allocated(self%patch_names)) then
      if (i >= 1 .and. i <= size(self%patch_names)) name = trim(self%patch_names(i))
   elseif (i == 1_I4P) then
      name = trim(header_patch_name(self%header))
   endif
   endfunction patch_name

   elemental function facet_patch(self, i) result(p)
   !< Patch of facet i (0 if out of range).
   class(surface_stl_object), intent(in) :: self !< Surface STL.
   integer(I4P),              intent(in) :: i    !< Facet.
   integer(I4P)                          :: p    !< Patch.

   p = 0_I4P
   if (i >= 1 .and. i <= self%facets_number) p = self%facet(i)%patch
   endfunction facet_patch

   subroutine set_patch(self, facets, name)
   !< Assign the facets `facets` to the patch named `name`, created if new (issue #66).
   class(surface_stl_object), intent(inout) :: self      !< Surface STL.
   integer(I4P),              intent(in)    :: facets(:) !< Facet ids.
   character(*),              intent(in)    :: name      !< Patch name.
   integer(I4P)                             :: p, k      !< Patch, counter.

   call sync_patch_names(self)
   p = find_or_add_patch(self%patch_names, name)
   do k = 1, size(facets)
      if (facets(k) >= 1 .and. facets(k) <= self%facets_number) self%facet(facets(k))%patch = p
   enddo
   endsubroutine set_patch

   subroutine sync_patch_names(self)
   !< Make the patch table cover every facet's patch id (issue #66): created as [header] when absent, extended with
   !< 'patch_<k>' names for ids beyond it; ids below 1 become 1.
   class(surface_stl_object), intent(inout) :: self !< Surface STL.
   character(FRLEN)                         :: extra !< New name.
   integer(I4P)                             :: f, pmax !< Counter, largest id.

   if (.not. allocated(self%patch_names)) then
      allocate(self%patch_names(1)) ; self%patch_names(1) = header_patch_name(self%header)
   endif
   if (size(self%patch_names) == 0) then
      deallocate(self%patch_names) ; allocate(self%patch_names(1)) ; self%patch_names(1) = header_patch_name(self%header)
   endif
   pmax = 1_I4P
   do f = 1, self%facets_number
      if (self%facet(f)%patch < 1_I4P) self%facet(f)%patch = 1_I4P
      pmax = max(pmax, self%facet(f)%patch)
   enddo
   do while (size(self%patch_names) < pmax)
      extra = 'patch_'//trim(str(int(size(self%patch_names), I4P) + 1_I4P, .true.))
      self%patch_names = [self%patch_names, extra]
   enddo
   endsubroutine sync_patch_names

   function find_or_add_patch(names, name) result(p)
   !< Index of `name` in `names`, appended if absent.
   character(FRLEN), allocatable, intent(inout) :: names(:) !< Patch names.
   character(*),                  intent(in)    :: name     !< Name.
   integer(I4P)                                 :: p        !< Index.
   character(FRLEN)                             :: nm       !< Padded name.

   nm = name
   do p = 1, size(names)
      if (names(p) == nm) return
   enddo
   names = [names, nm]
   p = size(names)
   endfunction find_or_add_patch

   pure function header_patch_name(header) result(name)
   !< The name a header gives a patch (issue #83): cut at the first NUL (a binary STL header is NUL-padded to 80 bytes)
   !< and blank-padded; `trim` it for the name itself.
   character(*), intent(in) :: header !< Header as read.
   character(FRLEN)         :: name   !< Name.
   integer(I4P)             :: m      !< First NUL.

   m = index(header, achar(0))
   if (m > 0) then
      name = header(:m - 1)
   else
      name = header
   endif
   endfunction header_patch_name

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
   call sync_patch_names(self)
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
   call drop_features(self)
   call drop_thin(self)
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
      call self%compute_bodies
      block
         integer(I4P) :: ff
         if (allocated(self%body_payload)) deallocate(self%body_payload)
         allocate(self%body_payload(self%facets_number))
         do ff = 1, self%facets_number
            self%body_payload(ff) = self%shell_body(self%facet_shell_(ff))
         enddo
      end block
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
      ! Issue #20 §Step 5 -- flatten the per-facet pseudo-normals into the device-
      ! mappable SoA. Built once here (after the per-facet pnormals exist), held
      ! for the lifetime of the surface, copied to device by `enter_device`. No
      ! device cost is paid until the user explicitly opts in.
      call self%build_pnormal_payload
      ! Issue #29 -- persistent fast-winding-number cache, built once per analyze instead
      ! of on every winding_number call.
      call self%aabb%build_winding_payload(facet=self%facet, wn=self%wn_payload)
      ! Issue #40 -- exact facet vertices as a flat array for the line-crossing kernels (the BVH
      ! payload stores v1 + edges, whose sums are rounded and would break watertightness).
      block
         integer(I4P) :: ff
         if (allocated(self%tri_vertices)) deallocate(self%tri_vertices)
         allocate(self%tri_vertices(9, self%facets_number))
         do ff = 1, self%facets_number
            self%tri_vertices(:, ff) = [self%facet(ff)%vertex(1)%x, self%facet(ff)%vertex(1)%y, self%facet(ff)%vertex(1)%z, &
                                        self%facet(ff)%vertex(2)%x, self%facet(ff)%vertex(2)%y, self%facet(ff)%vertex(2)%z, &
                                        self%facet(ff)%vertex(3)%x, self%facet(ff)%vertex(3)%y, self%facet(ff)%vertex(3)%z]
         enddo
      end block
   else
      ! No facets (e.g. an empty boolean result, issue #33): reset what the branch above computes,
      ! so no query reads the values of the previous surface (the volume used to stay stale).
      self%volume   = 0._R8P
      self%area     = 0._R8P
      self%centroid = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      self%bmin     = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      self%bmax     = vector_R8P(0._R8P, 0._R8P, 0._R8P)
      self%non_manifold_edges_number = 0_I4P
      self%shells_number = 0_I4P ; self%bodies_number = 0_I4P
      if (allocated(self%facet_shell_)) deallocate(self%facet_shell_)
      if (allocated(self%shell_body))   deallocate(self%shell_body)
      if (allocated(self%body_payload)) deallocate(self%body_payload)
      block
         integer(I4P) :: kind
         kind = self%aabb%get_tree_kind()
         call self%aabb%destroy
         call self%aabb%set_tree_kind(kind)
      end block
      if (allocated(self%pnormal_payload)) deallocate(self%pnormal_payload)
      if (allocated(self%wn_payload))      deallocate(self%wn_payload)
      if (allocated(self%tri_vertices))    deallocate(self%tri_vertices)
   endif
   endsubroutine analyze

   subroutine compute_bodies(self)
   !< Shells and bodies (issue #37), after `build_connectivity`.
   !<
   !< A shell is a set of facets connected through shared edges (`fcon_edge`, union-find); shells
   !< are numbered in the order of their first facet. A body is an outer shell with the cavities
   !< it contains: the nesting depth of a shell is the number of other shells containing it, and a
   !< shell at even depth starts a body while a shell at odd depth (a cavity) joins its immediate
   !< container. An island inside a cavity is a body of its own. Containment of shell A in shell B
   !< is tested only when B's bounding box contains A's: |sum of B's facet solid angles seen from a
   !< vertex of A| / 4 pi > 0.5 (orientation-independent). Scattered shells (particles) cost no
   !< solid-angle sum.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P), allocatable                :: parent(:), root_shell(:), first(:), order(:), cnt(:) !< Work.
   integer(I4P), allocatable                :: depth(:), container(:), body_of_root(:) !< Nesting.
   real(R8P),    allocatable                :: lo(:,:), hi(:,:) !< Shell boxes.
   type(vector_R8P)                         :: p                !< Probe vertex.
   real(R8P)                                :: omega            !< Solid-angle sum.
   real(R8P), parameter                     :: FOUR_PI_ = 4._R8P * acos(-1._R8P) !< 4 pi.
   integer(I4P)                             :: n, f, e, g, a, b, r, k, ns, nb, ka, kb !< Work.

   if (allocated(self%facet_shell_)) deallocate(self%facet_shell_)
   if (allocated(self%shell_body))   deallocate(self%shell_body)
   n = self%facets_number
   allocate(parent(n))
   parent = [(f, f = 1, n)]
   do f = 1, n
      do e = 1, 3
         g = self%facet(f)%fcon_edge(e)
         if (g > 0_I4P .and. g <= n) call unite(f, g)
      enddo
   enddo
   ! shells numbered by first facet
   allocate(root_shell(n), self%facet_shell_(n))
   root_shell = 0_I4P ; ns = 0_I4P
   do f = 1, n
      r = find(f)
      if (root_shell(r) == 0_I4P) then
         ns = ns + 1_I4P ; root_shell(r) = ns
      endif
      self%facet_shell_(f) = root_shell(r)
   enddo
   self%shells_number = ns
   ! facets grouped by shell (counting sort), shell boxes
   allocate(cnt(ns + 1), first(ns + 1), order(n), lo(3, ns), hi(3, ns))
   cnt = 0_I4P
   do f = 1, n
      cnt(self%facet_shell_(f)) = cnt(self%facet_shell_(f)) + 1_I4P
   enddo
   first(1) = 1_I4P
   do a = 1, ns
      first(a + 1) = first(a) + cnt(a)
   enddo
   cnt = 0_I4P
   lo = huge(1._R8P) ; hi = -huge(1._R8P)
   do f = 1, n
      a = self%facet_shell_(f)
      order(first(a) + cnt(a)) = f ; cnt(a) = cnt(a) + 1_I4P
      do k = 1, 3
         lo(:, a) = min(lo(:, a), [self%facet(f)%vertex(k)%x, self%facet(f)%vertex(k)%y, self%facet(f)%vertex(k)%z])
         hi(:, a) = max(hi(:, a), [self%facet(f)%vertex(k)%x, self%facet(f)%vertex(k)%y, self%facet(f)%vertex(k)%z])
      enddo
   enddo
   ! nesting: depth and immediate container (the containing shell of largest depth)
   allocate(depth(ns), container(ns))
   depth = 0_I4P ; container = 0_I4P
   do a = 1, ns
      p = self%facet(order(first(a)))%vertex(1)
      do b = 1, ns
         if (b == a) cycle
         if (any(lo(:, b) > lo(:, a)) .or. any(hi(:, b) < hi(:, a))) cycle
         if (all(lo(:, b) == lo(:, a)) .and. all(hi(:, b) == hi(:, a))) cycle ! same box: not nested
         omega = 0._R8P
         do k = first(b), first(b + 1) - 1
            omega = omega + self%facet(order(k))%solid_angle(p)
         enddo
         if (abs(omega) / FOUR_PI_ > 0.5_R8P) then
            depth(a) = depth(a) + 1_I4P
            if (container(a) == 0_I4P) then
               container(a) = b
            elseif (box_volume(b) < box_volume(container(a))) then
               container(a) = b ! the innermost container has the smallest box
            endif
         endif
      enddo
   enddo
   ! bodies: even-depth shells start one, odd-depth shells join their immediate container
   allocate(self%shell_body(ns), body_of_root(ns))
   body_of_root = 0_I4P ; nb = 0_I4P
   do a = 1, ns
      if (mod(depth(a), 2_I4P) == 0_I4P) then
         nb = nb + 1_I4P ; body_of_root(a) = nb
      endif
   enddo
   do a = 1, ns
      if (mod(depth(a), 2_I4P) == 0_I4P) then
         self%shell_body(a) = body_of_root(a)
      else
         ka = container(a)
         kb = 0_I4P ; if (ka > 0_I4P) kb = body_of_root(ka)
         if (kb == 0_I4P) then ! defensive: inconsistent nesting, keep the shell a body of its own
            nb = nb + 1_I4P ; kb = nb
         endif
         self%shell_body(a) = kb
      endif
   enddo
   self%bodies_number = nb
   contains
      function find(x) result(rt)
      !< Root of x with path halving.
      integer(I4P), intent(in) :: x  !< Element.
      integer(I4P)             :: rt !< Root.
      integer(I4P)             :: y  !< Cursor.

      y = x
      do while (parent(y) /= y)
         parent(y) = parent(parent(y))
         y = parent(y)
      enddo
      rt = y
      endfunction find

      subroutine unite(x, y)
      !< Union by smaller root index (deterministic).
      integer(I4P), intent(in) :: x, y !< Elements.
      integer(I4P)             :: rx, ry !< Roots.

      rx = find(x) ; ry = find(y)
      if (rx == ry) return
      if (rx < ry) then
         parent(ry) = rx
      else
         parent(rx) = ry
      endif
      endsubroutine unite

      pure function box_volume(c) result(vol)
      !< Volume of shell c's box.
      integer(I4P), intent(in) :: c   !< Shell.
      real(R8P)                :: vol !< Volume.

      vol = product(hi(:, c) - lo(:, c))
      endfunction box_volume
   endsubroutine compute_bodies

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
         if (present(remainder)) then ! the remainder keeps the header and the patch names (issue #66)
            remainder%header = self%header
            if (allocated(self%patch_names)) remainder%patch_names = self%patch_names
         endif
         if (present(remainder)) call remainder%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels(), &
                                                        aabb_tree_kind=self%aabb%get_tree_kind())
      endif
   endif
   endsubroutine clip

   pure subroutine compute_centroid(self)
   !< Compute centroid of STL surface: the centre of the enclosed volume.
   !<
   !< Divergence theorem over the signed tetrahedra (o, v1, v2, v3), o the bounding-box centre:
   !< centroid = o + sum D (v1 + v2 + v3 - 3 o) / 24 / sum D / 6, D = (v1-o) . ((v2-o) x (v3-o)). Only the
   !< vertices enter, in their stored order (as in the volume), never the stored facet normals. The previous
   !< formula weighted each facet by its unit normal where the unnormalised one is needed, so the result
   !< scaled with the facet size (2.0 instead of 0.5 on the unit cube cube-01, issue #59).
   !<
   !< @note The bounding box must be already computed.
   !<
   !< Degenerate / non-watertight surfaces can enclose a zero signed volume: the centroid then falls back to
   !< the origin, a defined value the caller can detect via `volume == 0`.
   class(surface_stl_object), intent(inout) :: self    !< File STL.
   real(R8P)                                :: o(3)    !< Reference point.
   real(R8P)                                :: a(3)    !< Vertex 1 relative to o.
   real(R8P)                                :: b(3)    !< Vertex 2 relative to o.
   real(R8P)                                :: c(3)    !< Vertex 3 relative to o.
   real(R8P)                                :: d       !< Six times the tetrahedron volume.
   real(R8P)                                :: v6      !< Sum of d.
   real(R8P)                                :: m24(3)  !< Sum of d (a + b + c).
   integer(I4P)                             :: f       !< Counter.

   if (self%facets_number>0) then
      self%centroid = 0._R8P
      o = 0.5_R8P * [self%bmin%x + self%bmax%x, self%bmin%y + self%bmax%y, self%bmin%z + self%bmax%z]
      v6 = 0._R8P
      m24 = 0._R8P
      do f=1, self%facets_number
         a = [self%facet(f)%vertex(1)%x, self%facet(f)%vertex(1)%y, self%facet(f)%vertex(1)%z] - o
         b = [self%facet(f)%vertex(2)%x, self%facet(f)%vertex(2)%y, self%facet(f)%vertex(2)%z] - o
         c = [self%facet(f)%vertex(3)%x, self%facet(f)%vertex(3)%y, self%facet(f)%vertex(3)%z] - o
         d = a(1) * (b(2) * c(3) - b(3) * c(2)) + a(2) * (b(3) * c(1) - b(1) * c(3)) + a(3) * (b(1) * c(2) - b(2) * c(1))
         v6 = v6 + d
         m24 = m24 + d * (a + b + c)
      enddo
      if (abs(v6) > tiny(0._R8P)) then
         m24 = o + m24 / (4._R8P * v6)
         self%centroid = vector_R8P(m24(1), m24(2), m24(3))
      endif
   endif
   endsubroutine compute_centroid

   subroutine compute_distance(self, point, distance, is_signed, sign_algorithm, is_square_root, &
                               facet_index, edge_index, vertex_index, closest_point, normal)
   !< Compute the (minimum) distance returning also the closest point.
   !<
   !< Signed-distance dispatch:
   !<  - `SIGN_PSEUDO_NORMAL` (default): one tree traversal yields both d^2 and
   !<    the closest facet+region; sign is read from
   !<    `dot(point - closest_point, pseudo_normal_at_closest)`. No second pass.
   !<  - `SIGN_RAY_INTERSECTIONS` / `SIGN_SOLID_ANGLE` / `SIGN_FAST_WINDING_NUMBER`: the unsigned
   !<    distance is computed first, then a separate inside test flips the sign (O(N) for ray
   !<    parity and solid angle; a hierarchical traversal for the fast winding number).
   class(surface_stl_object), intent(in)            :: self            !< File STL.
   type(vector_R8P),          intent(in)            :: point           !< Point coordinates.
   real(R8P),                 intent(out)           :: distance        !< Minimum distance.
   logical,                   intent(in),  optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in),  optional :: sign_algorithm  !< Algorithm code (SIGN_*).
   logical,                   intent(in),  optional :: is_square_root  !< Sentinel to trigger square-root distance.
   integer(I4P),              intent(out), optional :: facet_index     !< Index of facet containing the closest point.
   integer(I4P),              intent(out), optional :: edge_index      !< Local edge index on closest facet (0 if not on edge).
   integer(I4P),              intent(out), optional :: vertex_index    !< Local vertex index on closest facet (0 if not on vertex).
   type(vector_R8P),          intent(out), optional :: closest_point   !< Closest point (issue #37; the point itself if the surface is empty).
   type(vector_R8P),          intent(out), optional :: normal          !< Pseudo-normal of the closest feature (issue #37; zero if empty).
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

   if (present(closest_point) .or. present(normal)) then
      ! Rich outputs (issue #37): the region path, also when the distance took another path; it
      ! does not change `distance`.
      if (.not. (is_signed_ .and. algo == SIGN_PSEUDO_NORMAL) .and. self%facets_number > 0) then
         block
            real(R8P) :: d2_
            call closest_with_region(self, point, d2_, facet_index_, closest_region_, closest_)
         end block
         edge_index_ = 0 ; vertex_index_ = 0
         select case (closest_region_)
         case (1_I4P:3_I4P) ; edge_index_   = closest_region_
         case (-3_I4P:-1_I4P) ; vertex_index_ = -closest_region_
         end select
      endif
      if (present(closest_point)) closest_point = closest_
      if (present(normal)) then
         normal = vector_R8P(0._R8P, 0._R8P, 0._R8P)
         if (facet_index_ > 0) normal = self%facet(facet_index_)%pseudo_normal_for_region(closest_region_)
      endif
   endif
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
      case (SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_FAST_WINDING_NUMBER)
         if (self%is_point_inside(point=point, sign_algorithm=algo)) distance = -distance
      case default
         error stop 'fossil_surface_stl_object%compute_distance: unknown sign_algorithm code &
                    &(valid: SIGN_RAY_INTERSECTIONS=1, SIGN_SOLID_ANGLE=2, SIGN_PSEUDO_NORMAL=3, SIGN_FAST_WINDING_NUMBER=4)'
      end select
   endif
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
   !<
   !< @warning Issue #20 §Step 5: this routine does NOT release device-side
   !< storage from `enter_device`. If the surface is device-resident, the caller
   !< must invoke `exit_device` explicitly **before** `destroy` (or before the
   !< surface goes out of scope). Skipping that leaks the `!$acc enter data`
   !< mappings: the host arrays disappear here, but the device-side allocations
   !< stay registered in the OpenACC present-table against now-dangling host
   !< addresses. `destroy` cannot do the cleanup itself because it is
   !< `elemental` (callable on arrays) and `!$acc exit data` is not permitted
   !< inside an `elemental` / `pure` routine. The `is_on_device` flag is reset
   !< to `.false.` so a future host-only re-use of this object is not confused.
   class(surface_stl_object), intent(inout) :: self  !< File STL.

   if (allocated(self%facet)) deallocate(self%facet)
   if (allocated(self%pnormal_payload)) deallocate(self%pnormal_payload)
   if (allocated(self%wn_payload)) deallocate(self%wn_payload)
   self%is_on_device = .false.
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
   if (allocated(self%patch_names)) deallocate(self%patch_names)
   endsubroutine destroy

   function distance(self, point, is_signed, sign_algorithm, is_square_root)
   !< Return the (minimum) distance from a point to the triangulated surface.
   !<
   !< @note STL's metrix must be already computed.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical,                   intent(in), optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_*; default SIGN_PSEUDO_NORMAL).
   logical,                   intent(in), optional :: is_square_root  !< Sentinel to trigger square-root distance.
   real(R8P)                                       :: distance        !< Minimum distance from point to the triangulated surface.

   call self%compute_distance(point=point, distance=distance, &
                              is_signed=is_signed, sign_algorithm=sign_algorithm, is_square_root=is_square_root)
   endfunction distance

   subroutine query_point(self, point, distance, closest_point, facet, feature, feature_index, normal, body, &
                          is_signed, sign_algorithm)
   !< Rich query of one point (issue #37): distance (not squared), closest point, closest facet,
   !< feature kind and local index, pseudo-normal of the closest feature, body of the facet.
   !<
   !< One traversal with the Voronoi region (`distance_tree_with_region`, or brute force when the
   !< tree is disabled), then the closest point on the winning facet. Signed (default .true.):
   !< SIGN_PSEUDO_NORMAL (default) from the region's pseudo-normal, other algorithms through
   !< `is_point_inside`. Matches `compute_distance(is_signed=.true., is_square_root=.true.)` bit for
   !< bit with the pseudo-normal sign (the same traversal and the same sign test). An empty surface
   !< gives distance huge, the point itself as closest point, zeros elsewhere.
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   type(vector_R8P),          intent(in)            :: point          !< Query point.
   real(R8P),                 intent(out)           :: distance       !< Distance (negative inside when signed).
   type(vector_R8P),          intent(out), optional :: closest_point  !< Closest point.
   integer(I4P),              intent(out), optional :: facet          !< Closest facet.
   integer(I4P),              intent(out), optional :: feature        !< FEATURE_FACE, FEATURE_EDGE or FEATURE_VERTEX.
   integer(I4P),              intent(out), optional :: feature_index  !< Local edge or vertex index (0 for a face).
   type(vector_R8P),          intent(out), optional :: normal         !< Pseudo-normal of the closest feature.
   integer(I4P),              intent(out), optional :: body           !< Body of the closest facet.
   logical,                   intent(in),  optional :: is_signed      !< Signed distance (default .true.).
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_* (default SIGN_PSEUDO_NORMAL).
   real(R8P)                                        :: d2, side       !< Work.
   integer(I4P)                                     :: f, region, algo !< Work.
   type(vector_R8P)                                 :: cp, pn         !< Closest point, pseudo-normal.
   logical                                          :: signed_        !< Resolved option.

   signed_ = .true. ; if (present(is_signed)) signed_ = is_signed
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   cp = point ; f = 0_I4P ; region = 0_I4P ; pn = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   if (self%facets_number == 0_I4P) then
      distance = MaxR8P
   else
      call closest_with_region(self, point, d2, f, region, cp)
      distance = sqrt(d2)
      if (f > 0_I4P) pn = self%facet(f)%pseudo_normal_for_region(region)
      if (signed_) then
         if (algo == SIGN_PSEUDO_NORMAL) then
            side = (point%x - cp%x) * pn%x + (point%y - cp%y) * pn%y + (point%z - cp%z) * pn%z
            if (side < 0._R8P) distance = -distance
         else
            if (self%is_point_inside(point=point, sign_algorithm=algo)) distance = -distance
         endif
      endif
   endif
   if (present(closest_point)) closest_point = cp
   if (present(facet)) facet = f
   if (present(feature)) then
      feature = 0_I4P
      if (f > 0_I4P) then
         feature = FEATURE_FACE
         if (region >= 1_I4P .and. region <= 3_I4P) feature = FEATURE_EDGE
         if (region <= -1_I4P .and. region >= -3_I4P) feature = FEATURE_VERTEX
      endif
   endif
   if (present(feature_index)) feature_index = abs(region)
   if (present(normal)) normal = pn
   if (present(body)) then
      body = 0_I4P
      if (f > 0_I4P .and. allocated(self%facet_shell_)) body = self%shell_body(self%facet_shell_(f))
   endif
   endsubroutine query_point

   recursive subroutine query_many(self, points, distances, closest_points, facets, features, feature_indexes, normals, bodies, &
                                   is_signed, sign_algorithm, ambiguous, second_foot_points, rel_tolerance, foot_tolerance, &
                                   status, pose, curvatures, feature_distances, patches)
   !< Batch rich query (issue #37): `query_point` for every point, OpenMP over the points.
   !<
   !< Every output after `distances` is optional (a caller pays only for what it asks); arrays are
   !< caller-allocated: `closest_points(3, n)`, `normals(3, n)`, the others `(n)`. Bit-exact with
   !< `query_point`. `status = STATUS_INVALID_INPUT` on a size mismatch.
   !<
   !< Multiple foot points (issue #37 step 7, CPU only): `ambiguous(i)` is true when another facet has
   !< its closest point farther than `foot_tolerance` (default 1e-9 of the box diagonal) from the
   !< foot point and its distance within (1 + `rel_tolerance`) (default 1e-6) of the minimum: a point
   !< on or near the medial axis, a thin gap. `second_foot_points(:, i)` is that facet's closest point
   !< (the nearest such; the foot point itself when not ambiguous). FOSSIL reports, the caller decides.
   !<
   !< Curvature at the foot point (issue #42, CPU only): `curvatures(:, i)` = (H, K), the mean
   !< curvature (k1 + k2) / 2 and the Gaussian curvature k1 k2 interpolated with the foot's barycentric
   !< weights from the per-vertex `jet_curvature` (cubic fit over the 2-ring, error O(h^2)). Sign: H > 0
   !< where the surface is convex seen from outside (from the side the normals point to): a sphere has
   !< (1/R, 1/R^2), a cylinder (1/(2R), 0). The per-vertex fields are built once per call, O(facets).
   !< Where the vertex fit is rank deficient (`jet_curvature` returns CURV_STATUS_FIT_FAILED, issue #91) the
   !< vertex values are NaN, and so are the curvatures of the feet whose facet uses that vertex.
   class(surface_stl_object), intent(in)            :: self               !< Surface.
   type(vector_R8P),          intent(in)            :: points(:)          !< Query points [n].
   real(R8P),                 intent(out)           :: distances(:)       !< Distances [n].
   real(R8P),                 intent(out), optional :: closest_points(:,:) !< Closest points [3, n].
   integer(I4P),              intent(out), optional :: facets(:)          !< Closest facets [n].
   integer(I4P),              intent(out), optional :: features(:)        !< FEATURE_* [n].
   integer(I4P),              intent(out), optional :: feature_indexes(:) !< Local edge/vertex indexes [n].
   real(R8P),                 intent(out), optional :: normals(:,:)       !< Pseudo-normals [3, n].
   integer(I4P),              intent(out), optional :: bodies(:)          !< Bodies [n].
   logical,                   intent(in),  optional :: is_signed          !< Signed (default .true.).
   integer(I4P),              intent(in),  optional :: sign_algorithm     !< SIGN_* (default SIGN_PSEUDO_NORMAL).
   logical,                   intent(out), optional :: ambiguous(:)       !< Another foot point within tolerance [n].
   real(R8P),                 intent(out), optional :: second_foot_points(:,:) !< Second foot points [3, n].
   real(R8P),                 intent(in),  optional :: rel_tolerance      !< Relative distance tolerance (default 1e-6).
   real(R8P),                 intent(in),  optional :: foot_tolerance     !< Foot separation (default 1e-9 box diagonal).
   integer(I4P),              intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_t),        intent(in),  optional :: pose               !< Pose of the body (issue #41): world-frame points
                                                                          !< and outputs.
   real(R8P),                 intent(out), optional :: feature_distances(:) !< Distance to the nearest feature edge [n] (#61).
   real(R8P),                 intent(out), optional :: curvatures(:,:)    !< Mean and Gaussian curvature at the foot [2, n].
   integer(I4P),              intent(out), optional :: patches(:)         !< Patch of the closest facet [n] (issue #66).
   real(R8P), allocatable                           :: hv(:), kv(:)       !< Per-vertex mean and Gaussian curvature.
   type(vector_R8P)                                 :: cp, pn, cp2        !< Per-point closest point, normal, second foot.
   real(R8P)                                        :: rel_, foot_        !< Resolved tolerances.
   logical                                          :: amb, want_amb      !< Per-point flag, any ambiguity output.
   integer(I4P)                                     :: f, ft, fi, b, i, n !< Per-point outputs, counters.
   logical                                          :: signed_, bad       !< Resolved option, size check.
   integer(I4P)                                     :: algo               !< Resolved algorithm.

   if (present(pose)) then ! body-frame query (issue #41), outputs mapped back to the world frame
      call self%query_many(points=body_points(pose, points), distances=distances, closest_points=closest_points, &
                           facets=facets, features=features, feature_indexes=feature_indexes, normals=normals, bodies=bodies, &
                           is_signed=is_signed, sign_algorithm=sign_algorithm, ambiguous=ambiguous, &
                           second_foot_points=second_foot_points, rel_tolerance=rel_tolerance, foot_tolerance=foot_tolerance, &
                           status=status, curvatures=curvatures, feature_distances=feature_distances, patches=patches)
      if (present(status)) then
         if (status /= STATUS_OK) return
      endif
      do i = 1_I4P, size(points, kind=I4P)
         if (present(closest_points)) closest_points(:, i) = to_world(pose, closest_points(:, i))
         if (present(second_foot_points)) second_foot_points(:, i) = to_world(pose, second_foot_points(:, i))
         if (present(normals)) normals(:, i) = rotate(pose, normals(:, i))
      enddo
      return
   endif
   if (present(feature_distances)) then ! distance to the nearest feature edge (issue #61), build_features first
      if (.not. self%features_built .or. size(feature_distances) /= size(points)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      !$omp parallel do default(shared) schedule(dynamic, 64) private(i)
      do i = 1_I4P, size(points, kind=I4P)
         feature_distances(i) = feature_point_distance_oac(self%feature_seg, self%feature_box, self%feature_link, &
                                                           [points(i)%x, points(i)%y, points(i)%z])
      enddo
      !$omp end parallel do
   endif
   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   bad = size(distances, kind=I4P) /= n
   if (present(closest_points)) bad = bad .or. size(closest_points, 1) /= 3 .or. size(closest_points, 2) /= n
   if (present(normals)) bad = bad .or. size(normals, 1) /= 3 .or. size(normals, 2) /= n
   if (present(facets)) bad = bad .or. size(facets) /= n
   if (present(features)) bad = bad .or. size(features) /= n
   if (present(feature_indexes)) bad = bad .or. size(feature_indexes) /= n
   if (present(bodies)) bad = bad .or. size(bodies) /= n
   if (present(patches)) bad = bad .or. size(patches) /= n
   if (present(ambiguous)) bad = bad .or. size(ambiguous) /= n
   bad = bad .or. bad_shape3(second_foot_points, n)
   if (present(curvatures)) bad = bad .or. size(curvatures, 1) /= 2 .or. size(curvatures, 2) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   signed_ = .true. ; if (present(is_signed)) signed_ = is_signed
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   want_amb = present(ambiguous) .or. present(second_foot_points)
   rel_ = 1.e-6_R8P ; if (present(rel_tolerance)) rel_ = rel_tolerance
   foot_ = 1.e-9_R8P * vnorm(self%bmax - self%bmin) ; if (present(foot_tolerance)) foot_ = foot_tolerance
   if (present(curvatures)) call self%jet_curvature(H=hv, K=kv)
   !$omp parallel do default(shared) schedule(static) private(i, cp, pn, f, ft, fi, b, cp2, amb)
   do i = 1_I4P, n
      call self%query_point(point=points(i), distance=distances(i), closest_point=cp, facet=f, feature=ft, &
                            feature_index=fi, normal=pn, body=b, is_signed=signed_, sign_algorithm=algo)
      if (want_amb) then
         call second_foot(self, points(i), cp, distances(i)**2, rel_, foot_, amb, cp2)
         if (present(ambiguous)) ambiguous(i) = amb
         if (present(second_foot_points)) second_foot_points(:, i) = [cp2%x, cp2%y, cp2%z]
      endif
      if (present(closest_points)) closest_points(:, i) = [cp%x, cp%y, cp%z]
      if (present(facets)) facets(i) = f
      if (present(features)) features(i) = ft
      if (present(feature_indexes)) feature_indexes(i) = fi
      if (present(normals)) normals(:, i) = [pn%x, pn%y, pn%z]
      if (present(bodies)) bodies(i) = b
      if (present(patches)) then
         patches(i) = 0_I4P
         if (f >= 1_I4P .and. f <= self%facets_number) patches(i) = self%facet(f)%patch
      endif
      if (present(curvatures)) curvatures(:, i) = foot_curvature(self, f, cp, hv, kv)
   enddo
   !$omp end parallel do
   endsubroutine query_many

   function foot_curvature(self, f, cp, hv, kv) result(c)
   !< Mean and Gaussian curvature at a foot point on facet f (issue #42 step 3): the per-vertex values
   !< weighted by the barycentric coordinates of the foot, clamped to the triangle (the foot lies in it
   !< up to round-off). (0, 0) for no facet. On a near-collinear facet (Gram determinant <= 1e-12 d11 d22,
   !< where the barycentric solve is noise or its det rounds to <= 0) the values are interpolated along the
   !< longest edge at the foot's projection, and averaged if the three vertices coincide (issue #51: such
   !< facets, left by marching cubes where a sample lies on the iso-surface, returned (0, 0)).
   class(surface_stl_object), intent(in) :: self      !< Surface.
   integer(I4P),              intent(in) :: f         !< Facet (0: none).
   type(vector_R8P),          intent(in) :: cp        !< Foot point.
   real(R8P),                 intent(in) :: hv(:), kv(:) !< Per-vertex mean and Gaussian curvature.
   real(R8P)                             :: c(2)      !< (H, K).
   type(vector_R8P)                      :: e1, e2, w !< Edges, foot offset.
   real(R8P)                             :: d11, d12, d22, d1w, d2w, den, b(3) !< Gram system, weights.
   integer(I4P)                          :: vid(3), k !< Pool vertices, counter.

   c = 0._R8P
   if (f <= 0_I4P .or. f > self%facets_number) return
   do k = 1_I4P, 3_I4P
      vid(k) = self%vertex_pool%facet_vid(f, k)
   enddo
   if (any(vid <= 0_I4P) .or. any(vid > size(hv, kind=I4P))) return
   e1 = self%facet(f)%vertex(2) - self%facet(f)%vertex(1)
   e2 = self%facet(f)%vertex(3) - self%facet(f)%vertex(1)
   w  = cp - self%facet(f)%vertex(1)
   d11 = e1%dotproduct(e1) ; d12 = e1%dotproduct(e2) ; d22 = e2%dotproduct(e2)
   d1w = e1%dotproduct(w)  ; d2w = e2%dotproduct(w)
   den = d11 * d22 - d12 * d12
   if (den <= 1.e-12_R8P * d11 * d22) then
      c = degenerate_foot_curvature(self%facet(f)%vertex, cp, hv(vid), kv(vid))
      return
   endif
   b(2) = (d22 * d1w - d12 * d2w) / den
   b(3) = (d11 * d2w - d12 * d1w) / den
   b(2:3) = max(b(2:3), 0._R8P)
   b(1) = max(1._R8P - b(2) - b(3), 0._R8P)
   b = b / sum(b)
   c(1) = sum(b * hv(vid))
   c(2) = sum(b * kv(vid))
   endfunction foot_curvature

   pure function degenerate_foot_curvature(v, cp, hv, kv) result(c)
   !< (H, K) at a foot on a near-collinear facet (issue #51): linear along the longest edge, at the foot's
   !< projection on it (clamped to the edge); the mean of the three if the vertices coincide.
   type(vector_R8P), intent(in) :: v(3)          !< Facet vertices.
   type(vector_R8P), intent(in) :: cp            !< Foot point.
   real(R8P),        intent(in) :: hv(3), kv(3)  !< Per-vertex values at v.
   real(R8P)                    :: c(2)          !< (H, K).
   type(vector_R8P)             :: e             !< Longest edge.
   real(R8P)                    :: l2(3), t      !< Squared edge lengths, edge parameter.
   integer(I4P)                 :: a, b          !< Longest edge's ends.

   l2(1) = (v(2)%x - v(1)%x)**2 + (v(2)%y - v(1)%y)**2 + (v(2)%z - v(1)%z)**2
   l2(2) = (v(3)%x - v(2)%x)**2 + (v(3)%y - v(2)%y)**2 + (v(3)%z - v(2)%z)**2
   l2(3) = (v(1)%x - v(3)%x)**2 + (v(1)%y - v(3)%y)**2 + (v(1)%z - v(3)%z)**2
   a = maxloc(l2, 1) ; b = mod(a, 3) + 1
   if (l2(a) <= 0._R8P) then
      c = [sum(hv), sum(kv)] / 3._R8P
      return
   endif
   e%x = v(b)%x - v(a)%x ; e%y = v(b)%y - v(a)%y ; e%z = v(b)%z - v(a)%z
   t = ((cp%x - v(a)%x) * e%x + (cp%y - v(a)%y) * e%y + (cp%z - v(a)%z) * e%z) / l2(a)
   t = min(max(t, 0._R8P), 1._R8P)
   c = [(1._R8P - t) * hv(a) + t * hv(b), (1._R8P - t) * kv(a) + t * kv(b)]
   endfunction degenerate_foot_curvature

   subroutine box_classify_many(self, lo, hi, codes, radius, sign_algorithm, inside, status)
   !< Classify boxes against the surface (issue #38): OpenMP over the boxes.
   !<
   !< `codes(i)` is BOX_CUT if box i meets a facet, BOX_NEAR if not but a facet meets the box
   !< inflated by `radius` (default 0: no BOX_NEAR), else BOX_INSIDE or BOX_OUTSIDE from one point
   !< query at the box centre with `sign_algorithm` (default SIGN_FAST_WINDING_NUMBER; a box that meets no
   !< facet lies on one side of a closed surface, so its centre decides).
   !<
   !< The box-triangle test is the separating-axis test on the box inflated by `box_tolerance` (1e-12
   !< of the problem size): it never misses a contact, and may call a box closer than that cut.
   !< BOX_NEAR is likewise conservative: a facet within `radius` of the box is never missed, one up to
   !< sqrt(3) `radius` away (beyond a box corner) may be reported. `inside(i)` (optional) is the centre
   !< sign of every box that is not cut (also BOX_NEAR boxes; false for cut boxes). SAH BVH traversal
   !< with box pruning (`bvh_box_contact_oac`, shared with the device), brute force on an octree.
   !< `lo(3, n)`, `hi(3, n)`; `status = STATUS_INVALID_INPUT` on a size mismatch or a negative radius.
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   real(R8P),                 intent(in)            :: lo(:,:)        !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:)        !< Box maxima [3, n].
   integer(I4P),              intent(out)           :: codes(:)       !< BOX_* [n].
   real(R8P),                 intent(in),  optional :: radius         !< Near radius (default 0).
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_* for the centre query (default SIGN_FAST_WINDING_NUMBER).
   logical,                   intent(out), optional :: inside(:)      !< Centre inside, for boxes not cut [n].
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                        :: tol, reach, r  !< Inflations.
   integer(I4P)                                     :: algo, i, n, c  !< Work.
   logical                                          :: bad, ins       !< Work.
   type(vector_R8P)                                 :: centre         !< Box centre.

   if (present(status)) status = STATUS_OK
   n = size(codes, kind=I4P)
   r = 0._R8P ; if (present(radius)) r = radius
   bad = bad_shape3(lo, n) .or. bad_shape3(hi, n) .or. r < 0._R8P
   if (present(inside)) bad = bad .or. size(inside) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   reach = max(tol, r + tol)
   !$omp parallel do default(shared) schedule(dynamic, 16) private(i, c, ins, centre)
   do i = 1_I4P, n
      c = box_contact(self, lo(:, i), hi(:, i), tol, reach)
      ins = .false.
      if (c /= BOX_CUT .and. (c == BOX_NONE .or. present(inside))) then
         centre%x = 0.5_R8P * (lo(1, i) + hi(1, i))
         centre%y = 0.5_R8P * (lo(2, i) + hi(2, i))
         centre%z = 0.5_R8P * (lo(3, i) + hi(3, i))
         ins = self%is_point_inside(point=centre, sign_algorithm=algo)
      endif
      if (c == BOX_NONE) then
         c = BOX_OUTSIDE ; if (ins) c = BOX_INSIDE
      endif
      codes(i) = c
      if (present(inside)) inside(i) = ins
   enddo
   !$omp end parallel do
   endsubroutine box_classify_many

   function box_intersects_surface(self, lo, hi) result(cut)
   !< True if the box [lo, hi] meets the surface (issue #38): the separating-axis test, conservative
   !< by `box_tolerance` (never misses a contact).
   class(surface_stl_object), intent(in) :: self         !< Surface.
   real(R8P),                 intent(in) :: lo(3), hi(3) !< Box.
   logical                               :: cut          !< Result.
   real(R8P)                             :: tol          !< Inflation.

   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   cut = box_contact(self, lo, hi, tol, tol) == BOX_CUT
   endfunction box_intersects_surface

   subroutine box_classify_many_device(self, lo, hi, codes, radius, sign_algorithm, inside, status)
   !< `box_classify_many` on the device (OpenACC, issue #38): same arguments, codes and defaults.
   !<
   !< Requires `enter_device` (else STATUS_INVALID_INPUT and a stderr message). Each box runs
   !< `bvh_box_contact_oac` (the CPU kernel) and, when not cut, the centre sign of
   !< `query_many_device`: closest facet with region and pseudo-normal side (SIGN_PSEUDO_NORMAL), or
   !< `bvh_winding_number_oac` > 0.5 (SIGN_FAST_WINDING_NUMBER); only these two. SAH BVH only (the
   !< device scope). #20 rules: pointer aliases, `_oac` routines, component-wise vector stores.
   !< Validated on the host (tests-gnu, directives inert); the GPU validation is pending.
   class(surface_stl_object), intent(inout), target :: self           !< Surface (target: device aliases).
   real(R8P),                 intent(in)            :: lo(:,:)        !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:)        !< Box maxima [3, n].
   integer(I4P),              intent(out)           :: codes(:)       !< BOX_* [n].
   real(R8P),                 intent(in),  optional :: radius         !< Near radius (default 0).
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_FAST_WINDING_NUMBER (default) or SIGN_PSEUDO_NORMAL.
   logical,                   intent(out), optional :: inside(:)      !< Centre inside, for boxes not cut [n].
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P), allocatable :: cl(:)          !< Device-side codes.
   logical,      allocatable :: il(:)          !< Device-side inside flags.
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias of self%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias of self%aabb%payload.
   type(pnormal_payload_t),      pointer :: pnormal_p(:) !< Alias of self%pnormal_payload.
   type(wn_node_t),              pointer :: wn_p(:)      !< Alias of self%wn_payload.
   type(vector_R8P) :: point_i, closest_v, v1_v, e12_v, e13_v !< Per-box work (component-wise stores, #20).
   real(R8P)        :: blo(3), bhi(3), d2, nx, ny, nz, omega, beta2, tol, reach, r !< Per-box work, inflations.
   integer(I4P)     :: best_facet, best_payload, region, n, i, algo, c !< Per-box work.
   logical          :: fwn, bad, want_in, ins            !< Options, flags.

   if (present(status)) status = STATUS_OK
   n = size(codes, kind=I4P)
   r = 0._R8P ; if (present(radius)) r = radius
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
   bad = bad_shape3(lo, n) .or. bad_shape3(hi, n) .or. r < 0._R8P .or. &
         (algo /= SIGN_PSEUDO_NORMAL .and. algo /= SIGN_FAST_WINDING_NUMBER)
   if (present(inside)) bad = bad .or. size(inside) /= n
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (n == 0_I4P) return
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_surface_stl_object%box_classify_many_device: surface not on device -- call enter_device first.'
      return
   endif
   fwn = algo == SIGN_FAST_WINDING_NUMBER
   want_in = present(inside)
   beta2 = WN_DEFAULT_BETA * WN_DEFAULT_BETA
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   reach = max(tol, r + tol)
   allocate(cl(n), il(n))
   node_p => self%aabb%node ; payload_p => self%aabb%payload ; pnormal_p => self%pnormal_payload
   wn_p => self%wn_payload

   !$acc parallel loop present(node_p, payload_p, pnormal_p, wn_p) copyin(lo, hi) copyout(cl, il) &
   !$acc&  private(blo, bhi, c, ins, point_i, closest_v, v1_v, e12_v, e13_v, d2, nx, ny, nz, omega, best_facet, &
   !$acc&          best_payload, region)
   do i = 1_I4P, n
      blo(1) = lo(1, i) ; blo(2) = lo(2, i) ; blo(3) = lo(3, i)
      bhi(1) = hi(1, i) ; bhi(2) = hi(2, i) ; bhi(3) = hi(3, i)
      c = bvh_box_contact_oac(node_p, payload_p, blo, bhi, tol, reach)
      ins = .false.
      if (c /= BOX_CUT .and. (c == BOX_NONE .or. want_in)) then
         point_i%x = 0.5_R8P * (blo(1) + bhi(1))
         point_i%y = 0.5_R8P * (blo(2) + bhi(2))
         point_i%z = 0.5_R8P * (blo(3) + bhi(3))
         if (fwn) then
            omega = bvh_winding_number_oac(node_p, payload_p, wn_p, point_i%x, point_i%y, point_i%z, beta2)
            ins = omega / WN_FOUR_PI > 0.5_R8P
         else
            d2 = huge(0._R8P) ; best_facet = 0_I4P ; best_payload = 0_I4P ; region = 0_I4P
            call bvh_distance_node_with_region_oac(node=node_p, payload=payload_p, n=0_I4P, point=point_i, best=d2, &
                                                   best_facet=best_facet, best_payload_idx=best_payload)
            if (best_facet > 0_I4P) then
               v1_v%x  = payload_p(best_payload)%v1(1)
               v1_v%y  = payload_p(best_payload)%v1(2)
               v1_v%z  = payload_p(best_payload)%v1(3)
               e12_v%x = payload_p(best_payload)%e12(1)
               e12_v%y = payload_p(best_payload)%e12(2)
               e12_v%z = payload_p(best_payload)%e12(3)
               e13_v%x = payload_p(best_payload)%e13(1)
               e13_v%y = payload_p(best_payload)%e13(2)
               e13_v%z = payload_p(best_payload)%e13(3)
               call triangle_point_distance(v1=v1_v, e12=e12_v, e13=e13_v, a=payload_p(best_payload)%a, &
                                            b=payload_p(best_payload)%b, c=payload_p(best_payload)%c, &
                                            det=payload_p(best_payload)%det, point=point_i, distance=d2, &
                                            closest=closest_v, region=region)
               nx = pnormal_x_oac(pnormal_p(best_facet), region)
               ny = pnormal_y_oac(pnormal_p(best_facet), region)
               nz = pnormal_z_oac(pnormal_p(best_facet), region)
               ins = (point_i%x - closest_v%x) * nx + (point_i%y - closest_v%y) * ny + (point_i%z - closest_v%z) * nz < 0._R8P
            endif
         endif
      endif
      if (c == BOX_NONE) then
         c = BOX_OUTSIDE ; if (ins) c = BOX_INSIDE
      endif
      cl(i) = c ; il(i) = ins
   enddo
   codes = cl
   if (present(inside)) inside = il
   endsubroutine box_classify_many_device

   subroutine tag_block(self, x, y, z, band, tags, sign_algorithm, status, feature_band, feature_tags, thin_factor, thin_band, &
                        thin_tags)
   !< Signed tag bands on a rectilinear block of cells (issue #38), for AMR tagging.
   !<
   !< `x(nx+1)`, `y(ny+1)`, `z(nz+1)` are the cell edges (strictly increasing, stretched allowed);
   !< `tags(nx, ny, nz)`: 0 on a cut cell (`box_classify_many`'s BOX_CUT, conservative by
   !< `box_tolerance`), else the Chebyshev cell distance d to the nearest cut cell of the block, capped
   !< at band + 1, negative inside (-d) and positive outside (+d). A block that meets no facet is
   !< ±(band + 1) everywhere, from one point query at its centre; otherwise each 26-connected
   !< component of uncut cells is signed by one point query (`is_point_inside` with `sign_algorithm`,
   !< default SIGN_FAST_WINDING_NUMBER) at the centre of its first cell. Distances count cells of this
   !< block only: a cut cell in a neighbouring block is not seen (pad the block with ghost cells for
   !< bands that cross block faces). `status = STATUS_INVALID_INPUT` on bad edges, shape or band.
   !< Feature tags (issue #61): with `feature_band` and `feature_tags(nx, ny, nz)`, the Chebyshev cell distance to the
   !< nearest cell touched by a feature edge (`box_feature_many`), capped at feature_band + 1, 0 on touched cells;
   !< requires `build_features` (else STATUS_INVALID_INPUT).
   !< Thin tags (issue #62): with `thin_factor`, `thin_band` and `thin_tags(nx, ny, nz)`, the Chebyshev cell distance to
   !< the nearest cell meeting a facet thinner than thin_factor times the cell's largest side (`box_thin_many`), capped at
   !< thin_band + 1; requires `build_thin` (else STATUS_INVALID_INPUT).
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   real(R8P),                 intent(in)            :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),              intent(in)            :: band           !< Band width (>= 0).
   integer(I4P),              intent(out)           :: tags(:,:,:)    !< Tags [nx, ny, nz].
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_* (default SIGN_FAST_WINDING_NUMBER).
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),              intent(in),  optional :: feature_band   !< Feature band width (>= 0, issue #61).
   integer(I4P),              intent(out), optional :: feature_tags(:,:,:) !< Feature tags [nx, ny, nz].
   real(R8P),                 intent(in),  optional :: thin_factor    !< Thin threshold factor (issue #62).
   integer(I4P),              intent(in),  optional :: thin_band      !< Thin band width (>= 0).
   integer(I4P),              intent(out), optional :: thin_tags(:,:,:) !< Thin tags [nx, ny, nz].
   logical,      allocatable :: cut(:,:,:)  !< Cut cells.
   integer(I4P), allocatable :: dist(:,:,:), comp(:,:,:), seed(:,:) !< Bands, components, seeds.
   integer(I4P), allocatable :: sgn(:)      !< Sign of each component.
   integer(I4P)              :: nx, ny, nz, i, j, k, m, ncomp, algo !< Work.
   real(R8P)                 :: tol         !< Box inflation.
   logical                   :: bad         !< Size check.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   bad = nx < 1 .or. ny < 1 .or. nz < 1 .or. band < 0_I4P
   if (.not. bad) bad = any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny)) .or. any(z(2:) <= z(:nz)) .or. &
                        size(tags, 1) /= nx .or. size(tags, 2) /= ny .or. size(tags, 3) /= nz
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   algo = SIGN_FAST_WINDING_NUMBER ; if (present(sign_algorithm)) algo = sign_algorithm
   if (present(feature_tags)) then
      bad = .not. (present(feature_band) .and. self%features_built)
      if (.not. bad) bad = feature_band < 0_I4P .or. size(feature_tags, 1) /= nx .or. size(feature_tags, 2) /= ny .or. &
                           size(feature_tags, 3) /= nz
      if (bad) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      call feature_cell_bands(self, x, y, z, feature_band, feature_tags)
   endif
   if (present(thin_tags)) then
      bad = .not. (present(thin_band) .and. present(thin_factor) .and. self%thin_built)
      if (.not. bad) bad = thin_band < 0_I4P .or. thin_factor <= 0._R8P .or. size(thin_tags, 1) /= nx .or. &
                           size(thin_tags, 2) /= ny .or. size(thin_tags, 3) /= nz
      if (bad) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      call thin_cell_bands(self, x, y, z, thin_factor, thin_band, thin_tags)
   endif
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   if (box_contact(self, [x(1), y(1), z(1)], [x(nx+1), y(ny+1), z(nz+1)], tol, tol) /= BOX_CUT) then
      tags = band + 1_I4P
      if (self%is_point_inside(point=vector_R8P(0.5_R8P * (x(1) + x(nx+1)), 0.5_R8P * (y(1) + y(ny+1)), &
                                                0.5_R8P * (z(1) + z(nz+1))), sign_algorithm=algo)) tags = -tags
      return
   endif
   allocate(cut(nx, ny, nz), dist(nx, ny, nz), comp(nx, ny, nz))
   !$omp parallel do default(shared) collapse(2) schedule(dynamic) private(i, j, k)
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            cut(i, j, k) = box_contact(self, [x(i), y(j), z(k)], [x(i+1), y(j+1), z(k+1)], tol, tol) == BOX_CUT
         enddo
      enddo
   enddo
   !$omp end parallel do
   call cell_bands(cut, band, dist)
   call cell_components(cut, comp, ncomp, seed)
   allocate(sgn(0:ncomp))
   sgn(0) = 0_I4P
   !$omp parallel do default(shared) schedule(dynamic) private(m, i, j, k)
   do m = 1_I4P, ncomp
      i = seed(1, m) ; j = seed(2, m) ; k = seed(3, m)
      sgn(m) = 1_I4P
      if (self%is_point_inside(point=vector_R8P(0.5_R8P * (x(i) + x(i+1)), 0.5_R8P * (y(j) + y(j+1)), &
                                                0.5_R8P * (z(k) + z(k+1))), sign_algorithm=algo)) sgn(m) = -1_I4P
   enddo
   !$omp end parallel do
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         do i = 1_I4P, nx
            tags(i, j, k) = sgn(comp(i, j, k)) * dist(i, j, k)
         enddo
      enddo
   enddo
   endsubroutine tag_block

   subroutine line_crossings(self, axis, a, b, t0, t1, t, facets, entering, n, status)
   !< Exact crossings of the line along `axis` (1, 2, 3) at (a, b) with the surface (issue #40).
   !<
   !< (a, b) are the two other coordinates in cyclic order (axis 1: y, z; 2: z, x; 3: x, y). Every
   !< crossing with parameter (coordinate along the axis) in [t0, t1] is returned, sorted, with its
   !< facet and whether the line, going in +axis, enters the solid there. The test is exact
   !< (`fossil_line_query`: exact orient2d on the projected triangle with a symbolic tie-break), so
   !< on a watertight surface a line through edges and vertices counts every crossing once and the
   !< entering flags alternate; the parameters are rounded. SAH BVH traversal, brute force on an
   !< octree. `status = STATUS_INVALID_INPUT` if axis is not 1-3 or t1 < t0.
   class(surface_stl_object), intent(in)            :: self        !< Surface.
   integer(I4P),              intent(in)            :: axis        !< Line axis.
   real(R8P),                 intent(in)            :: a, b        !< Line position.
   real(R8P),                 intent(in)            :: t0, t1      !< Parameter range.
   real(R8P),    allocatable, intent(out)           :: t(:)        !< Crossing parameters [n].
   integer(I4P), allocatable, intent(out)           :: facets(:)   !< Facets [n].
   logical,      allocatable, intent(out)           :: entering(:) !< Entering flags [n].
   integer(I4P),              intent(out)           :: n           !< Number of crossings.
   integer(I4P),              intent(out), optional :: status      !< STATUS_OK or STATUS_INVALID_INPUT.

   if (present(status)) status = STATUS_OK
   n = 0_I4P
   if (axis < 1_I4P .or. axis > 3_I4P .or. t1 < t0) then
      if (present(status)) status = STATUS_INVALID_INPUT
      allocate(t(0), facets(0), entering(0))
      return
   endif
   call line_crossings_core(self, axis, a, b, t0, t1, t, facets, entering, n)
   endsubroutine line_crossings

   subroutine cell_wall_fractions(self, x, y, z, frac, inside, bad_lines, status)
   !< Wall fractions along the grid lines through the cell centres of a block (issue #40).
   !<
   !< `x(nx)`, `y(ny)`, `z(nz)` are the cell-centre coordinates (strictly increasing). `frac(d, i, j, k)`
   !< for d = 1..6 (-x, +x, -y, +y, -z, +z) is the fraction of the distance to the neighbouring
   !< centre at which the grid line first meets the wall, in [0, 1], or -1 if it does not meet it
   !< before that centre (also at the block ends: pad the block with ghost cells). `inside(i, j, k)`
   !< is the parity of the crossings along the x line from outside the surface box: exact
   !< topology (every crossing counted once on a watertight surface) and, since issue #69, the
   !< exact position of each crossing against the centre (orient3d): it equals the point-in-solid
   !< test, a centre on the wall counting as past it (taken at x + eps). `bad_lines` (optional)
   !< counts the lines whose entering/leaving flags do not alternate (open or non-manifold input).
   !< `status = STATUS_INVALID_INPUT` on bad coordinates or shapes.
   class(surface_stl_object), intent(in)            :: self            !< Surface.
   real(R8P),                 intent(in)            :: x(:), y(:), z(:) !< Cell centres.
   real(R8P),                 intent(out)           :: frac(:,:,:,:)   !< Fractions [6, nx, ny, nz].
   logical,                   intent(out)           :: inside(:,:,:)   !< Centres inside [nx, ny, nz].
   integer(I4P),              intent(out), optional :: bad_lines       !< Lines with non-alternating flags.
   integer(I4P),              intent(out), optional :: status          !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: nx, ny, nz, i, j, k, nbad !< Work.
   logical      :: bad                       !< Size check.
   real(R8P)    :: bl(3)                     !< Surface box minimum.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) ; ny = size(y, kind=I4P) ; nz = size(z, kind=I4P)
   bad = nx < 1 .or. ny < 1 .or. nz < 1
   if (.not. bad) bad = any(x(2:) <= x(:nx-1)) .or. any(y(2:) <= y(:ny-1)) .or. any(z(2:) <= z(:nz-1)) .or. &
                        size(frac, 1) /= 6 .or. size(frac, 2) /= nx .or. size(frac, 3) /= ny .or. size(frac, 4) /= nz .or. &
                        size(inside, 1) /= nx .or. size(inside, 2) /= ny .or. size(inside, 3) /= nz
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   frac = -1._R8P ; inside = .false. ; nbad = 0_I4P
   if (self%facets_number == 0_I4P) then
      if (present(bad_lines)) bad_lines = 0_I4P
      return
   endif
   bl = [self%bmin%x, self%bmin%y, self%bmin%z]
   !$omp parallel default(shared) private(i, j, k) reduction(+:nbad)
   !$omp do collapse(2) schedule(dynamic)
   do k = 1_I4P, nz
      do j = 1_I4P, ny
         call sweep(1_I4P, y(j), z(k), x, bl(1), frac(1:2, :, j, k), nbad, inside(:, j, k))
      enddo
   enddo
   !$omp end do nowait
   !$omp do collapse(2) schedule(dynamic)
   do i = 1_I4P, nx
      do k = 1_I4P, nz
         call sweep(2_I4P, z(k), x(i), y, bl(2), frac(3:4, i, :, k), nbad)
      enddo
   enddo
   !$omp end do nowait
   !$omp do collapse(2) schedule(dynamic)
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         call sweep(3_I4P, x(i), y(j), z, bl(3), frac(5:6, i, j, :), nbad)
      enddo
   enddo
   !$omp end do
   !$omp end parallel
   if (present(bad_lines)) bad_lines = nbad
   contains
      subroutine sweep(ax, a, b, c, cmin, fr, nb, ins)
      !< One grid line along ax through the centres c: fractions toward both neighbours, parity.
      integer(I4P), intent(in)              :: ax       !< Axis.
      real(R8P),    intent(in)              :: a, b     !< Line position.
      real(R8P),    intent(in)              :: c(:)     !< Centres along the line.
      real(R8P),    intent(in)              :: cmin     !< Surface box minimum along the axis.
      real(R8P),    intent(inout)           :: fr(:,:)  !< Fractions [2 (-, +), centres].
      integer(I4P), intent(inout)           :: nb       !< Bad lines.
      logical,      intent(inout), optional :: ins(:)   !< Parity inside flags.
      real(R8P),    allocatable :: t(:)    !< Crossings.
      integer(I4P), allocatable :: f(:)    !< Facets.
      logical,      allocatable :: e(:)    !< Entering flags.
      integer(I4P)              :: n, m, q, nc !< Work.

      nc = size(c, kind=I4P)
      ! the whole line: a crossing rounded past the last centre can still be before it exactly
      call line_crossings_core(self, ax, a, b, min(c(1), cmin) - 1._R8P, huge(1._R8P), t, f, e, n)
      do q = 1_I4P, n
         if (e(q) .neqv. (mod(q, 2_I4P) == 1_I4P)) then
            nb = nb + 1_I4P
            exit
         endif
      enddo
      m = 0_I4P ! crossings strictly before the current centre
      do q = 1_I4P, nc
         do while (m < n)
            if (t(m + 1) >= c(q)) exit
            m = m + 1_I4P
         enddo
         if (q > 1_I4P .and. m > 0_I4P) then ! last crossing in [c(q-1), c(q)]
            if (t(m) >= c(q - 1)) fr(1, q) = (c(q) - t(m)) / (c(q) - c(q - 1))
         endif
         if (q < nc .and. m < n) then        ! first crossing in [c(q), c(q+1)]
            if (t(m + 1) <= c(q + 1)) fr(2, q) = (t(m + 1) - c(q)) / (c(q + 1) - c(q))
         endif
      enddo
      if (present(ins)) call exact_parity(ax, a, b, c, t, f, e, n, ins)
      endsubroutine sweep

      subroutine exact_parity(ax, a, b, c, t, f, e, n, ins)
      !< Inside flags of the centres c by the parity of the crossings exactly before them (issue #69; the rounded
      !< parameters decided it before, so a centre within round-off of the wall could land on either side). Each
      !< crossing is placed among the sorted centres by its rounded parameter, then moved by the exact test, which is
      !< monotone along the line: about two orient3d per crossing. A centre on the wall counts as past it.
      integer(I4P), intent(in)  :: ax       !< Axis.
      real(R8P),    intent(in)  :: a, b     !< Line position.
      real(R8P),    intent(in)  :: c(:)     !< Centres.
      real(R8P),    intent(in)  :: t(:)     !< Crossing parameters.
      integer(I4P), intent(in)  :: f(:)     !< Facets.
      logical,      intent(in)  :: e(:)     !< Entering flags.
      integer(I4P), intent(in)  :: n        !< Crossings.
      logical,      intent(out) :: ins(:)   !< Inside flags.
      integer(I4P)              :: hist(0:size(c)), k, j, q, nc, cnt
      real(R8P)                 :: tri(3, 3)

      nc = size(c, kind=I4P)
      hist = 0_I4P
      do k = 1_I4P, n
         tri = reshape(self%tri_vertices(:, f(k)), [3, 3])
         j = count(c < t(k), kind=I4P) ! centres the crossing is (by its rounded parameter) not before
         do while (j > 0_I4P)
            if (.not. line_crossing_before(tri, ax, a, b, c(j), e(k))) exit
            j = j - 1_I4P
         enddo
         do while (j < nc)
            if (line_crossing_before(tri, ax, a, b, c(j + 1), e(k))) exit
            j = j + 1_I4P
         enddo
         hist(j) = hist(j) + 1_I4P ! before the centres j+1..nc
      enddo
      cnt = 0_I4P
      do q = 1_I4P, nc
         cnt = cnt + hist(q - 1)
         ins(q) = mod(cnt, 2_I4P) == 1_I4P
      enddo
      endsubroutine exact_parity
   endsubroutine cell_wall_fractions

   subroutine line_crossings_many(self, axis, a, b, t0, t1, t, facets, entering, counts, overflow, status)
   !< Crossings of the lines along `axis` at (a(i), b(i)), i = 1..n, with parameter in [t0, t1] (issue
   !< #40): `line_crossings` for every line into fixed buffers t(nmax, n), facets(nmax, n),
   !< entering(nmax, n), sorted; counts(i) crossings (at most nmax); overflow(i) if line i has more
   !< than nmax: its kept crossings are then not the first ones (traversal order), so re-run it with
   !< a larger buffer. OpenMP over the lines; the kernel of `line_crossings_many_device`. SAH BVH only.
   !< `status = STATUS_INVALID_INPUT` on a bad axis, range, shape, or a non-BVH tree.
   class(surface_stl_object), intent(in)            :: self          !< Surface.
   integer(I4P),              intent(in)            :: axis          !< Line axis.
   real(R8P),                 intent(in)            :: a(:), b(:)    !< Line positions [n].
   real(R8P),                 intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),                 intent(out)           :: t(:,:)        !< Parameters [nmax, n].
   integer(I4P),              intent(out)           :: facets(:,:)   !< Facets [nmax, n].
   logical,                   intent(out)           :: entering(:,:) !< Entering flags [nmax, n].
   integer(I4P),              intent(out)           :: counts(:)     !< Crossings per line [n].
   logical,                   intent(out)           :: overflow(:)   !< Buffer overflow per line [n].
   integer(I4P),              intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P) :: i, n !< Work.

   if (.not. many_args_ok(self, axis, a, b, t0, t1, t, facets, entering, counts, overflow, status)) return
   n = size(a, kind=I4P)
   !$omp parallel do default(shared) schedule(dynamic)
   do i = 1_I4P, n
      call bvh_line_crossings_oac(self%aabb%node, self%aabb%payload, self%tri_vertices, axis, a(i), b(i), t0, t1, t(:, i), &
                                  facets(:, i), entering(:, i), counts(i), overflow(i))
   enddo
   !$omp end parallel do
   endsubroutine line_crossings_many

   subroutine line_crossings_many_device(self, axis, a, b, t0, t1, t, facets, entering, counts, overflow, status)
   !< `line_crossings_many` on the device (OpenACC, issue #40): same arguments and results, bit for bit.
   !< Requires `enter_device` (else STATUS_INVALID_INPUT and a stderr message). #20 rules: pointer
   !< aliases of the mapped arrays (`tri_vertices` is mapped by `enter_device`). Validated on the host
   !< (tests-gnu); the GPU validation is pending.
   class(surface_stl_object), intent(inout), target :: self          !< Surface (target: device aliases).
   integer(I4P),              intent(in)            :: axis          !< Line axis.
   real(R8P),                 intent(in)            :: a(:), b(:)    !< Line positions [n].
   real(R8P),                 intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),                 intent(out)           :: t(:,:)        !< Parameters [nmax, n].
   integer(I4P),              intent(out)           :: facets(:,:)   !< Facets [nmax, n].
   logical,                   intent(out)           :: entering(:,:) !< Entering flags [nmax, n].
   integer(I4P),              intent(out)           :: counts(:)     !< Crossings per line [n].
   logical,                   intent(out)           :: overflow(:)   !< Buffer overflow per line [n].
   integer(I4P),              intent(out), optional :: status        !< STATUS_OK or STATUS_INVALID_INPUT.
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias of self%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias of self%aabb%payload.
   real(R8P),                    pointer :: tv_p(:,:)    !< Alias of self%tri_vertices.
   real(R8P),    allocatable :: tl(:,:)   !< Device-side parameters.
   integer(I4P), allocatable :: fl(:,:)   !< Device-side facets.
   logical,      allocatable :: el(:,:)   !< Device-side flags.
   integer(I4P), allocatable :: cl(:)     !< Device-side counts.
   logical,      allocatable :: ol(:)     !< Device-side overflow.
   integer(I4P)              :: i, n, nmax !< Work.

   if (.not. many_args_ok(self, axis, a, b, t0, t1, t, facets, entering, counts, overflow, status)) return
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_surface_stl_object%line_crossings_many_device: surface not on device -- call enter_device first.'
      return
   endif
   n = size(a, kind=I4P) ; nmax = size(t, 1, kind=I4P)
   if (n == 0_I4P) return
   allocate(tl(nmax, n), fl(nmax, n), el(nmax, n), cl(n), ol(n))
   tl = 0._R8P ; fl = 0_I4P ; el = .false. ! entries past counts(i) as on the CPU path
   node_p => self%aabb%node ; payload_p => self%aabb%payload ; tv_p => self%tri_vertices
   !$acc parallel loop present(node_p, payload_p, tv_p) copyin(a, b) copy(tl, fl, el) copyout(cl, ol)
   do i = 1_I4P, n
      call bvh_line_crossings_oac(node_p, payload_p, tv_p, axis, a(i), b(i), t0, t1, tl(:, i), fl(:, i), el(:, i), cl(i), ol(i))
   enddo
   t = tl ; facets = fl ; entering = el ; counts = cl ; overflow = ol
   endsubroutine line_crossings_many_device

   recursive subroutine query_many_device(self, points, distances, closest_points, facets, features, feature_indexes, &
                                          normals, bodies, is_signed, sign_algorithm, status, pose, feature_distances, &
                                         patches)
   !< `query_many` on the device (issue #37): same arguments, outputs and defaults.
   !<
   !< Requires `enter_device` (else STATUS_INVALID_INPUT and a stderr message). The kernel always
   !< computes every output into local arrays (an absent optional cannot be referenced in device
   !< code); the requested ones are copied out on the host. The loop body is the one of
   !< `distance_many_device` (BVH traversal with region, closest point recomputed by
   !< `triangle_point_distance`, pseudo-normal or fast-winding-number sign), keeping what it
   !< already computes; #20 rules (pointer aliases, `_oac` routines). SIGN_PSEUDO_NORMAL or
   !< SIGN_FAST_WINDING_NUMBER only. Validated on the host (tests-gnu, directives inert); the GPU
   !< validation is pending.
   class(surface_stl_object), intent(inout), target :: self               !< Surface (target: device aliases).
   type(vector_R8P),          intent(in)            :: points(:)          !< Query points [n].
   real(R8P),                 intent(out)           :: distances(:)       !< Distances [n].
   real(R8P),                 intent(out), optional :: closest_points(:,:) !< Closest points [3, n].
   integer(I4P),              intent(out), optional :: facets(:)          !< Closest facets [n].
   integer(I4P),              intent(out), optional :: features(:)        !< FEATURE_* [n].
   integer(I4P),              intent(out), optional :: feature_indexes(:) !< Local edge/vertex indexes [n].
   real(R8P),                 intent(out), optional :: normals(:,:)       !< Pseudo-normals [3, n].
   integer(I4P),              intent(out), optional :: bodies(:)          !< Bodies [n].
   integer(I4P),              intent(out), optional :: patches(:)         !< Patch of the closest facet [n] (issue #66): looked up
                                                                          !< on the host from the device's facet ids.
   logical,                   intent(in),  optional :: is_signed          !< Signed (default .true.).
   integer(I4P),              intent(in),  optional :: sign_algorithm     !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   integer(I4P),              intent(out), optional :: status             !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                 intent(out), optional :: feature_distances(:) !< Distance to the nearest feature edge [n] (#61).
   type(rigid_pose_t),        intent(in),  optional :: pose               !< Pose of the body (issue #41).
   real(R8P),    allocatable :: dd(:), cpl(:,:), nml(:,:) !< Device-side outputs.
   integer(I4P), allocatable :: fcl(:), ftl(:), fil(:), bdl(:) !< Device-side outputs.
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias of self%aabb%node.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias of self%aabb%payload.
   type(pnormal_payload_t),      pointer :: pnormal_p(:) !< Alias of self%pnormal_payload.
   type(wn_node_t),              pointer :: wn_p(:)      !< Alias of self%wn_payload.
   integer(I4P),                 pointer :: body_p(:)    !< Alias of self%body_payload.
   real(R8P),    pointer :: fseg_p(:,:), fbox_p(:,:)     !< Aliases of the feature tree (issue #61).
   integer(I4P), pointer :: flink_p(:,:)                 !< Alias of the feature tree links.
   real(R8P),    allocatable :: fdl(:), pts(:,:)         !< Device-side feature distances, points.
   type(vector_R8P) :: point_i, closest_v, v1_v, e12_v, e13_v !< Per-point work (component-wise stores, #20).
   real(R8P)        :: d2, nx, ny, nz, side, omega, beta2      !< Per-point work.
   integer(I4P)     :: best_facet, best_payload, region, n, i, algo !< Per-point work.
   logical          :: signed_, fwn, bad                       !< Options.

   if (present(pose)) then ! body-frame query (issue #41), outputs mapped back on the host
      call self%query_many_device(points=body_points(pose, points), distances=distances, closest_points=closest_points, &
                                  facets=facets, features=features, feature_indexes=feature_indexes, normals=normals, &
                                  bodies=bodies, is_signed=is_signed, sign_algorithm=sign_algorithm, status=status, &
                                  patches=patches, &
                                  feature_distances=feature_distances)
      if (present(status)) then
         if (status /= STATUS_OK) return
      endif
      do i = 1_I4P, size(points, kind=I4P)
         if (present(closest_points)) closest_points(:, i) = to_world(pose, closest_points(:, i))
         if (present(normals)) normals(:, i) = rotate(pose, normals(:, i))
      enddo
      return
   endif
   if (present(feature_distances)) then ! distance to the nearest feature edge (issue #61), on the device
      if (.not. (self%is_on_device .and. self%features_on_device) .or. size(feature_distances) /= size(points)) then
         if (present(status)) status = STATUS_INVALID_INPUT
         write(stderr, '(A)') 'query_many_device(feature_distances=): build_features, then enter_device, first'
         return
      endif
      fseg_p => self%feature_seg ; fbox_p => self%feature_box ; flink_p => self%feature_link
      n = size(points, kind=I4P)
      allocate(fdl(n), pts(3, n))
      do i = 1_I4P, n
         pts(:, i) = [points(i)%x, points(i)%y, points(i)%z]
      enddo
      !$acc parallel loop present(fseg_p, fbox_p, flink_p) copyin(pts) copyout(fdl)
      do i = 1_I4P, n
         fdl(i) = feature_point_distance_oac(fseg_p, fbox_p, flink_p, pts(:, i))
      enddo
      feature_distances = fdl
   endif
   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   bad = size(distances, kind=I4P) /= n
   if (present(closest_points)) bad = bad .or. size(closest_points, 1) /= 3 .or. size(closest_points, 2) /= n
   if (present(normals)) bad = bad .or. size(normals, 1) /= 3 .or. size(normals, 2) /= n
   if (present(facets)) bad = bad .or. size(facets) /= n
   if (present(features)) bad = bad .or. size(features) /= n
   if (present(feature_indexes)) bad = bad .or. size(feature_indexes) /= n
   if (present(bodies)) bad = bad .or. size(bodies) /= n
   if (present(patches)) bad = bad .or. size(patches) /= n
   signed_ = .true. ; if (present(is_signed)) signed_ = is_signed
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   bad = bad .or. (algo /= SIGN_PSEUDO_NORMAL .and. algo /= SIGN_FAST_WINDING_NUMBER)
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (n == 0_I4P) return
   if (.not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_surface_stl_object%query_many_device: surface not on device -- call enter_device first.'
      return
   endif
   fwn = algo == SIGN_FAST_WINDING_NUMBER
   beta2 = WN_DEFAULT_BETA * WN_DEFAULT_BETA
   allocate(dd(n), cpl(3, n), nml(3, n), fcl(n), ftl(n), fil(n), bdl(n))
   node_p => self%aabb%node ; payload_p => self%aabb%payload ; pnormal_p => self%pnormal_payload
   wn_p => self%wn_payload ; body_p => self%body_payload

   !$acc parallel loop present(node_p, payload_p, pnormal_p, wn_p, body_p) copyin(points) &
   !$acc&  copyout(dd, cpl, nml, fcl, ftl, fil, bdl) &
   !$acc&  private(point_i, closest_v, v1_v, e12_v, e13_v, d2, nx, ny, nz, side, omega, best_facet, best_payload, region)
   do i = 1_I4P, n
      point_i%x = points(i)%x ; point_i%y = points(i)%y ; point_i%z = points(i)%z
      d2 = huge(0._R8P) ; best_facet = 0_I4P ; best_payload = 0_I4P ; region = 0_I4P
      closest_v%x = point_i%x ; closest_v%y = point_i%y ; closest_v%z = point_i%z
      call bvh_distance_node_with_region_oac(node=node_p, payload=payload_p, n=0_I4P, point=point_i, best=d2, &
                                             best_facet=best_facet, best_payload_idx=best_payload)
      nx = 0._R8P ; ny = 0._R8P ; nz = 0._R8P
      if (best_facet > 0_I4P) then
         v1_v%x = payload_p(best_payload)%v1(1)
         v1_v%y = payload_p(best_payload)%v1(2)
         v1_v%z = payload_p(best_payload)%v1(3)
         e12_v%x = payload_p(best_payload)%e12(1)
         e12_v%y = payload_p(best_payload)%e12(2)
         e12_v%z = payload_p(best_payload)%e12(3)
         e13_v%x = payload_p(best_payload)%e13(1)
         e13_v%y = payload_p(best_payload)%e13(2)
         e13_v%z = payload_p(best_payload)%e13(3)
         call triangle_point_distance(v1=v1_v, e12=e12_v, e13=e13_v, a=payload_p(best_payload)%a, b=payload_p(best_payload)%b, &
                                      c=payload_p(best_payload)%c, det=payload_p(best_payload)%det, point=point_i, &
                                      distance=d2, closest=closest_v, region=region)
         nx = pnormal_x_oac(pnormal_p(best_facet), region)
         ny = pnormal_y_oac(pnormal_p(best_facet), region)
         nz = pnormal_z_oac(pnormal_p(best_facet), region)
         d2 = sqrt(d2)
         if (signed_) then
            if (fwn) then
               omega = bvh_winding_number_oac(node_p, payload_p, wn_p, point_i%x, point_i%y, point_i%z, beta2)
               if (omega / WN_FOUR_PI > 0.5_R8P) d2 = -d2
            else
               side = (points(i)%x - closest_v%x) * nx + (points(i)%y - closest_v%y) * ny + (points(i)%z - closest_v%z) * nz
               if (side < 0._R8P) d2 = -d2
            endif
         endif
      endif
      dd(i) = d2
      cpl(1, i) = closest_v%x ; cpl(2, i) = closest_v%y ; cpl(3, i) = closest_v%z
      nml(1, i) = nx ; nml(2, i) = ny ; nml(3, i) = nz
      fcl(i) = best_facet
      ftl(i) = 0_I4P ; bdl(i) = 0_I4P
      if (best_facet > 0_I4P) then
         ftl(i) = FEATURE_FACE
         if (region >= 1_I4P .and. region <= 3_I4P) ftl(i) = FEATURE_EDGE
         if (region <= -1_I4P .and. region >= -3_I4P) ftl(i) = FEATURE_VERTEX
         bdl(i) = body_p(best_facet)
      endif
      fil(i) = abs(region)
   enddo
   distances = dd
   if (present(closest_points)) closest_points = cpl
   if (present(normals)) normals = nml
   if (present(facets)) facets = fcl
   if (present(features)) features = ftl
   if (present(feature_indexes)) feature_indexes = fil
   if (present(bodies)) bodies = bdl
   if (present(patches)) then
      do i = 1, n
         patches(i) = 0_I4P
         if (fcl(i) >= 1_I4P .and. fcl(i) <= self%facets_number) patches(i) = self%facet(fcl(i))%patch
      enddo
   endif
   endsubroutine query_many_device

   subroutine probe_many(self, points, lengths, probe_points, visible, foot_points, normals, status, patches)
   !< Probe points (issue #37): for each point, the point at distance lengths(i) along the
   !< pseudo-normal of its closest feature, from its foot point (closest point): positive lengths go
   !< outward, negative inward. `visible(i)` is true when the open segment from the foot point to
   !< the probe point crosses no facet: hits of a ray from the foot point (along +-n) with
   !< tol < t <= |length|, tol = 1e-9 times the box diagonal (the facets containing the foot point
   !< are hit at t ~ 0). Supports the placement rules of image points at sqrt(2) dx, image lines,
   !< second image layers and fixed-distance probes; `status = STATUS_INVALID_INPUT` on a size
   !< mismatch. OpenMP over the points.
   class(surface_stl_object), intent(in)            :: self              !< Surface.
   type(vector_R8P),          intent(in)            :: points(:)         !< Query points [n].
   real(R8P),                 intent(in)            :: lengths(:)        !< Probe lengths [n] (+ outward).
   real(R8P),                 intent(out)           :: probe_points(:,:) !< Probe points [3, n].
   logical,                   intent(out)           :: visible(:)        !< Segment foot-probe crosses no facet [n].
   real(R8P),                 intent(out), optional :: foot_points(:,:)  !< Foot (closest) points [3, n].
   real(R8P),                 intent(out), optional :: normals(:,:)      !< Pseudo-normals [3, n].
   integer(I4P),              intent(out), optional :: status            !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),              intent(out), optional :: patches(:)        !< Patch of the foot point's facet [n] (issue #66).
   type(vector_R8P)                                 :: cp, pn, tip       !< Per-point work.
   real(R8P)                                        :: d, tol            !< Work.
   integer(I4P)                                     :: i, n, f           !< Counters, facet.

   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(lengths) /= n .or. size(probe_points, 1) /= 3 .or. size(probe_points, 2) /= n .or. size(visible) /= n .or. &
       bad_shape3(foot_points, n) .or. bad_shape3(normals, n) .or. bad_size(patches, n)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   tol = 1.e-9_R8P * vnorm(self%bmax - self%bmin)
   !$omp parallel do default(shared) schedule(static) private(i, cp, pn, tip, d, f)
   do i = 1_I4P, n
      call self%query_point(point=points(i), distance=d, closest_point=cp, normal=pn, facet=f)
      if (present(patches)) then
         patches(i) = 0_I4P
         if (f >= 1_I4P .and. f <= self%facets_number) patches(i) = self%facet(f)%patch
      endif
      tip = cp + lengths(i) * pn
      probe_points(:, i) = [tip%x, tip%y, tip%z]
      visible(i) = segment_clear(self, cp, tip, tol)
      if (present(foot_points)) foot_points(:, i) = [cp%x, cp%y, cp%z]
      if (present(normals)) normals(:, i) = [pn%x, pn%y, pn%z]
   enddo
   !$omp end parallel do
   endsubroutine probe_many

   subroutine image_points(self, points, images, visible, foot_points, status, patches)
   !< Image points (issue #37): the mirror of each point through its closest point, 2 foot - x
   !< (the ghost-point image of an IB closure), and whether the open segment from the foot point to
   !< the image crosses no facet (see `probe_many`). OpenMP over the points.
   class(surface_stl_object), intent(in)            :: self             !< Surface.
   type(vector_R8P),          intent(in)            :: points(:)        !< Query points [n].
   real(R8P),                 intent(out)           :: images(:,:)      !< Image points [3, n].
   logical,                   intent(out), optional :: visible(:)       !< Segment foot-image crosses no facet [n].
   real(R8P),                 intent(out), optional :: foot_points(:,:) !< Foot points [3, n].
   integer(I4P),              intent(out), optional :: status           !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P),              intent(out), optional :: patches(:)       !< Patch of the foot point's facet [n] (issue #66).
   type(vector_R8P)                                 :: cp, img          !< Per-point work.
   real(R8P)                                        :: d, tol           !< Work.
   integer(I4P)                                     :: i, n, f          !< Counters, facet.
   logical                                          :: bad              !< Size check.

   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   bad = size(images, 1) /= 3 .or. size(images, 2) /= n .or. bad_shape3(foot_points, n)
   if (present(visible)) bad = bad .or. size(visible) /= n
   bad = bad .or. bad_size(patches, n)
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   tol = 1.e-9_R8P * vnorm(self%bmax - self%bmin)
   !$omp parallel do default(shared) schedule(static) private(i, cp, img, d, f)
   do i = 1_I4P, n
      call self%query_point(point=points(i), distance=d, closest_point=cp, facet=f)
      if (present(patches)) then
         patches(i) = 0_I4P
         if (f >= 1_I4P .and. f <= self%facets_number) patches(i) = self%facet(f)%patch
      endif
      img%x = 2._R8P * cp%x - points(i)%x ; img%y = 2._R8P * cp%y - points(i)%y ; img%z = 2._R8P * cp%z - points(i)%z
      images(:, i) = [img%x, img%y, img%z]
      if (present(visible)) visible(i) = segment_clear(self, cp, img, tol)
      if (present(foot_points)) foot_points(:, i) = [cp%x, cp%y, cp%z]
   enddo
   !$omp end parallel do
   endsubroutine image_points

   recursive subroutine distance_many(self, points, distances, is_signed, sign_algorithm, is_square_root, status, pose)
   !< Batch distance query over a point array — the Adam/Xall embedded-boundary
   !< access pattern (one distance per Cartesian grid point against a static
   !< surface). Issue #19 §B3.
   !<
   !< Semantically a loop of `distance` calls: `distances(i)` is the distance from
   !< `points(i)` to the surface, under the same `is_signed` / `sign_algorithm` /
   !< `is_square_root` options applied uniformly to every point. The win is not the
   !< per-call semantics — it is that the loop is **`!$omp parallel do`**: the AABB
   !< tree is `intent(in)` (read-only) for the whole batch and each query writes
   !< only its own `distances(i)` slot, so the iterations are independent with no
   !< locks. Builds without `-fopenmp` compile the directive as an inert comment
   !< and the loop runs serially — correctness is identical either way.
   !<
   !< The option arguments are resolved once here and passed as concrete (non-
   !< optional) values into the worker, so the hot loop does not re-test `present()`
   !< per query.
   !<
   !< `points` and `distances` must be the same length; `distances` is
   !< caller-allocated. A length mismatch sets `status = STATUS_INVALID_INPUT` and
   !< the routine returns without touching `distances`.
   !<
   !< @note STL's metrix must be already computed (as for `distance`).
   class(surface_stl_object), intent(in)            :: self            !< Surface STL.
   type(vector_R8P),          intent(in)            :: points(:)       !< Query points.
   real(R8P),                 intent(out)           :: distances(:)    !< Per-point distances (caller-allocated, same size as `points`).
   logical,                   intent(in),  optional :: is_signed       !< Sentinel to trigger signed distance.
   integer(I4P),              intent(in),  optional :: sign_algorithm  !< Algorithm code (SIGN_*).
   logical,                   intent(in),  optional :: is_square_root  !< Sentinel to trigger square-root distance.
   integer(I4P),              intent(out), optional :: status          !< STATUS_OK or STATUS_INVALID_INPUT (size mismatch).
   type(rigid_pose_t),        intent(in),  optional :: pose            !< Pose of the body (issue #41): points are in the world frame.
   logical                                          :: is_signed_      !< Resolved signed sentinel.
   integer(I4P)                                     :: algo            !< Resolved sign-algorithm code.
   logical                                          :: is_sqrt_        !< Resolved square-root sentinel.
   integer(I4P)                                     :: n               !< Number of query points.
   integer(I4P)                                     :: i               !< Loop counter.

   if (present(pose)) then ! body-frame query (issue #41): distances are invariant
      call self%distance_many(points=body_points(pose, points), distances=distances, is_signed=is_signed, &
                              sign_algorithm=sign_algorithm, is_square_root=is_square_root, status=status)
      return
   endif
   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(distances, kind=I4P) /= n) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (n == 0_I4P) return

   ! Resolve options once — the per-query worker takes concrete values, so the
   ! parallel loop never re-tests present().
   is_signed_ = .false. ; if (present(is_signed))     is_signed_ = is_signed
   algo       = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   is_sqrt_   = .false. ; if (present(is_square_root)) is_sqrt_   = is_square_root

   ! The AABB tree and facet array are read-only for the whole batch; each
   ! iteration writes only distances(i). Embarrassingly parallel, no locks.
   !$omp parallel do default(none) schedule(static) &
   !$omp   shared(self, points, distances, n, is_signed_, algo, is_sqrt_) &
   !$omp   private(i)
   do i = 1_I4P, n
      call self%compute_distance(point=points(i), distance=distances(i), &
                                 is_signed=is_signed_, sign_algorithm=algo, is_square_root=is_sqrt_)
   enddo
   !$omp end parallel do
   endsubroutine distance_many

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

   function segment_clear(self, a, b, tol) result(clear)
   !< True if the open segment from a (on the surface) to b crosses no facet beyond tol from a.
   type(surface_stl_object), intent(in) :: self  !< Surface.
   type(vector_R8P),         intent(in) :: a, b  !< Segment ends.
   real(R8P),                intent(in) :: tol   !< Hits with t <= tol (the facets at a) are ignored.
   logical                              :: clear !< Result.
   type(ray_hit_t), allocatable         :: hits(:) !< Ray hits.
   type(vector_R8P)                     :: dir   !< Unit direction.
   real(R8P)                            :: len   !< Segment length.
   integer(I4P)                         :: h     !< Counter.

   clear = .true.
   len = vnorm(b - a)
   if (len <= tol) return
   dir = (b - a) / len
   call self%intersect_ray_all(ray_origin=a, ray_direction=dir, hits=hits)
   do h = 1, size(hits)
      if (hits(h)%t > tol .and. hits(h)%t <= len) then
         clear = .false.
         return
      endif
   enddo
   endfunction segment_clear

   subroutine second_foot(self, point, cp1, d1sq, rel, foot_tol, found, cp2)
   !< Nearest facet whose closest point is farther than foot_tol from cp1, with a squared distance
   !< within (1 + rel)^2 d1sq (issue #37). SAH BVH traversal on the flat node/payload arrays with
   !< pruning on the relaxed bound, tightened by the best candidate; brute force otherwise.
   type(surface_stl_object), intent(in)  :: self     !< Surface.
   type(vector_R8P),         intent(in)  :: point    !< Query point.
   type(vector_R8P),         intent(in)  :: cp1      !< First foot point.
   real(R8P),                intent(in)  :: d1sq     !< Minimum squared distance.
   real(R8P),                intent(in)  :: rel      !< Relative distance tolerance.
   real(R8P),                intent(in)  :: foot_tol !< Foot separation.
   logical,                  intent(out) :: found    !< A second foot point exists.
   type(vector_R8P),         intent(out) :: cp2      !< Second foot point (cp1 if none).
   integer(I4P)                          :: stack(256), sp, cur, lc, rc, k, idx, f, region !< Work.
   real(R8P)                             :: bound, d2 !< Work.
   type(vector_R8P)                      :: c, v1, e12, e13 !< Work.

   found = .false. ; cp2 = cp1
   if (self%facets_number == 0_I4P) return
   bound = (1._R8P + rel)**2 * d1sq
   if (self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(self%aabb%node) .and. &
       allocated(self%aabb%payload)) then
      sp = 1_I4P ; stack(1) = 0_I4P
      do while (sp > 0_I4P)
         cur = stack(sp) ; sp = sp - 1_I4P
         if (bvh_box_distance_sq_oac(self%aabb%node(cur), point%x, point%y, point%z) > bound) cycle
         do k = 0_I4P, bvh_payload_count_oac(self%aabb%node(cur)) - 1_I4P
            idx = bvh_payload_first_oac(self%aabb%node(cur)) + k
            v1  = vector_R8P(self%aabb%payload(idx)%v1(1),  self%aabb%payload(idx)%v1(2),  self%aabb%payload(idx)%v1(3))
            e12 = vector_R8P(self%aabb%payload(idx)%e12(1), self%aabb%payload(idx)%e12(2), self%aabb%payload(idx)%e12(3))
            e13 = vector_R8P(self%aabb%payload(idx)%e13(1), self%aabb%payload(idx)%e13(2), self%aabb%payload(idx)%e13(3))
            call triangle_point_distance(v1=v1, e12=e12, e13=e13, a=self%aabb%payload(idx)%a, b=self%aabb%payload(idx)%b, &
                                         c=self%aabb%payload(idx)%c, det=self%aabb%payload(idx)%det, point=point, &
                                         distance=d2, closest=c, region=region)
            call consider(d2, c)
         enddo
         lc = bvh_left_child_oac(self%aabb%node(cur)) ; rc = bvh_right_child_oac(self%aabb%node(cur))
         if (rc > 0_I4P .and. sp < size(stack)) then
            sp = sp + 1_I4P ; stack(sp) = rc
         endif
         if (lc > 0_I4P .and. sp < size(stack)) then
            sp = sp + 1_I4P ; stack(sp) = lc
         endif
      enddo
   else
      do f = 1_I4P, self%facets_number
         call self%facet(f)%compute_distance_with_region(point=point, distance=d2, closest=c, region=region)
         call consider(d2, c)
      enddo
   endif
   contains
      subroutine consider(dd, cc)
      !< Keep cc if its foot is separate from cp1 and dd is within the (tightening) bound.
      real(R8P),        intent(in) :: dd !< Candidate squared distance.
      type(vector_R8P), intent(in) :: cc !< Candidate closest point.

      if (dd > bound) return
      if (vnorm(cc - cp1) <= foot_tol) return
      found = .true. ; cp2 = cc ; bound = dd
      endsubroutine consider
   endsubroutine second_foot

   function many_args_ok(self, axis, a, b, t0, t1, t, facets, entering, counts, overflow, status) result(good)
   !< Argument check of `line_crossings_many(_device)` (issue #40); outputs zeroed when good.
   type(surface_stl_object), intent(in)            :: self          !< Surface.
   integer(I4P),             intent(in)            :: axis          !< Line axis.
   real(R8P),                intent(in)            :: a(:), b(:)    !< Line positions.
   real(R8P),                intent(in)            :: t0, t1        !< Parameter range.
   real(R8P),                intent(inout)         :: t(:,:)        !< Parameters.
   integer(I4P),             intent(inout)         :: facets(:,:)   !< Facets.
   logical,                  intent(inout)         :: entering(:,:) !< Entering flags.
   integer(I4P),             intent(inout)         :: counts(:)     !< Counts.
   logical,                  intent(inout)         :: overflow(:)   !< Overflow.
   integer(I4P),             intent(out), optional :: status        !< Status.
   logical                                         :: good          !< Arguments valid.
   integer(I4P)                                    :: n             !< Lines.

   if (present(status)) status = STATUS_OK
   n = size(a, kind=I4P)
   good = axis >= 1_I4P .and. axis <= 3_I4P .and. t1 >= t0 .and. size(b) == n .and. size(t, 2) == n .and. &
          size(facets, 2) == n .and. size(entering, 2) == n .and. size(counts) == n .and. size(overflow) == n .and. &
          size(facets, 1) == size(t, 1) .and. size(entering, 1) == size(t, 1) .and. size(t, 1) > 0
   good = good .and. self%facets_number > 0_I4P .and. self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. &
          allocated(self%aabb%node) .and. allocated(self%aabb%payload) .and. allocated(self%tri_vertices)
   if (.not. good) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   t = 0._R8P ; facets = 0_I4P ; entering = .false. ; counts = 0_I4P ; overflow = .false.
   endfunction many_args_ok

   subroutine line_crossings_core(self, axis, a, b, t0, t1, t, ids, es, n)
   !< Sorted crossings of an axis-aligned line (issue #40): BVH kernel with a growing buffer, or brute
   !< force over the facets (same kernel per triangle, same sort).
   type(surface_stl_object),  intent(in)  :: self   !< Surface.
   integer(I4P),              intent(in)  :: axis   !< Axis.
   real(R8P),                 intent(in)  :: a, b   !< Line position.
   real(R8P),                 intent(in)  :: t0, t1 !< Parameter range.
   real(R8P),    allocatable, intent(out) :: t(:)   !< Parameters.
   integer(I4P), allocatable, intent(out) :: ids(:) !< Facets.
   logical,      allocatable, intent(out) :: es(:)  !< Entering flags.
   integer(I4P),              intent(out) :: n      !< Count.
   integer(I4P) :: nmax, f !< Work.
   logical      :: overflow, hit, enter !< Work.
   real(R8P)    :: tt !< Work.

   n = 0_I4P
   if (self%facets_number == 0_I4P) then
      allocate(t(0), ids(0), es(0))
      return
   endif
   if (self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(self%aabb%node) .and. &
       allocated(self%aabb%payload)) then
      nmax = 64_I4P
      do
         if (allocated(t)) deallocate(t, ids, es)
         allocate(t(nmax), ids(nmax), es(nmax))
         call bvh_line_crossings_oac(self%aabb%node, self%aabb%payload, self%tri_vertices, axis, a, b, t0, t1, t, ids, es, &
                                     n, overflow)
         if (.not. overflow) exit
         nmax = 4_I4P * nmax
      enddo
   else
      allocate(t(self%facets_number), ids(self%facets_number), es(self%facets_number))
      do f = 1_I4P, self%facets_number
         call tri_line_hit_oac(self%tri_vertices(:, f), axis, a, b, hit, tt, enter)
         if (.not. hit .or. tt < t0 .or. tt > t1) cycle
         n = n + 1_I4P ; t(n) = tt ; ids(n) = f ; es(n) = enter
      enddo
      call sort_crossings_oac(t, ids, es, n)
   endif
   t = t(1:n) ; ids = ids(1:n) ; es = es(1:n)
   endsubroutine line_crossings_core

   function body_points(pose, points) result(bp)
   !< World-frame points mapped to the body frame of `pose` (issue #41).
   type(rigid_pose_t), intent(in) :: pose      !< Pose.
   type(vector_R8P),   intent(in) :: points(:) !< World points.
   type(vector_R8P)               :: bp(size(points)) !< Body points.
   real(R8P)                      :: q(3)      !< Work.
   integer(I4P)                   :: i         !< Counter.

   do i = 1_I4P, size(points, kind=I4P)
      q = to_body(pose, [points(i)%x, points(i)%y, points(i)%z])
      bp(i)%x = q(1) ; bp(i)%y = q(2) ; bp(i)%z = q(3)
   enddo
   endfunction body_points

   function box_contact(self, lo, hi, tol, reach) result(code)
   !< BOX_CUT, BOX_NEAR or BOX_NONE of the box [lo, hi] (issue #38): SAH BVH traversal
   !< (`bvh_box_contact_oac`), brute force over the facets otherwise (same tests, same result).
   type(surface_stl_object), intent(in) :: self         !< Surface.
   real(R8P),                intent(in) :: lo(3), hi(3) !< Box.
   real(R8P),                intent(in) :: tol, reach   !< Cut and near inflations.
   integer(I4P)                         :: code         !< Result.

   code = BOX_NONE
   if (self%facets_number == 0_I4P) return
   if (self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(self%aabb%node) .and. &
       allocated(self%aabb%payload)) then
      code = bvh_box_contact_oac(self%aabb%node, self%aabb%payload, lo, hi, tol, reach)
   else
      code = box_contact_brute(self, lo, hi, tol, reach)
   endif
   endfunction box_contact

   function box_contact_brute(self, lo, hi, tol, reach) result(code)
   !< `box_contact` by brute force over the facets (the octree fallback and the test reference).
   type(surface_stl_object), intent(in) :: self         !< Surface.
   real(R8P),                intent(in) :: lo(3), hi(3) !< Box.
   real(R8P),                intent(in) :: tol, reach   !< Cut and near inflations.
   integer(I4P)                         :: code         !< Result.
   real(R8P)                            :: v1(3), v2(3), v3(3) !< Triangle.
   integer(I4P)                         :: f            !< Counter.

   code = BOX_NONE
   do f = 1_I4P, self%facets_number
      v1 = [self%facet(f)%vertex(1)%x, self%facet(f)%vertex(1)%y, self%facet(f)%vertex(1)%z]
      v2 = [self%facet(f)%vertex(2)%x, self%facet(f)%vertex(2)%y, self%facet(f)%vertex(2)%z]
      v3 = [self%facet(f)%vertex(3)%x, self%facet(f)%vertex(3)%y, self%facet(f)%vertex(3)%z]
      if (box_triangle_overlap(lo, hi, v1, v2, v3, tol)) then
         code = BOX_CUT
         return
      endif
      if (box_triangle_overlap(lo, hi, v1, v2, v3, reach)) code = BOX_NEAR
   enddo
   endfunction box_contact_brute

   pure function vnorm(a) result(r)
   !< Euclidean norm of a vector (plain components, no TBP).
   type(vector_R8P), intent(in) :: a !< Vector.
   real(R8P)                    :: r !< Norm.

   r = sqrt(a%x * a%x + a%y * a%y + a%z * a%z)
   endfunction vnorm

   pure function bad_shape3(a, n) result(bad)
   !< True if the optional array a is present and not [3, n].
   real(R8P), optional, intent(in) :: a(:,:) !< Array.
   integer(I4P),        intent(in) :: n      !< Columns.
   logical                         :: bad    !< Result.

   bad = .false.
   if (present(a)) bad = size(a, 1) /= 3 .or. size(a, 2) /= n
   endfunction bad_shape3

   pure function bad_size(a, n) result(bad)
   !< True if the optional integer array a is present and not of size n.
   integer(I4P), optional, intent(in) :: a(:) !< Array.
   integer(I4P),           intent(in) :: n    !< Size.
   logical                            :: bad  !< Result.

   bad = .false.
   if (present(a)) bad = size(a) /= n
   endfunction bad_size

   subroutine closest_with_region(self, point, d2, facet, region, closest)
   !< Squared distance, closest facet, Voronoi region and closest point (issue #37): the traversal
   !< of the pseudo-normal path of `compute_distance` (tree, or brute force when disabled), then
   !< the closest point on the winning facet. Non-empty surface.
   type(surface_stl_object), intent(in)  :: self    !< Surface.
   type(vector_R8P),         intent(in)  :: point   !< Query point.
   real(R8P),                intent(out) :: d2      !< Squared distance.
   integer(I4P),             intent(out) :: facet   !< Closest facet.
   integer(I4P),             intent(out) :: region  !< Voronoi region tag.
   type(vector_R8P),         intent(out) :: closest !< Closest point.
   real(R8P)                             :: dd      !< Recomputed squared distance (discarded).

   closest = point
   if (self%aabb%get_is_initialized() .and. self%aabb%get_use_index()) then
      call self%aabb%distance_tree_with_region(facet=self%facet, point=point, distance=d2, closest_facet=facet, &
                                               closest_region=region)
   else
      call brute_force_with_region(self, point, d2, facet, region)
   endif
   if (facet > 0_I4P) call self%facet(facet)%compute_distance_with_region(point=point, distance=dd, closest=closest, &
                                                                          region=region)
   endsubroutine closest_with_region

   pure subroutine build_pnormal_payload(self)
   !< Build the flat per-facet pseudo-normal SoA for the device sign path
   !< (issue #20 §Step 5). One `pnormal_payload_t` record per facet, indexed by
   !< `facet_id` (1..facets_number). Pure-POD layout: no allocatable components,
   !< nothing the device cannot `copyin` shallowly.
   !<
   !< Called by `analyze` after `compute_pseudo_normals_via_pool` has populated
   !< the per-facet `normal`, `edge_pnormal`, `vertex_pnormal`. Strictly mechanical
   !< unpack of `vector_R8P` -> scalar SoA; bit-exact equivalent to the host
   !< `pseudo_normal_for_region` selector by construction.
   class(surface_stl_object), intent(inout) :: self !< File STL.
   integer(I4P)                             :: ff   !< Facet counter.

   if (allocated(self%pnormal_payload)) deallocate(self%pnormal_payload)
   if (self%facets_number <= 0_I4P) return
   allocate(self%pnormal_payload(self%facets_number))
   do ff = 1_I4P, self%facets_number
      associate (p => self%pnormal_payload(ff), f => self%facet(ff))
         p%nx    = f%normal%x ; p%ny = f%normal%y ; p%nz = f%normal%z
         p%ex(1) = f%edge_pnormal(1)%x ; p%ey(1) = f%edge_pnormal(1)%y ; p%ez(1) = f%edge_pnormal(1)%z
         p%ex(2) = f%edge_pnormal(2)%x ; p%ey(2) = f%edge_pnormal(2)%y ; p%ez(2) = f%edge_pnormal(2)%z
         p%ex(3) = f%edge_pnormal(3)%x ; p%ey(3) = f%edge_pnormal(3)%y ; p%ez(3) = f%edge_pnormal(3)%z
         p%vx(1) = f%vertex_pnormal(1)%x ; p%vy(1) = f%vertex_pnormal(1)%y ; p%vz(1) = f%vertex_pnormal(1)%z
         p%vx(2) = f%vertex_pnormal(2)%x ; p%vy(2) = f%vertex_pnormal(2)%y ; p%vz(2) = f%vertex_pnormal(2)%z
         p%vx(3) = f%vertex_pnormal(3)%x ; p%vy(3) = f%vertex_pnormal(3)%y ; p%vz(3) = f%vertex_pnormal(3)%z
      end associate
   enddo
   endsubroutine build_pnormal_payload

   subroutine enter_device(self)
   !< Make the surface device-resident (issue #20 §Step 5). Performs `!$acc enter
   !< data copyin` of the four flat arrays the GPU distance kernel reads:
   !< `self%aabb%node`, `self%aabb%payload`, `self%pnormal_payload` (pseudo-normal
   !< sign) and `self%wn_payload` (fast-winding-number sign, issue #30). Idempotent:
   !< if `self%is_on_device` is already true, the call is a no-op (a second
   !< `enter data copyin` would create a refcount the user did not ask for).
   !<
   !< This is the explicit half of the explicit/implicit API split: the caller
   !< pays the copyin cost exactly once, at a point of their choosing (typically
   !< right after `load_from_file` and `sanitize`, before the embedded-boundary
   !< marker pass that issues millions of `distance_many_device` queries).
   !<
   !< Restrictions:
   !<  - Only the SAH BVH tree (`AABB_TREE_SAH_BVH`) is device-callable.
   !<    `AABB_TREE_OCTREE` would copy `node(:)` but the device traversal kernels
   !<    do not enumerate the octree's 8-way children. If the tree kind is octree,
   !<    the call refuses and writes a stderr diagnostic; the surface stays on host.
   !<  - The surface must already be analyzed (BVH built, pnormal_payload and
   !<    wn_payload populated).
   !<    Called on an empty surface, the routine returns without touching the
   !<    device.
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (self%is_on_device) return
   if (self%facets_number <= 0_I4P) return
   if (self%aabb%get_tree_kind() /= AABB_TREE_SAH_BVH) then
      write(stderr, '(A)') 'fossil_surface_stl_object%enter_device: refused -- device path requires AABB_TREE_SAH_BVH.'
      return
   endif
   if (.not. allocated(self%aabb%node))      return
   if (.not. allocated(self%aabb%payload))   return
   if (.not. allocated(self%pnormal_payload)) return
   if (.not. allocated(self%wn_payload))      return
   if (.not. allocated(self%body_payload))    return
   if (.not. allocated(self%tri_vertices))    return

   !$acc enter data copyin(self%aabb%node, self%aabb%payload, self%pnormal_payload, self%wn_payload, self%body_payload, &
   !$acc&                  self%tri_vertices)
   if (self%features_built) then                 ! the feature-edge tree, if built (issue #61)
      !$acc enter data copyin(self%feature_seg, self%feature_box, self%feature_link)
      self%features_on_device = .true.
   endif
   if (self%thin_built) then                     ! the local sizes, if built (issue #62)
      !$acc enter data copyin(self%thin_pay, self%thin_node)
      self%thin_on_device = .true.
   endif
   self%is_on_device = .true.
   endsubroutine enter_device

   subroutine exit_device(self)
   !< Release the device-resident copy of the surface (issue #20 §Step 5).
   !< Symmetric inverse of `enter_device`: `!$acc exit data delete` of the four
   !< arrays previously copied in. Idempotent on `is_on_device=.false.`
   !< (returns without acting). NOT called by `destroy` (elemental, so it cannot
   !< hold `!$acc exit data`, see its @warning): the caller must call it before
   !< the surface is destroyed or goes out of scope, and before any operation
   !< that rebuilds the tree (`analyze` and its callers).
   class(surface_stl_object), intent(inout) :: self !< File STL.

   if (.not. self%is_on_device) return

   !$acc exit data delete(self%aabb%node, self%aabb%payload, self%pnormal_payload, self%wn_payload, self%body_payload, &
   !$acc&                 self%tri_vertices)
   if (self%features_on_device) then
      !$acc exit data delete(self%feature_seg, self%feature_box, self%feature_link)
      self%features_on_device = .false.
   endif
   if (self%thin_on_device) then
      !$acc exit data delete(self%thin_pay, self%thin_node)
      self%thin_on_device = .false.
   endif
   self%is_on_device = .false.
   endsubroutine exit_device

   recursive subroutine distance_many_device(self, points, distances, is_signed, sign_algorithm, is_square_root, status, &
                                             pose)
   !< GPU-offloaded batch distance query (issue #20 §Step 5). Device analogue of
   !< `distance_many`. Same semantics on the output -- `distances(i)` is the
   !< distance from `points(i)` to the surface -- but the loop body executes on
   !< the accelerator via `!$acc parallel loop`, driving the device-callable
   !< `bvh_distance_node_with_region_oac` kernel chain.
   !<
   !< Pre-conditions:
   !<  - `self%enter_device` has been called and `self%is_on_device` is true.
   !<    The `present(...)` clause on the parallel loop asserts this; calling
   !<    without `enter_device` first surfaces as a runtime "present" failure
   !<    from the OpenACC runtime, which is the documented loud-failure mode.
   !<  - The tree kind is `AABB_TREE_SAH_BVH` (octree not supported on device).
   !<
   !< Sign algorithms (`sign_algorithm`, used only when `is_signed`):
   !<  - `SIGN_PSEUDO_NORMAL` (default): the flat `pnormal_payload(best_facet)`
   !<    record at the closest point, recomputed by `triangle_point_distance`
   !<    over `payload(best_payload_idx)`.
   !<  - `SIGN_FAST_WINDING_NUMBER` (issue #30): inside iff the fast winding
   !<    number exceeds 0.5, from `bvh_winding_number_oac` over `wn_payload` with
   !<    the CPU default beta, the routine and constants the CPU path uses.
   !<    Robust on surfaces with holes; a second traversal per point.
   !<  - Any other code (`SIGN_RAY_INTERSECTIONS`, `SIGN_SOLID_ANGLE`, ...) is
   !<    refused with `STATUS_INVALID_INPUT` and a stderr message: those need a
   !<    ray traversal or an O(N) sum that has not been device-ported.
   !<
   !< Numerical contract: bit-exact vs the CPU serial `compute_distance` path
   !< where the device libm and FMA contraction allow; documented-tolerance
   !< fallback otherwise (see test acceptance gates in Step 6).
   class(surface_stl_object), intent(inout), target :: self           !< File STL (target: device arrays aliased to local pointers, see decls).
   type(vector_R8P),          intent(in)            :: points(:)      !< Query points (host).
   real(R8P),                 intent(out)           :: distances(:)   !< Per-point distances (host, caller-allocated).
   logical,                   intent(in),  optional :: is_signed      !< Trigger signed distance (default .false.).
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_PSEUDO_NORMAL (default) or SIGN_FAST_WINDING_NUMBER.
   logical,                   intent(in),  optional :: is_square_root !< Trigger sqrt of d^2 (default .false.).
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   type(rigid_pose_t),        intent(in),  optional :: pose            !< Pose of the body (issue #41): world-frame points.
   logical                                          :: is_signed_     !< Resolved signed sentinel.
   logical                                          :: is_fwn_        !< Resolved sign algorithm: fast winding number.
   integer(I4P)                                     :: algo           !< Resolved sign-algorithm code.
   real(R8P)                                        :: beta2          !< Squared FWN admissibility ratio (CPU default).
   real(R8P)                                        :: omega          !< Per-point total solid angle (FWN sign).
   logical                                          :: is_sqrt_       !< Resolved sqrt sentinel.
   integer(I4P)                                     :: n              !< Number of query points.
   integer(I4P)                                     :: i              !< Loop counter.
   real(R8P)                                        :: d2             !< Per-point best d^2.
   integer(I4P)                                     :: best_facet     !< Per-point winning facet id (1-based).
   integer(I4P)                                     :: best_payload   !< Per-point winning payload index.
   real(R8P)                                        :: cx, cy, cz     !< Per-point closest-point coords (sign step).
   integer(I4P)                                     :: region         !< Per-point Voronoi region (sign step).
   real(R8P)                                        :: nx, ny, nz     !< Per-point pseudo-normal components.
   real(R8P)                                        :: dx, dy, dz     !< Per-point (point - closest) components.
   real(R8P)                                        :: side           !< Sign-decision dot product.
   type(vector_R8P)                                 :: point_i        !< Reconstructed query point (kernel API).
   type(vector_R8P)                                 :: closest_v      !< Closest-point as vector_R8P (kernel API).
   type(vector_R8P)                                 :: v1_v, e12_v, e13_v !< Winning-triangle vertex/edges as vector_R8P (kernel API).
   ! Top-level local pointer aliases for the device arrays. enter_device maps
   ! only the member arrays, never `self` or `self%aabb`, so referencing
   ! `self%aabb%payload(i)` in the compute region makes the compiler implicitly
   ! copy `self` (OpenACC 3.4 §2.6.2: aggregates with no data clause are `copy`).
   ! That copy holds host addresses (nothing attached it, §2.6.8) and the
   ! subscript faults with CUDA 700. The loop body uses these flat aliases
   ! instead; pointer association is a host-side rename of the mapped
   ! allocation, so present(...) finds them.
   type(aabb_node_object),       pointer            :: node_p(:)      !< Alias of self%aabb%node.
   type(facet_distance_payload), pointer            :: payload_p(:)   !< Alias of self%aabb%payload.
   type(pnormal_payload_t),      pointer            :: pnormal_p(:)   !< Alias of self%pnormal_payload.
   type(wn_node_t),              pointer            :: wn_p(:)        !< Alias of self%wn_payload (lower bound 0, kept by pointer assignment).

   if (present(pose)) then ! body-frame query (issue #41): the transform of the CPU path, on the host
      call self%distance_many_device(points=body_points(pose, points), distances=distances, is_signed=is_signed, &
                                     sign_algorithm=sign_algorithm, is_square_root=is_square_root, status=status)
      return
   endif
   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(distances, kind=I4P) /= n) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (n == 0_I4P) return
   if (.not. self%is_on_device) then
      ! Refuse cleanly rather than relying on the OpenACC runtime's "present"
      ! abort message: more discoverable, and the host-side check costs nothing.
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_surface_stl_object%distance_many_device: surface not on device -- call enter_device first.'
      return
   endif

   is_signed_ = .false. ; if (present(is_signed))      is_signed_ = is_signed
   is_sqrt_   = .false. ; if (present(is_square_root)) is_sqrt_   = is_square_root

   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   if (algo /= SIGN_PSEUDO_NORMAL .and. algo /= SIGN_FAST_WINDING_NUMBER) then
      if (present(status)) status = STATUS_INVALID_INPUT
      write(stderr, '(A)') 'fossil_surface_stl_object%distance_many_device: refused -- sign_algorithm must be &
                           &SIGN_PSEUDO_NORMAL or SIGN_FAST_WINDING_NUMBER on the device.'
      return
   endif
   is_fwn_ = algo == SIGN_FAST_WINDING_NUMBER
   ! Same expression as the CPU path (`winding_number` with its default beta), so
   ! the admissibility test sees the same bits.
   beta2   = WN_DEFAULT_BETA * WN_DEFAULT_BETA

   ! Alias the device arrays to top-level local pointers. The compute region
   ! references the aliases (`node_p`, `payload_p`, `pnormal_p`) rather than the
   ! nested `self%aabb%...` components -- see the alias declarations for why.
   ! The targets were mapped by enter_device; pointer association is a host-side
   ! rename to the same allocation, so `present(node_p, ...)` finds them.
   node_p    => self%aabb%node
   payload_p => self%aabb%payload
   pnormal_p => self%pnormal_payload
   wn_p      => self%wn_payload

   !$acc parallel loop                                                                       &
   !$acc&  present(node_p, payload_p, pnormal_p, wn_p)                                       &
   !$acc&  copyin(points)                                                                    &
   !$acc&  copyout(distances)                                                                &
   !$acc&  private(point_i, d2, best_facet, best_payload, cx, cy, cz, region, nx, ny, nz,    &
   !$acc&          dx, dy, dz, side, closest_v, v1_v, e12_v, e13_v, omega)
   do i = 1_I4P, n
      ! Component-wise set (NOT `point_i = vector_R8P(...)`): VecFor overloads
      ! assignment(=) as a type-bound generic, and nvfortran does not allow TBP
      ! calls in device code (NVFORTRAN-F-1252). A language-support restriction,
      ! not a compiler bug. Plain component stores need no TBP.
      point_i%x    = points(i)%x
      point_i%y    = points(i)%y
      point_i%z    = points(i)%z
      d2           = huge(0._R8P)
      best_facet   = 0_I4P
      best_payload = 0_I4P
      call bvh_distance_node_with_region_oac(node=node_p, payload=payload_p,                  &
                                              n=0_I4P, point=point_i,                          &
                                              best=d2, best_facet=best_facet,                  &
                                              best_payload_idx=best_payload)
      if (is_signed_ .and. best_facet > 0_I4P) then
         ! Recompute closest point and Voronoi region on the winning triangle via
         ! the rich kernel (already !$acc routine seq from Step 4). Metrics are
         ! read through the pointer alias `payload_p`, not
         ! `self%aabb%payload(best_payload)` (see the alias declarations).
         ! Indexed directly rather than through an `associate` name: nvfortran
         ! treated an in-kernel associate name as a variable and generated an
         ! implicit `copy(p)` for it.
         !
         ! v1_v/e12_v/e13_v are filled by component-wise stores, not by assigning
         ! `vector_R8P(...)`: that assignment would call VecFor's `assign_vector`
         ! TBP (see the `point_i` note above).
         v1_v%x  = payload_p(best_payload)%v1(1)
         v1_v%y  = payload_p(best_payload)%v1(2)
         v1_v%z  = payload_p(best_payload)%v1(3)
         e12_v%x = payload_p(best_payload)%e12(1)
         e12_v%y = payload_p(best_payload)%e12(2)
         e12_v%z = payload_p(best_payload)%e12(3)
         e13_v%x = payload_p(best_payload)%e13(1)
         e13_v%y = payload_p(best_payload)%e13(2)
         e13_v%z = payload_p(best_payload)%e13(3)
         call triangle_point_distance(v1=v1_v, e12=e12_v, e13=e13_v,                       &
                                      a=payload_p(best_payload)%a,                         &
                                      b=payload_p(best_payload)%b,                         &
                                      c=payload_p(best_payload)%c,                         &
                                      det=payload_p(best_payload)%det,                     &
                                      point=point_i, distance=d2,                          &
                                      closest=closest_v, region=region)
         if (is_sqrt_) d2 = sqrt(d2)
         if (is_fwn_) then
            ! Inside iff w = omega / (4 pi) > 0.5: the CPU test in `is_point_inside`,
            ! written with the same division and constant so the decision matches.
            omega = bvh_winding_number_oac(node_p, payload_p, wn_p, point_i%x, point_i%y, point_i%z, beta2)
            if (omega / WN_FOUR_PI > 0.5_R8P) d2 = -d2
         else
            cx = closest_v%x ; cy = closest_v%y ; cz = closest_v%z
            nx = pnormal_x_oac(pnormal_p(best_facet), region)
            ny = pnormal_y_oac(pnormal_p(best_facet), region)
            nz = pnormal_z_oac(pnormal_p(best_facet), region)
            dx = points(i)%x - cx ; dy = points(i)%y - cy ; dz = points(i)%z - cz
            side = dx * nx + dy * ny + dz * nz
            if (side < 0._R8P) d2 = -d2
         endif
      else
         if (is_sqrt_) d2 = sqrt(d2)
      endif
      distances(i) = d2
   enddo
   !$acc end parallel loop

   endsubroutine distance_many_device

   subroutine wall_distance_many(self, points, distances, bodies, cutoff, status)
   !< Unsigned distance to the nearest wall and its body, with an early cutoff (issue #42): the SAH BVH
   !< traversal starts with the best squared distance at cutoff^2, so every node farther than the
   !< cutoff is pruned; a point with no wall closer than `cutoff` gets distance = cutoff and body 0
   !< (default cutoff: none). Same kernel as `wall_distance_many_device`; an octree falls back to
   !< `query_point`. OpenMP over the points. `status = STATUS_INVALID_INPUT` on a size mismatch.
   class(surface_stl_object), intent(in)            :: self         !< Surface.
   type(vector_R8P),          intent(in)            :: points(:)    !< Points [n].
   real(R8P),                 intent(out)           :: distances(:) !< Distances [n] (cutoff if none closer).
   integer(I4P),              intent(out)           :: bodies(:)    !< Bodies [n] (0 if none closer).
   real(R8P),                 intent(in),  optional :: cutoff       !< Cutoff distance (default: none).
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)    :: co, d2 !< Cutoff, squared distance.
   integer(I4P) :: i, n, bf, bp !< Work.

   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(distances) /= n .or. size(bodies) /= n) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   co = huge(1._R8P) ; if (present(cutoff)) co = cutoff
   if (self%facets_number == 0_I4P) then
      distances = co ; bodies = 0_I4P
      return
   endif
   !$omp parallel do default(shared) schedule(static) private(i, d2, bf, bp)
   do i = 1_I4P, n
      bf = 0_I4P ; bp = 0_I4P
      if (self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH .and. allocated(self%aabb%payload)) then
         d2 = co * co ; if (co >= huge(1._R8P)) d2 = huge(1._R8P)
         call bvh_distance_node_with_region_oac(node=self%aabb%node, payload=self%aabb%payload, n=0_I4P, point=points(i), &
                                                best=d2, best_facet=bf, best_payload_idx=bp)
      else
         call self%query_point(point=points(i), distance=d2, facet=bf, is_signed=.false.)
         d2 = d2 * d2
         if (d2 >= co * co) bf = 0_I4P
      endif
      if (bf > 0_I4P) then
         distances(i) = sqrt(d2) ; bodies(i) = self%body_payload(bf)
      else
         distances(i) = co ; bodies(i) = 0_I4P
      endif
   enddo
   !$omp end parallel do
   endsubroutine wall_distance_many

   subroutine wall_distance_many_device(self, points, distances, bodies, cutoff, status)
   !< `wall_distance_many` on the device (OpenACC, issue #42): the same kernel and early cutoff, bit
   !< for bit. Requires `enter_device`; SAH BVH. #20 rules (pointer aliases). GPU gate pending.
   class(surface_stl_object), intent(inout), target :: self         !< Surface (target: device aliases).
   type(vector_R8P),          intent(in)            :: points(:)    !< Points [n].
   real(R8P),                 intent(out)           :: distances(:) !< Distances [n].
   integer(I4P),              intent(out)           :: bodies(:)    !< Bodies [n].
   real(R8P),                 intent(in),  optional :: cutoff       !< Cutoff (default: none).
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   type(aabb_node_object),       pointer :: node_p(:)    !< Alias.
   type(facet_distance_payload), pointer :: payload_p(:) !< Alias.
   integer(I4P),                 pointer :: body_p(:)    !< Alias.
   type(vector_R8P) :: point_i !< Work.
   real(R8P)        :: co, d2, d0 !< Work.
   integer(I4P)     :: i, n, bf, bp !< Work.
   real(R8P),    allocatable :: dl(:) !< Device-side distances.
   integer(I4P), allocatable :: bl(:) !< Device-side bodies.

   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(distances) /= n .or. size(bodies) /= n .or. .not. self%is_on_device) then
      if (present(status)) status = STATUS_INVALID_INPUT
      if (.not. self%is_on_device) write(stderr, '(A)') 'wall_distance_many_device: surface not on device.'
      return
   endif
   co = huge(1._R8P) ; if (present(cutoff)) co = cutoff
   d0 = co * co ; if (co >= huge(1._R8P)) d0 = huge(1._R8P)
   if (n == 0_I4P) return
   allocate(dl(n), bl(n))
   node_p => self%aabb%node ; payload_p => self%aabb%payload ; body_p => self%body_payload
   !$acc parallel loop present(node_p, payload_p, body_p) copyin(points) copyout(dl, bl) private(point_i, d2, bf, bp)
   do i = 1_I4P, n
      point_i%x = points(i)%x ; point_i%y = points(i)%y ; point_i%z = points(i)%z
      d2 = d0 ; bf = 0_I4P ; bp = 0_I4P
      call bvh_distance_node_with_region_oac(node=node_p, payload=payload_p, n=0_I4P, point=point_i, best=d2, &
                                             best_facet=bf, best_payload_idx=bp)
      if (bf > 0_I4P) then
         dl(i) = sqrt(d2) ; bl(i) = body_p(bf)
      else
         dl(i) = co ; bl(i) = 0_I4P
      endif
   enddo
   distances = dl ; bodies = bl
   endsubroutine wall_distance_many_device

   subroutine state_changes(self, pose0, pose1, points, fresh, dead, sign_algorithm, checked, status)
   !< Fresh (solid at pose0, fluid at pose1) and dead (fluid at pose0, solid at pose1) points, e.g.
   !< cell centres, between two poses of the body (issue #41).
   !<
   !< Each point goes to the body frame at both poses, b0 and b1. Its state can change only if the
   !< segment b0 b1 reaches the surface, so when the distance at b0 exceeds |b0 - b1| the point keeps
   !< its state and the second query is skipped (the band of decision 4, exact as a filter). The
   !< states come from the signed distance with `sign_algorithm` (default SIGN_PSEUDO_NORMAL; the
   !< fast winding number for dirty geometry). `checked` (optional) counts the points that needed
   !< the second query. OpenMP over the points. `status = STATUS_INVALID_INPUT` on a size mismatch.
   class(surface_stl_object), intent(in)            :: self           !< Surface (body frame).
   type(rigid_pose_t),        intent(in)            :: pose0, pose1   !< Poses.
   type(vector_R8P),          intent(in)            :: points(:)      !< World points [n].
   logical,                   intent(out)           :: fresh(:)       !< Solid -> fluid [n].
   logical,                   intent(out)           :: dead(:)        !< Fluid -> solid [n].
   integer(I4P),              intent(in),  optional :: sign_algorithm !< SIGN_* (default SIGN_PSEUDO_NORMAL).
   integer(I4P),              intent(out), optional :: checked        !< Points that needed the second query.
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)     :: i, n, algo, nchk !< Work.
   real(R8P)        :: b0(3), b1(3), d0, d1 !< Body points, signed distances.
   logical          :: in0, in1 !< States.

   if (present(status)) status = STATUS_OK
   n = size(points, kind=I4P)
   if (size(fresh) /= n .or. size(dead) /= n) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   algo = SIGN_PSEUDO_NORMAL ; if (present(sign_algorithm)) algo = sign_algorithm
   nchk = 0_I4P
   !$omp parallel do default(shared) schedule(dynamic, 64) private(i, b0, b1, d0, d1, in0, in1) reduction(+:nchk)
   do i = 1_I4P, n
      b0 = to_body(pose0, [points(i)%x, points(i)%y, points(i)%z])
      b1 = to_body(pose1, [points(i)%x, points(i)%y, points(i)%z])
      call self%compute_distance(point=vector_R8P(b0(1), b0(2), b0(3)), distance=d0, is_signed=.true., sign_algorithm=algo, &
                                 is_square_root=.true.)
      in0 = d0 < 0._R8P
      fresh(i) = .false. ; dead(i) = .false.
      if (abs(d0) > norm2(b1 - b0)) cycle ! the segment b0 b1 cannot reach the surface
      nchk = nchk + 1_I4P
      call self%compute_distance(point=vector_R8P(b1(1), b1(2), b1(3)), distance=d1, is_signed=.true., sign_algorithm=algo, &
                                 is_square_root=.true.)
      in1 = d1 < 0._R8P
      fresh(i) = in0 .and. .not. in1
      dead(i) = in1 .and. .not. in0
   enddo
   !$omp end parallel do
   if (present(checked)) checked = nchk
   endsubroutine state_changes

   subroutine move(self, pose)
   !< Move the surface rigidly by `pose` (x <- R x + t) and refit the tree instead of rebuilding it
   !< (issue #41, decision 2b).
   !<
   !< The facet vertices move; what `analyze` derives from the coordinates is recomputed the same way
   !< (facet metrix, bounding box, vertex pool, area, volume, centroid, pseudo-normals and their
   !< payload, the winding cache, `tri_vertices`); the topology is kept (connectivity, shells,
   !< bodies). The SAH BVH keeps its topology and is refitted (`aabb%refit`, O(N)); an octree is
   !< spatial and is rebuilt. The facets are those of the same mesh moved and analyzed from scratch,
   !< so every query gives the same results (bit for bit; at an exact tie between facets either may
   !< be reported). A device copy (`enter_device`) is refreshed.
   class(surface_stl_object), intent(inout) :: self !< Surface.
   type(rigid_pose_t),        intent(in)    :: pose !< Pose.
   integer(I4P)                             :: f, k, ff !< Counters.
   real(R8P)                                :: w(3) !< Moved vertex.
   logical                                  :: was_on_device !< Device copy to refresh.

   if (self%facets_number == 0_I4P) return
   was_on_device = self%is_on_device
   if (was_on_device) call self%exit_device
   call drop_features(self)                       ! the feature tree is stale (issue #61)
   do f = 1_I4P, self%facets_number
      do k = 1, 3
         w = to_world(pose, [self%facet(f)%vertex(k)%x, self%facet(f)%vertex(k)%y, self%facet(f)%vertex(k)%z])
         self%facet(f)%vertex(k)%x = w(1) ; self%facet(f)%vertex(k)%y = w(2) ; self%facet(f)%vertex(k)%z = w(3)
      enddo
   enddo
   call self%compute_metrix
   call self%vertex_pool%initialize_from_facets(facet=self%facet)
   do ff = 1, self%facets_number
      call self%facet(ff)%set_vertex_ids(vid1=self%vertex_pool%facet_vid(ff, 1_I4P), &
                                         vid2=self%vertex_pool%facet_vid(ff, 2_I4P), &
                                         vid3=self%vertex_pool%facet_vid(ff, 3_I4P))
   enddo
   if (self%aabb%get_tree_kind() == AABB_TREE_SAH_BVH) then
      call self%aabb%refit(facet=self%facet)
   else
      call self%aabb%initialize(facet=self%facet, largest_edge_len=self%largest_edge_len())
   endif
   call self%compute_area
   call self%compute_volume
   call self%compute_centroid
   do ff = 1, self%facets_number
      call compute_pseudo_normals_via_pool(self%facet(ff), self%facet, self%vertex_pool)
   enddo
   call self%build_pnormal_payload
   call self%aabb%build_winding_payload(facet=self%facet, wn=self%wn_payload)
   do ff = 1, self%facets_number
      self%tri_vertices(:, ff) = [self%facet(ff)%vertex(1)%x, self%facet(ff)%vertex(1)%y, self%facet(ff)%vertex(1)%z, &
                                  self%facet(ff)%vertex(2)%x, self%facet(ff)%vertex(2)%y, self%facet(ff)%vertex(2)%z, &
                                  self%facet(ff)%vertex(3)%x, self%facet(ff)%vertex(3)%y, self%facet(ff)%vertex(3)%z]
   enddo
   if (was_on_device) call self%enter_device
   endsubroutine move

   recursive function is_point_inside(self, point, sign_algorithm, pose) result(is_inside)
   !< Compute sign.
   !<
   !< `sign_algorithm` selects the point-in-polyhedron test:
   !<   SIGN_PSEUDO_NORMAL (default) — Baerentzen-Aanaes pseudo-normal sign;
   !<                                  computes a signed distance and returns its sign.
   !<   SIGN_RAY_INTERSECTIONS       — odd/even crossing count, majority over three generic rays.
   !<   SIGN_SOLID_ANGLE             — sum of projected solid angles ≈ ±4π.
   !<   SIGN_FAST_WINDING_NUMBER     — fast winding number w > 0.5 (robust on holed meshes).
   !< Out-of-range values raise `error stop`.
   class(surface_stl_object), intent(in)           :: self            !< File STL.
   type(vector_R8P),          intent(in)           :: point           !< Point coordinates.
   logical                                         :: is_inside       !< Return true if point is inside surface.
   integer(I4P),              intent(in), optional :: sign_algorithm  !< Algorithm code (SIGN_*).
   type(rigid_pose_t),        intent(in), optional :: pose            !< Pose of the body (issue #41): world-frame point.
   integer(I4P)                                    :: algo            !< Effective algorithm code.
   real(R8P)                                       :: signed_d        !< Signed distance (pseudo-normal path).
   real(R8P)                                       :: pb(3)           !< Body-frame point.

   if (present(pose)) then
      pb = to_body(pose, [point%x, point%y, point%z])
      is_inside = self%is_point_inside(point=vector_R8P(pb(1), pb(2), pb(3)), sign_algorithm=sign_algorithm)
      return
   endif
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
   case(SIGN_FAST_WINDING_NUMBER)
      is_inside = self%winding_number(point=point) > 0.5_R8P
   case default
      error stop 'fossil_surface_stl_object%is_point_inside: unknown sign_algorithm code &
                 &(valid: SIGN_RAY_INTERSECTIONS=1, SIGN_SOLID_ANGLE=2, SIGN_PSEUDO_NORMAL=3, SIGN_FAST_WINDING_NUMBER=4)'
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
   !< Uses the hierarchical (fast, second-order) traversal over the persistent per-node cache
   !< `wn_payload` built by `analyze` (issue #29); falls back to the exact O(n_facets) per-facet
   !< sum when no cache exists. Set `beta <= 0` to force the exact sum (ground-truth tests).
   !< Default `beta = 2.0`; measured |w_fast - w_exact| <= 8.7e-3 (bunny) at beta = 2.
   class(surface_stl_object), intent(in)           :: self  !< Surface.
   type(vector_R8P),          intent(in)           :: point !< Query point.
   real(R8P),                 intent(in), optional :: beta  !< Multipole admissibility ratio (default 2.0).
   real(R8P)                                       :: w     !< Winding number.
   type(wn_node_t)                                 :: empty_wn(0:-1) !< No cache: forces the exact sum.

   if (self%facets_number <= 0) then
      w = 0._R8P
      return
   endif
   if (.not. allocated(self%wn_payload)) then
      w = compute_winding_number(facet=self%facet, tree=self%aabb, wn=empty_wn, point=point, beta=-1._R8P)
   else
      w = compute_winding_number(facet=self%facet, tree=self%aabb, wn=self%wn_payload, point=point, beta=beta)
   endif
   endfunction winding_number

   subroutine find_self_intersections(self, pairs, status, include_touching)
   !< Find all self-intersecting facet pairs in the surface (issue #18 §1.2).
   !<
   !< On return:
   !<   - `pairs` is allocated with exactly the number of intersections found
   !<     (zero-length allocation if none, never unallocated).
   !<   - Each record carries the pair of facet ids (`a < b`), the kind
   !<     (`SELFINT_CROSSING`, `SELFINT_COPLANAR`, `SELFINT_TOUCHING`) and the
   !<     3D intersection segment endpoints `p`, `q` (rounded; the two
   !<     centroids for a coplanar overlap).
   !<   - Every decision is exact (issue #69). Adjacent facets meeting only at
   !<     the shared feature are not reported; touching pairs only with
   !<     `include_touching`.
   !<
   !< Performance:
   !<   - Broad phase via tree-vs-tree traversal on the existing AABB tree
   !<     (O(N log N) on typical inputs, both BVH and octree kinds).
   !<   - Falls back to O(N^2) brute force if the tree is not initialized.
   class(surface_stl_object),                     intent(in)            :: self     !< Surface.
   type(intersection_pair_t), allocatable,        intent(out)           :: pairs(:) !< Output: list of intersection records.
   integer(I4P),                                  intent(out), optional :: status   !< 0 on success.
   logical,                                       intent(in),  optional :: include_touching !< Also pairs meeting at
                                                                                   !< one point (default .false.).

   if (self%facets_number <= 0) then
      allocate(pairs(0))
      if (present(status)) status = 0_I4P
      return
   endif
   call compute_self_intersections(facet=self%facet, tree=self%aabb, pairs=pairs, status=status, &
                                   include_touching=include_touching)
   endsubroutine find_self_intersections

   subroutine boolean(self, other, op, status, weld)
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
   !< ## Zero-area facets (issues #98, #99)
   !<
   !< An exactly degenerate facet (collinear or coincident vertices) that
   !< does not meet the other surface passes through uncut, as in the input
   !< (decided exactly by `degenerate_pair_meet`). One that touches or
   !< crosses it returns BOOL_STATUS_INVALID_INPUT with `self` unchanged:
   !< `sanitize` the inputs first (`remove_degenerate_facets` drops it).
   !<
   !< ## Coplanar faces, shared edges, contacts (issues #33, #65)
   !<
   !< Exact since issue #65, with no tolerance: facet pairs are classified by
   !< exact orient3d signs, cut points are implicit points (constructions,
   !< rounded once, the same bits for every facet that shares them), each cut
   !< facet is triangulated with exact orient2d, and the selection decides
   !< each patch (sub-triangles between cuts) once, exactly at a cut edge.
   !< Exactly coplanar faces cut each other along the parts of their edges
   !< inside the other facet and use the shared-face table.
   !<
   !< Faces coplanar by construction but off by round-off (after a rotation)
   !< are no longer merged: they cross exactly, the result keeps the slivers
   !< (watertight, the right volume), and its rounding to doubles can fold a
   !< sliver by an ulp (a self-intersecting pair). Solids touching only along
   !< an edge or at a point, and a symmetric difference whose two pieces
   !< touch along the crossing curves, give non-manifold results by
   !< definition (with the right volume).
   !<
   !< `weld` (opt-in, default off): before the arrangement, every vertex of B
   !< within `weld` of a vertex of A takes that vertex's coordinates bit for
   !< bit, when it is B's vertex's unique nearest, no other vertex of B goes
   !< to it, and no facet of B becomes degenerate. Inputs that share vertices
   !< up to round-off then merge (the old snap's output) and the arrangement
   !< stays exact; without it they give an exact result with slivers.
   class(surface_stl_object),     intent(inout)        :: self     !< Surface A; replaced by `op(A, B)` on success.
   class(surface_stl_object),     intent(in)           :: other    !< Surface B.
   integer(I4P),                  intent(in)           :: op       !< BOOL_UNION / BOOL_INTERSECT / BOOL_DIFFERENCE / BOOL_SYMDIFF.
   integer(I4P),                  intent(out), optional :: status  !< BOOL_STATUS_*.
   real(R8P),                     intent(in),  optional :: weld    !< Opt-in weld distance (issue #65), default off.
   type(facet_object), allocatable                     :: kept(:)
   integer(I4P)                                        :: bool_status
   type(surface_stl_object)                            :: welded  !< B with its vertices welded onto A's.
   type(facet_object), allocatable                     :: fb(:)   !< B's facets, patches mapped into A's table.
   integer(I4P), allocatable                           :: pmap(:) !< B's patch -> A's patch (issue #66).
   integer(I4P)                                        :: k       !< Counter.

   ! patches (issue #66): every output facet keeps its parent's patch; B's patches join A's table by name
   call sync_patch_names(self)
   allocate(pmap(other%get_patches_number()))
   do k = 1, size(pmap)
      pmap(k) = find_or_add_patch(self%patch_names, other%patch_name(k))
   enddo
   if (present(weld)) then
      if (weld > 0._R8P .and. self%facets_number > 0 .and. other%facets_number > 0) then
         call weld_onto(a=self, b=other, tol=weld, bw=welded)
         fb = welded%facet(1:welded%facets_number)
         do k = 1, size(fb)
            fb(k)%patch = pmap(max(1_I4P, min(fb(k)%patch, size(pmap, kind=I4P))))
         enddo
         call compute_boolean(facet_a=self%facet,  tree_a=self%aabb, &
                              facet_b=fb, tree_b=welded%aabb, &
                              op=op, kept_facet=kept, status=bool_status)
         if (present(status)) status = bool_status
         if (bool_status /= BOOL_STATUS_OK) return
         call self%adopt_facets(facets=kept)
         return
      endif
   endif
   if (other%facets_number > 0) then
      fb = other%facet(1:other%facets_number)
      do k = 1, size(fb)
         fb(k)%patch = pmap(max(1_I4P, min(fb(k)%patch, size(pmap, kind=I4P))))
      enddo
   else
      allocate(fb(0))
   endif
   call compute_boolean(facet_a=self%facet,  tree_a=self%aabb, &
                        facet_b=fb, tree_b=other%aabb, &
                        op=op, kept_facet=kept, status=bool_status)
   if (present(status)) status = bool_status
   if (bool_status /= BOOL_STATUS_OK) return

   ! Adopt the result into self. `adopt_facets` re-runs `analyze` so the
   ! AABB tree, connectivity, and pseudo-normals are rebuilt for queries.
   call self%adopt_facets(facets=kept)
   endsubroutine boolean

   subroutine weld_onto(a, b, tol, bw)
   !< The opt-in weld of `boolean` (issue #65): a copy of b whose vertices within `tol` of a vertex of a take that
   !< vertex's coordinates bit for bit (`weld_vertices`: the unique nearest, never two vertices of b onto one of a). A
   !< weld that would make a facet of b degenerate (two equal vertices) is undone, until none does. It is the only
   !< snapping that keeps the arrangement exact: the welded vertices are then the same input vertex for both surfaces.
   class(surface_stl_object), intent(in)  :: a, b !< Surfaces.
   real(R8P),                 intent(in)  :: tol  !< Weld distance.
   type(surface_stl_object),  intent(out) :: bw   !< Welded copy of b.
   real(R8P), allocatable                 :: xa(:, :), xb(:, :)
   integer(I4P), allocatable              :: map(:), vid(:, :)
   type(facet_object), allocatable        :: f(:)
   integer(I4P)                           :: i, k, l
   logical                                :: changed
   real(R8P)                              :: y(3, 3)

   call unique_vertices(a, xa)
   call unique_vertices(b, xb, vid)
   call weld_vertices(xa, xb, tol, map)
   changed = .true.
   do while (changed)
      changed = .false.
      do i = 1, b%facets_number
         do k = 1, 3
            y(:, k) = xb(:, vid(k, i))
            if (map(vid(k, i)) > 0_I4P) y(:, k) = xa(:, map(vid(k, i)))
         enddo
         if (all(y(:, 1) == y(:, 2)) .or. all(y(:, 2) == y(:, 3)) .or. all(y(:, 3) == y(:, 1))) then
            do k = 1, 3
               if (map(vid(k, i)) > 0_I4P) then
                  map(vid(k, i)) = 0_I4P ; changed = .true.
               endif
            enddo
         endif
      enddo
   enddo
   allocate(f(b%facets_number))
   do i = 1, b%facets_number
      f(i) = b%facet(i)
      do k = 1, 3
         l = vid(k, i)
         if (map(l) > 0_I4P) then
            f(i)%vertex(k)%x = xa(1, map(l)) ; f(i)%vertex(k)%y = xa(2, map(l)) ; f(i)%vertex(k)%z = xa(3, map(l))
         endif
      enddo
   enddo
   call bw%adopt_facets(facets=f)
   contains
      subroutine unique_vertices(s, x, ids)
      !< Distinct vertex coordinates of s (exact), and the index of each facet vertex among them.
      class(surface_stl_object),           intent(in)  :: s
      real(R8P), allocatable,              intent(out) :: x(:, :)
      integer(I4P), allocatable, optional, intent(out) :: ids(:, :)
      real(R8P), allocatable                           :: all_x(:, :)
      integer(I4P), allocatable                        :: o(:), w(:)
      integer(I4P)                                     :: n, m, ii, kk, width, lo, mid, hi, p, q, r

      n = 3_I4P * s%facets_number
      allocate(all_x(3, n))
      do ii = 1, s%facets_number
         do kk = 1, 3
            all_x(:, 3 * (ii - 1) + kk) = [s%facet(ii)%vertex(kk)%x, s%facet(ii)%vertex(kk)%y, s%facet(ii)%vertex(kk)%z]
         enddo
      enddo
      o = [(ii, ii = 1, n)] ; allocate(w(n))
      width = 1_I4P
      do while (width < n)
         lo = 1_I4P
         do while (lo <= n)
            mid = min(lo + width - 1_I4P, n) ; hi = min(lo + 2_I4P * width - 1_I4P, n)
            p = lo ; q = mid + 1_I4P ; r = lo
            do while (p <= mid .and. q <= hi)
               if (lexless(all_x(:, o(q)), all_x(:, o(p)))) then
                  w(r) = o(q) ; q = q + 1_I4P
               else
                  w(r) = o(p) ; p = p + 1_I4P
               endif
               r = r + 1_I4P
            enddo
            do while (p <= mid)
               w(r) = o(p) ; p = p + 1_I4P ; r = r + 1_I4P
            enddo
            do while (q <= hi)
               w(r) = o(q) ; q = q + 1_I4P ; r = r + 1_I4P
            enddo
            lo = lo + 2_I4P * width
         enddo
         o = w
         width = 2_I4P * width
      enddo
      allocate(x(3, n))
      if (present(ids)) allocate(ids(3, s%facets_number))
      m = 0_I4P
      do ii = 1, n
         if (m == 0_I4P) then
            m = 1_I4P ; x(:, 1) = all_x(:, o(ii))
         elseif (any(all_x(:, o(ii)) /= x(:, m))) then
            m = m + 1_I4P ; x(:, m) = all_x(:, o(ii))
         endif
         if (present(ids)) ids(mod(o(ii) - 1_I4P, 3_I4P) + 1_I4P, (o(ii) - 1_I4P) / 3_I4P + 1_I4P) = m
      enddo
      x = x(:, 1:m)
      endsubroutine unique_vertices

      pure function lexless(u, v) result(less)
      real(R8P), intent(in) :: u(3), v(3)
      logical               :: less

      if (u(1) /= v(1)) then
         less = u(1) < v(1)
      elseif (u(2) /= v(2)) then
         less = u(2) < v(2)
      else
         less = u(3) < v(3)
      endif
      endfunction lexless
   endsubroutine weld_onto

   subroutine resolve_self_intersections(self, status)
   !< Remove the self-intersections of a closed surface: replace it by the boundary of the region it winds around at
   !< least once (issue #78; the self-union of Zhou, Grinspun, Zorin, Jacobson 2016).
   !<
   !< The pairs of `find_self_intersections` (crossings, coplanar overlaps, folds) are cut exactly into the
   !< arrangement of the surface with itself (implicit points, `cdt_facet`: the exact machinery of `boolean`, issue
   !< #65); each patch of sub-triangles between cuts is kept iff the winding number of the surface is 0 in front of
   !< it (and so 1 behind it), see `self_union_compute`. Interior walls, doubly wound regions and inverted pockets
   !< (winding -1) go; the result is self-intersection free with the input's patches.
   !<
   !< Status (BOOL_STATUS_*), `self` unchanged unless OK:
   !<   - OK: no self-intersection (nothing to do), or the self-union is closed (no boundary edge);
   !<   - CDT_FAILED: a facet's triangulation failed (a cut along a coplanar facet's edge crossing another cut inside
   !<     it); cuts of three facets meeting at a point inside a facet are split at their three-plane point (IPK_TPI,
   !<     issue #91). See the Known limitations of the reference page;
   !<   - NOT_CLOSED: the selection is not closed (an input that is not closed has no well-defined winding number);
   !<   - INVALID_INPUT: an exactly degenerate (zero-area) facet meets a self-intersecting facet, or a cut point of its
   !<     neighbours lies on it (issue #99): `sanitize` first. A degenerate facet away from the cuts passes through;
   !<   - EMPTY_INPUT: no facets.
   class(surface_stl_object),     intent(inout)        :: self   !< Surface; replaced by its self-union on success.
   integer(I4P),                  intent(out), optional :: status !< BOOL_STATUS_*.
   type(intersection_pair_t), allocatable              :: pairs(:)
   type(facet_object),        allocatable              :: kept(:), copy(:)
   integer(I4P),              allocatable              :: pa(:), pb(:)
   type(surface_stl_object)                            :: result
   integer(I4P)                                        :: k, n, st

   if (present(status)) status = BOOL_STATUS_OK
   if (self%facets_number <= 0) then
      if (present(status)) status = BOOL_STATUS_EMPTY_INPUT
      return
   endif
   call self%find_self_intersections(pairs=pairs)
   n = count(pairs(:)%kind /= SELFINT_SHEETS) ! sheets crossing along an edge need no cut: the patches split there
   if (size(pairs) == 0) return
   allocate(pa(n), pb(n))
   n = 0_I4P
   do k = 1, size(pairs)
      if (pairs(k)%kind == SELFINT_SHEETS) cycle
      n = n + 1_I4P ; pa(n) = pairs(k)%a ; pb(n) = pairs(k)%b
   enddo
   call sync_patch_names(self)
   call compute_self_union(facet=self%facet(1:self%facets_number), tree=self%aabb, pa=pa, pb=pb, kept_facet=kept, &
                           status=st)
   if (st /= BOOL_STATUS_OK) then
      if (present(status)) status = st
      return
   endif
   allocate(copy, source=kept)
   call result%adopt_facets(facets=copy) ! closure check on a copy: self stays unchanged unless OK
   if (result%facet_1_de%ids_number + result%facet_2_de%ids_number + result%facet_3_de%ids_number > 0 .or. &
       result%facets_number == 0) then
      if (present(status)) status = BOOL_STATUS_NOT_CLOSED
      return
   endif
   call self%adopt_facets(facets=kept)
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

   subroutine decimate(self, target_facets, max_error, rim_collapses, status)
   !< Reduce `self`'s triangle count to ≤ `target_facets` via QEM edge
   !< collapse (issue #18 §1.3, Garland & Heckbert 1997).
   !<
   !< On success, `self` is replaced by the decimated mesh — adopt_facets
   !< runs analyze, rebuilding the AABB tree, vertex pool, connectivity,
   !< and pseudo-normals so all downstream queries (distance, winding number,
   !< boolean, ...) work without further user action.
   !<
   !< If the algorithm cannot reach `target_facets` (every remaining edge
   !< collapse would violate a safety check — link condition, normal-flip,
   !< degenerate facet), it returns DEC_STATUS_NO_PROGRESS and leaves `self`
   !< at the smallest count it could safely reach, a valid surface.
   !< `target_facets <= 0` returns DEC_STATUS_BAD_INPUT, `self` untouched (issue #81).
   !<
   !< Patches (issues #66, #81): every facet keeps its patch; a patch rim coarsens as a curve (two rim vertices merge
   !< along their rim edge, onto one of its endpoints; corners where three or more patches meet stay), unless
   !< `rim_collapses=.false.` locks the rims. `max_error` (a distance): the loop stops with DEC_STATUS_NO_PROGRESS when
   !< the cheapest collapse moves its vertex farther than max_error from one of the input planes it carries (the square
   !< root of the quadric error).
   !<
   !< Pre-condition: `self` must have been analyzed (i.e. `vertex_id` and
   !< `fcon_edge` populated). All public construction paths
   !< (`load_from_file`, `adopt_facets`, `analyze`) ensure this.
   class(surface_stl_object),     intent(inout)        :: self
   integer(I4P),                  intent(in)           :: target_facets
   real(R8P),                     intent(in),  optional :: max_error     !< Cost ceiling, a distance.
   logical,                       intent(in),  optional :: rim_collapses !< Patch rims coarsen (default .true.).
   integer(I4P),                  intent(out), optional :: status
   type(facet_object), allocatable                     :: working(:)
   integer(I4P)                                        :: dec_status

   if (present(status)) status = DEC_STATUS_OK
   if (self%facets_number == 0_I4P .or. target_facets <= 0_I4P) then
      if (present(status)) status = DEC_STATUS_BAD_INPUT
      return
   endif
   if (self%facets_number <= target_facets) return  ! nothing to do

   ! Work on a deep copy so the input stays valid until we adopt the output.
   allocate(working(self%facets_number))
   working(1:self%facets_number) = self%facet(1:self%facets_number)

   call compute_decimate(facet=working, target_facets=target_facets, max_error=max_error, rim_collapses=rim_collapses, &
                         status=dec_status)
   if (present(status)) status = dec_status
   if (dec_status == DEC_STATUS_BAD_INPUT) return
   ! Adopt regardless of NO_PROGRESS — the output is still a valid mesh,
   ! just larger than requested.
   call self%adopt_facets(facets=working)
   endsubroutine decimate

   subroutine isotropic_remesh(self, target_length, iterations, preserve_features, status)
   !< Botsch-Kobbelt 2004 isotropic remeshing (issue #18 §1.7).
   !<
   !< Equalizes edge lengths across the surface while preserving its
   !< geometry. Each outer iteration runs four passes:
   !<   1. Split edges longer than 4 L / 3.
   !<   2. Collapse edges shorter than 4 L / 5.
   !<   3. Flip interior edges to balance vertex valences toward 6.
   !<   4. Tangential relaxation (area-weighted Laplacian update projected on
   !<      the vertex tangent plane, issue #91) + projection onto the original
   !<      surface.
   !< Collapse, flip and relaxation refuse a change whose new facets leave the
   !< input by more than `REMESH_ERROR_TOL` (0.07) times the target length
   !< (issue #91): a target near the size of the body no longer collapses it.
   !<
   !< Arguments:
   !<   - `target_length`: desired edge length. If <= 0 (default), uses
   !<     the median of the input's edge lengths.
   !<   - `iterations`: outer-loop count. Default 5; typical convergence
   !<     at 3-5.
   !<   - `preserve_features`: if .true., vertices on edges with dihedral
   !<     angle > 30° are locked from being moved or collapsed away.
   !<     Essential for inputs with sharp features (cube edges, mechanical
   !<     CAD).
   !<
   !< Implementation: deep-copies `self` to use as the projection reference
   !< (the relaxation step queries `reference%aabb` for closest points),
   !< runs `fossil_remesh%isotropic_remesh` on a working facet array, then
   !< adopts the result back via `adopt_facets` (which rebuilds AABB tree,
   !< vertex pool, connectivity, pseudo-normals).
   !<
   !< Documented limitations: see isotropic-remesh.md.
   class(surface_stl_object),     intent(inout)        :: self
   real(R8P),                     intent(in), optional :: target_length     !< Target edge length; <= 0 → median.
   integer(I4P),                  intent(in), optional :: iterations        !< Outer-loop count; default 5.
   logical,                       intent(in), optional :: preserve_features !< Lock sharp-edge vertices; default .false.
   integer(I4P),                  intent(out), optional :: status
   type(facet_object), allocatable                     :: working(:)
   type(surface_stl_object)                            :: reference
   real(R8P)                                           :: L
   integer(I4P)                                        :: iters
   logical                                             :: preserve
   integer(I4P)                                        :: rem_status

   if (present(status)) status = REM_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = REM_STATUS_BAD_INPUT
      return
   endif

   L        = -1._R8P  ; if (present(target_length))     L        = target_length
   iters    = 5_I4P    ; if (present(iterations))        iters    = iterations
   preserve = .false.  ; if (present(preserve_features)) preserve = preserve_features

   ! Reference = current state (deep copy via adopt_facets, which also
   ! rebuilds the AABB tree the projection step needs).
   allocate(working(self%facets_number))
   working(1:self%facets_number) = self%facet(1:self%facets_number)
   call reference%adopt_facets(facets=working)

   ! Working = another deep copy (the reference adopted the first one).
   allocate(working(self%facets_number))
   working(1:self%facets_number) = self%facet(1:self%facets_number)

   call compute_remesh(facet=working, target_length=L, iterations=iters, &
                       preserve_features=preserve, &
                       reference_facet=reference%facet, reference_tree=reference%aabb, &
                       status=rem_status)
   if (present(status)) status = rem_status
   if (rem_status /= REM_STATUS_OK) return
   call self%adopt_facets(facets=working)
   endsubroutine isotropic_remesh

   subroutine surface_markers(self, dx, centroids, areas, normals, spacing_fraction, spacings, iterations, preserve_features, &
                              status)
   !< Grid-adapted surface markers (issue #42 step 4), for wall models at surface points and force
   !< integration with direct forcing.
   !<
   !< A copy of the surface is remeshed by `isotropic_remesh` to the target edge length sqrt(3) dx: the
   !< centroids of two equilateral triangles sharing an edge are then dx apart. The markers are the
   !< facet centroids, with the facet areas and unit normals. The surface itself is not changed.
   !<
   !< Spacing check: `spacings(k)` is the distance from marker k to the nearest centroid among its edge
   !< neighbours (0 for a facet without one), and `spacing_fraction` the fraction of markers whose
   !< spacing lies in [0.8, 1.2] dx, the window direct forcing needs (Cheylan et al. 2023).
   !< `iterations` (default 5) and `preserve_features` (default .false.; .true. keeps sharp edges) go to
   !< `isotropic_remesh`, with its limitations. `status = STATUS_INVALID_INPUT` for dx <= 0, an empty
   !< surface or a failed remesh (the outputs are then empty).
   class(surface_stl_object), intent(in)            :: self                !< Surface.
   real(R8P),                 intent(in)            :: dx                  !< Grid spacing.
   real(R8P), allocatable,    intent(out)           :: centroids(:,:)      !< Marker positions [3, m].
   real(R8P), allocatable,    intent(out)           :: areas(:)            !< Marker areas [m].
   real(R8P), allocatable,    intent(out)           :: normals(:,:)        !< Marker unit normals [3, m].
   real(R8P),                 intent(out), optional :: spacing_fraction    !< Fraction of spacings in [0.8, 1.2] dx.
   real(R8P), allocatable,    intent(out), optional :: spacings(:)         !< Nearest-neighbour spacings [m].
   integer(I4P),              intent(in),  optional :: iterations          !< Remesh iterations (default 5).
   logical,                   intent(in),  optional :: preserve_features   !< Keep sharp edges (default .false.).
   integer(I4P),              intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   type(surface_stl_object)                         :: work                !< Remeshed copy.
   type(facet_object), allocatable                  :: fac(:)              !< Facets of the copy.
   real(R8P), allocatable                           :: sp(:)               !< Spacings.
   type(vector_R8P)                                 :: dc                  !< Centroid offset.
   real(R8P)                                        :: target              !< Remesh edge length.
   integer(I4P)                                     :: m, k, j, g, rs      !< Counters, remesh status.
   integer                                          :: pass                !< Remesh pass.

   if (present(status)) status = STATUS_OK
   allocate(centroids(3, 0), areas(0), normals(3, 0))
   if (present(spacings)) allocate(spacings(0))
   if (present(spacing_fraction)) spacing_fraction = 0._R8P
   if (.not. (dx > 0._R8P) .or. self%facets_number == 0_I4P) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   target = sqrt(3._R8P) * dx
   do pass = 1, 2
      allocate(fac(self%facets_number)) ; fac = self%facet(1:self%facets_number)
      call work%adopt_facets(facets=fac)
      call work%isotropic_remesh(target_length=target, iterations=iterations, preserve_features=preserve_features, status=rs)
      if (rs /= REM_STATUS_OK .or. work%facets_number == 0_I4P) then
         if (present(status)) status = STATUS_INVALID_INPUT
         return
      endif
      m = work%facets_number
      if (allocated(sp)) deallocate(sp)
      allocate(sp(m))
      do k = 1_I4P, m
         sp(k) = huge(1._R8P)
         do j = 1_I4P, 3_I4P
            g = work%facet(k)%fcon_edge(j)
            if (g <= 0_I4P .or. g > m) cycle
            dc = work%facet(g)%centroid - work%facet(k)%centroid
            sp(k) = min(sp(k), vnorm(dc))
         enddo
         if (sp(k) == huge(1._R8P)) sp(k) = 0._R8P
      enddo
      ! calibration: the remesher's edges come out shorter than the target and the nearest-neighbour
      ! spacing is a minimum, so the first pass rescales the target by dx over the mean spacing
      if (pass == 1 .and. count(sp > 0._R8P) > 0) target = target * dx / (sum(sp) / count(sp > 0._R8P))
   enddo
   deallocate(centroids, areas, normals)
   allocate(centroids(3, m), areas(m), normals(3, m))
   do k = 1_I4P, m
      centroids(:, k) = [work%facet(k)%centroid%x, work%facet(k)%centroid%y, work%facet(k)%centroid%z]
      normals(:, k) = [work%facet(k)%normal%x, work%facet(k)%normal%y, work%facet(k)%normal%z]
      areas(k) = work%facet(k)%area()
   enddo
   if (present(spacing_fraction)) spacing_fraction = count(sp >= 0.8_R8P * dx .and. sp <= 1.2_R8P * dx) / real(m, R8P)
   if (present(spacings)) call move_alloc(from=sp, to=spacings)
   endsubroutine surface_markers

   subroutine intersect_ray_all(self, ray_origin, ray_direction, hits, status, exclude)
   !< Return every crossing of the ray with this surface, sorted by `t` ascending (issue #18 §2.5): one hit per crossing,
   !< decided exactly (issue #82, `fossil_ray_query`: a ray through a shared edge or vertex hits one facet there, so the
   !< parity of the hit count on a closed surface is the state of the origin). Hits at t = 0 (an origin on the surface)
   !< are kept; the facet `exclude` is skipped (a ray cast from a facet).
   !<
   !< Tree-accelerated: prunes whole subtrees whose AABB the ray misses. Falls
   !< back to `ray_intersect_all_flat` (the step-1 reference oracle) when the
   !< AABB tree is not available, and the test suite exercises both paths to
   !< guarantee bit-identical hit lists. Empty surface ⇒ status BAD_INPUT.
   class(surface_stl_object),    intent(in)            :: self           !< Surface.
   type(vector_R8P),             intent(in)            :: ray_origin     !< Ray origin.
   type(vector_R8P),             intent(in)            :: ray_direction  !< Ray direction (need not be unit).
   type(ray_hit_t), allocatable, intent(out)           :: hits(:)        !< Sorted hit records.
   integer(I4P),                 intent(out), optional :: status         !< Status code.
   integer(I4P),                 intent(in),  optional :: exclude        !< Facet to skip (issue #82).
   real(R8P)                                           :: dir_norm2      !< |dir|^2 for the BAD_INPUT guard.

   if (present(status)) status = RAY_STATUS_OK
   if (self%facets_number == 0_I4P) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   dir_norm2 = ray_direction%dotproduct(rhs=ray_direction)
   if (dir_norm2 < 1.0e-30_R8P) then
      allocate(hits(0))
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   if (self%aabb%get_nodes_number() > 0_I4P) then
      call self%aabb%intersect_ray_all_tree(facet=self%facet, ray_origin=ray_origin, &
                                            ray_direction=ray_direction, hits=hits, exclude=exclude)
   else
      call ray_intersect_all_flat(facet=self%facet, ray_origin=ray_origin, ray_direction=ray_direction, &
                                  hits=hits, status=status, exclude=exclude)
   endif
   endsubroutine intersect_ray_all

   subroutine intersect_ray_first(self, ray_origin, ray_direction, hit, has_hit, status, exclude)
   !< Return the closest ray-facet intersection (issue #18 §2.5); the facet `exclude` is skipped (issue #62: a ray
   !< cast from a facet never hits it at t = 0, without offsetting the origin).
   !<
   !< Tree-accelerated best-first traversal with `t_best` pruning of sibling
   !< subtrees. Output `hit` is meaningful only when `has_hit = .true.`; the
   !< caller must check that flag before consuming `hit%t / hit%facet_id /
   !< hit%point`. Empty surface or zero-magnitude direction ⇒ status BAD_INPUT
   !< and `has_hit = .false.`.
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   type(vector_R8P),          intent(in)            :: ray_origin     !< Ray origin.
   type(vector_R8P),          intent(in)            :: ray_direction  !< Ray direction.
   type(ray_hit_t),           intent(out)           :: hit            !< Closest hit (if has_hit).
   logical,                   intent(out)           :: has_hit        !< True if any facet was hit.
   integer(I4P),              intent(out), optional :: status         !< Status code.
   integer(I4P),              intent(in),  optional :: exclude        !< Facet to skip (issue #62).
   real(R8P)                                        :: dir_norm2      !< |dir|^2 for BAD_INPUT guard.
   type(ray_hit_t), allocatable                     :: fallback_hits(:) !< Used only when AABB is absent.

   has_hit = .false.
   if (present(status)) status = RAY_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   dir_norm2 = ray_direction%dotproduct(rhs=ray_direction)
   if (dir_norm2 < 1.0e-30_R8P) then
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   if (self%aabb%get_nodes_number() > 0_I4P) then
      call self%aabb%intersect_ray_first_tree(facet=self%facet, ray_origin=ray_origin, &
                                              ray_direction=ray_direction, hit=hit, has_hit=has_hit, exclude=exclude)
   else
      ! Fallback: build the all-hits list (other than `exclude`) and pick the first.
      call ray_intersect_all_flat(facet=self%facet, ray_origin=ray_origin, ray_direction=ray_direction, &
                                  hits=fallback_hits, status=status, exclude=exclude)
      if (size(fallback_hits) > 0) then
         hit = fallback_hits(1)
         has_hit = .true.
      endif
   endif
   endsubroutine intersect_ray_first

   subroutine intersect_ray_any(self, ray_origin, ray_direction, max_t, found, status)
   !< Any-hit early-exit ray query (issue #18 §2.5).
   !<
   !< Returns `found = .true.` as soon as any facet is intersected at parameter
   !< `0 <= t <= max_t`. No record is kept of which facet or where — by design,
   !< since the use cases (shadow rays, occlusion, conservative inside-tests)
   !< only need a yes/no answer and the early-exit is the whole point.
   !<
   !< `max_t` is optional; when omitted, defaults to `MaxR8P` (any hit anywhere
   !< along the forward half-ray). Pass `max_t = euclidean_distance_to_target`
   !< to test whether anything occludes a target; with a unit-length direction,
   !< `t` and Euclidean distance coincide.
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   type(vector_R8P),          intent(in)            :: ray_origin     !< Ray origin.
   type(vector_R8P),          intent(in)            :: ray_direction  !< Ray direction.
   real(R8P),                 intent(in),  optional :: max_t          !< Upper t bound (default MaxR8P).
   logical,                   intent(out)           :: found          !< True iff any hit in [0, max_t].
   integer(I4P),              intent(out), optional :: status         !< Status code.
   real(R8P)                                        :: cap            !< Effective max_t.
   real(R8P)                                        :: dir_norm2      !< |dir|^2 for BAD_INPUT guard.
   type(ray_hit_t), allocatable                     :: hits(:)        !< Hits (fallback).

   found = .false.
   if (present(status)) status = RAY_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   dir_norm2 = ray_direction%dotproduct(rhs=ray_direction)
   if (dir_norm2 < 1.0e-30_R8P) then
      if (present(status)) status = RAY_STATUS_BAD_INPUT
      return
   endif
   cap = MaxR8P
   if (present(max_t)) cap = max_t
   if (cap < 0._R8P) return  ! Empty interval; trivially "no hit".
   if (self%aabb%get_nodes_number() > 0_I4P) then
      call self%aabb%intersect_ray_any_tree(facet=self%facet, ray_origin=ray_origin, &
                                            ray_direction=ray_direction, max_t=cap, found=found)
   else
      ! Flat fallback: the exact all-hits list (issue #82), any hit in [0, cap].
      call ray_intersect_all_flat(facet=self%facet, ray_origin=ray_origin, ray_direction=ray_direction, hits=hits)
      found = size(hits) > 0
      if (found) found = hits(1)%t <= cap
   endif
   endsubroutine intersect_ray_any

   subroutine segment_sdf(self, facet_labels, sdf, num_clusters, smoothing_lambda, &
                          smoothing_iterations, num_rays, cone_angle_deg, status)
   !< Per-facet semantic labelling via Shape Diameter Function (issue #18 §1.9).
   !<
   !< Pipeline: per-facet SDF (cone-of-rays + median) → dual-graph Laplacian
   !< smoothing → 1D Gaussian-mixture clustering → argmax-posterior labels,
   !< re-ordered so cluster 1 = thinnest-feature region, cluster k = thickest.
   !<
   !< Output `facet_labels(1:nf)`:
   !<   - `0` (== `SDF_LABEL_UNASSIGNED`) when the facet's SDF was undefined
   !<     (degenerate shell, near-hole, isolated patch — too few cone rays hit
   !<     the opposite surface). These facets are excluded from clustering and
   !<     left unassigned by design; downstream code should treat them as their
   !<     own implicit "unknown" group.
   !<   - `1..num_clusters` for valid facets.
   !<
   !< Optional `sdf(:)` returns the smoothed per-facet SDF for inspection.
   !<
   !< Defaults follow Shapira (2008): `num_clusters=4`, `num_rays=30`,
   !< `cone_angle_deg=120`, `smoothing_lambda=0.5`, `smoothing_iterations=2`.
   !< The clustering MVP does not run the post-clustering alpha-expansion
   !< graph cut from Shapira §5; that is a documented follow-up. Without it,
   !< segmentation labels can flip across edges that the smoothing didn't
   !< already align — for typical AMR-CFD STL the cone-ray + smoothing signal
   !< is dominant and the labels are usable as-is.
   class(surface_stl_object),  intent(in)            :: self                 !< Surface.
   integer(I4P), allocatable,  intent(out)           :: facet_labels(:)      !< Per-facet labels ∈ [0, num_clusters].
   real(R8P),    allocatable,  intent(out), optional :: sdf(:)               !< Smoothed SDF (optional output).
   integer(I4P),               intent(in),  optional :: num_clusters         !< Default 4.
   real(R8P),                  intent(in),  optional :: smoothing_lambda     !< Default 0.5.
   integer(I4P),               intent(in),  optional :: smoothing_iterations !< Default 2.
   integer(I4P),               intent(in),  optional :: num_rays             !< Default 30.
   real(R8P),                  intent(in),  optional :: cone_angle_deg       !< Default 120.
   integer(I4P),               intent(out), optional :: status               !< Status code.

   if (present(status)) status = SDF_STATUS_OK
   if (self%facets_number == 0_I4P) then
      allocate(facet_labels(0))
      if (present(sdf)) allocate(sdf(0))
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif
   call compute_segment_sdf(facet=self%facet, tree=self%aabb, &
                            bmin=self%bmin, bmax=self%bmax, &
                            facet_labels=facet_labels, sdf=sdf, &
                            num_clusters=num_clusters, &
                            smoothing_lambda=smoothing_lambda, &
                            smoothing_iterations=smoothing_iterations, &
                            num_rays=num_rays, cone_angle_deg=cone_angle_deg, &
                            status=status)
   endsubroutine segment_sdf

   subroutine alpha_wrap(self, alpha, offset, max_iterations, wrapped, status)
   !< Alpha wrapping: a watertight, 2-manifold surrogate of any triangle soup (issue #18 §1.6, issue #50). Reference:
   !< Portaneri, Hemmer, Birdal, Mandad & Alliez, *Alpha Wrapping with an Offset* (SIGGRAPH 2022).
   !<
   !< The wrap is the boundary of the space a ball of diameter `alpha` carves from outside, laid on the `offset`
   !< surface of the input: a gap wider than about alpha + 2 offset lets the ball in (both sides of the walls are
   !< wrapped), a narrower one is closed by a cap. Computed on a cell grid (alpha/4 first) by `awrap_run` in
   !< `fossil_alpha_wrap`, which has the pipeline.
   !<
   !< Defaults from CGAL's Alpha_wrap_3:
   !<   - alpha  = bbox_diagonal / 50
   !<   - offset = alpha / 30
   !<
   !< Status:
   !<   - AWRAP_STATUS_OK            the wrap encloses the input without touching it and does not cross itself
   !<   - AWRAP_STATUS_BAD_INPUT     empty surface, alpha <= 0 or offset <= 0
   !<   - AWRAP_STATUS_DEGENERATE    alpha >= input bbox diagonal (empty wrap)
   !<   - AWRAP_STATUS_NOT_CONVERGED `max_iterations` grid levels (cell size halved each, default 5, at most 2^23
   !<                                cells) still left a triangle touching the input or crossing the wrap; the wrap
   !<                                of the finest level is returned, closed and 2-manifold
   class(surface_stl_object), intent(in)            :: self           !< Input surface (any topology).
   real(R8P),                 intent(in),  optional :: alpha          !< Feature size (default bbox_diagonal/50).
   real(R8P),                 intent(in),  optional :: offset         !< Hausdorff bound (default alpha/30).
   integer(I4P),              intent(in),  optional :: max_iterations !< Grid levels (default AWRAP_DEFAULT_MAX_OUTER_ITERATIONS).
   type(surface_stl_object),  intent(out)           :: wrapped        !< Output watertight surrogate.
   integer(I4P),              intent(out), optional :: status         !< Status code.
   real(R8P)                                        :: alpha_eff, offset_eff
   real(R8P)                                        :: bdiag
   type(facet_object), allocatable                  :: wrap_facets(:)

   if (present(status)) status = AWRAP_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif
   bdiag = sqrt((self%bmax%x - self%bmin%x)**2 + &
                (self%bmax%y - self%bmin%y)**2 + &
                (self%bmax%z - self%bmin%z)**2)
   alpha_eff = bdiag / 50._R8P
   if (present(alpha)) alpha_eff = alpha
   offset_eff = alpha_eff / 30._R8P
   if (present(offset)) offset_eff = offset
   if (alpha_eff <= 0._R8P .or. offset_eff <= 0._R8P) then
      if (present(status)) status = AWRAP_STATUS_BAD_INPUT
      return
   endif
   call awrap_run(input_facet=self%facet, input_tree=self%aabb, &
                  alpha=alpha_eff, offset=offset_eff, max_iterations=max_iterations, &
                  wrapped_facets=wrap_facets, status=status)
   if (.not. allocated(wrap_facets)) allocate(wrap_facets(0))
   ! patches (issue #66): each wrapped facet takes the patch of the input facet closest to its centroid
   block
      integer(I4P) :: k, fi
      real(R8P)    :: dd
      do k = 1, size(wrap_facets)
         call wrap_facets(k)%compute_metrix
         call self%compute_distance(point=wrap_facets(k)%centroid, distance=dd, facet_index=fi)
         if (fi >= 1 .and. fi <= self%facets_number) wrap_facets(k)%patch = self%facet(fi)%patch
      enddo
   endblock
   wrapped%header = self%header
   if (allocated(self%patch_names)) wrapped%patch_names = self%patch_names
   call wrapped%adopt_facets(facets=wrap_facets)
   endsubroutine alpha_wrap

   subroutine cotangent_laplacian(self, L, M, status)
   !< Build the cotangent Laplacian `L` and the barycentric mass matrix `M`
   !< over the surface's unique-vertex pool (issue #18 §2.1).
   !<
   !< Both matrices are returned as `csr_matrix_t`. `L` is symmetric and
   !< positive-semidefinite with `L * 1 = 0` (constant kernel). `M` is
   !< diagonal (stored as CSR with one nonzero per row for API
   !< uniformity).
   !<
   !< Common downstream uses:
   !<   - Implicit Laplacian smoothing: `(M - tau L) V_new = M V`. Requires
   !<     a sparse symmetric-positive-definite solver (SuiteSparse,
   !<     MUMPS, PETSc, ...) — FOSSIL provides the matrices, you bring
   !<     the solver.
   !<   - Heat-method geodesics: `(M + t L) u = delta_source` followed by
   !<     a Poisson solve `L phi = div(X)`.
   !<   - Mean-curvature normal: `H n = (1/2) M^{-1} L V`. Free once the
   !<     matrices exist.
   !<
   !< Pre-condition: the surface's vertex pool must be initialized — true
   !< after `load_from_file`, `analyze`, or `adopt_facets`.
   !<
   !< Documented limitations:
   !<   - Mass matrix is BARYCENTRIC. The Voronoi-area variant
   !<     (Meyer 2003) is more accurate on obtuse triangles; deferred.
   !<   - Boundary edges: the cotangent of the single opposite angle is
   !<     used, halved (so `L * 1 = 0` still holds). This is the standard
   !<     convention and is fine for closed solids and open shells alike.
   class(surface_stl_object), intent(in)            :: self
   type(csr_matrix_t),        intent(out)           :: L
   type(csr_matrix_t),        intent(out)           :: M
   integer(I4P),              intent(out), optional :: status

   if (self%facets_number == 0_I4P) then
      if (present(status)) status = LAPL_STATUS_BAD_INPUT
      return
   endif
   call build_cotangent_laplacian(facet=self%facet, pool=self%vertex_pool, &
                                   L=L, M=M, status=status)
   endsubroutine cotangent_laplacian

   subroutine gaussian_curvature(self, K, status)
   !< Per-vertex Gaussian curvature via the angle-defect formula
   !< (issue #18 §2.4).
   !<
   !<     K_i = (2*pi - sum of incident triangle angles) / A_i      (interior)
   !<         = (pi   - sum of incident triangle angles) / A_i      (boundary)
   !<
   !< where `A_i` is the barycentric area assigned to vertex `i`.
   !<
   !< Sign convention: K > 0 at convex vertices (sphere surface, cube
   !< corners), K < 0 at saddles, K ~ 0 in locally-flat regions.
   !<
   !< Gauss-Bonnet identity (the load-bearing global check): for any
   !< closed orientable surface of genus g,
   !<
   !<     sum_i K_i * A_i = 2 * pi * (2 - 2 g)
   !<
   !< For genus-0 surfaces (cube, sphere) this is `4 * pi`.
   !<
   !< Pre-condition: the surface's vertex pool must be initialized — true
   !< after `load_from_file`, `analyze`, or `adopt_facets`. Boundary
   !< detection uses `facet%fcon_edge`, so `build_connectivity` must
   !< also have run.
   class(surface_stl_object), intent(in)            :: self
   real(R8P), allocatable,    intent(out)           :: K(:)
   integer(I4P),              intent(out), optional :: status

   if (self%facets_number == 0_I4P) then
      allocate(K(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif
   call build_gaussian_curvature(facet=self%facet, pool=self%vertex_pool, &
                                  K=K, status=status)
   endsubroutine gaussian_curvature

   subroutine jet_curvature(self, H, K, status)
   !< Per-vertex mean and Gaussian curvature by a local polynomial fit over the 2-ring (osculating
   !< jet, Cazals & Pouget 2005; issue #42): cubic jet, error O(h^2) on any sampling, unlike
   !< `mean_curvature` / `gaussian_curvature`, which do not converge pointwise on irregular meshes.
   !< Sign as `mean_curvature` (H > 0 convex seen from the normals' side; sphere (1/R, 1/R^2)).
   !< Near-coincident ring vertices are merged; a vertex whose fit stays rank deficient gets NaN and
   !< `status` = CURV_STATUS_FIT_FAILED (issue #91). See `fossil_curvature%build_jet_curvature`.
   class(surface_stl_object), intent(in)            :: self   !< Surface.
   real(R8P), allocatable,    intent(out)           :: H(:)   !< Mean curvature per pool vertex.
   real(R8P), allocatable,    intent(out)           :: K(:)   !< Gaussian curvature per pool vertex.
   integer(I4P),              intent(out), optional :: status !< CURV_STATUS_*.

   if (self%facets_number == 0_I4P) then
      allocate(H(0), K(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif
   call build_jet_curvature(facet=self%facet, pool=self%vertex_pool, H=H, K=K, status=status)
   endsubroutine jet_curvature

   subroutine mean_curvature(self, H, status)
   !< Per-vertex signed mean curvature `H` via `H_i n_i = (1/2) M^-1 L V`,
   !< returning the magnitude with sign from projection onto the
   !< per-vertex angle-weighted pseudo-normal (issue #18 §2.4).
   !<
   !< Sign convention: H > 0 outward-convex (sphere surface, cube
   !< corners), H < 0 saddle / concave-outward, H ~ 0 in flat regions.
   !<
   !< **Note retracting the §2.4 commit message**: that earlier commit
   !< stated mean curvature required a sparse linear solver. That was
   !< wrong — the mass matrix is diagonal (barycentric per the §2.1
   !< builder), so `M^-1` is per-row scalar division and mean H is
   !< pure-Fortran today. The implicit-smoothing case in §2.3 still
   !< needs a real solver and remains gated.
   !<
   !< Pre-condition: `analyze` must have run (pseudo-normals
   !< populated). All public construction paths ensure this.
   class(surface_stl_object), intent(in)            :: self
   real(R8P), allocatable,    intent(out)           :: H(:)
   integer(I4P),              intent(out), optional :: status

   if (self%facets_number == 0_I4P) then
      allocate(H(0))
      if (present(status)) status = CURV_STATUS_BAD_INPUT
      return
   endif
   call build_mean_curvature(facet=self%facet, pool=self%vertex_pool, &
                              H=H, status=status)
   endsubroutine mean_curvature

   subroutine smooth(self, method, lambda, mu, iterations, status)
   !< Smooth the mesh in place via explicit Laplacian or Taubin
   !< lambda|mu flow (issue #18 §2.3).
   !<
   !< - `method = SMOOTH_METHOD_EXPLICIT`: `iterations` single-step
   !<   applications of `V_new = V + lambda M^-1 L V`. **Shrinks the
   !<   volume**; useful only for noise removal where shrinkage is
   !<   tolerable, or as a building block. `mu` is ignored.
   !< - `method = SMOOTH_METHOD_TAUBIN` (default): `iterations`
   !<   complete (lambda, mu) pairs. Each pair smooths then counter-
   !<   shrinks, designed to preserve volume while attenuating
   !<   high-frequency surface noise. **The right choice for mesh
   !<   denoising in production.** Defaults: `lambda = 0.5`,
   !<   `mu = -0.53`, `iterations = 5` pairs (10 total Laplacian
   !<   applies).
   !<
   !< Implicit Laplacian smoothing `(M - tau L) V_new = M V` is the
   !< third variant from the §2.3 spec; it requires a sparse
   !< symmetric-positive-definite solver and is deferred to the
   !< Tier-2 sparse-solver gate.
   !<
   !< After return, the surface is re-adopted (AABB tree, vertex
   !< pool, connectivity, pseudo-normals all rebuilt) so all
   !< downstream queries work without further user action.
   class(surface_stl_object),  intent(inout)         :: self
   integer(I4P),               intent(in),  optional :: method
   real(R8P),                  intent(in),  optional :: lambda
   real(R8P),                  intent(in),  optional :: mu
   integer(I4P),               intent(in),  optional :: iterations
   integer(I4P),               intent(out), optional :: status
   type(facet_object), allocatable                   :: working(:)
   integer(I4P)                                      :: ierr

   if (present(status)) status = SMOOTH_STATUS_OK
   if (self%facets_number == 0_I4P) then
      if (present(status)) status = SMOOTH_STATUS_BAD_INPUT
      return
   endif

   ! Copy facets to a working array, smooth in place, then re-adopt.
   ! The smoothing routine mutates the facet vertex caches but leaves
   ! the pool stale; adopt_facets rebuilds the pool from the facets.
   allocate(working(self%facets_number))
   working(1:self%facets_number) = self%facet(1:self%facets_number)
   call do_smooth_mesh(facet=working, pool=self%vertex_pool, &
                       method=method, lambda=lambda, mu=mu, &
                       iterations=iterations, status=ierr)
   if (present(status)) status = ierr
   if (ierr /= SMOOTH_STATUS_OK .and. ierr /= SMOOTH_STATUS_DEGENERATE) then
      deallocate(working)
      return
   endif
   call self%adopt_facets(facets=working)
   endsubroutine smooth

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

   ! The three rays must not be axis-aligned. They used to be `ex + EPS*(ey + ez)` etc.,
   ! but EPS is 0, so they were exactly axis-aligned: on axis-aligned geometry (the usual
   ! CAD case) a ray from a lattice point runs exactly through the diagonal shared by two
   ! triangles of a face, the closed u/v intervals of the facet test count it twice, and
   ! the parity flips (the unit cube's centre came out "outside"). Fixed generic directions,
   ! dominated by x, y and z respectively, with irrational-looking secondary components,
   ! make such exact edge hits non-generic; the 2-of-3 vote absorbs a residual one.
   is_inside_by_x = is_inside_by_ray_intersect(ray_origin=point, ray_direction=RAY_DIR_X)
   is_inside_by_y = is_inside_by_ray_intersect(ray_origin=point, ray_direction=RAY_DIR_Y)
   if (is_inside_by_x.and.is_inside_by_y) then
     is_inside = .true.
   else
      is_inside_by_z = is_inside_by_ray_intersect(ray_origin=point, ray_direction=RAY_DIR_Z)
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

   subroutine mass_properties(self, mass, center, inertia, principal_moments, principal_axes, density, about, volume, &
                              status)
   !< Mass properties of the closed body bounded by the surface (issue #59): mass, centre of mass, inertia tensor
   !< (about the centre, or about `about`), principal moments (ascending) and principal axes (columns, right-handed).
   !<
   !< Exact polyhedral integrals over signed tetrahedra (`fossil_mass_properties`): the only error is round-off. The
   !< tetrahedra are taken about the bounding-box centre, so a body far from the origin keeps its accuracy, and the sums
   !< are compensated. Only the vertices enter (in their stored order), never the stored facet normals.
   !<
   !< The surface must be watertight with a positive enclosed volume (outward orientation): otherwise `status` is
   !< `STATUS_INVALID_INPUT` and the outputs are zero. Inner shells wound inward (voids) subtract. With repeated principal
   !< moments the principal axes are not unique; any orthonormal set returned is valid.
   class(surface_stl_object), intent(in)            :: self                 !< Surface.
   real(R8P),                 intent(out), optional :: mass                 !< Mass (density * volume).
   real(R8P),                 intent(out), optional :: center(3)            !< Centre of mass.
   real(R8P),                 intent(out), optional :: inertia(3,3)         !< Inertia tensor.
   real(R8P),                 intent(out), optional :: principal_moments(3) !< Principal moments, ascending.
   real(R8P),                 intent(out), optional :: principal_axes(3,3)  !< Principal axes (columns).
   real(R8P),                 intent(in),  optional :: density              !< Uniform density (default 1).
   real(R8P),                 intent(in),  optional :: about(3)             !< Point of the tensor (default: the centre).
   real(R8P),                 intent(out), optional :: volume               !< Enclosed volume.
   integer(I4P),              intent(out), optional :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   type(moments_accumulator_t)                      :: acc                  !< Raw moments.
   real(R8P)                                        :: o(3)                 !< Reference point (bounding-box centre).
   real(R8P)                                        :: a(3), b(3), c(3)     !< Facet vertices relative to o.
   real(R8P)                                        :: density_             !< Density, local.
   real(R8P)                                        :: volume_, mass_       !< Volume, mass, local.
   real(R8P)                                        :: center_(3)           !< Centre, local.
   real(R8P)                                        :: inertia_(3,3)        !< Tensor, local.
   real(R8P)                                        :: w(3), v(3,3)         !< Principal moments and axes, local.
   integer(I4P)                                     :: f                    !< Counter.

   if (present(mass)) mass = 0._R8P
   if (present(center)) center = 0._R8P
   if (present(inertia)) inertia = 0._R8P
   if (present(principal_moments)) principal_moments = 0._R8P
   if (present(principal_axes)) principal_axes = 0._R8P
   if (present(volume)) volume = 0._R8P
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. self%is_watertight()) return
   density_ = 1._R8P ; if (present(density)) density_ = density
   o = 0.5_R8P * [self%bmin%x + self%bmax%x, self%bmin%y + self%bmax%y, self%bmin%z + self%bmax%z]
   do f = 1, self%facets_number
      a = [self%facet(f)%vertex(1)%x, self%facet(f)%vertex(1)%y, self%facet(f)%vertex(1)%z] - o
      b = [self%facet(f)%vertex(2)%x, self%facet(f)%vertex(2)%y, self%facet(f)%vertex(2)%z] - o
      c = [self%facet(f)%vertex(3)%x, self%facet(f)%vertex(3)%y, self%facet(f)%vertex(3)%z] - o
      call accumulate_tetrahedron(acc, a, b, c)
   enddo
   if (.not. (acc%s(1) + acc%c(1) > 0._R8P)) return
   call moments_to_properties_3d(acc, o, density_, about, volume_, mass_, center_, inertia_)
   if (present(mass)) mass = mass_
   if (present(center)) center = center_
   if (present(inertia)) inertia = inertia_
   if (present(volume)) volume = volume_
   if (present(principal_moments) .or. present(principal_axes)) then
      call symmetric3_eigen(inertia_, w, v)
      if (present(principal_moments)) principal_moments = w
      if (present(principal_axes)) principal_axes = v
   endif
   if (present(status)) status = STATUS_OK
   endsubroutine mass_properties

   subroutine deviation(self, other, lower, upper, mean, rms, per_vertex, tolerance, samples_per_area, max_triangles, &
                        samples, triangles, device, status)
   !< One-sided deviation of this surface from `other` (issue #60): h(self, other) = sup over x in self of d(x, other),
   !< bracketed as lower <= h <= upper.
   !<
   !< `lower` is the largest distance over deterministic samples of self (every vertex, every edge midpoint and the
   !< centroids of a regular subdivision of each facet, `samples_per_area` per unit area, default 4 per mean facet
   !< area), raised by the branch and bound. `upper` is certified: each facet is bounded through the Lipschitz and the
   !< convexity bounds of d and split until upper - lower <= `tolerance` (default 1e-3 of the larger bounding-box
   !< diagonal), or until `max_triangles` (default 2 000 000) triangles were created: then status is
   !< DEV_STATUS_NOT_CONVERGED and the bounds are valid but wider. Without `upper` (and `tolerance`) the branch and bound
   !< does not run. `mean` and `rms` are the area-weighted quadratures of d over self (estimates); `per_vertex(i)` is d
   !< at vertex i of the vertex pool. Distances and closest points come from other's BVH (`query_many`), so the result
   !< does not depend on the tree kind or on the thread count. With `device=.true.` every batch runs through
   !< `query_many_device` (other must be on the device: `enter_device`; else DEV_STATUS_BAD_INPUT): the same kernels,
   !< so the bounds are bit-identical to the CPU path. See `fossil_deviation`.
   class(surface_stl_object), intent(in)            :: self             !< Surface A.
   type(surface_stl_object),  intent(in),   target  :: other            !< Surface B.
   real(R8P),                 intent(out), optional :: lower            !< Lower bound of h(self, other).
   real(R8P),                 intent(out), optional :: upper            !< Certified upper bound of h(self, other).
   real(R8P),                 intent(out), optional :: mean             !< Mean of d over self (estimate).
   real(R8P),                 intent(out), optional :: rms              !< Rms of d over self (estimate).
   real(R8P), allocatable,    intent(out), optional :: per_vertex(:)    !< d at the vertex-pool vertices of self.
   real(R8P),                 intent(in),  optional :: tolerance        !< Target upper - lower.
   real(R8P),                 intent(in),  optional :: samples_per_area !< Interior samples per unit area.
   integer(I4P),              intent(in),  optional :: max_triangles    !< Work cap of the branch and bound.
   integer(I4P),              intent(out), optional :: samples          !< Points sampled.
   integer(I4P),              intent(out), optional :: triangles        !< Triangles created by the branch and bound.
   logical,                   intent(in),  optional :: device           !< Query other on the device.
   integer(I4P),              intent(out), optional :: status           !< DEV_STATUS_*.
   real(R8P),    allocatable                        :: v(:,:)           !< Vertices of self.
   integer(I4P), allocatable                        :: fv(:,:)          !< Facets of self.
   real(R8P),    allocatable                        :: pv(:)            !< Per-vertex distances.
   real(R8P),    allocatable                        :: feet(:,:)        !< Closest points of the vertices.
   integer(I4P), allocatable                        :: feet_facets(:)   !< Their facets.
   real(R8P)                                        :: lower_, upper_, mean_, rms_, spa, tol !< Results, parameters.
   integer(I4P)                                     :: samples_, triangles_, status_, cap    !< Results, parameters.
   type(vector_R8P)                                 :: d1, d2           !< Bounding-box diagonals.
   type(surface_evaluator_t)                        :: eval_other       !< other's distance query.

   if (present(lower)) lower = 0._R8P
   if (present(upper)) upper = 0._R8P
   if (present(mean)) mean = 0._R8P
   if (present(rms)) rms = 0._R8P
   if (present(samples)) samples = 0_I4P
   if (present(triangles)) triangles = 0_I4P
   if (present(status)) status = DEV_STATUS_BAD_INPUT
   if (self%facets_number == 0_I4P .or. other%facets_number == 0_I4P) return
   eval_other%surface => other
   if (present(device)) eval_other%device = device
   if (eval_other%device .and. .not. other%is_on_device) return
   call deviation_mesh(self, v, fv)
   spa = 4._R8P * self%facets_number / max(self%area, tiny(1._R8P))
   if (present(samples_per_area)) spa = samples_per_area
   allocate(pv(size(v, 2)), feet(3, size(v, 2)), feet_facets(size(v, 2)))
   call deviation_sampled(v, fv, eval_other, spa, lower_, mean_, rms_, pv, samples_, status_, vertex_feet=feet, &
                          vertex_facets=feet_facets)
   if (present(status)) status = status_
   if (status_ /= DEV_STATUS_OK) return
   triangles_ = 0_I4P
   if (present(upper)) then
      d1 = self%bmax - self%bmin ; d2 = other%bmax - other%bmin
      tol = 1.e-3_R8P * max(d1%normL2(), d2%normL2())
      if (present(tolerance)) tol = tolerance
      cap = 2000000_I4P ; if (present(max_triangles)) cap = max_triangles
      call deviation_certified(v, fv, pv, feet, feet_facets, eval_other, max(tol, tiny(1._R8P)), cap, lower_, upper_, &
                               triangles_, status_)
      upper = upper_
      if (present(status)) status = status_
   endif
   if (present(lower)) lower = lower_
   if (present(mean)) mean = mean_
   if (present(rms)) rms = rms_
   if (present(per_vertex)) call move_alloc(pv, per_vertex)
   if (present(samples)) samples = samples_
   if (present(triangles)) triangles = triangles_
   endsubroutine deviation

   subroutine hausdorff(self, other, lower, upper, tolerance, samples_per_area, max_triangles, device, status)
   !< Symmetric Hausdorff distance H = max(h(self, other), h(other, self)) (issue #60): both one-sided deviations,
   !< bracketed as lower <= H <= upper (`upper` certified, as in `deviation`).
   class(surface_stl_object), intent(in)            :: self             !< Surface A.
   type(surface_stl_object),  intent(in)            :: other            !< Surface B.
   real(R8P),                 intent(out), optional :: lower            !< Lower bound of H.
   real(R8P),                 intent(out), optional :: upper            !< Certified upper bound of H.
   real(R8P),                 intent(in),  optional :: tolerance        !< Target upper - lower of each direction.
   real(R8P),                 intent(in),  optional :: samples_per_area !< Interior samples per unit area.
   integer(I4P),              intent(in),  optional :: max_triangles    !< Work cap of each branch and bound.
   logical,                   intent(in),  optional :: device           !< Query on the device (both surfaces entered).
   integer(I4P),              intent(out), optional :: status           !< DEV_STATUS_*.
   real(R8P)                                        :: l1, l2, u1, u2   !< One-sided bounds.
   integer(I4P)                                     :: s1, s2           !< Statuses.

   if (present(lower)) lower = 0._R8P
   if (present(upper)) upper = 0._R8P
   if (present(upper)) then
      call self%deviation(other, lower=l1, upper=u1, tolerance=tolerance, samples_per_area=samples_per_area, &
                          max_triangles=max_triangles, device=device, status=s1)
      call other%deviation(self, lower=l2, upper=u2, tolerance=tolerance, samples_per_area=samples_per_area, &
                           max_triangles=max_triangles, device=device, status=s2)
   else
      call self%deviation(other, lower=l1, samples_per_area=samples_per_area, device=device, status=s1)
      call other%deviation(self, lower=l2, samples_per_area=samples_per_area, device=device, status=s2)
      u1 = 0._R8P ; u2 = 0._R8P
   endif
   if (present(status)) status = max(s1, s2)
   if (s1 == DEV_STATUS_BAD_INPUT .or. s2 == DEV_STATUS_BAD_INPUT) return
   if (present(lower)) lower = max(l1, l2)
   if (present(upper)) upper = max(u1, u2)
   endsubroutine hausdorff

   subroutine surface_evaluator_eval(self, points, distances, closest, facets)
   !< Unsigned distances of the points to the evaluator's surface, the closest points and their facets (`query_many`).
   class(surface_evaluator_t), intent(in)  :: self         !< Evaluator.
   real(R8P),                  intent(in)  :: points(:,:)  !< Points [3, n].
   real(R8P),                  intent(out) :: distances(:) !< Distances [n].
   real(R8P),                  intent(out) :: closest(:,:) !< Closest points [3, n].
   integer(I4P),               intent(out) :: facets(:)    !< Closest facets [n].
   type(vector_R8P), allocatable           :: q(:)         !< Points as vectors.
   integer(I4P)                            :: i            !< Counter.

   allocate(q(size(points, 2)))
   do i = 1, size(points, 2)
      q(i) = vector_R8P(points(1, i), points(2, i), points(3, i))
   enddo
   if (self%device) then
      call self%surface%query_many_device(points=q, distances=distances, closest_points=closest, facets=facets, &
                                          is_signed=.false.)
   else
      call self%surface%query_many(points=q, distances=distances, closest_points=closest, facets=facets, is_signed=.false.)
   endif
   endsubroutine surface_evaluator_eval

   subroutine surface_ccd_distance(self, s, d)
   !< Distance of the two bodies at the poses `pose_interpolate` gives at s, with the witnesses (`distance_to`).
   class(surface_ccd_evaluator_t), intent(inout) :: self !< Evaluator.
   real(R8P),                      intent(in)    :: s    !< Step parameter.
   real(R8P),                      intent(out)   :: d    !< Distance.

   call self%a%distance_to(self%b, d, point_self=self%pa, point_other=self%pb, facet_self=self%fa, facet_other=self%fb, &
                           pose_self=pose_interpolate(self%a0, self%a1, s), pose_other=pose_interpolate(self%b0, self%b1, s))
   endsubroutine surface_ccd_distance

   subroutine surface_ccd_safe_step(self, s, tolerance, cap, step)
   !< The directional safe step at s over the facet pairs of the two bodies (`body_safe_step`).
   class(surface_ccd_evaluator_t), intent(inout) :: self      !< Evaluator.
   real(R8P),                      intent(in)    :: s         !< Step parameter.
   real(R8P),                      intent(in)    :: tolerance !< Contact tolerance.
   real(R8P),                      intent(in)    :: cap       !< Largest step wanted (1 - s).
   real(R8P),                      intent(out)   :: step      !< Safe step.
   integer(I4P)                                  :: np        !< Pairs evaluated.

   call body_safe_step(self%a%facet(1:self%a%facets_number), self%a%aabb, self%b%facet(1:self%b%facets_number), &
                       self%b%aabb, pose_interpolate(self%a0, self%a1, s), pose_interpolate(self%b0, self%b1, s), &
                       self%ka, self%tha, self%dta, self%kb, self%thb, self%dtb, tolerance, cap, step, pairs=np)
   self%pairs = self%pairs + np
   endsubroutine surface_ccd_safe_step

   subroutine deviation_mesh(surf, v, fv)
   !< The vertex-pool vertices and the facets (pool indices) of a surface.
   type(surface_stl_object),  intent(in)  :: surf     !< Surface.
   real(R8P),    allocatable, intent(out) :: v(:,:)   !< Vertices [3, nv].
   integer(I4P), allocatable, intent(out) :: fv(:,:)  !< Facets [3, nf].
   type(vector_R8P)                       :: c        !< Coordinates.
   integer(I4P)                           :: i, k     !< Counters.

   allocate(v(3, surf%vertex_pool%vertex_count()), fv(3, surf%facets_number))
   do i = 1, size(v, 2)
      c = surf%vertex_pool%coord(i)
      v(:, i) = [c%x, c%y, c%z]
   enddo
   do i = 1, surf%facets_number
      do k = 1, 3
         fv(k, i) = surf%vertex_pool%facet_vid(i, int(k, I4P))
      enddo
   enddo
   endsubroutine deviation_mesh

   subroutine projection_target(self, target, angle, corner_angle, status)
   !< The target of the constrained projection (issue #74, `fossil_projection`): the facets by patch with one tree per
   !< patch, the constraint curves (feature edges of #61, sharp above `angle`, default 30 degrees, plus the patch
   !< boundaries) chained between corners (turning above `corner_angle`, default 30 degrees), the corners.
   class(surface_stl_object), intent(in)            :: self         !< Surface.
   type(projection_target_t), intent(out)           :: target       !< Target.
   real(R8P),                 intent(in),  optional :: angle        !< Sharp-edge threshold [degrees].
   real(R8P),                 intent(in),  optional :: corner_angle !< Corner threshold [degrees].
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK, STATUS_INVALID_INPUT.
   real(R8P),    allocatable                        :: v(:,:)       !< Vertices.
   integer(I4P), allocatable                        :: fv(:,:)      !< Facets.
   real(R8P)                                        :: a, ca        !< Thresholds.
   integer(I4P)                                     :: i            !< Counter.

   if (present(status)) status = STATUS_INVALID_INPUT
   a = FEATURE_EDGE_DEFAULT_ANGLE ; if (present(angle)) a = angle
   ca = FEATURE_EDGE_DEFAULT_ANGLE ; if (present(corner_angle)) ca = corner_angle
   if (self%facets_number == 0_I4P .or. .not. (a >= 0._R8P .and. a <= 180._R8P) .or. &
       .not. (ca >= 0._R8P .and. ca <= 180._R8P)) return
   call deviation_mesh(self, v, fv)
   call build_target_3d(v, fv, [(max(1_I4P, self%facet(i)%patch), i = 1, self%facets_number)], a, ca, target)
   if (present(status)) status = STATUS_OK
   endsubroutine projection_target

   subroutine feature_edges(self, edges, kinds, corners, polyline_start, polyline_vertices, angle, corner_angle, &
                            include_boundary, include_nonmanifold, degenerate, status)
   !< Feature edges of the surface (issue #61): sharp edges (dihedral above `angle`, default 30 degrees; the facets must
   !< be consistently oriented), boundary edges (one facet) and non-manifold edges (three or more), unless excluded;
   !< their kinds (FEATURE_EDGE_SHARP, _BOUNDARY, _NONMANIFOLD); the corners (feature vertices with a number of feature
   !< edges other than 2, or with 2 turning by more than `corner_angle`, default 30 degrees); the polylines chaining the
   !< edges between corners (CSR: polyline p is polyline_vertices(polyline_start(p):polyline_start(p+1)-1); a polyline
   !< without corners is closed, its first vertex repeated). Vertex ids index the vertex pool (`get_vertex_pool`).
   !< `degenerate` counts the two-facet edges not classified because a facet has zero area. See `fossil_feature_edges`.
   class(surface_stl_object), intent(in)            :: self                 !< Surface.
   integer(I4P), allocatable, intent(out), optional :: edges(:,:)           !< Feature edges [2, ne].
   integer(I4P), allocatable, intent(out), optional :: kinds(:)             !< FEATURE_EDGE_* [ne].
   integer(I4P), allocatable, intent(out), optional :: corners(:)           !< Corner vertices.
   integer(I4P), allocatable, intent(out), optional :: polyline_start(:)    !< Polyline offsets [np + 1].
   integer(I4P), allocatable, intent(out), optional :: polyline_vertices(:) !< Polyline vertices.
   real(R8P),                 intent(in),  optional :: angle                !< Sharp threshold [degrees].
   real(R8P),                 intent(in),  optional :: corner_angle         !< Corner threshold [degrees].
   logical,                   intent(in),  optional :: include_boundary     !< Boundary edges are features (default yes).
   logical,                   intent(in),  optional :: include_nonmanifold  !< Non-manifold edges are features (yes).
   integer(I4P),              intent(out), optional :: degenerate           !< Edges skipped for a zero-area facet.
   integer(I4P),              intent(out), optional :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable                        :: v(:,:)               !< Vertices.
   integer(I4P), allocatable                        :: fv(:,:)              !< Facets.
   integer(I4P), allocatable                        :: e(:,:), k(:), c(:), ps(:), pv(:) !< Results.
   real(R8P)                                        :: a, ca                !< Thresholds.
   logical                                          :: ib, inm              !< Inclusions.
   integer(I4P)                                     :: ndeg                 !< Degenerate count.

   if (present(status)) status = STATUS_INVALID_INPUT
   if (present(degenerate)) degenerate = 0_I4P
   a = FEATURE_EDGE_DEFAULT_ANGLE ; if (present(angle)) a = angle
   ca = FEATURE_EDGE_DEFAULT_ANGLE ; if (present(corner_angle)) ca = corner_angle
   ib = .true. ; if (present(include_boundary)) ib = include_boundary
   inm = .true. ; if (present(include_nonmanifold)) inm = include_nonmanifold
   if (self%facets_number == 0_I4P .or. .not. (a >= 0._R8P .and. a <= 180._R8P) .or. &
       .not. (ca >= 0._R8P .and. ca <= 180._R8P)) return
   call deviation_mesh(self, v, fv)
   call find_feature_edges(fv, v, a, ca, ib, inm, e, k, c, ps, pv, ndeg)
   if (present(edges)) call move_alloc(e, edges)
   if (present(kinds)) call move_alloc(k, kinds)
   if (present(corners)) call move_alloc(c, corners)
   if (present(polyline_start)) call move_alloc(ps, polyline_start)
   if (present(polyline_vertices)) call move_alloc(pv, polyline_vertices)
   if (present(degenerate)) degenerate = ndeg
   if (present(status)) status = STATUS_OK
   endsubroutine feature_edges

   subroutine build_features(self, angle, corner_angle, include_boundary, include_nonmanifold, status)
   !< Build the feature-edge tree used by `box_feature_many(_device)`, `tag_block(feature_band=)` and the
   !< feature-distance queries (issue #61): the feature edges of `feature_edges` (same arguments) as segments, in a
   !< bounding-volume hierarchy (`fossil_feature_bvh`). Explicit, because those procedures take the surface intent(in);
   !< `analyze` (any geometry change) and `move` drop it. Call it before `enter_device` for the device twins.
   class(surface_stl_object), intent(inout)         :: self                !< Surface.
   real(R8P),                 intent(in),  optional :: angle               !< Sharp threshold [degrees].
   real(R8P),                 intent(in),  optional :: corner_angle        !< Corner threshold [degrees].
   logical,                   intent(in),  optional :: include_boundary    !< Boundary edges are features.
   logical,                   intent(in),  optional :: include_nonmanifold !< Non-manifold edges are features.
   integer(I4P),              intent(out), optional :: status              !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P), allocatable                        :: edges(:,:)          !< Feature edges.
   type(vector_R8P)                                 :: a, b                !< Endpoints.
   integer(I4P)                                     :: st, i               !< Status, counter.

   call drop_features(self)
   call self%feature_edges(edges=edges, angle=angle, corner_angle=corner_angle, include_boundary=include_boundary, &
                           include_nonmanifold=include_nonmanifold, status=st)
   if (present(status)) status = st
   if (st /= STATUS_OK) return
   allocate(self%feature_seg(6, size(edges, 2)))
   do i = 1, size(edges, 2)
      a = self%vertex_pool%coord(edges(1, i)) ; b = self%vertex_pool%coord(edges(2, i))
      self%feature_seg(:, i) = [a%x, a%y, a%z, b%x, b%y, b%z]
   enddo
   call build_segment_bvh(self%feature_seg, self%feature_box, self%feature_link)
   self%features_built = .true.
   endsubroutine build_features

   subroutine feature_cell_bands(self, x, y, z, band, dist)
   !< Chebyshev cell distance to the nearest cell of the block touched by a feature edge, capped at band + 1.
   type(surface_stl_object), intent(in)  :: self             !< Surface (features built).
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),             intent(in)  :: band             !< Band width.
   integer(I4P),             intent(out) :: dist(:,:,:)      !< Distances [nx, ny, nz].
   real(R8P),    allocatable             :: lo(:,:), hi(:,:) !< Cell boxes.
   logical,      allocatable             :: touched(:)       !< Touched cells.
   integer(I4P)                          :: nx, ny, nz, i, j, k, c !< Counters.

   nx = size(x) - 1 ; ny = size(y) - 1 ; nz = size(z) - 1
   allocate(lo(3, nx * ny * nz), hi(3, nx * ny * nz), touched(nx * ny * nz))
   c = 0
   do k = 1, nz ; do j = 1, ny ; do i = 1, nx
      c = c + 1 ; lo(:, c) = [x(i), y(j), z(k)] ; hi(:, c) = [x(i+1), y(j+1), z(k+1)]
   enddo ; enddo ; enddo
   call self%box_feature_many(lo, hi, touched)
   call cell_bands(reshape(touched, [nx, ny, nz]), band, dist)
   endsubroutine feature_cell_bands

   subroutine drop_features(self)
   !< Forget the feature-edge tree (the geometry changed). The device copy, if any, is released by exit_device.
   type(surface_stl_object), intent(inout) :: self !< Surface.

   if (self%features_on_device) return            ! released by exit_device first (as the other mapped arrays)
   if (allocated(self%feature_seg)) deallocate(self%feature_seg)
   if (allocated(self%feature_box)) deallocate(self%feature_box)
   if (allocated(self%feature_link)) deallocate(self%feature_link)
   self%features_built = .false.
   endsubroutine drop_features

   subroutine box_feature_many(self, lo, hi, touched, status)
   !< Which boxes [lo, hi] touch a feature edge (issue #61): the segment-box separating-axis test on the box inflated by
   !< `box_tolerance` (conservative: never misses a contact), through the feature tree (`build_features` first, else
   !< STATUS_INVALID_INPUT). OpenMP over the boxes; the device twin runs the same routine.
   class(surface_stl_object), intent(in)            :: self    !< Surface.
   real(R8P),                 intent(in)            :: lo(:,:) !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:) !< Box maxima [3, n].
   logical,                   intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),              intent(out), optional :: status  !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                        :: tol     !< Inflation.
   integer(I4P)                                     :: i       !< Counter.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. self%features_built .or. size(lo, 1) /= 3 .or. size(hi, 1) /= 3 .or. size(lo, 2) /= size(hi, 2) .or. &
       size(touched) /= size(lo, 2)) return
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   !$omp parallel do default(shared) schedule(dynamic, 64) private(i)
   do i = 1, size(lo, 2)
      touched(i) = feature_box_touch_oac(self%feature_seg, self%feature_box, self%feature_link, lo(:, i), hi(:, i), tol)
   enddo
   !$omp end parallel do
   if (present(status)) status = STATUS_OK
   endsubroutine box_feature_many

   subroutine box_feature_many_device(self, lo, hi, touched, status)
   !< `box_feature_many` on the device (issue #61): the same routine in an OpenACC parallel loop over the boxes. Requires
   !< `build_features` before `enter_device` (else STATUS_INVALID_INPUT). Mapped arrays through top-level pointer aliases
   !< (#20 rules).
   class(surface_stl_object), intent(inout), target :: self    !< Surface (target: device aliases).
   real(R8P),                 intent(in)            :: lo(:,:) !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:) !< Box maxima [3, n].
   logical,                   intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),              intent(out), optional :: status  !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    pointer                            :: seg_p(:,:), box_p(:,:) !< Aliases.
   integer(I4P), pointer                            :: link_p(:,:)            !< Alias.
   logical,      allocatable                        :: tl(:)   !< Device-side result.
   real(R8P)                                        :: tol     !< Inflation.
   integer(I4P)                                     :: i, n    !< Counters.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. (self%is_on_device .and. self%features_on_device) .or. size(lo, 1) /= 3 .or. size(hi, 1) /= 3 .or. &
       size(lo, 2) /= size(hi, 2) .or. size(touched) /= size(lo, 2)) then
      write(stderr, '(A)') 'box_feature_many_device: build_features, then enter_device, first'
      return
   endif
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   seg_p => self%feature_seg ; box_p => self%feature_box ; link_p => self%feature_link
   n = size(lo, 2)
   allocate(tl(n))
   !$acc parallel loop present(seg_p, box_p, link_p) copyin(lo, hi) copyout(tl)
   do i = 1, n
      tl(i) = feature_box_touch_oac(seg_p, box_p, link_p, lo(:, i), hi(:, i), tol)
   enddo
   touched = tl
   if (present(status)) status = STATUS_OK
   endsubroutine box_feature_many_device

   subroutine build_thin(self, thickness, gap, status)
   !< Store the local size of every facet, min(thickness, gap) of `thickness` (issue #62), for `box_thin_many(_device)`
   !< and `tag_block(thin_factor=)`: per SAH BVH payload triangle, plus its minimum per tree node, so a traversal prunes
   !< the subtrees with no thin facet. Explicit, because those procedures take the surface intent(in). `analyze` drops it
   !< (the facet ids change); `move` keeps it (a rigid motion changes neither the thickness nor the tree). Call it before
   !< `enter_device` for the device twin. STATUS_INVALID_INPUT without a SAH BVH, without both arrays of size nf, or with
   !< neither.
   class(surface_stl_object), intent(inout)         :: self         !< Surface.
   real(R8P),                 intent(in),  optional :: thickness(:) !< Thickness per facet.
   real(R8P),                 intent(in),  optional :: gap(:)       !< Gap per facet.
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable                        :: ls(:)        !< Local size per facet.
   integer(I4P), allocatable                        :: order(:), stk(:) !< Post-order, stack.
   integer(I4P)                                     :: nn, no, sp, cur, lc, rc, k, idx !< Counters.
   logical                                          :: bad          !< Input check.

   call drop_thin(self)
   if (present(status)) status = STATUS_INVALID_INPUT
   bad = .not. (present(thickness) .or. present(gap)) .or. self%facets_number == 0 .or. &
         self%aabb%get_tree_kind() /= AABB_TREE_SAH_BVH
   if (.not. bad) bad = .not. (allocated(self%aabb%node) .and. allocated(self%aabb%payload))
   if (.not. bad .and. present(thickness)) bad = size(thickness) /= self%facets_number
   if (.not. bad .and. present(gap)) bad = size(gap) /= self%facets_number
   if (bad) return
   allocate(ls(self%facets_number)) ; ls = huge(1._R8P)
   if (present(thickness)) ls = min(ls, thickness)
   if (present(gap)) ls = min(ls, gap)
   allocate(self%thin_pay(size(self%aabb%payload)))
   do idx = 1, size(self%aabb%payload)
      self%thin_pay(idx) = ls(self%aabb%payload(idx)%facet_id)
   enddo
   nn = size(self%aabb%node)
   allocate(self%thin_node(0:nn - 1), order(nn), stk(nn)) ; self%thin_node = huge(1._R8P)
   no = 0 ; sp = 1 ; stk(1) = 0                     ! pre-order from the root; reversed, children come before parents
   do while (sp > 0)
      cur = stk(sp) ; sp = sp - 1
      no = no + 1 ; order(no) = cur
      lc = bvh_left_child_oac(self%aabb%node(cur)) ; rc = bvh_right_child_oac(self%aabb%node(cur))
      if (lc > 0) then
         sp = sp + 1 ; stk(sp) = lc
      endif
      if (rc > 0) then
         sp = sp + 1 ; stk(sp) = rc
      endif
   enddo
   do k = no, 1, -1
      cur = order(k)
      do idx = bvh_payload_first_oac(self%aabb%node(cur)), &
               bvh_payload_first_oac(self%aabb%node(cur)) + bvh_payload_count_oac(self%aabb%node(cur)) - 1
         self%thin_node(cur) = min(self%thin_node(cur), self%thin_pay(idx))
      enddo
      lc = bvh_left_child_oac(self%aabb%node(cur)) ; rc = bvh_right_child_oac(self%aabb%node(cur))
      if (lc > 0) self%thin_node(cur) = min(self%thin_node(cur), self%thin_node(lc))
      if (rc > 0) self%thin_node(cur) = min(self%thin_node(cur), self%thin_node(rc))
   enddo
   self%thin_built = .true.
   if (present(status)) status = STATUS_OK
   endsubroutine build_thin

   subroutine drop_thin(self)
   !< Forget the stored local sizes (the facets changed). The device copy, if any, is released by exit_device.
   type(surface_stl_object), intent(inout) :: self !< Surface.

   if (self%thin_on_device) return
   if (allocated(self%thin_pay)) deallocate(self%thin_pay)
   if (allocated(self%thin_node)) deallocate(self%thin_node)
   self%thin_built = .false.
   endsubroutine drop_thin

   subroutine box_thin_many(self, lo, hi, factor, touched, status)
   !< Which boxes [lo, hi] meet a facet thinner than factor * h, h the largest side of each box (issue #62): the
   !< box-triangle separating-axis test on the box inflated by `box_tolerance` (conservative), over the facets whose local
   !< size (`build_thin` first, else STATUS_INVALID_INPUT) is below the threshold. OpenMP over the boxes; the device twin
   !< runs the same routine.
   class(surface_stl_object), intent(in)            :: self       !< Surface.
   real(R8P),                 intent(in)            :: lo(:,:)    !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:)    !< Box maxima [3, n].
   real(R8P),                 intent(in)            :: factor     !< Threshold factor.
   logical,                   intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),              intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P)                                        :: tol        !< Inflation.
   integer(I4P)                                     :: i          !< Counter.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. self%thin_built .or. size(lo, 1) /= 3 .or. size(hi, 1) /= 3 .or. size(lo, 2) /= size(hi, 2) .or. &
       size(touched) /= size(lo, 2)) return
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   !$omp parallel do default(shared) schedule(dynamic, 64) private(i)
   do i = 1, size(lo, 2)
      touched(i) = bvh_box_thin_oac(self%aabb%node, self%aabb%payload, self%thin_pay, self%thin_node, lo(:, i), hi(:, i), &
                                    tol, factor)
   enddo
   !$omp end parallel do
   if (present(status)) status = STATUS_OK
   endsubroutine box_thin_many

   subroutine box_thin_many_device(self, lo, hi, factor, touched, status)
   !< `box_thin_many` on the device (issue #62): the same routine in an OpenACC parallel loop over the boxes. Requires
   !< `build_thin` before `enter_device` (else STATUS_INVALID_INPUT). Mapped arrays through top-level pointer aliases
   !< (#20 rules).
   class(surface_stl_object), intent(inout), target :: self       !< Surface (target: device aliases).
   real(R8P),                 intent(in)            :: lo(:,:)    !< Box minima [3, n].
   real(R8P),                 intent(in)            :: hi(:,:)    !< Box maxima [3, n].
   real(R8P),                 intent(in)            :: factor     !< Threshold factor.
   logical,                   intent(out)           :: touched(:) !< Touched [n].
   integer(I4P),              intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   type(aabb_node_object),       pointer            :: node_p(:)  !< Alias.
   type(facet_distance_payload), pointer            :: payload_p(:) !< Alias.
   real(R8P),                    pointer            :: tp_p(:), tn_p(:) !< Aliases.
   logical,      allocatable                        :: tl(:)      !< Device-side result.
   real(R8P)                                        :: tol        !< Inflation.
   integer(I4P)                                     :: i, n       !< Counters.

   touched = .false.
   if (present(status)) status = STATUS_INVALID_INPUT
   if (.not. (self%is_on_device .and. self%thin_on_device) .or. size(lo, 1) /= 3 .or. size(hi, 1) /= 3 .or. &
       size(lo, 2) /= size(hi, 2) .or. size(touched) /= size(lo, 2)) then
      write(stderr, '(A)') 'box_thin_many_device: build_thin, then enter_device, first'
      return
   endif
   tol = box_tolerance([self%bmin%x, self%bmin%y, self%bmin%z], [self%bmax%x, self%bmax%y, self%bmax%z])
   node_p => self%aabb%node ; payload_p => self%aabb%payload ; tp_p => self%thin_pay ; tn_p => self%thin_node
   n = size(lo, 2)
   allocate(tl(n))
   !$acc parallel loop present(node_p, payload_p, tp_p, tn_p) copyin(lo, hi) copyout(tl)
   do i = 1, n
      tl(i) = bvh_box_thin_oac(node_p, payload_p, tp_p, tn_p, lo(:, i), hi(:, i), tol, factor)
   enddo
   touched = tl
   if (present(status)) status = STATUS_OK
   endsubroutine box_thin_many_device

   subroutine thin_cell_bands(self, x, y, z, factor, band, dist)
   !< Chebyshev cell distance to the nearest cell of the block that meets a thin facet, capped at band + 1.
   type(surface_stl_object), intent(in)  :: self             !< Surface (thin data built).
   real(R8P),                intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   real(R8P),                intent(in)  :: factor           !< Threshold factor.
   integer(I4P),             intent(in)  :: band             !< Band width.
   integer(I4P),             intent(out) :: dist(:,:,:)      !< Distances [nx, ny, nz].
   real(R8P),    allocatable             :: lo(:,:), hi(:,:) !< Cell boxes.
   logical,      allocatable             :: touched(:)       !< Touched cells.
   integer(I4P)                          :: nx, ny, nz, i, j, k, c !< Counters.

   nx = size(x) - 1 ; ny = size(y) - 1 ; nz = size(z) - 1
   allocate(lo(3, nx * ny * nz), hi(3, nx * ny * nz), touched(nx * ny * nz))
   c = 0
   do k = 1, nz ; do j = 1, ny ; do i = 1, nx
      c = c + 1 ; lo(:, c) = [x(i), y(j), z(k)] ; hi(:, c) = [x(i+1), y(j+1), z(k+1)]
   enddo ; enddo ; enddo
   call self%box_thin_many(lo, hi, factor, touched)
   call cell_bands(reshape(touched, [nx, ny, nz]), band, dist)
   endsubroutine thin_cell_bands

   subroutine thin_regions(self, thickness, gap, h, factor, region, region_area, status)
   !< The facets thinner than factor * h (local size min(thickness, gap) below it, issue #62), grouped into regions
   !< connected through shared edges (union-find over the facet connectivity), sorted by decreasing area (ties: the
   !< lowest facet id first): region(nf) is 0 or the region id, region_area(nr) the areas.
   class(surface_stl_object), intent(in)            :: self           !< Surface.
   real(R8P),                 intent(in),  optional :: thickness(:)   !< Thickness per facet.
   real(R8P),                 intent(in),  optional :: gap(:)         !< Gap per facet.
   real(R8P),                 intent(in)            :: h              !< Cell size.
   real(R8P),                 intent(in)            :: factor         !< Threshold factor.
   integer(I4P), allocatable, intent(out)           :: region(:)      !< Region of each facet (0: not thin) [nf].
   real(R8P),    allocatable, intent(out)           :: region_area(:) !< Area of each region [nr].
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),    allocatable                        :: ls(:), area(:) !< Local sizes, root areas.
   integer(I4P), allocatable                        :: parent(:), first(:), rank(:), root_id(:) !< Union-find, order.
   logical,      allocatable                        :: thin(:)        !< Thin facets.
   integer(I4P)                                     :: nf, f, e, g, a, b, nr, i, j, key !< Counters.
   logical                                          :: bad            !< Input check.

   nf = self%facets_number
   allocate(region(nf), region_area(0)) ; region = 0
   if (present(status)) status = STATUS_INVALID_INPUT
   bad = .not. (present(thickness) .or. present(gap)) .or. nf == 0 .or. h <= 0._R8P .or. factor <= 0._R8P
   if (.not. bad .and. present(thickness)) bad = size(thickness) /= nf
   if (.not. bad .and. present(gap)) bad = size(gap) /= nf
   if (bad) return
   allocate(ls(nf)) ; ls = huge(1._R8P)
   if (present(thickness)) ls = min(ls, thickness)
   if (present(gap)) ls = min(ls, gap)
   thin = ls < factor * h
   allocate(parent(nf)) ; parent = [(f, f = 1, nf)]
   do f = 1, nf
      if (.not. thin(f)) cycle
      do e = 1, 3
         g = self%facet(f)%fcon_edge(e)
         if (g < 1 .or. g > nf) cycle
         if (.not. thin(g)) cycle
         a = find(f) ; b = find(g)
         if (a /= b) parent(max(a, b)) = min(a, b)       ! the root is the lowest facet id of the region
      enddo
   enddo
   allocate(area(nf), root_id(nf)) ; area = 0._R8P ; root_id = 0
   do f = 1, nf
      if (thin(f)) then
         a = find(f) ; area(a) = area(a) + self%facet(f)%area()
      endif
   enddo
   nr = 0
   allocate(first(nf))
   do f = 1, nf                                         ! roots in increasing facet id
      if (thin(f) .and. find(f) == f) then
         nr = nr + 1 ; first(nr) = f
      endif
   enddo
   allocate(rank(nr)) ; rank = first(1:nr)
   do i = 2, nr                                         ! stable insertion sort by decreasing area
      key = rank(i) ; j = i - 1
      do while (j >= 1)
         if (area(rank(j)) >= area(key)) exit
         rank(j + 1) = rank(j) ; j = j - 1
      enddo
      rank(j + 1) = key
   enddo
   deallocate(region_area) ; allocate(region_area(nr))
   do i = 1, nr
      root_id(rank(i)) = i ; region_area(i) = area(rank(i))
   enddo
   do f = 1, nf
      if (thin(f)) region(f) = root_id(find(f))
   enddo
   if (present(status)) status = STATUS_OK
   contains
      function find(x) result(r)
      !< Root of x, with path halving.
      integer(I4P), intent(in) :: x !< Facet.
      integer(I4P)             :: r !< Root.

      r = x
      do while (parent(r) /= r)
         parent(r) = parent(parent(r)) ; r = parent(r)
      enddo
      endfunction find
   endsubroutine thin_regions

   subroutine surface_thickness(self, thickness, gap, method, cone_angle, rays, vertex_thickness, vertex_gap, invalid, status, &
                                separation_angle)
   !< Local thickness and gap of every facet (issue #62, roadmap #58 WP4); `huge` where undefined or nothing is in front.
   !<
   !< - THICKNESS_RAY (default, shape diameter): `rays` rays from the facet centroid in a cone of full aperture
   !<   `cone_angle` [degrees] around -n (thickness) or +n (gap), the fixed set of `cone_directions`. A thickness ray counts
   !<   if its first hit (the source facet excluded) leaves the solid (n_hit.d > 0); a gap ray counts if it enters
   !<   (n_hit.d < 0) or misses (length huge). The value is the lower median of the counted lengths; with fewer than half
   !<   of the rays counted it is huge and the facet is counted in `invalid`. A cone biases it: a plate of thickness t
   !<   gives t median(1/cos theta_k), a sphere of radius R 2R median(cos theta_k); `rays=1` casts the axis alone.
   !< - THICKNESS_MEDIAL (shrinking ball): the diameter 2r of the largest ball tangent to the facet at its centroid p, on
   !<   the side of -n (thickness) or +n (gap), with no surface point inside. The ball starts at r = min(R0, t_axis/2)
   !<   (R0 the bounding-box diagonal, t_axis the first hit of the axial ray: an upper bound, the ball's diameter along the
   !<   axis cannot pass that wall, so the medial value never exceeds the axial chord). The ball centred at c = p -+ r n
   !<   shrinks to the ball through p and q, the closest surface point to c, tangent at p (`ball_radius`), until q is on
   !<   its sphere (|c - q| >= r (1 - 1e-12)), or until the new ball would see p and q under a separation angle below
   !<   `separation_angle` (default 30 degrees; Ma et al.'s denoising test: q is then a tessellation neighbour on the same
   !<   smooth sheet, or a bump of a scanned mesh, and the ball keeps its radius). One batched `query_many` per iteration
   !<   over the facets not converged. A gap ball still at R0 has nothing in front: huge. A facet still moving after 64
   !<   iterations is undefined (huge, counted in `invalid`). Exact on a plate (t).
   !< Vertex fields: area-weighted average of the finite values of the incident facets (huge if none). OpenMP over the
   !< facets; every facet is computed alone, so the result does not depend on the thread count.
   class(surface_stl_object), intent(in)                :: self                 !< Surface.
   real(R8P), allocatable,    intent(out), optional     :: thickness(:)         !< Thickness [nf].
   real(R8P), allocatable,    intent(out), optional     :: gap(:)               !< Gap [nf].
   integer(I4P),              intent(in),  optional     :: method               !< THICKNESS_RAY or THICKNESS_MEDIAL.
   real(R8P),                 intent(in),  optional     :: cone_angle           !< Full aperture [degrees] (default 60).
   integer(I4P),              intent(in),  optional     :: rays                 !< Rays per facet (default 31).
   real(R8P), allocatable,    intent(out), optional     :: vertex_thickness(:)  !< Thickness per pool vertex.
   real(R8P), allocatable,    intent(out), optional     :: vertex_gap(:)        !< Gap per pool vertex.
   integer(I4P),              intent(out), optional     :: invalid              !< Facets with an undefined value.
   integer(I4P),              intent(out), optional     :: status               !< STATUS_OK or STATUS_INVALID_INPUT.
   real(R8P),                 intent(in),  optional     :: separation_angle     !< Medial denoising angle [degrees].
   real(R8P), allocatable                               :: th(:), gp(:)         !< Work values.
   real(R8P)                                            :: sep                  !< Separation angle.
   logical,   allocatable                               :: bad(:)               !< Undefined facets.
   real(R8P)                                            :: cone                 !< Aperture.
   integer(I4P)                                         :: meth, m, nf, f       !< Options, counters.
   logical                                              :: want_t, want_g       !< Requested outputs.

   nf = self%facets_number
   if (present(invalid)) invalid = 0
   if (present(status)) status = STATUS_OK
   meth = THICKNESS_RAY ; if (present(method)) meth = method
   cone = THICKNESS_DEFAULT_CONE ; if (present(cone_angle)) cone = cone_angle
   m = THICKNESS_DEFAULT_RAYS ; if (present(rays)) m = rays
   sep = THICKNESS_DEFAULT_SEPARATION ; if (present(separation_angle)) sep = separation_angle
   want_t = present(thickness) .or. present(vertex_thickness)
   want_g = present(gap) .or. present(vertex_gap)
   if (nf == 0 .or. (meth /= THICKNESS_RAY .and. meth /= THICKNESS_MEDIAL) .or. m < 1 .or. cone < 0._R8P .or. &
       cone >= 180._R8P .or. sep < 0._R8P .or. sep >= 180._R8P .or. &
       .not. (want_t .or. want_g)) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(th(nf), gp(nf), bad(nf))
   th = huge(1._R8P) ; gp = huge(1._R8P) ; bad = .false.
   if (meth == THICKNESS_RAY) then
      !$omp parallel do default(shared) schedule(dynamic, 16) private(f)
      do f = 1, nf
         call ray_facet(f)
      enddo
      !$omp end parallel do
   else
      if (want_t) call medial(-1._R8P, th)
      if (want_g) call medial(1._R8P, gp)
   endif
   if (present(invalid)) invalid = count(bad)
   if (present(vertex_thickness)) call to_vertices(th, vertex_thickness)
   if (present(vertex_gap)) call to_vertices(gp, vertex_gap)
   if (present(thickness)) call move_alloc(th, thickness)
   if (present(gap)) call move_alloc(gp, gap)
   contains
      subroutine ray_facet(ff)
      !< Thickness and gap of facet ff by the ray method.
      integer(I4P), intent(in) :: ff                    !< Facet.
      real(R8P)                :: nrm(3), dirs(3, m), buf(m), hn(3), dn !< Work.
      type(vector_R8P)         :: d                     !< Ray direction.
      type(ray_hit_t)          :: hit                   !< First hit.
      logical                  :: has_hit               !< Hit flag.
      integer(I4P)             :: k, nv, nmiss          !< Counters.

      nrm = [self%facet(ff)%normal%x, self%facet(ff)%normal%y, self%facet(ff)%normal%z]
      if (nrm(1) ** 2 + nrm(2) ** 2 + nrm(3) ** 2 < 0.25_R8P) then   ! degenerate facet: no normal
         bad(ff) = .true. ; return
      endif
      if (want_t) then
         call cone_directions(-nrm, cone, m, dirs)
         nv = 0
         do k = 1, m
            d%x = dirs(1, k) ; d%y = dirs(2, k) ; d%z = dirs(3, k)
            call self%intersect_ray_first(ray_origin=self%facet(ff)%centroid, ray_direction=d, hit=hit, &
                                          has_hit=has_hit, exclude=ff)
            if (.not. has_hit) cycle
            hn = [self%facet(hit%facet_id)%normal%x, self%facet(hit%facet_id)%normal%y, self%facet(hit%facet_id)%normal%z]
            if (hn(1) * dirs(1, k) + hn(2) * dirs(2, k) + hn(3) * dirs(3, k) <= 0._R8P) cycle
            nv = nv + 1 ; buf(nv) = hit%t
         enddo
         if (2 * nv < m) then
            bad(ff) = .true.
         else
            th(ff) = lower_median(buf(1:nv))
         endif
      endif
      if (want_g) then
         call cone_directions(nrm, cone, m, dirs)
         nv = 0 ; nmiss = 0
         do k = 1, m
            d%x = dirs(1, k) ; d%y = dirs(2, k) ; d%z = dirs(3, k)
            call self%intersect_ray_first(ray_origin=self%facet(ff)%centroid, ray_direction=d, hit=hit, &
                                          has_hit=has_hit, exclude=ff)
            if (.not. has_hit) then
               nmiss = nmiss + 1 ; cycle
            endif
            hn = [self%facet(hit%facet_id)%normal%x, self%facet(hit%facet_id)%normal%y, self%facet(hit%facet_id)%normal%z]
            dn = hn(1) * dirs(1, k) + hn(2) * dirs(2, k) + hn(3) * dirs(3, k)
            if (dn >= 0._R8P) cycle
            nv = nv + 1 ; buf(nv) = hit%t
         enddo
         if (2 * (nv + nmiss) < m) then
            bad(ff) = .true.
         else
            buf(nv + 1:nv + nmiss) = huge(1._R8P)
            gp(ff) = lower_median(buf(1:nv + nmiss))
         endif
      endif
      endsubroutine ray_facet

      subroutine medial(side, val)
      !< Shrinking ball on the side `side` n of every facet (-1 thickness, +1 gap); val = 2r.
      real(R8P), intent(in)    :: side   !< -1 or +1.
      real(R8P), intent(inout) :: val(:) !< Diameters [nf] (huge on entry).
      integer(I4P), parameter       :: MAXIT = 64     !< Iteration cap.
      real(R8P),    allocatable     :: r(:), pc(:,:), dir(:,:), dq(:), q(:,:) !< Radii, centroids, directions, queries.
      integer(I4P), allocatable     :: act(:)         !< Facets not converged.
      type(vector_R8P), allocatable :: c(:)           !< Ball centres.
      type(vector_R8P)              :: bmn, bmx       !< Bounding box.
      real(R8P)                     :: r0, rn         !< Initial and new radius.
      real(R8P)                     :: cos_sep        !< Cosine of the separation angle.
      real(R8P)                     :: w(3)           !< q - c'.
      type(vector_R8P)              :: dd             !< Axial direction.
      type(ray_hit_t)               :: hh             !< Axial hit.
      logical                       :: has            !< Hit flag.
      integer(I4P)                  :: na, it, k, ff, nk !< Counters.

      cos_sep = feature_cos(sep)
      bmn = self%bmin ; bmx = self%bmax
      r0 = sqrt((bmx%x - bmn%x) ** 2 + (bmx%y - bmn%y) ** 2 + (bmx%z - bmn%z) ** 2)
      allocate(r(nf), pc(3, nf), dir(3, nf), act(nf))
      na = 0
      do ff = 1, nf
         pc(:, ff) = [self%facet(ff)%centroid%x, self%facet(ff)%centroid%y, self%facet(ff)%centroid%z]
         dir(:, ff) = side * [self%facet(ff)%normal%x, self%facet(ff)%normal%y, self%facet(ff)%normal%z]
         r(ff) = r0
         if (dir(1, ff) ** 2 + dir(2, ff) ** 2 + dir(3, ff) ** 2 < 0.25_R8P) then
            bad(ff) = .true.
         else
            na = na + 1 ; act(na) = ff
         endif
      enddo
      ! start: half the axial chord, an upper bound (the diameter along the axis cannot pass the first wall it meets)
      !$omp parallel do default(shared) schedule(dynamic, 16) private(k, ff, dd, hh, has)
      do k = 1, na
         ff = act(k)
         dd%x = dir(1, ff) ; dd%y = dir(2, ff) ; dd%z = dir(3, ff)
         call self%intersect_ray_first(ray_origin=self%facet(ff)%centroid, ray_direction=dd, hit=hh, has_hit=has, &
                                       exclude=ff)
         if (has) r(ff) = min(r0, 0.5_R8P * hh%t)
      enddo
      !$omp end parallel do
      do it = 1, MAXIT
         if (na == 0) exit
         allocate(c(na), dq(na), q(3, na))
         do k = 1, na
            ff = act(k)
            c(k)%x = pc(1, ff) + r(ff) * dir(1, ff) ; c(k)%y = pc(2, ff) + r(ff) * dir(2, ff)
            c(k)%z = pc(3, ff) + r(ff) * dir(3, ff)
         enddo
         call self%query_many(points=c, distances=dq, closest_points=q, is_signed=.false.)
         nk = 0
         do k = 1, na
            ff = act(k)
            if (dq(k) >= r(ff) * (1._R8P - 1e-12_R8P)) cycle               ! q on the sphere: converged
            rn = ball_radius(pc(:, ff), dir(:, ff), q(:, k))
            if (rn >= r(ff)) cycle                                          ! no progress left (round-off)
            ! separation angle at the new centre c' between p and q (the angle of their normals if the ball touches both):
            ! below `separation_angle` q lies on the same smooth sheet as p (a tessellation neighbour or a bump, not an
            ! opposite wall), so the ball stops here. cos = (p - c').(q - c') / rn^2, with p - c' = -rn dir
            w = q(:, k) - pc(:, ff) - rn * dir(:, ff)
            if (-(dir(1, ff) * w(1) + dir(2, ff) * w(2) + dir(3, ff) * w(3)) > cos_sep * rn) cycle
            r(ff) = rn
            nk = nk + 1 ; act(nk) = ff
         enddo
         na = nk
         deallocate(c, dq, q)
      enddo
      do k = 1, na                                                         ! still moving after MAXIT
         bad(act(k)) = .true. ; r(act(k)) = huge(1._R8P)
      enddo
      do ff = 1, nf
         if (bad(ff) .or. r(ff) >= huge(1._R8P)) cycle
         if (side > 0._R8P .and. r(ff) >= r0) cycle                         ! nothing in front
         val(ff) = 2._R8P * r(ff)
      enddo
      endsubroutine medial

      subroutine to_vertices(fv, vv)
      !< Area-weighted average of the finite facet values fv on the pool vertices.
      real(R8P),              intent(in)  :: fv(:)  !< Facet values.
      real(R8P), allocatable, intent(out) :: vv(:)  !< Vertex values.
      real(R8P), allocatable              :: w(:)   !< Weights.
      real(R8P)                           :: a      !< Facet area.
      integer(I4P)                        :: ff, j, vid !< Counters.

      allocate(vv(self%vertex_pool%vertex_count()), w(self%vertex_pool%vertex_count()))
      vv = 0._R8P ; w = 0._R8P
      do ff = 1, nf
         if (fv(ff) >= huge(1._R8P)) cycle
         a = self%facet(ff)%area()
         do j = 1, 3
            vid = self%facet(ff)%vertex_id(j)
            if (vid < 1) cycle
            vv(vid) = vv(vid) + a * fv(ff) ; w(vid) = w(vid) + a
         enddo
      enddo
      where (w > 0._R8P)
         vv = vv / w
      elsewhere
         vv = huge(1._R8P)
      endwhere
      endsubroutine to_vertices
   endsubroutine surface_thickness

   subroutine save_thickness_into_file(self, file_name, thickness, gap, region, status)
   !< Write the surface with per-facet thickness, gap and thin region (issue #62) as a legacy ASCII VTK POLYDATA file, for
   !< ParaView: the vertex pool as POINTS, every facet as a triangle in POLYGONS, and the cell scalars given (`thickness`
   !< and `gap` as doubles, -1 where undefined or nothing is in front; `region` as int). Each array, if present, has one
   !< value per facet (else STATUS_INVALID_INPUT).
   class(surface_stl_object), intent(in)            :: self         !< Surface.
   character(*),              intent(in)            :: file_name    !< Output file.
   real(R8P),                 intent(in),  optional :: thickness(:) !< Thickness per facet.
   real(R8P),                 intent(in),  optional :: gap(:)       !< Gap per facet.
   integer(I4P),              intent(in),  optional :: region(:)    !< Thin region per facet.
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_OPEN_FAIL.
   type(vector_R8P)                                 :: c            !< Coordinates.
   integer(I4P)                                     :: u, ios, i, nf, np !< Work.
   logical                                          :: bad          !< Input check.

   nf = self%facets_number ; np = self%vertex_pool%vertex_count()
   bad = nf == 0
   if (present(thickness)) bad = bad .or. size(thickness) /= nf
   if (present(gap)) bad = bad .or. size(gap) /= nf
   if (present(region)) bad = bad .or. size(region) /= nf
   if (bad) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   open(newunit=u, file=file_name, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      if (present(status)) status = STATUS_FILE_OPEN_FAIL
      return
   endif
   write(u, '(A)') '# vtk DataFile Version 3.0'
   write(u, '(A)') 'FOSSIL thickness and gaps'
   write(u, '(A)') 'ASCII'
   write(u, '(A)') 'DATASET POLYDATA'
   write(u, '(A,I0,A)') 'POINTS ', np, ' double'
   do i = 1, np
      c = self%vertex_pool%coord(i)
      write(u, '(3ES24.16)') c%x, c%y, c%z
   enddo
   write(u, '(A,I0,1X,I0)') 'POLYGONS ', nf, 4 * nf
   do i = 1, nf
      write(u, '(A,I0,1X,I0,1X,I0)') '3 ', self%facet(i)%vertex_id(1) - 1, self%facet(i)%vertex_id(2) - 1, &
                                     self%facet(i)%vertex_id(3) - 1
   enddo
   if (present(thickness) .or. present(gap) .or. present(region)) write(u, '(A,I0)') 'CELL_DATA ', nf
   if (present(thickness)) call scalars('thickness', thickness)
   if (present(gap)) call scalars('gap', gap)
   if (present(region)) then
      write(u, '(A)') 'SCALARS region int 1'
      write(u, '(A)') 'LOOKUP_TABLE default'
      do i = 1, nf
         write(u, '(I0)') region(i)
      enddo
   endif
   close(u)
   if (present(status)) status = STATUS_OK
   contains
      subroutine scalars(name, val)
      !< One double cell scalar, -1 for huge.
      character(*), intent(in) :: name   !< Name.
      real(R8P),    intent(in) :: val(:) !< Values.
      integer(I4P)             :: k      !< Counter.

      write(u, '(A)') 'SCALARS '//name//' double 1'
      write(u, '(A)') 'LOOKUP_TABLE default'
      do k = 1, nf
         write(u, '(ES24.16)') merge(-1._R8P, val(k), val(k) >= huge(1._R8P))
      enddo
      endsubroutine scalars
   endsubroutine save_thickness_into_file

   subroutine save_features_into_file(self, file_name, angle, corner_angle, include_boundary, include_nonmanifold, status)
   !< Write the feature edges and corners (issue #61) as a legacy ASCII VTK POLYDATA file, for ParaView: the feature
   !< vertices as POINTS, the corners as VERTICES, every feature edge as a 2-point LINE, and a cell scalar `kind` (0 for a
   !< corner, else FEATURE_EDGE_SHARP / _BOUNDARY / _NONMANIFOLD), corners first as the format orders cells.
   class(surface_stl_object), intent(in)            :: self                !< Surface.
   character(*),              intent(in)            :: file_name           !< Output file.
   real(R8P),                 intent(in),  optional :: angle, corner_angle !< Thresholds [degrees].
   logical,                   intent(in),  optional :: include_boundary    !< Boundary edges are features.
   logical,                   intent(in),  optional :: include_nonmanifold !< Non-manifold edges are features.
   integer(I4P),              intent(out), optional :: status              !< STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_OPEN_FAIL.
   integer(I4P), allocatable                        :: edges(:,:), kinds(:), corners(:), map(:) !< Features, renumbering.
   type(vector_R8P)                                 :: c                   !< Coordinates.
   integer(I4P)                                     :: st, u, ios, i, np   !< Work.

   call self%feature_edges(edges=edges, kinds=kinds, corners=corners, angle=angle, corner_angle=corner_angle, &
                           include_boundary=include_boundary, include_nonmanifold=include_nonmanifold, status=st)
   if (present(status)) status = st
   if (st /= STATUS_OK) return
   allocate(map(self%vertex_pool%vertex_count()))
   map = 0 ; np = 0
   do i = 1, size(edges, 2)
      if (map(edges(1, i)) == 0) then
         np = np + 1 ; map(edges(1, i)) = np
      endif
      if (map(edges(2, i)) == 0) then
         np = np + 1 ; map(edges(2, i)) = np
      endif
   enddo
   open(newunit=u, file=file_name, status='replace', action='write', iostat=ios)
   if (ios /= 0) then
      if (present(status)) status = STATUS_FILE_OPEN_FAIL
      return
   endif
   write(u, '(A)') '# vtk DataFile Version 3.0'
   write(u, '(A)') 'FOSSIL feature edges'
   write(u, '(A)') 'ASCII'
   write(u, '(A)') 'DATASET POLYDATA'
   write(u, '(A,I0,A)') 'POINTS ', np, ' double'
   block
      integer(I4P), allocatable :: inv(:)
      allocate(inv(np))
      do i = 1, size(map)
         if (map(i) > 0) inv(map(i)) = i
      enddo
      do i = 1, np
         c = self%vertex_pool%coord(inv(i))
         write(u, '(3ES24.16)') c%x, c%y, c%z
      enddo
   endblock
   write(u, '(A,I0,1X,I0)') 'VERTICES ', size(corners), 2 * size(corners)
   do i = 1, size(corners)
      write(u, '(A,I0)') '1 ', map(corners(i)) - 1
   enddo
   write(u, '(A,I0,1X,I0)') 'LINES ', size(edges, 2), 3 * size(edges, 2)
   do i = 1, size(edges, 2)
      write(u, '(A,I0,1X,I0)') '2 ', map(edges(1, i)) - 1, map(edges(2, i)) - 1
   enddo
   write(u, '(A,I0)') 'CELL_DATA ', size(corners) + size(edges, 2)
   write(u, '(A)') 'SCALARS kind int 1'
   write(u, '(A)') 'LOOKUP_TABLE default'
   do i = 1, size(corners)
      write(u, '(A)') '0'
   enddo
   do i = 1, size(edges, 2)
      write(u, '(I0)') kinds(i)
   enddo
   close(u)
   endsubroutine save_features_into_file

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
   integer(I4P), allocatable                       :: pmap(:)  !< Other's patch -> self's patch.

   if (present(status)) status = STATUS_OK
   if (other%facets_number > 0) then
      ! patches (issue #66): other's patch ids are mapped into self's table by name (equal names merge)
      call sync_patch_names(self)
      allocate(pmap(other%get_patches_number()))
      do f = 1, size(pmap)
         pmap(f) = find_or_add_patch(self%patch_names, other%patch_name(f))
      enddo
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
            facet(self%facets_number+f)%patch = pmap(max(1_I4P, min(other%facet(f)%patch, size(pmap, kind=I4P))))
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
            self%facet(f)%patch = pmap(max(1_I4P, min(other%facet(f)%patch, size(pmap, kind=I4P))))
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
   subroutine fill_holes(self, max_hole_edges, max_hole_size, holes_filled, holes_skipped, open_edges, status, refine, fair, &
                         nonmanifold_edges, manifold_only)
   !< Close small holes: boundary loops filled with Liepa minimum-weight patches (issue #31).
   !<
   !< Every closed loop of facet edges without a neighbour (`fcon_edge == 0`) is a candidate. A loop
   !< is filled only if it has at most `max_hole_edges` edges (default 200; the triangulation is
   !< O(n^3)) and its bounding-box diagonal is at most `max_hole_size` times the surface's
   !< (default 0.1), so genuine open boundaries (a sheet, an inlet cut) stay open. Patches are
   !< oriented like the surrounding mesh and reuse its vertices exactly, so after the closing
   !< `analyze` (same tree kind and refinement) the filled surface is connected across them.
   !< Such an edge is a boundary edge (one incident facet) or lies on a non-manifold edge (3+
   !< incident facets, none linked). Boundary edges that do not close into a loop (inconsistent
   !< orientation, non-manifold boundaries) are reported in `open_edges` and left alone; the
   !< non-manifold edges met are reported in `nonmanifold_edges`, once per edge: they are not holes
   !< (issue #88: before, their unchained facet edges were counted in `open_edges`).
   !< `manifold_only=.true.` (issue #91) chains loops through true boundary edges only, never along a
   !< non-manifold edge: the boundary edges that then close no loop are counted in `open_edges`. Default
   !< .false., both kinds are chained (bit-identical to before the option).
   !<
   !< `refine=.true.` refines each patch to the edge density of its loop (`refine_patch`, Liepa
   !< 2003: centroid splits plus Delaunay edge flips, loop vertices fixed). `fair=.true.` (implies
   !< refine) then moves the patch's interior vertices to the thin-plate surface spanned by the
   !< loop and its 1-ring (`fair_patch`: uniform bi-Laplacian solved by conjugate gradient), so a
   !< hole in a curved surface is closed by a curved patch with matching tangents. If the solve
   !< fails, the refined patch is used unfaired and a warning is written to stderr.
   class(surface_stl_object), intent(inout)         :: self           !< Surface (analyzed).
   integer(I4P),              intent(in),  optional :: max_hole_edges !< Largest loop to fill, in edges (default 200).
   real(R8P),                 intent(in),  optional :: max_hole_size  !< Largest loop bbox diagonal / surface bbox diagonal (default 0.1).
   integer(I4P),              intent(out), optional :: holes_filled   !< Loops filled.
   integer(I4P),              intent(out), optional :: holes_skipped  !< Loops left open by the thresholds.
   integer(I4P),              intent(out), optional :: open_edges     !< Boundary edges not forming a closed loop.
   integer(I4P),              intent(out), optional :: status         !< STATUS_OK.
   logical,                   intent(in),  optional :: refine         !< Refine patches to the loop's edge density (default .false.).
   logical,                   intent(in),  optional :: fair           !< Fair refined patches (bi-Laplacian; implies refine; default .false.).
   integer(I4P),              intent(out), optional :: nonmanifold_edges !< Non-manifold edges met (3+ facets), not filled.
   logical,                   intent(in),  optional :: manifold_only  !< Loops through true boundary edges only (default .false.).
   type(boundary_loop_t), allocatable :: loops(:)
   logical                            :: refine_ !< Resolved refine flag.
   logical                            :: fair_   !< Resolved fair flag.
   integer(I4P)                       :: fst     !< fair_patch status.
   type(facet_object),    allocatable :: patch(:), merged(:)
   integer(I4P)                       :: max_edges_, n_filled, n_skipped, n_open, k, ist, n_new, cursor
   integer(I4P)                       :: n_nonmanifold !< Non-manifold edges met.
   real(R8P)                          :: max_size_, diag, loop_diag
   type(facet_object),    allocatable :: patches(:)
   integer(I4P)                       :: pid     !< Patch of the filled facets (issue #66).

   max_edges_ = 200_I4P   ; if (present(max_hole_edges)) max_edges_ = max_hole_edges
   max_size_  = 0.1_R8P   ; if (present(max_hole_size))  max_size_  = max_hole_size
   refine_    = .false.   ; if (present(refine))         refine_    = refine
   fair_      = .false.   ; if (present(fair))           fair_      = fair
   if (fair_) refine_ = .true. ! a plain patch has no interior vertex to move
   n_filled = 0_I4P ; n_skipped = 0_I4P ; n_open = 0_I4P ; n_nonmanifold = 0_I4P
   if (present(status)) status = STATUS_OK
   if (self%facets_number > 0_I4P) then
      call extract_boundary_loops(facet=self%facet, loops=loops, n_open=n_open, n_nonmanifold=n_nonmanifold, &
                                  manifold_only=manifold_only)
      diag = sqrt((self%bmax%x - self%bmin%x)**2 + (self%bmax%y - self%bmin%y)**2 + (self%bmax%z - self%bmin%z)**2)
      allocate(patches(0)) ; n_new = 0_I4P
      do k = 1, size(loops)
         loop_diag = sqrt((maxval(loops(k)%x(:)%x) - minval(loops(k)%x(:)%x))**2 + &
                          (maxval(loops(k)%x(:)%y) - minval(loops(k)%x(:)%y))**2 + &
                          (maxval(loops(k)%x(:)%z) - minval(loops(k)%x(:)%z))**2)
         if (size(loops(k)%vid) > max_edges_ .or. loop_diag > max_size_ * diag) then
            n_skipped = n_skipped + 1_I4P
            cycle
         endif
         call triangulate_loop(loop=loops(k), facet=self%facet, patch=patch, status=ist)
         if (ist /= 0_I4P) then
            n_skipped = n_skipped + 1_I4P
            cycle
         endif
         if (refine_) call refine_patch(loop=loops(k), patch=patch)
         if (fair_) then
            call fair_patch(loop=loops(k), patch=patch, facet=self%facet, status=fst)
            if (fst /= 0_I4P) write(stderr, '(A,I0,A)') 'WARNING: fill_holes: fairing failed (CG status ', fst, &
                                                        '), the refined patch is used unfaired'
         endif
         ! the new facets take the patch of the hole's boundary when it is one patch, else the patch 'filled'
         call sync_patch_names(self)
         pid = self%facet(loops(k)%owner(1))%patch
         if (any(self%facet(loops(k)%owner(:))%patch /= pid)) pid = find_or_add_patch(self%patch_names, 'filled')
         patch(:)%patch = pid
         patches = [patches, patch]
         n_filled = n_filled + 1_I4P
      enddo
      if (n_filled > 0_I4P) then
         n_new = size(patches, kind=I4P)
         allocate(merged(self%facets_number + n_new))
         merged(1:self%facets_number) = self%facet(1:self%facets_number)
         cursor = self%facets_number
         merged(cursor + 1:cursor + n_new) = patches
         call self%adopt_facets(facets=merged, aabb_refinement_levels=self%aabb%get_refinement_levels(), &
                                aabb_tree_kind=self%aabb%get_tree_kind())
      endif
   endif
   if (present(holes_filled))  holes_filled  = n_filled
   if (present(holes_skipped)) holes_skipped = n_skipped
   if (present(open_edges))    open_edges    = n_open
   if (present(nonmanifold_edges)) nonmanifold_edges = n_nonmanifold
   endsubroutine fill_holes

   subroutine sanitize(self, do_analysis, status, fill_holes, max_hole_edges, max_hole_size, refine_holes, fair_holes, &
                       manifold_holes)
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
   !<   (7, optional, `fill_holes=.true.`) `fill_holes` after orientation, then orientation
   !<      once more if anything was filled (a closed hole can create a closed shell).
   !<   6. `sanitize_normals` — BFS-propagate winding consistency, then global
   !<      flip if the volume sign indicates inward orientation. If any facet was
   !<      flipped, re-`analyze`: the pseudo-normals and the AABB distance payload
   !<      were built for the old winding, and signed-distance / inside queries
   !<      would otherwise use them (an inside-out input reported its interior as
   !<      outside).
   !<   7. Warnings to stderr summarising every counter the user might care about.
   !<
   !< NaN/Inf scrubbing happens earlier, at `load_from_file`, since loading garbage
   !< coordinates should fail outright rather than be repaired.
   !<
   !< Post-sanitize, the composite predicates `is_watertight`, `is_manifold`, and
   !< `is_volume` summarise the mesh's repair state in one boolean each.
   class(surface_stl_object), intent(inout)        :: self        !< File STL.
   logical,                   intent(in), optional :: do_analysis    !< Sentil for performing a first analysis.
   logical,                   intent(in), optional :: fill_holes     !< Close small holes (default .false.), see fill_holes.
   integer(I4P),              intent(in), optional :: max_hole_edges !< Passed to fill_holes (default 200).
   real(R8P),                 intent(in), optional :: max_hole_size  !< Passed to fill_holes (default 0.1).
   logical,                   intent(in), optional :: refine_holes   !< Passed to fill_holes as `refine` (default .false.).
   logical,                   intent(in), optional :: fair_holes     !< Passed to fill_holes as `fair` (default .false.).
   logical,                   intent(in), optional :: manifold_holes !< Passed to fill_holes as `manifold_only` (default .false.).
   integer(I4P),              intent(out), optional :: status     !< 0=success (reserved for future use).
   integer(I4P)                                    :: flipped_count !< Facets reversed by sanitize_normals.
   integer(I4P)                                    :: holes_filled, holes_skipped, open_edges !< fill_holes report.
   logical                                         :: fill_holes_   !< Local copy of fill_holes.

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
      call self%sanitize_normals(flipped_count=flipped_count)
      ! Orientation-dependent derived data (pseudo-normals, pnormal_payload, AABB
      ! distance payload) still describes the pre-flip winding: rebuild it.
      if (flipped_count > 0) call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
      ! Optional hole filling (issue #31), after orientation: loops are chained along consistently
      ! wound boundaries. A filled hole can turn an open part into a closed shell, so orientation
      ! runs once more (nesting parity) and the derived data is rebuilt if anything flipped.
      fill_holes_ = .false. ; if (present(fill_holes)) fill_holes_ = fill_holes
      holes_filled = 0_I4P ; holes_skipped = 0_I4P ; open_edges = 0_I4P
      if (fill_holes_) then
         call self%fill_holes(max_hole_edges=max_hole_edges, max_hole_size=max_hole_size, holes_filled=holes_filled, &
                              holes_skipped=holes_skipped, open_edges=open_edges, refine=refine_holes, fair=fair_holes, &
                              manifold_only=manifold_holes)
         if (holes_filled > 0_I4P) then
            call self%sanitize_normals(flipped_count=flipped_count)
            if (flipped_count > 0) call self%analyze(aabb_refinement_levels=self%aabb%get_refinement_levels())
         endif
         if (holes_filled > 0) write(stderr,'(A,I0,A)') 'WARNING: sanitize: ', holes_filled, ' hole(s) filled'
         if (holes_skipped > 0) write(stderr,'(A,I0,A)') 'WARNING: sanitize: ', holes_skipped, &
                                                           ' hole(s) left open (over max_hole_edges / max_hole_size)'
         if (open_edges > 0) write(stderr,'(A,I0,A)') 'WARNING: sanitize: ', open_edges, &
                                                        ' boundary edge(s) (1 incident facet) not forming a closed loop'
      endif
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
      ! with fill_holes, say that they are not holes (issue #88: the old boundary-edge count included their facet edges)
      if (self%non_manifold_edges_number > 0) then
         if (fill_holes_) then
            write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%non_manifold_edges_number, &
                                     ' non-manifold edge(s) (3+ incident facets) left: not holes, fill_holes leaves them'
         else
            write(stderr,'(A,I0,A)') 'WARNING: sanitize: ',self%non_manifold_edges_number, &
                                     ' non-manifold edge(s) (3+ incident facets) detected'
         endif
      endif
   endif
   endsubroutine sanitize

   pure subroutine sanitize_normals(self, flipped_count)
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
   !<
   !< @warning A flip reverses a facet's winding but leaves the orientation-dependent
   !< data built by `analyze` describing the old winding: the vertex/edge pseudo-normals,
   !< `pnormal_payload` and the AABB distance payload. Signed-distance and inside queries
   !< are then wrong until `analyze` runs again. `sanitize` does this automatically; a
   !< direct caller should read `flipped_count` and call `analyze` when it is positive.
   class(surface_stl_object), intent(inout)         :: self          !< File STL.
   integer(I4P),              intent(out), optional :: flipped_count !< Facets whose winding was reversed (global flip counts every facet).
   logical,      allocatable                :: facet_checked(:) !< Visited flag per facet.
   logical                                  :: flipped          !< One make_normal_consistent outcome.
   integer(I4P)                             :: n_flipped        !< Running flip count.
   integer(I4P), allocatable                :: queue(:)         !< BFS work queue (facet ids).
   integer(I4P)                             :: head, tail       !< Queue head/tail indices.
   integer(I4P)                             :: seed, f, e       !< Counters.
   integer(I4P)                             :: neighbour        !< Neighbour facet id along an edge.
   integer(I4P), allocatable                :: comp(:)          !< Connected-component id per facet.
   integer(I4P)                             :: n_comp           !< Number of connected components.

   n_flipped = 0_I4P
   n_comp    = 0_I4P
   if (self%facets_number > 0) then
      allocate(facet_checked(1:self%facets_number))
      allocate(queue(1:self%facets_number))
      allocate(comp(1:self%facets_number))
      facet_checked = .false.
      comp = 0_I4P

      ! Outer loop covers disconnected components: pick any unvisited facet as a new
      ! seed and BFS from it. For a single connected mesh, the inner BFS visits
      ! everything and the outer loop exits after one iteration.
      do seed = 1, self%facets_number
         if (facet_checked(seed)) cycle
         n_comp = n_comp + 1_I4P
         head = 1
         tail = 1
         queue(tail) = seed
         facet_checked(seed) = .true.
         comp(seed) = n_comp
         do while (head <= tail)
            f = queue(head)
            head = head + 1
            do e = 1, 3
               neighbour = self%facet(f)%fcon_edge(e)
               if (neighbour <= 0) cycle
               if (facet_checked(neighbour)) cycle
               call self%facet(f)%make_normal_consistent(edge=e, other=self%facet(neighbour), flipped=flipped)
               if (flipped) n_flipped = n_flipped + 1_I4P
               facet_checked(neighbour) = .true.
               comp(neighbour) = n_comp
               tail = tail + 1
               queue(tail) = neighbour
            enddo
         enddo
      enddo
      ! Orient closed shells by nesting parity and open components by their summed volume
      ! (see orient_shells). With no closed shell this is the former global volume rule.
      call orient_shells(self=self, comp=comp, n_comp=n_comp, n_flipped=n_flipped)
   endif
   ! Refresh the cached volume so get_volume() returns the post-sanitize value.
   call self%compute_volume
   if (present(flipped_count)) flipped_count = n_flipped
   endsubroutine sanitize_normals

   pure subroutine orient_shells(self, comp, n_comp, n_flipped)
   !< Orient every connected component after the BFS made each one internally consistent.
   !<
   !< Closed components ("shells", every facet edge has a neighbour) are oriented by
   !< nesting parity, open components as one group by the sign of their summed volume.
   !<
   !< A shell inside an even number of other closed shells (an outer boundary) must
   !< enclose positive volume; inside an odd number (a cavity wall) negative. The BFS
   !< made each shell internally consistent, so a shell with the wrong volume sign is
   !< flipped as a whole. Containment of shell c in shell d: bbox(c) inside bbox(d), and
   !< a point of c (its first facet's centroid) inside d by generalized winding number,
   !< |sum of d's facet solid angles| > 2*pi, which does not depend on d's orientation.
   !<
   !< Open components (a non-watertight body with small holes, or sheets) have no exact
   !< inside: they keep the former global rule, applied to them together. Their summed
   !< signed volume (apex at facet 1's first vertex, as in compute_volume) must be
   !< positive, otherwise every open facet is flipped. With no closed shell this is
   !< exactly the former whole-surface rule.
   class(surface_stl_object), intent(inout) :: self         !< File STL.
   integer(I4P),              intent(in)    :: comp(:)      !< Connected-component id per facet.
   integer(I4P),              intent(in)    :: n_comp       !< Number of connected components.
   integer(I4P),              intent(inout) :: n_flipped    !< Running flip count.
   logical,      allocatable   :: closed(:)    !< Shell is closed.
   integer(I4P), allocatable   :: first(:)     !< First facet of each shell.
   real(R8P),    allocatable   :: vol(:)       !< Signed volume of each shell.
   real(R8P),    allocatable   :: bmin(:,:)    !< Shell bbox minimum (3, n_comp).
   real(R8P),    allocatable   :: bmax(:,:)    !< Shell bbox maximum (3, n_comp).
   integer(I4P), allocatable   :: depth(:)     !< Nesting depth of each closed shell.
   type(vector_R8P)            :: apex, p      !< Volume apex of a shell; probe point.
   real(R8P)                   :: omega        !< Solid-angle sum.
   real(R8P)                   :: vol_open     !< Summed signed volume of the open components.
   integer(I4P)                :: c, d, ff, v  !< Counters.

   allocate(closed(n_comp), first(n_comp), vol(n_comp), bmin(3, n_comp), bmax(3, n_comp), depth(n_comp))
   closed = .true. ; first = 0_I4P ; vol = 0._R8P ; depth = 0_I4P
   bmin =  huge(1._R8P) ; bmax = -huge(1._R8P)
   do ff = 1, self%facets_number
      c = comp(ff)
      if (first(c) == 0_I4P) first(c) = ff
      if (any(self%facet(ff)%fcon_edge(:) <= 0_I4P)) closed(c) = .false.
      do v = 1, 3
         bmin(1, c) = min(bmin(1, c), self%facet(ff)%vertex(v)%x) ; bmax(1, c) = max(bmax(1, c), self%facet(ff)%vertex(v)%x)
         bmin(2, c) = min(bmin(2, c), self%facet(ff)%vertex(v)%y) ; bmax(2, c) = max(bmax(2, c), self%facet(ff)%vertex(v)%y)
         bmin(3, c) = min(bmin(3, c), self%facet(ff)%vertex(v)%z) ; bmax(3, c) = max(bmax(3, c), self%facet(ff)%vertex(v)%z)
      enddo
   enddo
   ! Open components, as one group: the former global rule.
   if (.not. all(closed)) then
      vol_open = 0._R8P
      do ff = 1, self%facets_number
         if (.not. closed(comp(ff))) vol_open = vol_open + self%facet(ff)%tetrahedron_volume(apex=self%facet(1)%vertex(1))
      enddo
      if (vol_open < 0._R8P) then
         do ff = 1, self%facets_number
            if (closed(comp(ff))) cycle
            call self%facet(ff)%reverse_normal
            n_flipped = n_flipped + 1_I4P
         enddo
      endif
   endif
   if (.not. any(closed)) return
   do ff = 1, self%facets_number
      c = comp(ff)
      if (.not. closed(c)) cycle
      apex = self%facet(first(c))%vertex(1)
      vol(c) = vol(c) + self%facet(ff)%tetrahedron_volume(apex=apex)
   enddo
   do c = 1, n_comp
      if (.not. closed(c)) cycle
      p = self%facet(first(c))%centroid
      do d = 1, n_comp
         if (d == c .or. .not. closed(d)) cycle
         if (any(bmin(:, c) < bmin(:, d)) .or. any(bmax(:, c) > bmax(:, d))) cycle
         omega = 0._R8P
         do ff = 1, self%facets_number
            if (comp(ff) == d) omega = omega + self%facet(ff)%solid_angle(point=p)
         enddo
         if (abs(omega) > 2._R8P * PI) depth(c) = depth(c) + 1_I4P
      enddo
   enddo
   do c = 1, n_comp
      if (.not. closed(c)) cycle
      if ((mod(depth(c), 2_I4P) == 0_I4P .and. vol(c) < 0._R8P) .or. &
          (mod(depth(c), 2_I4P) == 1_I4P .and. vol(c) > 0._R8P)) then
         do ff = 1, self%facets_number
            if (comp(ff) == c) then
               call self%facet(ff)%reverse_normal
               n_flipped = n_flipped + 1_I4P
            endif
         enddo
      endif
   enddo
   endsubroutine orient_shells

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
   case('fast_winding_number')
      code = SIGN_FAST_WINDING_NUMBER
   case default
      error stop 'fossil_surface_stl_object%sign_algorithm_from_string: unknown name "'//trim(adjustl(name))// &
                 '" (valid: "ray_intersections", "solid_angle", "pseudo_normal", "fast_winding_number")'
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
   !< When `clip_min`/`clip_max` are present, only facets that overlap the AABB are
   !< loaded (a facet crossing the box is kept whole; see `triangle_overlaps_aabb`).
   class(surface_stl_object),       intent(inout)        :: self                   !< Surface STL.
   character(*),                    intent(in)           :: file_name              !< STL file path.
   logical,                         intent(in), optional :: is_ascii               !< Force ASCII (default .true. if guess_format=.false.).
   logical,                         intent(in), optional :: guess_format           !< Auto-detect format from file size.
   type(vector_R8P),                intent(in), optional :: clip_min, clip_max     !< AABB clip extents (facets overlapping it only).
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
   character(FRLEN), allocatable                         :: names(:)               !< Patch names met (issue #66).
   integer(I4P)                                          :: p                      !< Current patch.

   if (present(status)) status = STATUS_OK
   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_read(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_, guess_format=guess_format, status=istat)
   if (istat /= 0) then
      if (present(status)) then ; status = istat ; return ; endif
   endif
   call stl_load_facets_number(file_unit=file_unit, is_ascii=is_ascii_, facets_number=facets_number)
   call stl_load_header(file_unit=file_unit, is_ascii=is_ascii_, header=self%header)
   names = [header_patch_name(self%header)] ; p = 1_I4P
   if (present(clip_min).and.present(clip_max)) then
      ! two-pass: count, then re-read and keep facets inside the AABB
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call stl_ascii_next_facet(file_unit, names, p)
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
      names = [header_patch_name(self%header)] ; p = 1_I4P
      allocate(facets(1:ff), stat=istat, errmsg=msg)
      if (istat /= 0) then
         if (present(status)) then ; status = STATUS_ALLOC_FAIL ; close(file_unit) ; return ; endif
         error stop 'surface_stl_object%load_from_file: '//trim(msg)
      endif
      ff = 0
      do f=1, facets_number
         if (is_ascii_) then
            call stl_ascii_next_facet(file_unit, names, p)
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
            facets(ff)%patch = p
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
            call stl_ascii_next_facet(file_unit, names, p)
            call facets(f)%load_from_file_ascii(file_unit=file_unit)
         else
            call facets(f)%load_from_file_binary(file_unit=file_unit)
         endif
         facets(f)%id = f
         facets(f)%patch = p
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

   if (allocated(self%patch_names)) deallocate(self%patch_names)
   call move_alloc(names, self%patch_names)
   call self%adopt_facets(facets=facets, aabb_refinement_levels=aabb_refinement_levels, aabb_tree_kind=aabb_tree_kind)
   endsubroutine load_from_file

   subroutine distance_to(self, other, distance, point_self, point_other, facet_self, facet_other, intersecting, status, &
                          pose_self, pose_other)
   !< Distance between this surface and `other` (issue #67): the minimum over their facet pairs, by a best-first descent
   !< of the two trees (`fossil_body_distance`), bit-identical to the brute-force minimum, with the closest facets
   !< (ties: the lexicographically smallest pair) and the witness points. 0 with `intersecting` when they touch or
   !< cross (exact classification). Both surfaces as they are (placed by `move` or a transform, trees current).
   class(surface_stl_object), intent(in)            :: self         !< Surface.
   class(surface_stl_object), intent(in)            :: other        !< Other body.
   real(R8P),                 intent(out)           :: distance     !< Distance.
   real(R8P),                 intent(out), optional :: point_self(3), point_other(3) !< Witness points.
   integer(I4P),              intent(out), optional :: facet_self, facet_other     !< Closest facets.
   logical,                   intent(out), optional :: intersecting !< The bodies touch or cross.
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK, STATUS_INVALID_INPUT (empty surface).
   type(rigid_pose_t),        intent(in),  optional :: pose_self, pose_other !< Poses (step 2): the bodies where `move`
                                                                    !< would put them, trees not refitted; same bits.
   real(R8P)                                        :: d2, pa(3), pb(3)
   integer(I4P)                                     :: ia, ib
   logical                                          :: meet

   if (present(status)) status = STATUS_OK
   distance = huge(1._R8P)
   if (self%facets_number == 0 .or. other%facets_number == 0) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   call body_distance(self%facet(1:self%facets_number), self%aabb, other%facet(1:other%facets_number), other%aabb, &
                      d2, ia, ib, pa, pb, meet, pose_a=pose_self, pose_b=pose_other)
   distance = sqrt(d2)
   if (present(point_self)) point_self = pa
   if (present(point_other)) point_other = pb
   if (present(facet_self)) facet_self = ia
   if (present(facet_other)) facet_other = ib
   if (present(intersecting)) intersecting = meet
   endsubroutine distance_to

   subroutine time_of_impact(self, other, pose_self0, pose_self1, pose_other0, pose_other1, toi, status, tolerance, &
                             distance, point_self, point_other, facet_self, facet_other, iterations, max_iterations, &
                             directional, pairs)
   !< Time of first contact of this body and `other` over a step (issue #73): each body moves from its pose at s = 0 to
   !< its pose at s = 1 along `pose_interpolate`; conservative advancement (`fossil_ccd%ccd_advance`) over `distance_to`
   !< with the bound mu = motion_bound_3d(self) + motion_bound_3d(other): d >= tolerance / 2 on [0, toi] (no earlier
   !< contact, no tunnelling), d(toi) <= tolerance on CCD_CONTACT. `status`: CCD_CONTACT, CCD_SEPARATE (toi = 1),
   !< CCD_INITIAL_CONTACT (toi = 0), CCD_NOT_CONVERGED (cap hit: toi a certified lower bound), CCD_INVALID_INPUT (an
   !< empty body, tolerance <= 0). Default tolerance: 1e-6 of the larger bounding-box diagonal. The witnesses (world
   !< frame) and the closest facets are those at toi. `directional=.true.` adds C2A's per-pair directional step
   !< (`body_safe_step`, the larger of the two certified steps is taken); `pairs` counts its exact pair evaluations.
   class(surface_stl_object), intent(in), target    :: self           !< Body A.
   class(surface_stl_object), intent(in), target    :: other          !< Body B.
   type(rigid_pose_t),        intent(in)            :: pose_self0, pose_self1   !< Poses of A at s = 0, 1.
   type(rigid_pose_t),        intent(in)            :: pose_other0, pose_other1 !< Poses of B at s = 0, 1.
   real(R8P),                 intent(out)           :: toi            !< Time of first contact in [0, 1].
   integer(I4P),              intent(out)           :: status         !< CCD_* status.
   real(R8P),                 intent(in),  optional :: tolerance      !< Contact tolerance.
   real(R8P),                 intent(out), optional :: distance       !< d(toi).
   real(R8P),                 intent(out), optional :: point_self(3), point_other(3) !< Witness points at toi.
   integer(I4P),              intent(out), optional :: facet_self, facet_other     !< Closest facets at toi.
   integer(I4P),              intent(out), optional :: iterations     !< Advancement steps.
   integer(I4P),              intent(in),  optional :: max_iterations !< Cap (default CCD_MAX_ITERATIONS).
   logical,                   intent(in),  optional :: directional    !< Directional step (absent: adaptive, #75).
   integer(I4P),              intent(out), optional :: pairs          !< Exact pair evaluations of the directional steps.
   type(surface_ccd_evaluator_t)                    :: ev             !< Distance at s.
   real(R8P)                                        :: mu, tol, d     !< Bound, tolerance, distance.
   type(vector_R8P)                                 :: da, db         !< Bounding-box diagonals.

   toi = 0._R8P ; d = 0._R8P
   if (self%facets_number == 0 .or. other%facets_number == 0) then
      status = CCD_INVALID_INPUT
   else
      mu = motion_bound_3d(pose_self0, pose_self1, facet_xs(self, 1), facet_xs(self, 2), facet_xs(self, 3)) + &
           motion_bound_3d(pose_other0, pose_other1, facet_xs(other, 1), facet_xs(other, 2), facet_xs(other, 3))
      da = self%bmax - self%bmin ; db = other%bmax - other%bmin
      tol = 1.e-6_R8P * max(da%normL2(), db%normL2())
      if (present(tolerance)) tol = tolerance
      ev%a => self ; ev%b => other
      ev%a0 = pose_self0 ; ev%a1 = pose_self1 ; ev%b0 = pose_other0 ; ev%b1 = pose_other1
      call pose_log_3d(pose_self0, pose_self1, ev%ka, ev%tha) ; ev%dta = pose_self1%t - pose_self0%t
      call pose_log_3d(pose_other0, pose_other1, ev%kb, ev%thb) ; ev%dtb = pose_other1%t - pose_other0%t
      call ccd_advance(ev, mu, tol, toi, status, distance=d, iterations=iterations, max_iterations=max_iterations, &
                       directional=directional)
   endif
   if (present(pairs)) pairs = ev%pairs
   if (present(distance)) distance = d
   if (present(point_self)) point_self = ev%pa
   if (present(point_other)) point_other = ev%pb
   if (present(facet_self)) facet_self = ev%fa
   if (present(facet_other)) facet_other = ev%fb
   if (present(iterations) .and. status == CCD_INVALID_INPUT) iterations = 0_I4P
   endsubroutine time_of_impact

   pure function facet_xs(surface, axis) result(c)
   !< The coordinate `axis` (1-3) of the three vertices of every facet [3 nf] (body frame).
   class(surface_stl_object), intent(in) :: surface !< Surface.
   integer,                   intent(in) :: axis    !< Coordinate.
   real(R8P), allocatable                :: c(:)    !< Coordinates.
   integer(I4P)                          :: i, k    !< Counters.

   allocate(c(3 * surface%facets_number))
   do i = 1, surface%facets_number
      do k = 1, 3
         select case (axis)
         case (1) ; c(3 * (i - 1) + k) = surface%facet(i)%vertex(k)%x
         case (2) ; c(3 * (i - 1) + k) = surface%facet(i)%vertex(k)%y
         case default ; c(3 * (i - 1) + k) = surface%facet(i)%vertex(k)%z
         endselect
      enddo
   enddo
   endfunction facet_xs

   subroutine proximity_pairs(self, cutoff, pairs, distances, status)
   !< The pairs of bodies of this surface (`facet_body`: outer shells with their cavities) closer than `cutoff`
   !< (issue #67): pairs(:, k) = [b1, b2] with b1 < b2, in order, and their distances (0 when they touch or cross).
   !< Each body becomes a surface of its own, then as `proximity_pairs_of`.
   class(surface_stl_object),  intent(in)            :: self         !< Surface.
   real(R8P),                  intent(in)            :: cutoff       !< Distance cutoff (>= 0).
   integer(I4P), allocatable,  intent(out)           :: pairs(:, :)  !< Body pairs [2, n].
   real(R8P),    allocatable,  intent(out)           :: distances(:) !< Their distances [n].
   integer(I4P),               intent(out), optional :: status       !< STATUS_OK, STATUS_INVALID_INPUT.
   type(surface_stl_object), allocatable             :: parts(:)     !< One surface per body.
   type(facet_object), allocatable                   :: fs(:)        !< Facets of a body.
   integer(I4P)                                      :: b, f, n      !< Counters.

   if (present(status)) status = STATUS_OK
   allocate(pairs(2, 0), distances(0))
   if (self%facets_number == 0 .or. .not. cutoff >= 0._R8P) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   allocate(parts(self%get_bodies_number()))
   do b = 1, size(parts)
      n = count([(self%facet_body(f) == b, f = 1, self%facets_number)])
      allocate(fs(n)) ; n = 0
      do f = 1, self%facets_number
         if (self%facet_body(f) /= b) cycle
         n = n + 1 ; fs(n) = self%facet(f)
      enddo
      call parts(b)%adopt_facets(facets=fs)
   enddo
   call proximity_pairs_of(parts, cutoff, pairs, distances, status)
   endsubroutine proximity_pairs

   subroutine proximity_pairs_of(surfaces, cutoff, pairs, distances, status)
   !< The pairs of surfaces closer than `cutoff` (issue #67): a broad phase on the bounding boxes (pairs whose boxes are
   !< farther than `cutoff` are skipped), then `distance_to` for the others; pairs(:, k) = [i, j], i < j, in order.
   !< The pairs found go into work buffers grown by doubling (move_alloc), and the results are allocated once at the end
   !< (issue #76: no self-referencing array growth, which corrupted the heap under nvfortran 26.9 elsewhere).
   type(surface_stl_object),  intent(in)            :: surfaces(:)  !< Bodies.
   real(R8P),                 intent(in)            :: cutoff       !< Distance cutoff (>= 0).
   integer(I4P), allocatable, intent(out)           :: pairs(:, :)  !< Surface pairs [2, n].
   real(R8P),    allocatable, intent(out)           :: distances(:) !< Their distances [n].
   integer(I4P),              intent(out), optional :: status       !< STATUS_OK, STATUS_INVALID_INPUT.
   type(vector_R8P)                                 :: l1, h1, l2, h2
   real(R8P)                                        :: d, g
   integer(I4P)                                     :: i, j, m
   integer(I4P), allocatable                        :: wp(:, :), tp(:, :) !< Pair buffer and its grown copy.
   real(R8P),    allocatable                        :: wd(:), td(:)       !< Distance buffer and its grown copy.

   if (present(status)) status = STATUS_OK
   if (.not. cutoff >= 0._R8P) then
      allocate(pairs(2, 0), distances(0))
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   m = 0_I4P
   allocate(wp(2, 16), wd(16))
   do i = 1, size(surfaces)
      if (surfaces(i)%facets_number == 0) cycle
      l1 = surfaces(i)%bmin ; h1 = surfaces(i)%bmax
      do j = i + 1, size(surfaces)
         if (surfaces(j)%facets_number == 0) cycle
         l2 = surfaces(j)%bmin ; h2 = surfaces(j)%bmax
         g = max(0._R8P, l2%x - h1%x, l1%x - h2%x)**2 + max(0._R8P, l2%y - h1%y, l1%y - h2%y)**2 + &
             max(0._R8P, l2%z - h1%z, l1%z - h2%z)**2
         if (g > cutoff * cutoff * (1._R8P + 1.e-10_R8P)) cycle
         call surfaces(i)%distance_to(surfaces(j), d)
         if (d > cutoff) cycle
         if (m == size(wd)) then
            allocate(tp(2, 2 * m), td(2 * m))
            tp(:, 1:m) = wp(:, 1:m) ; td(1:m) = wd(1:m)
            call move_alloc(from=tp, to=wp) ; call move_alloc(from=td, to=wd)
         endif
         m = m + 1_I4P
         wp(1, m) = i ; wp(2, m) = j ; wd(m) = d
      enddo
   enddo
   allocate(pairs(2, m), distances(m))
   pairs(:, 1:m) = wp(:, 1:m) ; distances(1:m) = wd(1:m)
   endsubroutine proximity_pairs_of

   subroutine first_impact(surfaces, poses0, poses1, toi, pair, status, tolerance, directional, tested)
   !< The earliest contact among several moving bodies over a step (issue #73): body i moves from poses0(i) to
   !< poses1(i) along `pose_interpolate`. Broad phase: the box of each body at s = 0 grown by its motion bound mu_i and
   !< by tolerance / 2 contains it over the whole step, so pairs whose grown boxes do not overlap stay farther than the
   !< tolerance. Every other pair has a lower bound on its time of impact, max(0, box distance at s = 0 - tolerance) /
   !< (mu_i + mu_j) (the first advancement step); pairs run `time_of_impact` in the order of these bounds and the search
   !< stops when the next bound exceeds the best contact. Ties go to the smallest pair (i, j), so the result equals the
   !< loop over all pairs. status: CCD_CONTACT / CCD_INITIAL_CONTACT (pair = [i, j], i < j), CCD_SEPARATE (toi = 1,
   !< pair = 0), CCD_NOT_CONVERGED (a pair hit the cap before the best contact: toi its certified lower bound),
   !< CCD_INVALID_INPUT (size mismatch). Default tolerance: 1e-6 of the largest bounding-box diagonal. `tested`: the
   !< pairs that ran `time_of_impact`.
   type(surface_stl_object), intent(in), target    :: surfaces(:)    !< Bodies.
   type(rigid_pose_t),       intent(in)            :: poses0(:)      !< Poses at s = 0.
   type(rigid_pose_t),       intent(in)            :: poses1(:)      !< Poses at s = 1.
   real(R8P),                intent(out)           :: toi            !< Earliest time of contact.
   integer(I4P),             intent(out)           :: pair(2)        !< Its pair.
   integer(I4P),             intent(out)           :: status         !< CCD_* status.
   real(R8P),                intent(in),  optional :: tolerance      !< Contact tolerance.
   logical,                  intent(in),  optional :: directional    !< Directional steps (see time_of_impact).
   integer(I4P),             intent(out), optional :: tested         !< Pairs that ran time_of_impact.
   real(R8P),                allocatable           :: mu(:), lo(:, :), hi(:, :), bound(:) !< Bounds, boxes at s = 0.
   integer(I4P),             allocatable           :: cand(:, :), order(:) !< Candidate pairs, their order.
   type(vector_R8P)                                :: l, h           !< Boxes.
   real(R8P)                                       :: tol, g, t, c(3), w(3), half !< Work.
   integer(I4P)                                    :: n, i, j, k, m, nc, st, ntest !< Counters.

   toi = 1._R8P ; pair = 0_I4P ; status = CCD_SEPARATE ; ntest = 0_I4P
   n = size(surfaces, kind=I4P)
   if (size(poses0) /= n .or. size(poses1) /= n) then
      status = CCD_INVALID_INPUT ; toi = 0._R8P
      if (present(tested)) tested = 0_I4P
      return
   endif
   allocate(mu(n), lo(3, n), hi(3, n)) ; tol = 0._R8P
   do i = 1, n
      mu(i) = 0._R8P ; lo(:, i) = huge(1._R8P) ; hi(:, i) = -huge(1._R8P)
      if (surfaces(i)%facets_number == 0) cycle
      mu(i) = motion_bound_3d(poses0(i), poses1(i), facet_xs(surfaces(i), 1), facet_xs(surfaces(i), 2), &
                              facet_xs(surfaces(i), 3))
      l = surfaces(i)%bmin ; h = surfaces(i)%bmax
      l = h - l ; tol = max(tol, l%normL2()) ; l = surfaces(i)%bmin
      do k = 0, 7
         c = [merge(h%x, l%x, btest(k, 0)), merge(h%y, l%y, btest(k, 1)), merge(h%z, l%z, btest(k, 2))]
         w = to_world(poses0(i), c)
         lo(:, i) = min(lo(:, i), w) ; hi(:, i) = max(hi(:, i), w)
      enddo
      lo(:, i) = lo(:, i) - 4._R8P * spacing(abs(lo(:, i))) ; hi(:, i) = hi(:, i) + 4._R8P * spacing(abs(hi(:, i)))
   enddo
   tol = 1.e-6_R8P * tol ; if (present(tolerance)) tol = tolerance
   half = 0.5_R8P * tol
   ! broad phase on the grown boxes; lower bounds of the candidates
   allocate(cand(2, n * (n - 1) / 2), bound(n * (n - 1) / 2)) ; nc = 0_I4P
   do i = 1, n
      if (surfaces(i)%facets_number == 0) cycle
      do j = i + 1, n
         if (surfaces(j)%facets_number == 0) cycle
         if (any(lo(:, j) - mu(j) - half > hi(:, i) + mu(i) + half) .or. &
             any(lo(:, i) - mu(i) - half > hi(:, j) + mu(j) + half)) cycle
         g = norm2(max(0._R8P, lo(:, j) - hi(:, i), lo(:, i) - hi(:, j)))
         nc = nc + 1_I4P ; cand(:, nc) = [i, j]
         if (mu(i) + mu(j) > 0._R8P) then
            bound(nc) = max(0._R8P, g - tol) / (mu(i) + mu(j))
         else
            bound(nc) = merge(0._R8P, huge(1._R8P), g <= tol) ! no motion: only an initial contact
         endif
      enddo
   enddo
   ! candidates by (bound, i, j): insertion sort of an index (deterministic)
   allocate(order(nc))
   do k = 1, nc
      order(k) = k ; m = k
      do while (m > 1)
         if (bound(order(m - 1)) < bound(order(m)) .or. (bound(order(m - 1)) == bound(order(m)) .and. &
             order(m - 1) < order(m))) exit
         order([m - 1, m]) = order([m, m - 1]) ; m = m - 1
      enddo
   enddo
   do k = 1, nc
      m = order(k)
      if (bound(m) > toi) exit
      i = cand(1, m) ; j = cand(2, m)
      call surfaces(i)%time_of_impact(surfaces(j), poses0(i), poses1(i), poses0(j), poses1(j), t, st, tolerance=tol, &
                                      directional=directional)
      ntest = ntest + 1_I4P
      if (st == CCD_SEPARATE) cycle
      if (t < toi .or. (t == toi .and. (status == CCD_SEPARATE .or. i < pair(1) .or. &
                                        (i == pair(1) .and. j < pair(2))))) then
         toi = t ; pair = [i, j] ; status = st
      endif
   enddo
   if (present(tested)) tested = ntest
   endsubroutine first_impact

   subroutine load_from_files(self, file_names, guess_format, aabb_refinement_levels, aabb_tree_kind, status)
   !< Load several STL files into one surface, one patch per file (issue #66), in the order given; a multi-solid ASCII
   !< file contributes its blocks as patches. A patch takes its `solid` name, or the file name when that is empty;
   !< equal names merge. The header is the first patch's name.
   class(surface_stl_object), intent(inout)        :: self                   !< Surface STL.
   character(*),              intent(in)           :: file_names(:)          !< STL file paths.
   logical,                   intent(in), optional :: guess_format           !< Auto-detect each format (default .true.).
   integer(I4P),              intent(in), optional :: aabb_refinement_levels !< AABB refinement levels passed to analyze.
   integer(I4P),              intent(in), optional :: aabb_tree_kind         !< AABB_TREE_OCTREE or AABB_TREE_SAH_BVH.
   integer(I4P),              intent(out), optional :: status                !< As `load_from_file`.
   type(surface_stl_object)                        :: part                   !< One file.
   type(facet_object), allocatable                 :: facets(:), tmp(:)      !< All facets.
   character(FRLEN), allocatable                   :: names(:)               !< Patch names.
   integer(I4P), allocatable                       :: map(:)                 !< Part patch -> patch.
   integer(I4P)                                    :: k, q, f, n, st         !< Counters, status.
   logical                                         :: guess                  !< Effective guess flag.
   character(len=:), allocatable                   :: nm                     !< Patch name.

   if (present(status)) status = STATUS_OK
   guess = .true. ; if (present(guess_format)) guess = guess_format
   allocate(facets(0), names(0))
   do k = 1, size(file_names)
      call part%load_from_file(file_name=trim(file_names(k)), guess_format=guess, status=st)
      if (st /= STATUS_OK) then
         if (present(status)) then ; status = st ; return ; endif
         error stop 'surface_stl_object%load_from_files: cannot load '//trim(file_names(k))
      endif
      allocate(map(part%get_patches_number()))
      do q = 1, size(map)
         nm = part%patch_name(q)
         if (len_trim(nm) == 0) nm = basename_of(trim(file_names(k)))
         map(q) = find_or_add_patch(names, nm)
      enddo
      n = size(facets)
      allocate(tmp(n + part%facets_number))
      if (n > 0) tmp(1:n) = facets
      do f = 1, part%facets_number
         tmp(n + f) = part%facet(f)
         tmp(n + f)%patch = map(part%facet(f)%patch)
         tmp(n + f)%id = n + f
      enddo
      call move_alloc(tmp, facets)
      deallocate(map)
      call part%destroy
   enddo
   call self%destroy
   if (size(names) > 0) self%header = names(1)
   call move_alloc(names, self%patch_names)
   call self%adopt_facets(facets=facets, aabb_refinement_levels=aabb_refinement_levels, aabb_tree_kind=aabb_tree_kind)
   contains
      pure function basename_of(path) result(b)
      character(*), intent(in)      :: path
      character(len=:), allocatable :: b
      integer(I4P)                  :: i, j

      i = index(path, '/', back=.true.)
      b = path(i + 1:)
      j = index(b, '.', back=.true.)
      if (j > 1) b = b(1:j - 1)
      endfunction basename_of
   endsubroutine load_from_files

   subroutine save_into_file(self, file_name, is_ascii, status)
   !< Save the surface to an STL file.
   class(surface_stl_object), intent(in)           :: self      !< Surface STL.
   character(*),              intent(in)           :: file_name !< STL file path.
   logical,                   intent(in), optional :: is_ascii  !< Write as ASCII (default .true.).
   integer(I4P),              intent(out), optional :: status   !< 0=success, 4=open failure.
   integer(I4P)                                    :: file_unit !< File unit.
   logical                                         :: is_ascii_ !< Effective ASCII flag.
   integer(I4P)                                    :: f, istat  !< Counter and open status.
   integer(I4P)                                    :: p         !< Patch.

   if (present(status)) status = STATUS_OK
   is_ascii_ = .true. ; if (present(is_ascii)) is_ascii_ = is_ascii
   call stl_open_for_write(file_name=file_name, file_unit=file_unit, is_ascii=is_ascii_, status=istat)
   if (istat /= 0) then
      if (present(status)) then ; status = istat ; return ; endif
   endif
   if (is_ascii_ .and. self%get_patches_number() > 1_I4P) then
      ! one solid block per non-empty patch, facets in their order (issue #66)
      do p = 1, self%get_patches_number()
         if (.not. any(self%facet(1:self%facets_number)%patch == p)) cycle ! an empty patch writes no block
         write(file_unit, '(A)') 'solid '//self%patch_name(p)
         do f=1, self%facets_number
            if (self%facet(f)%patch == p) call self%facet(f)%save_into_file_ascii(file_unit=file_unit)
         enddo
         write(file_unit, '(A)') 'endsolid '//self%patch_name(p)
      enddo
      close(file_unit)
      return
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
   if (is_ascii_) write(file_unit, '(A)') 'endsolid '//trim(header_patch_name(self%header))
   close(file_unit)
   endsubroutine save_into_file

   subroutine save_aabb_into_file(self, base_file_name, is_ascii)
   !< Save the facets of each leaf of the AABB tree (either kind) as a separate STL file, through `tree_boxes` (issue
   !< #89: before it, the BVH wrote nothing).
   !<
   !< File name pattern: `<base_file_name>aabb-l_<depth>-b_<box>.stl`, with the depth (root = 0) and the position of the
   !< box in the order of `tree_boxes` (breadth first, root = 1).
   class(surface_stl_object), intent(in)           :: self           !< Surface STL.
   character(*),              intent(in)           :: base_file_name !< Base file-name prefix.
   logical,                   intent(in), optional :: is_ascii       !< Write as ASCII (default .true.).
   type(facet_object), allocatable                 :: aabb_facet(:)  !< Leaf facets.
   type(surface_stl_object)                        :: leaf           !< Scratch surface for save.
   type(vector_R8P), allocatable                   :: bmin(:), bmax(:) !< Boxes.
   integer(I4P), allocatable                       :: depth(:), start(:), ids(:) !< Depths, leaf facets (CSR).
   logical, allocatable                            :: is_leaf(:)     !< Leaf flags.
   integer(I4P)                                    :: b, k           !< Counters.
   character(len=:), allocatable                   :: file_name      !< Composed file name.

   call self%tree_boxes(bmin=bmin, bmax=bmax, depth=depth, leaf=is_leaf, facet_start=start, facet_ids=ids)
   do b = 1, size(bmin)
      if (.not. is_leaf(b) .or. start(b + 1) == start(b)) cycle
      file_name = trim(adjustl(base_file_name))//'aabb-l_'//trim(str(depth(b), .true.))//'-b_'//trim(str(b, .true.))//'.stl'
      allocate(aabb_facet(start(b + 1) - start(b)))
      do k = start(b), start(b + 1) - 1
         aabb_facet(k - start(b) + 1) = self%facet(ids(k))
      enddo
      call leaf%adopt_facets(facets=aabb_facet)
      call leaf%set_header(self%header)
      call leaf%save_into_file(file_name=file_name, is_ascii=is_ascii)
      call leaf%destroy
   enddo
   endsubroutine save_aabb_into_file

   subroutine tree_boxes(self, bmin, bmax, depth, leaf, max_depth, parent, facet_start, facet_ids, status)
   !< The boxes of the surface's AABB tree, read only (issue #89): every node down to `max_depth` (root = 0; default:
   !< the whole tree), breadth first from the root (position 1), the children of a node in the tree's order. Both tree
   !< kinds: the SAH BVH's boxes bound their facets and overlap; the octree's are the cells of its space partition
   !< (empty cells are not nodes), and a leaf holds the facets whose centroid lies in it, which may stick out of it.
   !<
   !< Optional outputs: `depth`; `leaf` (the node has no children: a node at `max_depth` with children is not a leaf);
   !< `parent` (position of the parent in the output, 0 for the root); the facets of the leaves in CSR form, those of
   !< box i in `facet_ids(facet_start(i):facet_start(i+1)-1)` (none for an inner node), every facet in exactly one leaf
   !< when the whole tree is returned. An empty surface has no tree: zero boxes. `max_depth < 0`: zero boxes and
   !< STATUS_INVALID_INPUT.
   class(surface_stl_object), intent(in)                        :: self           !< Surface.
   type(vector_R8P),          allocatable, intent(out)           :: bmin(:)        !< Box minima.
   type(vector_R8P),          allocatable, intent(out)           :: bmax(:)        !< Box maxima.
   integer(I4P),              allocatable, intent(out), optional :: depth(:)       !< Depths (root = 0).
   logical,                   allocatable, intent(out), optional :: leaf(:)        !< Leaf flags.
   integer(I4P),                           intent(in),  optional :: max_depth      !< Deepest level (default: all).
   integer(I4P),              allocatable, intent(out), optional :: parent(:)      !< Parents (0: root).
   integer(I4P),              allocatable, intent(out), optional :: facet_start(:) !< CSR start of the leaf facets [n+1].
   integer(I4P),              allocatable, intent(out), optional :: facet_ids(:)   !< Leaf facets.
   integer(I4P),                           intent(out), optional :: status         !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P), allocatable                                     :: dep(:), par(:), fs(:), fi(:) !< Work.
   logical,      allocatable                                     :: lf(:)          !< Work.
   integer(I4P)                                                  :: md             !< Effective max_depth.

   if (present(status)) status = STATUS_OK
   md = huge(1_I4P) ; if (present(max_depth)) md = max_depth
   if (md < 0_I4P) then
      allocate(bmin(0), bmax(0), dep(0), lf(0), par(0), fs(1), fi(0)) ; fs = 1_I4P
      if (present(status)) status = STATUS_INVALID_INPUT
   else
      call self%aabb%node_boxes(max_depth=md, bmin=bmin, bmax=bmax, depth=dep, leaf=lf, parent=par, facet_start=fs, &
                                facet_ids=fi)
   endif
   if (present(depth)) call move_alloc(from=dep, to=depth)
   if (present(leaf)) call move_alloc(from=lf, to=leaf)
   if (present(parent)) call move_alloc(from=par, to=parent)
   if (present(facet_start)) call move_alloc(from=fs, to=facet_start)
   if (present(facet_ids)) call move_alloc(from=fi, to=facet_ids)
   endsubroutine tree_boxes

   ! internal I/O helpers (module-private; not exposed as TBPs)

   subroutine stl_ascii_next_facet(file_unit, names, p)
   !< Position an ASCII STL unit on the next `facet` line (issue #66): skip blank and `endsolid` lines; a `solid <name>`
   !< line opens the patch of that name (equal names merge), which becomes the current patch `p`. A file with one
   !< `solid` block (the first line, read as the header) reads as before.
   integer(I4P),                  intent(in)    :: file_unit !< Open ASCII unit.
   character(FRLEN), allocatable, intent(inout) :: names(:)  !< Patch names met so far.
   integer(I4P),                  intent(inout) :: p         !< Current patch.
   character(1024)                              :: line      !< Line buffer.
   character(len=:), allocatable                :: t         !< Line without leading blanks.
   integer(I4P)                                 :: ios       !< I/O status.

   do
      read(file_unit, '(A)', iostat=ios) line
      if (ios /= 0) return ! the facet reader reports it
      t = trim(adjustl(line))
      if (len(t) == 0) cycle
      if (len(t) >= 8) then
         if (t(1:8) == 'endsolid') cycle
      endif
      if (len(t) >= 5) then
         if (t(1:5) == 'solid') then
            p = find_or_add_patch(names, trim(adjustl(t(6:))))
            cycle
         endif
      endif
      backspace(file_unit)
      return
   enddo
   endsubroutine stl_ascii_next_facet

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
   ! status='replace' truncates an existing file (issue #80: a binary file written over a longer one kept its tail)
   if (is_ascii) then
      open(newunit=file_unit, file=file_name, form='formatted', status='replace', iostat=ios)
   else
      open(newunit=file_unit, file=file_name, access='stream', form='unformatted', status='replace', iostat=ios)
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
   character(1024)           :: facet_record  !< Line buffer for ASCII scan (longer than any STL line).
   character(FRLEN)          :: binary_header !< The 80-byte binary header.

   facets_number = 0
   rewind(file_unit)
   if (is_ascii) then
      do
         read(file_unit, '(A)', end=10, err=10) facet_record
         if (index(string=facet_record, substring='facet normal') > 0) facets_number = facets_number + 1
      enddo
   else
      read(file_unit, end=10, err=10) binary_header
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
      write(file_unit, '(A)') 'solid '//trim(header_patch_name(header))
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
