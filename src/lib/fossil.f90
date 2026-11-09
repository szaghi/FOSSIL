!< FOSSIL, FOrtran Surfaces, Solids, Immersed boundaries and Level sets.

module fossil
!< FOSSIL, FOrtran Surfaces, Solids, Immersed boundaries and Level sets.

use fossil_aabb_tree_object, only : AABB_AUTO_REFINEMENT, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use fossil_decimate, only : DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS
use fossil_ray_query, only : ray_hit_t, RAY_STATUS_OK, RAY_STATUS_BAD_INPUT
use fossil_remesh, only : REM_STATUS_OK, REM_STATUS_BAD_INPUT
use fossil_sdf, only : SDF_STATUS_OK, SDF_STATUS_BAD_INPUT, &
                       SDF_LABEL_UNASSIGNED, SDF_SENTINEL
use fossil_alpha_wrap, only : AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, &
                              AWRAP_STATUS_DEGENERATE, AWRAP_STATUS_NOT_CONVERGED
use fossil_csr_matrix, only : csr_matrix_t, &
                              CSR_STATUS_OK, CSR_STATUS_BAD_INPUT, CSR_STATUS_OUT_OF_RANGE
use fossil_laplacian, only : LAPL_STATUS_OK, LAPL_STATUS_BAD_INPUT, &
                             LAPL_STATUS_DEGENERATE_TRIANGLE
use fossil_curvature, only : CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, &
                             CURV_STATUS_DEGENERATE_TRIANGLE, CURV_STATUS_FIT_FAILED
use fossil_smoothing, only : SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, &
                             SMOOTH_STATUS_DEGENERATE, &
                             SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN, &
                             SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU, &
                             SMOOTH_DEFAULT_ITERATIONS
use fossil_facet_object
use fossil_marching_cubes, only : extract_isosurface, MC_STATUS_OK, MC_STATUS_BAD_DIMENSIONS
use fossil_surface_stl_object
use fossil_distance_field, only : distance_field_object, FAR_VALUE
use fossil_curve_slice, only : slice_surface
use fossil_cut_cells_3d, only : cut_cells_3d_t, cut_cells_3d
use fossil_cut_coarsen, only : cut_coarsen_3d, cut_coarsen_2d
use fossil_cut_moments, only : MOMENTS_MAX, mono_count, mono_exponents
use fossil_cut_graph, only : cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, &
                            SRD_NO_VOLUME, frd_t, frd_neighbourhoods, merge_t, merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, &
                            MERGE_CYCLE, MERGE_ALONE
use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, detect_topology, &
                                     detect_topology_2d, record_map, TOPO_BOUNDARY, TOPO_ABUT, TOPO_PERIODIC, TOPO_SINGULAR, &
                                     BLOCK_STATUS_OK, BLOCK_STATUS_BAD_INPUT, locate_in_blocks, locate_in_blocks_2d
use fossil_block_locate, only : LOC_INSIDE, LOC_OUTSIDE, LOC_WARPED, LOC_DEGENERATE, LOC_NOT_CONVERGED
use fossil_block_stencil, only : STENCIL_CENTRED, STENCIL_SHIFTED, STENCIL_REDUCED, STENCIL_FAILED, STENCIL_MAX
use fossil_block_geometry, only : CELL_VALID, CELL_SINGULAR, CELL_INVALID, CELL_UNDETERMINED
use fossil_overset, only : wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, block_real_t, &
                           hole_map, hole_map_2d, xray_flags, xray_flags_2d, donor_set_t, overlap_matrix, &
                           overlap_matrix_2d, donor_candidates, donor_candidates_2d, PAIR_NONE, PAIR_OVERLAP, PAIR_SHARED, &
                           overset_mask_t, classify, classify_2d, SELECT_EXPLICIT, SELECT_IMPLICIT, fringe_levels, &
                           fringe_levels_2d, overlap_width, overlap_width_2d, ORPHAN_NONE, ORPHAN_NO_CANDIDATE, &
                           ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR, ORPHAN_BOUNDARY, ORPHAN_WALL, FRINGE_INNER_1, &
                           FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2, wall_projection, wall_projection_2d, &
                           exposed_points, connectivity_t, build_connectivity, save_connectivity, &
                           cell_centre_block, cell_centre_block_2d
use fossil_ccd, only : motion_bound_3d, motion_bound_2d, ccd_evaluator_t, ccd_advance, CCD_CONTACT, CCD_SEPARATE, &
                       CCD_INITIAL_CONTACT, CCD_NOT_CONVERGED, CCD_INVALID_INPUT, CCD_MAX_ITERATIONS, CCD_STALL_STEPS, &
                       CCD_STALL_FRACTION
use fossil_projection, only : projection_target_t, PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER, project_points, &
                              project_points_2d, classify_points, classify_points_2d, PROJ_STATUS_OK, PROJ_STATUS_BAD_INPUT, &
                              project_block_face, project_block_side
use fossil_block_report, only : block_report_t, block_report, block_report_2d, save_block_report, utau_flat_plate
use fossil_hyperbolic, only : hyperbolic_grid_2d, geometric_spacing
use fossil_rigid_pose, only : rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, pose_from_quaternion, &
                             pose_2d, pose_compose, pose_inverse, to_body, to_world, rotate, pose_interpolate
use fossil_swept_cells, only : swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, swept_cells_3d, gauss_legendre_01
use fossil_mass_properties, only : inertia_to_world, principal_frame
use fossil_feature_edges, only : FEATURE_EDGE_SHARP, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD, &
                                 FEATURE_EDGE_DEFAULT_ANGLE
use fossil_thickness, only : THICKNESS_RAY, THICKNESS_MEDIAL, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS, &
                            THICKNESS_DEFAULT_SEPARATION
use fossil_deviation, only : DEV_STATUS_OK, DEV_STATUS_BAD_INPUT, DEV_STATUS_NOT_CONVERGED
use fossil_box_query, only : BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
use fossil_curve_object, only : curve_object, first_impact_2d, cut_cells_2d_t, CURVE_VALID, CURVE_TOO_FEW_VERTICES, &
                               CURVE_ZERO_LENGTH_SEGMENT, &
                               CURVE_ZERO_AREA_LOOP, CURVE_WRONG_ORIENTATION, CURVE_OPEN_LOOP, CURVE_SELF_INTERSECTING

implicit none
private
public :: facet_object
public :: surface_stl_object, proximity_pairs_of, first_impact
public :: distance_field_object, FAR_VALUE
public :: slice_surface
public :: cut_cells_3d_t, cut_cells_3d
public :: cut_coarsen_3d, cut_coarsen_2d, MOMENTS_MAX, mono_count, mono_exponents
public :: cut_links_t, cut_links_2d, srd_t, srd_neighbourhoods, SRD_REGULAR, SRD_COMPLETE, SRD_INCOMPLETE, SRD_NO_VOLUME
public :: frd_t, frd_neighbourhoods, merge_t, merge_targets, MERGE_NONE, MERGE_OK, MERGE_CHAIN, MERGE_CYCLE, MERGE_ALONE
public :: curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, detect_topology, detect_topology_2d, record_map
public :: TOPO_BOUNDARY, TOPO_ABUT, TOPO_PERIODIC, TOPO_SINGULAR, BLOCK_STATUS_OK, BLOCK_STATUS_BAD_INPUT
public :: CELL_VALID, CELL_SINGULAR, CELL_INVALID, CELL_UNDETERMINED
public :: locate_in_blocks, locate_in_blocks_2d, LOC_INSIDE, LOC_OUTSIDE, LOC_WARPED, LOC_DEGENERATE, LOC_NOT_CONVERGED
public :: STENCIL_CENTRED, STENCIL_SHIFTED, STENCIL_REDUCED, STENCIL_FAILED, STENCIL_MAX
public :: wall_cutter, wall_cutter_2d, hole_flags, hole_flags_2d, block_logical_t, block_real_t
public :: hole_map, hole_map_2d, xray_flags, xray_flags_2d
public :: donor_set_t, overlap_matrix, overlap_matrix_2d, donor_candidates, donor_candidates_2d, PAIR_NONE, PAIR_OVERLAP
public :: PAIR_SHARED, overset_mask_t, classify, classify_2d, SELECT_EXPLICIT, SELECT_IMPLICIT, fringe_levels
public :: fringe_levels_2d, overlap_width, overlap_width_2d, ORPHAN_NONE, ORPHAN_NO_CANDIDATE, ORPHAN_STENCIL_HOLE
public :: ORPHAN_RECEPTOR_DONOR, ORPHAN_BOUNDARY, ORPHAN_WALL, FRINGE_INNER_1, FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
public :: wall_projection, wall_projection_2d, exposed_points, connectivity_t, build_connectivity, save_connectivity
public :: cell_centre_block, cell_centre_block_2d
public :: rigid_pose_t, rigid_pose_2d_t, pose_identity, pose_from_axis_angle, pose_from_quaternion, pose_2d, pose_compose
public :: pose_inverse, to_body, to_world, rotate, pose_interpolate
public :: motion_bound_3d, motion_bound_2d, ccd_evaluator_t, ccd_advance, CCD_CONTACT, CCD_SEPARATE, CCD_INITIAL_CONTACT
public :: CCD_NOT_CONVERGED, CCD_INVALID_INPUT, CCD_MAX_ITERATIONS, CCD_STALL_STEPS, CCD_STALL_FRACTION
public :: projection_target_t, PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER, project_points, project_points_2d
public :: classify_points, classify_points_2d, PROJ_STATUS_OK, PROJ_STATUS_BAD_INPUT, project_block_face, project_block_side
public :: block_report_t, block_report, block_report_2d, save_block_report, utau_flat_plate
public :: hyperbolic_grid_2d, geometric_spacing
public :: swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, swept_cells_3d, gauss_legendre_01
public :: inertia_to_world, principal_frame
public :: FEATURE_EDGE_SHARP, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD, FEATURE_EDGE_DEFAULT_ANGLE
public :: DEV_STATUS_OK, DEV_STATUS_BAD_INPUT, DEV_STATUS_NOT_CONVERGED
public :: THICKNESS_RAY, THICKNESS_MEDIAL, THICKNESS_DEFAULT_CONE, THICKNESS_DEFAULT_RAYS, THICKNESS_DEFAULT_SEPARATION
public :: curve_object, first_impact_2d, cut_cells_2d_t, CURVE_VALID, CURVE_TOO_FEW_VERTICES, CURVE_ZERO_LENGTH_SEGMENT, &
          CURVE_ZERO_AREA_LOOP, &
          CURVE_WRONG_ORIENTATION, CURVE_OPEN_LOOP, CURVE_SELF_INTERSECTING
public :: intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING, SELFINT_COPLANAR, SELFINT_EDGE_FOLD, SELFINT_SHEETS
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
public :: FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX
public :: BOX_CUT, BOX_NEAR, BOX_INSIDE, BOX_OUTSIDE
public :: sign_algorithm_from_string
public :: STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
public :: STATUS_INVALID_INPUT
public :: AABB_AUTO_REFINEMENT
public :: AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
public :: BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF
public :: BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT
public :: BOOL_STATUS_NOT_CLOSED, BOOL_STATUS_INVALID_INPUT
public :: extract_isosurface
public :: MC_STATUS_OK, MC_STATUS_BAD_DIMENSIONS
public :: DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS
public :: REM_STATUS_OK, REM_STATUS_BAD_INPUT
public :: ray_hit_t
public :: RAY_STATUS_OK, RAY_STATUS_BAD_INPUT
public :: SDF_STATUS_OK, SDF_STATUS_BAD_INPUT
public :: SDF_LABEL_UNASSIGNED, SDF_SENTINEL
public :: AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT
public :: AWRAP_STATUS_DEGENERATE, AWRAP_STATUS_NOT_CONVERGED
public :: csr_matrix_t
public :: CSR_STATUS_OK, CSR_STATUS_BAD_INPUT, CSR_STATUS_OUT_OF_RANGE
public :: LAPL_STATUS_OK, LAPL_STATUS_BAD_INPUT, LAPL_STATUS_DEGENERATE_TRIANGLE
public :: CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, CURV_STATUS_DEGENERATE_TRIANGLE, CURV_STATUS_FIT_FAILED
public :: SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, SMOOTH_STATUS_DEGENERATE
public :: SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN
public :: SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU, SMOOTH_DEFAULT_ITERATIONS
endmodule fossil
