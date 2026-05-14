!< FOSSIL, FOrtran Stereo (si) Litography parser.

module fossil
!< FOSSIL, FOrtran Stereo (si) Litography parser.

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
                             CURV_STATUS_DEGENERATE_TRIANGLE
use fossil_smoothing, only : SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, &
                             SMOOTH_STATUS_DEGENERATE, &
                             SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN, &
                             SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU, &
                             SMOOTH_DEFAULT_ITERATIONS
use fossil_facet_object
use fossil_marching_cubes, only : extract_isosurface, MC_STATUS_OK, MC_STATUS_BAD_DIMENSIONS
use fossil_surface_stl_object

implicit none
private
public :: facet_object
public :: surface_stl_object
public :: intersection_pair_t
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, SIGN_PSEUDO_NORMAL
public :: sign_algorithm_from_string
public :: STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
public :: STATUS_INVALID_INPUT
public :: AABB_AUTO_REFINEMENT
public :: AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
public :: BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF
public :: BOOL_STATUS_OK, BOOL_STATUS_CDT_FAILED, BOOL_STATUS_NOT_IMPLEMENTED, BOOL_STATUS_EMPTY_INPUT
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
public :: CURV_STATUS_OK, CURV_STATUS_BAD_INPUT, CURV_STATUS_DEGENERATE_TRIANGLE
public :: SMOOTH_STATUS_OK, SMOOTH_STATUS_BAD_INPUT, SMOOTH_STATUS_DEGENERATE
public :: SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN
public :: SMOOTH_DEFAULT_LAMBDA, SMOOTH_DEFAULT_MU, SMOOTH_DEFAULT_ITERATIONS
endmodule fossil
