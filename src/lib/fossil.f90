!< FOSSIL, FOrtran Stereo (si) Litography parser.

module fossil
!< FOSSIL, FOrtran Stereo (si) Litography parser.

use fossil_aabb_tree_object, only : AABB_AUTO_REFINEMENT, AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use fossil_decimate, only : DEC_STATUS_OK, DEC_STATUS_BAD_INPUT, DEC_STATUS_NO_PROGRESS
use fossil_ray_query, only : ray_hit_t, RAY_STATUS_OK, RAY_STATUS_BAD_INPUT
use fossil_remesh, only : REM_STATUS_OK, REM_STATUS_BAD_INPUT
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
endmodule fossil
