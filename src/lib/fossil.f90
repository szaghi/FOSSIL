!< FOSSIL, FOrtran Stereo (si) Litography parser.

module fossil
!< FOSSIL, FOrtran Stereo (si) Litography parser.

use fossil_facet_object
use fossil_surface_stl_object

implicit none
private
public :: facet_object
public :: surface_stl_object
public :: SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE
public :: sign_algorithm_from_string
public :: STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
endmodule fossil
