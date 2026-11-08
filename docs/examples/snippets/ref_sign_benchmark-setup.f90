use fossil,              only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                SIGN_FAST_WINDING_NUMBER
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = 4._R8P * atan(1._R8P), HALF = 15._R8P * PI / 180._R8P  ! half opening
integer(I4P), parameter         :: ND = 12                                           ! decades 1e-1 .. 1e-12
integer(I4P), parameter         :: ALGO(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                              SIGN_FAST_WINDING_NUMBER]
character(9), parameter         :: VARIANT(3) = ['clean    ', 'hole     ', 'flipped  ']
type(surface_stl_object)        :: wedge
type(facet_object), allocatable :: f(:)
type(vector_R8P),   allocatable :: p(:)
logical,            allocatable :: truth(:)
integer(I4P)                    :: v, a, i, wrong(5), ne, per_decade(ND)

call make_points(p, truth, ne)
print '(A,I0,A,I0,A)', 'wedge of 30 degrees, ', size(p), ' points (', ne, ' at the apex edge)'
print '(A)', 'variant    pseudo-normal  ray  solid angle  winding  naive normal'
