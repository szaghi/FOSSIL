!< FOSSIL, expected-failure test: invalid sign_algorithm triggers `error stop`.
!<
!< Regression test for audit S6 + D4-rest. The `case default` in `is_point_inside`
!< previously returned an undefined logical; it now raises `error stop`. The API also
!< switched from a string algorithm name to an integer enum — `sign_algorithm` is
!< `integer(I4P)` with valid codes `SIGN_RAY_INTERSECTIONS=1` and `SIGN_SOLID_ANGLE=2`.
!< This program passes an out-of-range integer (99) and expects a non-zero exit code.
!<
!< The surface needs no facets for this dispatch: the algorithm-code check happens
!< before any geometric work.

program fossil_xfail_bad_sign_algorithm

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface     !< Empty surface; sufficient for the dispatch test.
logical                  :: is_inside   !< Sink for the function result.

is_inside = surface%is_point_inside(point=0._R8P * ex_R8P, sign_algorithm=99_I4P)

! Unreachable: `is_point_inside` must `error stop` before returning.
print '(A)', 'FAIL: is_point_inside did not error_stop on bad sign_algorithm code'
print '(A,L1)', 'returned is_inside = ', is_inside
endprogram fossil_xfail_bad_sign_algorithm
