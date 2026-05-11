!< FOSSIL, expected-failure test: invalid sign_algorithm triggers `error stop`.
!<
!< Regression test for audit S6. Before the fix, the `case default` in
!< `is_point_inside` was a comment — passing an unknown algorithm name returned an
!< undefined logical. After the fix, an `error stop` is raised. This program is
!< expected to terminate with a non-zero exit code; the test harness inverts the
!< exit status for `fossil_xfail_*` programs.
!<
!< The surface must have facets before `is_point_inside` is meaningful, but the
!< algorithm-name check happens before any geometric work — an empty surface
!< suffices to exercise the dispatch.

program fossil_xfail_bad_sign_algorithm

use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : ex_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface     !< Empty surface; sufficient for the dispatch test.
logical                  :: is_inside   !< Sink for the function result.

is_inside = surface%is_point_inside(point=0._R8P * ex_R8P, sign_algorithm='not_a_real_algorithm')

! Unreachable: `is_point_inside` must `error stop` before returning.
print '(A)', 'FAIL: is_point_inside did not error_stop on bad sign_algorithm'
print '(A,L1)', 'returned is_inside = ', is_inside
endprogram fossil_xfail_bad_sign_algorithm
