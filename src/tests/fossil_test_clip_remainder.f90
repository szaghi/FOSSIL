!< FOSSIL test: clip remainder contains exactly the facets excluded from the clip (audit #14 S3).
!<
!< Verifies three invariants:
!<   1. facets_in + facets_out == original total
!<   2. remainder is non-empty
!<   3. remainder bmax%x > clip bmax%x (at least one facet is truly outside)

program fossil_test_clip_remainder

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: surface
type(surface_stl_object) :: remainder
type(vector_R8P)         :: bmin, bmax
integer(I4P)             :: n_orig, n_in, n_out
logical                  :: are_tests_passed(3)
type(vector_R8P)         :: rem_bmax

are_tests_passed = .false.

bmin = vector_R8P(-15._R8P, -5._R8P,  0._R8P)
bmax = vector_R8P(  0._R8P,  5._R8P, 20._R8P)

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
n_orig = surface%get_facets_number()

call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)
n_in  = surface%get_facets_number()
n_out = remainder%get_facets_number()

! 1. partition invariant: in + out == original
are_tests_passed(1) = (n_in + n_out == n_orig)

! 2. remainder is non-empty
are_tests_passed(2) = (n_out > 0)

! 3. remainder has at least one facet beyond the clip bmax on x
rem_bmax = remainder%get_bmax()
are_tests_passed(3) = (rem_bmax%x > bmax%x)

print '(A,I0,A,I0,A,I0)', 'facets: original=', n_orig, '  in=', n_in, '  out=', n_out
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_clip_remainder
