!< FOSSIL test: report per-facet storage and pool footprint after stage 3c.
!<
!< Issue #5 motivation check. Prints sizeof(facet_object) and the per-facet
!< amortized pool overhead on a real STL (dragon). Asserts only that totals are
!< finite and that the pool covers exactly the facet vertices it should.

program fossil_test_facet_size

use fossil,                    only : surface_stl_object
use fossil_facet_object,       only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none

type(facet_object)                :: probe
type(surface_stl_object)          :: surface
type(vertex_pool_object), pointer :: pool
integer                           :: facet_bits
integer(I4P)                      :: nf, nv
real(R8P)                         :: facet_bytes, pool_bytes, total_bytes_per_facet
logical                           :: are_tests_passed(2)

are_tests_passed = .false.

facet_bits = storage_size(probe)
facet_bytes = real(facet_bits, R8P) / 8._R8P

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
nf = surface%get_facets_number()
pool => surface%get_vertex_pool()
nv = pool%vertex_count()

! Pool footprint estimate (dominant terms):
!   coord_:        nv * 24 B   (three R8P per unique vertex)
!   facet_to_pool: 3 * nf * 4 B
!   at_offset:     (nv + 1) * 4 B
!   at_pairs:      2 * 3 * nf * 4 B
pool_bytes = real(nv, R8P) * 24._R8P              &
           + real(3 * nf, R8P) * 4._R8P           &
           + real(nv + 1, R8P) * 4._R8P           &
           + real(2 * 3 * nf, R8P) * 4._R8P

total_bytes_per_facet = facet_bytes + pool_bytes / real(nf, R8P)

print '(A,I0)',          'dragon.stl facets:                    ', nf
print '(A,I0)',          'dragon.stl pool unique vertices:      ', nv
print '(A,F8.4)',        'V / F ratio:                          ', real(nv, R8P) / real(nf, R8P)
print '(A,F8.1,A)',      'sizeof(facet_object):                 ', facet_bytes, ' B'
print '(A,F8.1,A)',      'pool overhead (dominant):             ', pool_bytes, ' B total'
print '(A,F8.1,A)',      'pool overhead amortized per facet:    ', pool_bytes / real(nf, R8P), ' B'
print '(A,F8.1,A)',      'effective bytes per facet (facet + amortized pool): ', total_bytes_per_facet, ' B'

are_tests_passed(1) = facet_bytes > 0._R8P .and. facet_bytes < 1000._R8P
are_tests_passed(2) = nv > 0 .and. nv < 3 * nf

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_facet_size
