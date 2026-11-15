!< FOSSIL test: decimate's max_error is a certified bound of the distance from the result to the input (roadmap #92
!< WP5 b).
!<
!< Before, max_error stopped the loop when the cheapest collapse moved its vertex farther than max_error from one of the
!< input planes it carries (the square root of the quadric error): a bound on planes, not on the surface. On dragon.stl
!< with max_error = 3e-4 of its diagonal (8.0e-3), the collapse that takes it from 6546 to 6544 facets puts its vertex
!< within 0.04 max_error of the input but leaves the midpoint of its edge to a ring vertex 1.15 max_error from it (a chord
!< across a curved patch); 25 collapses on, at NO_PROGRESS, the result is 1.7 max_error from the input. Now a collapse
!< is accepted only if every facet of its new one-ring is certified within max_error of the input (a branch and bound
!< against the input's tree), else skipped and the loop goes on.
!<   1. dragon.stl, max_error = 3e-4 of the diagonal, target 6545 (that collapse): the certified distance of the result
!<      from the input (tolerance max_error / 20) <= max_error: [0.83, 0.87] max_error (before: [1.18, 1.22]).
!<   2. The same to target 4 (the whole run, NO_PROGRESS at 6538 facets): [0.83, 0.87] (before: [1.70, 1.74] at 6534).
!< The distance is certified on the facets that are not input facets (an input facet is at distance 0): the dragon's own
!< folds keep the certificate of the whole surface from converging (2e6 triangles, upper 6.3 max_error).
!< Without max_error nothing changes (the perf hashes of fossil_val_perf guard that bit for bit).

program fossil_test_decimate_certified

use fossil, only : surface_stl_object, facet_object, DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: dragon, dec
type(vector_R8P)         :: diagonal
real(R8P)                :: max_error, lower, upper
integer(I4P)             :: status, dev_status
logical                  :: are_tests_passed(2)

are_tests_passed = .false.
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
diagonal = dragon%get_bmax() - dragon%get_bmin()
max_error = 3.e-4_R8P * diagonal%normL2()

! ----- 1. The collapse 6546 -> 6544.
dec = dragon
call dec%decimate(target_facets=6545_I4P, max_error=max_error, status=status)
call distance_from_input(dec, dragon, max_error, lower, upper, dev_status)
print '(A,I0,A,I0,A,F7.4,A,F7.4,A)', 'dragon to 6545: status ', status, ', facets ', dec%get_facets_number(), &
      ', distance from the input in [', lower / max_error, ', ', upper / max_error, '] max_error'
are_tests_passed(1) = status == DEC_STATUS_OK .and. dev_status == DEV_STATUS_OK .and. upper <= max_error

! ----- 2. The whole run.
dec = dragon
call dec%decimate(target_facets=4_I4P, max_error=max_error, status=status)
call distance_from_input(dec, dragon, max_error, lower, upper, dev_status)
print '(A,I0,A,I0,A,F7.4,A,F7.4,A)', 'dragon to 4: status ', status, ', facets ', dec%get_facets_number(), &
      ', distance from the input in [', lower / max_error, ', ', upper / max_error, '] max_error'
are_tests_passed(2) = status == DEC_STATUS_NO_PROGRESS .and. dev_status == DEV_STATUS_OK .and. &
                      upper <= max_error .and. dec%get_facets_number() < dragon%get_facets_number()

print '(A,2L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   subroutine distance_from_input(s, s0, max_error, lower, upper, status)
   !< Certified bounds of h(s, s0) (tolerance max_error / 20) over the facets of s that are not facets of s0: decimate
   !< keeps the surviving facets in their input order, so facet k of s, if an input facet, is one of the facets
   !< k .. k + (deleted facets) of s0. Every facet of s that is an input facet is at distance 0 from s0.
   type(surface_stl_object), intent(in)  :: s, s0         !< Result, input.
   real(R8P),                intent(in)  :: max_error     !< max_error of the decimation.
   real(R8P),                intent(out) :: lower, upper  !< Bounds of h(s, s0).
   integer(I4P),             intent(out) :: status        !< DEV_STATUS_*.
   type(surface_stl_object)              :: part          !< Facets of s that are not facets of s0.
   type(facet_object), allocatable       :: f(:), g(:)    !< Their copies.
   integer(I4P)                          :: k, j, n, gap  !< Counters.
   logical                               :: old           !< An input facet.

   allocate(f(s%get_facets_number()))
   gap = s0%get_facets_number() - s%get_facets_number()
   n = 0
   do k = 1, s%get_facets_number()
      old = .false.
      do j = k, min(k + gap, s0%get_facets_number())
         old = same(s%facet_at(k), s0%facet_at(j))
         if (old) exit
      enddo
      if (old) cycle
      n = n + 1 ; f(n) = s%facet_at(k)
   enddo
   lower = 0._R8P ; upper = 0._R8P ; status = DEV_STATUS_OK
   if (n == 0) return
   allocate(g(n)) ; g = f(1:n)
   call part%adopt_facets(facets=g)
   call part%deviation(s0, lower=lower, upper=upper, tolerance=max_error / 20._R8P, status=status)
   endsubroutine distance_from_input

   logical function same(a, b)
   !< The same three vertices in the same order, bit for bit.
   type(facet_object), intent(in) :: a, b !< Facets.
   integer(I4P)                   :: j    !< Counter.

   same = .true.
   do j = 1, 3
      same = same .and. a%vertex(j)%x == b%vertex(j)%x .and. a%vertex(j)%y == b%vertex(j)%y .and. &
             a%vertex(j)%z == b%vertex(j)%z
   enddo
   endfunction same
endprogram fossil_test_decimate_certified
