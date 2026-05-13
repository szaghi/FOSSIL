!< FOSSIL test: self-intersection detection (issue #18 §1.2).
!<
!< Three invariants:
!<   1. Hand-crafted intersecting pair — build two facets that cross along a
!<      known segment, wrap them as a 2-facet surface, and assert the detector
!<      finds exactly one pair with segment endpoints matching the analytic
!<      intersection within tolerance.
!<   2. Clean meshes return zero intersections — cube.stl (12 facets, exact)
!<      and bunny.stl (69k facets, real-world but verified clean) must report
!<      zero pairs. This is the no-false-positive guarantee against the
!<      adjacency filter (every edge in a closed mesh is shared by two
!<      facets; if the filter were too narrow, every shared edge would
!<      register as a hit).
!<   3. Dirty mesh is detected non-empty — the dragon.stl fixture in this
!<      repo has known real self-intersections (~759 pairs at the time of
!<      writing). Asserting > 0 catches silent regressions where the broad
!<      phase or narrow phase falsely prunes everything. The exact count is
!<      not asserted because it depends on EPS-sensitive behavior at the
!<      narrow phase boundary; only its non-zero-ness is load-bearing.
!<
!< Note on dragon.stl: a future PR may either replace this fixture with a
!< cleaned version or use it as ground truth for §1.2 resolution. Either way,
!< invariant 3 documents the current state explicitly so it doesn't surprise.

program fossil_test_self_intersection

use fossil, only : surface_stl_object, intersection_pair_t
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: pair_surface, cube, dragon, bunny
type(intersection_pair_t), allocatable :: pairs(:)
type(facet_object), allocatable :: tri(:)
type(vector_R8P) :: p_expected_lo, p_expected_hi
type(vector_R8P) :: p_lo, p_hi
real(R8P) :: err_lo, err_hi
integer(I4P) :: status
integer(I4P) :: n_dragon
logical :: are_tests_passed(3)

real(R8P), parameter :: TOL_SEGMENT = 1.0e-12_R8P  ! analytic intersection — full FP precision

are_tests_passed = .false.

! ----- 1. Hand-crafted intersecting pair.
!         Triangle A: in z=0 plane, vertices (0,0,0), (2,0,0), (1,2,0).
!         Triangle B: in y=0.5 plane, vertices (0.5,0.5,-1), (1.5,0.5,-1), (1,0.5,1).
!         They cross along the segment y=0.5, z=0, x in [0.75, 1.25].
allocate(tri(2))
tri(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(3) = vector_R8P(1._R8P, 2._R8P, 0._R8P)
tri(1)%id        = 1_I4P
tri(2)%vertex(1) = vector_R8P(0.5_R8P, 0.5_R8P, -1._R8P)
tri(2)%vertex(2) = vector_R8P(1.5_R8P, 0.5_R8P, -1._R8P)
tri(2)%vertex(3) = vector_R8P(1.0_R8P, 0.5_R8P,  1._R8P)
tri(2)%id        = 2_I4P
call pair_surface%adopt_facets(facets=tri)

call pair_surface%find_self_intersections(pairs=pairs, status=status)

! Expected segment endpoints (order-agnostic).
p_expected_lo = vector_R8P(0.75_R8P, 0.5_R8P, 0._R8P)
p_expected_hi = vector_R8P(1.25_R8P, 0.5_R8P, 0._R8P)

if (status == 0_I4P .and. size(pairs) == 1) then
   ! Order the returned segment so its low-x endpoint comes first.
   if (pairs(1)%p%x <= pairs(1)%q%x) then
      p_lo = pairs(1)%p ; p_hi = pairs(1)%q
   else
      p_lo = pairs(1)%q ; p_hi = pairs(1)%p
   endif
   err_lo = sqrt((p_lo%x - p_expected_lo%x)**2 + (p_lo%y - p_expected_lo%y)**2 + (p_lo%z - p_expected_lo%z)**2)
   err_hi = sqrt((p_hi%x - p_expected_hi%x)**2 + (p_hi%y - p_expected_hi%y)**2 + (p_hi%z - p_expected_hi%z)**2)
   are_tests_passed(1) = pairs(1)%a == 1 .and. pairs(1)%b == 2 .and. &
                         err_lo <= TOL_SEGMENT .and. err_hi <= TOL_SEGMENT
endif

! ----- 2. Clean meshes return zero intersections.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%find_self_intersections(pairs=pairs, status=status)
if (.not. (status == 0_I4P .and. size(pairs) == 0)) goto 100

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny%find_self_intersections(pairs=pairs, status=status)
if (.not. (status == 0_I4P .and. size(pairs) == 0)) goto 100

are_tests_passed(2) = .true.

100 continue

! ----- 3. Dirty mesh (dragon) reports > 0 real intersections.
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%find_self_intersections(pairs=pairs, status=status)
n_dragon = size(pairs)
are_tests_passed(3) = (status == 0_I4P .and. n_dragon > 0_I4P)

print '(A,I0)', 'dragon: nfound=', n_dragon
print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_self_intersection
