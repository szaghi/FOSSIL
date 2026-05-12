!< FOSSIL test: facet%area() and surface%get_area() (issue #7).
!<
!< Four invariants:
!<   1. The reporter's exact triangle [0,0,0]/[1,0,0]/[0.5,1,0] has area 0.5
!<      and corresponding det = (2*area)^2 = 1.
!<   2. On a real mesh, surface%get_area() equals the sum of per-facet
!<      facet%area() to within accumulated FP round-off.
!<   3. Surface area is invariant under rigid motions (translate + rotate).
!<   4. Surface area scales as k^2 under uniform scaling by k.

program fossil_test_area

use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : ey_R8P, vector_R8P

implicit none

type(surface_stl_object) :: surface
type(facet_object)       :: single
real(R8P)                :: area_facet_sum, area_surface, area_baseline, area_scaled
real(R8P), parameter     :: HALF_PI = 1.5707963267948966_R8P
real(R8P), parameter     :: SCALE_K = 3.0_R8P
real(R8P), parameter     :: TOL_TRIANGLE = 1.0e-15_R8P
real(R8P), parameter     :: TOL_FP       = 1.0e-9_R8P
real(R8P), parameter     :: TOL_INVARIANT = 1.0e-6_R8P
integer(I4P)             :: f, nf
logical                  :: are_tests_passed(4)

are_tests_passed = .false.

! --- 1. Single-triangle exact area. The Gram det for this triangle is 1,
!         so facet%area() = 1/2 = 0.5 exactly in IEEE doubles.
single%vertex(1) = vector_R8P(0.0_R8P, 0.0_R8P, 0.0_R8P)
single%vertex(2) = vector_R8P(1.0_R8P, 0.0_R8P, 0.0_R8P)
single%vertex(3) = vector_R8P(0.5_R8P, 1.0_R8P, 0.0_R8P)
call single%compute_metrix
are_tests_passed(1) = (abs(single%area() - 0.5_R8P) <= TOL_TRIANGLE)

! --- 2. Sum of facet areas == surface area, on dragon.
call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
nf = surface%get_facets_number()
area_facet_sum = 0._R8P
do f = 1, nf
   associate (fp => surface%facet_at(f))
      area_facet_sum = area_facet_sum + fp%area()
   end associate
enddo
area_surface = surface%get_area()
are_tests_passed(2) = (abs(area_facet_sum - area_surface) <= TOL_FP * area_surface)

! --- 3. Area invariant under translate + rotate.
area_baseline = surface%get_area()
call surface%translate(delta=vector_R8P(100._R8P, -50._R8P, 25._R8P))
call surface%rotate(axis=ey_R8P, angle=HALF_PI, center=surface%get_centroid())
call surface%analyze
are_tests_passed(3) = (abs(surface%get_area() - area_baseline) <= TOL_INVARIANT * area_baseline)

! --- 4. Area scales as k^2 under uniform scaling. Reload to reset state.
call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
area_baseline = surface%get_area()
call surface%resize(factor=vector_R8P(SCALE_K, SCALE_K, SCALE_K), recompute_metrix=.true.)
call surface%analyze
area_scaled = surface%get_area()
are_tests_passed(4) = (abs(area_scaled - SCALE_K**2 * area_baseline) <= TOL_INVARIANT * area_scaled)

print '(A,F8.4,A,F8.4)',     'single triangle: area=', single%area(), '  det=', single%det
print '(A,ES14.6,A,ES14.6)', 'dragon: sum_facets=', area_facet_sum, '  surface=', area_surface
print '(A,4L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_area
