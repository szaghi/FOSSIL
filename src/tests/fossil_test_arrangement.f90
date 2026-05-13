!< FOSSIL test: mesh-arrangement scaffolding (issue #18 §1.1 stage 1).
!<
!< Three invariants:
!<   1. Two non-overlapping cubes: arrangement reports zero cuts on every
!<      facet. The cross-mesh broad phase + narrow phase must produce no
!<      false positives on disjoint inputs.
!<   2. Two unit cubes offset along +x by 0.5: the four side faces of each
!<      cube are crossed by the other cube's side faces. The total cut count
!<      on each cube must be > 0 (qualitative invariant — exact counts depend
!<      on the per-facet cube tessellation, which is 2 triangles per face).
!<   3. project_to_plane / lift_from_plane round-trip: lifting then projecting
!<      a known point reproduces it within floating-point tolerance, AND
!<      projecting then lifting a 3D point that lies on the facet plane
!<      reproduces it within tolerance.

program fossil_test_arrangement

use fossil, only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_arrangement, only : arrangement_t, cut_list_t, &
                               arrangement_initialize, &
                               arrangement_collect_intersections, &
                               project_to_plane, lift_from_plane
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube_a, cube_b
type(arrangement_t)      :: arr
type(facet_object)       :: f
type(facet_object), pointer :: fp_a(:), fp_b(:)
type(vector_R8P) :: p_orig, p_back, p_in
real(R8P)        :: u, v
integer(I4P)     :: i, total_a_cuts, total_b_cuts
logical          :: are_tests_passed(3)

real(R8P), parameter :: TOL_ROUNDTRIP = 1.0e-12_R8P

are_tests_passed = .false.

! ----- 1. Two non-overlapping cubes (B translated by +10 in x).
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(10._R8P, 0._R8P, 0._R8P))
call cube_b%analyze
fp_a => cube_a%facets_ref()
fp_b => cube_b%facets_ref()
call arrangement_initialize(arr=arr, facet_a=fp_a, tree_a=cube_a%aabb, &
                                     facet_b=fp_b, tree_b=cube_b%aabb)
call arrangement_collect_intersections(arr=arr)
total_a_cuts = 0 ; total_b_cuts = 0
do i = 1, arr%n_a ; total_a_cuts = total_a_cuts + arr%cut(i)%n_segments ; enddo
do i = 1, arr%n_b ; total_b_cuts = total_b_cuts + arr%cut(arr%n_a + i)%n_segments ; enddo
are_tests_passed(1) = (total_a_cuts == 0 .and. total_b_cuts == 0)
print '(A,I0,A,I0)', 'disjoint cubes: A cuts=', total_a_cuts, '  B cuts=', total_b_cuts

! ----- 2. Two unit cubes offset by (0.5, 0, 0): they overlap along x in [0.5, 1].
call cube_a%destroy ; call cube_b%destroy
call cube_a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_b%translate(delta=vector_R8P(0.5_R8P, 0._R8P, 0._R8P))
call cube_b%analyze
fp_a => cube_a%facets_ref()
fp_b => cube_b%facets_ref()
call arrangement_initialize(arr=arr, facet_a=fp_a, tree_a=cube_a%aabb, &
                                     facet_b=fp_b, tree_b=cube_b%aabb)
call arrangement_collect_intersections(arr=arr)
total_a_cuts = 0 ; total_b_cuts = 0
do i = 1, arr%n_a ; total_a_cuts = total_a_cuts + arr%cut(i)%n_segments ; enddo
do i = 1, arr%n_b ; total_b_cuts = total_b_cuts + arr%cut(arr%n_a + i)%n_segments ; enddo
are_tests_passed(2) = (total_a_cuts > 0 .and. total_b_cuts > 0 .and. &
                       total_a_cuts == total_b_cuts)  ! every segment lands on both meshes
print '(A,I0,A,I0)', 'overlapping cubes: A cuts=', total_a_cuts, '  B cuts=', total_b_cuts

! ----- 3. project / lift round-trip on a synthetic facet.
f%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
f%vertex(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)
f%vertex(3) = vector_R8P(0._R8P, 3._R8P, 0._R8P)
call f%compute_metrix
! 3a: project a known in-plane point, then lift back.
p_orig = vector_R8P(0.5_R8P, 1.0_R8P, 0._R8P)  ! in z=0 plane through f's vertices
call project_to_plane(facet=f, p3d=p_orig, u=u, v=v)
call lift_from_plane(facet=f, u=u, v=v, p3d=p_back)
! 3b: lift a known (u, v), then project back.
call lift_from_plane(facet=f, u=1.5_R8P, v=2.0_R8P, p3d=p_in)
call project_to_plane(facet=f, p3d=p_in, u=u, v=v)
are_tests_passed(3) = abs(p_back%x - p_orig%x) <= TOL_ROUNDTRIP .and. &
                      abs(p_back%y - p_orig%y) <= TOL_ROUNDTRIP .and. &
                      abs(p_back%z - p_orig%z) <= TOL_ROUNDTRIP .and. &
                      abs(u - 1.5_R8P) <= TOL_ROUNDTRIP .and. &
                      abs(v - 2.0_R8P) <= TOL_ROUNDTRIP

print '(A,3F8.4,A,3F8.4)', 'roundtrip: orig=', p_orig%x, p_orig%y, p_orig%z, &
                           '  back=', p_back%x, p_back%y, p_back%z
print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_arrangement
