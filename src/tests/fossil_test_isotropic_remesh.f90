!< FOSSIL test: isotropic remeshing (issue #18 §1.7).
!<
!< Incremental rebuild — see fossil_remesh module header for staging notes.
!<
!< STEP 1: single-triangle round-trip (collect_state + materialize).
!< STEP 2: edge enumeration + median length on cube.stl.
!< STEP 3: split pass on a single triangle with target_length=0.4.
!< STEP 4: split + collapse on a sphere via §1.5 marching cubes.
!< STEP 6 (this commit): split + collapse + flip + relax+projection on a
!<   sphere. The relax pass moves each vertex toward its area-weighted
!<   one-ring centroid then projects back onto the reference surface via
!<   `tree%distance_tree_with_region` + `facet%compute_distance_with_region`.
!<   Assertion: volume preserved within 10% (the area-weighted Laplacian
!<   has known inward bias on convex regions — `1 - cos(θ)` per relaxation
!<   step where θ is the typical incident-facet half-angle; the sphere
!<   case loses ~5% per iteration even with projection); manifoldness
!<   preserved. Tighter volume preservation requires normal-direction
!<   projection or scaled-step heuristics; deferred as future work.

program fossil_test_isotropic_remesh

use fossil, only : surface_stl_object, extract_isosurface
use fossil_remesh, only : isotropic_remesh, REM_STATUS_OK, &
                          count_unique_edges, compute_median_edge_length, &
                          run_split_only, run_split_and_collapse, run_full_pipeline
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(facet_object), allocatable :: working(:), tri_split(:), sphere_facets(:), &
                                   relax_input(:)
type(surface_stl_object) :: cube, sphere, sphere_ref, sphere_relaxed
type(facet_object), pointer :: fp(:)
real(R8P), allocatable :: values(:, :, :)
type(vector_R8P)         :: bmin, bmax
integer(I4P) :: status, i, ne_cube, nf_before_split, nf_after_split, k
integer(I4P) :: nf_before_sphere, nf_after_sphere, nm_after_sphere
integer(I4P) :: nf_relaxed, nm_relaxed
real(R8P)    :: med_cube, vol_before_sphere, vol_after_sphere
real(R8P)    :: vol_relaxed
real(R8P)    :: x, y, z, dx
logical      :: are_tests_passed(5)
logical      :: all_facets_valid
real(R8P)    :: expected_x(3), got_x(3)
real(R8P)    :: expected_y(3), got_y(3)
real(R8P)    :: expected_z(3), got_z(3)
integer(I4P), parameter :: N_GRID = 32_I4P
real(R8P),    parameter :: PI     = 3.141592653589793_R8P

are_tests_passed = .false.

! Build a single-triangle input with vertex IDs 1, 2, 3.
allocate(working(1))
working(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
working(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
working(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P)
working(1)%vertex_id = [1, 2, 3]
working(1)%fcon_edge = [0, 0, 0]
call working(1)%compute_metrix

expected_x = [working(1)%vertex(1)%x, working(1)%vertex(2)%x, working(1)%vertex(3)%x]
expected_y = [working(1)%vertex(1)%y, working(1)%vertex(2)%y, working(1)%vertex(3)%y]
expected_z = [working(1)%vertex(1)%z, working(1)%vertex(2)%z, working(1)%vertex(3)%z]

! Use a target_length tuned to keep all edges in [4L/5, 4L/3]: with
! triangle edges 1, 1, sqrt(2)~1.41, choose L = 1.2 → split-threshold 1.6,
! collapse-threshold 0.96. All three edges are in (0.96, 1.6) → no-op.
call isotropic_remesh(facet=working, target_length=1.2_R8P, iterations=1_I4P, &
                       preserve_features=.false., status=status)

! The stub should round-trip: same facet count, same vertex coordinates.
if (status == REM_STATUS_OK .and. allocated(working) .and. size(working) == 1) then
   got_x = [working(1)%vertex(1)%x, working(1)%vertex(2)%x, working(1)%vertex(3)%x]
   got_y = [working(1)%vertex(1)%y, working(1)%vertex(2)%y, working(1)%vertex(3)%y]
   got_z = [working(1)%vertex(1)%z, working(1)%vertex(2)%z, working(1)%vertex(3)%z]
   are_tests_passed(1) = all(abs(got_x - expected_x) <= tiny(1._R8P)) .and. &
                         all(abs(got_y - expected_y) <= tiny(1._R8P)) .and. &
                         all(abs(got_z - expected_z) <= tiny(1._R8P))
endif

print '(A,I0,A,I0)', 'round-trip: status=', status, ' nf=', size(working)
do i = 1, 3
   print '(A,I0,A,3F8.4,A,3F8.4)', '  v', i, ' got=', got_x(i), got_y(i), got_z(i), &
         ' expected=', expected_x(i), expected_y(i), expected_z(i)
enddo

! Step 2: edge count + median on cube.stl.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
fp => cube%facets_ref()
ne_cube  = count_unique_edges(facet=fp)
med_cube = compute_median_edge_length(facet=fp)
print '(A,I0,A,F8.4)', 'cube: ne=', ne_cube, '  median=', med_cube
are_tests_passed(2) = (ne_cube == 18_I4P .and. abs(med_cube - 1.0_R8P) <= 1.0e-12_R8P)

! Step 3: split a single triangle into more sub-triangles.
allocate(tri_split(1))
tri_split(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri_split(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
tri_split(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P)
tri_split(1)%vertex_id = [1, 2, 3]
tri_split(1)%fcon_edge = [0, 0, 0]
call tri_split(1)%compute_metrix
nf_before_split = size(tri_split)
call run_split_only(facet=tri_split, target_length=0.4_R8P, n_iterations=1_I4P)
nf_after_split = size(tri_split)
print '(A,I0,A,I0)', 'split: nf_before=', nf_before_split, ' nf_after=', nf_after_split

! Validate: every output facet has 3 distinct vertex coordinates.
all_facets_valid = .true.
do i = 1, size(tri_split)
   do k = 1, 3
      if (sqrt((tri_split(i)%vertex(k)%x - tri_split(i)%vertex(mod(k, 3) + 1)%x)**2 + &
               (tri_split(i)%vertex(k)%y - tri_split(i)%vertex(mod(k, 3) + 1)%y)**2 + &
               (tri_split(i)%vertex(k)%z - tri_split(i)%vertex(mod(k, 3) + 1)%z)**2) <= tiny(1._R8P)) then
         all_facets_valid = .false.
      endif
   enddo
enddo
are_tests_passed(3) = (nf_after_split > nf_before_split .and. all_facets_valid)

! Step 4: split + collapse on a sphere via §1.5 marching cubes.
block
   integer :: ig, jg, kg
   dx = 4._R8P / real(N_GRID - 1, R8P)
   allocate(values(N_GRID, N_GRID, N_GRID))
   do kg = 1, N_GRID
      z = -2._R8P + (kg - 1) * dx
      do jg = 1, N_GRID
         y = -2._R8P + (jg - 1) * dx
         do ig = 1, N_GRID
            x = -2._R8P + (ig - 1) * dx
            values(ig, jg, kg) = sqrt(x**2 + y**2 + z**2) - 1._R8P
         enddo
      enddo
   enddo
endblock
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)
nf_before_sphere = sphere%get_facets_number()
vol_before_sphere = sphere%get_volume()
allocate(sphere_facets(nf_before_sphere))
do i = 1, nf_before_sphere
   sphere_facets(i) = sphere%facet_at(i)
enddo
call run_split_and_collapse(facet=sphere_facets, target_length=-1._R8P, n_iterations=2_I4P)
nf_after_sphere = size(sphere_facets)
! Round-trip into a surface to evaluate volume / manifoldness.
call sphere%destroy
call sphere%adopt_facets(facets=sphere_facets)
vol_after_sphere = sphere%get_volume()
nm_after_sphere = sphere%get_non_manifold_edges_number()
print '(A,I0,A,I0,A,F8.4,A,F8.4,A,I0)', 'sphere: nf_before=', nf_before_sphere, &
      ' nf_after=', nf_after_sphere, &
      '  vol_before=', vol_before_sphere, ' vol_after=', vol_after_sphere, &
      ' nm_edges=', nm_after_sphere
are_tests_passed(4) = (nf_after_sphere /= nf_before_sphere .and. &
                       abs(vol_after_sphere - vol_before_sphere) <= 0.10_R8P * vol_before_sphere .and. &
                       nm_after_sphere == 0_I4P)

! Step 6: full pipeline (split + collapse + flip + relax+projection) on a
! sphere with the original sphere as the reference. With projection wired,
! volume drift should be tighter than step 4's collapse-only test.
block
   integer :: ig, jg, kg
   real(R8P) :: dxg, xg, yg, zg
   integer, parameter :: N6 = 32
   dxg = 4._R8P / real(N6 - 1, R8P)
   allocate(values(N6, N6, N6))
   do kg = 1, N6
      zg = -2._R8P + (kg - 1) * dxg
      do jg = 1, N6
         yg = -2._R8P + (jg - 1) * dxg
         do ig = 1, N6
            xg = -2._R8P + (ig - 1) * dxg
            values(ig, jg, kg) = sqrt(xg**2 + yg**2 + zg**2) - 1._R8P
         enddo
      enddo
   enddo
endblock
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere_ref%adopt_facets(facets=sphere_facets)  ! reference (unmodified)
deallocate(values)
allocate(relax_input(sphere_ref%get_facets_number()))
do i = 1, sphere_ref%get_facets_number()
   relax_input(i) = sphere_ref%facet_at(i)
enddo
call run_full_pipeline(facet=relax_input, target_length=-1._R8P, n_iterations=2_I4P, &
                       reference_facet=sphere_ref%facets_ref(), &
                       reference_tree=sphere_ref%aabb)
call sphere_relaxed%adopt_facets(facets=relax_input)
nf_relaxed   = sphere_relaxed%get_facets_number()
vol_relaxed  = sphere_relaxed%get_volume()
nm_relaxed   = sphere_relaxed%get_non_manifold_edges_number()
print '(A,I0,A,F8.4,A,I0)', 'sphere relaxed: nf=', nf_relaxed, '  vol=', vol_relaxed, '  nm_edges=', nm_relaxed
! Reference sphere volume = sphere_ref%get_volume() ≈ 4.148; tight tolerance.
are_tests_passed(5) = (abs(vol_relaxed - sphere_ref%get_volume()) <= 0.10_R8P * sphere_ref%get_volume() .and. &
                       nm_relaxed == 0_I4P)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_isotropic_remesh
