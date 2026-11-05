!< VALIDATION version (issue #45): the full, slow test, run by scripts/run_validation.sh at -O2.
!< The fast regression variant run per push is src/tests/fossil_test_alpha_wrap.f90 (no invariants 11 and
!< 14 on the sphere wrap, sphere sampled on a 20^3 grid).
!<
!< FOSSIL test: alpha-wrap octree (issue #18 §1.6, step 1).
!<
!< Step 1: validate the leaf-size-α octree built by `awrap_build_octree`.
!< Step 2: validate the inside/outside flood fill from `awrap_classify_leaves`.
!< Step 3: validate the dual-contour boundary extraction
!<         from `awrap_extract_surface`.
!< Step 4: validate vertex projection from `awrap_project_vertices`.
!< Step 5: validate the capstone `surface%alpha_wrap` TBP — adaptive
!<         refinement loop integrating steps 1-4.
!<
!< 3 invariants for step 1, 3 for step 2, 3 for step 3, 3 for step 4, 7 for the wrap of awrap_run:
!<   1. Cube — refinement structure. cube.stl is a closed unit cube. Build
!<      octree with α = 0.1. Assert: n_boundary_leaves > 0 (geometry is
!<      captured), every boundary leaf actually overlaps at least one input
!<      facet (no false positives), every non-boundary leaf does NOT overlap
!<      any input facet (no false negatives — the load-bearing property
!<      step 2 will rely on for its flood-fill barriers), max depth ≤
!<      AWRAP_MAX_DEPTH, all boundary leaves are at the target size (size ≤
!<      α) since they reached the splitting stop condition.
!<   2. Sphere — refinement structure on a curved boundary. Sphere built
!<      from §1.5 marching cubes. Same property checks as inv 1, plus
!<      verify the boundary leaves form a "shell" — no boundary leaf is
!<      strictly interior (centroid distance from the surface > α).
!<   3. Degenerate input — α larger than bbox diagonal. Status =
!<      AWRAP_STATUS_DEGENERATE; the octree still has 1 root leaf marked
!<      BOUNDARY. Empty facet array → AWRAP_STATUS_BAD_INPUT.
!<   4. Cube — flood fill correctness on a closed solid. After
!<      `awrap_classify_leaves`: every empty leaf is now INSIDE or OUTSIDE
!<      (no EMPTY remaining), n_inside + n_outside + n_boundary = n_leaves,
!<      every INSIDE leaf's centroid lies inside [0,1]^3 (geometric truth
!<      for the cube fixture).
!<   5. Sphere — flood fill on a curved closed surface. Every INSIDE leaf
!<      centroid has distance < R from origin (R = 1 for the unit MC
!<      sphere); every OUTSIDE centroid has distance > R - α (some
!<      conservatism for boundary thickness). No leaks.
!<   6. Cube with one facet deleted — the EXPECTED-LEAK test. The hole
!<      lets the flood fill leak through, so n_inside drops drastically
!<      compared to a closed cube (assert ≤ 30% of closed-cube n_inside).
!<      This is NOT a bug — step 3+ closes the wrap surface across the
!<      hole, NOT by repairing the input. Without this expected-leak the
!<      algorithm is broken.
!<   7. Cube wrap surface — extract via awrap_extract_surface, adopt as a
!<      surface, assert is_watertight() AND volume within 30% of unit cube
!<      (axis-aligned voxel-y wrap is bigger than the input by up to one
!<      α-thick layer). The watertightness is the load-bearing property:
!<      a 2-manifold output by construction is the whole point of this step.
!<   8. Sphere wrap surface — same on a curved input. Watertight + volume
!<      within 30% of (4/3)π. Voxel-y at α=0.1 → expect ~10% inflation.
!<   9. Holed cube wrap surface — n_facets > 0 (something is emitted) but
!<      step 3 alone CANNOT close the hole (no INSIDE leaves means the
!<      wrap emits twin quads on both sides of the boundary band; topology
!<      is degenerate). Step 5's adaptive refinement closes it. Asserts
!<      facet count > 0 and skips watertightness — documented gap.
!<  10. Cube wrap + projection — volume drops closer to the analytic value
!<      (1.0 for unit cube) compared to the unprojected step-3 wrap.
!<      Watertightness preserved (vertex dedup is doing its job).
!<      offset = α/3 per CGAL's recommended ratio.
!<  11. Sphere wrap + projection — same on a curved input. Volume drops
!<      closer to (4/3)π. Watertightness preserved.
!<  12. Projection manifold preservation. Compute non_manifold_edges_number
!<      for the projected wrap of the cube; must be zero (vertex dedup
!<      ensures shared vertices stay shared after projection, so no
!<      seams open up).
!<  13-19: the wrap of awrap_run (issue #50); every wrap must be status OK,
!<      watertight, manifold and free of self-crossings (find_self_intersections).
!<  13. Cube, alpha 0.1: volume within 0.5% of the offset cube (Minkowski sum).
!<  14. Marching-cubes sphere: >= 90% of the vertices on the offset surface, none closer than offset / 2.
!<  15. Holed cube, hole wider than the ball: a hollow shell lining both sides.
!<  16. Holed cube, hole narrower than the ball: capped.
!<  17. One triangle: a closed shell around both sides.
!<  18. The issue #50 table: cube-01 closed and minus one facet, alpha 0.1, 0.3, 0.7.
!<  19. dragon.stl with the defaults; the wrap encloses it, vertices as in 14.

program fossil_val_alpha_wrap

use fossil, only : surface_stl_object, extract_isosurface, intersection_pair_t
use fossil_alpha_wrap, only : awrap_octree_t, awrap_build_octree, awrap_classify_leaves, &
                              awrap_extract_surface, awrap_project_vertices, &
                              AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE, &
                              AWRAP_LEAF_FLAG_BOUNDARY, AWRAP_LEAF_FLAG_INTERIOR, &
                              AWRAP_LEAF_FLAG_EMPTY, AWRAP_LEAF_FLAG_INSIDE, AWRAP_LEAF_FLAG_OUTSIDE, &
                              AWRAP_MAX_DEPTH
use fossil_facet_object, only : facet_object
use fossil_utils, only : PI, triangle_overlaps_aabb
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere, cube_holed, wrap_surface, wrapped_tbp, plate, cube01, cube01_holed, dragon
type(facet_object), allocatable :: plate_facets(:)
real(R8P), parameter :: TABLE_ALPHA(3) = [0.1_R8P, 0.3_R8P, 0.7_R8P]
type(facet_object), allocatable :: sphere_facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
type(facet_object), allocatable :: wrap_facets(:)
type(facet_object), pointer :: facets(:)
type(awrap_octree_t) :: octree
real(R8P), allocatable :: values(:,:,:)
type(vector_R8P) :: bmin, bmax, centroid
real(R8P) :: x, y, z, dx, alpha
real(R8P) :: max_leaf_size, leaf_sx, leaf_sy, leaf_sz
real(R8P) :: wrap_volume, expected_volume
real(R8P) :: vol_before, vol_after, offset
integer(I4P) :: status, i, j, k, n, max_depth
integer(I4P) :: n_false_pos, n_false_neg
integer(I4P) :: n_in_geometry, n_outside_geometry, n_empty_remaining
integer(I4P) :: n_inside_closed_cube
integer(I4P) :: n_wrap_facets, n_nm_edges
logical :: any_overlap, prop_no_false_pos, prop_no_false_neg, prop_size_ok
logical :: are_tests_passed(19)

integer(I4P), parameter :: N_GRID = 32_I4P

are_tests_passed = .false.

! ---- Invariant 1: cube — refinement structure properties.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
facets => cube%facets_ref()
alpha = 0.1_R8P
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call check_node_properties(octree=octree, facets=facets, alpha=alpha, &
                           max_depth=max_depth, max_leaf_size=max_leaf_size, &
                           n_false_pos=n_false_pos, n_false_neg=n_false_neg)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', &
      'inv 1 (cube): n_nodes=', octree%n_nodes, ' n_leaves=', octree%n_leaves, &
      ' n_boundary=', octree%n_boundary_leaves, ' max_depth=', max_depth, &
      ' n_false_pos=', n_false_pos
print '(A,F8.4,A,I0)', '              max_leaf_size=', max_leaf_size, ' n_false_neg=', n_false_neg
are_tests_passed(1) = (status == AWRAP_STATUS_OK)              .and. &
                      (octree%n_boundary_leaves > 0_I4P)       .and. &
                      (n_false_pos == 0_I4P)                   .and. &
                      (n_false_neg == 0_I4P)                   .and. &
                      (max_depth <= AWRAP_MAX_DEPTH)           .and. &
                      (max_leaf_size <= alpha + 1.0e-9_R8P)

! ---- Invariant 2: sphere — same properties on a curved surface.
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)
do k = 1_I4P, N_GRID
   z = -2._R8P + (k - 1) * dx
   do j = 1_I4P, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1_I4P, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1._R8P  ! unit sphere
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)
facets => sphere%facets_ref()
alpha = 0.1_R8P
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call check_node_properties(octree=octree, facets=facets, alpha=alpha, &
                           max_depth=max_depth, max_leaf_size=max_leaf_size, &
                           n_false_pos=n_false_pos, n_false_neg=n_false_neg)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', &
      'inv 2 (sphere): n_nodes=', octree%n_nodes, ' n_leaves=', octree%n_leaves, &
      ' n_boundary=', octree%n_boundary_leaves, ' max_depth=', max_depth, &
      ' n_false_pos=', n_false_pos
print '(A,F8.4,A,I0)', '                max_leaf_size=', max_leaf_size, ' n_false_neg=', n_false_neg
are_tests_passed(2) = (status == AWRAP_STATUS_OK)              .and. &
                      (octree%n_boundary_leaves > 0_I4P)       .and. &
                      (n_false_pos == 0_I4P)                   .and. &
                      (n_false_neg == 0_I4P)                   .and. &
                      (max_depth <= AWRAP_MAX_DEPTH)           .and. &
                      (max_leaf_size <= alpha + 1.0e-9_R8P)

! ---- Invariant 3: degenerate inputs.
!      3a: alpha > bbox diagonal → AWRAP_STATUS_DEGENERATE, 1-leaf octree.
facets => cube%facets_ref()
call awrap_build_octree(facet=facets, alpha=100._R8P, octree=octree, status=status)
print '(A,I0,A,I0,A,I0)', 'inv 3a (alpha > bbox): status=', status, &
      ' n_nodes=', octree%n_nodes, ' n_boundary=', octree%n_boundary_leaves
are_tests_passed(3) = (status == AWRAP_STATUS_DEGENERATE) .and. &
                      (octree%n_nodes == 1_I4P)           .and. &
                      (octree%n_boundary_leaves == 1_I4P)
!      3b: empty facet input → AWRAP_STATUS_BAD_INPUT.
block
   type(facet_object), allocatable, target :: empty_facets(:)
   allocate(empty_facets(0))
   call awrap_build_octree(facet=empty_facets, alpha=0.1_R8P, octree=octree, status=status)
   print '(A,I0)', 'inv 3b (empty input): status=', status
   are_tests_passed(3) = are_tests_passed(3) .and. (status == AWRAP_STATUS_BAD_INPUT)
   deallocate(empty_facets)
endblock

! ---- Invariant 4: cube — flood fill correctness on a closed solid.
facets => cube%facets_ref()
alpha = 0.1_R8P
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
n_in_geometry      = 0_I4P
n_outside_geometry = 0_I4P
n_empty_remaining  = 0_I4P
do i = 1_I4P, octree%n_nodes
   if (octree%node(i)%first_child /= -1_I4P) cycle
   if (octree%node(i)%leaf_flag /= AWRAP_LEAF_FLAG_INSIDE) cycle
   centroid = 0.5_R8P * (octree%node(i)%bmin + octree%node(i)%bmax)
   if (centroid%x > 0._R8P .and. centroid%x < 1._R8P .and. &
       centroid%y > 0._R8P .and. centroid%y < 1._R8P .and. &
       centroid%z > 0._R8P .and. centroid%z < 1._R8P) then
      n_in_geometry = n_in_geometry + 1_I4P
   else
      n_outside_geometry = n_outside_geometry + 1_I4P
   endif
enddo
do i = 1_I4P, octree%n_nodes
   if (octree%node(i)%first_child /= -1_I4P) cycle
   if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_EMPTY) n_empty_remaining = n_empty_remaining + 1_I4P
enddo
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', &
      'inv 4 (cube classify): n_inside=', octree%n_inside_leaves, &
      ' n_outside=', octree%n_outside_leaves, &
      ' n_empty_remaining=', n_empty_remaining, &
      ' INSIDE_in_cube=', n_in_geometry, ' INSIDE_outside_cube=', n_outside_geometry
are_tests_passed(4) = (status == AWRAP_STATUS_OK)                    .and. &
                      (n_empty_remaining == 0_I4P)                   .and. &
                      (octree%n_inside_leaves > 0_I4P)               .and. &
                      (octree%n_outside_leaves > 0_I4P)              .and. &
                      (octree%n_inside_leaves + octree%n_outside_leaves &
                          + octree%n_boundary_leaves == octree%n_leaves) .and. &
                      (n_outside_geometry == 0_I4P)
n_inside_closed_cube = octree%n_inside_leaves   ! save for inv 6 comparison

! ---- Invariant 5: sphere — flood fill on a curved closed surface.
facets => sphere%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
n_in_geometry      = 0_I4P
n_outside_geometry = 0_I4P
do i = 1_I4P, octree%n_nodes
   if (octree%node(i)%first_child /= -1_I4P) cycle
   if (octree%node(i)%leaf_flag /= AWRAP_LEAF_FLAG_INSIDE) cycle
   centroid = 0.5_R8P * (octree%node(i)%bmin + octree%node(i)%bmax)
   if (sqrt(centroid%x**2 + centroid%y**2 + centroid%z**2) < 1._R8P) then
      n_in_geometry = n_in_geometry + 1_I4P
   else
      n_outside_geometry = n_outside_geometry + 1_I4P
   endif
enddo
print '(A,I0,A,I0,A,I0,A,I0)', &
      'inv 5 (sphere classify): n_inside=', octree%n_inside_leaves, &
      ' n_outside=', octree%n_outside_leaves, &
      ' INSIDE_in_sphere=', n_in_geometry, ' INSIDE_outside_sphere=', n_outside_geometry
are_tests_passed(5) = (status == AWRAP_STATUS_OK)                    .and. &
                      (octree%n_inside_leaves > 0_I4P)               .and. &
                      (octree%n_outside_leaves > 0_I4P)              .and. &
                      (n_outside_geometry == 0_I4P)

! ---- Invariant 6: cube minus one facet — expected leak through the hole.
!      Drop the first facet and re-build a surface; flood fill should leak.
facets => cube%facets_ref()  ! re-fetch (block above may have invalidated)
allocate(holed_facets(size(facets, kind=I4P) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call cube_holed%adopt_facets(facets=holed_facets)
facets => cube_holed%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
print '(A,I0,A,I0,A,I0)', &
      'inv 6 (holed cube): n_inside=', octree%n_inside_leaves, &
      ' (closed cube was=', n_inside_closed_cube, &
      ');  n_outside=', octree%n_outside_leaves
are_tests_passed(6) = (status == AWRAP_STATUS_OK) .and. &
                      (octree%n_inside_leaves < n_inside_closed_cube * 30_I4P / 100_I4P)

! ---- Invariant 7: cube wrap surface — watertight, volume reasonable.
facets => cube%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
n_wrap_facets = size(wrap_facets, kind=I4P)
call wrap_surface%adopt_facets(facets=wrap_facets)
wrap_volume = wrap_surface%get_volume()
expected_volume = 1._R8P  ! unit cube
print '(A,I0,A,L1,A,F8.4,A,F8.4)', &
      'inv 7 (cube wrap): n_facets=', n_wrap_facets, ' watertight=', wrap_surface%is_watertight(), &
      ' volume=', wrap_volume, ' expected=', expected_volume
are_tests_passed(7) = (status == AWRAP_STATUS_OK)                   .and. &
                      (n_wrap_facets > 0_I4P)                       .and. &
                      wrap_surface%is_watertight()                  .and. &
                      (abs(wrap_volume - expected_volume) <= 0.30_R8P * expected_volume)

! ---- Invariant 8: sphere wrap surface — watertight, volume reasonable.
facets => sphere%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
n_wrap_facets = size(wrap_facets, kind=I4P)
call wrap_surface%adopt_facets(facets=wrap_facets)
wrap_volume = wrap_surface%get_volume()
expected_volume = 4._R8P / 3._R8P * 3.14159265358979_R8P  ! unit sphere
print '(A,I0,A,L1,A,F8.4,A,F8.4)', &
      'inv 8 (sphere wrap): n_facets=', n_wrap_facets, ' watertight=', wrap_surface%is_watertight(), &
      ' volume=', wrap_volume, ' expected=', expected_volume
are_tests_passed(8) = (status == AWRAP_STATUS_OK)                   .and. &
                      (n_wrap_facets > 0_I4P)                       .and. &
                      wrap_surface%is_watertight()                  .and. &
                      (abs(wrap_volume - expected_volume) <= 0.30_R8P * expected_volume)

! ---- Invariant 9: holed cube wrap — emits something but not closed at this step.
facets => cube_holed%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
n_wrap_facets = size(wrap_facets, kind=I4P)
print '(A,I0,A,A)', 'inv 9 (holed cube wrap): n_facets=', n_wrap_facets, &
      ' (watertightness deferred to step 5 adaptive refinement)'
are_tests_passed(9) = (status == AWRAP_STATUS_OK) .and. (n_wrap_facets > 0_I4P)

! ---- Invariant 10: cube wrap + projection — volume snaps closer to input.
facets => cube%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
call wrap_surface%adopt_facets(facets=wrap_facets)
vol_before = wrap_surface%get_volume()
! Re-extract for projection (adopt_facets consumed wrap_facets via move_alloc).
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
offset = alpha / 3._R8P
call awrap_project_vertices(wrap_facets=wrap_facets, input_facet=facets, input_tree=cube%aabb, &
                            offset=offset, status=status)
call wrap_surface%adopt_facets(facets=wrap_facets)
vol_after = wrap_surface%get_volume()
print '(A,F8.4,A,F8.4,A,L1)', &
      'inv 10 (cube projected): vol_before=', vol_before, ' vol_after=', vol_after, &
      ' watertight=', wrap_surface%is_watertight()
are_tests_passed(10) = (status == AWRAP_STATUS_OK)              .and. &
                       (abs(vol_after - 1._R8P) <= abs(vol_before - 1._R8P)) .and. &
                       wrap_surface%is_watertight()

! ---- Invariant 11: sphere wrap + projection — volume snaps closer to input.
facets => sphere%facets_ref()
call awrap_build_octree(facet=facets, alpha=alpha, octree=octree, status=status)
call awrap_classify_leaves(octree=octree, status=status)
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
call wrap_surface%adopt_facets(facets=wrap_facets)
vol_before = wrap_surface%get_volume()
call awrap_extract_surface(octree=octree, wrapped_facets=wrap_facets, status=status)
call awrap_project_vertices(wrap_facets=wrap_facets, input_facet=facets, input_tree=sphere%aabb, &
                            offset=offset, status=status)
call wrap_surface%adopt_facets(facets=wrap_facets)
vol_after = wrap_surface%get_volume()
expected_volume = 4._R8P / 3._R8P * 3.14159265358979_R8P
print '(A,F8.4,A,F8.4,A,L1)', &
      'inv 11 (sphere projected): vol_before=', vol_before, ' vol_after=', vol_after, &
      ' watertight=', wrap_surface%is_watertight()
are_tests_passed(11) = (status == AWRAP_STATUS_OK)              .and. &
                       (abs(vol_after - expected_volume) <= abs(vol_before - expected_volume)) .and. &
                       wrap_surface%is_watertight()

! ---- Invariant 12: projection preserves manifoldness (vertex dedup works).
n_nm_edges = wrap_surface%get_non_manifold_edges_number()
print '(A,I0)', 'inv 12 (projected manifold): non_manifold_edges=', n_nm_edges
are_tests_passed(12) = (n_nm_edges == 0_I4P)

! ---- Invariant 13: capstone TBP on the cube (issue #50): status OK, closed, manifold, no self-crossing, and its
!      volume is the cube grown by the offset (Minkowski sum with a ball: 1 + 6e + 3 pi e^2 + 4/3 pi e^3).
offset = 0.033_R8P
call cube%alpha_wrap(alpha=0.1_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
are_tests_passed(13) = wrap_ok('inv 13 (cube TBP)', wrapped_tbp, status) .and. &
                       (abs(wrap_volume - expected_volume) <= 0.005_R8P * expected_volume)
print '(A,F8.5,A,F8.5)', '        volume=', wrap_volume, ' offset cube=', expected_volume

! ---- Invariant 14: capstone TBP on the marching-cubes sphere: at least 90% of the wrap vertices on the offset
!      surface, none closer to the input than half the offset (the repair pass moves a few towards their cell faces).
call sphere%alpha_wrap(alpha=0.1_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
are_tests_passed(14) = wrap_ok('inv 14 (sphere TBP)', wrapped_tbp, status) .and. offset_ok(sphere, wrapped_tbp, offset)

! ---- Invariant 15: the holed cube (a hole of inradius 0.29) at alpha = 0.15: the ball gets through the hole, the
!      wrap lines both sides of every wall, two offset layers the barrier keeps apart.
call cube_holed%alpha_wrap(alpha=0.15_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
are_tests_passed(15) = wrap_ok('inv 15 (holed cube, hole open)', wrapped_tbp, status) .and. &
                       (wrap_volume >= 11._R8P * offset) .and. (wrap_volume <= 0.7_R8P)
print '(A,F8.5)', '        volume=', wrap_volume

! ---- Invariant 16: the same holed cube at alpha = 0.6: the hole is capped, volume just below the offset cube's.
offset = 0.02_R8P
call cube_holed%alpha_wrap(alpha=0.6_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
expected_volume = 1._R8P + 6._R8P * offset + 3._R8P * PI * offset**2 + 4._R8P / 3._R8P * PI * offset**3
are_tests_passed(16) = wrap_ok('inv 16 (holed cube, hole capped)', wrapped_tbp, status) .and. &
                       (wrap_volume < expected_volume) .and. (wrap_volume > 0.98_R8P * expected_volume)
print '(A,F8.5,A,F8.5)', '        volume=', wrap_volume, ' offset cube=', expected_volume

! ---- Invariant 17: one triangle (zero thickness, area 1/2): a closed shell around both sides, 2 offset x area + rim.
allocate(plate_facets(1))
plate_facets(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
plate_facets(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P)
call plate_facets(1)%compute_metrix
call plate%adopt_facets(facets=plate_facets)
call plate%alpha_wrap(alpha=0.2_R8P, offset=offset, wrapped=wrapped_tbp, status=status)
wrap_volume = wrapped_tbp%get_volume()
are_tests_passed(17) = wrap_ok('inv 17 (one triangle)', wrapped_tbp, status) .and. &
                       (wrap_volume >= 0.02_R8P) .and. (wrap_volume <= 0.0225_R8P)
print '(A,F8.5)', '        volume=', wrap_volume

! ---- Invariant 18: the measured table of issue #50: cube-01.stl closed and minus one facet, alpha 0.1, 0.3, 0.7,
!      default offset (alpha / 30): every wrap OK, closed, manifold, not self-crossing.
call cube01%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
facets => cube01%facets_ref()
if (allocated(holed_facets)) deallocate(holed_facets) ! adopt_facets may have taken it
allocate(holed_facets(size(facets, kind=I4P) - 1_I4P))
do i = 1_I4P, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1_I4P)
enddo
call cube01_holed%adopt_facets(facets=holed_facets)
are_tests_passed(18) = .true.
do i = 1_I4P, 3_I4P
   alpha = TABLE_ALPHA(i)
   call cube01%alpha_wrap(alpha=alpha, wrapped=wrapped_tbp, status=status)
   are_tests_passed(18) = wrap_ok('inv 18 (cube-01, closed)', wrapped_tbp, status) .and. are_tests_passed(18)
   print '(A,F5.2,A,F8.5)', '        alpha=', alpha, ' volume=', wrapped_tbp%get_volume()
   call cube01_holed%alpha_wrap(alpha=alpha, wrapped=wrapped_tbp, status=status)
   are_tests_passed(18) = wrap_ok('inv 18 (cube-01, holed)', wrapped_tbp, status) .and. are_tests_passed(18)
   print '(A,F5.2,A,F8.5)', '        alpha=', alpha, ' volume=', wrapped_tbp%get_volume()
enddo

! ---- Invariant 19: dragon.stl with the defaults (it was neither watertight nor manifold before issue #50): OK, closed,
!      manifold, not self-crossing, and it encloses the dragon (larger volume).
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%alpha_wrap(wrapped=wrapped_tbp, status=status)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
offset = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2) / 50._R8P / 30._R8P
are_tests_passed(19) = wrap_ok('inv 19 (dragon, defaults)', wrapped_tbp, status) .and. &
                       offset_ok(dragon, wrapped_tbp, offset) .and. (wrapped_tbp%get_volume() > dragon%get_volume())
print '(A,F10.3,A,F10.3)', '        volume=', wrapped_tbp%get_volume(), ' dragon=', dragon%get_volume()

print '(A,19L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. (all(are_tests_passed))) error stop 1

contains

   function wrap_ok(label, w, st) result(ok)
   !< Status OK, watertight, manifold, no self-crossing (find_self_intersections); prints the line.
   !< A reported pair whose facets share a vertex, have normals parallel within 1e-8 rad and do not overlap in their
   !< plane is not counted (printed apart): find_self_intersections skips such pairs only within 1e-12 rad (issue
   !< #46), above it pair_classify's on-plane snap (1e-10) can build a segment that is not there. Checked in exact
   !< rational arithmetic on the 46 such pairs of the sphere wrap and the one of the dragon's (issue #50 night log):
   !< every one meets at the shared vertex only.
   character(*),             intent(in) :: label
   type(surface_stl_object), intent(in) :: w
   integer(I4P),             intent(in) :: st
   logical                              :: ok
   type(intersection_pair_t), allocatable :: pairs(:)
   type(facet_object), pointer          :: fa, fb
   integer(I4P)                         :: p, n_real, n_excused

   call w%find_self_intersections(pairs=pairs)
   n_real = 0 ; n_excused = 0
   do p = 1, size(pairs)
      fa => w%facet_at(pairs(p)%a) ; fb => w%facet_at(pairs(p)%b)
      if (spurious_coplanar(fa, fb)) then
         n_excused = n_excused + 1
      else
         n_real = n_real + 1
      endif
   enddo
   print '(A,A,I0,A,I0,2(A,L1),A,I0,A,I0)', label, ': status=', st, ' n_facets=', w%get_facets_number(), &
         ' watertight=', w%is_watertight(), ' manifold=', w%is_manifold(), ' self-crossings=', n_real, &
         ' near-coplanar touches=', n_excused
   ok = (st == AWRAP_STATUS_OK) .and. w%is_watertight() .and. w%is_manifold() .and. (n_real == 0)
   endfunction wrap_ok

   function offset_ok(input, w, eps) result(ok)
   !< At least 90% of the vertices of w at eps from input (relative 1e-6), none closer than eps / 2 (the repair pass
   !< moves a vertex towards its cell face, which can be closer than eps: 0.85 eps on the sphere, issue #50).
   type(surface_stl_object), intent(in) :: input, w
   real(R8P),                intent(in) :: eps
   logical                              :: ok
   type(facet_object), pointer          :: f
   real(R8P)                            :: d, dmin
   integer(I4P)                         :: i, v, n_on, n_all

   dmin = huge(1._R8P) ; n_on = 0 ; n_all = 0
   do i = 1, w%get_facets_number()
      f => w%facet_at(i)
      do v = 1, 3
         d = input%distance(point=f%vertex(v), is_signed=.false., is_square_root=.true.)
         dmin = min(dmin, d) ; n_all = n_all + 1
         if (abs(d - eps) <= 1.e-6_R8P * eps) n_on = n_on + 1
      enddo
   enddo
   print '(A,F12.9,A,F6.3)', '        min vertex distance / offset=', dmin / eps, ' on the offset surface: ', &
         real(n_on, R8P) / real(n_all, R8P)
   ok = (dmin >= 0.5_R8P * eps) .and. (real(n_on, R8P) >= 0.9_R8P * real(n_all, R8P))
   endfunction offset_ok

   function spurious_coplanar(fa, fb) result(yes)
   !< Two facets sharing a vertex, normals parallel within 1e-8 rad, no overlap in the plane (2-D separating axes
   !< on the edges of either, projected along the dominant normal component).
   type(facet_object), intent(in) :: fa, fb
   logical                        :: yes
   real(R8P)                      :: na(3), nb(3), a(2, 3), b(2, 3)
   integer(I4P)                   :: i, j, shared, drop

   yes = .false.
   shared = 0
   do i = 1, 3
      do j = 1, 3
         if (fa%vertex(i)%x == fb%vertex(j)%x .and. fa%vertex(i)%y == fb%vertex(j)%y .and. &
             fa%vertex(i)%z == fb%vertex(j)%z) shared = shared + 1
      enddo
   enddo
   if (shared /= 1) return
   na = [fa%normal%x, fa%normal%y, fa%normal%z] ; nb = [fb%normal%x, fb%normal%y, fb%normal%z]
   if (norm2([na(2) * nb(3) - na(3) * nb(2), na(3) * nb(1) - na(1) * nb(3), na(1) * nb(2) - na(2) * nb(1)]) > 1.e-8_R8P) return
   drop = maxloc(abs(na), dim=1)
   do i = 1, 3
      a(:, i) = pack([fa%vertex(i)%x, fa%vertex(i)%y, fa%vertex(i)%z], [1, 2, 3] /= drop)
      b(:, i) = pack([fb%vertex(i)%x, fb%vertex(i)%y, fb%vertex(i)%z], [1, 2, 3] /= drop)
   enddo
   yes = separated(a, b) .or. separated(b, a)
   endfunction spurious_coplanar

   pure function separated(p, q) result(sep)
   !< Is there an edge of p with all of q on its outer side or on it (relative 1e-9)?
   real(R8P), intent(in) :: p(2, 3), q(2, 3)
   logical               :: sep
   real(R8P)             :: e(2), sq(3), inner
   integer(I4P)          :: i, j, k, m

   sep = .false.
   do i = 1, 3
      j = mod(i, 3) + 1 ; k = mod(j, 3) + 1
      e = p(:, j) - p(:, i)
      inner = e(1) * (p(2, k) - p(2, i)) - e(2) * (p(1, k) - p(1, i))
      if (inner == 0._R8P) cycle
      do m = 1, 3
         sq(m) = sign(1._R8P, inner) * (e(1) * (q(2, m) - p(2, i)) - e(2) * (q(1, m) - p(1, i)))
      enddo
      if (all(sq <= 1.e-9_R8P * abs(inner))) then
         sep = .true.
         return
      endif
   enddo
   endfunction separated

   subroutine check_node_properties(octree, facets, alpha, max_depth, max_leaf_size, &
                                    n_false_pos, n_false_neg)
   !< For every leaf in `octree`, verify the load-bearing property step 2
   !< depends on: leaf is BOUNDARY iff its bbox overlaps at least one facet.
   !<
   !< `n_false_pos` counts BOUNDARY leaves that overlap NO facet (algorithm
   !< error: should never happen). `n_false_neg` counts INTERIOR leaves that
   !< DO overlap at least one facet (algorithm error: should never happen
   !< for input that hasn't been moved between build and check).
   !< Also returns max recursion depth and max leaf bbox-side reached.
   type(awrap_octree_t), intent(in)  :: octree
   type(facet_object),   intent(in)  :: facets(:)
   real(R8P),            intent(in)  :: alpha
   integer(I4P),         intent(out) :: max_depth
   real(R8P),            intent(out) :: max_leaf_size
   integer(I4P),         intent(out) :: n_false_pos, n_false_neg
   integer(I4P)                      :: i, f, nf
   real(R8P)                         :: sx, sy, sz, max_side
   logical                           :: any_overlap

   max_depth = 0_I4P
   max_leaf_size = 0._R8P
   n_false_pos = 0_I4P
   n_false_neg = 0_I4P
   nf = size(facets, kind=I4P)
   do i = 1_I4P, octree%n_nodes
      if (octree%node(i)%first_child /= -1_I4P) cycle  ! not a leaf
      if (octree%node(i)%depth > max_depth) max_depth = octree%node(i)%depth
      ! Track leaf side for boundary leaves only — these are the ones that should be ≤ α.
      if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_BOUNDARY) then
         sx = octree%node(i)%bmax%x - octree%node(i)%bmin%x
         sy = octree%node(i)%bmax%y - octree%node(i)%bmin%y
         sz = octree%node(i)%bmax%z - octree%node(i)%bmin%z
         max_side = max(sx, sy, sz)
         if (max_side > max_leaf_size) max_leaf_size = max_side
      endif
      ! Cross-check the leaf flag against actual facet overlap.
      any_overlap = .false.
      do f = 1_I4P, nf
         if (triangle_overlaps_aabb(bmin=octree%node(i)%bmin, bmax=octree%node(i)%bmax, &
                                    v1=facets(f)%vertex(1), &
                                    v2=facets(f)%vertex(2), &
                                    v3=facets(f)%vertex(3))) then
            any_overlap = .true.
            exit
         endif
      enddo
      if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_BOUNDARY .and. .not. any_overlap) &
         n_false_pos = n_false_pos + 1_I4P
      if (octree%node(i)%leaf_flag == AWRAP_LEAF_FLAG_INTERIOR .and. any_overlap) &
         n_false_neg = n_false_neg + 1_I4P
   enddo
   endsubroutine check_node_properties

endprogram fossil_val_alpha_wrap
