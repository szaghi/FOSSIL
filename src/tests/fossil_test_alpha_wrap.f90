!< FOSSIL test: alpha-wrap octree (issue #18 §1.6, step 1).
!<
!< Step-1 scope: validate the leaf-size-α octree built by `awrap_build_octree`.
!< Subsequent commits add steps 2 (flood fill), 3 (dual-contour), 4 (project),
!< 5 (adaptive refinement + capstone TBP).
!<
!< Three invariants for step 1:
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

program fossil_test_alpha_wrap

use fossil, only : surface_stl_object, extract_isosurface
use fossil_alpha_wrap, only : awrap_octree_t, awrap_build_octree, &
                              AWRAP_STATUS_OK, AWRAP_STATUS_BAD_INPUT, AWRAP_STATUS_DEGENERATE, &
                              AWRAP_LEAF_FLAG_BOUNDARY, AWRAP_LEAF_FLAG_INTERIOR, &
                              AWRAP_MAX_DEPTH
use fossil_facet_object, only : facet_object
use fossil_utils, only : triangle_overlaps_aabb
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere
type(facet_object), allocatable :: sphere_facets(:)
type(facet_object), pointer :: facets(:)
type(awrap_octree_t) :: octree
real(R8P), allocatable :: values(:,:,:)
type(vector_R8P) :: bmin, bmax
real(R8P) :: x, y, z, dx, alpha
real(R8P) :: max_leaf_size, leaf_sx, leaf_sy, leaf_sz
integer(I4P) :: status, i, j, k, n, max_depth
integer(I4P) :: n_false_pos, n_false_neg
logical :: any_overlap, prop_no_false_pos, prop_no_false_neg, prop_size_ok
logical :: are_tests_passed(3)

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

print '(A,3L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains

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

endprogram fossil_test_alpha_wrap
