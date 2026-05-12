!< FOSSIL, assert AABB_AUTO_REFINEMENT picks sensible depths per facet count.
!<
!< Confirms the heuristic produces levels that match the empirical sweet spot
!< (within the [AUTO_MIN_LEVELS, AUTO_MAX_LEVELS] range) and that the resulting
!< distance queries remain bit-exact against brute force.
!<
!< Reference numbers from the benchmark sweep (see commit log for step 5.A):
!<   cube    (12 facets)    -> level 1 from heuristic, may be trimmed to 0 by the
!<                              "largest_edge_len > octant median" geometric safety
!<                              check in aabb_tree%initialize (legitimate; a cube
!<                              spans the whole bbox so deeper octree boxes are
!<                              smaller than facet edges).
!<   naca    (188 facets)   -> level 1 from heuristic.
!<   dragon  (6588 facets)  -> level 3 from heuristic, geometry-trimmed to 2.

program fossil_test_auto_refinement

use fossil, only : surface_stl_object, AABB_AUTO_REFINEMENT, AABB_TREE_OCTREE
use fossil_aabb_tree_object, only : AABB_USE_INDEX, AABB_USE_BRUTE_FORCE
use penf, only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

logical :: all_passed

all_passed = .true.

call check_one('src/tests/cube.stl',           'cube',   expected_min=0, expected_max=1, passed=all_passed)
call check_one('src/tests/naca0012-ascii.stl', 'naca',   expected_min=0, expected_max=2, passed=all_passed)
call check_one('src/tests/dragon.stl',         'dragon', expected_min=1, expected_max=4, passed=all_passed)

print '(A,L1)', 'Are all tests passed? ', all_passed
if (.not. all_passed) error stop 1

contains

   subroutine check_one(file_name, label, expected_min, expected_max, passed)
   !< Load with AABB_AUTO_REFINEMENT, confirm the chosen depth falls in the expected
   !< range, then run a single distance query both with the tree and brute force.
   !< Distances must match bit-exactly (squared d^2).
   character(*), intent(in)    :: file_name
   character(*), intent(in)    :: label
   integer(I4P), intent(in)    :: expected_min, expected_max
   logical,      intent(inout) :: passed
   type(surface_stl_object)    :: surface
   type(vector_R8P)            :: p
   real(R8P)                   :: d_tree, d_brute
   integer(I4P)                :: levels
   logical                     :: this_ok

   ! Explicitly request the octree tree kind — this test exercises the octree's
   ! auto-tune heuristic. Since the library default is now AABB_TREE_SAH_BVH,
   ! refinement_levels is meaningless unless we opt back into the octree.
   call surface%load_from_file(file_name=file_name, guess_format=.true.,    &
                               aabb_refinement_levels=AABB_AUTO_REFINEMENT, &
                               aabb_tree_kind=AABB_TREE_OCTREE)
   call surface%sanitize

   levels = surface%aabb%get_refinement_levels()

   associate(bmin => surface%get_bmin(), bmax => surface%get_bmax())
      p = (bmin + bmax) * 0.5_R8P + 0.1_R8P * (bmax - bmin)
   end associate

   call surface%aabb%set_use_index(AABB_USE_BRUTE_FORCE)
   d_brute = surface%distance(point=p, is_signed=.false.)
   call surface%aabb%set_use_index(AABB_USE_INDEX)
   d_tree  = surface%distance(point=p, is_signed=.false.)

   this_ok = (levels >= expected_min) .and. (levels <= expected_max) .and. (d_tree == d_brute)
   passed  = passed .and. this_ok
   print '(A,A,A,I0,A,I0,A,I0,A,L1)', '  ', label, ': auto_levels=', levels, &
       '  (expected [', expected_min, ',', expected_max, '])  passed=', this_ok
   endsubroutine check_one

endprogram fossil_test_auto_refinement
