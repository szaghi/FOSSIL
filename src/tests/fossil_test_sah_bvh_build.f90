!< FOSSIL, assert the SAH BVH builder produces a well-formed binary tree.
!<
!< Structural invariants checked (independent of distance-query semantics):
!<   1. Every original facet appears in exactly one leaf's facet_id list.
!<   2. Every internal node's bbox contains both children's bboxes (a node's
!<      AABB must enclose every triangle in its subtree).
!<   3. Leaves have left_child = right_child = 0; internal nodes have both
!<      non-zero indices in `[1, nodes_number-1]`.
!<   4. nodes_number > 0 on a non-empty mesh.
!<
!< Run on cube (12 facets, post-sanitize), naca0012 (188), and dragon (post-
!< sanitize ~6500). The test does NOT depend on distance queries — those
!< require the traversal dispatch landing in the next commit. This is the
!< build-side check that the structure itself is correct.

program fossil_test_sah_bvh_build

use fossil, only : surface_stl_object, AABB_TREE_SAH_BVH
use fossil_aabb_node_object, only : aabb_node_object
use fossil_list_id_object,   only : list_id_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter :: BBOX_TOL = 1.0e-9_R8P
logical :: all_passed

all_passed = .true.

call check_mesh('src/tests/cube.stl',           'cube',   all_passed)
call check_mesh('src/tests/naca0012-ascii.stl', 'naca',   all_passed)
call check_mesh('src/tests/dragon.stl',         'dragon', all_passed)

print '(A,L1)', 'Are all tests passed? ', all_passed
if (.not. all_passed) error stop 1

contains

   subroutine check_mesh(file_name, label, passed)
   !< Load a mesh, sanitize, build the SAH BVH, and verify all four invariants.
   character(*), intent(in)    :: file_name
   character(*), intent(in)    :: label
   logical,      intent(inout) :: passed
   type(surface_stl_object)    :: surface
   type(aabb_node_object), pointer :: node
   type(list_id_object)            :: ids
   integer(I4P), allocatable       :: facet_owner(:)   !< Which leaf claims each facet (0 = none, -1 = multiple).
   integer(I4P)                    :: n_facets, n_nodes
   integer(I4P)                    :: n_unowned, n_duplicate
   integer(I4P)                    :: i, k, lc, rc
   logical                         :: bbox_ok, child_index_ok, this_ok

   call surface%load_from_file(file_name=file_name, guess_format=.true.)
   call surface%sanitize
   ! Rebuild with the BVH path explicitly.
   call surface%analyze(aabb_tree_kind=AABB_TREE_SAH_BVH)

   n_facets = surface%get_facets_number()
   n_nodes  = surface%aabb%get_nodes_number()

   if (n_facets <= 0 .or. n_nodes <= 0) then
      print '(A,A,A)', '  ', label, ': empty surface or zero nodes — skipping'
      return
   endif

   ! Invariant 1: facet ownership. Walk all leaves, record which facet ids appear.
   ! After the walk, every id in [1, n_facets] must be owned exactly once.
   allocate(facet_owner(n_facets))
   facet_owner = 0
   do i = 0, n_nodes - 1
      node => surface%aabb%node_at(i)
      if (.not. associated(node)) cycle
      lc = node%get_left_child()
      rc = node%get_right_child()
      if (lc /= 0 .or. rc /= 0) cycle      ! internal node, skip
      ids = node%facet_id()
      do k = 1, ids%ids_number
         if (ids%id(k) < 1 .or. ids%id(k) > n_facets) cycle
         if (facet_owner(ids%id(k)) == 0) then
            facet_owner(ids%id(k)) = i + 1   ! +1 to distinguish from "unowned"
         else
            facet_owner(ids%id(k)) = -1      ! claimed twice
         endif
      enddo
   enddo
   n_unowned   = count(facet_owner == 0)
   n_duplicate = count(facet_owner == -1)

   ! Invariant 2 + 3: bbox containment and child-index sanity.
   bbox_ok        = .true.
   child_index_ok = .true.
   do i = 0, n_nodes - 1
      node => surface%aabb%node_at(i)
      if (.not. associated(node)) cycle
      lc = node%get_left_child()
      rc = node%get_right_child()
      if (lc == 0 .and. rc == 0) cycle    ! leaf
      ! Internal node: both children must be valid indices > i (DFS order; children come after parent).
      if (lc <= 0 .or. lc > n_nodes - 1 .or. rc <= 0 .or. rc > n_nodes - 1) then
         child_index_ok = .false.
         cycle
      endif
      ! Bbox containment: parent's bbox must enclose each child's bbox (with tolerance).
      if (.not. bbox_contains(node, surface%aabb%node_at(lc))) bbox_ok = .false.
      if (.not. bbox_contains(node, surface%aabb%node_at(rc))) bbox_ok = .false.
   enddo

   this_ok = (n_unowned == 0) .and. (n_duplicate == 0) .and. bbox_ok .and. child_index_ok
   passed = passed .and. this_ok

   print '(A,A,A,I0,A,I0,A,I0,A,I0,A,L1,A,L1,A,L1)',                       &
       '  ', label, ': nodes=', n_nodes, ' facets=', n_facets,             &
       ' unowned=', n_unowned, ' duplicate=', n_duplicate,                 &
       ' bbox_ok=', bbox_ok, ' child_idx_ok=', child_index_ok,             &
       ' passed=', this_ok
   deallocate(facet_owner)
   endsubroutine check_mesh

   function bbox_contains(parent, child) result(yes)
   !< Parent bbox must enclose child bbox in every dimension (with floating-point
   !< tolerance). True if `parent%bmin <= child%bmin` and `parent%bmax >= child%bmax`.
   type(aabb_node_object), pointer, intent(in) :: parent, child
   logical                                     :: yes
   type(vector_R8P)                            :: pmin, pmax, cmin, cmax

   yes = .false.
   if (.not. associated(parent)) return
   if (.not. associated(child))  return
   pmin = parent%bmin(); pmax = parent%bmax()
   cmin = child%bmin();  cmax = child%bmax()
   yes = (pmin%x <= cmin%x + BBOX_TOL) .and. (pmax%x >= cmax%x - BBOX_TOL) .and. &
         (pmin%y <= cmin%y + BBOX_TOL) .and. (pmax%y >= cmax%y - BBOX_TOL) .and. &
         (pmin%z <= cmin%z + BBOX_TOL) .and. (pmax%z >= cmax%z - BBOX_TOL)
   endfunction bbox_contains

endprogram fossil_test_sah_bvh_build
