!< FOSSIL, assert facet connectivity is symmetric and non-manifold edges are
!< explicit, not silently collapsed.
!<
!< Invariant 1 (symmetry): for every facet `f1` and every local edge `e` with
!<   `f1.fcon_edge(e) = f2 > 0`, some edge of `f2` references `f1`.
!<   When this fails, downstream code that walks connectivity must add ad-hoc
!<   guards or `error_stop` — exactly the band-aid we removed by switching to
!<   sort-and-pair edge matching.
!<
!< Invariant 2 (non-manifold detection): on a mesh that contains a deliberate
!<   3-way edge (three triangles sharing the same edge), the surface reports
!<   at least one non-manifold edge and does NOT pick an arbitrary "winner"
!<   for that edge — all three facets have `fcon_edge` = 0 across the shared
!<   edge. trimesh silently drops such edges; we follow Open3D/libigl and
!<   surface them.
!<
!< Run on:
!<   - cube.stl   : clean closed manifold, expect 0 asymmetric, 0 non-manifold.
!<   - dragon.stl : real mesh that historically exposed 58 asymmetric edges.
!<   - synthetic 3-triangle fixture built in-memory.

program fossil_test_connectivity_symmetry

use fossil, only : surface_stl_object, facet_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

logical :: all_passed

all_passed = .true.

call check_mesh_symmetry(file_name='src/tests/cube.stl',   label='cube',   passed=all_passed)
call check_mesh_symmetry(file_name='src/tests/dragon.stl', label='dragon', passed=all_passed)
call check_non_manifold_fixture(passed=all_passed)

print '(A,L1)', 'Are all tests passed? ', all_passed
if (.not. all_passed) error stop 1

contains

   subroutine check_mesh_symmetry(file_name, label, passed)
   !< Load a mesh, run sanitize, then count asymmetric fcon_edge links.
   character(*), intent(in)    :: file_name
   character(*), intent(in)    :: label
   logical,      intent(inout) :: passed
   type(surface_stl_object)    :: surface
   type(facet_object), pointer :: f
   integer(I4P)                :: i, e, neighbour, n_facets, asymmetric
   logical                     :: this_ok

   call surface%load_from_file(file_name=file_name, guess_format=.true.)
   call surface%sanitize

   n_facets = surface%get_facets_number()
   asymmetric = 0
   do i = 1, n_facets
      f => surface%facet_at(i)
      if (.not. associated(f)) cycle
      do e = 1, 3
         neighbour = f%fcon_edge(e)
         if (neighbour <= 0) cycle
         if (.not. has_back_link(surface, neighbour, i)) asymmetric = asymmetric + 1
      enddo
   enddo

   this_ok = (asymmetric == 0)
   passed = passed .and. this_ok
   print '(A,A,A,I0,A,I0,A,L1)', '  ', label, ': facets=', n_facets, &
       ' asymmetric_edges=', asymmetric, ' passed=', this_ok
   endsubroutine check_mesh_symmetry

   function has_back_link(surf, neighbour_id, self_id) result(yes)
   !< Look up `neighbour_id`'s edges and check that one of them references `self_id`.
   type(surface_stl_object),     intent(in) :: surf
   integer(I4P),                 intent(in) :: neighbour_id
   integer(I4P),                 intent(in) :: self_id
   logical                                  :: yes
   type(facet_object), pointer              :: g
   integer(I4P)                             :: k

   yes = .false.
   g => surf%facet_at(neighbour_id)
   if (.not. associated(g)) return
   do k = 1, 3
      if (g%fcon_edge(k) == self_id) then
         yes = .true.
         return
      endif
   enddo
   endfunction has_back_link

   subroutine check_non_manifold_fixture(passed)
   !< Build a deliberate non-manifold mesh in-memory: three triangles share a
   !< single edge along the z-axis. The shared edge is V1=(0,0,0)-V2=(0,0,1).
   !< Triangles fan out at +x, at +y rotated 120 degrees, and at -x rotated 240.
   !< All three triangles list the shared edge as their first edge (vertices 1-2).
   !<
   !< Correct connectivity: the shared edge has multiplicity 3 -> non-manifold.
   !< Every facet's fcon_edge(1) (the shared edge) should be 0 (not linked).
   !< The surface must report at least one non-manifold edge.
   logical, intent(inout) :: passed
   type(surface_stl_object) :: surface
   type(facet_object), allocatable :: facets(:)
   type(facet_object), pointer     :: f
   integer(I4P)                    :: i, link_count, n_nm
   logical                         :: this_ok

   allocate(facets(3))
   ! shared edge V1=(0,0,0) to V2=(0,0,1) is vertex(1)-vertex(2) on each facet
   facets(1)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
   facets(1)%vertex(2) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 1._R8P * ez_R8P
   facets(1)%vertex(3) = 1._R8P * ex_R8P + 0._R8P * ey_R8P + 0.5_R8P * ez_R8P

   facets(2)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
   facets(2)%vertex(2) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 1._R8P * ez_R8P
   facets(2)%vertex(3) = -0.5_R8P * ex_R8P + 0.866_R8P * ey_R8P + 0.5_R8P * ez_R8P

   facets(3)%vertex(1) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 0._R8P * ez_R8P
   facets(3)%vertex(2) = 0._R8P * ex_R8P + 0._R8P * ey_R8P + 1._R8P * ez_R8P
   facets(3)%vertex(3) = -0.5_R8P * ex_R8P - 0.866_R8P * ey_R8P + 0.5_R8P * ez_R8P

   call surface%adopt_facets(facets)

   link_count = 0
   do i = 1, 3
      f => surface%facet_at(i)
      if (.not. associated(f)) cycle
      if (f%fcon_edge(1) > 0) link_count = link_count + 1
   enddo
   n_nm = surface%get_non_manifold_edges_number()

   this_ok = (link_count == 0) .and. (n_nm >= 1)
   passed = passed .and. this_ok
   print '(A,I0,A,I0,A,L1)', '  3-facet fan: linked_across_shared_edge=', link_count, &
       ' non_manifold_edges=', n_nm, ' passed=', this_ok
   endsubroutine check_non_manifold_fixture

endprogram fossil_test_connectivity_symmetry
