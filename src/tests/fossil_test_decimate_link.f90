!< FOSSIL test: decimate keeps a valid surface (issue #81).
!<
!<   1. cube-01 (48 facets, flat faces, watertight) to 36 facets: before #81 the result had 2 non-manifold edges and
!<      2 shells with DEC_STATUS_OK: a collapse left a zero-area facet (two pool vertices at one place, merged into a
!<      pinch by the re-analysis), which the normal-flip check let through. Now the result is watertight, one shell,
!<      the unit volume kept, with DEC_STATUS_OK (or NO_PROGRESS if the target is unreachable).
!<   2. cube-01 to as few facets as it goes (target 4): watertight and manifold whatever the status.
!<   3. target_facets <= 0: DEC_STATUS_BAD_INPUT, the surface untouched.
!<   4. An open box (cube-01 without its top face; boundary variant of the link condition): manifold, one shell, its
!<      boundary a single loop.
!<   5. A tetrahedron to 2 facets: every edge fails the link condition (the edge (c, d) opposite to (a, b) is in both
!<      links), the only check that sees it; before #81 the collapse left two facets on the same three vertices.
!<      Expected: NO_PROGRESS, the tetrahedron unchanged.

program fossil_test_decimate_link

use fossil, only : surface_stl_object, facet_object, DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS, DEC_STATUS_BAD_INPUT
use penf,   only : I4P, R8P

implicit none

type(surface_stl_object)        :: cube, box, tet
type(facet_object), allocatable :: f(:)
integer(I4P)                    :: status, k, n
logical                         :: are_tests_passed(5)

are_tests_passed = .false.

! ----- 1. cube-01 to 36 facets.
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call cube%decimate(target_facets=36_I4P, status=status)
call report('cube-01 to 36', cube, status)
are_tests_passed(1) = (status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. valid_closed(cube) .and. &
                      abs(cube%get_volume() - 1._R8P) <= 1.e-12_R8P
if (status == DEC_STATUS_OK) are_tests_passed(1) = are_tests_passed(1) .and. cube%get_facets_number() <= 36_I4P

! ----- 2. cube-01 to 4 facets: unreachable, the surface stays valid.
call cube%destroy
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call cube%decimate(target_facets=4_I4P, status=status)
call report('cube-01 to 4', cube, status)
are_tests_passed(2) = (status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. valid_closed(cube)

! ----- 3. target_facets <= 0.
call cube%destroy
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call cube%decimate(target_facets=0_I4P, status=status)
are_tests_passed(3) = status == DEC_STATUS_BAD_INPUT .and. cube%get_facets_number() == 48_I4P
call cube%decimate(target_facets=-3_I4P, status=status)
are_tests_passed(3) = are_tests_passed(3) .and. status == DEC_STATUS_BAD_INPUT .and. cube%get_facets_number() == 48_I4P
print '(A,I0,A,I0)', 'target 0 and -3: status ', status, ', facets ', cube%get_facets_number()

! ----- 4. An open box: the facets of cube-01 off its top face z = 1.
call cube%destroy
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
n = 0
allocate(f(cube%get_facets_number()))
do k = 1, cube%get_facets_number()
   associate(fk => cube%facet_at(k))
      if (fk%vertex(1)%z > 0.999_R8P .and. fk%vertex(2)%z > 0.999_R8P .and. fk%vertex(3)%z > 0.999_R8P) cycle
      n = n + 1 ; f(n) = fk
   endassociate
enddo
f = f(1:n)
call box%adopt_facets(facets=f)
print '(A,I0,A,L1)', 'open box: ', box%get_facets_number(), ' facets, manifold ', box%is_manifold()
call box%decimate(target_facets=8_I4P, status=status)
call report('open box to 8', box, status)
are_tests_passed(4) = (status == DEC_STATUS_OK .or. status == DEC_STATUS_NO_PROGRESS) .and. box%is_manifold() .and. &
                      .not. box%is_watertight() .and. box%get_shells_number() == 1_I4P
if (status == DEC_STATUS_OK) are_tests_passed(4) = are_tests_passed(4) .and. box%get_facets_number() <= 8_I4P

! ----- 5. A tetrahedron to 2 facets.
if (allocated(f)) deallocate(f)
allocate(f(4))
call tri(f(1), [0._R8P, 0._R8P, 0._R8P], [0._R8P, 1._R8P, 0._R8P], [1._R8P, 0._R8P, 0._R8P])
call tri(f(2), [0._R8P, 0._R8P, 0._R8P], [1._R8P, 0._R8P, 0._R8P], [0._R8P, 0._R8P, 1._R8P])
call tri(f(3), [0._R8P, 0._R8P, 0._R8P], [0._R8P, 0._R8P, 1._R8P], [0._R8P, 1._R8P, 0._R8P])
call tri(f(4), [1._R8P, 0._R8P, 0._R8P], [0._R8P, 1._R8P, 0._R8P], [0._R8P, 0._R8P, 1._R8P])
call tet%adopt_facets(facets=f)
call tet%decimate(target_facets=2_I4P, status=status)
call report('tetrahedron to 2', tet, status)
are_tests_passed(5) = status == DEC_STATUS_NO_PROGRESS .and. tet%get_facets_number() == 4_I4P .and. valid_closed(tet) &
                      .and. abs(tet%get_volume() - 1._R8P / 6._R8P) <= 1.e-15_R8P

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1)%x = p1(1) ; f%vertex(1)%y = p1(2) ; f%vertex(1)%z = p1(3)
   f%vertex(2)%x = p2(1) ; f%vertex(2)%y = p2(2) ; f%vertex(2)%z = p2(3)
   f%vertex(3)%x = p3(1) ; f%vertex(3)%y = p3(2) ; f%vertex(3)%z = p3(3)
   call f%compute_metrix
   endsubroutine tri

   function valid_closed(s) result(yes)
   !< Watertight, manifold, one shell.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes

   yes = s%is_watertight() .and. s%is_manifold() .and. s%get_shells_number() == 1_I4P
   endfunction valid_closed

   subroutine report(what, s, status)
   !< One line per case.
   character(*),             intent(in) :: what
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: status

   print '(2A,I0,A,I0,A,F16.12,A,I0,A,I0,A,L1)', what, ': status ', status, ', facets ', s%get_facets_number(), &
         ', volume ', s%get_volume(), ', non-manifold edges ', s%get_non_manifold_edges_number(), ', shells ', &
         s%get_shells_number(), ', watertight ', s%is_watertight()
   endsubroutine report
endprogram fossil_test_decimate_link
