!< FOSSIL, test load and write ASCII file.

program fossil_test_load_save_ascii
!< FOSSIL, test load and write ASCII file.

use fossil, only : file_stl_object, surface_stl_object
use penf, only : R8P

implicit none

type(file_stl_object)    :: file_1              !< STL file.
type(surface_stl_object) :: surface_1           !< STL surface.
type(file_stl_object)    :: file_2              !< STL file.
type(surface_stl_object) :: surface_2           !< STL surface.
integer                  :: file_unit           !< File unit.
logical                  :: are_tests_passed(8) !< Result of tests check.

are_tests_passed = .false.

call file_1%load_from_file(facet=surface_1%facet, file_name='src/tests/naca0012-ascii.stl', is_ascii=.true.)
call surface_1%analize
are_tests_passed(1) = surface_1%facets_number == 188
are_tests_passed(2) = surface_1%facet(5)%vertex(2)%x == 0.683601_R8P
are_tests_passed(3) = surface_1%facet(5)%vertex(2)%y == -0.00763869_R8P
are_tests_passed(4) = surface_1%facet(5)%vertex(2)%z == 0._R8P

call file_1%save_into_file(facet=surface_1%facet, file_name='fossil_test_load_save-naca0012-ascii.stl')

call file_2%load_from_file(facet=surface_2%facet, file_name='fossil_test_load_save-naca0012-ascii.stl', is_ascii=.true.)
call surface_2%analize
are_tests_passed(5) = surface_2%facets_number == 188
are_tests_passed(6) = surface_2%facet(5)%vertex(2)%x == 0.683601_R8P
are_tests_passed(7) = surface_2%facet(5)%vertex(2)%y == -0.00763869_R8P
are_tests_passed(8) = surface_2%facet(5)%vertex(2)%z == 0._R8P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

open(newunit=file_unit, file='fossil_test_load_save-naca0012-ascii.stl')
close(unit=file_unit, status='delete')
endprogram fossil_test_load_save_ascii
