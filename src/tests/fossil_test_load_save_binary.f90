!< FOSSIL, test load and write binary file.

program fossil_test_load_save_binary
!< FOSSIL, test load and write binary file.

use fossil, only : file_stl_object
use penf, only : R4P

implicit none

type(file_stl_object) :: file_1              !< STL file.
type(file_stl_object) :: file_2              !< STL file.
integer               :: file_unit           !< File unit.
logical               :: are_tests_passed(8) !< Result of tests check.

are_tests_passed = .false.

call file_1%initialize(file_name='src/tests/naca0012-binary.stl', is_ascii=.false.)
call file_1%load_from_file
are_tests_passed(1) = file_1%facets_number == 188
are_tests_passed(2) = real(file_1%facet(5)%vertex(2)%x, R4P) == 0.683601_R4P
are_tests_passed(3) = real(file_1%facet(5)%vertex(2)%y, R4P) == -0.00763869_R4P
are_tests_passed(4) = real(file_1%facet(5)%vertex(2)%z, R4P) == 0._R4P
call file_1%save_into_file(file_name='fossil_test_load_save-naca0012-binary.stl')

call file_2%initialize(file_name='fossil_test_load_save-naca0012-binary.stl', is_ascii=.false.)
call file_2%load_from_file
are_tests_passed(5) = file_2%facets_number == 188
are_tests_passed(6) = real(file_2%facet(5)%vertex(2)%x, R4P) == 0.683601_R4P
are_tests_passed(7) = real(file_2%facet(5)%vertex(2)%y, R4P) == -0.00763869_R4P
are_tests_passed(8) = real(file_2%facet(5)%vertex(2)%z, R4P) == 0._R4P

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

open(newunit=file_unit, file='fossil_test_load_save-naca0012-binary.stl')
close(unit=file_unit, status='delete')
endprogram fossil_test_load_save_binary
