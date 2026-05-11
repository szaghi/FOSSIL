!< FOSSIL, test load and write binary file.

program fossil_test_load_save_binary
!< FOSSIL, test load and write binary file.

use fossil, only : surface_stl_object
use penf, only : R4P

implicit none

type(surface_stl_object) :: surface_1           !< STL surface.
type(surface_stl_object) :: surface_2           !< STL surface.
integer                  :: file_unit           !< File unit.
logical                  :: are_tests_passed(8) !< Result of tests check.

are_tests_passed = .false.

call surface_1%load_from_file(file_name='src/tests/naca0012-binary.stl', is_ascii=.false.)
are_tests_passed(1) = surface_1%get_facets_number() == 188
associate (f5 => surface_1%facet_at(5))
   are_tests_passed(2) = real(f5%vertex(2)%x, R4P) == 0.683601_R4P
   are_tests_passed(3) = real(f5%vertex(2)%y, R4P) == -0.00763869_R4P
   are_tests_passed(4) = real(f5%vertex(2)%z, R4P) == 0._R4P
end associate
call surface_1%save_into_file(file_name='fossil_test_load_save-naca0012-binary.stl', is_ascii=.false.)

call surface_2%load_from_file(file_name='fossil_test_load_save-naca0012-binary.stl', is_ascii=.false.)
are_tests_passed(5) = surface_2%get_facets_number() == 188
associate (f5 => surface_2%facet_at(5))
   are_tests_passed(6) = real(f5%vertex(2)%x, R4P) == 0.683601_R4P
   are_tests_passed(7) = real(f5%vertex(2)%y, R4P) == -0.00763869_R4P
   are_tests_passed(8) = real(f5%vertex(2)%z, R4P) == 0._R4P
end associate

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

open(newunit=file_unit, file='fossil_test_load_save-naca0012-binary.stl')
close(unit=file_unit, status='delete')
endprogram fossil_test_load_save_binary
