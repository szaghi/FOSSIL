!< FOSSIL test: optional status argument on mutating procedures (audit #14 S1).
!<
!< Verifies that each procedure that can fail returns a non-zero status
!< instead of error-stopping when the status argument is present.

program fossil_test_status_arg

use fossil, only : surface_stl_object, &
                   STATUS_OK, STATUS_ALLOC_FAIL, STATUS_AMBIGUOUS_ARGS, &
                   STATUS_FILE_NOT_FOUND, STATUS_FILE_OPEN_FAIL
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: surface
integer(I4P)             :: st
logical                  :: are_tests_passed(4)

are_tests_passed = .false.

! 1. load_from_file: missing file -> STATUS_FILE_NOT_FOUND
call surface%load_from_file(file_name='no_such_file.stl', guess_format=.true., status=st)
are_tests_passed(1) = (st == STATUS_FILE_NOT_FOUND)

! 2. translate: ambiguous args -> STATUS_AMBIGUOUS_ARGS; load a real surface first
call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true., status=st)
call surface%translate(x=1._R8P, delta=vector_R8P(1._R8P, 0._R8P, 0._R8P), status=st)
are_tests_passed(2) = (st == STATUS_AMBIGUOUS_ARGS)

! 3. resize: ambiguous args -> STATUS_AMBIGUOUS_ARGS
call surface%resize(x=2._R8P, factor=vector_R8P(2._R8P, 2._R8P, 2._R8P), status=st)
are_tests_passed(3) = (st == STATUS_AMBIGUOUS_ARGS)

! 4. save_into_file: unwritable path -> STATUS_FILE_OPEN_FAIL
call surface%save_into_file(file_name='/no_permission/out.stl', status=st)
are_tests_passed(4) = (st == STATUS_FILE_OPEN_FAIL)

print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_status_arg
