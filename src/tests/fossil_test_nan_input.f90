!< FOSSIL, assert load_from_file refuses NaN/Inf coordinates via STATUS_INVALID_INPUT.
!<
!< Without this guard, NaN coordinates propagate into facet normals, AABB extents,
!< pseudo-normal sign, and every downstream geometric query. Detecting at load time
!< gives the caller a chance to fail early with a clear status code, rather than
!< chase NaNs through a fully-populated surface state.
!<
!< Strategy: write a tiny ASCII STL containing a `nan` vertex coordinate (the ASCII
!< STL format allows arbitrary float-formatted strings, and gfortran's list-directed
!< read accepts 'nan'/'NaN' / 'inf' / 'Inf' as IEEE values). Try to load it; the
!< status must be STATUS_INVALID_INPUT and the surface must remain empty.

program fossil_test_nan_input

use fossil, only : surface_stl_object, STATUS_INVALID_INPUT, STATUS_OK
use penf, only : I4P

implicit none

character(*), parameter  :: TMP_PATH = '/tmp/fossil_test_nan_input.stl'
type(surface_stl_object) :: surface
integer(I4P)             :: file_unit, st
logical                  :: tests_passed(2)

! Write a minimal one-facet ASCII STL with a NaN in vertex 2's y coordinate.
open(newunit=file_unit, file=TMP_PATH, status='replace', action='write')
write(file_unit, '(A)') 'solid bad'
write(file_unit, '(A)') '  facet normal 0.0 0.0 1.0'
write(file_unit, '(A)') '    outer loop'
write(file_unit, '(A)') '      vertex 0.0 0.0 0.0'
write(file_unit, '(A)') '      vertex 1.0 nan 0.0'
write(file_unit, '(A)') '      vertex 0.0 1.0 0.0'
write(file_unit, '(A)') '    endloop'
write(file_unit, '(A)') '  endfacet'
write(file_unit, '(A)') 'endsolid bad'
close(file_unit)

! Attempt to load — must yield STATUS_INVALID_INPUT, surface stays empty.
call surface%load_from_file(file_name=TMP_PATH, is_ascii=.true., status=st)
tests_passed(1) = (st == STATUS_INVALID_INPUT)
tests_passed(2) = (surface%get_facets_number() == 0)

print '(A,I0)',     'status returned:       ', st
print '(A,I0)',     'facets after load:     ', surface%get_facets_number()
print '(A,L1)', 'Are all tests passed? ', all(tests_passed)

! Tidy up the temp file (best-effort).
open(newunit=file_unit, file=TMP_PATH, status='old', action='read', iostat=st)
if (st == 0) close(file_unit, status='delete')

if (.not. all(tests_passed)) error stop 1

endprogram fossil_test_nan_input
