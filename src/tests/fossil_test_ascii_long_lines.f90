!< FOSSIL test: ASCII STL lines longer than 80 characters (issue #66).
!<
!< A vertex line with three coordinates of 17 significant digits ('(A,3ES26.17)', 90 characters) must read back the
!< exact doubles. The line buffer was 80 characters: the last coordinate lost its exponent ('2.5E-01' read as 2.5),
!< silently.
program fossil_test_ascii_long_lines

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P

implicit none

character(*), parameter     :: FN = 'fossil_test_ascii_long_lines.stl'
real(R8P),    parameter     :: X(3, 3) = reshape([0.1_R8P, 0.25_R8P, 0.75_R8P, 1.3_R8P, 2.5e-7_R8P, -0.123456789012345_R8P, &
                                                  0.2_R8P, 1.75_R8P, 3.e-3_R8P], [3, 3])
type(surface_stl_object)    :: s
type(facet_object), pointer :: f
integer(I4P)                :: u, k
logical                     :: ok

open(newunit=u, file=FN, action='write', status='replace')
write(u, '(A)') 'solid long'
write(u, '(A)') '  facet normal 0 0 0'
write(u, '(A)') '    outer loop'
do k = 1, 3
   write(u, '(A,3ES26.17)') '      vertex', X(:, k)
enddo
write(u, '(A)') '    endloop'
write(u, '(A)') '  endfacet'
write(u, '(A)') 'endsolid long'
close(u)
call s%load_from_file(file_name=FN, guess_format=.true.)
f => s%facet_at(1)
ok = s%get_facets_number() == 1
do k = 1, 3
   ok = ok .and. f%vertex(k)%x == X(1, k) .and. f%vertex(k)%y == X(2, k) .and. f%vertex(k)%z == X(3, k)
enddo
open(newunit=u, file=FN) ; close(u, status='delete')
print '(A,L1)', 'exact coordinates from 90-character lines: ', ok
print '(A,L1)', 'Are all tests passed? ', ok
endprogram fossil_test_ascii_long_lines
