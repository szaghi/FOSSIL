!< FOSSIL test: fossilizer --compare (issue #60, roadmap #58 WP2, step 6).
!<
!< `fossilizer -i cube-02.stl --compare <cube-02 translated by 0.3 along x> --tolerance 1e-4`: both one-sided
!< deviations and the Hausdorff distance printed as certified brackets of 0.3, read back from the lines.
!< Runs exe/fossilizer from the repository root, like every FOSSIL test.
program fossil_test_deviation_cli

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P

implicit none

character(*), parameter  :: MOVED = 'fossil_test_deviation_cli_moved.stl'
character(*), parameter  :: OUT   = 'fossil_test_deviation_cli.txt'
type(surface_stl_object) :: cube
character(256)           :: line
real(R8P)                :: b(2, 3)
integer(I4P)             :: u, exitstat, ios, k
logical                  :: ok

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call cube%translate(delta=0.3_R8P * ex_R8P)
call cube%analyze
call cube%save_into_file(file_name=MOVED)
call execute_command_line('exe/fossilizer -i src/tests/cube-02.stl --compare '//MOVED//' --tolerance 1e-4 > '//OUT, &
                          exitstat=exitstat)
b = -1._R8P ; k = 0
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, 'lower, upper') > 0) then
      k = k + 1
      read(line(index(line, 'upper') + 6:), *) b(:, k)
   endif
enddo
close(u, status='delete')
open(newunit=u, file=MOVED, status='old') ; close(u, status='delete')
ok = exitstat == 0 .and. k == 3 .and. all(b(1, :) <= 0.3_R8P + 1e-14_R8P) .and. all(b(2, :) >= 0.3_R8P) .and. &
     all(b(2, :) - b(1, :) <= 1e-4_R8P)
print '(A,I0,A,6ES12.4,A,L1)', 'fossilizer --compare: exit ', exitstat, '  brackets ', b, '  ok ', ok
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_deviation_cli
