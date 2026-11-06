!< FOSSIL test: fossilizer --distance (issue #67, step 5).
!<
!< `fossilizer -i cube-02.stl --distance <cube-02 translated by its width + 0.25 along x>`: the printed distance is
!< 0.25 and the bodies do not touch; against the cube itself the distance is 0 and they touch. Read back from the lines.
!< Runs exe/fossilizer from the repository root, like every FOSSIL test.
program fossil_test_distance_cli

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, vector_R8P

implicit none

character(*), parameter  :: MOVED = 'fossil_test_distance_cli_moved.stl'
character(*), parameter  :: OUT   = 'fossil_test_distance_cli.txt'
type(surface_stl_object) :: cube
type(vector_R8P)         :: lo, hi
real(R8P)                :: d(2)
logical                  :: hit(2)
integer(I4P)             :: exitstat(2), u
logical                  :: ok

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
lo = cube%get_bmin() ; hi = cube%get_bmax()
call cube%translate(delta=(hi%x - lo%x + 0.25_R8P) * ex_R8P)
call cube%analyze
call cube%save_into_file(file_name=MOVED)
call run(MOVED, 1)
call run('src/tests/cube-02.stl', 2)
open(newunit=u, file=MOVED, status='old') ; close(u, status='delete')
ok = all(exitstat == 0) .and. abs(d(1) - 0.25_R8P) < 1e-14_R8P .and. .not. hit(1) .and. d(2) == 0._R8P .and. hit(2)
print '(A,2ES12.4,A,2L2,A,L1)', 'fossilizer --distance: distances ', d, '  touching ', hit, '  ok ', ok
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
contains
   subroutine run(other, k)
   !< Run fossilizer --distance other, read the distance line into d(k), hit(k).
   character(*), intent(in) :: other !< The other body's file.
   integer(I4P), intent(in) :: k     !< Slot.
   character(256)           :: line
   integer(I4P)             :: u, ios, es

   call execute_command_line('exe/fossilizer -i src/tests/cube-02.stl --distance '//other//' > '//OUT, exitstat=es)
   exitstat(k) = es
   open(newunit=u, file=OUT, status='old')
   do
      read(u, '(A)', iostat=ios) line
      if (ios /= 0) exit
      if (index(line, '  distance ') == 1) then
         read(line(11:index(line, 'touching') - 1), *) d(k)
         hit(k) = index(line, 'crossing T') > 0
      endif
   enddo
   close(u, status='delete')
   endsubroutine run
endprogram fossil_test_distance_cli
