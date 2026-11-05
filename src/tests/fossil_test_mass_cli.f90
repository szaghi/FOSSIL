!< FOSSIL test: fossilizer --mass_properties (issue #59, roadmap #58 WP1, step 4).
!<
!< 1. `fossilizer -i src/tests/cube-02.stl --mass_properties --density 2`: mass 2, centre (0.5, 0.5, 0.5), principal
!<    moments m/6 = 1/3, inertia tensor diagonal (1e-14), read back from the printed lines.
!< 2. `fossilizer -i src/tests/dragon.stl --mass_properties`: the dragon is not watertight, the report says so.
!< Runs exe/fossilizer from the repository root, like every FOSSIL test.
program fossil_test_mass_cli

use penf, only : I4P, R8P

implicit none

character(*), parameter :: OUT = 'fossil_test_mass_cli.txt'
character(256)          :: line
real(R8P)               :: mass, center(3), row(3), moments(3), inertia(3,3)
integer(I4P)            :: u, exitstat, ios, k
logical                 :: ok, rejected

call execute_command_line('exe/fossilizer -i src/tests/cube-02.stl --mass_properties --density 2 > '//OUT, &
                          exitstat=exitstat)
mass = -1._R8P ; center = -1._R8P ; moments = -1._R8P ; inertia = -1._R8P
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, '  mass  ') == 1) read(line(23:), *) mass
   if (index(line, '  centre of mass') == 1) read(line(23:), *) center
   if (index(line, '  principal moments') == 1) read(line(23:), *) moments
   if (index(line, '  inertia tensor') == 1) then
      do k = 1, 3
         read(u, '(A)') line ; read(line(23:), *) row ; inertia(k,:) = row
      enddo
   endif
enddo
close(u, status='delete')
ok = exitstat == 0 .and. abs(mass - 2._R8P) <= 1e-14_R8P .and. maxval(abs(center - 0.5_R8P)) <= 1e-14_R8P .and. &
     maxval(abs(moments - 1._R8P / 3._R8P)) <= 1e-14_R8P .and. abs(inertia(1,2)) <= 1e-14_R8P .and. &
     abs(inertia(1,3)) <= 1e-14_R8P .and. abs(inertia(2,3)) <= 1e-14_R8P
print '(A,I0,A,ES23.15,A,3F9.5,A,L1)', '1. cube-02: exit ', exitstat, '  mass ', mass, '  moments', moments, '  ok ', ok

call execute_command_line('exe/fossilizer -i src/tests/dragon.stl --mass_properties > '//OUT//' 2>/dev/null', &
                          exitstat=exitstat)
rejected = .false.
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, 'mass properties: not a closed body') == 1) rejected = .true.
enddo
close(u, status='delete')
print '(A,I0,A,L1)', '2. dragon: exit ', exitstat, '  reported not closed ', rejected
ok = ok .and. exitstat == 0 .and. rejected

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_mass_cli
