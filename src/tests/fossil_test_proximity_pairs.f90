!< FOSSIL test: pairs of bodies closer than a cutoff (issue #67, step 3).
!<
!< A 3 x 3 x 1 grid of unit cubes at spacing 1.5 (gaps 0.5 between neighbours along x and y, 0.5 sqrt(2) on the
!< diagonals): proximity_pairs_of with cutoff 0.6 gives exactly the 12 axis neighbours at 0.5; with cutoff 0.8 also the
!< 8 diagonal pairs at 0.7071; with 0.4 none. The same cubes merged into one surface (9 bodies): proximity_pairs gives
!< the same pairs. Each distance equals the brute-force distance of its pair.
program fossil_test_proximity_pairs

use fossil, only : surface_stl_object, proximity_pairs_of
use penf,   only : I4P, R8P

implicit none

type(surface_stl_object)  :: cubes(9), merged
integer(I4P), allocatable :: pairs(:, :), pairs2(:, :)
real(R8P),    allocatable :: d(:), d2(:)
integer(I4P)              :: i, j, k
logical                   :: ok1, ok2, ok3

k = 0
do j = 0, 2
   do i = 0, 2
      k = k + 1
      call cubes(k)%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
      call cubes(k)%translate(x=1.5_R8P * i, y=1.5_R8P * j, z=0._R8P) ; call cubes(k)%analyze
   enddo
enddo
call proximity_pairs_of(cubes, 0.6_R8P, pairs, d)
ok1 = size(d) == 12 .and. all(abs(d - 0.5_R8P) < 1e-15_R8P)
call proximity_pairs_of(cubes, 0.8_R8P, pairs, d)
ok1 = ok1 .and. size(d) == 20 .and. count(abs(d - sqrt(0.5_R8P)) < 1e-15_R8P) == 8
call proximity_pairs_of(cubes, 0.4_R8P, pairs, d)
ok1 = ok1 .and. size(d) == 0
print '(A,L1)', '1. list of bodies: ', ok1

merged = cubes(1)
do k = 2, 9
   call merged%merge_solids(cubes(k))
enddo
call proximity_pairs_of(cubes, 0.8_R8P, pairs, d)
call merged%proximity_pairs(0.8_R8P, pairs2, d2)
ok2 = merged%get_bodies_number() == 9 .and. size(d2) == 20 .and. all(abs(d2 - d) < 1e-15_R8P)
print '(A,L1)', '2. bodies of one surface: ', ok2

ok3 = .true.
do k = 1, size(d)
   call cubes(pairs(1, k))%distance_to(cubes(pairs(2, k)), d2(1))
   ok3 = ok3 .and. d2(1) == d(k)
enddo
print '(A,L1)', '3. distances: ', ok3
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3
endprogram fossil_test_proximity_pairs
