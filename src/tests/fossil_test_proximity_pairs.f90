!< FOSSIL test: pairs of bodies closer than a cutoff (issue #67, step 3).
!<
!< A 3 x 3 x 1 grid of unit cubes at spacing 1.5 (gaps 0.5 between neighbours along x and y, 0.5 sqrt(2) on the
!< diagonals): proximity_pairs_of with cutoff 0.6 gives exactly the 12 axis neighbours at 0.5; with cutoff 0.8 also the
!< 8 diagonal pairs at 0.7071; with 0.4 none. The same cubes merged into one surface (9 bodies): proximity_pairs gives
!< the same pairs. Each distance equals the brute-force distance of its pair.
!< 4. Issue #76 (the result is now counted, then allocated once): pairs and distances equal, bit for bit and in order, the
!<    loop over all pairs i < j with distance_to <= cutoff, on cubes at irregular offsets (cutoffs 0.3, 0.8, 5).
program fossil_test_proximity_pairs

use fossil, only : surface_stl_object, proximity_pairs_of
use penf,   only : I4P, R8P

implicit none

type(surface_stl_object)  :: cubes(9), merged
integer(I4P), allocatable :: pairs(:, :), pairs2(:, :)
real(R8P),    allocatable :: d(:), d2(:)
integer(I4P)              :: i, j, k
logical                   :: ok1, ok2, ok3, ok4
real(R8P)                 :: dd, cut(3)
integer(I4P)              :: m, c

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
ok4 = .true.
do k = 1, 9
   call cubes(k)%translate(x=0.13_R8P * mod(k, 4), y=-0.07_R8P * mod(k, 3), z=0.05_R8P * k) ; call cubes(k)%analyze
enddo
cut = [0.3_R8P, 0.8_R8P, 5._R8P]
do c = 1, 3
   call proximity_pairs_of(cubes, cut(c), pairs, d)
   m = 0
   do i = 1, 9
      do j = i + 1, 9
         call cubes(i)%distance_to(cubes(j), dd)
         if (dd > cut(c)) cycle
         m = m + 1
         if (m > size(d)) then ; ok4 = .false. ; cycle ; endif
         ok4 = ok4 .and. pairs(1, m) == i .and. pairs(2, m) == j .and. d(m) == dd
      enddo
   enddo
   ok4 = ok4 .and. m == size(d) .and. size(pairs, 1) == 2 .and. size(pairs, 2) == m
   print '(A,F4.1,A,I0)', '   cutoff ', cut(c), ': pairs ', m
enddo
print '(A,L1)', '4. equal to the loop over all pairs: ', ok4
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4
endprogram fossil_test_proximity_pairs
