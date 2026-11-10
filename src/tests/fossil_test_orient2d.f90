!< FOSSIL test: exact orient2d (issue #36 step 1).
program fossil_test_orient2d
!< FOSSIL test: exact orient2d (issue #36 step 1).
!<
!< 1. Near-collinear grid (Shewchuk's classic failure case), at two magnitudes T = 0 and
!<    T = 2^40: a = (T + 0.5 + i u, T + 0.5 + j u), u = spacing(T + 0.5), i, j = 0..31,
!<    b = (T + 12, T + 12), c = (T + 24, T + 24). Every coordinate is exact, and
!<    (a - c) x (b - c) = 12 (a_y - a_x), so the exact sign is sign(j - i). orient2d must
!<    return it for all 2 x 1024 points, and stay consistent under the cyclic permutation
!<    (b, c, a) and flip under (a, c, b). The naive formula gets many of them wrong (printed).
!< 2. Exactly collinear points (t, 3t) with 3t exact: orient2d is exactly 0.
!< 3. Well-separated random triples: same sign as the naive formula.
!< 4. orient2d (vector_R8P) and orient2d_xy (reals) agree bit for bit.
!< 5. Coincident points (issue #100 step 3a early exit): +0 exactly (the bits of +0, as the exact stage gives), in every
!<    position and at magnitudes 2^-40 .. 2^40; one ulp apart, the exact nonzero sign (u > 0 by construction).

use fossil_predicates, only : orient2d, orient2d_xy
use penf,              only : I4P, I8P, R8P
use vecfor,            only : vector_R8P

implicit none

real(R8P)    :: t, u, ax, ay, o, onaive, ref, r(6)
integer(I4P) :: m, i, j, k, wrong, naive_wrong, perm_wrong
logical      :: ok(5)
type(vector_R8P) :: va, vb, vc

! 1. near-collinear grid
wrong = 0 ; naive_wrong = 0 ; perm_wrong = 0
do m = 1, 2
   t = 0._R8P ; if (m == 2) t = 2._R8P**40
   u = spacing(t + 0.5_R8P)
   do j = 0, 31
      do i = 0, 31
         ax = t + 0.5_R8P + i * u
         ay = t + 0.5_R8P + j * u
         ref = real(j - i, R8P)
         o = orient2d_xy(ax, ay, t + 12._R8P, t + 12._R8P, t + 24._R8P, t + 24._R8P)
         if (.not. same_sign(o, ref)) wrong = wrong + 1
         onaive = (ax - (t + 24._R8P)) * (t + 12._R8P - (t + 24._R8P)) - (ay - (t + 24._R8P)) * (t + 12._R8P - (t + 24._R8P))
         if (.not. same_sign(onaive, ref)) naive_wrong = naive_wrong + 1
         if (.not. same_sign(orient2d_xy(t + 12._R8P, t + 12._R8P, t + 24._R8P, t + 24._R8P, ax, ay), ref)) &
            perm_wrong = perm_wrong + 1
         if (.not. same_sign(orient2d_xy(ax, ay, t + 24._R8P, t + 24._R8P, t + 12._R8P, t + 12._R8P), -ref)) &
            perm_wrong = perm_wrong + 1
      enddo
   enddo
enddo
ok(1) = wrong == 0 .and. perm_wrong == 0
print '(A,I0,A,I0,A,I0,A)', 'near-collinear grid: wrong ', wrong, '  permutation wrong ', perm_wrong, &
      '  (naive formula wrong ', naive_wrong, ' of 2048)'

! 2. exactly collinear
ok(2) = orient2d_xy(0._R8P, 0._R8P, 1.25_R8P, 3.75_R8P, 1024.75_R8P, 3074.25_R8P) == 0._R8P .and. &
        orient2d_xy(0.5_R8P, 1.5_R8P, 2._R8P**40 + 0.5_R8P, 3._R8P * (2._R8P**40 + 0.5_R8P), 1.25_R8P, 3.75_R8P) == 0._R8P
print '(A,L1)', 'exactly collinear -> 0: ', ok(2)

! 3. well separated random triples
call random_seed(size=k)
call random_seed(put=[(12345 + i, i = 1, k)])
wrong = 0
do i = 1, 2000
   call random_number(r)
   r = 2._R8P * r - 1._R8P
   o = orient2d_xy(r(1), r(2), r(3), r(4), r(5), r(6))
   onaive = (r(1) - r(5)) * (r(4) - r(6)) - (r(2) - r(6)) * (r(3) - r(5))
   if (abs(onaive) > 1.e-6_R8P .and. .not. same_sign(o, onaive)) wrong = wrong + 1
enddo
ok(3) = wrong == 0
print '(A,I0)', 'random triples, sign disagreements: ', wrong

! 4. vector and real interfaces agree
va%x = 0.5_R8P + 3 * spacing(0.5_R8P) ; va%y = 0.5_R8P + 7 * spacing(0.5_R8P) ; va%z = 99._R8P
vb%x = 12._R8P ; vb%y = 12._R8P ; vb%z = -1._R8P
vc%x = 24._R8P ; vc%y = 24._R8P ; vc%z = 3._R8P
ok(4) = orient2d(va, vb, vc) == orient2d_xy(va%x, va%y, vb%x, vb%y, vc%x, vc%y) .and. orient2d(va, vb, vc) > 0._R8P
print '(A,L1)', 'vector and real interfaces agree: ', ok(4)

! 5. coincident and one-ulp-apart points
wrong = 0
do m = -40, 40, 8
   t = 0.75_R8P * 2._R8P**m ; u = spacing(t)
   if (transfer(orient2d_xy(t, t, t, t, t, 2._R8P * t), 0_I8P) /= 0_I8P) wrong = wrong + 1 ! a == b
   if (transfer(orient2d_xy(t, t, 3._R8P * t, -t, t, t), 0_I8P) /= 0_I8P) wrong = wrong + 1 ! a == c
   if (transfer(orient2d_xy(-t, 5._R8P * t, t, t, t, t), 0_I8P) /= 0_I8P) wrong = wrong + 1 ! b == c
   if (transfer(orient2d_xy(t, t, t, t, t, t), 0_I8P) /= 0_I8P) wrong = wrong + 1          ! all three
   if (.not. same_sign(orient2d_xy(t, t, t + u, t, t, t + 1._R8P), 1._R8P)) wrong = wrong + 1 ! (a-c) x (b-c) = u
   if (.not. same_sign(orient2d_xy(t + u, t, t, t, t, t + 1._R8P), -1._R8P)) wrong = wrong + 1
enddo
ok(5) = wrong == 0
print '(A,I0)', 'coincident / one-ulp-apart points, wrong: ', wrong

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   pure function same_sign(a, b) result(same)
   !< True if a and b have the same sign (both zero counts as the same).
   real(R8P), intent(in) :: a, b !< Values.
   logical               :: same !< Result.

   same = (a > 0._R8P .and. b > 0._R8P) .or. (a < 0._R8P .and. b < 0._R8P) .or. (a == 0._R8P .and. b == 0._R8P)
   endfunction same_sign
endprogram fossil_test_orient2d
