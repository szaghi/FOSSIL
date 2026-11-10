!< FOSSIL test: exact orient3d (issue #40).
program fossil_test_orient3d
!< FOSSIL test: exact orient3d (issue #40).
!<
!< 1. Near-coplanar grid, at two magnitudes T = 0 and T = 2^40: b, c, d lie exactly on the plane
!<    z = x - y; a = (T + 0.5 + i u, T + 0.5 + j u, (i - j) u + k u), u = spacing(T + 0.5),
!<    i, j = 0..31, k = -1, 0, 1. Every coordinate is exact and a is on the plane for k = 0, so
!<    the exact sign is sign(k) times a fixed sign (taken from a well-separated point). orient3d
!<    must return it for all 2 x 3 x 1024 points (exactly 0 for k = 0), keep it under the even
!<    permutation (b, c, a, d) and flip it under the odd permutations (b, a, c, d) and
!<    (b, c, d, a). The naive determinant gets many of them wrong (printed).
!< 2. Exactly coplanar points with dyadic coordinates on z = 2 x - 3 y + 0.5: exactly 0.
!< 3. Well-separated random quadruples: same sign as the naive determinant; the filter's double is
!<    returned unchanged.
!< 4. orient3d (vector_R8P) and orient3d_xyz (arrays) agree bit for bit.
!< 5. The device cases of `fossil_orient3d_cases` (issue #100 step 3c; the distance bench runs them on the device) through
!<    their batched loops as host code: the loop equals the scalar calls bit for bit, the known signs and exact zeros
!<    hold, the sign flips under the swap of a and b, and `line_crossing_before` gives the exact answer on every
!<    perturbation case (on the plane, an ulp and a unit off, entering and leaving, degenerate triangles).

use fossil_line_query,     only : line_crossing_before
use fossil_orient3d_cases, only : orient3d_quadruples, lcb_cases, orient3d_many, lcb_many
use fossil_predicates,     only : orient3d, orient3d_xyz
use penf,              only : I4P, I8P, R8P
use vecfor,            only : vector_R8P

implicit none

real(R8P)    :: t, u, a(3), b(3), c(3), d(3), o, ref, s0, r(12)
integer(I4P) :: m, i, j, k, n, wrong, naive_wrong, perm_wrong, zero_wrong, kk
logical      :: ok(5)
real(R8P),    allocatable :: q(:,:), oq(:), ltri(:,:), lab(:,:), lc(:) !< Device cases (5).
integer(I4P), allocatable :: expect(:), lax(:)                         !< Their expected signs, line axes.
logical,      allocatable :: len_(:), lexp(:), lbef(:)                 !< Entering flags, answers.
integer(I4P)              :: loop_bad, sign_bad, swap_bad, lcb_bad      !< Failures (5).

! 1. near-coplanar grid
wrong = 0 ; naive_wrong = 0 ; perm_wrong = 0 ; zero_wrong = 0
do m = 1, 2
   t = 0._R8P ; if (m == 2) t = 2._R8P**40
   u = spacing(t + 0.5_R8P)
   b = [t + 12._R8P, t + 3._R8P, 9._R8P]
   c = [t + 24._R8P, t + 7._R8P, 17._R8P]
   d = [t + 5._R8P, t + 30._R8P, -25._R8P]
   a = [t + 0.5_R8P, t + 0.5_R8P, 1._R8P] ! well separated: the fixed sign of a displacement in +z
   s0 = sign(1._R8P, orient3d_xyz(a, b, c, d))
   do k = -1, 1
      do j = 0, 31
         do i = 0, 31
            a = [t + 0.5_R8P + i * u, t + 0.5_R8P + j * u, (i - j) * u + k * u]
            ref = s0 * k
            o = orient3d_xyz(a, b, c, d)
            if (.not. same_sign(o, ref)) wrong = wrong + 1
            if (k == 0 .and. o /= 0._R8P) zero_wrong = zero_wrong + 1
            if (.not. same_sign(naive(a, b, c, d), ref)) naive_wrong = naive_wrong + 1
            if (.not. same_sign(orient3d_xyz(b, c, a, d), ref)) perm_wrong = perm_wrong + 1
            if (.not. same_sign(orient3d_xyz(b, a, c, d), -ref)) perm_wrong = perm_wrong + 1
            if (.not. same_sign(orient3d_xyz(b, c, d, a), -ref)) perm_wrong = perm_wrong + 1
         enddo
      enddo
   enddo
enddo
ok(1) = wrong == 0 .and. perm_wrong == 0 .and. zero_wrong == 0
print '(A,I0,A,I0,A,I0,A,I0,A)', 'near-coplanar grid: wrong ', wrong, '  coplanar not 0 ', zero_wrong, &
      '  permutation wrong ', perm_wrong, '  (naive determinant wrong ', naive_wrong, ' of 6144)'

! 2. exactly coplanar, dyadic coordinates on z = 2 x - 3 y + 0.5
call random_seed(size=n)
call random_seed(put=[(4242 + i, i = 1, n)])
wrong = 0
do kk = 1, 500
   call random_number(r)
   r = real(nint(r * 2._R8P**20), R8P) / 2._R8P**10 - 512._R8P ! dyadic, 10 fractional bits
   a = plane_point(r(1), r(2)) ; b = plane_point(r(3), r(4)) ; c = plane_point(r(5), r(6)) ; d = plane_point(r(7), r(8))
   if (orient3d_xyz(a, b, c, d) /= 0._R8P) wrong = wrong + 1
enddo
ok(2) = wrong == 0
print '(A,I0,A)', 'exactly coplanar -> 0: ', 500 - wrong, ' of 500'

! 3. well separated random quadruples
wrong = 0
do kk = 1, 2000
   call random_number(r)
   r = 2._R8P * r - 1._R8P
   a = r(1:3) ; b = r(4:6) ; c = r(7:9) ; d = r(10:12)
   o = orient3d_xyz(a, b, c, d)
   if (abs(naive(a, b, c, d)) > 1.e-6_R8P .and. o /= naive(a, b, c, d)) wrong = wrong + 1
enddo
ok(3) = wrong == 0
print '(A,I0)', 'random quadruples, filter double changed: ', wrong

! 4. vector and array interfaces agree
a = [0.5_R8P + 3 * spacing(0.5_R8P), 0.5_R8P + 7 * spacing(0.5_R8P), 2 * spacing(0.5_R8P)]
b = [12._R8P, 3._R8P, 9._R8P] ; c = [24._R8P, 7._R8P, 17._R8P] ; d = [5._R8P, 30._R8P, -25._R8P]
ok(4) = orient3d(vector_R8P(a(1), a(2), a(3)), vector_R8P(b(1), b(2), b(3)), vector_R8P(c(1), c(2), c(3)), &
                 vector_R8P(d(1), d(2), d(3))) == orient3d_xyz(a, b, c, d) .and. orient3d_xyz(a, b, c, d) /= 0._R8P
print '(A,L1)', 'vector and array interfaces agree: ', ok(4)

! 5. the device cases as host code
call orient3d_quadruples(q, expect)
allocate(oq(size(expect)))
call orient3d_many(size(expect, kind=I4P), q, oq)
loop_bad = 0 ; sign_bad = 0 ; swap_bad = 0
do kk = 1, size(expect)
   o = orient3d_xyz(q(1:3, kk), q(4:6, kk), q(7:9, kk), q(10:12, kk))
   if (transfer(o, 0_I8P) /= transfer(oq(kk), 0_I8P)) loop_bad = loop_bad + 1
   if (expect(kk) /= 2 .and. .not. same_sign(o, real(expect(kk), R8P))) sign_bad = sign_bad + 1
   if (.not. same_sign(orient3d_xyz(q(4:6, kk), q(1:3, kk), q(7:9, kk), q(10:12, kk)), -o)) swap_bad = swap_bad + 1
enddo
call lcb_cases(ltri, lax, lab, lc, len_, lexp)
allocate(lbef(size(lexp)))
call lcb_many(size(lexp, kind=I4P), ltri, lax, lab, lc, len_, lbef)
lcb_bad = 0
do kk = 1, size(lexp)
   if ((lbef(kk) .neqv. lexp(kk)) .or. (lbef(kk) .neqv. line_crossing_before(reshape(ltri(:, kk), [3, 3]), lax(kk), &
                                                                            lab(1, kk), lab(2, kk), lc(kk), len_(kk)))) &
      lcb_bad = lcb_bad + 1
enddo
ok(5) = loop_bad == 0 .and. sign_bad == 0 .and. swap_bad == 0 .and. lcb_bad == 0
print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', 'device cases on the host: ', size(expect), ' quadruples, loop vs scalar ', &
      loop_bad, ', known signs wrong ', sign_bad, ', swap not flipped ', swap_bad, '; ', size(lexp), &
      ' line_crossing_before cases wrong ', lcb_bad

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   pure function plane_point(x, y) result(p)
   !< The point (x, y) on z = 2 x - 3 y + 0.5 (exact for the dyadic inputs).
   real(R8P), intent(in) :: x, y !< Coordinates.
   real(R8P)             :: p(3) !< Point.

   p = [x, y, 2._R8P * x - 3._R8P * y + 0.5_R8P]
   endfunction plane_point

   pure function naive(a, b, c, d) result(det)
   !< Naive determinant det[a - d; b - d; c - d] (the filter's formula).
   real(R8P), intent(in) :: a(3), b(3), c(3), d(3) !< Points.
   real(R8P)             :: det                    !< Determinant.
   real(R8P)             :: p(3), q(3), w(3)       !< Rows.

   p = a - d ; q = b - d ; w = c - d
   det = p(3) * (q(1) * w(2) - w(1) * q(2)) + q(3) * (w(1) * p(2) - p(1) * w(2)) + w(3) * (p(1) * q(2) - q(1) * p(2))
   endfunction naive

   pure function same_sign(a, b) result(same)
   !< True if a and b have the same sign (both zero counts as the same).
   real(R8P), intent(in) :: a, b !< Values.
   logical               :: same !< Result.

   same = (a > 0._R8P .and. b > 0._R8P) .or. (a < 0._R8P .and. b < 0._R8P) .or. (a == 0._R8P .and. b == 0._R8P)
   endfunction same_sign
endprogram fossil_test_orient3d
