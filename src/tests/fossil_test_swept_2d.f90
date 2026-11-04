!< FOSSIL test: 2-D swept geometry between two poses (issue #41 step 5).
program fossil_test_swept_2d
!< FOSSIL test: 2-D swept geometry between two poses (issue #41 step 5).
!<
!< Translation (step 5a, exact path):
!< - a square translating across a uniform block: every time-averaged face fluid length equals an
!<   independent analytic integration (breakpoints from the square's sides, two-point Gauss on every
!<   piece) to 1e-14;
!< - the square, a non-convex heptagon, and the heptagon translating at a common rotated pose: the
!<   uniform-flow residual (swept area + time integral of the outflow of the body velocity) is
!<   round-off in every cell; the projected walls sum to 0 over the block (closed body inside it);
!<   cleared and covered cells occur and match the fractions;
!< - negative control: the same motions with the Gauss path (`exact=.false.`, 4 nodes, no
!<   breakpoints) leave residuals far above round-off;
!< - refusals: decreasing edges, nq = 0, a poses size other than nq.

use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, STATUS_INVALID_INPUT, STATUS_OK
use penf,   only : I4P, R8P

implicit none

type(curve_object)     :: c
type(swept_cells_2d_t) :: sw
real(R8P)              :: x(9), y(7), err, res, d(2), sq(4)
integer(I4P)           :: i, j, st, st2, st3
logical                :: ok(5)

x = [(0.25_R8P * i, i = 0, 8)] ; y = [(0.25_R8P * i, i = 0, 6)]
! square [0.3, 1.1] x [0.2, 0.9] translating by d
sq = [0.3_R8P, 1.1_R8P, 0.2_R8P, 0.9_R8P] ; d = [0.37_R8P, 0.21_R8P]
call c%adopt_loops(u=[sq(1), sq(2), sq(2), sq(1)], v=[sq(3), sq(3), sq(4), sq(4)], loop_sizes=[4_I4P])
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P, d), x, y, sw, status=st)
err = 0._R8P
do j = 1, 6
   do i = 1, 9
      err = max(err, abs(sw%face_u(i, j) - face_ref(1, x(i), y(j), y(j + 1))))
   enddo
enddo
do j = 1, 7
   do i = 1, 8
      err = max(err, abs(sw%face_v(i, j) - face_ref(2, y(j), x(i), x(i + 1))))
   enddo
enddo
ok(1) = st == STATUS_OK .and. sw%exact .and. err < 1.e-14_R8P
print '(A,ES9.1,A,I0,A,L1)', 'square, translation: max face error vs analytic ', err, '  evaluations/face ', sw%evaluations, &
      '  ok ', ok(1)
ok(2) = balanced(sw, 'square', .true.) .and. count(sw%cleared) > 0 .and. count(sw%covered) > 0
! non-convex heptagon, and the same at a common rotated pose
call c%adopt_loops(u=[0.4_R8P, 1.3_R8P, 1.0_R8P, 1.45_R8P, 0.9_R8P, 0.62_R8P, 0.35_R8P], &
                   v=[0.3_R8P, 0.45_R8P, 0.7_R8P, 1.1_R8P, 1.2_R8P, 0.8_R8P, 0.95_R8P], loop_sizes=[7_I4P])
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P, [0.29_R8P, -0.13_R8P]), x, y, sw, status=st)
ok(3) = st == STATUS_OK .and. balanced(sw, 'heptagon', .true.)
call swept_cells_2d(c, pose_2d(0.3_R8P, [0.2_R8P, -0.25_R8P]), pose_2d(0.3_R8P, [0.53_R8P, -0.08_R8P]), x, y, sw, status=st)
ok(3) = ok(3) .and. st == STATUS_OK .and. sw%exact .and. balanced(sw, 'heptagon, rotated pose', .true.)
! negative control: Gauss without breakpoints
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P, [0.29_R8P, -0.13_R8P]), x, y, sw, exact=.false., status=st)
res = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
ok(4) = st == STATUS_OK .and. .not. sw%exact .and. res > 1.e-6_R8P
print '(A,ES9.1,A,L1)', 'negative control (Gauss, 4 nodes, no breakpoints): max residual ', res, '  ok ', ok(4)
! refusals
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P), x(9:1:-1), y, sw, status=st)
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P), x, y, sw, nq=0_I4P, status=st2)
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P), x, y, sw, nq=3_I4P, poses=[pose_2d(0._R8P), pose_2d(0._R8P)], &
                    status=st3)
ok(5) = st == STATUS_INVALID_INPUT .and. st2 == STATUS_INVALID_INPUT .and. st3 == STATUS_INVALID_INPUT
print '(A,L1)', 'refusals: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function balanced(s, label, exact_expected) result(good)
   !< Residual round-off in every closed cell, projected walls summing to 0, cleared/covered consistent.
   type(swept_cells_2d_t), intent(in) :: s              !< Swept geometry.
   character(*),           intent(in) :: label          !< Label.
   logical,                intent(in) :: exact_expected !< Exact path expected.
   logical                            :: good           !< Result.
   real(R8P)                          :: r, w(2)        !< Work.

   r = maxval(abs(s%residual), mask=s%alpha0 >= 0._R8P .and. s%alpha1 >= 0._R8P)
   w = [sum(s%wall(1, :, :)), sum(s%wall(2, :, :))]
   good = (s%exact .eqv. exact_expected) .and. r < 1.e-14_R8P .and. maxval(abs(w)) < 1.e-13_R8P .and. &
          all(.not. s%cleared .or. (s%alpha0 == 0._R8P .and. s%alpha1 > 0._R8P)) .and. &
          all(.not. s%covered .or. (s%alpha0 > 0._R8P .and. s%alpha1 == 0._R8P)) .and. all(s%alpha0 >= 0._R8P)
   print '(A,A,ES9.1,A,2ES9.1,A,2I3,A,L1)', label, ': max uniform-flow residual ', r, '  sum of walls ', w, &
         '  cleared/covered ', count(s%cleared), count(s%covered), '  ok ', good
   endfunction balanced

   function face_ref(ax, a, lo, hi) result(avg)
   !< Analytic time average of the fluid length of the face (ax = 1: u = a, v in [lo, hi]; ax = 2:
   !< v = a, u in [lo, hi]) for the square sq translating by s d: breakpoints where a side crosses the
   !< face's line or its ends, two-point Gauss on every piece (piecewise linear, jumps at the sides).
   integer,   intent(in) :: ax          !< Face normal axis.
   real(R8P), intent(in) :: a, lo, hi   !< Line, span.
   real(R8P)             :: avg         !< Average.
   real(R8P)             :: g(8), cand(6), sa, sb !< Geometry, candidate breakpoints, piece.
   real(R8P)             :: b(6)        !< Breakpoints.
   integer               :: n, q, p     !< Work.

   if (ax == 1) then
      g = [sq(1), sq(2), sq(3), sq(4), d(1), d(2), lo, hi] ! n0, n1, t0, t1, dn, dt, lo, hi
   else
      g = [sq(3), sq(4), sq(1), sq(2), d(2), d(1), lo, hi]
   endif
   cand = [(a - g(1)) / g(5), (a - g(2)) / g(5), (lo - g(3)) / g(6), (lo - g(4)) / g(6), (hi - g(3)) / g(6), &
           (hi - g(4)) / g(6)]
   n = 0
   do q = 1, 6
      if (cand(q) > 0._R8P .and. cand(q) < 1._R8P) then
         n = n + 1 ; b(n) = cand(q)
      endif
   enddo
   do q = 2, n ! sort
      do p = q, 2, -1
         if (b(p) < b(p - 1)) b(p - 1:p) = [b(p), b(p - 1)]
      enddo
   enddo
   avg = 0._R8P
   do q = 0, n
      sa = 0._R8P ; if (q > 0) sa = b(q)
      sb = 1._R8P ; if (q < n) sb = b(q + 1)
      ! two-point Gauss (interior nodes: the length jumps where the face line leaves the square)
      avg = avg + 0.5_R8P * (sb - sa) * (flen(g, a, 0.5_R8P * (sa + sb) - (sb - sa) / (2._R8P * sqrt(3._R8P))) + &
                                         flen(g, a, 0.5_R8P * (sa + sb) + (sb - sa) / (2._R8P * sqrt(3._R8P))))
   enddo
   endfunction face_ref

   pure function flen(g, a, sv) result(l)
   !< Fluid length of the face at instant sv (g: square extent normal and along, motion, span).
   real(R8P), intent(in) :: g(8), a, sv !< Geometry, line, instant.
   real(R8P)             :: l           !< Fluid length.

   l = g(8) - g(7)
   if (a > g(1) + sv * g(5) .and. a < g(2) + sv * g(5)) &
      l = l - max(0._R8P, min(g(8), g(4) + sv * g(6)) - max(g(7), g(3) + sv * g(6)))
   endfunction flen
endprogram fossil_test_swept_2d
