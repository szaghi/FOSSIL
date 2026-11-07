!< FOSSIL test: exact 2-D swept geometry under rotation (issue #69, step 9).
!<
!< The screw motion of `pose_interpolate` with a rotation takes the exact path by default: certified breakpoints per
!< face (a wall vertex crossing the face's line, an end node crossing a wall segment), Gauss-Legendre doubled to
!< round-off on every piece. Cases: a square rotating by 0.4 rad about an off-centre point while translating; the
!< same square turning 1.2 rad; a NACA 0012 (200 vertices) pitching 0.15 rad about its quarter chord. Each:
!< - `sw%exact`, and the rigid-field residual at round-off (<= 1e-14) in every closed cell;
!< - the projected walls sum to zero over the block;
!< - the Gauss path (`exact=.false.`) converges to the exact faces: nq = 16 at least 4 times closer than nq = 4.
!< Negative control (base 04e5bbf: rotation by Gauss only, 4 nodes): residuals 3.8e-3, 1.4e-2, 1.7e-4.
program fossil_test_swept_rotation_2d

use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, STATUS_OK
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter   :: PI = acos(-1._R8P)
type(curve_object)     :: c
type(rigid_pose_2d_t)  :: p0, p1
real(R8P)              :: x(9), y(7), u(200), v(200), sq(4)
integer(I4P)           :: i
logical                :: ok(3)

x = [(0.25_R8P * i, i = 0, 8)] ; y = [(0.25_R8P * i, i = 0, 6)]
sq = [0.3_R8P, 1.1_R8P, 0.2_R8P, 0.9_R8P]
call c%adopt_loops(u=[sq(1), sq(2), sq(2), sq(1)], v=[sq(3), sq(3), sq(4), sq(4)], loop_sizes=[4_I4P])
p0 = pose_2d(0._R8P)
p1 = about(0.4_R8P, [1._R8P, 0.75_R8P], [0.1_R8P, 0.05_R8P])
ok(1) = check('square, 0.4 rad off-centre + shift')
p1 = about(1.2_R8P, [0.7_R8P, 0.55_R8P], [0.05_R8P, -0.03_R8P])
ok(2) = check('square, 1.2 rad + shift')
do i = 0, 99
   u(i + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - i) / 100)) ; v(i + 1) = thickness(u(i + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do i = 1, 99
   u(101 + i) = 0.5_R8P * (1._R8P - cos(PI * i / 100)) ; v(101 + i) = -thickness(u(101 + i))
enddo
u = u + 0.4_R8P ; v = v + 0.75_R8P
call c%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
p1 = about(0.15_R8P, [0.65_R8P, 0.75_R8P], [0._R8P, 0._R8P])
ok(3) = check('NACA 0012 pitching 0.15 rad')
print '(A,L1)', 'Are all tests passed? ', all(ok)

contains
   function check(label) result(good)
   !< Exact path between p0 and p1: residual, walls, and the Gauss path converging to it.
   character(*), intent(in) :: label
   logical                  :: good
   type(swept_cells_2d_t)   :: se, g4, g16
   real(R8P)                :: r, w(2), e4, e16
   integer(I4P)             :: st, s4, s16

   call swept_cells_2d(c, p0, p1, x, y, se, status=st)
   call swept_cells_2d(c, p0, p1, x, y, g4, nq=4_I4P, exact=.false., status=s4)
   call swept_cells_2d(c, p0, p1, x, y, g16, nq=16_I4P, exact=.false., status=s16)
   r = maxval(abs(se%residual), mask=se%alpha0 >= 0._R8P .and. se%alpha1 >= 0._R8P)
   w = [sum(se%wall(1, :, :)), sum(se%wall(2, :, :))]
   e4 = max(maxval(abs(g4%face_u - se%face_u)), maxval(abs(g4%face_v - se%face_v)))
   e16 = max(maxval(abs(g16%face_u - se%face_u)), maxval(abs(g16%face_v - se%face_v)))
   good = st == STATUS_OK .and. s4 == STATUS_OK .and. s16 == STATUS_OK .and. se%exact .and. r <= 1e-14_R8P .and. &
          all(abs(w) <= 1e-14_R8P) .and. e16 <= 0.25_R8P * e4
   print '(A,A,ES9.1,A,2ES9.1,A,2ES9.1,A,I0,A,L1)', label, ': residual ', r, '  walls ', w, &
         '  |Gauss 4, 16 - exact| ', e4, e16, '  evaluations/face ', se%evaluations, '  ok ', good
   endfunction check

   pure function about(angle, centre, shift) result(p)
   !< Rotation by angle about centre, then a translation.
   real(R8P), intent(in) :: angle, centre(2), shift(2)
   type(rigid_pose_2d_t) :: p

   p = pose_2d(angle)
   p%t = centre - matmul(p%r, centre) + shift
   endfunction about

   pure function thickness(xx) result(yy)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xx
   real(R8P)             :: yy

   yy = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.126_R8P * xx - 0.3516_R8P * xx**2 + 0.2843_R8P * xx**3 - 0.1036_R8P * xx**4)
   endfunction thickness
endprogram fossil_test_swept_rotation_2d
