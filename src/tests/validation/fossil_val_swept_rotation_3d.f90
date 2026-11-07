!< FOSSIL validation: exact 3-D swept geometry under rotation (issue #69, step 10), full-size copy of
!< fossil_test_swept_rotation_3d: the three bodies on 6^3 blocks.
!<
!< The screw motion of `pose_interpolate` with a rotation takes the exact path by default: per grid plane, certified
!< breakpoints (a mesh vertex on the plane, a mesh edge meeting a grid line, a grid node on a triangle) and Gauss-
!< Legendre doubled to round-off on every piece. Cases on a 6^3 block:
!< - the unit cube turning 0.3 rad about (1, 1, 1) through its centre while translating (its faces start on grid
!<   planes: zero zones at s = 0);
!< - the L-prism (non-convex) turning 0.5 rad about a tilted off-centre axis while translating;
!< - the 32-face cone pitching 0.4 rad about a horizontal axis.
!< Each: `sw%exact`, the rigid-field residual at round-off (<= 1e-14) in every closed cell, the projected walls
!< summing to zero over the block, and the Gauss path (`exact=.false.`) converging to the exact faces (nq = 8 at
!< least 4 times closer than nq = 2).
!< Negative control (base 04e5bbf: rotation by Gauss only, 4 nodes): residuals 1.0e-3, 9.3e-3, 2.8e-3.
program fossil_val_swept_rotation_3d

use fossil,                   only : surface_stl_object, rigid_pose_t, pose_identity, pose_from_axis_angle, &
                                     swept_cells_3d_t, swept_cells_3d, STATUS_OK
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(rigid_pose_t)        :: p0, p1
real(R8P)                 :: e(7), c(3)
integer(I4P)              :: i
logical                   :: ok(3)

call bench_bodies(b)
p0 = pose_identity()
call to_surface(b(1)%soup, s) ! unit cube [0, 1]^3
e = [(-0.25_R8P + 0.4_R8P * i / 1.5_R8P, i = 0, 6)]
c = 0.5_R8P
p1 = pose_from_axis_angle([1._R8P, 1._R8P, 1._R8P], 0.3_R8P)
p1%t = c - matmul(p1%r, c) + [0.13_R8P, -0.07_R8P, 0.05_R8P]
ok(1) = check('cube 0.3 rad about (1,1,1) + shift')
call to_surface(b(8)%soup, s) ! L-prism [0, 2]^2 x [0, 1]; the block holds the whole motion
e = [(-1.1_R8P + 0.6_R8P * i, i = 0, 6)]
c = [0.7_R8P, 0.9_R8P, 0.4_R8P]
p1 = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0.5_R8P)
p1%t = c - matmul(p1%r, c) + [0.05_R8P, 0.1_R8P, -0.08_R8P]
ok(2) = check('L-prism 0.5 rad, tilted axis + shift')
call to_surface(b(7)%soup, s) ! cone, 32 faces
c = 0.5_R8P
p1 = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0.4_R8P)
p1%t = c - matmul(p1%r, c)
ok(3) = check('cone pitching 0.4 rad')
print '(A,L1)', 'Are all tests passed? ', all(ok)

contains
   function check(label) result(good)
   !< Exact path between p0 and p1: residual, walls, the Gauss path converging to it.
   character(*), intent(in) :: label
   logical                  :: good
   type(swept_cells_3d_t)   :: se, g4, g16
   real(R8P)                :: r, w(3), e4, e16
   integer(I4P)             :: st, s4, s16

   call swept_cells_3d(s, p0, p1, e, e, e, se, status=st)
   call swept_cells_3d(s, p0, p1, e, e, e, g4, nq=2_I4P, exact=.false., status=s4)
   call swept_cells_3d(s, p0, p1, e, e, e, g16, nq=8_I4P, exact=.false., status=s16)
   r = maxval(abs(se%residual), mask=se%alpha0 >= 0._R8P .and. se%alpha1 >= 0._R8P)
   w = [sum(se%wall(1, :, :, :)), sum(se%wall(2, :, :, :)), sum(se%wall(3, :, :, :))]
   e4 = max(maxval(abs(g4%face_x - se%face_x)), maxval(abs(g4%face_y - se%face_y)), maxval(abs(g4%face_z - se%face_z)))
   e16 = max(maxval(abs(g16%face_x - se%face_x)), maxval(abs(g16%face_y - se%face_y)), &
             maxval(abs(g16%face_z - se%face_z)))
   good = st == STATUS_OK .and. s4 == STATUS_OK .and. s16 == STATUS_OK .and. se%exact .and. r <= 1e-14_R8P .and. &
          all(abs(w) <= 1e-14_R8P) .and. e16 <= 0.25_R8P * e4
   print '(A,A,ES9.1,A,3ES9.1,A,2ES9.1,A,I0,A,L1)', label, ': residual ', r, '  walls ', w, &
         '  |Gauss 2, 8 - exact| ', e4, e16, '  evaluations/plane ', se%evaluations, '  ok ', good
   endfunction check
endprogram fossil_val_swept_rotation_3d
