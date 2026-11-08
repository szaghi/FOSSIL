!run figref-moving-bodies-screw fig_ref_moving_bodies_screw
!fig ref-moving-bodies-screw
!< Figure of the reference page "Moving bodies", section "Poses": the screw motion of `pose_interpolate` from the
!< identity to p1 (a turn of 120 degrees about z, translation (1.8, 0.6, 1.2)). The unit cube centred at the body-frame
!< origin is drawn by its edges at s = 0, 1/5, ..., 1 (dark to yellow); the body-frame origin moves on the straight line
!< s t1 (teal) while a corner turns about it with the angle linear in s (its path in ember). The cube of
!< pose_compose(half, half) (half = the pose at s = 1/2) is drawn in teal: it has the end rotation but not the end
!< position, because the screw motion is not a power of a pose when t1 is off the rotation axis.
program fig_ref_moving_bodies_screw
use fossil,     only : rigid_pose_t, pose_identity, pose_from_axis_angle, pose_interpolate, pose_compose, to_world
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter :: PI = acos(-1._R8P)
integer(I4P), parameter :: NS = 5, NT = 60
type(rigid_pose_t)      :: p0, p1, q, half, twice
type(vector_R8P)        :: a(12), b(12), path_a(NT), path_b(NT)
real(R8P)               :: corner(3), x(3)
integer(I4P)            :: k, u
character(120)          :: caption

p0 = pose_identity()
p1 = pose_from_axis_angle(axis=[0._R8P, 0._R8P, 1._R8P], angle=2 * PI / 3, translation=[1.8_R8P, 0.6_R8P, 1.2_R8P])
half = pose_interpolate(p0, p1, 0.5_R8P)
twice = pose_compose(half, half)
print '(A,ES9.2,A,F7.4)', 'compose(half, half): rotation off p1 by ', maxval(abs(twice%r - p1%r)), &
                          ', translation off p1 by ', norm2(twice%t - p1%t)
q = pose_interpolate(p0, p1, 1._R8P)
print '(A,L1)', 'pose_interpolate at s = 1 is p1 bit for bit: ', all(q%r == p1%r) .and. all(q%t == p1%t)
corner = [0.5_R8P, -0.5_R8P, 0.5_R8P]
write(caption, '(A,F5.3,A)') 'screw motion, s = 0 to 1 by 1/5; compose(half, half) (teal) misses p1 by ', &
                             norm2(twice%t - p1%t), ' in position'
call fig_open(u, 'ref-moving-bodies-screw', trim(caption), -62._R8P, 28._R8P, width=760_I4P, height=480_I4P)
do k = 0, NS
   q = pose_interpolate(p0, p1, real(k, R8P) / NS)
   call cube_edges(q)
   call fig_segments(u, merge('accent', 'box   ', k == NS), a, b)
enddo
call cube_edges(twice)
call fig_segments(u, 'probe', a, b)
do k = 1, NT
   q = pose_interpolate(p0, p1, real(k - 1, R8P) / (NT - 1))
   x = to_world(q, corner) ; path_a(k) = vector_R8P(x(1), x(2), x(3))
   x = to_world(q, [0._R8P, 0._R8P, 0._R8P]) ; path_b(k) = vector_R8P(x(1), x(2), x(3))
enddo
call fig_segments(u, 'hole', path_a(1:NT - 1), path_a(2:NT))
call fig_segment(u, 'probe', path_b(1), path_b(NT))
call fig_points(u, 'probe', [path_b(1), path_b(NT)])
call fig_points(u, 'hole', [path_a(1), path_a(NT)])
call fig_label(u, 'accent', vector_R8P(0._R8P, -0.75_R8P, -0.6_R8P), 's = 0')
call fig_label(u, 'accent', vector_R8P(1.8_R8P, 0.6_R8P, 1.95_R8P), 's = 1')
call fig_close(u)
contains
   subroutine cube_edges(pose)
   !< The 12 edges of the cube [-1/2, 1/2]^3 at a pose, into a, b.
   type(rigid_pose_t), intent(in) :: pose
   real(R8P)                      :: c(3, 8), w(3)
   integer(I4P)                   :: i, j, m

   do i = 0, 7
      w = to_world(pose, [merge(0.5_R8P, -0.5_R8P, btest(i, 0)), merge(0.5_R8P, -0.5_R8P, btest(i, 1)), &
                          merge(0.5_R8P, -0.5_R8P, btest(i, 2))])
      c(:, i + 1) = w
   enddo
   m = 0
   do i = 0, 7
      do j = 0, 2
         if (btest(i, j)) cycle
         m = m + 1
         a(m) = vector_R8P(c(1, i + 1), c(2, i + 1), c(3, i + 1))
         b(m) = vector_R8P(c(1, ibset(i, j) + 1), c(2, ibset(i, j) + 1), c(3, ibset(i, j) + 1))
      enddo
   enddo
   endsubroutine cube_edges
endprogram fig_ref_moving_bodies_screw
