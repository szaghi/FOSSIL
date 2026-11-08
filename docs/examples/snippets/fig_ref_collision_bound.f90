!< Figure of the reference page "Continuous collision detection": one rectangle, the same placements at s = 0 and s = 1
!< (a quarter turn and a shift), two body frames. Left: the frame origin at the rectangle's centre; right: 1.2 to its
!< left. pose_interpolate turns about the frame origin while the origin moves on a straight line, so the two motions
!< differ in between: the vertex paths (gold), the origin's path (teal), the rectangle at s = 0, 1/4, ..., 1, and the
!< vertex with the largest rho (ember), which sets the speed bound mu.
program fig_ref_collision_bound
use fossil,     only : rigid_pose_2d_t, pose_2d, pose_interpolate, to_world, motion_bound_2d
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter  :: PI = acos(-1._R8P)
real(R8P), parameter  :: UR(4) = [-0.5_R8P, 0.5_R8P, 0.5_R8P, -0.5_R8P], VR(4) = [-0.2_R8P, -0.2_R8P, 0.2_R8P, 0.2_R8P]
type(rigid_pose_2d_t) :: p0, p1
real(R8P)             :: c(2), rho(4), mu
integer(I4P)          :: run, u
character(40)         :: text

call fig_open(u, 'ref-collision-bound', 'same end poses, two frames: the origin (teal) sets the motion and mu', &
              -90._R8P, 90._R8P, height=360_I4P)
do run = 1, 2
   c = merge([0._R8P, 0._R8P], [1.2_R8P, 0._R8P], run == 1) ! the rectangle's centre in the body frame
   ! world placements: centre at (0, 0) turned by 0 at s = 0, centre at (2, 0.6) turned by 90 degrees at s = 1
   p0 = pose_2d(angle=0._R8P, translation=[0._R8P, 0._R8P] - c + [3.6_R8P * (run - 1), 0._R8P])
   p1 = pose_2d(angle=PI / 2, translation=[2._R8P + 3.6_R8P * (run - 1), 0.6_R8P] - [-c(2), c(1)])
   rho = hypot(UR + c(1), VR + c(2))
   mu = motion_bound_2d(p0, p1, UR + c(1), VR + c(2))
   print '(A,F4.1,A,4F8.4,A,F8.4)', 'centre at x = ', c(1), ' in the body frame: rho ', rho, ', mu ', mu
   call draw()
   write(text, '(A,F6.3)') merge('origin at the centre, mu ', 'origin 1.2 left, mu      ', run == 1), mu
   call fig_label(u, 'accent', vector_R8P(3.6_R8P * (run - 1) - 0.6_R8P, -0.9_R8P, 0._R8P), trim(text))
enddo
call fig_close(u)
contains
   subroutine draw()
   !< The rectangle at five instants, the vertex paths, the origin's path, the vertex of largest rho.
   type(vector_R8P)         :: a(64), b(64), pa(4 * 40), pb(4 * 40)
   integer(I4P)             :: k, i, m

   m = 0
   do k = 0, 4
      do i = 1, 4
         m = m + 1 ; a(m) = w(0.25_R8P * k, i) ; b(m) = w(0.25_R8P * k, mod(i, 4_I4P) + 1_I4P)
      enddo
   enddo
   call fig_segments(u, 'ghost', a(5:16), b(5:16))
   call fig_quad(u, 'body', w(0._R8P, 1_I4P), w(0._R8P, 2_I4P), w(0._R8P, 3_I4P), w(0._R8P, 4_I4P))
   call fig_quad(u, 'body', w(1._R8P, 1_I4P), w(1._R8P, 2_I4P), w(1._R8P, 3_I4P), w(1._R8P, 4_I4P))
   m = 0
   do i = 1, 4
      do k = 1, 40
         m = m + 1 ; pa(m) = w((k - 1) / 40._R8P, i) ; pb(m) = w(k / 40._R8P, i)
      enddo
   enddo
   call fig_segments(u, 'accent', pa(1:m), pb(1:m))
   call fig_segments(u, 'probe', [o(0._R8P)], [o(1._R8P)])
   call fig_points(u, 'probe', [o(0._R8P), o(1._R8P)])
   i = maxloc(rho, 1)
   call fig_points(u, 'hole', [(w(0.25_R8P * k, i), k = 0, 4)])
   endsubroutine draw

   function w(s, i) result(q)
   !< Vertex i of the rectangle at s, world frame.
   real(R8P),    intent(in) :: s !< Instant.
   integer(I4P), intent(in) :: i !< Vertex.
   type(vector_R8P)         :: q !< Point.
   real(R8P)                :: y(2)

   y = to_world(pose_interpolate(p0, p1, s), [UR(i) + c(1), VR(i) + c(2)])
   q = vector_R8P(y(1), y(2), 0._R8P)
   endfunction w

   function o(s) result(q)
   !< The body-frame origin at s, world frame.
   real(R8P), intent(in) :: s !< Instant.
   type(vector_R8P)      :: q !< Point.
   real(R8P)             :: y(2)

   y = to_world(pose_interpolate(p0, p1, s), [0._R8P, 0._R8P])
   q = vector_R8P(y(1), y(2), 0._R8P)
   endfunction o
endprogram fig_ref_collision_bound
