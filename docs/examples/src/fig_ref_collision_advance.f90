!run figref-collision-advance fig_ref_collision_advance
!fig ref-collision-advance
!< Figure of the reference page "Continuous collision detection", section "Conservative advancement": a square (side
!< 0.4, frame at its centre) moves by (1.9, 2) and turns a quarter over the step towards a wall at x = 2.2. Top:
!< the distance d(s) to the wall (gold, 400 samples of distance_to at the interpolated poses), the iterates s_k of
!< curve%time_of_impact (ember) and from each the Lipschitz line d_k - mu (s - s_k) (teal), which stays under d(s)
!< and reaches tol / 2 at s_(k+1). Bottom: the square at s_0 to s_3 and at contact (ember).
program fig_ref_collision_advance
use fossil,     only : curve_object, rigid_pose_2d_t, pose_2d, pose_interpolate, to_world, motion_bound_2d, CCD_CONTACT
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter :: NS = 400, KMAX = 60
real(R8P),    parameter :: PI = acos(-1._R8P), XS = 5._R8P, YS = 1.6_R8P, TOL = 1.e-6_R8P
real(R8P),    parameter :: UA(4) = [-0.2_R8P, 0.2_R8P, 0.2_R8P, -0.2_R8P], VA(4) = [-0.2_R8P, -0.2_R8P, 0.2_R8P, 0.2_R8P]
real(R8P),    parameter :: UB(4) = [2.2_R8P, 2.5_R8P, 2.5_R8P, 2.2_R8P], VB(4) = [-1.4_R8P, -1.4_R8P, 1.4_R8P, 1.4_R8P]
type(curve_object)      :: a, b
type(rigid_pose_2d_t)   :: a0, a1, rest
type(vector_R8P)        :: c0(NS), c1(NS)
real(R8P)               :: d(0:NS), s(0:KMAX), dk(0:KMAX), mu, toi, y0
integer(I4P)            :: i, k, n, u, status, it, m
character(100)          :: caption

call a%adopt_loops(u=UA, v=VA, loop_sizes=[4])
call b%adopt_loops(u=UB, v=VB, loop_sizes=[4])
rest = pose_2d(angle=0._R8P)
a0 = pose_2d(angle=0._R8P, translation=[0.3_R8P, -1._R8P])
a1 = pose_2d(angle=PI / 2, translation=[2.2_R8P, 1._R8P])
mu = motion_bound_2d(a0, a1, UA, VA)
do i = 0, NS
   call a%distance_to(b, d(i), pose_self=pose_interpolate(a0, a1, real(i, R8P) / NS), pose_other=rest)
enddo
s(0) = 0._R8P ; dk(0) = d(0) ; n = KMAX
do k = 1, KMAX                                       ! the k-th iterate: the loop capped at k steps
   call a%time_of_impact(b, a0, a1, rest, rest, toi, status, tolerance=TOL, distance=dk(k), iterations=it, &
                         max_iterations=k)
   s(k) = toi
   if (status == CCD_CONTACT) then
      n = it ; exit
   endif
enddo
print '(A,F8.5,A,F8.5,A,I0,A,F9.6)', 'mu ', mu, ', d(0) ', d(0), ', contact after ', n, ' steps at s = ', s(n)
do k = 0, min(n, 6_I4P)
   print '(A,I2,A,F9.6,A,ES11.4)', 'k ', k, ': s ', s(k), ', d ', dk(k)
enddo

write(caption, '(A,F5.3,A,I0,A,F6.4)') 'd(s) (gold) and the bounds d_k - mu (s - s_k), mu = ', mu, ': contact in ', &
                                       n, ' steps at s = ', s(n)
call fig_open(u, 'ref-collision-advance', trim(caption), -90._R8P, 90._R8P, width=760_I4P, height=560_I4P)
y0 = 0._R8P
call fig_segments(u, 'box', [g(0._R8P, 0._R8P), g(0._R8P, 0._R8P)], [g(1._R8P, 0._R8P), g(0._R8P, 1.1_R8P * d(0))])
do i = 1, NS
   c0(i) = g(real(i - 1, R8P) / NS, d(i - 1)) ; c1(i) = g(real(i, R8P) / NS, d(i))
enddo
call fig_segments(u, 'accent', c0, c1)
do k = 0, n - 1                                      ! the Lipschitz line from (s_k, d_k) down to (s_k+1, tol / 2)
   call fig_segments(u, 'probe', [g(s(k), dk(k))], [g(s(k + 1), dk(k) - mu * (s(k + 1) - s(k)))])
   call fig_segments(u, 'box', [g(s(k + 1), 0._R8P)], [g(s(k + 1), dk(k + 1))])
enddo
call fig_points(u, 'hole', [(g(s(k), dk(k)), k = 0, n)])
call fig_label(u, 'box', g(1._R8P, 0._R8P), 's = 1')
call fig_label(u, 'box', g(0._R8P, 1.1_R8P * d(0)), 'd')
call fig_label(u, 'probe', g(s(0) + 0.06_R8P, dk(0) - mu * 0.06_R8P), 'slope -mu')
do k = 0, min(n, 3_I4P)                              ! the bodies at the first iterates and at contact, below the plot
   call fig_segments(u, 'ghost', [(w(s(k), m), m = 1, 4)], [(w(s(k), mod(m, 4) + 1), m = 1, 4)])
   call fig_label(u, 'ghost', w(s(k), 1_I4P), 's'//achar(48 + k))
enddo
call fig_segments(u, 'patch', [(w(s(n), m), m = 1, 4)], [(w(s(n), mod(m, 4) + 1), m = 1, 4)])
call fig_segments(u, 'body', [(wb(m), m = 1, 4)], [(wb(mod(m, 4) + 1), m = 1, 4)])
call fig_label(u, 'accent', wb(1), 'wall')
call fig_close(u)
contains
   function g(ss, dd) result(q)
   !< A point of the graph: s across, d up.
   real(R8P), intent(in) :: ss, dd !< Instant, distance.
   type(vector_R8P)      :: q      !< Point.

   q = vector_R8P(XS * ss, YS * dd, 0._R8P)
   endfunction g

   function w(ss, m) result(q)
   !< Vertex m of the square at instant ss, placed below the graph.
   real(R8P),    intent(in) :: ss !< Instant.
   integer(I4P), intent(in) :: m  !< Vertex.
   type(vector_R8P)         :: q  !< Point.
   real(R8P)                :: p(2)

   p = to_world(pose_interpolate(a0, a1, ss), [UA(m), VA(m)])
   q = vector_R8P(p(1) + 1.2_R8P, p(2) - 2.4_R8P, 0._R8P)
   endfunction w

   function wb(m) result(q)
   !< Vertex m of the wall, placed below the graph.
   integer(I4P), intent(in) :: m !< Vertex.
   type(vector_R8P)         :: q !< Point.

   q = vector_R8P(UB(m) + 1.2_R8P, VB(m) - 2.4_R8P, 0._R8P)
   endfunction wb
endprogram fig_ref_collision_advance
