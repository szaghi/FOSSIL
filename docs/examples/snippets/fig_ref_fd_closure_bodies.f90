!< Figure of the reference page "Finite-difference wall closure", several bodies (`fd_closure_2d_bodies`): two circles of
!< radius 0.25, 0.03 apart, on cell-centred lines of spacing 0.02, reach 2: the ghosts (teal), the shared ones (ember),
!< each shared ghost's foot (gold) joined to the nearest point of the other circle, its second foot (grey segment).
program fig_ref_fd_closure_bodies
use fossil,     only : curve_object, ib_closure_t, fd_closure_2d_bodies, rigid_pose_2d_t, pose_2d, to_world
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, H = 0.02_R8P
type(curve_object)            :: discs(2)
type(ib_closure_t)            :: cl
type(rigid_pose_2d_t)         :: poses(2)
type(vector_R8P), allocatable :: a(:), b(:), p1(:), p2(:), p3(:)
real(R8P)                     :: u(200), v(200), xc(60), yc(30), w(2)
integer(I4P)                  :: un, i, j, g, d

u = [(0.25_R8P * cos(2 * PI * (i - 1) / 200), i = 1, 200)] ; v = [(0.25_R8P * sin(2 * PI * (i - 1) / 200), i = 1, 200)]
call discs(1)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P]) ; call discs(2)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
xc = [(-0.6_R8P + (i - 0.5_R8P) * H, i = 1, 60)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * H, j = 1, 30)]
poses = [pose_2d(0._R8P, [-0.265_R8P, 0._R8P]), pose_2d(0._R8P, [0.265_R8P, 0._R8P])]
call fd_closure_2d_bodies(discs, xc, yc, [2_I4P, 2_I4P], cl, poses=poses)
call fig_open(un, 'ref-fd-closure-bodies', 'two circles 0.03 apart: ghosts (teal), shared ghosts (ember), feet (gold) '// &
              'joined to the other body', -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0._R8P, 0._R8P, 0._R8P), 0.16_R8P)
allocate(a(0))
do j = 1, size(yc) ; do i = 1, size(xc)
   a = [a, vector_R8P(xc(i), yc(j), 0._R8P)]
enddo ; enddo
call fig_points(un, 'box', a)
deallocate(a) ; allocate(a(0), b(0))
do d = 1, 2
   do i = 1, 200
      w = to_world(poses(d), [u(i), v(i)]) ; a = [a, vector_R8P(w(1), w(2), 0._R8P)]
      w = to_world(poses(d), [u(mod(i, 200) + 1), v(mod(i, 200) + 1)]) ; b = [b, vector_R8P(w(1), w(2), 0._R8P)]
   enddo
enddo
call fig_segments(un, 'accent', a, b)
deallocate(a, b) ; allocate(a(0), b(0), p1(0), p2(0), p3(0))
do g = 1, cl%n
   if (cl%shared(g)) then
      p2 = [p2, vector_R8P(cl%points(1, g), cl%points(2, g), 2 * LIFT)]
      p3 = [p3, vector_R8P(cl%foot(1, g), cl%foot(2, g), 2 * LIFT)]
      a = [a, vector_R8P(cl%foot(1, g), cl%foot(2, g), LIFT)]
      b = [b, vector_R8P(cl%second_foot(1, g), cl%second_foot(2, g), LIFT)]
   else
      p1 = [p1, vector_R8P(cl%points(1, g), cl%points(2, g), 2 * LIFT)]
   endif
enddo
call fig_segments(un, 'ghost', a, b)
call fig_points(un, 'probe', p1) ; call fig_points(un, 'patch', p2) ; call fig_points(un, 'foot', p3)
call fig_close(un)
print '(A,I0,A,I0,A,F7.4,A,F7.4)', 'bodies: ghosts ', cl%n, ', shared ', count(cl%shared), ', clearance of the shared ', &
      minval(cl%clearance, mask=cl%shared), ' to ', maxval(cl%clearance, mask=cl%shared)
endprogram fig_ref_fd_closure_bodies
