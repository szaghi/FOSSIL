!run figref-surface-loads-bodies fig_ref_surface_loads_bodies
!fig ref-surface-loads-bodies
!< Figure of the reference page "Surface loads", several bodies (`load_samples_2d_bodies`): two circles of radius 0.25,
!< 0.03 apart, two probes 0.02 and 0.04 out of every sample: the samples (gold), the probes visible from their sample
!< (teal segments), the probes the other circle blocks (ember segments), each gap sample's clearance to the other circle
!< (grey segment to its second foot).
program fig_ref_surface_loads_bodies
use fossil,     only : curve_object, load_samples_t, load_samples_2d_bodies, rigid_pose_2d_t, pose_2d, to_world
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter          :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, H = 0.02_R8P
type(curve_object)            :: discs(2)
type(load_samples_t)          :: ls
type(rigid_pose_2d_t)         :: poses(2)
type(vector_R8P), allocatable :: a(:), b(:), c(:), d(:), e(:), f(:), p1(:)
real(R8P)                     :: u(200), v(200), w(2)
integer(I4P)                  :: un, i, s, k, dd

u = [(0.25_R8P * cos(2 * PI * (i - 1) / 200), i = 1, 200)] ; v = [(0.25_R8P * sin(2 * PI * (i - 1) / 200), i = 1, 200)]
call discs(1)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P]) ; call discs(2)%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
poses = [pose_2d(0._R8P, [-0.265_R8P, 0._R8P]), pose_2d(0._R8P, [0.265_R8P, 0._R8P])]
call load_samples_2d_bodies(discs, samples=ls, poses=poses, order=1_I4P, lengths=[H, 2 * H])
call fig_open(un, 'ref-surface-loads-bodies', 'two circles 0.03 apart: samples (gold), visible probes (teal), probes the '// &
              'other body blocks (ember), clearance (grey)', -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0._R8P, 0._R8P, 0._R8P), 0.12_R8P)
allocate(a(0), b(0))
do dd = 1, 2
   do i = 1, 200
      w = to_world(poses(dd), [u(i), v(i)]) ; a = [a, vector_R8P(w(1), w(2), 0._R8P)]
      w = to_world(poses(dd), [u(mod(i, 200) + 1), v(mod(i, 200) + 1)]) ; b = [b, vector_R8P(w(1), w(2), 0._R8P)]
   enddo
enddo
call fig_segments(un, 'accent', a, b)
deallocate(a, b) ; allocate(a(0), b(0), c(0), d(0), e(0), f(0), p1(0))
do s = 1, ls%n
   p1 = [p1, vector_R8P(ls%points(1, s), ls%points(2, s), 3 * LIFT)]
   do k = ls%m, 1, -1 ! the longer probe first: its segment covers the shorter one's
      if (ls%visible(k, s)) then
         a = [a, vector_R8P(ls%points(1, s), ls%points(2, s), LIFT)]
         b = [b, vector_R8P(ls%probe(1, k, s), ls%probe(2, k, s), LIFT)]
      else
         c = [c, vector_R8P(ls%points(1, s), ls%points(2, s), 2 * LIFT)]
         d = [d, vector_R8P(ls%probe(1, k, s), ls%probe(2, k, s), 2 * LIFT)]
      endif
   enddo
   if (ls%clearance(s) < 2 * H) then
      e = [e, vector_R8P(ls%points(1, s), ls%points(2, s), LIFT)]
      f = [f, vector_R8P(ls%second_foot(1, s), ls%second_foot(2, s), LIFT)]
   endif
enddo
call fig_segments(un, 'probe', a, b) ; call fig_segments(un, 'patch', c, d) ; call fig_segments(un, 'ghost', e, f)
call fig_points(un, 'foot', p1)
call fig_close(un)
print '(A,I0,A,I0,A,I0,A,F7.4)', 'bodies: samples ', ls%n, ', probes blocked ', count(.not. ls%visible), &
      ', within 2 h of the other body ', count(ls%clearance < 2 * H), ', least clearance ', minval(ls%clearance)
endprogram fig_ref_surface_loads_bodies
