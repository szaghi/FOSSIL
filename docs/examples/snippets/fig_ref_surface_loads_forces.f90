!< Figure of the reference page "Surface loads", section "Forces and moments" (`integrate_loads`): a NACA 0012 of unit
!< chord (60 vertices per side, the upper and the lower side two patches) at 20 degrees nose up about its leading edge c,
!< in the pressure p = -v of a fluid at rest under gravity -v (zero at v = 0): the contributions -p n w of the samples
!< (ember, drawn ending on the wall), the force on each side drawn from c (teal), the total from the centroid (gold), all
!< three at one scale: the buoyancy (0, A), the moment about c the buoyancy's A (x_c - c) x (0, 1).
program fig_ref_surface_loads_forces
use fossil,     only : curve_object, load_samples_t, load_samples_2d, integrate_loads, rigid_pose_2d_t, pose_2d, to_world
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter         :: M = 60
real(R8P),    parameter         :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, SC = 2._R8P, SF = 2._R8P, SW = 6._R8P
type(curve_object)              :: naca
type(load_samples_t)            :: ls
type(rigid_pose_2d_t)           :: pose
type(vector_R8P), allocatable   :: a(:), b(:)
real(R8P),        allocatable   :: p(:), ppf(:,:), ppm(:)
real(R8P)                       :: u(2 * M), v(2 * M), w(2), q(2), force(2), moment, area, xc(2), xi
integer(I4P)                    :: i, k, s, un

k = 0
do i = M, 0, -1                                              ! trailing edge, upper side, leading edge
   xi = 0.5_R8P * (1._R8P - cos(PI * i / M)) ; k = k + 1 ; u(k) = xi ; v(k) = yt(xi)
enddo
do i = 1, M - 1                                              ! lower side
   xi = 0.5_R8P * (1._R8P - cos(PI * i / M)) ; k = k + 1 ; u(k) = xi ; v(k) = -yt(xi)
enddo
call naca%adopt_loops(u=u, v=v, loop_sizes=[2 * M], patches=[(1_I4P, i = 1, M), (2_I4P, i = M + 1, 2 * M)], &
                      patch_names=['upper', 'lower'])
pose = pose_2d(-20._R8P * PI / 180._R8P, [0._R8P, 0._R8P])   ! nose up about the leading edge
call load_samples_2d(naca, samples=ls, order=2_I4P, pose=pose)
allocate(p(ls%n))
do s = 1, ls%n
   p(s) = -ls%points(2, s)                                   ! a fluid at rest under gravity -v, p = 0 at v = 0
enddo
call integrate_loads(ls, force, moment, pressure=p, about=[0._R8P, 0._R8P], per_patch_force=ppf, &
                     per_patch_moment=ppm)
call naca%mass_properties(area=area, center=xc) ; xc = to_world(pose, xc)

call fig_open(un, 'ref-surface-loads-moment', 'NACA 0012 in p = -v: -p n w (ember), each side (teal), the total (gold)', &
              -90._R8P, 90._R8P)
allocate(a(2 * M), b(2 * M))
do i = 1, 2 * M
   w = to_world(pose, [u(i), v(i)]) ; q = to_world(pose, [u(mod(i, 2 * M) + 1), v(mod(i, 2 * M) + 1)])
   a(i) = vec(w(1), w(2)) ; b(i) = vec(q(1), q(2))
enddo
call fig_segments(un, 'body', a, b)
deallocate(a, b) ; allocate(a(ls%n), b(ls%n))
do s = 1, ls%n
   a(s) = vec(ls%points(1, s) + SW * p(s) * ls%normal(1, s) * ls%weight(s), &   ! drawn outside, ending on the wall
              ls%points(2, s) + SW * p(s) * ls%normal(2, s) * ls%weight(s), LIFT)
   b(s) = vec(ls%points(1, s), ls%points(2, s), LIFT)
enddo
call fig_segments(un, 'patch', a, b)
call fig_segments(un, 'probe', [vec(0._R8P, 0._R8P, 2 * LIFT), vec(0._R8P, 0._R8P, 2 * LIFT)], &
                  [vec(SC * ppf(1, 1), SC * ppf(2, 1), 2 * LIFT), vec(SC * ppf(1, 2), SC * ppf(2, 2), 2 * LIFT)])
call fig_segments(un, 'foot', [vec(xc(1), xc(2), 3 * LIFT)], [vec(xc(1) + SF * force(1), xc(2) + SF * force(2), 3 * LIFT)])
call fig_points(un, 'foot', [vec(0._R8P, 0._R8P, 3 * LIFT), vec(xc(1), xc(2), 3 * LIFT)])
call fig_label(un, 'foot', vec(-0.06_R8P, 0.02_R8P, 3 * LIFT), 'c')
call fig_label(un, 'foot', vec(xc(1) + 0.02_R8P, xc(2) - 0.05_R8P, 3 * LIFT), 'x_c')
call fig_label(un, 'probe', vec(SC * ppf(1, 1) + 0.02_R8P, SC * ppf(2, 1) - 0.03_R8P, 3 * LIFT), 'upper')
call fig_label(un, 'probe', vec(SC * ppf(1, 2) + 0.02_R8P, SC * ppf(2, 2) + 0.01_R8P, 3 * LIFT), 'lower')
call fig_label(un, 'foot', vec(xc(1) + SF * force(1) + 0.02_R8P, xc(2) + SF * force(2), 3 * LIFT), 'F')
call fig_close(un)
print '(A,I0,A,I0,A,F9.6)', 'NACA 0012 at 20 degrees: segments ', ls%ne, ', samples ', ls%n, ', area A ', area
print '(A,2F10.6,A,2F10.6,A,2F10.6)', 'force on the upper side', ppf(:, 1), ', lower side', ppf(:, 2), ', total', &
      anint(force * 1.e6_R8P) / 1.e6_R8P + 0._R8P
print '(A,F10.6,A,F10.6,A,F10.6)', 'moment about c: upper', ppm(1), ', lower', ppm(2), ', total', moment
print '(A,F10.6,A,F10.6,A,L1)', 'buoyancy A = ', area, ', A (x_c - c) x (0, 1) = ', area * xc(1), &
      ';  equal to 1e-14: ', abs(force(1)) < 1.e-14_R8P .and. abs(force(2) - area) < 1.e-14_R8P .and. &
      abs(moment - area * xc(1)) < 1.e-14_R8P

contains
   function vec(x, y, z) result(pt)
   !< A point of the drawing plane.
   real(R8P), intent(in)           :: x, y !< Coordinates.
   real(R8P), intent(in), optional :: z    !< Lift.
   type(vector_R8P)                :: pt   !< Point.

   pt = vector_R8P(x, y, 0._R8P)
   if (present(z)) pt%z = z
   endfunction vec

   pure real(R8P) function yt(x)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x ** 2 + 0.2843_R8P * x ** 3 - 0.1036_R8P * x ** 4)
   endfunction yt
endprogram fig_ref_surface_loads_forces
