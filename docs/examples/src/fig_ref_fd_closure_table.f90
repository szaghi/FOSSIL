!run figref-fd-closure-table fig_ref_fd_closure_table
!fig ref-fd-closure-table
!fig ref-fd-closure-ilw
!fig ref-fd-closure-velocity
!fig ref-fd-closure-image-line
!fig ref-fd-closure-resident
!< Figures of the reference page "Finite-difference wall closure", the closure table (`fd_closure_2d`, `fd_closure_3d`):
!< - table: the right end of an ellipse (semi-axes 0.4, 0.2) on cell-centred lines of spacing 0.02, reach 2: each ghost
!<   joined to its foot, the image line of three probes along the normal from the foot;
!< - ilw: the same end, each ghost joined along its grid lines to the ILW crossing of every direction of its mask;
!< - velocity: the ellipse at a pose, spinning about the pose's origin and moving along u: the wall velocity at the feet;
!< - image line: one ghost, its image line, and the solver's cells that enclose the three probes (`probe_cell`);
!< - resident: a marching-cubes sphere on node lines, the table built resident on the device: the feet and normals of the
!<   ghosts of its lower half.
program fig_ref_fd_closure_table
use fossil,     only : surface_stl_object, curve_object, facet_object, extract_isosurface, ib_closure_t, fd_closure_2d, &
                       fd_closure_3d, ib_closure_reset, rigid_pose_2d_t, pose_2d, to_world
use fossil_fig, only : fig_open, fig_zoom, fig_quad, fig_segments, fig_points, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P, H = 0.02_R8P
type(curve_object)              :: ellipse
type(surface_stl_object)        :: sphere
type(ib_closure_t)              :: cl
type(rigid_pose_2d_t)           :: pose
type(facet_object), allocatable :: facets(:)
type(vector_R8P),   allocatable :: a(:), b(:), p1(:), p2(:), p3(:)
real(R8P)                       :: u(400), v(400), xc(50), yc(30), xn(25), f(33, 33, 33), w(2), c(2)
integer(I4P)                    :: un, i, j, k, g, d, ax, ii, jj, nvis, nblk

u = [(0.4_R8P * cos(2 * PI * (i - 1) / 400), i = 1, 400)] ; v = [(0.2_R8P * sin(2 * PI * (i - 1) / 400), i = 1, 400)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[400_I4P])
call ellipse%build_curvature
xc = [(-0.5_R8P + (i - 0.5_R8P) * H, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * H, j = 1, 30)]

! ---- the table: feet and image lines
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], cl, lengths=[(k * sqrt(2._R8P) * H, k = 1, 3)], ilw=.true.)
call fig_open(un, 'ref-fd-closure-table', 'ghosts (teal), feet (gold), image lines of three probes (ember)', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vec(0.33_R8P, 0._R8P), 0.12_R8P)
call grid_points(a) ; call fig_points(un, 'box', a)
call curve_segments(u, v, a, b, 0._R8P) ; call fig_segments(un, 'accent', a, b)
allocate(p1(0), p2(0), p3(0)) ; deallocate(a, b) ; allocate(a(0), b(0))
nvis = 0 ; nblk = 0
do g = 1, cl%n
   if (cl%points(1, g) < 0.2_R8P) cycle
   p1 = [p1, vec(cl%points(1, g), cl%points(2, g), 2 * LIFT)] ; p2 = [p2, vec(cl%foot(1, g), cl%foot(2, g), 2 * LIFT)]
   a = [a, vec(cl%points(1, g), cl%points(2, g), LIFT)] ; b = [b, vec(cl%foot(1, g), cl%foot(2, g), LIFT)]
   a = [a, vec(cl%foot(1, g), cl%foot(2, g), LIFT)] ; b = [b, vec(cl%probe(1, cl%m, g), cl%probe(2, cl%m, g), LIFT)]
   do k = 1, cl%m
      p3 = [p3, vec(cl%probe(1, k, g), cl%probe(2, k, g), 2 * LIFT)]
      if (cl%visible(k, g)) then ; nvis = nvis + 1 ; else ; nblk = nblk + 1 ; endif
   enddo
enddo
call fig_segments(un, 'ghost', a, b)
call fig_points(un, 'probe', p1) ; call fig_points(un, 'foot', p2) ; call fig_points(un, 'patch', p3)
call fig_close(un)
print '(A,I0,A,I0,A,I0)', 'table: ghosts ', cl%n, ', probes drawn visible ', nvis, ', blocked ', nblk

! ---- ILW: along the grid lines from each ghost to its crossing
call fig_open(un, 'ref-fd-closure-ilw', 'ghosts (teal), walks to the ILW crossings (gold), normals (ember)', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vec(0.36_R8P, 0.06_R8P), 0.07_R8P)
call grid_points(a) ; call fig_points(un, 'box', a)
call curve_segments(u, v, a, b, 0._R8P) ; call fig_segments(un, 'accent', a, b)
deallocate(a, b, p1, p2, p3) ; allocate(a(0), b(0), p1(0), p2(0), p3(0))
do g = 1, cl%n
   if (cl%points(1, g) < 0.28_R8P .or. cl%points(2, g) < 0._R8P) cycle
   p1 = [p1, vec(cl%points(1, g), cl%points(2, g), 2 * LIFT)]
   do d = 1, 4
      if (cl%ilw_element(d, g) < 1) cycle
      ax = (d + 1) / 2
      if (ax == 1) then
         c = [cl%ilw_t(d, g), cl%points(2, g)]
      else
         c = [cl%points(1, g), cl%ilw_t(d, g)]
      endif
      a = [a, vec(cl%points(1, g), cl%points(2, g), LIFT)] ; b = [b, vec(c(1), c(2), LIFT)]
      p2 = [p2, vec(c(1), c(2), 2 * LIFT)]
      p3 = [p3, vec(c(1), c(2), LIFT)]
   enddo
enddo
call fig_segments(un, 'ghost', a, b)
deallocate(a, b) ; allocate(a(0), b(0))
do g = 1, cl%n
   if (cl%points(1, g) < 0.28_R8P .or. cl%points(2, g) < 0._R8P) cycle
   do d = 1, 4
      if (cl%ilw_element(d, g) < 1) cycle
      ax = (d + 1) / 2
      if (ax == 1) then
         c = [cl%ilw_t(d, g), cl%points(2, g)]
      else
         c = [cl%points(1, g), cl%ilw_t(d, g)]
      endif
      a = [a, vec(c(1), c(2), LIFT)] ; b = [b, vec(c(1) + 0.4_R8P * H * cl%ilw_normal(1, d, g), &
                                                   c(2) + 0.4_R8P * H * cl%ilw_normal(2, d, g), LIFT)]
   enddo
enddo
call fig_segments(un, 'patch', a, b)
call fig_points(un, 'probe', p1) ; call fig_points(un, 'foot', p2)
call fig_close(un)
print '(A,I0,A,I0)', 'ILW: slots with a crossing ', count(cl%ilw_element > 0), ', ghosts ', cl%n

! ---- velocity: the ellipse at a pose, spinning
pose = pose_2d(0.3_R8P, [0.05_R8P, 0._R8P])
call fd_closure_2d(ellipse, xc, yc, [1_I4P, 1_I4P], cl, pose=pose, linear_velocity=[0.5_R8P, 0._R8P], &
                   angular_velocity=2._R8P)
call fig_open(un, 'ref-fd-closure-velocity', 'the ellipse at its pose: wall velocity at the feet (ember)', &
              -90._R8P, 90._R8P)
deallocate(a, b) ; allocate(a(0), b(0))
do i = 1, 400
   w = to_world(pose, [u(i), v(i)]) ; a = [a, vec(w(1), w(2), 0._R8P)]
   w = to_world(pose, [u(mod(i, 400) + 1), v(mod(i, 400) + 1)]) ; b = [b, vec(w(1), w(2), 0._R8P)]
enddo
call fig_segments(un, 'accent', a, b)
deallocate(a, b, p2) ; allocate(a(0), b(0), p2(0))
do g = 1, cl%n, 2
   a = [a, vec(cl%foot(1, g), cl%foot(2, g), LIFT)]
   b = [b, vec(cl%foot(1, g) + 0.08_R8P * cl%velocity(1, g), cl%foot(2, g) + 0.08_R8P * cl%velocity(2, g), LIFT)]
   p2 = [p2, vec(cl%foot(1, g), cl%foot(2, g), 2 * LIFT)]
enddo
call fig_segments(un, 'patch', a, b) ; call fig_points(un, 'foot', p2)
call fig_points(un, 'probe', [vec(pose%t(1), pose%t(2), 2 * LIFT)])
call fig_close(un)
print '(A,I0,A,F7.4,A,F7.4)', 'velocity: ghosts ', cl%n, ', speed from ', minval(norm2(cl%velocity, 1)), ' to ', &
      maxval(norm2(cl%velocity, 1))

! ---- image line: one ghost, its probes and the cells enclosing them
call fd_closure_2d(ellipse, xc, yc, [2_I4P, 2_I4P], cl, lengths=[(k * sqrt(2._R8P) * H, k = 1, 3)])
g = minloc(abs(cl%points(1, :) - 0.31_R8P) + abs(cl%points(2, :) - 0.11_R8P), 1)
call fig_open(un, 'ref-fd-closure-image-line', 'a ghost (teal), its foot (gold), image line (ember), probe cells', -90._R8P, 90._R8P)
call fig_zoom(un, vec(cl%probe(1, 2, g), cl%probe(2, 2, g)), 0.07_R8P)
do k = 1, cl%m
   ii = cl%probe_cell(1, k, g) ; jj = cl%probe_cell(2, k, g)
   if (ii < 1 .or. jj < 1) cycle
   call fig_quad(un, 'frac', vec(xc(ii), yc(jj), -LIFT), vec(xc(ii + 1), yc(jj), -LIFT), vec(xc(ii + 1), yc(jj + 1), -LIFT), &
                 vec(xc(ii), yc(jj + 1), -LIFT), 0.6_R8P)
   print '(A,I0,A,2F9.5,A,2I3)', 'probe ', k, ' at', cl%probe(:, k, g), ' in cell', ii, jj
enddo
call grid_points(a) ; call fig_points(un, 'box', a)
call curve_segments(u, v, a, b, 0._R8P) ; call fig_segments(un, 'accent', a, b)
call fig_segments(un, 'ghost', [vec(cl%points(1, g), cl%points(2, g), LIFT)], &
                  [vec(cl%probe(1, cl%m, g), cl%probe(2, cl%m, g), LIFT)])
call fig_points(un, 'probe', [vec(cl%points(1, g), cl%points(2, g), 2 * LIFT)])
call fig_points(un, 'foot', [vec(cl%foot(1, g), cl%foot(2, g), 2 * LIFT)])
deallocate(p3) ; allocate(p3(0))
do k = 1, cl%m
   p3 = [p3, vec(cl%probe(1, k, g), cl%probe(2, k, g), 2 * LIFT)]
enddo
call fig_points(un, 'patch', p3)
call fig_close(un)

! ---- resident: a marching-cubes sphere, the table built on the device
do k = 1, 33 ; do j = 1, 33 ; do i = 1, 33
   f(i, j, k) = norm2([i - 17._R8P, j - 17._R8P, k - 17._R8P] / 16._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)]
call sphere%enter_device
call fd_closure_3d(sphere, xn, xn, xn, [2_I4P, 2_I4P, 2_I4P], cl, lengths=[0.05_R8P], resident=.true.)
call cl%update_host
call fig_open(un, 'ref-fd-closure-resident', 'resident table: feet (gold) and normals (ember), lower half', -60._R8P, 25._R8P)
call fig_surface(un, 'ghost', sphere)
deallocate(a, b, p2) ; allocate(a(0), b(0), p2(0))
do g = 1, cl%n
   if (cl%foot(3, g) > 0._R8P) cycle
   p2 = [p2, vector_R8P(cl%foot(1, g), cl%foot(2, g), cl%foot(3, g))]
   a = [a, vector_R8P(cl%foot(1, g), cl%foot(2, g), cl%foot(3, g))]
   b = [b, vector_R8P(cl%probe(1, 1, g), cl%probe(2, 1, g), cl%probe(3, 1, g))]
enddo
call fig_segments(un, 'patch', a, b) ; call fig_points(un, 'foot', p2)
call fig_close(un)
print '(A,I0,A,L1)', 'sphere: ghosts ', cl%n, ', resident ', cl%on_device
call ib_closure_reset(cl)
call sphere%exit_device

contains
   function vec(x, y, z) result(p)
   !< A point of the drawing plane.
   real(R8P), intent(in)           :: x, y !< Coordinates.
   real(R8P), intent(in), optional :: z    !< Lift.
   type(vector_R8P)                :: p    !< Point.

   p = vector_R8P(x, y, 0._R8P)
   if (present(z)) p%z = z
   endfunction vec

   subroutine grid_points(a)
   !< The solver's cell centres.
   type(vector_R8P), allocatable, intent(out) :: a(:) !< Points.
   integer(I4P) :: i, j                                !< Counters.

   allocate(a(0))
   do j = 1, size(yc) ; do i = 1, size(xc)
      a = [a, vec(xc(i), yc(j))]
   enddo ; enddo
   endsubroutine grid_points

   subroutine curve_segments(u, v, a, b, z)
   !< The segments of a closed polyline.
   real(R8P),                     intent(in)  :: u(:), v(:) !< Vertices.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Segment ends.
   real(R8P),                     intent(in)  :: z          !< Lift.
   integer(I4P) :: i, n                                     !< Counters.

   n = size(u) ; allocate(a(n), b(n))
   do i = 1, n
      a(i) = vec(u(i), v(i), z) ; b(i) = vec(u(mod(i, n) + 1), v(mod(i, n) + 1), z)
   enddo
   endsubroutine curve_segments
endprogram fig_ref_fd_closure_table
