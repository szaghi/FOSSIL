!< Run and figure of the reference page "Finite-difference wall closure", section "Wall curvature":
!< - 2-D: an ellipse (semi-axes 0.4, 0.2, 400 vertices) on cell-centred lines of spacing 0.02, reach 2: at the foot of
!<   every ghost point the curvature of `query_many(principal_curvatures=)` against the exact curvature of the ellipse
!<   there; the figure zooms on the right end: the ghosts (teal) joined to their feet (gold), and at three feet the
!<   osculating circle of radius 1 / kappa (ember);
!< - 3-D: a spheroid (semi-axes 0.5, 0.5, 0.3, a unit icosphere of 5 subdivisions scaled) spinning about its axis: v = Omega x x is a rigid motion,
!<   stress-free, so at the wall d_n v_t = kappa_n v_t with kappa_n the normal curvature along v_t (the parallel); at the
!<   feet of the ghost points on node lines of spacing 0.05 (reach 2), `normal_curvature` along v_t against the exact
!<   n_rho / rho; the device query returns the same bits.
program fig_ref_fd_closure_curvature
use fossil,     only : surface_stl_object, curve_object, facet_object, ghost_points_t, &
                       ghost_points_2d, ghost_points_3d, normal_curvature
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P), A = 0.4_R8P, B = 0.2_R8P
type(curve_object)              :: ellipse
type(surface_stl_object)        :: spheroid
type(ghost_points_t)            :: gp
type(facet_object), allocatable :: facets(:)
type(vector_R8P),   allocatable :: q(:), s1(:), s2(:), pg(:), pf(:), pa(:)
real(R8P),          allocatable :: d(:), foot(:,:), nrm(:,:), kap(:,:), dirs(:,:,:), kap2(:,:), dirs2(:,:,:)
real(R8P)                       :: u(400), v(400), xc(50), yc(30), xn(25), zn(17), t, kx, err, rho, n(3), &
                                   vt(3), kn, c(3), r
integer(I4P)                    :: i, j, k, m, un, sel(3)

! ---- 2-D: the curvature at the feet of the ghost points of an ellipse
u = [(A * cos(2 * PI * (i - 1) / 400), i = 1, 400)] ; v = [(B * sin(2 * PI * (i - 1) / 400), i = 1, 400)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[400_I4P])
xc = [(-0.5_R8P + (i - 0.5_R8P) * 0.02_R8P, i = 1, 50)] ; yc = [(-0.3_R8P + (j - 0.5_R8P) * 0.02_R8P, j = 1, 30)]
call ellipse%build_curvature                                      ! the vertex curvatures, once per geometry
call ghost_points_2d(ellipse, xc, yc, [2_I4P, 2_I4P], gp)         ! the ghost points of a stencil of half-width 2
allocate(d(gp%n), foot(2, gp%n), nrm(2, gp%n), kap(1, gp%n), dirs(2, 1, gp%n))
call ellipse%query_many(points=gp%points, distances=d, closest_points=foot, normals=nrm, principal_curvatures=kap, &
                        principal_directions=dirs)                ! kap(1, g): the wall curvature at the foot of ghost g
err = 0._R8P
do i = 1, gp%n
   t = atan2(foot(2, i) / B, foot(1, i) / A)
   kx = A * B / (A**2 * sin(t)**2 + B**2 * cos(t)**2)**1.5_R8P
   err = max(err, abs(kap(1, i) - kx) / kx)
enddo
print '(A,I0,A,F7.3,A,F7.3,A,ES9.2)', 'ellipse: ghosts ', gp%n, ', curvature at the feet from ', minval(kap), ' to ', &
      maxval(kap), ', max relative error ', err
allocate(pg(0), pf(0), s1(0), s2(0), pa(0))
do i = 1, gp%n
   if (gp%points(1, i) < 0.22_R8P) cycle
   pg = [pg, vector_R8P(gp%points(1, i), gp%points(2, i), 1.e-3_R8P)]
   pf = [pf, vector_R8P(foot(1, i), foot(2, i), 2.e-3_R8P)]
enddo
do i = 1, 400
   s1 = [s1, vector_R8P(u(i), v(i), 0._R8P)] ; s2 = [s2, vector_R8P(u(mod(i, 400) + 1), v(mod(i, 400) + 1), 0._R8P)]
enddo
call fig_open(un, 'ref-fd-closure-curvature', 'ghosts (teal), feet (gold), osculating circles of radius 1/kappa (ember)', &
              -90._R8P, 90._R8P)
call fig_zoom(un, vector_R8P(0.34_R8P, 0._R8P, 0._R8P), 0.13_R8P)
deallocate(pa) ; allocate(pa(0))
do j = 1, 30 ; do i = 1, 50
   if (xc(i) > 0.2_R8P) pa = [pa, vector_R8P(xc(i), yc(j), 0._R8P)]
enddo ; enddo
call fig_points(un, 'box', pa)
call fig_segments(un, 'accent', s1, s2)
m = size(pg)
call fig_segments(un, 'box', pg, pf)
call fig_points(un, 'probe', pg)
call fig_points(un, 'foot', pf)
deallocate(s1, s2) ; allocate(s1(0), s2(0))
sel = [maxloc(kap(1, :), 1, mask=foot(1, :) > 0._R8P), minloc(abs(foot(1, :) - 0.33_R8P), 1, mask=foot(2, :) > 0._R8P), &
       minloc(abs(foot(1, :) - 0.27_R8P), 1, mask=foot(2, :) < 0._R8P)]
do m = 1, 3                                                       ! the osculating circles at three feet
   i = sel(m)
   r = 1._R8P / kap(1, i) ; c(1:2) = foot(:, i) - r * nrm(:, i)   ! centre inside, at 1 / kappa along -n
   do k = 0, 95
      s1 = [s1, vector_R8P(c(1) + r * cos(2 * PI * k / 96), c(2) + r * sin(2 * PI * k / 96), 3.e-3_R8P)]
      s2 = [s2, vector_R8P(c(1) + r * cos(2 * PI * (k + 1) / 96), c(2) + r * sin(2 * PI * (k + 1) / 96), 3.e-3_R8P)]
   enddo
   print '(A,2F8.4,A,F7.4)', 'osculating circle at the foot ', foot(:, i), ': radius ', r
enddo
call fig_segments(un, 'patch', s1, s2)
call fig_close(un)

! ---- 3-D: a spinning spheroid, the stress-free identity d_n v_t = kappa_n v_t at the feet of the ghost points
call icospheroid(5_I4P, [0.5_R8P, 0.5_R8P, 0.3_R8P], facets)    ! vertices on the spheroid: an icosphere scaled
call spheroid%adopt_facets(facets=facets)
xn = [(-0.6_R8P + 0.05_R8P * i, i = 0, 24)] ; zn = [(-0.4_R8P + 0.05_R8P * i, i = 0, 16)]
call ghost_points_3d(spheroid, xn, xn, zn, [2_I4P, 2_I4P, 2_I4P], gp)
deallocate(d, foot, kap, dirs) ; allocate(d(gp%n), q(gp%n), foot(3, gp%n), kap(2, gp%n), dirs(3, 2, gp%n), kap2(2, gp%n), dirs2(3, 2, gp%n))
q = [(vector_R8P(gp%points(1, i), gp%points(2, i), gp%points(3, i)), i = 1, gp%n)]
call spheroid%build_curvature                                     ! before enter_device: mapped with the surface
call spheroid%query_many(points=q, distances=d, closest_points=foot, principal_curvatures=kap, principal_directions=dirs)
do i = 1, gp%n
   vt = [-foot(2, i), foot(1, i), 0._R8P]                          ! the wall velocity of a spin about z, Omega = 1
   if (norm2(vt) > 0._R8P) kn = normal_curvature(kap(:, i), dirs(:, :, i), vt) ! then d_n v_t = kn v_t along v_t
enddo
err = 0._R8P
do i = 1, gp%n
   rho = hypot(foot(1, i), foot(2, i))
   if (rho < 0.05_R8P) cycle
   n = [foot(1, i) / 0.25_R8P, foot(2, i) / 0.25_R8P, foot(3, i) / 0.09_R8P] ; n = n / norm2(n) ! exact normal
   kn = normal_curvature(kap(:, i), dirs(:, :, i), [-foot(2, i), foot(1, i), 0._R8P])
   err = max(err, abs(kn - hypot(n(1), n(2)) / rho) / (hypot(n(1), n(2)) / rho))
enddo
call spheroid%enter_device
call spheroid%query_many_device(points=q, distances=d, principal_curvatures=kap2, principal_directions=dirs2)
call spheroid%exit_device
print '(A,I0,A,I0,A,ES9.2)', 'spheroid: facets ', spheroid%get_facets_number(), ', ghosts ', gp%n, &
      ', max relative error of kappa_n along v_t ', err
print '(A,L1)', 'device = host bit for bit: ', all(kap2 == kap) .and. all(dirs2 == dirs)

contains
   subroutine icospheroid(levels, semi, fo)
   !< A unit icosahedron subdivided `levels` times (midpoints pushed onto the unit sphere), scaled by the semi-axes: every
   !< vertex lies on the spheroid.
   integer(I4P),                    intent(in)  :: levels  !< Subdivisions.
   real(R8P),                       intent(in)  :: semi(3) !< Semi-axes.
   type(facet_object), allocatable, intent(out) :: fo(:)   !< Facets.
   integer(I4P), parameter :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                    2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                    4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                    5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   real(R8P), allocatable :: tri(:,:,:), nxt(:,:,:)
   real(R8P)              :: vv(3, 12), g
   integer(I4P)           :: k, l, m

   g = (1._R8P + sqrt(5._R8P)) / 2._R8P
   vv = reshape([-1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, 0._R8P, &
                 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, &
                 g, 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, -g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P], [3, 12])
   do k = 1, 12
      vv(:, k) = vv(:, k) / norm2(vv(:, k))
   enddo
   allocate(tri(3, 3, 20))
   do k = 1, 20
      tri(:, :, k) = vv(:, IDX(:, k))
   enddo
   do l = 1, levels
      m = size(tri, dim=3)
      allocate(nxt(3, 3, 4 * m))
      do k = 1, m
         nxt(:, :, 4 * k - 3) = reshape([tri(:, 1, k), mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
         nxt(:, :, 4 * k - 2) = reshape([mid(tri(:, 1, k), tri(:, 2, k)), tri(:, 2, k), mid(tri(:, 2, k), tri(:, 3, k))], [3, 3])
         nxt(:, :, 4 * k - 1) = reshape([mid(tri(:, 3, k), tri(:, 1, k)), mid(tri(:, 2, k), tri(:, 3, k)), tri(:, 3, k)], [3, 3])
         nxt(:, :, 4 * k)     = reshape([mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 2, k), tri(:, 3, k)), &
                                         mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   m = size(tri, dim=3)
   allocate(fo(m))
   do k = 1, m
      do l = 1, 3
         fo(k)%vertex(l) = vector_R8P(semi(1) * tri(1, l, k), semi(2) * tri(2, l, k), semi(3) * tri(3, l, k))
      enddo
      call fo(k)%compute_metrix
   enddo
   endsubroutine icospheroid

   pure function mid(a, b) result(m)
   !< Unit midpoint.
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: m(3)

   m = (a + b) / norm2(a + b)
   endfunction mid
endprogram fig_ref_fd_closure_curvature
