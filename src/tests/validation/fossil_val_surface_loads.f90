!< FOSSIL validation: what accuracy the surface loads give (#109 WP2 step 06), through both paths: the load samples
!< (`load_samples_3d/_2d` + `integrate_loads`, the finite-difference path) and the cut cells' wall terms
!< (`cut_cells_3d(walls=)` / `cut_walls_2d` + `cut_loads_3d/_2d`, the finite-volume path). A quadrature study, no flow
!< solve: the pressure is given in closed form.
program fossil_val_surface_loads
!< FOSSIL validation: what accuracy the surface loads give (#109 WP2 step 06), through both paths: the load samples
!< (`load_samples_3d/_2d` + `integrate_loads`, the finite-difference path) and the cut cells' wall terms
!< (`cut_cells_3d(walls=)` / `cut_walls_2d` + `cut_loads_3d/_2d`, the finite-volume path). A quadrature study, no flow
!< solve: the pressure is given in closed form.
!<
!< 1. Buoyancy on large closed meshes. In p = p0 + g . x (tau = 0) the divergence theorem gives F = -g V and
!<    M = V (x_c - c) x (-g) on any closed faceted body (V, x_c from `mass_properties`; 2-D per unit depth, A and the
!<    scalar moment). Bodies: the STL fixtures bunny, dragon, dragon-fine if watertight (checked: none of them is, so they
!<    are reported and skipped), an icosphere of 6 subdivisions (81 920 facets, its vertices moved along the sphere
!<    towards a point so that it has no symmetry) at a pose, a marching-cubes torus at a pose; 2-D a NACA 0012 of 4 000
!<    segments turned by 10 degrees and an ellipse of 8 192 segments. FD: rule orders 1, 2, 3. FV: 64^3 cells (2-D 512^2)
!<    around the body, p at the walls' centroids, with the gradient and the walls' gyration for the exact moment. The
!<    error is printed in units of eps times the sum of the absolute terms: S_F = sum w |p| + V |g|, S_M = sum w |p|
!<    |x - c| + V |g| |x_c - c| (FV: + sum |K_w g| |A_w|).
!< 2. Potential flow (incompressible, inviscid, rho = U = 1, R = 0.5). The pressure is p = p_inf + rho (U^2 - |u|^2) / 2
!<    with u the exact potential-flow velocity ON THE TRUE WALL at the point nearest the sample (the radial projection
!<    c + R (x - c) / |x - c|), the convention of `fossil_val_fd_closure` (the solver imposes the true wall's data at the
!<    faceted wall): the error is then that of the faceted wall and of the rule alone. The alternative, u of the flow
!<    continued analytically to the sample point inside the true body, is printed at the finest level for comparison (it
!<    has the same order; that field does not exist in a solver, whose flow stops at its wall).
!<    - 2-D circular cylinder with circulation Gamma = -0.8 (2 pi U R): dw/dz = U (1 - R^2 / z^2) - i Gamma / (2 pi z),
!<      drag 0 (d'Alembert), lift -rho U Gamma (Kutta-Joukowski, Gamma counter-clockwise positive). Polygons of N = 32 ..
!<      4 096 segments, vertices on the circle, moved along it towards a point (the radial projection of a regular polygon
!<      from 0.3 R off the centre): no symmetry, so the drag is not zero by symmetry. Drag error |F_x| / (rho U^2 R) and
!<      lift error |F_y - L| / |L|, against the longest segment h_s, rule orders q = 1, 2, 3.
!<    - 3-D sphere in uniform flow along x: |u| = 3/2 U sin(theta) on the wall, zero force. The icospheres of 1 .. 6
!<      subdivisions, moved along the sphere as above. Error |F| / (rho U^2 pi R^2 / 2) against the longest edge h_s.
!<    - FV: the same polygons and icospheres refined with the cells (2-D N = 2 n on n^2 cells, n = 16 .. 1 024; 3-D
!<      2 .. 5 subdivisions on 8^3 .. 64^3 cells; longest element h_s about 1.4 h in 2-D, 1.15 h in 3-D), p at the
!<      walls' centroids: errors against h.
!< 3. FD vs FV: the same p on the same faceted wall by `integrate_loads` (order 3) and by the cut walls; their difference
!<    is the cut walls' midpoint rule on the wall fragments, of order 2.
!<
!< Pass criteria (`Are all tests passed? T` iff all hold; exit status 1 otherwise):
!<   1. buoyancy: every call STATUS_OK; the force within 1 (eps S_F) for every body, path and rule order; the moment
!<      within 1 (eps S_M) where it is exact: FD from order 2, FV with the gradient and the gyration (FD order 1 and the
!<      FV centroid-only moment are printed, not checked: their integrand is quadratic);
!<   2. potential flow: the d'Alembert and Kutta-Joukowski errors (3-D: the force) decrease at an observed order >= 1.8
!<      over the two finest levels, for every rule order (FD) and on the cut walls (FV): the faceted wall's second order;
!<   3. FD vs FV: their difference decreases at an observed order >= 1.8 over the two finest levels.
use fossil,                  only : surface_stl_object, curve_object, facet_object, load_samples_t, load_samples_3d, &
                                    load_samples_2d, integrate_loads, cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t, &
                                    cut_loads_3d, cut_walls_gyration, cut_cells_2d_t, cut_walls_2d_t, cut_walls_2d, &
                                    cut_loads_2d, extract_isosurface, rigid_pose_t, pose_from_axis_angle, STATUS_OK
use fossil_thickness_bodies, only : make_sphere
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P),    parameter :: PI = acos(-1._R8P), EPS = epsilon(1._R8P), ORDER_MIN = 1.8_R8P
real(R8P),    parameter :: G3(3) = [0.3_R8P, -0.7_R8P, 1.1_R8P], C3(3) = [0.1_R8P, -0.2_R8P, 0.3_R8P], P0 = 2._R8P
real(R8P),    parameter :: G2(2) = [0.4_R8P, -0.9_R8P], C2(2) = [0.15_R8P, -0.25_R8P]
real(R8P),    parameter :: R = 0.5_R8P, RHO = 1._R8P, UINF = 1._R8P, PINF = 0.7_R8P
real(R8P),    parameter :: GAM = -0.8_R8P * 2._R8P * PI * UINF * R, LIFT = -RHO * UINF * GAM
real(R8P),    parameter :: BETA = 0.3_R8P, ZC(2) = [0.13_R8P, -0.07_R8P], XC(3) = [0.11_R8P, -0.06_R8P, 0.08_R8P]
real(R8P),    parameter :: FREF2 = RHO * UINF**2 * R, FREF3 = 0.5_R8P * RHO * UINF**2 * PI * R**2
integer(I4P), parameter :: N2 = 8, N3 = 6, NV2 = 7, NV3 = 4
real(R8P)               :: a2(2), a3(3)
real(R8P)               :: h2(N2), ed2(N2, 3), el2(N2, 3), h3(N3), ef3(N3, 3)
real(R8P)               :: hv2(NV2), edv2(NV2), elv2(NV2), dd2(NV2), hv3(NV3), efv3(NV3), dd3(NV3)
real(R8P)               :: worst(7) = 0._R8P
integer(I4P)            :: nstat = 0_I4P
logical                 :: ok = .true.

a2 = [cos(0.7_R8P), sin(0.7_R8P)] ; a3 = [0.3_R8P, -0.5_R8P, 0.8_R8P] ; a3 = a3 / norm2(a3)
call buoyancy
call flow_2d
call flow_3d
call flow_fv_2d
call flow_fv_3d
call report
call check('every call STATUS_OK', nstat == 0_I4P)
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine buoyancy()
   !< Part 1: buoyancy through both paths on the large closed bodies.
   type(surface_stl_object)        :: s
   type(curve_object)              :: cv
   type(facet_object), allocatable :: fc(:)
   type(rigid_pose_t)              :: pose
   real(R8P), allocatable          :: f(:,:,:), u(:), v(:)
   real(R8P)                       :: ca, sa
   integer(I4P)                    :: i, j, k, st

   print '(A)', '1. Buoyancy in p = p0 + g . x: F = -g V, M = V (x_c - c) x (-g); error / (eps S), S = sum of |terms|'
   print '(A)', '   (q1 M and FV M c-only: relative error of the moment / S_M, not checked: quadratic integrand)'
   print '(A)', '   body                              elements  FD q1 F  q2 F  q2 M  q3 F  q3 M   FV F  FV M    q1 M rel  '// &
                'FV M c-only'
   call try_stl('src/tests/bunny.stl', 'bunny')
   call try_stl('src/tests/dragon.stl', 'dragon')
   call try_stl('src/tests/dragon-fine.stl', 'dragon-fine')
   pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.7_R8P, -0.3_R8P, 1.2_R8P])
   call stretched_sphere(6_I4P, XC, s) ; call s%move(pose)
   call buoyancy_3d(s, 'icosphere 6, stretched, posed')
   allocate(f(97, 97, 97))                                        ! a torus R 0.35, r 0.15 by marching cubes
   do k = 1, 97 ; do j = 1, 97 ; do i = 1, 97
      f(i, j, k) = torus_sdf(-0.6_R8P + 1.2_R8P * [i - 1, j - 1, k - 1] / 96._R8P)
   enddo ; enddo ; enddo
   call extract_isosurface(values=f, bmin=vector_R8P(-0.6_R8P, -0.6_R8P, -0.6_R8P), &
                           bmax=vector_R8P(0.6_R8P, 0.6_R8P, 0.6_R8P), iso=0._R8P, surface=fc, status=st)
   if (st /= 0) nstat = nstat + 1_I4P
   call s%adopt_facets(facets=fc)
   call s%move(pose_from_axis_angle([0.3_R8P, -1._R8P, 0.4_R8P], 0.9_R8P, [-0.2_R8P, 0.4_R8P, 0.1_R8P]))
   call buoyancy_3d(s, 'marching-cubes torus, posed')
   call naca0012(2000_I4P, u, v)
   ca = cos(10 * PI / 180) ; sa = sin(10 * PI / 180)
   call cv%adopt_loops(u=0.2_R8P + ca * u + sa * v, v=-0.1_R8P - sa * u + ca * v, loop_sizes=[size(u, kind=I4P)])
   call buoyancy_2d(cv, 'NACA 0012, 10 degrees')
   ca = cos(0.3_R8P) ; sa = sin(0.3_R8P)
   u = [(0.6_R8P * cos(2 * PI * i / 8192), i = 0, 8191)] ; v = [(0.25_R8P * sin(2 * PI * i / 8192), i = 0, 8191)]
   call cv%adopt_loops(u=0.31_R8P + ca * u - sa * v, v=0.17_R8P + sa * u + ca * v, loop_sizes=[8192_I4P])
   call buoyancy_2d(cv, 'ellipse 0.6 x 0.25, turned')
   endsubroutine buoyancy

   subroutine try_stl(file, name)
   !< A fixture, if watertight.
   character(*), intent(in) :: file, name !< File, label.
   type(surface_stl_object) :: s          !< Surface.

   call s%load_from_file(file_name=file, guess_format=.true.)
   if (s%is_watertight()) then
      call buoyancy_3d(s, name)
   else
      print '(3X,A,T38,I9,A)', name, s%get_facets_number(), '  not watertight: skipped (no volume to compare with)'
   endif
   endsubroutine try_stl

   subroutine buoyancy_3d(s, name)
   !< FD (orders 1-3) and FV (64^3 cells) buoyancy of one closed 3-D body.
   type(surface_stl_object), intent(in) :: s    !< Body.
   character(*),             intent(in) :: name !< Label.
   type(load_samples_t)   :: ls
   type(cut_cells_3d_t)   :: cut
   type(cut_walls_3d_t)   :: walls
   real(R8P), allocatable :: p(:), gy(:,:), x(:), y(:), z(:)
   real(R8P)              :: v, x0(3), fr(3), mr(3), f(3), m(3), sf, sm, ef(3), em(3), m1, efv, emv, m0, kg(3)
   integer(I4P)           :: q, i, st
   type(vector_R8P)       :: lo, hi

   call s%mass_properties(volume=v, center=x0, status=st)
   if (st /= STATUS_OK) then
      nstat = nstat + 1_I4P ; print '(3X,A,A)', name, ': mass_properties refused' ; return
   endif
   fr = -G3 * v ; mr = cross(v * (x0 - C3), -G3)
   do q = 1, 3
      call load_samples_3d(s, samples=ls, order=q, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(P0 + dot_product(G3, ls%points(:, i)), i = 1, ls%n)]
      call integrate_loads(ls, f, m, status=st, pressure=p, about=C3) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      sf = v * norm2(G3) ; sm = v * norm2(G3) * norm2(x0 - C3)
      do i = 1, ls%n
         sf = sf + ls%weight(i) * abs(p(i)) ; sm = sm + ls%weight(i) * abs(p(i)) * norm2(ls%points(:, i) - C3)
      enddo
      ef(q) = maxval(abs(f - fr)) / (EPS * sf) ; em(q) = maxval(abs(m - mr)) / (EPS * sm)
      if (q == 1) m1 = maxval(abs(m - mr)) / sm
   enddo
   lo = s%get_bmin() ; hi = s%get_bmax()
   x = lines(lo%x, hi%x, 64_I4P) ; y = lines(lo%y, hi%y, 64_I4P) ; z = lines(lo%z, hi%z, 64_I4P)
   call cut_cells_3d(s, x, y, z, cut, status=st, walls=walls) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
   p = [(P0 + dot_product(G3, walls%centroid(:, i)), i = 1, walls%n)]
   call cut_walls_gyration(walls, s, x, y, z, gy, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
   call cut_loads_3d(walls, f, m, status=st, pressure=p, gradient=spread(G3, 2, walls%n), gyration=gy, about=C3)
   if (st /= STATUS_OK) nstat = nstat + 1_I4P
   sf = v * norm2(G3) ; sm = v * norm2(G3) * norm2(x0 - C3)
   do i = 1, walls%n
      kg = [gy(1, i) * G3(1) + gy(4, i) * G3(2) + gy(6, i) * G3(3), gy(4, i) * G3(1) + gy(2, i) * G3(2) + gy(5, i) * G3(3), &
            gy(6, i) * G3(1) + gy(5, i) * G3(2) + gy(3, i) * G3(3)]
      sf = sf + abs(p(i)) * walls%wet_area(i)
      sm = sm + abs(p(i)) * walls%wet_area(i) * norm2(walls%centroid(:, i) - C3) + norm2(kg) * walls%wet_area(i)
   enddo
   efv = maxval(abs(f - fr)) / (EPS * sf) ; emv = maxval(abs(m - mr)) / (EPS * sm)
   call cut_loads_3d(walls, f, m, status=st, pressure=p, about=C3)
   m0 = maxval(abs(m - mr)) / sm
   print '(3X,A,T38,I9,2X,7F6.2,2ES11.2)', name, s%get_facets_number(), ef(1), ef(2), em(2), ef(3), em(3), efv, emv, m1, m0
   worst(1) = max(worst(1), maxval(ef), efv) ; worst(2) = max(worst(2), em(2), em(3), emv)
   call check(name//': buoyancy force within 1 eps S_F (FD q = 1, 2, 3, FV)', all(ef <= 1._R8P) .and. efv <= 1._R8P)
   call check(name//': buoyancy moment within 1 eps S_M (FD q = 2, 3, FV with gyration)', &
              em(2) <= 1._R8P .and. em(3) <= 1._R8P .and. emv <= 1._R8P)
   endsubroutine buoyancy_3d

   subroutine buoyancy_2d(cv, name)
   !< FD (orders 1-3) and FV (512^2 cells) buoyancy of one closed curve.
   type(curve_object), intent(in) :: cv   !< Curve.
   character(*),       intent(in) :: name !< Label.
   type(load_samples_t)   :: ls
   type(cut_cells_2d_t)   :: cut
   type(cut_walls_2d_t)   :: walls
   real(R8P), allocatable :: p(:), x(:), y(:)
   real(R8P)              :: aa, x0(2), fr(2), mr, f(2), m, sf, sm, ef(3), em(3), m1, efv, emv, m0, d(2), lo(2), hi(2)
   integer(I4P)           :: q, i, st

   call cv%mass_properties(area=aa, center=x0, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
   fr = -G2 * aa ; mr = aa * cross2(x0 - C2, -G2)
   do q = 1, 3
      call load_samples_2d(cv, samples=ls, order=q, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(P0 + dot_product(G2, ls%points(:, i)), i = 1, ls%n)]
      call integrate_loads(ls, f, m, status=st, pressure=p, about=C2) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      sf = aa * norm2(G2) ; sm = aa * norm2(G2) * norm2(x0 - C2)
      do i = 1, ls%n
         sf = sf + ls%weight(i) * abs(p(i)) ; sm = sm + ls%weight(i) * abs(p(i)) * norm2(ls%points(:, i) - C2)
      enddo
      ef(q) = maxval(abs(f - fr)) / (EPS * sf) ; em(q) = abs(m - mr) / (EPS * sm)
      if (q == 1) m1 = abs(m - mr) / sm
   enddo
   lo = cv%get_bmin() ; hi = cv%get_bmax()
   x = lines(lo(1), hi(1), 512_I4P) ; y = lines(lo(2), hi(2), 512_I4P)
   call cv%cut_cells(x, y, cut, status=st, boundaries=.true.) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
   call cut_walls_2d(cut, walls, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
   p = [(P0 + dot_product(G2, walls%centroid(:, i)), i = 1, walls%n)]
   call cut_loads_2d(walls, f, m, status=st, pressure=p, gradient=spread(G2, 2, walls%n), about=C2)
   if (st /= STATUS_OK) nstat = nstat + 1_I4P
   sf = aa * norm2(G2) ; sm = aa * norm2(G2) * norm2(x0 - C2)
   do i = 1, walls%n
      d = walls%ends(3:4, i) - walls%ends(1:2, i)
      sf = sf + abs(p(i)) * walls%length(i)
      sm = sm + abs(p(i)) * walls%length(i) * norm2(walls%centroid(:, i) - C2) + &
           norm2(d) * abs(dot_product(d, G2)) / 12._R8P * walls%length(i)
   enddo
   efv = maxval(abs(f - fr)) / (EPS * sf) ; emv = abs(m - mr) / (EPS * sm)
   call cut_loads_2d(walls, f, m, status=st, pressure=p, about=C2)
   m0 = abs(m - mr) / sm
   print '(3X,A,T38,I9,2X,7F6.2,2ES11.2)', name, cv%get_vertices_number(), ef(1), ef(2), em(2), ef(3), em(3), efv, emv, &
                                          m1, m0
   worst(1) = max(worst(1), maxval(ef), efv) ; worst(2) = max(worst(2), em(2), em(3), emv)
   call check(name//': buoyancy force within 1 eps S_F (FD q = 1, 2, 3, FV)', all(ef <= 1._R8P) .and. efv <= 1._R8P)
   call check(name//': buoyancy moment within 1 eps S_M (FD q = 2, 3, FV with the gradient)', &
              em(2) <= 1._R8P .and. em(3) <= 1._R8P .and. emv <= 1._R8P)
   endsubroutine buoyancy_2d

   subroutine flow_2d()
   !< Part 2, FD: the cylinder with circulation on stretched polygons, rule orders 1-3.
   type(curve_object)     :: cv
   type(load_samples_t)   :: ls
   real(R8P), allocatable :: p(:)
   real(R8P)              :: f(2), m
   integer(I4P)           :: k, q, i, n, st

   do k = 1, N2
      n = 32_I4P * 2_I4P**(k - 1)
      call stretched_circle(n, cv, h2(k))
      do q = 1, 3
         call load_samples_2d(cv, samples=ls, order=q, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
         p = [(pressure_2d(ls%points(:, i), .true.), i = 1, ls%n)]
         call integrate_loads(ls, f, m, status=st, pressure=p) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
         ed2(k, q) = abs(f(1)) / FREF2 ; el2(k, q) = abs(f(2) - LIFT) / abs(LIFT)
      enddo
   enddo
   p = [(pressure_2d(ls%points(:, i), .false.), i = 1, ls%n)]     ! the flow continued to the samples, finest, q = 3
   call integrate_loads(ls, f, m, status=st, pressure=p)
   print '(A)', ''
   print '(A)', '2. Potential flow, 2-D cylinder R = 0.5 with circulation: drag 0, lift -rho U Gamma (FD, load samples)'
   print '(A)', '   p from u on the true circle at the radial projection of the sample; errors |F_x| / (rho U^2 R), '// &
                '|F_y - L| / |L|'
   print '(A)', '      N    h_s       drag q=1   q=2        q=3        lift q=1   q=2        q=3'
   do k = 1, N2
      print '(I7,ES10.2,1X,6ES11.2)', 32 * 2**(k - 1), h2(k), ed2(k, :), el2(k, :)
   enddo
   print '(A17,1X,6F11.2)', 'order', (order(ed2(:, q), h2), q = 1, 3), (order(el2(:, q), h2), q = 1, 3)
   print '(A,2ES10.2)', '   u continued to the sample point instead (N = 4096, q = 3): drag, lift errors ', &
                        abs(f(1)) / FREF2, abs(f(2) - LIFT) / abs(LIFT)
   do q = 1, 3
      call check('2-D FD q = '//char(48 + q)//": d'Alembert error of order >= 1.8", order(ed2(:, q), h2) >= ORDER_MIN)
      call check('2-D FD q = '//char(48 + q)//': Kutta-Joukowski error of order >= 1.8', order(el2(:, q), h2) >= ORDER_MIN)
   enddo
   endsubroutine flow_2d

   subroutine flow_3d()
   !< Part 2, FD: the sphere in uniform flow on stretched icospheres, rule orders 1-3.
   type(surface_stl_object) :: s
   type(load_samples_t)     :: ls
   real(R8P), allocatable   :: p(:)
   real(R8P)                :: f(3), m(3), fa(3)
   integer(I4P)             :: k, q, i, st

   do k = 1, N3
      call stretched_sphere(k, XC, s, h3(k))
      do q = 1, 3
         call load_samples_3d(s, samples=ls, order=q, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
         p = [(pressure_3d(ls%points(:, i), .true.), i = 1, ls%n)]
         call integrate_loads(ls, f, m, status=st, pressure=p) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
         ef3(k, q) = norm2(f) / FREF3
      enddo
   enddo
   p = [(pressure_3d(ls%points(:, i), .false.), i = 1, ls%n)]
   call integrate_loads(ls, fa, m, status=st, pressure=p)
   print '(A)', ''
   print '(A)', '   3-D sphere R = 0.5 in uniform flow: F = 0 (FD, load samples); error |F| / (rho U^2 pi R^2 / 2)'
   print '(A)', '   level  facets    h_s       q=1        q=2        q=3'
   do k = 1, N3
      print '(I8,I8,ES10.2,1X,3ES11.2)', k, 20 * 4**k, h3(k), ef3(k, :)
   enddo
   print '(A26,1X,3F11.2)', 'order', (order(ef3(:, q), h3), q = 1, 3)
   print '(A,ES10.2)', '   u continued to the sample point instead (level 6, q = 3): ', norm2(fa) / FREF3
   do q = 1, 3
      call check('3-D FD q = '//char(48 + q)//": d'Alembert error of order >= 1.8", order(ef3(:, q), h3) >= ORDER_MIN)
   enddo
   endsubroutine flow_3d

   subroutine flow_fv_2d()
   !< Part 2 and 3, FV: the cylinder on n^2 cells, polygons of 2 n segments; FD (order 3) on the same polygon.
   type(curve_object)     :: cv
   type(cut_cells_2d_t)   :: cut
   type(cut_walls_2d_t)   :: walls
   type(load_samples_t)   :: ls
   real(R8P), allocatable :: p(:), x(:), y(:)
   real(R8P)              :: f(2), m, ffd(2), hs
   integer(I4P)           :: k, i, n, st

   do k = 1, NV2
      n = 16_I4P * 2_I4P**(k - 1) ; hv2(k) = 1.6_R8P / n
      call stretched_circle(2_I4P * n, cv, hs)
      x = [(ZC(1) - 0.8_R8P + 0.31_R8P * hv2(k) + i * hv2(k), i = 0, n)]
      y = [(ZC(2) - 0.8_R8P + 0.17_R8P * hv2(k) + i * hv2(k), i = 0, n)]
      call cv%cut_cells(x, y, cut, status=st, boundaries=.true.) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      call cut_walls_2d(cut, walls, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(pressure_2d(walls%centroid(:, i), .true.), i = 1, walls%n)]
      call cut_loads_2d(walls, f, m, status=st, pressure=p) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      edv2(k) = abs(f(1)) / FREF2 ; elv2(k) = abs(f(2) - LIFT) / abs(LIFT)
      call load_samples_2d(cv, samples=ls, order=3_I4P, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(pressure_2d(ls%points(:, i), .true.), i = 1, ls%n)]
      call integrate_loads(ls, ffd, m, status=st, pressure=p)
      dd2(k) = maxval(abs(ffd - f)) / FREF2
   enddo
   print '(A)', ''
   print '(A)', '   2-D cylinder on the cut walls (FV): n^2 cells of h = 1.6 / n, polygon of N = 2 n segments, p at the '// &
                'wall centroids'
   print '(A)', '      n      N    h         drag       lift       |F_FD(q=3) - F_FV| / (rho U^2 R)'
   do k = 1, NV2
      print '(2I7,ES10.2,3ES11.2)', 16 * 2**(k - 1), 32 * 2**(k - 1), hv2(k), edv2(k), elv2(k), dd2(k)
   enddo
   print '(A24,3F11.2)', 'order', order(edv2, hv2), order(elv2, hv2), order(dd2, hv2)
   call check("2-D FV: d'Alembert error of order >= 1.8", order(edv2, hv2) >= ORDER_MIN)
   call check('2-D FV: Kutta-Joukowski error of order >= 1.8', order(elv2, hv2) >= ORDER_MIN)
   call check('2-D FD vs FV: difference of order >= 1.8', order(dd2, hv2) >= ORDER_MIN)
   endsubroutine flow_fv_2d

   subroutine flow_fv_3d()
   !< Part 2 and 3, FV: the sphere on n^3 cells, icospheres of level 1 + log2(n / 4); FD (order 3) on the same mesh.
   type(surface_stl_object) :: s
   type(cut_cells_3d_t)     :: cut
   type(cut_walls_3d_t)     :: walls
   type(load_samples_t)     :: ls
   real(R8P), allocatable   :: p(:), x(:), y(:), z(:)
   real(R8P)                :: f(3), m(3), ffd(3), hs
   integer(I4P)             :: k, i, n, st

   print '(A)', ''
   print '(A)', '   3-D sphere on the cut walls (FV): n^3 cells of h = 1.6 / n, icosphere of level 1 + log2(n / 4), '// &
                'p at the wall centroids'
   print '(A)', '      n  level  facets   walls    h_s       h         |F| err    |F_FD(q=3) - F_FV| / '// &
                '(rho U^2 pi R^2 / 2)'
   do k = 1, NV3
      n = 8_I4P * 2_I4P**(k - 1) ; hv3(k) = 1.6_R8P / n
      call stretched_sphere(k + 1_I4P, XC, s, hs)
      x = [(XC(1) - 0.8_R8P + 0.31_R8P * hv3(k) + i * hv3(k), i = 0, n)]
      y = [(XC(2) - 0.8_R8P + 0.17_R8P * hv3(k) + i * hv3(k), i = 0, n)]
      z = [(XC(3) - 0.8_R8P + 0.23_R8P * hv3(k) + i * hv3(k), i = 0, n)]
      call cut_cells_3d(s, x, y, z, cut, status=st, walls=walls) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(pressure_3d(walls%centroid(:, i), .true.), i = 1, walls%n)]
      call cut_loads_3d(walls, f, m, status=st, pressure=p) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      efv3(k) = norm2(f) / FREF3
      call load_samples_3d(s, samples=ls, order=3_I4P, status=st) ; if (st /= STATUS_OK) nstat = nstat + 1_I4P
      p = [(pressure_3d(ls%points(:, i), .true.), i = 1, ls%n)]
      call integrate_loads(ls, ffd, m, status=st, pressure=p)
      dd3(k) = norm2(ffd - f) / FREF3
      print '(I7,I7,I8,I8,2ES10.2,2ES11.2)', n, k + 1, s%get_facets_number(), walls%n, hs, hv3(k), efv3(k), dd3(k)
   enddo
   print '(A49,2F11.2)', 'order', order(efv3, hv3), order(dd3, hv3)
   call check("3-D FV: d'Alembert error of order >= 1.8", order(efv3, hv3) >= ORDER_MIN)
   call check('3-D FD vs FV: difference of order >= 1.8', order(dd3, hv3) >= ORDER_MIN)
   endsubroutine flow_fv_3d

   subroutine report()
   !< Summary.
   print '(A)', ''
   print '(A,F6.2,A,F6.2,A)', 'buoyancy, largest error over the bodies and paths: force ', worst(1), ', exact moment ', &
                              worst(2), ' (eps S)'
   endsubroutine report

   real(R8P) function pressure_2d(x, project) result(pp)
   !< The cylinder's pressure p_inf + rho (U^2 - |u|^2) / 2, u at the radial projection of x on the circle (project) or
   !< at x itself (the flow continued inside the body).
   real(R8P), intent(in) :: x(2)    !< Point.
   logical,   intent(in) :: project !< Project on the circle.
   complex(R8P)          :: zz, dw  !< Point, complex velocity u - i v.
   real(R8P)             :: d(2)    !< Offset from the centre.

   d = x - ZC ; if (project) d = R * d / norm2(d)
   zz = cmplx(d(1), d(2), R8P)
   dw = UINF * (1._R8P - R**2 / zz**2) - cmplx(0._R8P, 1._R8P, R8P) * GAM / (2._R8P * PI * zz)
   pp = PINF + 0.5_R8P * RHO * (UINF**2 - abs(dw)**2)
   endfunction pressure_2d

   real(R8P) function pressure_3d(x, project) result(pp)
   !< The sphere's pressure in uniform flow along x: u = U (e + R^3 / 2 (e / r^3 - 3 (e . d) d / r^5)), d = x - c (or its
   !< radial projection on the sphere).
   real(R8P), intent(in) :: x(3)       !< Point.
   logical,   intent(in) :: project    !< Project on the sphere.
   real(R8P)             :: d(3), u(3), rr, e(3)

   e = [1._R8P, 0._R8P, 0._R8P]
   d = x - XC ; if (project) d = R * d / norm2(d)
   rr = norm2(d)
   u = UINF * (e + 0.5_R8P * R**3 * (e / rr**3 - 3._R8P * dot_product(e, d) * d / rr**5))
   pp = PINF + 0.5_R8P * RHO * (UINF**2 - dot_product(u, u))
   endfunction pressure_3d

   subroutine stretched_circle(n, cv, hs)
   !< The circle of radius R about ZC as n segments, vertices on it: a regular polygon projected radially from BETA R a2
   !< off the centre (counter-clockwise, no symmetry).
   integer(I4P),       intent(in)    :: n  !< Segments.
   type(curve_object), intent(inout) :: cv !< Curve.
   real(R8P),          intent(out)   :: hs !< Longest segment.
   real(R8P)    :: u(n), v(n), w(2), t
   integer(I4P) :: i

   do i = 1, n
      t = 0.1_R8P + 2._R8P * PI * (i - 1) / n
      w = [cos(t), sin(t)] + BETA * a2 ; w = w / norm2(w)
      u(i) = ZC(1) + R * w(1) ; v(i) = ZC(2) + R * w(2)
   enddo
   hs = 0._R8P
   do i = 1, n
      hs = max(hs, norm2([u(mod(i, n) + 1) - u(i), v(mod(i, n) + 1) - v(i)]))
   enddo
   call cv%adopt_loops(u=u, v=v, loop_sizes=[n])
   endsubroutine stretched_circle

   subroutine stretched_sphere(levels, cc, s, hs)
   !< The icosphere of radius R about cc, its unit vertices projected radially from BETA a3 (outward, no symmetry).
   integer(I4P),             intent(in)            :: levels !< Subdivisions.
   real(R8P),                intent(in)            :: cc(3)  !< Centre.
   type(surface_stl_object), intent(inout)         :: s      !< Surface.
   real(R8P),                intent(out), optional :: hs     !< Longest edge.
   type(facet_object), allocatable :: fc(:)
   real(R8P)                       :: w(3), vv(3, 3)
   integer(I4P)                    :: j, l

   call make_sphere(levels, fc)
   if (present(hs)) hs = 0._R8P
   do j = 1, size(fc)
      do l = 1, 3
         w = [fc(j)%vertex(l)%x, fc(j)%vertex(l)%y, fc(j)%vertex(l)%z] + BETA * a3 ; w = cc + R * w / norm2(w)
         fc(j)%vertex(l) = vector_R8P(w(1), w(2), w(3)) ; vv(:, l) = w
      enddo
      call fc(j)%compute_metrix
      if (present(hs)) hs = max(hs, norm2(vv(:, 2) - vv(:, 1)), norm2(vv(:, 3) - vv(:, 2)), norm2(vv(:, 1) - vv(:, 3)))
   enddo
   call s%adopt_facets(facets=fc)
   endsubroutine stretched_sphere

   real(R8P) function torus_sdf(x) result(d)
   !< Signed distance to the torus of radii 0.35, 0.15 about the z axis.
   real(R8P), intent(in) :: x(3) !< Point.

   d = sqrt((sqrt(x(1)**2 + x(2)**2) - 0.35_R8P)**2 + x(3)**2) - 0.15_R8P
   endfunction torus_sdf

   subroutine naca0012(mh, xx, yy)
   !< NACA 0012 of unit chord, closed trailing edge, 2 mh vertices, cosine spacing, counter-clockwise.
   integer(I4P),           intent(in)  :: mh           !< Points per side.
   real(R8P), allocatable, intent(out) :: xx(:), yy(:) !< Vertices.
   real(R8P)    :: xi
   integer(I4P) :: q, j

   allocate(xx(2 * mh), yy(2 * mh)) ; j = 0
   do q = mh, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = yt(xi)
   enddo
   do q = 1, mh - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * q / mh)) ; j = j + 1 ; xx(j) = xi ; yy(j) = -yt(xi)
   enddo
   endsubroutine naca0012

   pure real(R8P) function yt(xx)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: xx !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.1260_R8P * xx - 0.3516_R8P * xx ** 2 + 0.2843_R8P * xx ** 3 - &
                   0.1036_R8P * xx ** 4)
   endfunction yt

   function lines(lo, hi, n) result(l)
   !< n cells around [lo, hi] with uneven margins.
   real(R8P),    intent(in) :: lo, hi !< Extent.
   integer(I4P), intent(in) :: n      !< Cells.
   real(R8P)                :: l(n + 1)
   real(R8P)                :: a, b
   integer(I4P)             :: i

   a = lo - 0.0137_R8P - 0.031_R8P * (hi - lo) ; b = hi + 0.0119_R8P + 0.027_R8P * (hi - lo)
   l = [(a + (b - a) * i / n, i = 0, n)]
   endfunction lines

   real(R8P) function order(e, h) result(o)
   !< Observed order of the two finest levels.
   real(R8P), intent(in) :: e(:), h(:) !< Errors, sizes.
   integer(I4P)          :: n

   n = size(e, kind=I4P)
   o = log(e(n - 1) / e(n)) / log(h(n - 1) / h(n))
   endfunction order

   pure function cross(u, v) result(w)
   !< Cross product.
   real(R8P), intent(in) :: u(3), v(3) !< Vectors.
   real(R8P)             :: w(3)

   w = [u(2) * v(3) - u(3) * v(2), u(3) * v(1) - u(1) * v(3), u(1) * v(2) - u(2) * v(1)]
   endfunction cross

   pure real(R8P) function cross2(u, v)
   !< z component of u x v.
   real(R8P), intent(in) :: u(2), v(2) !< Vectors.

   cross2 = u(1) * v(2) - u(2) * v(1)
   endfunction cross2

   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (cond) then
      print '(A)', 'ok: '//label
   else
      print '(A)', 'FAILED: '//label ; ok = .false.
   endif
   endsubroutine check
endprogram fossil_val_surface_loads
