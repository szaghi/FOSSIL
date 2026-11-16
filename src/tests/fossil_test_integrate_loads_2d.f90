!< FOSSIL test: the loads of a curve from its load samples (#109 WP2 step 02).
program fossil_test_integrate_loads_2d
!< FOSSIL test: the loads of a curve from its load samples (#109 WP2 step 02).
!<
!< 1. Buoyancy per unit depth: in p = p0 + g . x (tau = 0) the force on a closed polygon is F = -g A and the moment about
!<    c is M = A ((x_c - c) x (-g))_z, A and x_c from the curve's `mass_properties`: a 64-gon circle, a disc with a hole
!<    (clockwise), a NACA 0012 polygon (200 vertices) unposed and at a pose, two discs (per body too). The force is exact
!<    for every order, the moment from order 2; errors checked against TOL eps times the sums of the absolute terms
!<    (S_F = sum w |p| + A |g|, S_M = sum w |p| |x - c| + A |g| |x_c - c|), the largest ratio printed; the order-1 moment
!<    error printed, not checked.
!< 2. Uniform shear: F = tau L, the length output = the perimeter.
!< 3. Per patch and per body: three patches over two discs add up to the totals.
!< 4. Twin: `device=.true.` = host bit for bit (p and tau, p only, tau only, resident samples); refusals; the empty curve.
use fossil, only : curve_object, surface_stl_object, load_samples_t, load_samples_2d, load_samples_3d, load_samples_reset, &
                   integrate_loads, rigid_pose_2d_t, pose_2d, to_world, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object, only : facet_object
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter    :: PI = acos(-1._R8P), EPS = epsilon(1._R8P), TOL = 1._R8P
real(R8P), parameter    :: G(2) = [0.4_R8P, -0.9_R8P], P0 = 1.5_R8P, C(2) = [0.15_R8P, -0.25_R8P]
real(R8P), parameter    :: TAU(2) = [0.3_R8P, -0.6_R8P]
type(curve_object)      :: circle, holed, naca, discs, empty
type(surface_stl_object):: tri
type(facet_object), allocatable :: fc(:)
type(load_samples_t)    :: ls, lr, l3
type(rigid_pose_2d_t)   :: pose
real(R8P), allocatable  :: u(:), v(:), p(:), t(:,:), ppf(:,:), ppm(:), pbf(:,:), pbm(:), dpf(:,:), dpm(:), dbf(:,:), &
                           dbm(:)
real(R8P)               :: f(2), m, fd(2), md, a, ad, rf, rm, rm1, sf, sm, area
integer(I4P)            :: i, o, s, st
logical                 :: ok

ok = .true. ; rf = 0._R8P ; rm = 0._R8P ; rm1 = 0._R8P
u = [(0.5_R8P + 0.4_R8P * cos(2 * PI * i / 64), i = 0, 63)] ; v = [(0.5_R8P + 0.4_R8P * sin(2 * PI * i / 64), i = 0, 63)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
u = [u, [(0.5_R8P + 0.15_R8P * cos(-2 * PI * i / 32), i = 0, 31)]]
v = [v, [(0.6_R8P + 0.15_R8P * sin(-2 * PI * i / 32), i = 0, 31)]]
call holed%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 32_I4P])
call naca0012(100_I4P, u, v)
call naca%adopt_loops(u=u, v=v, loop_sizes=[200_I4P])
u = [[(0.3_R8P + 0.25_R8P * cos(2 * PI * i / 48), i = 0, 47)], [(1.1_R8P + 0.2_R8P * cos(2 * PI * i / 40), i = 0, 39)]]
v = [[(0.4_R8P + 0.25_R8P * sin(2 * PI * i / 48), i = 0, 47)], [(0.2_R8P + 0.2_R8P * sin(2 * PI * i / 40), i = 0, 39)]]
call discs%adopt_loops(u=u, v=v, loop_sizes=[48_I4P, 40_I4P], &
                       patches=[(1_I4P, i = 1, 20), (2_I4P, i = 21, 48), (3_I4P, i = 1, 40)], &
                       patch_names=['one  ', 'two  ', 'three'])
pose = pose_2d(0.35_R8P, [0.4_R8P, -0.7_R8P])
! 1. buoyancy
do o = 1, 3
   call buoyancy(circle, 'circle 64-gon', o, .false.)
   call buoyancy(holed, 'disc with a hole', o, .false.)
   call buoyancy(naca, 'NACA 0012', o, .false.)
   call buoyancy(naca, 'NACA 0012 posed', o, .true.)
   call buoyancy(discs, 'two discs', o, .false.)
enddo
print '(A,F6.2,A,F6.2,A,ES9.2)', 'largest error / (eps S): force ', rf, ', moment (orders 2, 3) ', rm, &
      ';  order-1 moment error / S_M ', rm1
! per body: each disc its own buoyancy
call load_samples_2d(discs, samples=ls, order=2_I4P)
p = pressure(ls)
call integrate_loads(ls, f, m, status=st, pressure=p, about=C, per_body_force=pbf, per_body_moment=pbm)
call check('per body: two bodies', st == STATUS_OK .and. size(pbf, 2) == 2 .and. size(pbm) == 2)
call body_check(1_I4P, 0_I4P, 48_I4P) ; call body_check(2_I4P, 48_I4P, 40_I4P)
! 2. uniform shear
call load_samples_2d(naca, samples=ls, order=3_I4P)
allocate(t(2, ls%n))
do s = 1, ls%n
   t(:, s) = TAU
enddo
call integrate_loads(ls, f, m, status=st, shear=t, area=a)
area = naca%get_perimeter()
call check('uniform shear: F = tau L (4 eps)', st == STATUS_OK .and. all(abs(f - TAU * a) <= 4._R8P * EPS * abs(TAU) * a))
call check('uniform shear: length = perimeter (1e-14)', abs(a - area) <= 1.e-14_R8P * area)
! 3. per patch and per body add up
call load_samples_2d(discs, samples=ls, order=3_I4P)
p = pressure(ls) ; t = shear(ls)
call integrate_loads(ls, f, m, status=st, pressure=p, shear=t, about=C, per_patch_force=ppf, per_patch_moment=ppm, &
                     per_body_force=pbf, per_body_moment=pbm, area=a)
call scales(ls, p, 0._R8P, C, sf, sm)
call check('per patch: three patches, two bodies', st == STATUS_OK .and. size(ppf, 2) == 3 .and. size(ppm) == 3 .and. &
           size(pbf, 2) == 2)
call check('per patch: add up to the totals (4 eps S)', all(abs(sum(ppf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           abs(sum(ppm) - m) <= 4._R8P * EPS * (sm + 4 * a))
call check('per body: add up to the totals (4 eps S)', all(abs(sum(pbf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           abs(sum(pbm) - m) <= 4._R8P * EPS * (sm + 4 * a))
call check('length: the two perimeters', abs(a - discs%get_perimeter()) <= 1.e-14_R8P * a)
! 4. twins
call discs%enter_device
call integrate_loads(ls, fd, md, status=st, pressure=p, shear=t, about=C, per_patch_force=dpf, per_patch_moment=dpm, &
                     per_body_force=dbf, per_body_moment=dbm, area=ad, device=.true.)
call check('twin: device = host, bit for bit', st == STATUS_OK .and. same(fd, f) .and. same([md], [m]) .and. &
           same([ad], [a]) .and. same2(dpf, ppf) .and. same(dpm, ppm) .and. same2(dbf, pbf) .and. same(dbm, pbm))
call load_samples_2d(discs, samples=lr, order=3_I4P, resident=.true.)
call integrate_loads(lr, fd, md, status=st, pressure=p, shear=t, about=C, device=.true.)
call check('twin: resident samples on the device = host', st == STATUS_OK .and. lr%on_device .and. same(fd, f) .and. &
           same([md], [m]))
call load_samples_reset(lr)
call discs%exit_device
call load_samples_2d(naca, samples=ls, order=2_I4P, pose=pose)
p = pressure(ls) ; t = shear(ls)
call integrate_loads(ls, f, m, pressure=p, shear=t) ; call integrate_loads(ls, fd, md, pressure=p, shear=t, device=.true.)
call check('twin: posed NACA, p and tau, default about', same(fd, f) .and. same([md], [m]))
call integrate_loads(ls, f, m, pressure=p) ; call integrate_loads(ls, fd, md, pressure=p, device=.true.)
call check('twin: pressure only', same(fd, f) .and. same([md], [m]))
call integrate_loads(ls, f, m, shear=t) ; call integrate_loads(ls, fd, md, shear=t, device=.true.)
call check('twin: shear only', same(fd, f) .and. same([md], [m]))
! refusals, the empty curve
call integrate_loads(ls, f, m, status=st)
call check('refused: neither p nor tau', st == STATUS_INVALID_INPUT .and. all(f == 0._R8P) .and. m == 0._R8P)
call integrate_loads(ls, f, m, status=st, pressure=p(1:ls%n - 1))
call check('refused: p of the wrong size', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f, m, status=st, shear=t(1:1, :))
call check('refused: tau of the wrong shape', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f(1:1), m, status=st, pressure=p)
call check('refused: force of size 1', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f, m, status=st, pressure=p, about=[0._R8P, 0._R8P, 0._R8P], per_body_moment=pbm)
call check('refused: about of size 3, per body of size 0', st == STATUS_INVALID_INPUT .and. size(pbm) == 0)
allocate(fc(1)) ; fc(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; fc(1)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
fc(1)%vertex(3) = vector_R8P(0._R8P, 1._R8P, 0._R8P) ; call fc(1)%compute_metrix
call tri%adopt_facets(facets=fc)
call load_samples_3d(tri, samples=l3, order=1_I4P)
call integrate_loads(l3, f, m, status=st, pressure=[1._R8P])
call check('refused: 3-D samples to the 2-D routine', st == STATUS_INVALID_INPUT)
call load_samples_2d(empty, samples=ls)
deallocate(p) ; allocate(p(0))
call integrate_loads(ls, f, m, status=st, pressure=p, per_patch_force=pbf, area=a, device=.true.)
call check('empty curve: zero loads', st == STATUS_OK .and. all(f == 0._R8P) .and. m == 0._R8P .and. a == 0._R8P .and. &
           size(pbf, 2) == 0)

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine buoyancy(cur, name, order, posed)
   !< The buoyancy of one curve at one order: the force always, the moment from order 2.
   type(curve_object), intent(in) :: cur   !< Curve.
   character(*),       intent(in) :: name  !< Label.
   integer(I4P),       intent(in) :: order !< Rule order.
   logical,            intent(in) :: posed !< At the pose.
   real(R8P)    :: ff(2), mm, aa, x0(2), fr(2), mr, ssf, ssm, e1, e2 !< Work.
   integer(I4P) :: sst                                               !< Status.

   call cur%mass_properties(area=aa, center=x0, status=sst)
   if (posed) then
      call load_samples_2d(cur, samples=ls, order=order, pose=pose)
      x0 = to_world(pose, x0)
   else
      call load_samples_2d(cur, samples=ls, order=order)
   endif
   p = pressure(ls)
   call integrate_loads(ls, ff, mm, status=sst, pressure=p, about=C)
   fr = -G * aa ; mr = aa * cross2(x0 - C, -G)
   call scales(ls, p, aa, x0, ssf, ssm)
   e1 = maxval(abs(ff - fr)) / (EPS * ssf) ; e2 = abs(mm - mr) / (EPS * ssm)
   rf = max(rf, e1)
   if (order >= 2) then
      rm = max(rm, e2)
   else
      rm1 = max(rm1, e2 * EPS)
   endif
   call check('buoyancy '//name//' order '//char(48 + order)//': F = -g A', sst == STATUS_OK .and. e1 <= TOL)
   if (order >= 2) call check('buoyancy '//name//' order '//char(48 + order)//': M = A (xc - c) x (-g)', e2 <= TOL)
   endsubroutine buoyancy

   subroutine body_check(b, i0, nv)
   !< Body b (vertices i0 + 1 ... i0 + nv): its own buoyancy, from a curve of that loop alone.
   integer(I4P), intent(in) :: b, i0, nv !< Body, first vertex - 1, vertices.
   type(curve_object) :: one              !< The body alone.
   real(R8P)          :: aa, x0(2), ssf, ssm !< Area, centroid, scales.

   call one%adopt_loops(u=u(i0 + 1:i0 + nv), v=v(i0 + 1:i0 + nv), loop_sizes=[nv])
   call one%mass_properties(area=aa, center=x0)
   call scales(ls, p, aa, x0, ssf, ssm)
   call check('per body: disc '//char(48 + b)//' buoyancy', all(abs(pbf(:, b) + G * aa) <= TOL * EPS * ssf) .and. &
              abs(pbm(b) - aa * cross2(x0 - C, -G)) <= TOL * EPS * ssm)
   endsubroutine body_check

   subroutine scales(smp, pp, aa, x0, ssf, ssm)
   !< The sums of the absolute terms: S_F = sum w |p| + A |g|, S_M = sum w |p| |x - c| + A |g| |x0 - c|.
   type(load_samples_t), intent(in)  :: smp       !< Samples.
   real(R8P),            intent(in)  :: pp(:)     !< Pressure.
   real(R8P),            intent(in)  :: aa, x0(2) !< Area, centroid.
   real(R8P),            intent(out) :: ssf, ssm  !< Scales.
   integer(I4P) :: q                              !< Counter.

   ssf = aa * norm2(G) ; ssm = aa * norm2(G) * norm2(x0 - C)
   do q = 1, smp%n
      ssf = ssf + smp%weight(q) * abs(pp(q)) ; ssm = ssm + smp%weight(q) * abs(pp(q)) * norm2(smp%points(:, q) - C)
   enddo
   endsubroutine scales

   function pressure(smp) result(pp)
   !< p = p0 + g . x at the samples.
   type(load_samples_t), intent(in) :: smp   !< Samples.
   real(R8P), allocatable           :: pp(:) !< Pressure.
   integer(I4P) :: q                         !< Counter.

   allocate(pp(smp%n))
   do q = 1, smp%n
      pp(q) = P0 + (G(1) * smp%points(1, q) + G(2) * smp%points(2, q))
   enddo
   endfunction pressure

   function shear(smp) result(tt)
   !< A varying traction at the samples (the twins and sums only).
   type(load_samples_t), intent(in) :: smp     !< Samples.
   real(R8P), allocatable           :: tt(:,:) !< Traction.
   integer(I4P) :: q                           !< Counter.

   allocate(tt(2, smp%n))
   do q = 1, smp%n
      tt(:, q) = [sin(3 * smp%points(2, q)), 0.4_R8P * cos(2 * smp%points(1, q))]
   enddo
   endfunction shear

   subroutine naca0012(mh, x, y)
   !< NACA 0012 of unit chord, closed trailing edge, 2 mh vertices, cosine spacing, counter-clockwise.
   integer(I4P),           intent(in)  :: mh   !< Points per side.
   real(R8P), allocatable, intent(out) :: x(:), y(:) !< Vertices.
   real(R8P)    :: xi                          !< Abscissa.
   integer(I4P) :: k, j                        !< Counters.

   allocate(x(2 * mh), y(2 * mh)) ; j = 0
   do k = mh, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * k / mh)) ; j = j + 1 ; x(j) = xi ; y(j) = yt(xi)
   enddo
   do k = 1, mh - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * k / mh)) ; j = j + 1 ; x(j) = xi ; y(j) = -yt(xi)
   enddo
   endsubroutine naca0012

   pure real(R8P) function yt(x)
   !< NACA 0012 half thickness, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.

   yt = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x ** 2 + 0.2843_R8P * x ** 3 - 0.1036_R8P * x ** 4)
   endfunction yt

   pure real(R8P) function cross2(x, y)
   !< z component of x cross y.
   real(R8P), intent(in) :: x(2), y(2) !< Vectors.

   cross2 = x(1) * y(2) - x(2) * y(1)
   endfunction cross2

   logical function same(x, y)
   !< Equal bits.
   real(R8P), intent(in) :: x(:), y(:) !< Arrays.

   same = size(x) == size(y)
   if (same) same = all(transfer(x, 0_I8P, size(x)) == transfer(y, 0_I8P, size(y)))
   endfunction same

   logical function same2(x, y)
   !< Equal bits and shapes.
   real(R8P), intent(in) :: x(:,:), y(:,:) !< Arrays.

   same2 = all(shape(x) == shape(y))
   if (same2) same2 = same(reshape(x, [size(x)]), reshape(y, [size(y)]))
   endfunction same2

   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (.not. cond) then
      print '(A)', 'FAILED: '//label
      ok = .false.
   else
      print '(A)', 'ok: '//label
   endif
   endsubroutine check
endprogram fossil_test_integrate_loads_2d
