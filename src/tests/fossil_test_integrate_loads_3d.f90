!< FOSSIL test: the loads of a surface from its load samples (#109 WP2 step 02).
program fossil_test_integrate_loads_3d
!< FOSSIL test: the loads of a surface from its load samples (#109 WP2 step 02).
!<
!< 1. Buoyancy: in p = p0 + g . x (tau = 0) the force on a closed faceted body is F = -g V and the moment about c is
!<    M = V (x_c - c) x (-g), V and x_c from `mass_properties`: on cube-01, an icosphere (1280 facets) at a pose, two boxes
!<    (per body too), an irregular tetrahedron. The force is exact for every order (its integrand is linear on a flat facet), the moment from order 2
!<    (quadratic integrand); the error of each component is checked against TOL eps times the sum of the absolute terms
!<    (S_F = sum w |p| + V |g|, S_M = sum w |p| |x - c| + V |g| |x_c - c|); the largest ratio is printed. The order-1
!<    moment error is printed, not checked (not round-off on the tetrahedron; the symmetric bodies cancel it).
!< 2. Uniform shear: F = tau A, the area output = the sum of the facet areas.
!< 3. Per patch and per body: three patches over two boxes, the totals add up to the whole to round-off.
!< 4. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path bit for bit, p and tau
!<    together, resident samples too.
!< 5. Refusals (no p or tau, sizes, 2-D samples), the empty surface.
use fossil,                  only : surface_stl_object, load_samples_t, load_samples_3d, load_samples_reset, &
                                    integrate_loads, rigid_pose_t, pose_from_axis_angle, to_world, STATUS_OK, &
                                    STATUS_INVALID_INPUT, curve_object, load_samples_2d
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_boxes, make_sphere
use penf,                    only : I4P, I8P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: EPS = epsilon(1._R8P), TOL = 1._R8P
real(R8P), parameter            :: G(3) = [0.3_R8P, -0.7_R8P, 1.1_R8P], P0 = 2._R8P, C(3) = [0.1_R8P, -0.2_R8P, 0.3_R8P]
real(R8P), parameter            :: TAU(3) = [0.2_R8P, -0.5_R8P, 0.7_R8P]
type(surface_stl_object)        :: cube, sphere, boxes, tet, empty
type(facet_object), allocatable :: fc(:)
type(load_samples_t)            :: ls, lr, l2
type(rigid_pose_t)              :: pose
type(curve_object)              :: disc
real(R8P), allocatable          :: p(:), t(:,:), ppf(:,:), ppm(:,:), pbf(:,:), pbm(:,:), dpf(:,:), dpm(:,:), dbf(:,:), &
                                   dbm(:,:)
real(R8P)                       :: f(3), m(3), fd(3), md(3), a, ad, rf, rm, rm1, sf, sm, area
integer(I4P)                    :: o, s, st, k
logical                         :: ok

ok = .true. ; rf = 0._R8P ; rm = 0._R8P ; rm1 = 0._R8P
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call make_sphere(3_I4P, fc) ; call sphere%adopt_facets(facets=fc)
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, 0.25_R8P, fc) ; call boxes%adopt_facets(facets=fc)
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.7_R8P, -0.3_R8P, 1.2_R8P])
call make_tet([0._R8P, 0._R8P, 0._R8P], [1.3_R8P, 0.1_R8P, 0.2_R8P], [0.2_R8P, 0.9_R8P, -0.1_R8P], [0.4_R8P, 0.3_R8P, 1.1_R8P], fc)
call tet%adopt_facets(facets=fc)
! 1. buoyancy
do o = 1, 3
   call buoyancy(cube, 'cube-01', o, .false.)
   call buoyancy(sphere, 'icosphere posed', o, .true.)
   call buoyancy(boxes, 'two boxes', o, .false.)
   call buoyancy(tet, 'irregular tetrahedron', o, .false.)
enddo
print '(A,F6.2,A,F6.2,A,ES9.2)', 'largest error / (eps S): force ', rf, ', moment (orders 2, 3) ', rm, &
      ';  order-1 moment error / S_M ', rm1
! per body: each box -g V_b, V_b (x_b - c) x (-g)
call load_samples_3d(boxes, samples=ls, order=2_I4P)
p = pressure(ls)
call integrate_loads(ls, f, m, status=st, pressure=p, about=C, per_body_force=pbf, per_body_moment=pbm)
call scales(ls, p, 1._R8P, [0.5_R8P, 0.5_R8P, 0.5_R8P], sf, sm)
call check('per body: two bodies', st == STATUS_OK .and. size(pbf, 2) == 2 .and. size(pbm, 2) == 2)
call check('per body: box 1 buoyancy', all(abs(pbf(:, 1) + G) <= TOL * EPS * sf) .and. &
           all(abs(pbm(:, 1) - cross([0.5_R8P, 0.5_R8P, 0.5_R8P] - C, -G)) <= TOL * EPS * sm))
call check('per body: box 2 buoyancy', all(abs(pbf(:, 2) + G) <= TOL * EPS * sf) .and. &
           all(abs(pbm(:, 2) - cross([1.75_R8P, 0.5_R8P, 0.5_R8P] - C, -G)) <= TOL * EPS * sm))
! 2. uniform shear
call load_samples_3d(sphere, samples=ls, order=3_I4P)
allocate(t(3, ls%n))
do s = 1, ls%n
   t(:, s) = TAU
enddo
call integrate_loads(ls, f, m, status=st, shear=t, area=a)
area = surface_area(sphere)
call check('uniform shear: F = tau A (4 eps)', st == STATUS_OK .and. all(abs(f - TAU * a) <= 4._R8P * EPS * abs(TAU) * a))
call check('uniform shear: area = sum of the facet areas (1e-14)', abs(a - area) <= 1.e-14_R8P * area)
! 3. per patch and per body add up
call boxes%set_patch([1_I4P, 2_I4P, 3_I4P, 4_I4P, 5_I4P], 'two')
call boxes%set_patch([(k, k = 13_I4P, 20_I4P)], 'three')
call load_samples_3d(boxes, samples=ls, order=3_I4P)
p = pressure(ls) ; t = shear(ls)
call integrate_loads(ls, f, m, status=st, pressure=p, shear=t, about=C, per_patch_force=ppf, per_patch_moment=ppm, &
                     per_body_force=pbf, per_body_moment=pbm, area=a)
call scales(ls, p, 0._R8P, C, sf, sm)
call check('per patch: three patches, two bodies', st == STATUS_OK .and. size(ppf, 2) == 3 .and. size(pbf, 2) == 2)
call check('per patch: add up to the totals (4 eps S)', all(abs(sum(ppf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           all(abs(sum(ppm, 2) - m) <= 4._R8P * EPS * (sm + 4 * a)))
call check('per body: add up to the totals (4 eps S)', all(abs(sum(pbf, 2) - f) <= 4._R8P * EPS * (sf + 2 * a)) .and. &
           all(abs(sum(pbm, 2) - m) <= 4._R8P * EPS * (sm + 4 * a)))
call check('area: 12 unit faces', abs(a - 12._R8P) <= 8._R8P * EPS)
! 4. twins
call boxes%enter_device
call integrate_loads(ls, fd, md, status=st, pressure=p, shear=t, about=C, per_patch_force=dpf, per_patch_moment=dpm, &
                     per_body_force=dbf, per_body_moment=dbm, area=ad, device=.true.)
call check('twin: device = host, bit for bit', st == STATUS_OK .and. same(fd, f) .and. same(md, m) .and. same([ad], [a]) &
           .and. same2(dpf, ppf) .and. same2(dpm, ppm) .and. same2(dbf, pbf) .and. same2(dbm, pbm))
call load_samples_3d(boxes, samples=lr, order=3_I4P, resident=.true.)
call integrate_loads(lr, fd, md, status=st, pressure=p, shear=t, about=C, device=.true.)
call check('twin: resident samples on the device = host', st == STATUS_OK .and. lr%on_device .and. same(fd, f) .and. &
           same(md, m))
call load_samples_reset(lr)
call boxes%exit_device
call load_samples_3d(sphere, samples=ls, order=1_I4P, pose=pose)
p = pressure(ls) ; t = shear(ls)
call integrate_loads(ls, f, m, pressure=p, shear=t)
call integrate_loads(ls, fd, md, pressure=p, shear=t, device=.true.)
call check('twin: posed sphere, p and tau, default about', same(fd, f) .and. same(md, m))
call integrate_loads(ls, f, m, pressure=p) ; call integrate_loads(ls, fd, md, pressure=p, device=.true.)
call check('twin: pressure only', same(fd, f) .and. same(md, m))
call integrate_loads(ls, f, m, shear=t) ; call integrate_loads(ls, fd, md, shear=t, device=.true.)
call check('twin: shear only', same(fd, f) .and. same(md, m))
! 5. refusals, the empty surface
call integrate_loads(ls, f, m, status=st)
call check('refused: neither p nor tau', st == STATUS_INVALID_INPUT .and. all(f == 0._R8P) .and. all(m == 0._R8P))
call integrate_loads(ls, f, m, status=st, pressure=p(1:ls%n - 1))
call check('refused: p of the wrong size', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f, m, status=st, shear=t(1:2, :))
call check('refused: tau of the wrong shape', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f(1:2), m, status=st, pressure=p)
call check('refused: force of size 2', st == STATUS_INVALID_INPUT)
call integrate_loads(ls, f, m, status=st, pressure=p, about=C(1:2), per_patch_force=ppf)
call check('refused: about of size 2, per patch of size 0', st == STATUS_INVALID_INPUT .and. size(ppf, 2) == 0)
call disc%adopt_loops(u=[0._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P], loop_sizes=[3_I4P])
call load_samples_2d(disc, samples=l2)
call integrate_loads(l2, f, m, status=st, pressure=[1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P])
call check('refused: 2-D samples to the 3-D routine', st == STATUS_INVALID_INPUT)
call load_samples_3d(empty, samples=ls)
deallocate(p) ; allocate(p(0))
call integrate_loads(ls, f, m, status=st, pressure=p, per_body_force=pbf, area=a, device=.true.)
call check('empty surface: zero loads', st == STATUS_OK .and. all(f == 0._R8P) .and. all(m == 0._R8P) .and. a == 0._R8P &
           .and. size(pbf, 2) == 0)

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine buoyancy(surf, name, order, posed)
   !< The buoyancy of one body at one order: the force always, the moment from order 2.
   type(surface_stl_object), intent(in) :: surf  !< Body.
   character(*),             intent(in) :: name  !< Label.
   integer(I4P),             intent(in) :: order !< Rule order.
   logical,                  intent(in) :: posed !< At the pose.
   real(R8P) :: ff(3), mm(3), v, x0(3), fr(3), mr(3), ssf, ssm, e1, e2 !< Work.
   integer(I4P) :: sst                          !< Status.

   call surf%mass_properties(volume=v, center=x0, status=sst)
   if (posed) then
      call load_samples_3d(surf, samples=ls, order=order, pose=pose)
      x0 = to_world(pose, x0)
   else
      call load_samples_3d(surf, samples=ls, order=order)
   endif
   p = pressure(ls)
   call integrate_loads(ls, ff, mm, status=sst, pressure=p, about=C)
   fr = -G * v ; mr = cross(v * (x0 - C), -G)
   call scales(ls, p, v, x0, ssf, ssm)
   e1 = maxval(abs(ff - fr)) / (EPS * ssf) ; e2 = maxval(abs(mm - mr)) / (EPS * ssm)
   rf = max(rf, e1)
   if (order >= 2) then
      rm = max(rm, e2)
   else
      rm1 = max(rm1, e2 * EPS)
   endif
   call check('buoyancy '//name//' order '//char(48 + order)//': F = -g V', sst == STATUS_OK .and. e1 <= TOL)
   if (order >= 2) call check('buoyancy '//name//' order '//char(48 + order)//': M = V (xc - c) x (-g)', e2 <= TOL)
   endsubroutine buoyancy

   subroutine scales(smp, pp, v, x0, ssf, ssm)
   !< The sums of the absolute terms: S_F = sum w |p| + V |g|, S_M = sum w |p| |x - c| + V |g| |x0 - c|.
   type(load_samples_t), intent(in)  :: smp      !< Samples.
   real(R8P),            intent(in)  :: pp(:)    !< Pressure.
   real(R8P),            intent(in)  :: v, x0(3) !< Volume, centroid.
   real(R8P),            intent(out) :: ssf, ssm !< Scales.
   integer(I4P) :: q                             !< Counter.

   ssf = v * norm2(G) ; ssm = v * norm2(G) * norm2(x0 - C)
   do q = 1, smp%n
      ssf = ssf + smp%weight(q) * abs(pp(q)) ; ssm = ssm + smp%weight(q) * abs(pp(q)) * norm2(smp%points(:, q) - C)
   enddo
   endsubroutine scales

   function pressure(smp) result(pp)
   !< p = p0 + g . x at the samples.
   type(load_samples_t), intent(in) :: smp     !< Samples.
   real(R8P), allocatable           :: pp(:)   !< Pressure.
   integer(I4P) :: q                           !< Counter.

   allocate(pp(smp%n))
   do q = 1, smp%n
      pp(q) = P0 + ((G(1) * smp%points(1, q) + G(2) * smp%points(2, q)) + G(3) * smp%points(3, q))
   enddo
   endfunction pressure

   function shear(smp) result(tt)
   !< A varying traction at the samples (the twins and sums only).
   type(load_samples_t), intent(in) :: smp     !< Samples.
   real(R8P), allocatable           :: tt(:,:) !< Traction.
   integer(I4P) :: q                           !< Counter.

   allocate(tt(3, smp%n))
   do q = 1, smp%n
      tt(:, q) = [sin(3 * smp%points(2, q)), cos(2 * smp%points(3, q)), 0.3_R8P * smp%points(1, q)]
   enddo
   endfunction shear

   real(R8P) function surface_area(surf) result(aa)
   !< Sum of the facet areas.
   type(surface_stl_object), intent(in) :: surf  !< Surface.
   real(R8P)    :: vv(3, 3)                      !< Vertices.
   integer(I4P) :: ff                            !< Counter.

   aa = 0._R8P
   do ff = 1, surf%get_facets_number()
      vv = reshape(surf%tri_vertices(:, ff), [3, 3])
      aa = aa + 0.5_R8P * norm2(cross(vv(:, 2) - vv(:, 1), vv(:, 3) - vv(:, 1)))
   enddo
   endfunction surface_area

   subroutine make_tet(a, b, cc, d, ff)
   !< The tetrahedron a b cc d (d above the plane of a b cc, counter-clockwise from d), facets outward.
   real(R8P),                       intent(in)  :: a(3), b(3), cc(3), d(3) !< Vertices.
   type(facet_object), allocatable, intent(out) :: ff(:)                   !< Facets.

   allocate(ff(4))
   call put(ff(1), a, cc, b) ; call put(ff(2), a, b, d) ; call put(ff(3), b, cc, d) ; call put(ff(4), cc, a, d)
   endsubroutine make_tet

   subroutine put(fo, a, b, cc)
   !< One facet.
   type(facet_object), intent(inout) :: fo           !< Facet.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo%vertex(1) = vector_R8P(a(1), a(2), a(3)) ; fo%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo%vertex(3) = vector_R8P(cc(1), cc(2), cc(3)) ; call fo%compute_metrix
   endsubroutine put

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

   pure function cross(x, y) result(z)
   !< Cross product.
   real(R8P), intent(in) :: x(3), y(3) !< Vectors.
   real(R8P)             :: z(3)       !< x cross y.

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross

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
endprogram fossil_test_integrate_loads_3d
