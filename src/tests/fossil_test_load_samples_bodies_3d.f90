!< FOSSIL test: the load samples of several bodies (#109 WP3 step 02).
program fossil_test_load_samples_bodies_3d
!< FOSSIL test: the load samples of several bodies (#109 WP3 step 02).
!<
!< 1. One body: `load_samples_3d_bodies` equals `load_samples_3d` bit for bit in every single-body field (probes, lines,
!<    order 3; posed; the device path), body = 1, sub_body = the facet_body, no second body (clearance huge); an empty
!<    second surface is ignored.
!< 2. Two unit spheres far apart (posed): the samples of body b are the slots of its single-body set, bit for bit, in
!<    body order (body = b, sub_body its facet_body, probes shorter than the gap all as the single-body set).
!< 3. Clearance: every sample's clearance is |x - c_other| - 1 to the facet sag, the second body the other sphere, its
!<    second foot on the other sphere; the least clearance is the gap 0.6 to the sag.
!< 4. Visibility: two boxes 1/8 apart: on the facing faces a probe shorter than the gap is visible, a longer one enters
!<    the other box and is not (the single-body set sees both), clearance 1/8 and the second foot across to round-off.
!< 5. Buoyancy: a sphere and a box, posed: `integrate_loads` with p = p0 + g . x gives per body F_b = -g V_b and (order
!<    2) M_b = V_b (x_c,b - about) x (-g) to round-off; per_body means the surfaces; the totals are their sums.
!< 6. Twin: posed bodies 0.15 apart with probes and lines: `device=.true.` (built without OpenACC its loops run on the
!<    host) equals the host path in every field, bit for bit; resident, then `update_host`.
!< 7. Refusals: an octree surface, size(poses) /= size(surfaces), no surface, the device path of unmapped surfaces.
use fossil,                  only : surface_stl_object, load_samples_t, load_samples_3d, load_samples_3d_bodies, &
                                    load_samples_reset, load_samples_mismatches, integrate_loads, rigid_pose_t, &
                                    pose_from_axis_angle, to_world, AABB_TREE_OCTREE, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_boxes, make_sphere
use iso_fortran_env,         only : int64
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: G8 = 0.125_R8P             !< Box gap.
type(surface_stl_object)        :: sphere, boxa, boxb, octs, empty
type(surface_stl_object)        :: one(1), two(2), boxes(2), mix(2)
type(facet_object), allocatable :: fc(:)
type(load_samples_t)            :: sh, sb, sd, sa, sc, sr
type(rigid_pose_t)              :: pa, pb
type(vector_R8P)                :: cw
real(R8P), allocatable          :: x(:), p(:), pf(:,:), pm(:,:)
real(R8P)                       :: e, emin, gv(3), c(3), vol, fv(3), mv(3), ftot(3), mtot(3), sc1, xc(3)
integer(I4P)                    :: s, st, nbad, nf, b, n1
logical                         :: ok

ok = .true.
call make_sphere(2_I4P, fc) ; call sphere%adopt_facets(facets=fc)
x = [(-2.45_R8P + 0.2_R8P * s, s = 0, 24)]
! 1. one body
one(1) = sphere
call load_samples_3d(sphere, x, x, x, sh, status=st, order=3_I4P, lengths=[0.1_R8P, 0.4_R8P])
call load_samples_3d_bodies(one, x, x, x, sb, status=st, order=3_I4P, lengths=[0.1_R8P, 0.4_R8P])
call check('one body: status, flags', st == STATUS_OK .and. sb%n == sh%n .and. sh%n == 1920_I4P .and. sb%with_bodies .and. &
           sb%bodies == 1_I4P)
call check('one body: = load_samples_3d bit for bit in every single-body field', load_samples_mismatches(sh, strip(sb)) == 0)
call check('one body: body 1, sub_body the facet_body, no second body', all(sb%body == 1_I4P) .and. &
           all(sb%sub_body == sh%body) .and. all(sb%second_body == 0_I4P) .and. all(sb%clearance == huge(1._R8P)))
pa = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.3_R8P, [0.05_R8P, -0.02_R8P, 0.03_R8P])
call load_samples_3d(sphere, x, x, x, sh, status=st, lengths=[-0.1_R8P, 0.3_R8P], pose=pa)
call load_samples_3d_bodies(one, x, x, x, sb, status=st, poses=[pa], lengths=[-0.1_R8P, 0.3_R8P])
call check('one body posed: = load_samples_3d bit for bit', st == STATUS_OK .and. sb%n > 0_I4P .and. &
           load_samples_mismatches(sh, strip(sb)) == 0_I4P)
call sphere%enter_device ; call one(1)%enter_device ! map the array elements themselves (a copy is not mapped)
call load_samples_3d(sphere, x, x, x, sd, status=st, lengths=[-0.1_R8P, 0.3_R8P], pose=pa, device=.true.)
call load_samples_3d_bodies(one, x, x, x, sc, status=st, poses=[pa], lengths=[-0.1_R8P, 0.3_R8P], device=.true.)
call check('one body posed, device path: = load_samples_3d(device) and the host bodies set, bit for bit', &
           st == STATUS_OK .and. load_samples_mismatches(sd, strip(sc)) == 0_I4P .and. &
           load_samples_mismatches(sb, sc) == 0_I4P)
call sphere%exit_device ; call one(1)%exit_device
two(1) = sphere ; two(2) = empty
call load_samples_3d_bodies(two, x, x, x, sc, status=st, poses=[pa, pa], lengths=[-0.1_R8P, 0.3_R8P])
call check('an empty second surface is ignored', st == STATUS_OK .and. sc%bodies == 2_I4P .and. &
           load_samples_mismatches(sh, strip(sc)) == 0_I4P .and. all(sc%second_body == 0_I4P))
! 2. two spheres far apart (gap 0.6), probes shorter than the gap
pa = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.2_R8P, [-1.3_R8P, 0.03_R8P, 0.01_R8P])
pb = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], -0.3_R8P, [1.3_R8P, -0.02_R8P, 0.04_R8P])
two(1) = sphere ; two(2) = sphere
call load_samples_3d_bodies(two, x, x, x, sb, status=st, poses=[pa, pb], lengths=[0.1_R8P, 0.25_R8P])
call load_samples_3d(sphere, x, x, x, sa, lengths=[0.1_R8P, 0.25_R8P], pose=pa)
call load_samples_3d(sphere, x, x, x, sc, lengths=[0.1_R8P, 0.25_R8P], pose=pb)
call check('far apart: the two sets in body order', st == STATUS_OK .and. sb%n == sa%n + sc%n .and. &
           sb%ne == sa%ne + sc%ne .and. all(sb%body(1:sa%n) == 1_I4P) .and. all(sb%body(sa%n + 1:) == 2_I4P))
call check('far apart: every body''s slots = its single-body set, bit for bit', same_slice(sb, 0_I4P, sa) .and. &
           same_slice(sb, sa%n, sc))
! 3. clearance
e = 0._R8P ; nbad = 0_I4P ; emin = huge(1._R8P)
do s = 1_I4P, sb%n
   c = pb%t ; if (sb%body(s) == 2_I4P) c = pa%t
   e = max(e, abs(sb%clearance(s) - (norm2(sb%points(:, s) - c) - 1._R8P)))
   emin = min(emin, sb%clearance(s))
   if (sb%second_body(s) /= 3_I4P - sb%body(s)) nbad = nbad + 1_I4P
   if (abs(norm2(sb%second_foot(:, s) - c) - 1._R8P) > 0.02_R8P) nbad = nbad + 1_I4P
   if (abs(norm2(sb%second_foot(:, s) - sb%points(:, s)) - sb%clearance(s)) > 1.e-12_R8P) nbad = nbad + 1_I4P
enddo
call check('clearance: |x - c_other| - 1 (facet sag), second body and foot', e < 0.03_R8P .and. nbad == 0_I4P)
call check('clearance: the least one is the gap 0.6 (facet sag)', abs(emin - (norm2(pb%t - pa%t) - 2._R8P)) < 0.03_R8P)
! 4. visibility across a narrow gap: two boxes 1/8 apart
call make_boxes([-1._R8P, -0.6_R8P, -0.6_R8P], [-0.5_R8P * G8, 0.6_R8P, 0.6_R8P], 1_I4P, 0._R8P, fc)
call boxa%adopt_facets(facets=fc)
call make_boxes([0.5_R8P * G8, -0.6_R8P, -0.6_R8P], [1._R8P, 0.6_R8P, 0.6_R8P], 1_I4P, 0._R8P, fc)
call boxb%adopt_facets(facets=fc)
boxes(1) = boxa ; boxes(2) = boxb
call load_samples_3d_bodies(boxes, samples=sb, status=st, lengths=[0.1_R8P, 0.2_R8P])
call load_samples_3d(boxa, samples=sa, lengths=[0.1_R8P, 0.2_R8P])
nf = 0_I4P ; nbad = 0_I4P
do s = 1_I4P, sb%n
   if (sb%weight(s) <= 0._R8P) cycle
   if (abs(sb%points(1, s) + sb%normal(1, s) * 0.5_R8P * G8) > 1.e-14_R8P .or. abs(sb%normal(1, s)) < 0.5_R8P) then
      if (.not. all(sb%visible(:, s))) nbad = nbad + 1_I4P ! no other face is shadowed
      cycle
   endif
   nf = nf + 1_I4P ! a sample on a facing face
   if (.not. sb%visible(1, s) .or. sb%visible(2, s)) nbad = nbad + 1_I4P
   if (abs(sb%clearance(s) - G8) > 1.e-14_R8P .or. sb%second_body(s) /= 3_I4P - sb%body(s)) nbad = nbad + 1_I4P
   if (abs(sb%second_foot(1, s) + sb%points(1, s)) > 1.e-14_R8P .or. &
       any(abs(sb%second_foot(2:3, s) - sb%points(2:3, s)) > 1.e-14_R8P)) nbad = nbad + 1_I4P
   if (sb%body(s) == 1_I4P) then
      if (.not. all(sa%visible(:, s))) nbad = nbad + 1_I4P ! alone, box a does not block it
   endif
enddo
call check('narrow gap: facing faces, the long probe enters the other box, clearance 1/8, second foot across', &
           st == STATUS_OK .and. nf == 12_I4P .and. nbad == 0_I4P)
! 5. buoyancy per body: a sphere and a box, posed
pa = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], 0.4_R8P, [-1.2_R8P, 0.1_R8P, 0.2_R8P])
pb = pose_from_axis_angle([1._R8P, 0._R8P, -1._R8P], -0.2_R8P, [1.1_R8P, -0.3_R8P, 0.1_R8P])
mix(1) = sphere ; mix(2) = boxa
gv = [0.3_R8P, -0.2_R8P, 9.81_R8P] ; c = [0.1_R8P, 0.2_R8P, -0.3_R8P]
call load_samples_3d_bodies(mix, samples=sb, status=st, poses=[pa, pb])
allocate(p(sb%n))
do s = 1_I4P, sb%n
   p(s) = 2._R8P + dot_product(gv, sb%points(:, s))
enddo
call integrate_loads(sb, ftot, mtot, status=st, pressure=p, about=c, per_body_force=pf, per_body_moment=pm)
nbad = 0_I4P ; e = 0._R8P
do b = 1_I4P, 2_I4P
   vol = mix(b)%get_volume() ; cw = mix(b)%get_centroid() ; xc = [cw%x, cw%y, cw%z]
   if (b == 1_I4P) then
      xc = to_world(pa, xc)
   else
      xc = to_world(pb, xc)
   endif
   fv = -gv * vol
   mv = cross(xc - c, -gv) * vol
   sc1 = norm2(gv) * vol * 3._R8P
   e = max(e, maxval(abs(pf(:, b) - fv)) / sc1, maxval(abs(pm(:, b) - mv)) / sc1)
enddo
call check('buoyancy: per body F_b = -g V_b, M_b = V_b (x_c,b - c) x (-g) to round-off', st == STATUS_OK .and. &
           size(pf, 2) == 2 .and. e < 1.e-13_R8P)
call check('buoyancy: the totals are the sums of the bodies', &
           maxval(abs(ftot - (pf(:, 1) + pf(:, 2)))) < 1.e-13_R8P * norm2(ftot) .and. &
           maxval(abs(mtot - (pm(:, 1) + pm(:, 2)))) < 1.e-12_R8P * (norm2(mtot) + norm2(ftot)))
! 6. twin: posed spheres 0.15 apart with probes and lines
pa = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], 0.4_R8P, [-1.075_R8P, 0.02_R8P, -0.01_R8P])
pb = pose_from_axis_angle([1._R8P, 0._R8P, -1._R8P], -0.2_R8P, [1.075_R8P, -0.03_R8P, 0.02_R8P])
two(1) = sphere ; two(2) = sphere ; call two(1)%enter_device ; call two(2)%enter_device
call load_samples_3d_bodies(two, x, x, x, sh, status=st, poses=[pa, pb], lengths=[0.1_R8P, 0.25_R8P])
call load_samples_3d_bodies(two, x, x, x, sd, status=st, poses=[pa, pb], lengths=[0.1_R8P, 0.25_R8P], device=.true.)
n1 = count(.not. sh%visible)
call check('twin: device = host, every field', st == STATUS_OK .and. sh%n > 0_I4P .and. n1 > 0_I4P .and. &
           load_samples_mismatches(sh, sd) == 0_I4P)
call load_samples_3d_bodies(two, x, x, x, sr, status=st, poses=[pa, pb], lengths=[0.1_R8P, 0.25_R8P], resident=.true.)
call sr%update_host
call check('twin: resident, update_host = host', st == STATUS_OK .and. sr%on_device .and. &
           load_samples_mismatches(sh, sr) == 0_I4P)
call load_samples_reset(sr)
call two(1)%exit_device ; call two(2)%exit_device
! 7. refusals
call make_sphere(2_I4P, fc) ; call octs%adopt_facets(facets=fc, aabb_tree_kind=AABB_TREE_OCTREE)
two(1) = sphere ; two(2) = octs
call load_samples_3d_bodies(two, samples=sb, status=st)
call check('refused: an octree surface (stderr message)', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)
two(2) = sphere
call load_samples_3d_bodies(two, samples=sb, status=st, poses=[pa])
call check('refused: size(poses) /= size(surfaces)', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)
call load_samples_3d_bodies(two(1:0), samples=sb, status=st)
call check('refused: no surface', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)
call load_samples_3d_bodies(two, samples=sb, status=st, device=.true.)
call check('refused: device path of unmapped surfaces (stderr message)', st == STATUS_INVALID_INPUT .and. sb%n == 0_I4P)

print '(A,L1)', 'Are all tests passed? ', ok

contains
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

   function strip(c) result(s)
   !< A bodies set without its several-bodies fields (a single-body set's shapes and flags).
   type(load_samples_t), intent(in) :: c !< Bodies set.
   type(load_samples_t)             :: s !< Stripped copy.

   s = c ; s%on_device = .false. ; s%device_key = 0
   s%with_bodies = .false. ; s%bodies = 0_I4P
   deallocate(s%sub_body, s%clearance, s%second_body, s%second_foot)
   allocate(s%sub_body(0), s%clearance(0), s%second_body(0), s%second_foot(c%dims, 0))
   endfunction strip

   logical function same_slice(c, off, a) result(same)
   !< The slots off+1 ... off+a%n of the bodies set c equal the single-body set a, every field by bits (body = the
   !< position of a's body in c, sub_body = a's body).
   type(load_samples_t), intent(in) :: c, a !< Bodies set, single-body set.
   integer(I4P),         intent(in) :: off  !< Offset.
   integer(I4P) :: n                        !< Samples.

   n = a%n
   same = eqr(reshape(c%points(:, off + 1:off + n), [3 * n]), reshape(a%points, [3 * n])) .and. &
          eqr(c%weight(off + 1:off + n), a%weight) .and. &
          eqr(reshape(c%normal(:, off + 1:off + n), [3 * n]), reshape(a%normal, [3 * n])) .and. &
          all(c%element(off + 1:off + n) == a%element) .and. all(c%patch(off + 1:off + n) == a%patch) .and. &
          all(c%sub_body(off + 1:off + n) == a%body) .and. &
          eqr(reshape(c%probe(:, :, off + 1:off + n), [size(a%probe)]), reshape(a%probe, [size(a%probe)])) .and. &
          all(c%visible(:, off + 1:off + n) .eqv. a%visible) .and. all(c%probe_cell(:, :, off + 1:off + n) == a%probe_cell)
   endfunction same_slice

   logical function eqr(a, b)
   !< Equal bits.
   real(R8P), intent(in) :: a(:), b(:) !< Reals.

   eqr = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
   endfunction eqr

   pure function cross(a, b) result(r)
   !< a x b.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: r(3)       !< Product.

   r = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross
endprogram fossil_test_load_samples_bodies_3d
