!< FOSSIL test: the closure table of several bodies on a 3-D block (#109 WP3 step 01).
program fossil_test_fd_closure_bodies_3d
!< FOSSIL test: the closure table of several bodies on a 3-D block (#109 WP3 step 01).
!<
!< 1. One body: `fd_closure_3d_bodies` equals `fd_closure_3d` bit for bit in every single-body field (with curvature,
!<    probes, velocities and ILW; posed; on the device path), body = 1, sub_body = the facet_body, no second body, no
!<    shared or overlapping ghost; an empty second surface is ignored.
!< 2. Two unit spheres far apart (posed, each with its own velocity): the ghosts are the union of each body's, and every
!<    ghost equals the single-body table of its body, field by field, bit for bit; the clearance is |foot - c_other| - 1
!<    to the facet sag, the second body the other sphere.
!< 3. Narrow gap: two boxes 1/8 apart (one grid point in the gap): the feet on the facing faces have clearance 1/8 to
!<    round-off, the second foot on the other face, shared (default shared distance 2 h = 0.2); a probe longer than the
!<    gap enters the other box and is not visible, a shorter one is; shared_distance = 0.1 shares none.
!< 4. Overlap: two spheres 0.8 apart: points inside both are counted, the overlapping ghosts lie inside both, each
!<    ghost's |distance| is the smaller of the two.
!< 5. Twin: posed bodies with velocities, probes and ILW: `device=.true.` (built without OpenACC its loops run on the
!<    host) equals the host path in every field, bit for bit; resident, then `update_host`; the velocity of every ghost
!<    is V_b + omega_b x (foot - c_b) of its owner b, bit for bit.
!< 6. Refusals: an octree surface, size(poses) /= size(surfaces), a velocity array of the wrong shape, no surface, the
!<    device path of an unmapped surface: STATUS_INVALID_INPUT, n = 0.
use fossil,                  only : surface_stl_object, ib_closure_t, fd_closure_3d, fd_closure_3d_bodies, ib_closure_reset, &
                                    ib_closure_mismatches, rigid_pose_t, pose_from_axis_angle, AABB_TREE_OCTREE, STATUS_OK, &
                                    STATUS_INVALID_INPUT
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_boxes, make_sphere
use iso_fortran_env,         only : int64
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: G8 = 0.125_R8P           !< Box gap.
type(surface_stl_object)        :: sphere, boxa, boxb, octs, empty
type(surface_stl_object)        :: two(2), boxes(2), one(1)
type(facet_object), allocatable :: fc(:)
type(ib_closure_t)              :: ch, cb, cd, ca, cc, cr
type(rigid_pose_t)              :: pa, pb, pp(2)
type(vector_R8P), allocatable   :: qp(:)
real(R8P), allocatable          :: x(:), xb(:), yb(:), qd(:)
real(R8P)                       :: e, v(3), c(3), vl(3, 2), wl(3, 2), cl(3, 2)
integer(I4P)                    :: g, k, st, nbad, nf, ns, b
logical                         :: ok

ok = .true.
call make_sphere(2_I4P, fc) ; call sphere%adopt_facets(facets=fc) ; call sphere%build_curvature(status=st)
call check('sphere: 320 facets, curvature built', st == STATUS_OK .and. sphere%get_facets_number() == 320_I4P)
x = [(-2.45_R8P + 0.2_R8P * g, g = 0, 24)]
vl = reshape([0.1_R8P, -0.2_R8P, 0.3_R8P, -0.3_R8P, 0.05_R8P, 0.2_R8P], [3, 2])
wl = reshape([0.5_R8P, -1._R8P, 2._R8P, -0.7_R8P, 0.4_R8P, 0.1_R8P], [3, 2])
cl = reshape([0.01_R8P, 0.02_R8P, -0.03_R8P, 0.2_R8P, -0.1_R8P, 0.05_R8P], [3, 2])
! 1. one body
one(1) = sphere
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], ch, status=st, lengths=[0.2_R8P, 0.5_R8P], &
                   linear_velocity=vl(:, 1), angular_velocity=wl(:, 1), centre=cl(:, 1), ilw=.true.)
call fd_closure_3d_bodies(one, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, lengths=[0.2_R8P, 0.5_R8P], &
                          linear_velocities=vl(:, 1:1), angular_velocities=wl(:, 1:1), centres=cl(:, 1:1), ilw=.true.)
call check('one body: status, flags', st == STATUS_OK .and. cb%n == ch%n .and. ch%n > 100_I4P .and. cb%with_bodies .and. &
           cb%bodies == 1_I4P .and. cb%with_curvature .and. cb%with_velocity .and. cb%with_ilw)
call check('one body: = fd_closure_3d bit for bit in every single-body field', &
           ib_closure_mismatches(ch, strip(cb)) == 0_I4P)
call check('one body: body 1, sub_body the facet_body, no second body, nothing shared', all(cb%body == 1_I4P) .and. &
           all(cb%sub_body == ch%body) .and. all(cb%second_body == 0_I4P) .and. all(cb%clearance == huge(1._R8P)) .and. &
           .not. any(cb%shared) .and. .not. any(cb%overlap) .and. cb%overlap_points == 0_I4P .and. &
           all((cb%ilw_body == 1_I4P) .eqv. (cb%ilw_element > 0_I4P)) .and. all(cb%ilw_body >= 0_I4P))
pa = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.3_R8P, [0.05_R8P, -0.02_R8P, 0.03_R8P])
call fd_closure_3d(sphere, x, x, x, [1_I4P, 2_I4P, 3_I4P], ch, status=st, lengths=[-0.1_R8P, 0.3_R8P], pose=pa, &
                   linear_velocity=vl(:, 1), angular_velocity=wl(:, 1), ilw=.true.)
call fd_closure_3d_bodies(one, x, x, x, [1_I4P, 2_I4P, 3_I4P], cb, status=st, lengths=[-0.1_R8P, 0.3_R8P], poses=[pa], &
                          linear_velocities=vl(:, 1:1), angular_velocities=wl(:, 1:1), ilw=.true.)
call check('one body posed: = fd_closure_3d bit for bit', st == STATUS_OK .and. cb%n > 0_I4P .and. &
           ib_closure_mismatches(ch, strip(cb)) == 0_I4P)
call sphere%enter_device ; call one(1)%enter_device ! map the array elements themselves (a copy is not mapped)
call fd_closure_3d(sphere, x, x, x, [1_I4P, 2_I4P, 3_I4P], cd, status=st, lengths=[-0.1_R8P, 0.3_R8P], pose=pa, &
                   linear_velocity=vl(:, 1), angular_velocity=wl(:, 1), ilw=.true., device=.true.)
call fd_closure_3d_bodies(one, x, x, x, [1_I4P, 2_I4P, 3_I4P], cc, status=st, lengths=[-0.1_R8P, 0.3_R8P], poses=[pa], &
                          linear_velocities=vl(:, 1:1), angular_velocities=wl(:, 1:1), ilw=.true., device=.true.)
call check('one body posed, device path: = fd_closure_3d(device) and the host bodies table, bit for bit', &
           st == STATUS_OK .and. ib_closure_mismatches(cd, strip(cc)) == 0_I4P .and. ib_closure_mismatches(cb, cc) == 0_I4P)
call sphere%exit_device ; call one(1)%exit_device
two(1) = sphere ; two(2) = empty
call fd_closure_3d_bodies(two, x, x, x, [1_I4P, 2_I4P, 3_I4P], cc, status=st, lengths=[-0.1_R8P, 0.3_R8P], &
                          poses=[pa, pa], linear_velocities=vl, angular_velocities=wl, ilw=.true.)
call check('an empty second surface is ignored', st == STATUS_OK .and. cc%bodies == 2_I4P .and. &
           ib_closure_mismatches(ch, strip(cc)) == 0_I4P .and. all(cc%second_body == 0_I4P))
! 2. two spheres far apart (gap 0.6 > 2 h)
pa = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.2_R8P, [-1.3_R8P, 0.03_R8P, 0.01_R8P])
pb = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], -0.3_R8P, [1.3_R8P, -0.02_R8P, 0.04_R8P])
two(1) = sphere ; two(2) = sphere
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, lengths=[0.2_R8P, 0.5_R8P], poses=[pa, pb], &
                          linear_velocities=vl, angular_velocities=wl, centres=cl, ilw=.true.)
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], ca, lengths=[0.2_R8P, 0.5_R8P], pose=pa, &
                   linear_velocity=vl(:, 1), angular_velocity=wl(:, 1), centre=cl(:, 1), ilw=.true.)
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], cc, lengths=[0.2_R8P, 0.5_R8P], pose=pb, &
                   linear_velocity=vl(:, 2), angular_velocity=wl(:, 2), centre=cl(:, 2), ilw=.true.)
call check('far apart: the ghosts are the union of each body''s', st == STATUS_OK .and. cb%n == ca%n + cc%n .and. &
           count(cb%body == 1_I4P) == ca%n .and. count(cb%body == 2_I4P) == cc%n .and. cb%overlap_points == 0_I4P)
nbad = 0_I4P
do g = 1_I4P, cb%n
   if (cb%body(g) == 1_I4P) then
      if (.not. same_ghost(cb, g, ca, 1_I4P)) nbad = nbad + 1_I4P
   else
      if (.not. same_ghost(cb, g, cc, 2_I4P)) nbad = nbad + 1_I4P
   endif
enddo
call check('far apart: every ghost = the single-body table of its body, bit for bit', nbad == 0_I4P)
e = 0._R8P ; nbad = 0_I4P
do g = 1_I4P, cb%n
   c = pb%t ; if (cb%body(g) == 2_I4P) c = pa%t
   e = max(e, abs(cb%clearance(g) - (norm2(cb%foot(:, g) - c) - 1._R8P)))
   if (cb%second_body(g) /= 3_I4P - cb%body(g)) nbad = nbad + 1_I4P
   if (abs(norm2(cb%second_foot(:, g) - c) - 1._R8P) > 0.02_R8P) nbad = nbad + 1_I4P
enddo
call check('far apart: clearance |foot - c_other| - 1 (facet sag), second body and foot', e < 0.03_R8P .and. &
           nbad == 0_I4P .and. .not. any(cb%shared) .and. abs(cb%shared_distance - 0.4_R8P) < 1.e-12_R8P)
! 3. narrow gap: two boxes 1/8 apart
call make_boxes([-1._R8P, -0.6_R8P, -0.6_R8P], [-0.5_R8P * G8, 0.6_R8P, 0.6_R8P], 1_I4P, 0._R8P, fc)
call boxa%adopt_facets(facets=fc)
call make_boxes([0.5_R8P * G8, -0.6_R8P, -0.6_R8P], [1._R8P, 0.6_R8P, 0.6_R8P], 1_I4P, 0._R8P, fc)
call boxb%adopt_facets(facets=fc)
boxes(1) = boxa ; boxes(2) = boxb
xb = [(-1.25_R8P + 0.1_R8P * g, g = 0, 25)] ; yb = [(-0.83_R8P + 0.1_R8P * g, g = 0, 16)]
call fd_closure_3d_bodies(boxes, xb, yb, yb, [2_I4P, 2_I4P, 2_I4P], cb, status=st, lengths=[0.1_R8P, 0.2_R8P])
call fd_closure_3d(boxa, xb, yb, yb, [2_I4P, 2_I4P, 2_I4P], ca, lengths=[0.1_R8P, 0.2_R8P])
nf = 0_I4P ; nbad = 0_I4P ; ns = 0_I4P
do g = 1_I4P, cb%n
   if (cb%shared(g) .neqv. (cb%clearance(g) < cb%shared_distance)) nbad = nbad + 1_I4P
   if (cb%shared(g)) ns = ns + 1_I4P
   if (cb%body(g) /= 1_I4P .or. abs(cb%foot(1, g) + 0.5_R8P * G8) > 1.e-14_R8P) cycle
   if (abs(cb%foot(2, g)) > 0.55_R8P .or. abs(cb%foot(3, g)) > 0.55_R8P) cycle
   nf = nf + 1_I4P
   if (abs(cb%clearance(g) - G8) > 1.e-14_R8P .or. .not. cb%shared(g) .or. cb%second_body(g) /= 2_I4P) nbad = nbad + 1_I4P
   if (abs(cb%second_foot(1, g) - 0.5_R8P * G8) > 1.e-14_R8P .or. &
       any(abs(cb%second_foot(2:3, g) - cb%foot(2:3, g)) > 1.e-14_R8P)) nbad = nbad + 1_I4P
   if (.not. cb%visible(1, g) .or. cb%visible(2, g)) nbad = nbad + 1_I4P
   k = ghost_of(ca, cb%ijk(:, g))
   if (k < 1_I4P) then
      nbad = nbad + 1_I4P
   elseif (.not. all(ca%visible(:, k))) then
      nbad = nbad + 1_I4P
   endif
enddo
call check('narrow gap: facing feet, clearance 1/8, second foot, shared, the long probe enters the other box', &
           st == STATUS_OK .and. nf > 20_I4P .and. nbad == 0_I4P .and. ns > nf .and. ns < cb%n)
call fd_closure_3d_bodies(boxes, xb, yb, yb, [2_I4P, 2_I4P, 2_I4P], cb, status=st, shared_distance=0.1_R8P)
call check('narrow gap: shared_distance 0.1 < 1/8 shares none', st == STATUS_OK .and. .not. any(cb%shared))
! 4. overlap: two spheres 0.8 apart
pp(1) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [-0.4_R8P, 0._R8P, 0._R8P])
pp(2) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [0.4_R8P, 0._R8P, 0._R8P])
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, poses=pp)
nbad = 0_I4P
allocate(qp(cb%n), qd(cb%n))
do g = 1_I4P, cb%n
   qp(g) = vector_R8P(cb%points(1, g), cb%points(2, g), cb%points(3, g))
   if (cb%overlap(g)) then
      if (norm2(cb%points(:, g) - pp(1)%t) > 1.02_R8P .or. norm2(cb%points(:, g) - pp(2)%t) > 1.02_R8P) nbad = nbad + 1_I4P
   else
      if (norm2(cb%points(:, g) - pp(1)%t) < 0.98_R8P .and. norm2(cb%points(:, g) - pp(2)%t) < 0.98_R8P) nbad = nbad + 1_I4P
   endif
enddo
do b = 1_I4P, 2_I4P
   call sphere%query_many(points=qp, distances=qd, pose=pp(b))
   do g = 1_I4P, cb%n
      if (abs(cb%distance(g)) > abs(qd(g))) nbad = nbad + 1_I4P
   enddo
enddo
call check('overlap: points inside both counted, overlapping ghosts inside both, |d| the smaller', st == STATUS_OK .and. &
           cb%overlap_points > count(cb%overlap) .and. count(cb%overlap) > 10 .and. nbad == 0_I4P)
! 5. twin: posed bodies 0.15 apart with velocities, probes and ILW
pa = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], 0.4_R8P, [-1.075_R8P, 0.02_R8P, -0.01_R8P])
pb = pose_from_axis_angle([1._R8P, 0._R8P, -1._R8P], -0.2_R8P, [1.075_R8P, -0.03_R8P, 0.02_R8P])
two(1) = sphere ; two(2) = sphere ; call two(1)%enter_device ; call two(2)%enter_device
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 1_I4P, 2_I4P], ch, status=st, lengths=[0.1_R8P, 0.25_R8P], &
                          poses=[pa, pb], linear_velocities=vl, angular_velocities=wl, ilw=.true.)
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 1_I4P, 2_I4P], cd, status=st, lengths=[0.1_R8P, 0.25_R8P], &
                          poses=[pa, pb], linear_velocities=vl, angular_velocities=wl, ilw=.true., device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. ch%n > 0_I4P .and. any(ch%shared) .and. &
           .not. all(ch%visible) .and. any(ch%ilw_body == 2_I4P) .and. ib_closure_mismatches(ch, cd) == 0_I4P)
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 1_I4P, 2_I4P], cr, status=st, lengths=[0.1_R8P, 0.25_R8P], &
                          poses=[pa, pb], linear_velocities=vl, angular_velocities=wl, ilw=.true., resident=.true.)
call cr%update_host
call check('twin: resident, update_host = host', st == STATUS_OK .and. cr%on_device .and. &
           ib_closure_mismatches(ch, cr) == 0_I4P)
call ib_closure_reset(cr)
nbad = 0_I4P
do g = 1_I4P, ch%n
   b = ch%body(g) ; c = pa%t ; if (b == 2_I4P) c = pb%t
   v(1) = vl(1, b) + (wl(2, b) * (ch%foot(3, g) - c(3)) - wl(3, b) * (ch%foot(2, g) - c(2)))
   v(2) = vl(2, b) + (wl(3, b) * (ch%foot(1, g) - c(1)) - wl(1, b) * (ch%foot(3, g) - c(3)))
   v(3) = vl(3, b) + (wl(1, b) * (ch%foot(2, g) - c(2)) - wl(2, b) * (ch%foot(1, g) - c(1)))
   if (any(v /= ch%velocity(:, g))) nbad = nbad + 1_I4P
enddo
call check('twin: velocity V_b + omega_b x (foot - c_b) of the owner, bit for bit', nbad == 0_I4P)
call two(1)%exit_device ; call two(2)%exit_device
! 6. refusals
call make_sphere(2_I4P, fc) ; call octs%adopt_facets(facets=fc, aabb_tree_kind=AABB_TREE_OCTREE)
two(1) = sphere ; two(2) = octs
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st)
call check('refused: an octree surface (stderr message)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
two(2) = sphere
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, poses=[pa])
call check('refused: size(poses) /= size(surfaces)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, linear_velocities=vl(:, 1:1))
call check('refused: linear_velocities of shape [3, 1] for two bodies', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_3d_bodies(two(1:0), x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st)
call check('refused: no surface', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)
call fd_closure_3d_bodies(two, x, x, x, [2_I4P, 2_I4P, 2_I4P], cb, status=st, device=.true.)
call check('refused: device path of unmapped surfaces (stderr message)', st == STATUS_INVALID_INPUT .and. cb%n == 0_I4P)

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
   !< A bodies table without its several-bodies fields (a single-body table's shapes and flags).
   type(ib_closure_t), intent(in) :: c !< Bodies table.
   type(ib_closure_t)             :: s !< Stripped copy.

   s = c ; s%on_device = .false. ; s%device_key = 0
   s%with_bodies = .false. ; s%bodies = 0_I4P ; s%overlap_points = 0_I4P ; s%shared_distance = 0._R8P
   deallocate(s%sub_body, s%clearance, s%second_body, s%second_foot, s%shared, s%overlap, s%ilw_body)
   allocate(s%sub_body(0), s%clearance(0), s%second_body(0), s%second_foot(c%dims, 0), s%shared(0), s%overlap(0), &
            s%ilw_body(2 * c%dims, 0))
   endfunction strip

   integer(I4P) function ghost_of(c, ijk) result(k)
   !< The ghost of table c at indexes ijk (0: none).
   type(ib_closure_t), intent(in) :: c      !< Table.
   integer(I4P),       intent(in) :: ijk(:) !< Indexes.
   integer(I4P) :: g                        !< Counter.

   k = 0_I4P
   do g = 1_I4P, c%n
      if (all(c%ijk(:, g) == ijk)) then
         k = g ; return
      endif
   enddo
   endfunction ghost_of

   logical function same_ghost(cb, g, cs, b) result(same)
   !< Ghost g of the bodies table cb equals the ghost at the same indexes of body b's single-body table cs, every field
   !< by bits (body = b, sub_body = cs's body; ilw_body = b where a crossing).
   type(ib_closure_t), intent(in) :: cb, cs !< Bodies and single-body tables.
   integer(I4P),       intent(in) :: g, b   !< Ghost, body.
   integer(I4P) :: s                        !< Single-body ghost.

   same = .false.
   s = ghost_of(cs, cb%ijk(:, g))
   if (s < 1_I4P) return
   same = eqr(cb%points(:, g), cs%points(:, s)) .and. cb%mask(g) == cs%mask(s) .and. all(cb%steps(:, g) == cs%steps(:, s))
   same = same .and. eqr([cb%distance(g)], [cs%distance(s)]) .and. eqr(cb%foot(:, g), cs%foot(:, s)) .and. &
          eqr(cb%normal(:, g), cs%normal(:, s)) .and. cb%element(g) == cs%element(s) .and. &
          cb%feature(g) == cs%feature(s) .and. cb%patch(g) == cs%patch(s) .and. cb%body(g) == b .and. &
          cb%sub_body(g) == cs%body(s)
   same = same .and. eqr(reshape(cb%probe(:, :, g), [size(cb%probe(:, :, g))]), &
                         reshape(cs%probe(:, :, s), [size(cs%probe(:, :, s))])) .and. &
          all(cb%visible(:, g) .eqv. cs%visible(:, s)) .and. all(cb%probe_cell(:, :, g) == cs%probe_cell(:, :, s))
   same = same .and. eqr(cb%kappa(:, g), cs%kappa(:, s)) .and. &
          eqr(reshape(cb%direction(:, :, g), [6]), reshape(cs%direction(:, :, s), [6])) .and. &
          eqr(cb%velocity(:, g), cs%velocity(:, s))
   same = same .and. eqr(cb%ilw_t(:, g), cs%ilw_t(:, s)) .and. all(cb%ilw_element(:, g) == cs%ilw_element(:, s)) .and. &
          eqr(reshape(cb%ilw_normal(:, :, g), [18]), reshape(cs%ilw_normal(:, :, s), [18])) .and. &
          eqr(reshape(cb%ilw_kappa(:, :, g), [12]), reshape(cs%ilw_kappa(:, :, s), [12])) .and. &
          eqr(reshape(cb%ilw_direction(:, :, :, g), [36]), reshape(cs%ilw_direction(:, :, :, s), [36])) .and. &
          all((cb%ilw_body(:, g) == b) .eqv. (cs%ilw_element(:, s) > 0_I4P))
   endfunction same_ghost

   logical function eqr(a, b)
   !< Equal bits.
   real(R8P), intent(in) :: a(:), b(:) !< Reals.

   eqr = all(transfer(a, 1_int64, size(a)) == transfer(b, 1_int64, size(b)))
   endfunction eqr
endprogram fossil_test_fd_closure_bodies_3d
