!< FOSSIL test: the closure table of a finite-difference wall closure on a 3-D block (#109 WP1 step 04).
program fossil_test_fd_closure_3d
!< FOSSIL test: the closure table of a finite-difference wall closure on a 3-D block (#109 WP1 step 04).
!<
!< 1. Hand checks on an icosphere of radius R = 1/2 (3 subdivisions, curvature built): every foot on the sphere and
!<    along the normal (|foot| within the facet sag of R, the normal within the facet tilt of foot/|foot|), the signed
!<    distance |p| - R to the sag, k1, k2 = 1/R within 4 %, the probes foot + l n (visible: a convex body), their cells by
!<    the low-side rule, the velocity V + omega x (foot - c) bit for bit, a pure rotation tangent to the wall to the
!<    normal tilt, every ILW slot with a crossing between its two points, the normal there leaving the solid along the
!<    walk.
!< 2. ILW against the cut cells: on the cut-cell nodes, the crossing of every slot lies on the cell edge between the
!<    slot's points and gives the edge's fluid fraction of `cut_cells_3d(edges=.true.)` bit for bit.
!< 3. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path in every field, bit for
!<    bit (`ib_closure_mismatches`): with curvature, probes, velocities and ILW; posed; without curvature and probes. The
!<    posed ghosts and ILW crossings (refitted tree, no moved copy since step 06) equal those of a surface moved by
!<    `move`, bit for bit.
!< 4. Resident lifecycle: `resident=.true.` flags the table, `update_host` gives the host path's bits, a copy is not
!<    mapped (its `update_host` warns and drops the flag), `exit_device` / `enter_device` / `exit_device`, a new call
!<    into a resident table releases it first, `ib_closure_reset`.
!< 5. Refusals: the device path of an unmapped surface or of a curvature built after `enter_device`, a bad reach size, a
!<    bad velocity size: STATUS_INVALID_INPUT, n = 0. An empty surface has no ghost.
use fossil,              only : surface_stl_object, ib_closure_t, fd_closure_3d, ib_closure_reset, ib_closure_mismatches, &
                                cut_cells_3d_t, cut_cells_3d, rigid_pose_t, pose_from_axis_angle, STATUS_OK, STATUS_INVALID_INPUT
use fossil_facet_object, only : facet_object
use fossil_thickness_bodies, only : make_boxes
use iso_fortran_env,     only : int64
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

real(R8P), parameter            :: R = 0.5_R8P
type(surface_stl_object)        :: sphere, plain, late, empty, moved, slabs
type(facet_object), allocatable :: fc(:)
type(ib_closure_t)              :: ch, cd, cr, cc
type(cut_cells_3d_t)            :: cut
type(rigid_pose_t)              :: pose
real(R8P), allocatable          :: x(:), y(:), xs(:), ys(:)
real(R8P)                       :: e, ev, vx, vy, vz, rr, elo, ehi, frac, ef
integer(I4P)                    :: i, g, k, d, a, st, nbad, nslot, lo, hi, s, sgn, idx(3)
logical                         :: ok

ok = .true.
call icosphere(3_I4P, R, fc) ; call sphere%adopt_facets(facets=fc) ; call sphere%build_curvature(status=st)
call check('sphere: curvature built', st == STATUS_OK)
call icosphere(3_I4P, R, fc) ; call plain%adopt_facets(facets=fc)
call check('plain sphere: 1280 facets', plain%get_facets_number() == 1280_I4P)
x = [(-0.73_R8P + 0.1_R8P * i, i = 0, 15)]
! 1. hand checks
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], ch, status=st, lengths=[0.1_R8P, 0.25_R8P], &
                   linear_velocity=[0.1_R8P, -0.2_R8P, 0.3_R8P], angular_velocity=[0.5_R8P, -1._R8P, 2._R8P], &
                   centre=[0.01_R8P, 0.02_R8P, -0.03_R8P], ilw=.true.)
call check('sphere: status, ghosts, flags', st == STATUS_OK .and. ch%n > 100_I4P .and. ch%dims == 3_I4P .and. ch%m == 2_I4P &
           .and. ch%with_curvature .and. ch%with_velocity .and. ch%with_ilw)
e = 0._R8P ; ev = 0._R8P
do g = 1_I4P, ch%n
   rr = norm2(ch%foot(:, g))
   e = max(e, abs(rr - R) / R)
   ev = max(ev, acos(min(1._R8P, dot_product(ch%normal(:, g), ch%foot(:, g)) / rr)))
enddo
call check('sphere: feet on the sphere (sag < 1 %)', e < 0.01_R8P)
call check('sphere: normal along the foot (tilt < 0.1 rad)', ev < 0.1_R8P)
e = 0._R8P
do g = 1_I4P, ch%n
   e = max(e, abs(ch%distance(g) - (norm2(ch%points(:, g)) - R)) / R)
enddo
call check('sphere: signed distance |p| - R (< 1 %)', e < 0.01_R8P .and. all(ch%distance < 0._R8P))
call check('sphere: k1, k2 = 1/R within 4 %', maxval(abs(ch%kappa * R - 1._R8P)) < 0.04_R8P)
nbad = 0_I4P
do g = 1_I4P, ch%n
   do k = 1_I4P, ch%m
      if (any(ch%probe(:, k, g) /= ch%foot(:, g) + ch%lengths(k) * ch%normal(:, g))) nbad = nbad + 1_I4P
      do a = 1_I4P, 3_I4P
         i = ch%probe_cell(a, k, g)
         if (i == 0_I4P) then
            if (ch%probe(a, k, g) > x(1) .and. ch%probe(a, k, g) <= x(size(x))) nbad = nbad + 1_I4P
         elseif (.not. (x(i) < ch%probe(a, k, g) .and. ch%probe(a, k, g) <= x(i + 1))) then
            nbad = nbad + 1_I4P
         endif
      enddo
   enddo
enddo
call check('sphere: probes foot + l n, cells by the low-side rule', nbad == 0_I4P)
call check('sphere: every probe visible (convex)', all(ch%visible))
nbad = 0_I4P
do g = 1_I4P, ch%n
   vx = 0.1_R8P + (-1._R8P * (ch%foot(3, g) - (-0.03_R8P)) - 2._R8P * (ch%foot(2, g) - 0.02_R8P))
   vy = -0.2_R8P + (2._R8P * (ch%foot(1, g) - 0.01_R8P) - 0.5_R8P * (ch%foot(3, g) - (-0.03_R8P)))
   vz = 0.3_R8P + (0.5_R8P * (ch%foot(2, g) - 0.02_R8P) - (-1._R8P) * (ch%foot(1, g) - 0.01_R8P))
   if (vx /= ch%velocity(1, g) .or. vy /= ch%velocity(2, g) .or. vz /= ch%velocity(3, g)) nbad = nbad + 1_I4P
enddo
call check('sphere: velocity V + omega x (foot - c), bit for bit', nbad == 0_I4P)
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, angular_velocity=[0.5_R8P, -1._R8P, 2._R8P])
e = 0._R8P
do g = 1_I4P, cd%n
   e = max(e, abs(dot_product(cd%velocity(:, g), cd%normal(:, g))) / norm2(cd%velocity(:, g)))
enddo
call check('sphere: a pure rotation about the centre is tangent (to the facet tilt)', st == STATUS_OK .and. e < 0.1_R8P)
nbad = 0_I4P ; nslot = 0_I4P
do g = 1_I4P, ch%n
   do d = 1_I4P, 6_I4P
      s = ch%steps(d, g)
      if (s < 1_I4P) then
         if (ch%ilw_element(d, g) /= 0_I4P) nbad = nbad + 1_I4P
         cycle
      endif
      nslot = nslot + 1_I4P
      call slot(ch, g, d, lo, hi, a, sgn)
      elo = x(lo) ; ehi = x(hi)
      if (ch%ilw_element(d, g) < 1_I4P) nbad = nbad + 1_I4P
      if (ch%ilw_t(d, g) < elo .or. ch%ilw_t(d, g) > ehi) nbad = nbad + 1_I4P
      if (ch%ilw_normal(a, d, g) * sgn <= 0._R8P) nbad = nbad + 1_I4P
      if (abs(ch%ilw_kappa(1, d, g) * R - 1._R8P) > 0.04_R8P) nbad = nbad + 1_I4P
   enddo
enddo
call check('sphere: ILW crossing in every slot, between its points, normal leaving along the walk, k = 1/R', &
           nbad == 0_I4P .and. nslot > ch%n)
! 2. ILW against the cut cells: the edge fractions of the slots' edges
call cut_cells_3d(sphere, x, x, x, cut, edges=.true.)
nbad = 0_I4P
do g = 1_I4P, ch%n
   do d = 1_I4P, 6_I4P
      if (ch%steps(d, g) < 1_I4P) cycle
      call slot(ch, g, d, lo, hi, a, sgn)
      elo = x(lo) ; ehi = x(hi)
      idx = ch%ijk(:, g) ; idx(a) = lo
      select case (a)
      case (1) ; frac = cut%edge_x(idx(1), idx(2), idx(3))
      case (2) ; frac = cut%edge_y(idx(1), idx(2), idx(3))
      case default ; frac = cut%edge_z(idx(1), idx(2), idx(3))
      endselect
      if (sgn > 0_I4P) then
         ef = (ehi - ch%ilw_t(d, g)) / (ehi - elo)
      else
         ef = (ch%ilw_t(d, g) - elo) / (ehi - elo)
      endif
      if (ef /= frac) nbad = nbad + 1_I4P
   enddo
enddo
call check('sphere: ILW crossings give the cut cells'' edge fractions bit for bit', nbad == 0_I4P)
! 3. twins
call sphere%enter_device
call fd_closure_3d(sphere, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, lengths=[0.1_R8P, 0.25_R8P], &
                   linear_velocity=[0.1_R8P, -0.2_R8P, 0.3_R8P], angular_velocity=[0.5_R8P, -1._R8P, 2._R8P], &
                   centre=[0.01_R8P, 0.02_R8P, -0.03_R8P], ilw=.true., device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. ib_closure_mismatches(ch, cd) == 0_I4P)
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.3_R8P, [0.05_R8P, -0.02_R8P, 0.03_R8P])
call fd_closure_3d(sphere, x, x, x, [1_I4P, 2_I4P, 3_I4P], ch, status=st, lengths=[-0.05_R8P, 0.1_R8P, 0.3_R8P], &
                   pose=pose, linear_velocity=[0.1_R8P, 0._R8P, 0._R8P], angular_velocity=[0._R8P, 0._R8P, 1._R8P], &
                   ilw=.true.)
call fd_closure_3d(sphere, x, x, x, [1_I4P, 2_I4P, 3_I4P], cd, status=st, lengths=[-0.05_R8P, 0.1_R8P, 0.3_R8P], &
                   pose=pose, linear_velocity=[0.1_R8P, 0._R8P, 0._R8P], angular_velocity=[0._R8P, 0._R8P, 1._R8P], &
                   ilw=.true., device=.true.)
call check('twin posed: device = host, every field', st == STATUS_OK .and. ch%n > 0_I4P .and. &
           ib_closure_mismatches(ch, cd) == 0_I4P)
e = 0._R8P
do g = 1_I4P, ch%n
   e = max(e, abs(norm2(ch%foot(:, g) - pose%t) - R) / R)
enddo
call check('posed: feet on the moved sphere, velocity about the pose''s origin', e < 0.01_R8P .and. &
           abs(ch%velocity(1, 1) - (0.1_R8P - (ch%foot(2, 1) - pose%t(2)))) < 1.e-15_R8P)
! the posed ghosts and crossings are those of the moved surface (step 06: refitted tree, no moved copy), bit for bit
call icosphere(3_I4P, R, fc) ; call moved%adopt_facets(facets=fc) ; call moved%move(pose)
call fd_closure_3d(moved, x, x, x, [1_I4P, 2_I4P, 3_I4P], cc, status=st, ilw=.true.)
call check('posed: ghosts and ILW crossings = those of a moved copy, bit for bit', cc%n == ch%n .and. &
           all(cc%ijk == ch%ijk) .and. all(cc%mask == ch%mask) .and. all(cc%steps == ch%steps) .and. &
           all(cc%ilw_element == ch%ilw_element) .and. &
           all(transfer(cc%points, 1_int64, size(cc%points)) == transfer(ch%points, 1_int64, size(ch%points))) .and. &
           all(transfer(cc%ilw_t, 1_int64, size(cc%ilw_t)) == transfer(ch%ilw_t, 1_int64, size(ch%ilw_t))))
call ib_closure_reset(cc)
! more crossings in an ILW window than the device's private buffer (16): twelve slabs 0.004 thick, 0.004 apart, nodes
! inside the first; the overflowed slots take the exact-size device path (step 06), the same table bit for bit
call make_boxes([0.4_R8P, -0.3_R8P, -0.3_R8P], [0.404_R8P, 0.3_R8P, 0.3_R8P], 12_I4P, 0.004_R8P, fc)
call slabs%adopt_facets(facets=fc) ; call slabs%enter_device
xs = [(0.002_R8P + 0.1_R8P * i, i = 0, 7)] ; ys = [(-0.25_R8P + 0.1_R8P * i, i = 0, 5)]
call fd_closure_3d(slabs, xs, ys, ys, [1_I4P, 1_I4P, 1_I4P], ch, status=st, ilw=.true.)
call fd_closure_3d(slabs, xs, ys, ys, [1_I4P, 1_I4P, 1_I4P], cd, status=st, ilw=.true., device=.true.)
call check('ILW: 24 crossings in a window (buffer overflow), device = host, every field', st == STATUS_OK .and. &
           ch%n == 36_I4P .and. all(ch%ilw_element(1:2, :) > 0_I4P) .and. ib_closure_mismatches(ch, cd) == 0_I4P)
call slabs%exit_device
call plain%enter_device
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P, 2_I4P], ch, status=st, ilw=.true.)
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, ilw=.true., device=.true.)
call check('twin without curvature and probes: device = host', st == STATUS_OK .and. .not. cd%with_curvature .and. &
           cd%m == 0_I4P .and. size(cd%kappa) == 0 .and. ib_closure_mismatches(ch, cd) == 0_I4P)
! 4. resident lifecycle
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P, 2_I4P], cr, status=st, ilw=.true., resident=.true.)
call check('resident: flagged', st == STATUS_OK .and. cr%on_device)
call cr%update_host
call check('resident: update_host gives the host bits', ib_closure_mismatches(ch, cr) == 0_I4P .and. cr%on_device)
cc = cr
call cc%update_host                                     ! a copy: warns on stderr, drops the flag
call check('resident: a copy is not mapped', .not. cc%on_device .and. ib_closure_mismatches(ch, cc) == 0_I4P)
call cr%exit_device
call check('resident: exit_device', .not. cr%on_device)
call cr%enter_device
call check('resident: enter_device', cr%on_device)
call cr%enter_device
call check('resident: enter_device twice is a no-op', cr%on_device)
call fd_closure_3d(plain, x, x, x, [1_I4P, 1_I4P, 1_I4P], cr, status=st, device=.true.)
call check('resident: a new call releases it first', st == STATUS_OK .and. .not. cr%on_device .and. cr%n > 0_I4P)
call fd_closure_3d(plain, x, x, x, [1_I4P, 1_I4P, 1_I4P], cr, status=st, resident=.true.)
call ib_closure_reset(cr)
call check('resident: ib_closure_reset', .not. cr%on_device .and. cr%n == 0_I4P .and. .not. allocated(cr%points))
! 5. refusals and the empty surface
call icosphere(3_I4P, R, fc) ; call late%adopt_facets(facets=fc) ; call late%enter_device ; call late%build_curvature(status=st)
call fd_closure_3d(late, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, device=.true.)
call check('refused: curvature built after enter_device (stderr message)', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call plain%exit_device
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, device=.true.)
call check('refused: device path of an unmapped surface (stderr message)', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P], cd, status=st)
call check('refused: size(reach) /= 3', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_3d(plain, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, linear_velocity=[1._R8P, 0._R8P])
call check('refused: size(linear_velocity) /= 3', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
y = [0._R8P]
call fd_closure_3d(plain, x, y, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st)
call check('refused: one point on an axis', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_3d(empty, x, x, x, [2_I4P, 2_I4P, 2_I4P], cd, status=st, device=.true., lengths=[0.1_R8P], ilw=.true.)
call check('empty surface: no ghost, both paths', st == STATUS_OK .and. cd%n == 0_I4P .and. allocated(cd%points))
call sphere%exit_device

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

   subroutine slot(c, g, d, lo, hi, a, sgn)
   !< The last solid and the first fluid point of slot d of ghost g: their indexes lo < hi along axis a, walk direction.
   type(ib_closure_t), intent(in)  :: c          !< Table.
   integer(I4P),       intent(in)  :: g, d       !< Ghost, slot.
   integer(I4P),       intent(out) :: lo, hi, a, sgn !< Indexes, axis, direction.
   integer(I4P) :: p0, p1                        !< Points.

   a = (d + 1_I4P) / 2_I4P ; sgn = 1_I4P ; if (mod(d, 2_I4P) == 1_I4P) sgn = -1_I4P
   p0 = c%ijk(a, g) + sgn * (c%steps(d, g) - 1_I4P) ; p1 = c%ijk(a, g) + sgn * c%steps(d, g)
   lo = min(p0, p1) ; hi = max(p0, p1)
   endsubroutine slot

   subroutine put(fo, k, a, b, cc)
   !< Facet k from three points.
   type(facet_object), intent(inout) :: fo(:)             !< Facets.
   integer(I4P),       intent(in)    :: k                 !< Index.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo(k)%vertex(1) = vector_R8P(a(1), a(2), a(3))
   fo(k)%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo(k)%vertex(3) = vector_R8P(cc(1), cc(2), cc(3))
   call fo(k)%compute_metrix
   endsubroutine put

   subroutine icosphere(levels, radius, fo)
   !< Icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere.
   integer(I4P),                    intent(in)  :: levels !< Subdivisions.
   real(R8P),                       intent(in)  :: radius !< Radius.
   type(facet_object), allocatable, intent(out) :: fo(:)  !< Facets.
   real(R8P), allocatable  :: tri(:,:,:), nxt(:,:,:)
   real(R8P)               :: vv(3, 12), gr
   integer(I4P), parameter :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                    2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                    4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                    5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   integer(I4P)            :: k, l, m

   gr = (1._R8P + sqrt(5._R8P)) / 2._R8P
   vv = reshape([-1._R8P, gr, 0._R8P, 1._R8P, gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P, -gr, 0._R8P, &
                 0._R8P, -1._R8P, gr, 0._R8P, 1._R8P, gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P, -gr, &
                 gr, 0._R8P, -1._R8P, gr, 0._R8P, 1._R8P, -gr, 0._R8P, -1._R8P, -gr, 0._R8P, 1._R8P], [3, 12])
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
      call put(fo, k, radius * tri(:, 1, k), radius * tri(:, 2, k), radius * tri(:, 3, k))
   enddo
   endsubroutine icosphere

   pure function mid(a, b) result(m)
   !< Unit midpoint.
   real(R8P), intent(in) :: a(3), b(3) !< Points.
   real(R8P)             :: m(3)       !< Midpoint on the unit sphere.

   m = (a + b) / norm2(a + b)
   endfunction mid
endprogram fossil_test_fd_closure_3d
