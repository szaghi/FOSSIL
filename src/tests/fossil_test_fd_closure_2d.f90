!< FOSSIL test: the closure table of a finite-difference wall closure on a 2-D block (#109 WP1 step 04).
program fossil_test_fd_closure_2d
!< FOSSIL test: the closure table of a finite-difference wall closure on a 2-D block (#109 WP1 step 04).
!<
!< 1. Hand checks on a circle of radius R = 1/2 (128 uniform vertices, curvature built): feet on the circle to the chord
!<    sag, the normal along foot/|foot|, the signed distance |p| - R, the curvature 1/R to round-off (uniform circle),
!<    the tangent (-n_v, n_u), probes foot + l n with their cells, the velocity V + omega (-r_v, r_u) bit for bit, a pure
!<    rotation tangent to the wall within h / (2R), every ILW slot with a crossing between its two points and the normal
!<    leaving the solid along the walk.
!< 2. ILW against the cut cells: on the cut-cell nodes, the crossing of every slot gives the fluid fraction of the cell
!<    edge between the slot's points (`curve%cut_cells` apertures) to round-off.
!< 3. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path in every field, bit for
!<    bit, on the circle (all options), an ellipse posed with velocities, and the circle without curvature.
!< 4. Resident lifecycle as in 3-D; refusals; an empty curve has no ghost.
use fossil, only : curve_object, cut_cells_2d_t, ib_closure_t, fd_closure_2d, ib_closure_reset, ib_closure_mismatches, &
                   rigid_pose_2d_t, pose_2d, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P), R = 0.5_R8P
integer(I4P), parameter :: NV = 128_I4P
type(curve_object)    :: circle, ellipse, plain, empty
type(ib_closure_t)    :: ch, cd, cr, cc
type(cut_cells_2d_t)  :: cut
type(rigid_pose_2d_t) :: pose
real(R8P)             :: u(NV), v(NV), e, h, vu, vv, elo, ehi, ef, ap
real(R8P), allocatable :: x(:), y(:)
integer(I4P)          :: i, g, k, d, a, st, nbad, nslot, lo, hi, sgn, ic, jc
logical               :: ok

ok = .true.
u = [(R * cos(2 * PI * i / NV), i = 0, NV - 1)] ; v = [(R * sin(2 * PI * i / NV), i = 0, NV - 1)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[NV]) ; call circle%build_curvature(status=st)
call check('circle: curvature built', st == STATUS_OK)
call plain%adopt_loops(u=u, v=v, loop_sizes=[NV])
u = [(0.02_R8P + 0.6_R8P * cos(2 * PI * i / NV), i = 0, NV - 1)] ; v = [(-0.01_R8P + 0.3_R8P * sin(2 * PI * i / NV), i = 0, NV - 1)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[NV]) ; call ellipse%build_curvature(status=st)
h = 2._R8P * R * sin(PI / NV)
x = [(-0.73_R8P + 0.05_R8P * i, i = 0, 30)] ; y = [(-0.71_R8P + 0.05_R8P * i, i = 0, 29)]
! 1. hand checks
call fd_closure_2d(circle, x, y, [2_I4P, 3_I4P], ch, status=st, lengths=[0.05_R8P, 0.1_R8P, 0.15_R8P], &
                   linear_velocity=[0.1_R8P, -0.2_R8P], angular_velocity=2._R8P, centre=[0.01_R8P, -0.02_R8P], ilw=.true.)
call check('circle: status, ghosts, flags', st == STATUS_OK .and. ch%n > 50_I4P .and. ch%dims == 2_I4P .and. &
           ch%nk == 1_I4P .and. ch%m == 3_I4P .and. ch%with_curvature .and. ch%with_velocity .and. ch%with_ilw)
e = 0._R8P
do g = 1_I4P, ch%n
   e = max(e, abs(norm2(ch%foot(:, g)) - R))
enddo
call check('circle: feet on the circle (chord sag)', e <= R * (1._R8P - cos(PI / NV)) * 1.0001_R8P)
e = 0._R8P
do g = 1_I4P, ch%n
   e = max(e, abs(ch%distance(g) - (norm2(ch%points(:, g)) - R)))
enddo
call check('circle: signed distance |p| - R (to the sag)', e <= 1.0001_R8P * R * (1._R8P - cos(PI / NV)) .and. &
           all(ch%distance < 0._R8P))
call check('circle: curvature 1/R', maxval(abs(ch%kappa(1, :) * R - 1._R8P)) < 1.e-12_R8P)
call check('circle: tangent (-n_v, n_u)', all(ch%direction(1, 1, :) == -ch%normal(2, :)) .and. &
           all(ch%direction(2, 1, :) == ch%normal(1, :)))
nbad = 0_I4P
do g = 1_I4P, ch%n
   do k = 1_I4P, ch%m
      if (any(ch%probe(:, k, g) /= ch%foot(:, g) + ch%lengths(k) * ch%normal(:, g))) nbad = nbad + 1_I4P
      i = ch%probe_cell(1, k, g)
      if (i > 0_I4P) then
         if (.not. (x(i) < ch%probe(1, k, g) .and. ch%probe(1, k, g) <= x(i + 1))) nbad = nbad + 1_I4P
      elseif (ch%probe(1, k, g) > x(1) .and. ch%probe(1, k, g) <= x(size(x))) then
         nbad = nbad + 1_I4P
      endif
      i = ch%probe_cell(2, k, g)
      if (i > 0_I4P) then
         if (.not. (y(i) < ch%probe(2, k, g) .and. ch%probe(2, k, g) <= y(i + 1))) nbad = nbad + 1_I4P
      elseif (ch%probe(2, k, g) > y(1) .and. ch%probe(2, k, g) <= y(size(y))) then
         nbad = nbad + 1_I4P
      endif
   enddo
enddo
call check('circle: probes foot + l n, cells by the low-side rule', nbad == 0_I4P .and. all(ch%visible))
nbad = 0_I4P
do g = 1_I4P, ch%n
   vu = 0.1_R8P - 2._R8P * (ch%foot(2, g) - (-0.02_R8P))
   vv = -0.2_R8P + 2._R8P * (ch%foot(1, g) - 0.01_R8P)
   if (vu /= ch%velocity(1, g) .or. vv /= ch%velocity(2, g)) nbad = nbad + 1_I4P
enddo
call check('circle: velocity V + omega (-r_v, r_u), bit for bit', nbad == 0_I4P)
call fd_closure_2d(circle, x, y, [2_I4P, 3_I4P], cd, status=st, angular_velocity=2._R8P)
e = 0._R8P
do g = 1_I4P, cd%n
   e = max(e, abs(dot_product(cd%velocity(:, g), cd%normal(:, g))) / norm2(cd%velocity(:, g)))
enddo
call check('circle: a pure rotation is tangent within h/(2R)', st == STATUS_OK .and. e <= 0.5_R8P * h / R * 1.0001_R8P)
nbad = 0_I4P ; nslot = 0_I4P
do g = 1_I4P, ch%n
   do d = 1_I4P, 4_I4P
      if (ch%steps(d, g) < 1_I4P) then
         if (ch%ilw_element(d, g) /= 0_I4P) nbad = nbad + 1_I4P
         cycle
      endif
      nslot = nslot + 1_I4P
      call slot(ch, g, d, lo, hi, a, sgn)
      if (a == 1_I4P) then
         elo = x(lo) ; ehi = x(hi)
      else
         elo = y(lo) ; ehi = y(hi)
      endif
      if (ch%ilw_element(d, g) < 1_I4P) nbad = nbad + 1_I4P
      if (ch%ilw_t(d, g) < elo .or. ch%ilw_t(d, g) > ehi) nbad = nbad + 1_I4P
      if (ch%ilw_normal(a, d, g) * sgn <= 0._R8P) nbad = nbad + 1_I4P
      if (abs(ch%ilw_kappa(1, d, g) * R - 1._R8P) > 1.e-12_R8P) nbad = nbad + 1_I4P
   enddo
enddo
call check('circle: ILW crossing in every slot, between its points, normal leaving along the walk, k = 1/R', &
           nbad == 0_I4P .and. nslot > ch%n)
! 2. ILW against the cut cells (lines = cut-cell nodes)
call circle%cut_cells(x, y, cut, status=st)
nbad = 0_I4P
do g = 1_I4P, ch%n
   do d = 1_I4P, 4_I4P
      if (ch%steps(d, g) < 1_I4P) cycle
      call slot(ch, g, d, lo, hi, a, sgn)
      if (a == 1_I4P) then ! the edge along u at row j, from x(lo): -v of cell (lo, j) or +v of (lo, j - 1)
         elo = x(lo) ; ehi = x(hi) ; ic = lo ; jc = ch%ijk(2, g)
         if (jc <= size(y) - 1) then
            ap = cut%aperture(3, ic, jc)
         else
            ap = cut%aperture(4, ic, jc - 1)
         endif
      else                 ! the edge along v at column i, from y(lo): -u of cell (i, lo) or +u of (i - 1, lo)
         elo = y(lo) ; ehi = y(hi) ; ic = ch%ijk(1, g) ; jc = lo
         if (ic <= size(x) - 1) then
            ap = cut%aperture(1, ic, jc)
         else
            ap = cut%aperture(2, ic - 1, jc)
         endif
      endif
      if (sgn > 0_I4P) then
         ef = (ehi - ch%ilw_t(d, g)) / (ehi - elo)
      else
         ef = (ch%ilw_t(d, g) - elo) / (ehi - elo)
      endif
      if (abs(ef - ap) > 1.e-12_R8P) nbad = nbad + 1_I4P
   enddo
enddo
call check('circle: ILW crossings give the cut cells'' edge apertures', st == STATUS_OK .and. nbad == 0_I4P)
! 3. twins
call circle%enter_device ; call ellipse%enter_device ; call plain%enter_device
call fd_closure_2d(circle, x, y, [2_I4P, 3_I4P], cd, status=st, lengths=[0.05_R8P, 0.1_R8P, 0.15_R8P], &
                   linear_velocity=[0.1_R8P, -0.2_R8P], angular_velocity=2._R8P, centre=[0.01_R8P, -0.02_R8P], ilw=.true., &
                   device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. ib_closure_mismatches(ch, cd) == 0_I4P)
pose = pose_2d(0.4_R8P, [0.03_R8P, -0.05_R8P])
call fd_closure_2d(ellipse, x, y, [3_I4P, 1_I4P], ch, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pose, &
                   linear_velocity=[0.1_R8P, 0._R8P], angular_velocity=-1._R8P, ilw=.true.)
call fd_closure_2d(ellipse, x, y, [3_I4P, 1_I4P], cd, status=st, lengths=[-0.02_R8P, 0.07_R8P], pose=pose, &
                   linear_velocity=[0.1_R8P, 0._R8P], angular_velocity=-1._R8P, ilw=.true., device=.true.)
call check('twin posed ellipse: device = host, every field', st == STATUS_OK .and. ch%n > 0_I4P .and. &
           ib_closure_mismatches(ch, cd) == 0_I4P)
nbad = 0_I4P
do g = 1_I4P, ch%n
   do d = 1_I4P, 4_I4P
      if (ch%steps(d, g) > 0_I4P .and. ch%ilw_element(d, g) < 1_I4P) nbad = nbad + 1_I4P
   enddo
enddo
call check('posed ellipse: an ILW crossing in every slot', nbad == 0_I4P)
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P], ch, status=st, ilw=.true.)
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P], cd, status=st, ilw=.true., device=.true.)
call check('twin without curvature and probes: device = host', st == STATUS_OK .and. .not. cd%with_curvature .and. &
           ib_closure_mismatches(ch, cd) == 0_I4P)
! 4. resident lifecycle, refusals, empty curve
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P], cr, status=st, ilw=.true., resident=.true.)
call check('resident: flagged', st == STATUS_OK .and. cr%on_device)
call cr%update_host
call check('resident: update_host gives the host bits', ib_closure_mismatches(ch, cr) == 0_I4P)
cc = cr
call cc%exit_device                                     ! a copy: only the flag drops
call check('resident: a copy is not mapped', .not. cc%on_device .and. cr%on_device)
call fd_closure_2d(plain, x, y, [1_I4P, 1_I4P], cr, status=st)
call check('resident: a new call releases it first', st == STATUS_OK .and. .not. cr%on_device)
call ib_closure_reset(cr)
call plain%exit_device
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P], cd, status=st, device=.true.)
call check('refused: device path of an unmapped curve', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P], cd, status=st, centre=[0._R8P, 0._R8P, 0._R8P])
call check('refused: size(centre) /= 2', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_2d(plain, x, y, [2_I4P, 2_I4P, 2_I4P], cd, status=st)
call check('refused: size(reach) /= 2', st == STATUS_INVALID_INPUT .and. cd%n == 0_I4P)
call fd_closure_2d(empty, x, y, [2_I4P, 2_I4P], cd, status=st, device=.true., lengths=[0.1_R8P], ilw=.true.)
call check('empty curve: no ghost', st == STATUS_OK .and. cd%n == 0_I4P)
call circle%exit_device ; call ellipse%exit_device

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
   type(ib_closure_t), intent(in)  :: c              !< Table.
   integer(I4P),       intent(in)  :: g, d           !< Ghost, slot.
   integer(I4P),       intent(out) :: lo, hi, a, sgn !< Indexes, axis, direction.
   integer(I4P) :: p0, p1                            !< Points.

   a = (d + 1_I4P) / 2_I4P ; sgn = 1_I4P ; if (mod(d, 2_I4P) == 1_I4P) sgn = -1_I4P
   p0 = c%ijk(a, g) + sgn * (c%steps(d, g) - 1_I4P) ; p1 = c%ijk(a, g) + sgn * c%steps(d, g)
   lo = min(p0, p1) ; hi = max(p0, p1)
   endsubroutine slot
endprogram fossil_test_fd_closure_2d
