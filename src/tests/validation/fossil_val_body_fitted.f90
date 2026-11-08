!< FOSSIL validation: body-fitted grid support at size (issue #74, step 7).
!<
!< 1. The NACA0012 wing (src/tests/naca0012-ascii.stl: chord x in [-0.25, 0.75], span z in [0, 1], blunt trailing
!<    edge): a 41 x 21 structured grid 0.005 above its upper surface (cosine chordwise) projected as a block face with
!<    tolerance 0.01: its trailing-edge row lands on the upper trailing-edge line (x = 0.75, y = 0.00126) exactly, its
!<    root and tip rows on z = 0 and z = 1, the root and tip trailing-edge corners on the wing's corners; projected
!<    again every vertex moves by less than 1e-12 (measured 3.5e-13: feet on edges of the long spanwise triangles).
!<    Negative control: free projection moves trailing-edge vertices off that line.
!< 2. An O-grid extruded from the NACA0012 section (129 points, 40 layers from h1 = nu / u_tau, ratio 1.2) and its
!<    report: every cell certified valid, y1+ within 2 % of 1 away from the trailing edge (reported in full).
!< 3. The bunny: a 101 x 101 grid on a plane through it projected free (tolerance 0): time, largest move, idempotence.
!< Exit status 0 if every check holds.
program fossil_val_body_fitted

use fossil, only : surface_stl_object, curve_object, projection_target_t, curvilinear_block_t, curvilinear_block_2d_t, &
                   topo_record_t, TOPO_BOUNDARY, project_block_face, project_points, PROJ_FREE, PROJ_CURVE, PROJ_CORNER, &
                   hyperbolic_grid_2d, geometric_spacing, wall_cutter_2d, block_report_t, block_report_2d
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter         :: PI = acos(-1._R8P), NU = 1.5e-5_R8P, UTAU = 0.4_R8P
type(surface_stl_object)     :: wing, bunny
type(projection_target_t)    :: t
type(curvilinear_block_t)    :: b, bf
type(curvilinear_block_2d_t) :: g(1)
type(curve_object)           :: wall
type(topo_record_t)          :: rec
type(block_report_t)         :: rep
real(R8P), allocatable       :: x(:,:,:), y(:,:,:), z(:,:,:), moved(:,:), fx(:), fy(:), fz(:)
integer(I4P), allocatable    :: kinds(:,:)
real(R8P)                    :: xc, yt, u(129), v(129), t0, t1, w, dev
type(vector_R8P)             :: lo, hi
integer(I4P)                 :: i, j, n, invalid, st, noff
logical                      :: all_ok, ok

all_ok = .true.

! 1. NACA0012 wing
call wing%load_from_file(file_name='src/tests/naca0012-ascii.stl', guess_format=.true.)
call wing%projection_target(t)
allocate(x(41, 21, 2), y(41, 21, 2), z(41, 21, 2))
do j = 1, 21
   do i = 1, 41
      xc = 0.5_R8P * (1._R8P - cos(PI * (i - 1) / 40._R8P))
      yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
      x(i, j, :) = xc - 0.25_R8P ; y(i, j, :) = [yt + 0.005_R8P, yt + 0.1_R8P] ; z(i, j, :) = (j - 1) / 20._R8P
   enddo
enddo
call b%initialize(x, y, z) ; bf = b
call project_block_face(b, 5_I4P, t, 0.01_R8P, moved=moved, kinds=kinds, status=st)
ok = st == 0
ok = ok .and. all(b%x(41, :, 1) == 0.75_R8P) .and. all(abs(b%y(41, :, 1) - 0.00126_R8P) <= 1e-15_R8P)
ok = ok .and. all(b%z(:, 1, 1) == 0._R8P) .and. all(b%z(:, 21, 1) == 1._R8P)
ok = ok .and. kinds(41, 1) == PROJ_CORNER .and. kinds(41, 21) == PROJ_CORNER
print '(A,I0,A,I0,A,ES12.4)', 'wing: corners ', count(kinds == PROJ_CORNER), ', on curves ', count(kinds == PROJ_CURVE), &
                               ', largest move ', maxval(moved)
call check(ok, 'wing: trailing edge on its line, root and tip rows on z = 0, 1, corners on corners', maxval(moved))
allocate(fx(41 * 21), fy(41 * 21), fz(41 * 21))
call project_points(t, pack(b%x(:, :, 1), .true.), pack(b%y(:, :, 1), .true.), pack(b%z(:, :, 1), .true.), &
                    [(PROJ_FREE, i = 1, 41 * 21)], [(0_I4P, i = 1, 41 * 21)], fx, fy, fz)
w = maxval(abs([fx - pack(b%x(:, :, 1), .true.), fy - pack(b%y(:, :, 1), .true.), fz - pack(b%z(:, :, 1), .true.)]))
call check(w <= 1e-12_R8P, 'wing: the projected grid is on the surface (re-projection moves by)', w)
call project_points(t, bf%x(41, :, 1), bf%y(41, :, 1), bf%z(41, :, 1), [(PROJ_FREE, i = 1, 21)], [(0_I4P, i = 1, 21)], &
                    fx(1:21), fy(1:21), fz(1:21))
noff = count(abs(fy(1:21) - 0.00126_R8P) > 1e-6_R8P .or. fx(1:21) /= 0.75_R8P)
call check(noff > 0, 'negative control: free projection moves trailing-edge vertices off the line (count)', real(noff, R8P))

! 2. O-grid and report
do i = 1, 129
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / 129._R8P))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= 65)
enddo
call hyperbolic_grid_2d(u, v, geometric_spacing(NU / UTAU, 1.2_R8P, 40_I4P), g(1), invalid=invalid, status=st)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [g(1)%ni, 1, 1] ; rec%name = 'wall'
call g(1)%add_record(rec)
call wall_cutter_2d(g, wall)
call block_report_2d(g(1), wall, NU, rep, u_tau=UTAU)
dev = maxval(abs(rep%y1p(4:g(1)%ni - 3) - 1._R8P))
print '(A,3F9.5,A,F8.4,A,F8.4)', 'O-grid y1+ min mean max ', rep%y1p_min, rep%y1p_mean, rep%y1p_max, &
                                 ', growth max ', rep%growth_max, ', skewness max ', rep%skew_max
call check(st == 0 .and. invalid == 0 .and. dev <= 0.02_R8P, 'O-grid: every cell valid, y1+ within 2 % of 1 (max dev)', dev)

! 3. bunny
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call cpu_time(t0) ; call bunny%projection_target(t) ; call cpu_time(t1)
print '(A,I0,A,F8.3,A)', 'bunny: ', bunny%get_facets_number(), ' facets, target built in ', t1 - t0, ' s'
lo = bunny%get_bmin() ; hi = bunny%get_bmax()
deallocate(x, y, z) ; allocate(x(101, 101, 2), y(101, 101, 2), z(101, 101, 2))
do j = 1, 101
   do i = 1, 101
      x(i, j, :) = lo%x + (hi%x - lo%x) * (i - 1) / 100 ; y(i, j, :) = lo%y + (hi%y - lo%y) * (j - 1) / 100
      z(i, j, :) = [0.5_R8P * (lo%z + hi%z), hi%z]
   enddo
enddo
call b%initialize(x, y, z)
call cpu_time(t0) ; call project_block_face(b, 5_I4P, t, 0._R8P, moved=moved, status=st) ; call cpu_time(t1)
n = 101 * 101 ; deallocate(fx, fy, fz) ; allocate(fx(n), fy(n), fz(n))
call project_points(t, pack(b%x(:, :, 1), .true.), pack(b%y(:, :, 1), .true.), pack(b%z(:, :, 1), .true.), &
                    [(PROJ_FREE, i = 1, n)], [(0_I4P, i = 1, n)], fx, fy, fz)
w = maxval(abs([fx - pack(b%x(:, :, 1), .true.), fy - pack(b%y(:, :, 1), .true.), fz - pack(b%z(:, :, 1), .true.)]))
print '(A,F8.3,A,ES12.4)', 'bunny: 10201 vertices projected in ', t1 - t0, ' s, largest move ', maxval(moved)
call check(st == 0 .and. w <= 1e-12_R8P, 'bunny: on the surface (re-projection moves by)', w)

print '(A,L1)', 'Are all tests passed? ', all_ok
if (.not. all_ok) error stop 1

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES23.15)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check
endprogram fossil_val_body_fitted
