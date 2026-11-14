!< FOSSIL test: the 2-D piece map of `swept_cells_2d_device` equals that of `swept_cells_2d` bit for bit (issue #105 step 2).
program fossil_test_piece_map_2d_device_twin
!< FOSSIL test: the 2-D piece map of `swept_cells_2d_device` equals that of `swept_cells_2d` bit for bit (issue #105 step 2).
!<
!< Built without OpenACC (this CI build) the device loops of the map (`piece_map_2d_build(device=)`: the evaluations'
!< sizes, count, evaluation and fill, the tracks per cell, the nodes of every measure round with their sums, the round's
!< test per pair, the assembly) run on the host: the test catches the driver bugs of the device twin (the flat arrays
!< handed to the launch sets, the padded slices, the launch sets themselves, the per-cell sums), not device codegen (the
!< curve bench does that on the GPU, `piece map 2-D device vs serial mismatches`). Every map compared by its bits
!< (`piece_map_2d_mismatches`), the swept geometry and the status too:
!< - exact translation: an ellipse on 24^2 cells; the NACA 0012's thin aft region moving 2 cells across (multi-track
!<   cells measured per event interval); a plate 2e-11 wide crossing a cell (states between events 1e-11 apart);
!< - exact rotation: the ellipse, a lattice of squares and bars on 12^2 cells (many multi-track cells, doubling rounds),
!<   a plate crossing a cell about a far pivot;
!< - the Gauss path (nq = 4) of the lattice rotating, and with `refine` (every evaluated cell measured);
!< - the lattice cases again in launch sets of 4 KB of workspace (`chunk_bytes`: many sets per pass);
!< - an empty curve, and the refusal of a map with node poses: the same statuses;
!< - negative control: a pose one ulp off gives mismatches.
use fossil,                   only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, &
                                     swept_cells_2d_device, piece_map_2d_t
use fossil_piece_map_compare, only : piece_map_2d_mismatches
use fossil_swept_compare,     only : swept_cells_2d_mismatches
use penf,                     only : I4P, I8P, R8P

implicit none

real(R8P), parameter   :: PI = acos(-1._R8P) !< Pi.
type(curve_object)     :: ce, cl, cn, cp, em !< Ellipse, lattice, NACA 0012, plates, empty curve.
type(piece_map_2d_t)   :: a, b               !< Maps (negative control).
type(swept_cells_2d_t) :: sa, sb             !< Swept geometry (negative control).
real(R8P), allocatable :: u(:), v(:), xg(:), xl(:) !< Ellipse, grids.
real(R8P)              :: x3(4)              !< The 3 x 3 block [-1, 2]^2.
integer(I4P)           :: k, m, total        !< Work.
logical                :: ok                 !< All cases.

ok = .true. ; total = 0_I4P
allocate(u(64), v(64))
do k = 1, 64
   u(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 64) ; v(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 64)
enddo
call ce%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
call lattice(cl) ; call naca(cn)
xg = [(0.0031_R8P + 0.9938_R8P * k / 24._R8P, k = 0, 24)]
xl = [(-0.25_R8P + 0.125_R8P * k, k = 0, 12)]
x3 = [-1._R8P, 0._R8P, 1._R8P, 2._R8P]
! exact translation
call check('ellipse, translation', twin(ce, pose_2d(0._R8P), pose_2d(0._R8P, [1.7_R8P / 24._R8P, 0.6_R8P / 24._R8P]), xg, &
                                        .true.))
call check('NACA 0012 aft region, 2 cells across', twin(cn, pose_2d(0._R8P), &
                                                        pose_2d(0._R8P, [0.3_R8P, 2._R8P] * (xg(2) - xg(1))), xg, .true.))
call rect(cp, [-0.55_R8P, -0.55_R8P + 2e-11_R8P, -3._R8P, 4._R8P])
call check('plate 2e-11 wide crossing, translation', twin(cp, pose_2d(0._R8P), pose_2d(0._R8P, [2._R8P, 0._R8P]), x3, &
                                                          .true.))
call check('plate 2e-11 wide crossing, Gauss', twin(cp, pose_2d(0._R8P), pose_2d(0._R8P, [2._R8P, 0._R8P]), x3, .false.))
! exact rotation
call check('ellipse, rotation', twin(ce, about(0._R8P), about(PI / 18._R8P), xg, .true.))
call check('lattice, rotation', twin(cl, about(0._R8P), about(0.35_R8P), xl, .true.))
call rect(cp, [-0.05_R8P, 0.05_R8P, 3._R8P, 8._R8P])
call check('plate crossing, rotation about a far pivot', twin(cp, pose_2d(0.2_R8P, [0.5_R8P, -5._R8P]), &
                                                              pose_2d(-0.2_R8P, [0.5_R8P, -5._R8P]), x3, .true.))
! Gauss path, refine
call check('lattice, rotation, Gauss', twin(cl, about(0._R8P), about(0.35_R8P), xl, .false.))
call check('lattice, rotation, Gauss + refine', twin(cl, about(0._R8P), about(0.35_R8P), xl, .false., refine=.true.))
! launch sets of 4 KB
call check('lattice, rotation, 4 KB sets', twin(cl, about(0._R8P), about(0.35_R8P), xl, .true., chunk=4096_I8P))
call check('lattice, Gauss + refine, 4 KB sets', twin(cl, about(0._R8P), about(0.35_R8P), xl, .false., refine=.true., &
                                                      chunk=4096_I8P))
call check('lattice, translation, 4 KB sets', twin(cl, pose_2d(0._R8P), pose_2d(0._R8P, [0.31_R8P, 0.17_R8P]), xl, .true., &
                                                   chunk=4096_I8P))
! empty curve, refusal
call check('empty curve, translation', twin(em, pose_2d(0._R8P), pose_2d(0._R8P, [0.1_R8P, 0._R8P]), x3, .true.))
call swept_cells_2d(ce, pose_2d(0._R8P), pose_2d(0._R8P), xg, xg, sa, nq=1_I4P, poses=[pose_2d(0._R8P)], status=k, map=a)
call swept_cells_2d_device(ce, pose_2d(0._R8P), pose_2d(0._R8P), xg, xg, sb, nq=1_I4P, poses=[pose_2d(0._R8P)], status=m, &
                           map=b)
call check('refusal: map with node poses', abs(k - m) + merge(0_I4P, 1_I4P, k /= 0_I4P))
! negative control
call swept_cells_2d(ce, about(0._R8P), about(PI / 18._R8P), xg, xg, sa, map=a)
call swept_cells_2d_device(ce, about(0._R8P), about(nearest(PI / 18._R8P, 1._R8P)), xg, xg, sb, map=b)
m = piece_map_2d_mismatches(a, b)
print '(A,I0)', 'negative control (one ulp): mismatches ', m
ok = ok .and. m > 0_I4P
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   function twin(c, p0, p1, x, exact, refine, chunk) result(n)
   !< Mismatches of the device twin against the host call (same arguments, nq = 4): map, swept geometry, status.
   type(curve_object),    intent(in)           :: c      !< Curve.
   type(rigid_pose_2d_t), intent(in)           :: p0, p1 !< End poses.
   real(R8P),             intent(in)           :: x(:)   !< Cell edges (both directions).
   logical,               intent(in)           :: exact  !< Exact paths.
   logical,               intent(in), optional :: refine !< Gauss path refined per event interval.
   integer(I8P),          intent(in), optional :: chunk  !< Workspace bytes per launch set.
   integer(I4P)                                :: n      !< Mismatches.
   type(piece_map_2d_t)                        :: mh, md !< Maps.
   type(swept_cells_2d_t)                      :: wh, wd !< Swept geometry.
   integer(I4P)                                :: sh, sd !< Statuses.

   call swept_cells_2d(c, p0, p1, x, x, wh, nq=4_I4P, exact=exact, status=sh, map=mh, refine=refine)
   call swept_cells_2d_device(c, p0, p1, x, x, wd, nq=4_I4P, exact=exact, status=sd, map=md, refine=refine, chunk_bytes=chunk)
   n = piece_map_2d_mismatches(mh, md, verbose=.true.) + swept_cells_2d_mismatches(wh, wd) + abs(sh - sd)
   if (sh == 0_I4P .and. mh%nt > 0_I4P) print '(A,I0,A,I0,A,I0,A,I0)', '   tracks ', mh%nt, ', evaluated cells ', &
      mh%evaluated_cells, ', multi-track cells ', mh%multi_track_cells, ', narrow cells ', mh%narrow_cells
   endfunction twin

   subroutine check(label, n)
   !< One case.
   character(*), intent(in) :: label !< Case.
   integer(I4P), intent(in) :: n     !< Mismatches.

   print '(A,I0)', label//': mismatches ', n
   total = total + n ; ok = ok .and. n == 0_I4P
   endsubroutine check

   function about(angle) result(p)
   !< A rotation by angle about (0.5, 0.5).
   real(R8P), intent(in) :: angle !< Angle.
   type(rigid_pose_2d_t) :: p     !< Pose.

   p = pose_2d(angle)
   p%t = [0.5_R8P, 0.5_R8P] - matmul(p%r, [0.5_R8P, 0.5_R8P])
   endfunction about

   subroutine rect(cv, b)
   !< The rectangle [b(1), b(2)] x [b(3), b(4)], counter-clockwise.
   type(curve_object), intent(inout) :: cv   !< Curve.
   real(R8P),          intent(in)    :: b(4) !< Bounds.

   call cv%adopt_loops(u=[b(1), b(2), b(2), b(1)], v=[b(3), b(3), b(4), b(4)], loop_sizes=[4_I4P])
   endsubroutine rect

   subroutine lattice(cv)
   !< 4 x 4 squares of side 0.12 at pitch 0.25 from (0.065, 0.065), and three bars 0.02 thick (two horizontal, one
   !< vertical) between them; all counter-clockwise (the lattice of `fossil_test_piece_map_2d`).
   type(curve_object), intent(inout) :: cv          !< Curve.
   real(R8P)                         :: uu(76), vv(76), o(2) !< Points.
   integer(I4P)                      :: i, j, n     !< Counters.

   n = 0_I4P
   do i = 0, 3
      do j = 0, 3
         o = [0.065_R8P + 0.25_R8P * j, 0.065_R8P + 0.25_R8P * i]
         uu(n + 1:n + 4) = o(1) + [0._R8P, 0.12_R8P, 0.12_R8P, 0._R8P]
         vv(n + 1:n + 4) = o(2) + [0._R8P, 0._R8P, 0.12_R8P, 0.12_R8P]
         n = n + 4_I4P
      enddo
   enddo
   uu(65:68) = [0.02_R8P, 0.98_R8P, 0.98_R8P, 0.02_R8P] ; vv(65:68) = [0.24_R8P, 0.24_R8P, 0.26_R8P, 0.26_R8P]
   uu(69:72) = [0.02_R8P, 0.98_R8P, 0.98_R8P, 0.02_R8P] ; vv(69:72) = [0.74_R8P, 0.74_R8P, 0.76_R8P, 0.76_R8P]
   uu(73:76) = [0.49_R8P, 0.51_R8P, 0.51_R8P, 0.49_R8P] ; vv(73:76) = [0.27_R8P, 0.27_R8P, 0.73_R8P, 0.73_R8P]
   call cv%adopt_loops(u=uu, v=vv, loop_sizes=[(4_I4P, i = 1, 19)])
   endsubroutine lattice

   subroutine naca(cv)
   !< NACA 0012 of chord 0.7 from (0.15, 0.5), 2 x 40 cosine-spaced points, closed trailing edge.
   type(curve_object), intent(inout) :: cv !< Curve.
   real(R8P)                         :: uu(80), vv(80), xc !< Points.
   integer(I4P)                      :: i  !< Counter.

   do i = 0, 39 ! upper side from the trailing edge, then lower side from the leading edge
      xc = 0.5_R8P * (1._R8P + cos(PI * i / 40._R8P))
      uu(i + 1) = xc ; vv(i + 1) = half(xc)
      xc = 0.5_R8P * (1._R8P - cos(PI * i / 40._R8P))
      uu(41 + i) = xc ; vv(41 + i) = -half(xc)
   enddo
   uu = 0.15_R8P + 0.7_R8P * uu ; vv = 0.5_R8P + 0.7_R8P * vv
   call cv%adopt_loops(u=uu, v=vv, loop_sizes=[80_I4P])
   endsubroutine naca

   pure function half(xc) result(yc)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xc !< Chord position.
   real(R8P)             :: yc !< Half thickness.

   yc = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1036_R8P * xc**4)
   endfunction half
endprogram fossil_test_piece_map_2d_device_twin
