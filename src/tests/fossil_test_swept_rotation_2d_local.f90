!< FOSSIL test: the 2-D exact rotation by local integration against the per-face reference (issue #104 decision 2).
program fossil_test_swept_rotation_2d_local
!< FOSSIL test: the 2-D exact rotation by local integration against the per-face reference (issue #104 decision 2).
!<
!< `swept_cells_2d` on the exact rotation path integrates every cell the wall can meet per event interval on its own
!< wall segments, gives the faces of no such cell in closed form and sends the faces of cells the piece map cannot
!< integrate to the per-face integration (`faces_rotation_2d`, counted in `sw%fallback_faces`).
!< `swept_cells_2d_reference` runs the per-face integration on every face. Per case:
!< - both `exact`, the same alpha and swept bit for bit;
!< - every face length within 16 eps xm of the reference's and every face outflow within 16 eps xm vmax (absolute),
!<   xm = max |grid| + |t0| + |dt| the coordinate magnitude (the round-off of one crossing), vmax = |dth| xm + |dt| the
!<   largest rigid-field speed: both integrations are converged to 64 eps max(h, 1e-6) of the face side times its scale,
!<   so they agree to round-off (measured: below 1 eps xm in step 0 of the decision);
!< - the residual of every closed cell at round-off: |r| <= 256 eps xm (vmax + h) h (measured up to 43 on the half-cell
!<   disc in the debug build, the reference's own floor: the end-pose areas' round-off), against 1e-3 to 1e-1 of the
!<   cell area for the Gauss path at 4 nodes on such steps (issue #104), so a quadrature-level error fails it; and per
!<   cell within 64 eps xm vmax of the reference's (the same swept area plus four outflows within the bound above);
!< - the outflows through `outflow_u` / `outflow_v` balance the residual: r = swept + net outflow, bit for bit;
!< - with `map=`: the faces, outflows and residual of the call without map, bit for bit.
!< Cases on the local path (no fallback face): a disc (64-gon) turning half a cell about an off-centre point; a screw
!< step (a square turning 0.3 rad while translating); a disc whose body frame lies 20 units from its origin (large |p|:
!< the arc's curvature bound dth^2 |p| of the candidate boxes); a hexagon translating 70 cells while turning (more than
!< the 64 sub-steps of the candidate boxes: each sub-step box spans more than one cell).
!< Cases reaching the fallback (fallback faces > 0, asserted): the corner (1, 1) of the 3 x 3 block grazing a block's
!< face at the turning point of the rotation (an exact tangential contact: the event is not certified, the cell is
!< untrusted), and an open polyline (no cut cell is closed: every face falls back).
use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, swept_cells_2d_reference, &
                   piece_map_2d_t, STATUS_OK
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter   :: PI = acos(-1._R8P)
type(curve_object)     :: c
real(R8P), allocatable :: x(:), y(:), u(:), v(:)
integer(I4P)           :: i
logical                :: ok(6)

! 1. a 64-gon disc (radius 0.2) turning half a cell (h = 1/16) at its far side about (0.3, 0.45), 16 x 16 cells
x = [(real(i, R8P) / 16._R8P, i = 0, 16)] ; y = x
call polygon(c, 64_I4P, [0.55_R8P, 0.5_R8P], 0.2_R8P, 0._R8P)
ok(1) = check('disc, half a cell', about(0.1_R8P, [0.3_R8P, 0.45_R8P], [0._R8P, 0._R8P]), &
              about(0.1_R8P + 0.5_R8P / 16._R8P / 0.45_R8P, [0.3_R8P, 0.45_R8P], [0._R8P, 0._R8P]), .false.)
! 2. a screw step: a square turning 0.3 rad about (0.5, 0.5) while translating (0.07, -0.04)
call c%adopt_loops(u=[0.35_R8P, 0.62_R8P, 0.62_R8P, 0.35_R8P], v=[0.33_R8P, 0.33_R8P, 0.61_R8P, 0.61_R8P], &
                   loop_sizes=[4_I4P])
ok(2) = check('square, screw step', about(0._R8P, [0.5_R8P, 0.5_R8P], [0._R8P, 0._R8P]), &
              about(0.3_R8P, [0.5_R8P, 0.5_R8P], [0.07_R8P, -0.04_R8P]), .false.)
! 3. a disc 20 units from its body origin: world pose = rotation about the body origin (20 units away) plus a shift
!    placing it near (0.5, 0.5); the step turns 0.02 rad (an arc of 0.4 = 6.4 cells), from 0.004 rad: the extreme
!    vertices' turning point (angle 0) stays out of the step
call polygon(c, 48_I4P, [20._R8P, 0._R8P], 0.15_R8P, 0._R8P)
ok(3) = check('disc far from its body origin', far_pose(0.004_R8P), far_pose(0.024_R8P), .false.)
! 4. a hexagon (radius 0.02) translating 0.7 (70 cells of 0.01) while turning 0.5 rad, on 100 x 8 cells
x = [(0.01_R8P * i, i = 0, 100)] ; y = [(0.01_R8P * i, i = 0, 8)]
call polygon(c, 6_I4P, [0._R8P, 0._R8P], 0.02_R8P, 0.1_R8P)
ok(4) = check('hexagon, 70 cells per step', pose_2d(0._R8P, [0.1_R8P, 0.041_R8P]), &
              pose_2d(0.5_R8P, [0.8_R8P, 0.037_R8P]), .false.)
! 5. fallback: the corner (1, 1) of the 3 x 3 block grazing the face of a block rotating about a pivot beyond it, at
!    the turning point (exact tangential contact: the middle cell untrusted)
x = [-1._R8P, 0._R8P, 1._R8P, 2._R8P] ; y = x
call c%adopt_loops(u=[-22._R8P, -2._R8P, -2._R8P, -22._R8P], v=[-10._R8P, -10._R8P, 13._R8P, 13._R8P], loop_sizes=[4_I4P])
ok(5) = check('corner grazing a face (contact)', graze_pose(-1._R8P), graze_pose(1._R8P), .true.)
! 6. fallback: an open polyline (no closed cut cell)
x = [(0.125_R8P * i, i = 0, 8)] ; y = x
u = [0.2_R8P, 0.45_R8P, 0.7_R8P, 0.8_R8P] ; v = [0.3_R8P, 0.55_R8P, 0.5_R8P, 0.75_R8P]
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P], closed=[.false.])
ok(6) = check('open polyline', about(0._R8P, [0.5_R8P, 0.5_R8P], [0._R8P, 0._R8P]), &
              about(0.25_R8P, [0.5_R8P, 0.5_R8P], [0.03_R8P, 0._R8P]), .true., open=.true.)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(label, p0, p1, fallback, open) result(good)
   !< The local path against the reference between p0 and p1 on the grid x, y (see the program notes); `fallback`: the
   !< fallback must be reached (else it must not); `open`: an open curve, whose cut-cell fractions and residuals have no
   !< meaning (no solid side): the residual and `exact` equal to the reference's instead of at round-off and true (the
   !< per-face integration does not converge to round-off on every face of this polyline).
   character(*),          intent(in) :: label    !< Label.
   type(rigid_pose_2d_t), intent(in) :: p0, p1   !< Poses.
   logical,               intent(in) :: fallback !< Fallback expected.
   logical, optional,     intent(in) :: open     !< An open curve.
   logical                           :: good     !< Result.
   type(swept_cells_2d_t)            :: sl, sr, sm !< Local, reference, local with the map.
   type(piece_map_2d_t)              :: mp       !< Map.
   real(R8P), allocatable            :: lu(:,:), lv(:,:), ru(:,:), rv(:,:), mu(:,:), mv(:,:) !< Outflows.
   real(R8P)                         :: xm, vmax, dth, h, ef, eo, er, rb !< Scales, errors.
   integer(I4P)                      :: nx, ny, s1, s2, s3, i_, j_, nbad !< Work.
   logical                           :: closed   !< Closed cell.
   logical                           :: ex       !< Both exact (an open curve: the same exact).

   nx = size(x) - 1 ; ny = size(y) - 1
   allocate(lu(nx + 1, ny), lv(nx, ny + 1), ru(nx + 1, ny), rv(nx, ny + 1), mu(nx + 1, ny), mv(nx, ny + 1))
   call swept_cells_2d(c, p0, p1, x, y, sl, status=s1, outflow_u=lu, outflow_v=lv)
   call swept_cells_2d_reference(c, p0, p1, x, y, sr, status=s2, outflow_u=ru, outflow_v=rv)
   call swept_cells_2d(c, p0, p1, x, y, sm, status=s3, map=mp, outflow_u=mu, outflow_v=mv)
   dth = atan2(p1%r(2, 1), p1%r(1, 1)) - atan2(p0%r(2, 1), p0%r(1, 1)) ; dth = atan2(sin(dth), cos(dth))
   xm = max(maxval(abs(x)), maxval(abs(y))) + norm2(p0%t) + norm2(p1%t - p0%t)
   vmax = abs(dth) * xm + norm2(p1%t - p0%t)
   h = min(minval(x(2:) - x(:nx)), minval(y(2:) - y(:ny)))
   ef = max(maxval(abs(sl%face_u - sr%face_u)), maxval(abs(sl%face_v - sr%face_v))) / (epsilon(1._R8P) * xm)
   eo = max(maxval(abs(lu - ru)), maxval(abs(lv - rv))) / (epsilon(1._R8P) * xm * vmax)
   er = 0._R8P ; rb = 0._R8P ; nbad = 0_I4P
   do j_ = 1, ny
      do i_ = 1, nx
         closed = sl%alpha0(i_, j_) >= 0._R8P .and. sl%alpha1(i_, j_) >= 0._R8P
         if (closed) er = max(er, abs(sl%residual(i_, j_)))
         if (closed) rb = max(rb, abs(sr%residual(i_, j_)))
         if (closed .and. abs(sl%residual(i_, j_) - sr%residual(i_, j_)) > 64._R8P * epsilon(1._R8P) * xm * vmax) &
            nbad = nbad + 1_I4P
         ! the residual from the outflows, as the library forms it
         if (sl%residual(i_, j_) /= sl%swept(i_, j_) + (lu(i_ + 1, j_) - lu(i_, j_) + lv(i_, j_ + 1) - lv(i_, j_))) &
            nbad = nbad + 1_I4P
      enddo
   enddo
   ex = sl%exact .and. sr%exact
   er = er / (epsilon(1._R8P) * xm * (vmax + h) * h) ; rb = rb / (epsilon(1._R8P) * xm * (vmax + h) * h)
   if (present(open)) then
      if (open) then
         if (any(sl%residual /= sr%residual) .or. (sl%exact .neqv. sr%exact)) nbad = nbad + 1_I4P
         er = 0._R8P ; rb = 0._R8P ; ex = .true.
      endif
   endif
   if (any(sl%alpha0 /= sr%alpha0) .or. any(sl%alpha1 /= sr%alpha1) .or. any(sl%swept /= sr%swept)) nbad = nbad + 1_I4P
   if (any(sm%face_u /= sl%face_u) .or. any(sm%face_v /= sl%face_v) .or. any(sm%residual /= sl%residual) .or. &
       any(mu /= lu) .or. any(mv /= lv)) nbad = nbad + 1_I4P
   good = s1 == STATUS_OK .and. s2 == STATUS_OK .and. s3 == STATUS_OK .and. ex .and. &
          nbad == 0_I4P .and. ef <= 16._R8P .and. eo <= 16._R8P .and. er <= 256._R8P
   if (fallback) then
      good = good .and. sl%fallback_faces > 0_I4P
   else
      good = good .and. sl%fallback_faces == 0_I4P
   endif
   print '(A,A,ES9.2,A,ES9.2,A,ES9.2,A,ES9.2,A,I0,A,I0,A,I0,A,I0,A,L1,A,L1,A,L1)', label, ': |dface| / eps xm ', ef, &
         ', |doutflow| / eps xm vmax ', eo, ', max |residual| / eps xm (vmax + h) h ', er, ' (reference ', rb, &
         '), fallback faces ', sl%fallback_faces, ', evaluations ', sl%evaluations, ' (reference ', sr%evaluations, &
         '), bad ', nbad, ', exact ', sl%exact, ' (reference ', sr%exact, '), ok ', good
   endfunction check

   subroutine polygon(cv, n, cen, r, ph)
   !< A regular n-gon of circumradius r about cen, counter-clockwise, its first vertex at angle ph.
   type(curve_object), intent(inout) :: cv       !< Curve.
   integer(I4P),       intent(in)    :: n        !< Vertices.
   real(R8P),          intent(in)    :: cen(2), r, ph !< Centre, radius, phase.
   real(R8P)                         :: pu(n), pv(n) !< Vertices.
   integer(I4P)                      :: k        !< Counter.

   do k = 1, n
      pu(k) = cen(1) + r * cos(ph + 2._R8P * PI * (k - 1) / n) ; pv(k) = cen(2) + r * sin(ph + 2._R8P * PI * (k - 1) / n)
   enddo
   call cv%adopt_loops(u=pu, v=pv, loop_sizes=[n])
   endsubroutine polygon

   pure function about(angle, centre, shift) result(p)
   !< Rotation by angle about centre, then a translation.
   real(R8P), intent(in) :: angle, centre(2), shift(2) !< Motion.
   type(rigid_pose_2d_t) :: p                          !< Pose.

   p = pose_2d(angle)
   p%t = centre - matmul(p%r, centre) + shift
   endfunction about

   pure function far_pose(angle) result(p)
   !< The body origin 20 units from the disc (body frame (20, 0)): rotation by angle about it, the disc's centre placed
   !< at (0.5, 0.5) for angle 0.
   real(R8P), intent(in) :: angle !< Angle.
   type(rigid_pose_2d_t) :: p     !< Pose.

   p = pose_2d(angle, [-19.5_R8P, 0.5_R8P])
   endfunction far_pose

   function graze_pose(sgn) result(p)
   !< The body origin on the pivot (1, 1) + 2 (1, 1) / sqrt(2), its +u axis along the diagonal turned by sgn 0.3: the
   !< block's face u = -2 passes the corner (1, 1) at the middle of the step, touching it.
   real(R8P), intent(in) :: sgn !< -1 (pose0) or 1 (pose1).
   type(rigid_pose_2d_t) :: p   !< Pose.

   p = pose_2d(PI / 4._R8P + sgn * 0.3_R8P, [1._R8P, 1._R8P] + sqrt(2._R8P) * [1._R8P, 1._R8P])
   endfunction graze_pose
endprogram fossil_test_swept_rotation_2d_local
