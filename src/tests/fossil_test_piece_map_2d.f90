!< FOSSIL test: the piece map of a time step in 2-D, the tracks of every cell (issue #102 step 1).
program fossil_test_piece_map_2d
!< FOSSIL test: the piece map of a time step in 2-D, the tracks of every cell (issue #102 step 1).
!<
!< Hand-built cases on the 3 x 3 block [-1, 2]^2, each a curve and a motion whose tracks in the middle cell [0, 1]^2
!< are known (the end-pose overlap, which links an old and a new piece iff they overlap, is given for contrast):
!< - a thin plate crossing the whole cell within the step (translation; rotation about a far pivot; the Gauss path
!<   of both): one dead and one fresh track (overlap: continuing);
!< - a narrow gap moving more than its width inside the cell: continuing (overlap: dead + fresh);
!< - a plate entering from the top and crossing the cell: split; leaving: merge;
!< - a block uncovering the cell: fresh (cleared); covering it: dead (covered);
!< - a ring whose hole is a fluid island moving inside the cell: two continuing tracks;
!< - a strip sweeping over a cell corner and stopping past it: the main piece continues and the corner piece is fresh
!<   (overlap: a split);
!< - a fluid wedge whose apex enters the cell through a face: continuing. The two segments at the apex cross the face
!<   line at the event and persist through it; a rule calling them volatile would declare the piece dead + fresh.
!< Every case also passes the map checks (`check_map`): every piece at pose0 and at pose1 in exactly one track of its
!< cell, numbered as the cut cells of a copy moved to the pose (and the track volumes equal to their pieces' areas bit
!< for bit), per-cell sums over the tracks equal to the swept geometry (faces, wall, residual, swept = V1 - V0), every
!< per-track residual at round-off on the exact paths, residual = V1 - V0 + sum of the track's face outflows, and each face's
!< outflow summed over the tracks equal from its two cells (and to d times its aperture on a translation). Then the same checks on an ellipse and a NACA 0012 at a modest
!< grid, translating and rotating, and on a cell of the NACA's thin aft region crossed within the step (one piece at
!< each end, dead + fresh).
use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, piece_map_2d_t, &
                   cut_cells_2d_t, STATUS_OK, STATUS_INVALID_INPUT, TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, &
                   TRACK_SPLIT, TRACK_MERGE, TRACK_UNTRUSTED
use penf,   only : I4P, R8P

implicit none

type(curve_object)     :: c
type(swept_cells_2d_t) :: sw
type(piece_map_2d_t)   :: map
real(R8P)              :: x(4)
real(R8P), allocatable :: u(:), v(:), xg(:)
real(R8P), parameter   :: PI = acos(-1._R8P)
integer(I4P)           :: i, st
logical                :: ok(20)

ok = .true.
x = [-1._R8P, 0._R8P, 1._R8P, 2._R8P]
! 1-2. thin plate crossing the cell within the step: dead + fresh (exact translation, then the Gauss path)
call rect(c, [-0.55_R8P, -0.45_R8P, -3._R8P, 4._R8P])
ok(1) = case_2d('plate crossing, translation', pose_2d(0._R8P), pose_2d(0._R8P, [2._R8P, 0._R8P]), .true., &
                [TRACK_DEAD, TRACK_FRESH])
ok(2) = case_2d('plate crossing, translation, Gauss', pose_2d(0._R8P), pose_2d(0._R8P, [2._R8P, 0._R8P]), .false., &
                [TRACK_DEAD, TRACK_FRESH])
! 3-4. the same plate rotating about a far pivot through the cell: exact rotation, Gauss
call rect(c, [-0.05_R8P, 0.05_R8P, 3._R8P, 8._R8P])
ok(3) = case_2d('plate crossing, rotation', pose_2d(0.2_R8P, [0.5_R8P, -5._R8P]), pose_2d(-0.2_R8P, [0.5_R8P, -5._R8P]), &
                .true., [TRACK_DEAD, TRACK_FRESH])
ok(4) = case_2d('plate crossing, rotation, Gauss', pose_2d(0.2_R8P, [0.5_R8P, -5._R8P]), &
                pose_2d(-0.2_R8P, [0.5_R8P, -5._R8P]), .false., [TRACK_DEAD, TRACK_FRESH])
! 5. a gap of width 0.05 moving by 0.2 inside the cell: continuing
call c%adopt_loops(u=[-3._R8P, 4._R8P, 4._R8P, -3._R8P, -3._R8P, 4._R8P, 4._R8P, -3._R8P], &
                   v=[-3._R8P, -3._R8P, 0.40_R8P, 0.40_R8P, 0.45_R8P, 0.45_R8P, 4._R8P, 4._R8P], loop_sizes=[4_I4P, 4_I4P])
ok(5) = case_2d('narrow gap moving past its width', pose_2d(0._R8P), pose_2d(0._R8P, [0._R8P, 0.2_R8P]), .true., &
                [TRACK_CONTINUING])
! 6-7. a plate entering from the top and crossing the cell: split; the reverse motion: merge
call rect(c, [0.45_R8P, 0.55_R8P, 1.2_R8P, 5._R8P])
ok(6) = case_2d('plate entering: split', pose_2d(0._R8P), pose_2d(0._R8P, [0._R8P, -1.6_R8P]), .true., [TRACK_SPLIT])
ok(7) = case_2d('plate leaving: merge', pose_2d(0._R8P, [0._R8P, -1.6_R8P]), pose_2d(0._R8P), .true., [TRACK_MERGE])
! 8-9. a block uncovering the cell: cleared (fresh); covering it: covered (dead)
call rect(c, [-3._R8P, 1.2_R8P, -3._R8P, 4._R8P])
ok(8) = case_2d('block uncovering: cleared', pose_2d(0._R8P), pose_2d(0._R8P, [-1.5_R8P, 0._R8P]), .true., [TRACK_FRESH])
ok(9) = case_2d('block covering: covered', pose_2d(0._R8P, [-1.5_R8P, 0._R8P]), pose_2d(0._R8P), .true., [TRACK_DEAD])
ok(9) = ok(9) .and. sw%covered(2, 2) .and. map%tracks(2, 2) == 1_I4P
! 10. a ring inside the cell, its hole a fluid island: two continuing tracks
call c%adopt_loops(u=[0.2_R8P, 0.8_R8P, 0.8_R8P, 0.2_R8P, 0.35_R8P, 0.35_R8P, 0.65_R8P, 0.65_R8P], &
                   v=[0.2_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, 0.35_R8P, 0.65_R8P, 0.65_R8P, 0.35_R8P], loop_sizes=[4_I4P, 4_I4P])
ok(10) = case_2d('fluid island', pose_2d(0._R8P), pose_2d(0._R8P, [0.05_R8P, 0.03_R8P]), .true., &
                 [TRACK_CONTINUING, TRACK_CONTINUING])
! 11. a strip at 45 degrees sweeping over the corner (1, 1) and stopping past it: continuing + fresh
call c%adopt_loops(u=[0.8_R8P + 3._R8P, 0.85_R8P + 3._R8P, 0.85_R8P - 3._R8P, 0.8_R8P - 3._R8P], &
                   v=[0.8_R8P - 3._R8P, 0.85_R8P - 3._R8P, 0.85_R8P + 3._R8P, 0.8_R8P + 3._R8P], loop_sizes=[4_I4P])
ok(11) = case_2d('strip over a corner', pose_2d(0._R8P, [0.3_R8P, 0.3_R8P]), pose_2d(0._R8P), .true., &
                 [TRACK_CONTINUING, TRACK_FRESH])
! 12. a fluid wedge whose apex enters the cell through its left face: continuing
call c%adopt_loops(u=[-2._R8P, 3._R8P, 3._R8P, -0.1_R8P, 3._R8P, 3._R8P, -2._R8P], &
                   v=[-2._R8P, -2._R8P, 0.19_R8P, 0.5_R8P, 0.81_R8P, 3._R8P, 3._R8P], loop_sizes=[7_I4P])
ok(12) = case_2d('wedge apex entering', pose_2d(0._R8P), pose_2d(0._R8P, [0.3_R8P, 0._R8P]), .true., [TRACK_CONTINUING])
! 13. refusal: a map needs the screw motion, not node poses
call swept_cells_2d(c, pose_2d(0._R8P), pose_2d(0._R8P), x, x, sw, nq=1_I4P, poses=[pose_2d(0._R8P)], status=st, map=map)
ok(13) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'refusal (map with node poses): ', ok(13)
! 14-17. ellipse and NACA 0012, translating and rotating, every cell checked
allocate(u(64), v(64))
do i = 1, 64
   u(i) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (i - 1) / 64) ; v(i) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (i - 1) / 64)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
xg = [(0.0031_R8P + 0.9938_R8P * i / 24._R8P, i = 0, 24)]
ok(14) = body_2d('ellipse, translation', pose_2d(0._R8P), pose_2d(0._R8P, [1.7_R8P / 24._R8P, 0.6_R8P / 24._R8P]))
ok(15) = body_2d('ellipse, rotation', about(0._R8P), about(PI / 18._R8P))
call naca(c)
ok(16) = body_2d('NACA 0012, translation', pose_2d(0._R8P), pose_2d(0._R8P, [-1.9_R8P / 24._R8P, 0.7_R8P / 24._R8P]))
ok(17) = body_2d('NACA 0012, rotation', about(0._R8P), about(0.12_R8P))
! 20. the NACA 0012's thin aft region moving across by 2 cells: cell (20, 13) has one piece at each end and yet a dead
! and a fresh track (a cell with one piece at each end is trivial only without events; overlap: continuing)
ok(20) = body_2d('NACA 0012 aft region, 2 cells across', pose_2d(0._R8P), &
                 pose_2d(0._R8P, [0.3_R8P, 2._R8P] * (xg(2) - xg(1))))
i = map%first_track(20, 13)
ok(20) = ok(20) .and. map%tracks(20, 13) == 2_I4P .and. map%status(i) == TRACK_DEAD .and. &
         map%status(i + 1) == TRACK_FRESH .and. map%old_start(i + 2) - map%old_start(i) == 1_I4P .and. &
         map%new_start(i + 2) - map%new_start(i) == 1_I4P
print '(A,L1)', '  NACA 0012 aft cell (20, 13): one piece at each end, dead + fresh: ', ok(20)

! 18-19. a lattice of 4 x 4 squares (sides 0.12, pitch 0.25) and three thin bars on 12 x 12 cells of 0.125: many
! multi-piece cells, splits, merges, fresh and dead pieces; translating, then rotating about the lattice centre
call lattice(c)
xg = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)]
ok(18) = body_2d('lattice, translation', pose_2d(0._R8P), pose_2d(0._R8P, [0.31_R8P, 0.17_R8P]))
ok(19) = body_2d('lattice, rotation', about(0._R8P), about(0.35_R8P))

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine rect(cv, b)
   !< The rectangle [b(1), b(2)] x [b(3), b(4)], counter-clockwise.
   type(curve_object), intent(inout) :: cv   !< Curve.
   real(R8P),          intent(in)    :: b(4) !< Extent.

   call cv%adopt_loops(u=[b(1), b(2), b(2), b(1)], v=[b(3), b(3), b(4), b(4)], loop_sizes=[4_I4P])
   endsubroutine rect

   function case_2d(label, p0, p1, exact, expected) result(good)
   !< Swept geometry and map of c between p0 and p1 on the block x^2: the map checks, and the statuses of the tracks of
   !< the middle cell equal `expected`.
   character(*),          intent(in) :: label       !< Label.
   type(rigid_pose_2d_t), intent(in) :: p0, p1      !< Poses.
   logical,               intent(in) :: exact       !< Exact path.
   integer(I4P),          intent(in) :: expected(:) !< Statuses of the middle cell's tracks.
   logical                           :: good        !< Result.
   integer(I4P)                      :: s_, t0, nt  !< Work.

   call swept_cells_2d(c, p0, p1, x, x, sw, exact=exact, status=s_, map=map)
   good = s_ == STATUS_OK .and. check_map(p0, p1, x, x, exact)
   t0 = map%first_track(2, 2) ; nt = map%tracks(2, 2)
   good = good .and. nt == size(expected, kind=I4P)
   if (good) good = all(map%status(t0:t0 + nt - 1) == expected)
   print '(A,A,A,*(I3))', '  ', label, ': middle cell statuses', map%status(t0:t0 + nt - 1)
   print '(A,A,A,L1)', '  ', label, ': ok ', good
   endfunction case_2d

   function body_2d(label, p0, p1) result(good)
   !< Map checks of c between p0 and p1 on the grid xg^2 (exact path).
   character(*),          intent(in) :: label  !< Label.
   type(rigid_pose_2d_t), intent(in) :: p0, p1 !< Poses.
   logical                           :: good   !< Result.
   integer(I4P)                      :: s_     !< Work.

   call swept_cells_2d(c, p0, p1, xg, xg, sw, status=s_, map=map)
   good = s_ == STATUS_OK .and. check_map(p0, p1, xg, xg, .true.) .and. map%untrusted_cells == 0_I4P
   print '(A,A,A,L1)', '  ', label, ': ok ', good
   endfunction body_2d

   function check_map(p0, p1, gx, gy, exact) result(good)
   !< The map checks (see the program notes) on sw and map.
   type(rigid_pose_2d_t), intent(in) :: p0, p1    !< Poses.
   real(R8P),             intent(in) :: gx(:), gy(:) !< Grid.
   logical,               intent(in) :: exact     !< Exact path: per-track residuals at round-off.
   logical                           :: good      !< Result.
   type(curve_object)                :: cm        !< Moved copy.
   type(cut_cells_2d_t)              :: k0, k1    !< Cuts of the moved copies.
   integer(I4P), allocatable         :: seen0(:), seen1(:) !< Pieces listed.
   real(R8P)                         :: fsum(4), wsum(2), rsum, dv, efa, ewa, era, esw, eres, vol !< Work.
   real(R8P), allocatable            :: ou(:,:), ov(:,:) !< Face outflows (+u, +v) summed over the tracks, per side.
   real(R8P)                         :: eid, eof, etr, d(2) !< Identity, face sides, translation errors; displacement.
   integer(I4P)                      :: i_, j_, t, k, p, nbad !< Work.

   cm = c ; call cm%move(p0) ; call cm%cut_cells(gx, gy, k0)
   cm = c ; call cm%move(p1) ; call cm%cut_cells(gx, gy, k1)
   nbad = 0_I4P
   if (map%np0 /= k0%np .or. map%np1 /= k1%np .or. map%nx /= size(gx) - 1 .or. map%ny /= size(gy) - 1) nbad = nbad + 1_I4P
   if (nbad > 0_I4P) then
      good = .false. ; print '(A)', '  sizes differ' ; return
   endif
   allocate(seen0(k0%np), seen1(k1%np)) ; seen0 = 0_I4P ; seen1 = 0_I4P
   efa = 0._R8P ; ewa = 0._R8P ; era = 0._R8P ; esw = 0._R8P ; eres = 0._R8P ; eid = 0._R8P ; eof = 0._R8P ; etr = 0._R8P
   allocate(ou(2 * size(gx), size(gy)), ov(size(gx), 2 * size(gy))) ; ou = 0._R8P ; ov = 0._R8P
   do j_ = 1, size(gy) - 1
      do i_ = 1, size(gx) - 1
         fsum = 0._R8P ; wsum = 0._R8P ; rsum = 0._R8P ; dv = 0._R8P
         do t = map%first_track(i_, j_), map%first_track(i_, j_) + map%tracks(i_, j_) - 1
            if (any(map%track_cell(:, t) /= [i_, j_])) nbad = nbad + 1_I4P
            vol = 0._R8P
            do k = map%old_start(t), map%old_start(t + 1) - 1
               p = map%old_piece(k) ; seen0(p) = seen0(p) + 1_I4P ; vol = vol + k0%piece_area(p)
               if (map%old_track(p) /= t .or. any(k0%piece_cell(:, p) /= [i_, j_])) nbad = nbad + 1_I4P
            enddo
            if (vol /= map%volume0(t)) nbad = nbad + 1_I4P
            vol = 0._R8P
            do k = map%new_start(t), map%new_start(t + 1) - 1
               p = map%new_piece(k) ; seen1(p) = seen1(p) + 1_I4P ; vol = vol + k1%piece_area(p)
               if (map%new_track(p) /= t .or. any(k1%piece_cell(:, p) /= [i_, j_])) nbad = nbad + 1_I4P
            enddo
            if (vol /= map%volume1(t)) nbad = nbad + 1_I4P
            fsum = fsum + map%faces(:, t) ; wsum = wsum + map%wall(:, t) ; rsum = rsum + map%residual(t)
            ! outflow: residual = V1 - V0 + sum_f outflow; face F from the cell above it at 2F - 1, below it at 2F
            eid = max(eid, abs(map%residual(t) - (map%volume1(t) - map%volume0(t) + sum(map%outflow(:, t)))))
            ou(2 * i_ - 1, j_) = ou(2 * i_ - 1, j_) - map%outflow(1, t)
            ou(2 * i_ + 2, j_) = ou(2 * i_ + 2, j_) + map%outflow(2, t)
            ov(i_, 2 * j_ - 1) = ov(i_, 2 * j_ - 1) - map%outflow(3, t)
            ov(i_, 2 * j_ + 2) = ov(i_, 2 * j_ + 2) + map%outflow(4, t)
            dv = dv + (map%volume1(t) - map%volume0(t))
            if (exact .and. map%status(t) /= TRACK_UNTRUSTED) eres = max(eres, abs(map%residual(t)))
         enddo
         if (map%tracks(i_, j_) == 0_I4P) cycle
         efa = max(efa, maxval(abs(fsum - [sw%face_u(i_, j_), sw%face_u(i_ + 1, j_), sw%face_v(i_, j_), &
                                          sw%face_v(i_, j_ + 1)])))
         ewa = max(ewa, maxval(abs(wsum - sw%wall(:, i_, j_))))
         era = max(era, abs(rsum - sw%residual(i_, j_)))
         esw = max(esw, abs(dv - sw%swept(i_, j_)))
      enddo
   enddo
   nbad = nbad + count(seen0 /= 1_I4P, kind=I4P) + count(seen1 /= 1_I4P, kind=I4P)
   ! a face's outflow from its two cells (interior faces) and, on the translation path, against d times the aperture
   d = p1%t - p0%t
   do j_ = 1, size(gy) - 1
      do i_ = 2, size(gx) - 1
         eof = max(eof, abs(ou(2 * i_ - 1, j_) - ou(2 * i_, j_)))
      enddo
   enddo
   do j_ = 2, size(gy) - 1
      do i_ = 1, size(gx) - 1
         eof = max(eof, abs(ov(i_, 2 * j_ - 1) - ov(i_, 2 * j_)))
      enddo
   enddo
   if (all(p0%r == p1%r) .and. exact) then
      do j_ = 1, size(gy) - 1
         do i_ = 1, size(gx)
            if (i_ > 1) etr = max(etr, abs(ou(2 * i_, j_) - d(1) * sw%face_u(i_, j_)))
         enddo
      enddo
   endif
   good = nbad == 0_I4P .and. max(efa, ewa, era, esw) < 1.e-13_R8P .and. eres < 1.e-13_R8P .and. &
          max(eid, eof, etr) < 1.e-13_R8P
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', '  tracks ', map%nt, ', multi-track cells ', map%multi_track_cells, &
         ', untrusted cells ', map%untrusted_cells, ', evaluated cells ', map%evaluated_cells, ', evaluations ', &
         map%evaluations, ', skipped nodes ', map%skipped_nodes, ', pieces ', map%np0, ' -> ', map%np1
   print '(A,I0,A,4ES9.1,A,ES9.1,A,L1)', '  bad links ', nbad, ', per-cell sums - swept (faces, wall, residual, swept) ', &
         efa, ewa, era, esw, ', max per-track residual ', eres, ', exact ', map%exact
   print '(A,3ES9.1)', '  outflow: residual - (V1 - V0 + sum outflow), face from its two sides, translation vs d L: ', &
         eid, eof, etr
   endfunction check_map

   pure function about(angle) result(p)
   !< Rotation by angle about the point (0.5, 0.5).
   real(R8P), intent(in) :: angle !< Angle.
   type(rigid_pose_2d_t) :: p     !< Pose.

   p = pose_2d(angle)
   p%t = [0.5_R8P, 0.5_R8P] - matmul(p%r, [0.5_R8P, 0.5_R8P])
   endfunction about

   subroutine lattice(cv)
   !< 4 x 4 squares of side 0.12 at pitch 0.25 from (0.065, 0.065), and three bars 0.02 thick (two horizontal, one
   !< vertical) between them; all counter-clockwise.
   type(curve_object), intent(inout) :: cv          !< Curve.
   real(R8P)                         :: uu(76), vv(76), o(2) !< Points.
   integer(I4P)                      :: a, k, n      !< Counters.

   n = 0_I4P
   do a = 0, 3
      do k = 0, 3
         o = [0.065_R8P + 0.25_R8P * k, 0.065_R8P + 0.25_R8P * a]
         uu(n + 1:n + 4) = o(1) + [0._R8P, 0.12_R8P, 0.12_R8P, 0._R8P]
         vv(n + 1:n + 4) = o(2) + [0._R8P, 0._R8P, 0.12_R8P, 0.12_R8P]
         n = n + 4_I4P
      enddo
   enddo
   uu(65:68) = [0.02_R8P, 0.98_R8P, 0.98_R8P, 0.02_R8P] ; vv(65:68) = [0.24_R8P, 0.24_R8P, 0.26_R8P, 0.26_R8P]
   uu(69:72) = [0.02_R8P, 0.98_R8P, 0.98_R8P, 0.02_R8P] ; vv(69:72) = [0.74_R8P, 0.74_R8P, 0.76_R8P, 0.76_R8P]
   uu(73:76) = [0.49_R8P, 0.51_R8P, 0.51_R8P, 0.49_R8P] ; vv(73:76) = [0.27_R8P, 0.27_R8P, 0.73_R8P, 0.73_R8P]
   call cv%adopt_loops(u=uu, v=vv, loop_sizes=[(4_I4P, k = 1, 19)])
   endsubroutine lattice

   subroutine naca(cv)
   !< NACA 0012 of chord 0.7 from (0.15, 0.5), 2 x 40 cosine-spaced points, closed trailing edge.
   type(curve_object), intent(inout) :: cv !< Curve.
   real(R8P)                         :: uu(80), vv(80), xc !< Points.
   integer(I4P)                      :: k  !< Counter.

   do k = 0, 39 ! upper side from the trailing edge, then lower side from the leading edge
      xc = 0.5_R8P * (1._R8P + cos(PI * k / 40._R8P))
      uu(k + 1) = xc ; vv(k + 1) = half(xc)
      xc = 0.5_R8P * (1._R8P - cos(PI * k / 40._R8P))
      uu(41 + k) = xc ; vv(41 + k) = -half(xc)
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
endprogram fossil_test_piece_map_2d
