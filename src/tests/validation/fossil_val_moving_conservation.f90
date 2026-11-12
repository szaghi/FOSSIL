!< FOSSIL validation: a conservative toy update across moving cut cells through the piece map (issue #102 steps 3, 2).
program fossil_val_moving_conservation
!< FOSSIL validation: a conservative toy update across moving cut cells through the piece map (issue #102 steps 3, 2).
!<
!< A passive scalar is carried by the rigid field of the body's motion on a fixed grid, one mass per cut piece, over
!< 100 steps of a periodic motion (the body returns to its first pose). Per step P(k-1) -> P(k), with the swept
!< geometry and its piece map (`swept_cells_2d(..., map=)`):
!< - one flux per grid face, F_f = q_up O_f: O_f the face's time-integrated outflow of the rigid field (the map's
!<   `outflow`, summed over the tracks of the cell below the face), q_up the mass-weighted value of the upwind cell
!<   (a cell without fluid at the start: the mean of its neighbours'; the domain boundary: zero gradient);
!< - each side's tracks take F_f in proportion to their own outflow through f (aperture weights where it sums to 0),
!<   so both sides exchange exactly F_f;
!< - a track's mass is its old pieces' mass minus its outflows; a fresh track fills by its inflow only; a track without
!<   new piece (dead, or transient) hands its remainder to a receiver: the new piece of the neighbour across the face of
!<   the track's cell with the largest aperture of the track that has the largest aperture on that face (else the
!<   largest new piece of the 3 x 3, then 5 x 5 neighbourhood);
!< - the new pieces of a track share its mass in proportion to their areas.
!< Checks, per case (gate: the first two, no untrusted cell, and no mass lost beyond 1e-14 of the total). The face
!< outflows of a multi-track cell on the Gauss path are exact per event interval while a single-track neighbour copies
!< the Gauss value: the two sides of a face then differ by the quadrature error, which the split spreads by aperture
!< (both sides still exchange F_f: conservation stays exact).
!< 1. mass: a smooth q0 = 1 + sin(2 pi x) cos(2 pi y) / 2; the total over the pieces plus the domain-boundary outflow
!<    stays constant, drift <= 1e-12 relative per step accumulated (printed);
!< 2. uniform flow moving with the body: q = 1 stays 1, max |M - V| / |cell| <= 1e-11 on the exact paths, printed with
!<    the map's largest per-track |residual| / |cell| of a step: the floor of the cut-cell areas (Green sums of unit
!<    coordinates, ~1e-16 absolute, a few 1e-12 of a 1/100 cell); |q - 1| printed over pieces of at least |cell| / 10;
!<    on the Gauss path (gated) the floor is the swept geometry's own per-cell Gauss residual, which single-track cells
!<    copy: max |M - V| / |cell| <= 2 times the largest per-track |residual| / |cell| of a step, and the tracks of
!<    multi-track cells (integrated per event interval on every path) at <= 1e-11 of a cell;
!< 3. counts over the steps: tracks per status, untrusted, dead and transient remainders (count, total |mass|);
!< 4. control, reported: two scalars a = 1 above the body's chord line (body frame, extended both ways), b = 1 below;
!<    the exact solution is the initial field carried with the body. After the period, the mass of a in pieces whose
!<    body-frame centroid is below the chord, and of b above, with the tracks; with the end-pose overlap (tracks merged
!<    when an old piece of one and a new piece of the other share a fluid interval of a cell edge); with every cell
!<    lumped (all its pieces one).
!< Bodies of 256 vertices on 96^2 cells (the exact rotation on 40^2: the swept faces dominate); argument `full`: the
!< chain bodies (1024 / 1022 vertices) at 256^2 (106^2 for the exact rotation).
!< 3-D (`swept_cells_3d(..., map=)`, the same update and checks, six faces per cell): the marching-cubes sphere and the
!< cone of `fossil_val_device_chain`, on a cube holding the whole path (24^3, 32^3), one period per case: translation
!< at 0.5 and 2 cells per step (cone 20 steps at both sizes; sphere at 24^3 6 and 3 steps, its exact swept faces
!< being the costly part; `full`: 20 steps and the sphere at 32^3 too) and the Gauss screw at 0.5 cells (24^3). Pieces
!< with fluid on a common cell face stand for the overlap. Untrusted map cells are reported, not gated (one track with
!< the cell's values, which the update handles like any other); empty and dead are keyed on pieces, never on alpha.
use fossil, only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, piece_map_2d_t, &
                   cut_cells_2d_t, to_body, TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, TRACK_SPLIT, TRACK_MERGE, &
                   TRACK_SPLIT_MERGE, TRACK_TRANSIENT, TRACK_UNTRUSTED
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, swept_cells_3d, &
                   piece_map_3d_t, cut_cells_3d_t, cut_cells_3d
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,   only : I4P, I8P, R8P
use vecfor, only : vector_R8P

implicit none

integer(I4P), parameter :: PATH_TRANSLATION = 1_I4P, PATH_SCREW = 2_I4P, PATH_GAUSS = 3_I4P !< Motions and paths.
integer(I4P), parameter :: NSTEP = 100_I4P   !< Steps (one period), 2-D.
integer(I4P), parameter :: NF = 8_I4P        !< Fields: mass, uniform, a and b with tracks, overlap, lumped.
integer(I4P), parameter :: M_TRACKS = 1_I4P, M_OVERLAP = 2_I4P, M_LUMPED = 3_I4P !< Groupings.
integer(I4P), parameter :: GSEL(NF) = [M_TRACKS, M_TRACKS, M_TRACKS, M_TRACKS, M_OVERLAP, M_OVERLAP, M_LUMPED, &
                                        M_LUMPED] !< Grouping of each field.
real(R8P),    parameter :: PI = acos(-1._R8P)
character(10), parameter :: PATHS(3) = [character(10) :: 'transl.', 'screw', 'screw, G4']
type(curve_object)       :: cb(2)            !< Bodies (body frame).
type(surface_stl_object) :: sb(2)            !< 3-D bodies (body frame).
character(12)            :: names(2)         !< Body names.
character(32)            :: arg              !< Argument.
logical                  :: full, ok         !< Full size, gate.
integer(I4P)             :: n, nv, m3(2)     !< Cells per axis, body vertices, 3-D cells per axis.
integer(I8P)             :: c0, c1, rate     !< Clock.

full = .false.
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg) ; full = trim(arg) == 'full'
endif
n = 96_I4P ; nv = 256_I4P
if (full) then
   n = 256_I4P ; nv = 1024_I4P
endif
call bodies(nv, cb)
names = [character(12) :: 'ellipse', 'NACA 0012']
call system_clock(c0, rate)
print '(A,I0,A,I0,A)', 'moving conservation (issue #102 step 3): ', n, '^2 cells, ', NSTEP, ' steps per case'
ok = .true.
ok = run(names(1), cb(1), n, PATH_TRANSLATION, 0.5_R8P) .and. ok
ok = run(names(2), cb(2), n, PATH_TRANSLATION, 0.5_R8P) .and. ok
ok = run(names(1), cb(1), n, PATH_TRANSLATION, 2._R8P) .and. ok
ok = run(names(2), cb(2), n, PATH_TRANSLATION, 2._R8P) .and. ok
ok = run(names(2), cb(2), (5_I4P * n) / 12_I4P, PATH_SCREW, 0.5_R8P) .and. ok ! exact rotation: the swept faces' cost
ok = run(names(2), cb(2), (5_I4P * n) / 12_I4P, PATH_SCREW, 2._R8P) .and. ok
ok = run(names(2), cb(2), n, PATH_GAUSS, 2._R8P) .and. ok
call bodies3(sb)
names = [character(12) :: 'mc sphere', 'cone']
m3 = [24_I4P, 32_I4P]
print '(A)', '3-D (steps per case: one period of the motion; the sphere''s exact translation is the costly part)'
ok = run3(names(2), sb(2), m3(1), PATH_TRANSLATION, 0.5_R8P, 20_I4P) .and. ok
ok = run3(names(2), sb(2), m3(1), PATH_TRANSLATION, 2._R8P, 20_I4P) .and. ok
ok = run3(names(2), sb(2), m3(2), PATH_TRANSLATION, 0.5_R8P, 20_I4P) .and. ok
ok = run3(names(2), sb(2), m3(2), PATH_TRANSLATION, 2._R8P, 20_I4P) .and. ok
ok = run3(names(1), sb(1), m3(1), PATH_TRANSLATION, 0.5_R8P, merge(20_I4P, 6_I4P, full)) .and. ok
ok = run3(names(1), sb(1), m3(1), PATH_TRANSLATION, 2._R8P, merge(20_I4P, 3_I4P, full)) .and. ok
if (full) then ! the sphere at 32^3: 42 s (0.5 cells) and 117 s (2 cells) for 20 steps
   ok = run3(names(1), sb(1), m3(2), PATH_TRANSLATION, 0.5_R8P, 20_I4P) .and. ok
   ok = run3(names(1), sb(1), m3(2), PATH_TRANSLATION, 2._R8P, 20_I4P) .and. ok
endif
ok = run3(names(1), sb(1), m3(1), PATH_GAUSS, 0.5_R8P, merge(20_I4P, 10_I4P, full)) .and. ok
ok = run3(names(2), sb(2), m3(1), PATH_GAUSS, 0.5_R8P, 20_I4P) .and. ok
call system_clock(c1)
print '(A,F8.1,A)', 'total time ', real(c1 - c0, R8P) / rate, ' s'
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   function run(label, body, n, path, delta) result(good)
   !< One case: NSTEP steps of the periodic motion with a largest per-step displacement of about delta cells.
   character(*),       intent(in) :: label   !< Body name.
   type(curve_object), intent(in) :: body    !< Body (body frame).
   integer(I4P),       intent(in) :: n       !< Cells per axis.
   integer(I4P),       intent(in) :: path    !< PATH_*.
   real(R8P),          intent(in) :: delta   !< Displacement per step [cells].
   logical                        :: good    !< Gate.
   type(rigid_pose_2d_t)          :: p(0:NSTEP) !< Poses.
   type(curve_object)             :: c       !< Moved copy.
   type(cut_cells_2d_t)           :: k0, k1  !< Cuts at the step's ends.
   type(swept_cells_2d_t)         :: sw      !< Swept geometry.
   type(piece_map_2d_t)           :: map     !< Piece map.
   real(R8P),    allocatable      :: x(:), m0(:,:), m1(:,:) !< Grid, masses per field [np, NF].
   integer(I4P), allocatable      :: grp(:,:)                !< Group of each track per grouping [nt, 3].
   real(R8P)                      :: h, tot0(NF), bout(NF), lost(NF), rmass(NF), drift, efree, qfree, eest, dmax, w(2)
   real(R8P)                      :: leak(2, 3), cel, floor, fmulti
   integer(I4P)                   :: k, st, f, nst(-1:7), nrem(NF), nmulti, nuntr, nskip, ndead_nt, nsw
   integer(I8P)                   :: t0, t1, rate

   call system_clock(t0, rate)
   x = [(0.0031_R8P - 0.02_R8P * 0.9938_R8P + 1.04_R8P * 0.9938_R8P * k / real(n, R8P), k = 0, n)] ; h = x(2) - x(1)
   cel = h * h
   call poses(body, h, path, delta, p, dmax)
   c = body ; call c%move(p(0)) ; call c%cut_cells(x, x, k0)
   allocate(m0(k0%np, NF))
   call initial(k0, p(0), m0)
   do f = 1_I4P, NF
      tot0(f) = sum(m0(:, f))
   enddo
   bout = 0._R8P ; lost = 0._R8P ; rmass = 0._R8P ; nrem = 0_I4P ; drift = 0._R8P ; efree = 0._R8P ; qfree = 0._R8P
   eest = 0._R8P ; floor = 0._R8P ; fmulti = 0._R8P ; nsw = 0_I4P ; nst = 0_I4P ; nmulti = 0_I4P ; nuntr = 0_I4P ; nskip = 0_I4P
   good = .true. ; ndead_nt = 0_I4P
   do k = 1_I4P, NSTEP
      call swept_cells_2d(body, p(k - 1), p(k), x, x, sw, exact=path /= PATH_GAUSS, status=st, map=map)
      c = body ; call c%move(p(k)) ; call c%cut_cells(x, x, k1)
      if (st /= 0_I4P .or. map%np0 /= k0%np .or. map%np1 /= k1%np) then
         print '(A,I0)', '  map and cuts disagree at step ', k ; good = .false. ; return
      endif
      nst = nst + map%n_status ; nmulti = nmulti + map%multi_track_cells ; nuntr = nuntr + map%untrusted_cells
      nskip = nskip + map%skipped_nodes
      if (path == PATH_GAUSS) eest = max(eest, sw%error_estimate)
      call groupings(map, k0, k1, grp)
      allocate(m1(k1%np, NF))
      do f = 1_I4P, NF
         call advance(map, k0, k1, x, f, p(k), p(0), grp(:, GSEL(f)), m0(:, f), m1(:, f), bout(f), lost(f), nrem(f), &
                      rmass(f), nsw)
      enddo
      drift = max(drift, abs(sum(m1(:, 1)) + bout(1) - tot0(1)) / tot0(1))
      efree = max(efree, maxval(abs(m1(:, 2) - k1%piece_area)) / cel)
      qfree = max(qfree, maxval(abs(m1(:, 2) / k1%piece_area - 1._R8P), mask=k1%piece_area >= 0.1_R8P * cel))
      floor = max(floor, maxval(abs(map%residual)) / cel)
      do f = 1_I4P, map%nt ! the floor of the multi-track cells' tracks (integrated per event interval on every path)
         if (map%tracks(map%track_cell(1, f), map%track_cell(2, f)) > 1_I4P) fmulti = max(fmulti, abs(map%residual(f)) / cel)
      enddo
      call move_alloc(m1, m0) ; k0 = k1
      ndead_nt = ndead_nt + map%n_status(TRACK_DEAD) + map%n_status(TRACK_TRANSIENT)
   enddo
   call leakage(k0, p(NSTEP), m0, leak)
   call system_clock(t1)
   print '(A,A,1X,I0,A,A,A,F4.1,A,F5.2,A,F6.1,A)', '- ', trim(label), n, '^2, ', trim(PATHS(path)), ', delta ', delta, &
         ' cells (measured ', dmax / h, '), ', real(t1 - t0, R8P) / rate, ' s'
   print '(A,ES9.2,A,ES9.2,A)', '  1. mass drift (max over steps, relative) ', drift, ' (bound ', 1e-12_R8P * NSTEP, ')'
   print '(A,ES9.2,A,ES9.2)', '  2. uniform flow: max |M - V| / |cell| ', efree, ', max |q - 1| (pieces >= |cell| / 10) ', &
         qfree
   print '(A,ES9.2,A,ES9.2,A,ES9.2)', '     largest per-track |residual| / |cell| of a step ', floor, &
         ' (tracks of multi-track cells ', fmulti, '), nq-vs-2nq face estimate ', eest
   print '(A,8I8)', '  3. tracks continuing/fresh/dead/split/merge/split-merge/transient/untrusted ', &
         nst([TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, TRACK_SPLIT, TRACK_MERGE, TRACK_SPLIT_MERGE, TRACK_TRANSIENT, &
              TRACK_UNTRUSTED])
   print '(A,I0,A,I0,A,I0,A,I0)', '     multi-track cells ', nmulti, ', untrusted cells ', nuntr, ', skipped nodes ', nskip, &
         ', Gauss-Seidel sweeps (largest) ', nsw
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2,A,ES9.2)', '     dead + transient tracks ', ndead_nt, ', remainders handed on: ', &
         nrem(1), ' (mass field, total |mass| ', rmass(1), '; uniform field ', rmass(2), '), lost ', maxval(abs(lost))
   w = [sum(m0(:, 3)), sum(m0(:, 4))]
   print '(A,3(ES10.3,A,ES10.3,A))', '  4. leak a below / b above: tracks ', leak(1, 1), ' / ', leak(2, 1), &
         ', overlap ', leak(1, 2), ' / ', leak(2, 2), ', lumped ', leak(1, 3), ' / ', leak(2, 3), ''
   good = drift <= 1e-12_R8P * NSTEP .and. all(lost <= 1e-14_R8P * abs(tot0)) .and. nuntr == 0_I4P
   if (path /= PATH_GAUSS) then
      good = good .and. efree <= 1e-11_R8P
   else ! the floor is the cells' own Gauss residual (single-track cells copy it); multi-track cells' tracks are exact
      good = good .and. efree <= 2._R8P * floor .and. fmulti <= 1e-11_R8P
   endif
   print '(A,L1)', '  ok ', good
   endfunction run

   subroutine poses(body, h, path, delta, p, dmax)
   !< P_k, k = 0..NSTEP, one period: translation A (sin wk, sin 2wk / 2) and, on the screw paths, a rotation by
   !< T sin wk about (0.5, 0.5); A and T give a largest vertex displacement per step of about delta h (35 % from the
   !< translation, 65 % from the rotation on the screw paths). dmax: the measured largest displacement.
   type(curve_object),    intent(in)  :: body          !< Body.
   real(R8P),             intent(in)  :: h, delta      !< Cell size, displacement per step [cells].
   integer(I4P),          intent(in)  :: path          !< PATH_*.
   type(rigid_pose_2d_t), intent(out) :: p(0:NSTEP)    !< Poses.
   real(R8P),             intent(out) :: dmax          !< Largest vertex displacement per step.
   real(R8P)                          :: wk, at, th, r, uv(2), a(2), b(2) !< Work.
   integer(I4P)                       :: k, g          !< Counters.

   wk = 2._R8P * PI / NSTEP ; r = 0._R8P
   do g = 1_I4P, body%get_vertices_number()
      uv = body%get_vertex(g) ; r = max(r, sqrt((uv(1) - 0.5_R8P)**2 + (uv(2) - 0.5_R8P)**2))
   enddo
   if (path == PATH_TRANSLATION) then
      at = delta * h / (wk * sqrt(1.25_R8P)) ; th = 0._R8P
   else
      at = 0.35_R8P * delta * h / (wk * sqrt(1.25_R8P)) ; th = 0.65_R8P * delta * h / (wk * r)
   endif
   do k = 0_I4P, NSTEP
      p(k) = pose_2d(th * sin(wk * k))
      p(k)%t = [0.5_R8P, 0.5_R8P] - matmul(p(k)%r, [0.5_R8P, 0.5_R8P]) + at * [sin(wk * k), 0.5_R8P * sin(2._R8P * wk * k)]
   enddo
   dmax = 0._R8P
   do k = 1_I4P, NSTEP
      do g = 1_I4P, body%get_vertices_number()
         uv = body%get_vertex(g)
         a = matmul(p(k)%r, uv) + p(k)%t ; b = matmul(p(k - 1)%r, uv) + p(k - 1)%t
         dmax = max(dmax, sqrt(sum((a - b)**2)))
      enddo
   enddo
   endsubroutine poses

   subroutine initial(k0, p0, m)
   !< Initial masses of the fields per piece: q0 (mass), 1 (uniform), and a / b (three copies each) by the side of the
   !< chord line v = 0.5 (body frame) of the piece's centroid.
   type(cut_cells_2d_t),  intent(in)  :: k0     !< Cut.
   type(rigid_pose_2d_t), intent(in)  :: p0     !< Pose.
   real(R8P),             intent(out) :: m(:,:) !< Masses [np, NF].
   integer(I4P)                       :: q, f   !< Piece, field.

   do q = 1_I4P, k0%np
      do f = 1_I4P, NF
         m(q, f) = k0%piece_area(q) * exact_q(f, p0, p0, k0%piece_centroid(:, q))
      enddo
   enddo
   endsubroutine initial

   function exact_q(fid, pk, p0, w) result(qv)
   !< The exact solution of field fid at the world point w at pose pk: the initial field carried with the body (the
   !< rigid field leaves the fluid at rest in the body frame). Fields: q0 = 1 + sin(2 pi x) cos(2 pi y) / 2 at the
   !< point's position at p0; 1; a = 1 above the chord line v = 0.5 of the body frame, b = 1 below (three copies).
   integer(I4P),          intent(in) :: fid    !< Field.
   type(rigid_pose_2d_t), intent(in) :: pk, p0 !< Pose, initial pose.
   real(R8P),             intent(in) :: w(2)   !< World point.
   real(R8P)                         :: qv     !< Value.
   real(R8P)                         :: b(2), w0(2) !< Body point, its position at p0.

   b = to_body(pk, w)
   select case (fid)
   case (1_I4P)
      w0 = matmul(p0%r, b) + p0%t
      qv = 1._R8P + 0.5_R8P * sin(2._R8P * PI * w0(1)) * cos(2._R8P * PI * w0(2))
   case (2_I4P)
      qv = 1._R8P
   case (3_I4P, 5_I4P, 7_I4P)
      qv = merge(1._R8P, 0._R8P, b(2) >= 0.5_R8P)
   case default
      qv = merge(0._R8P, 1._R8P, b(2) >= 0.5_R8P)
   endselect
   endfunction exact_q

   subroutine groupings(map, k0, k1, grp)
   !< The group of every track per grouping: tracks (itself), overlap (tracks of a cell united when an old piece of one
   !< and a new piece of the other share a fluid interval of the same cell edge, with positive length), lumped (the
   !< cell's first track). Group ids are track ids.
   type(piece_map_2d_t),      intent(in)  :: map    !< Map.
   type(cut_cells_2d_t),      intent(in)  :: k0, k1 !< Cuts.
   integer(I4P), allocatable, intent(out) :: grp(:,:) !< Groups [nt, 3].
   integer(I4P), allocatable :: par(:)              !< Union-find over tracks.
   integer(I4P)              :: i, j, p, q, ra, rb, t !< Work.

   allocate(grp(map%nt, 3), par(map%nt))
   do t = 1_I4P, map%nt
      par(t) = t ; grp(t, M_TRACKS) = t
      grp(t, M_LUMPED) = map%first_track(map%track_cell(1, t), map%track_cell(2, t))
   enddo
   do j = 1_I4P, map%ny
      do i = 1_I4P, map%nx
         do p = k0%first_piece(i, j), k0%first_piece(i, j) + k0%pieces(i, j) - 1_I4P
            do q = k1%first_piece(i, j), k1%first_piece(i, j) + k1%pieces(i, j) - 1_I4P
               if (.not. overlap(k0, p, k1, q)) cycle
               ra = root(par, map%old_track(p)) ; rb = root(par, map%new_track(q))
               if (ra /= rb) par(max(ra, rb)) = min(ra, rb)
            enddo
         enddo
      enddo
   enddo
   do t = 1_I4P, map%nt
      grp(t, M_OVERLAP) = root(par, t)
   enddo
   endsubroutine groupings

   pure integer(I4P) function root(par, t) result(r)
   !< Union-find root.
   integer(I4P), intent(in) :: par(:) !< Parents.
   integer(I4P), intent(in) :: t      !< Element.

   r = t
   do while (par(r) /= r)
      r = par(r)
   enddo
   endfunction root

   pure logical function overlap(k0, p, k1, q)
   !< Old piece p and new piece q share a fluid interval of a cell edge with positive length.
   type(cut_cells_2d_t), intent(in) :: k0, k1 !< Cuts.
   integer(I4P),         intent(in) :: p, q   !< Pieces.
   integer(I4P)                     :: a, b   !< Intervals.

   overlap = .false.
   do a = k0%piece_iv_start(p), k0%piece_iv_start(p + 1) - 1_I4P
      do b = k1%piece_iv_start(q), k1%piece_iv_start(q + 1) - 1_I4P
         if (k0%piece_iv_edge(a) /= k1%piece_iv_edge(b)) cycle
         if (min(k0%piece_iv_hi(a), k1%piece_iv_hi(b)) > max(k0%piece_iv_lo(a), k1%piece_iv_lo(b))) then
            overlap = .true. ; return
         endif
      enddo
   enddo
   endfunction overlap

   subroutine advance(map, k0, k1, x, fid, pk, p0, grp, m0, m1, bout, lost, nrem, rmass, nsweep)
   !< One step of one field with one grouping of the tracks (see the program notes). The upwind values are the cells'
   !< mass-weighted values at the end of the step (backward Euler, so stable at any displacement per step): per cell,
   !< V1 q + sum_out O q = M0 + sum_in O q_neighbour, solved by Gauss-Seidel sweeps in the four directions to 1e-15.
   !< Inflow through the domain boundary carries the exact solution at the face's middle (`exact_q`).
   type(piece_map_2d_t), intent(in)    :: map    !< Map.
   type(cut_cells_2d_t), intent(in)    :: k0, k1 !< Cuts.
   real(R8P),            intent(in)    :: x(:)   !< Grid (both axes).
   integer(I4P),         intent(in)    :: fid    !< Field.
   type(rigid_pose_2d_t), intent(in)   :: pk, p0 !< Pose at the step's end, initial pose.
   integer(I4P),         intent(in)    :: grp(:) !< Group of each track.
   real(R8P),            intent(in)    :: m0(:)  !< Old masses.
   real(R8P),            intent(out)   :: m1(:)  !< New masses.
   real(R8P),            intent(inout) :: bout   !< Domain-boundary outflow, accumulated.
   real(R8P),            intent(inout) :: lost   !< Mass with nowhere to go, accumulated.
   integer(I4P),         intent(inout) :: nrem   !< Remainders handed on, accumulated.
   real(R8P),            intent(inout) :: rmass  !< Their total |mass|, accumulated.
   integer(I4P),         intent(inout) :: nsweep !< Largest number of sweeps.
   integer(I4P), parameter   :: DI(4) = [-1, 1, 0, 0], DJ(4) = [0, 0, -1, 1] !< Neighbour across faces -u, +u, -v, +v.
   real(R8P),    allocatable :: q(:,:), mc(:,:), v1(:,:), so(:,:,:), sa(:,:,:), ob(:,:,:), fu(:,:), fv(:,:), mg(:), &
                                vg(:), ag(:,:), qb(:,:,:) !< Work.
   integer(I4P), allocatable :: ng(:)             !< New pieces per group.
   integer(I4P)              :: nx, ny, i, j, t, p, f, g, sw, ii, jj, ia, ib, is, ja, jb, js !< Work.
   real(R8P)                 :: o, fo, sh, a, b, dq, qn !< Work.

   nx = map%nx ; ny = map%ny
   allocate(q(nx, ny), mc(nx, ny), v1(nx, ny), so(4, nx, ny), sa(4, nx, ny), ob(4, nx, ny), fu(nx + 1, ny), &
            fv(nx, ny + 1), mg(map%nt), vg(map%nt), ng(map%nt), ag(4, map%nt), qb(4, nx, ny))
   qb = 0._R8P ! the exact values at the middles of the domain-boundary faces
   do j = 1_I4P, ny
      qb(1, 1, j) = exact_q(fid, pk, p0, [x(1), 0.5_R8P * (x(j) + x(j + 1))])
      qb(2, nx, j) = exact_q(fid, pk, p0, [x(nx + 1), 0.5_R8P * (x(j) + x(j + 1))])
   enddo
   do i = 1_I4P, nx
      qb(3, i, 1) = exact_q(fid, pk, p0, [0.5_R8P * (x(i) + x(i + 1)), x(1)])
      qb(4, i, ny) = exact_q(fid, pk, p0, [0.5_R8P * (x(i) + x(i + 1)), x(ny + 1)])
   enddo
   mc = 0._R8P ; v1 = 0._R8P
   do p = 1_I4P, k0%np
      mc(k0%piece_cell(1, p), k0%piece_cell(2, p)) = mc(k0%piece_cell(1, p), k0%piece_cell(2, p)) + m0(p)
   enddo
   do p = 1_I4P, k1%np
      v1(k1%piece_cell(1, p), k1%piece_cell(2, p)) = v1(k1%piece_cell(1, p), k1%piece_cell(2, p)) + k1%piece_area(p)
   enddo
   ! per cell and face: the tracks' outflows and apertures; one outflow per face (from the cell below it, above it on
   ! the domain's low side), seen outward from each of its cells (ob)
   so = 0._R8P ; sa = 0._R8P
   do t = 1_I4P, map%nt
      i = map%track_cell(1, t) ; j = map%track_cell(2, t)
      so(:, i, j) = so(:, i, j) + map%outflow(:, t) ; sa(:, i, j) = sa(:, i, j) + map%faces(:, t)
   enddo
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         ob(1, i, j) = so(1, i, j) ; ob(3, i, j) = so(3, i, j)
         if (i > 1_I4P) ob(1, i, j) = -so(2, i - 1, j)
         if (j > 1_I4P) ob(3, i, j) = -so(4, i, j - 1)
         ob(2, i, j) = so(2, i, j) ; ob(4, i, j) = so(4, i, j)
      enddo
   enddo
   ! the end-of-step cell values
   q = 0._R8P
   do j = 1_I4P, ny
      do i = 1_I4P, nx
         if (v1(i, j) > 0._R8P) q(i, j) = mc(i, j) / v1(i, j)
      enddo
   enddo
   do sw = 1_I4P, 4000_I4P
      dq = 0._R8P
      ia = 1_I4P ; ib = nx ; is = 1_I4P ; ja = 1_I4P ; jb = ny ; js = 1_I4P
      if (mod(sw, 2_I4P) == 0_I4P) then
         ia = nx ; ib = 1_I4P ; is = -1_I4P
      endif
      if (mod((sw - 1_I4P) / 2_I4P, 2_I4P) == 1_I4P) then
         ja = ny ; jb = 1_I4P ; js = -1_I4P
      endif
      do j = ja, jb, js
         do i = ia, ib, is
            a = v1(i, j) ; b = mc(i, j)
            do f = 1_I4P, 4_I4P
               o = ob(f, i, j)
               if (o > 0._R8P) then
                  a = a + o
               elseif (o < 0._R8P) then
                  ii = i + DI(f) ; jj = j + DJ(f)
                  if (ii < 1 .or. ii > nx .or. jj < 1 .or. jj > ny) then
                     b = b - o * qb(f, i, j) ! inflow through the domain boundary: the exact value
                  else
                     b = b - o * q(ii, jj)
                  endif
               endif
            enddo
            if (a > 0._R8P) then
               qn = b / a ; dq = max(dq, abs(qn - q(i, j)) / max(1._R8P, abs(qn))) ; q(i, j) = qn
            endif
         enddo
      enddo
      if (dq <= 1e-15_R8P) exit
   enddo
   nsweep = max(nsweep, sw)
   ! face fluxes (+u, +v), upwind at the end of the step
   do j = 1_I4P, ny
      do i = 1_I4P, nx + 1_I4P
         if (i > 1_I4P) then
            o = so(2, i - 1, j)
         else
            o = -so(1, 1, j)
         endif
         if (o > 0._R8P .and. i > 1_I4P) then
            fu(i, j) = o * q(i - 1, j)
         elseif (o > 0._R8P) then
            fu(i, j) = o * qb(1, 1, j)
         elseif (i <= nx) then
            fu(i, j) = o * q(i, j)
         else
            fu(i, j) = o * qb(2, nx, j)
         endif
      enddo
   enddo
   do j = 1_I4P, ny + 1_I4P
      do i = 1_I4P, nx
         if (j > 1_I4P) then
            o = so(4, i, j - 1)
         else
            o = -so(3, i, 1)
         endif
         if (o > 0._R8P .and. j > 1_I4P) then
            fv(i, j) = o * q(i, j - 1)
         elseif (o > 0._R8P) then
            fv(i, j) = o * qb(3, i, 1)
         elseif (j <= ny) then
            fv(i, j) = o * q(i, j)
         else
            fv(i, j) = o * qb(4, i, ny)
         endif
      enddo
   enddo
   do j = 1_I4P, ny
      bout = bout + fu(nx + 1, j) - fu(1, j)
   enddo
   do i = 1_I4P, nx
      bout = bout + fv(i, ny + 1) - fv(i, 1)
   enddo
   do j = 1_I4P, ny ! faces of cells without tracks: their flux has nowhere to go
      do i = 1_I4P, nx
         if (map%tracks(i, j) > 0_I4P) cycle
         lost = lost + abs(fu(i, j)) + abs(fu(i + 1, j)) + abs(fv(i, j)) + abs(fv(i, j + 1))
      enddo
   enddo
   ! groups: old mass minus the tracks' parts of the face fluxes. A track takes q times its own outflow; the face's
   ! flux and the sum over the cell's tracks differ by round-off (the two sides' outflow sums), which goes by aperture.
   mg = 0._R8P ; vg = 0._R8P ; ng = 0_I4P ; ag = 0._R8P
   do p = 1_I4P, k0%np
      g = grp(map%old_track(p)) ; mg(g) = mg(g) + m0(p)
   enddo
   do t = 1_I4P, map%nt
      i = map%track_cell(1, t) ; j = map%track_cell(2, t) ; g = grp(t)
      do f = 1_I4P, 4_I4P
         select case (f)
         case (1_I4P) ; fo = -fu(i, j)
         case (2_I4P) ; fo = fu(i + 1, j)
         case (3_I4P) ; fo = -fv(i, j)
         case default ; fo = fv(i, j + 1)
         endselect
         if (ob(f, i, j) /= 0._R8P) then ! the value carried: fo / ob (the upwind q)
            sh = map%outflow(f, t) * (fo / ob(f, i, j))
            if (sa(f, i, j) > 0._R8P) then
               sh = sh + (fo - so(f, i, j) * (fo / ob(f, i, j))) * (map%faces(f, t) / sa(f, i, j))
            else
               sh = sh + (fo - so(f, i, j) * (fo / ob(f, i, j))) / map%tracks(i, j)
            endif
         elseif (sa(f, i, j) > 0._R8P) then
            sh = fo * (map%faces(f, t) / sa(f, i, j))
         else
            sh = fo / map%tracks(i, j)
         endif
         mg(g) = mg(g) - sh
      enddo
      ag(:, g) = ag(:, g) + map%faces(:, t)
   enddo
   do p = 1_I4P, k1%np
      g = grp(map%new_track(p)) ; vg(g) = vg(g) + k1%piece_area(p) ; ng(g) = ng(g) + 1_I4P
   enddo
   m1 = 0._R8P
   do p = 1_I4P, k1%np
      g = grp(map%new_track(p))
      if (vg(g) > 0._R8P) then
         m1(p) = mg(g) * (k1%piece_area(p) / vg(g))
      else
         m1(p) = mg(g) / ng(g)
      endif
   enddo
   ! groups without new pieces: the remainder to a receiver
   do t = 1_I4P, map%nt
      if (grp(t) /= t .or. ng(t) > 0_I4P) cycle
      if (mg(t) == 0._R8P) cycle
      p = receiver(map, k1, t, ag(:, t))
      if (p == 0_I4P) then
         lost = lost + abs(mg(t))
      else
         m1(p) = m1(p) + mg(t) ; nrem = nrem + 1_I4P ; rmass = rmass + abs(mg(t))
      endif
   enddo
   endsubroutine advance

   function receiver(map, k1, t, ap) result(r)
   !< The new piece receiving the remainder of group t (cell of track t): across the faces of the cell in decreasing
   !< order of the group's aperture, the neighbour's new piece with the largest aperture on that face; else the largest
   !< new piece of the 3 x 3, then 5 x 5 neighbourhood; 0 if none.
   type(piece_map_2d_t), intent(in) :: map    !< Map.
   type(cut_cells_2d_t), intent(in) :: k1     !< New cut.
   integer(I4P),         intent(in) :: t      !< Group (its track id).
   real(R8P),            intent(in) :: ap(4)  !< The group's apertures (-u, +u, -v, +v).
   integer(I4P)                     :: r      !< Receiver piece.
   real(R8P)                        :: a(4), best !< Work.
   integer(I4P)                     :: i, j, f, k, ii, jj, side, q, rad !< Work.

   r = 0_I4P
   i = map%track_cell(1, t) ; j = map%track_cell(2, t) ; a = ap
   do k = 1_I4P, 4_I4P
      f = maxloc(a, dim=1) ; if (a(f) <= 0._R8P) exit
      a(f) = -1._R8P
      select case (f)
      case (1_I4P) ; ii = i - 1_I4P ; jj = j ; side = 2_I4P
      case (2_I4P) ; ii = i + 1_I4P ; jj = j ; side = 1_I4P
      case (3_I4P) ; ii = i ; jj = j - 1_I4P ; side = 4_I4P
      case default ; ii = i ; jj = j + 1_I4P ; side = 3_I4P
      endselect
      if (ii < 1 .or. ii > map%nx .or. jj < 1 .or. jj > map%ny) cycle
      best = 0._R8P
      do q = k1%first_piece(ii, jj), k1%first_piece(ii, jj) + k1%pieces(ii, jj) - 1_I4P
         if (k1%piece_aperture(side, q) > best) then
            best = k1%piece_aperture(side, q) ; r = q
         endif
      enddo
      if (r > 0_I4P) return
   enddo
   do rad = 1_I4P, 2_I4P
      best = 0._R8P
      do jj = max(1_I4P, j - rad), min(map%ny, j + rad)
         do ii = max(1_I4P, i - rad), min(map%nx, i + rad)
            do q = k1%first_piece(ii, jj), k1%first_piece(ii, jj) + k1%pieces(ii, jj) - 1_I4P
               if (k1%piece_area(q) > best) then
                  best = k1%piece_area(q) ; r = q
               endif
            enddo
         enddo
      enddo
      if (r > 0_I4P) return
   enddo
   endfunction receiver

   subroutine leakage(k, p, m, leak)
   !< Mass of a in pieces whose body-frame centroid is below the chord, and of b above, per grouping.
   type(cut_cells_2d_t),  intent(in)  :: k         !< Final cut.
   type(rigid_pose_2d_t), intent(in)  :: p         !< Final pose.
   real(R8P),             intent(in)  :: m(:,:)    !< Masses.
   real(R8P),             intent(out) :: leak(2, 3) !< (a below, b above) x (tracks, overlap, lumped).
   real(R8P)                          :: cb(2)     !< Body-frame centroid.
   integer(I4P)                       :: q, g      !< Work.

   leak = 0._R8P
   do q = 1_I4P, k%np
      cb = to_body(p, k%piece_centroid(:, q))
      do g = 1_I4P, 3_I4P
         if (cb(2) < 0.5_R8P) then
            leak(1, g) = leak(1, g) + m(q, 1 + 2 * g)
         else
            leak(2, g) = leak(2, g) + m(q, 2 + 2 * g)
         endif
      enddo
   enddo
   endsubroutine leakage

   function run3(label, body, m, path, delta, ns) result(good)
   !< One 3-D case: ns steps of the periodic motion (one period) with a largest per-step displacement of about delta
   !< cells, on m^3 cells over a cube holding the body's whole path. Untrusted map cells are reported, not gated: the map
   !< gives them one track with the cell's values, which the update handles like any other.
   character(*),             intent(in) :: label   !< Body name.
   type(surface_stl_object), intent(in) :: body    !< Body (body frame).
   integer(I4P),             intent(in) :: m       !< Cells per axis.
   integer(I4P),             intent(in) :: path    !< PATH_*.
   real(R8P),                intent(in) :: delta   !< Displacement per step [cells].
   integer(I4P),             intent(in) :: ns      !< Steps (one period).
   logical                              :: good    !< Gate.
   type(rigid_pose_t), allocatable      :: p(:)    !< Poses [0:ns].
   type(surface_stl_object)             :: s       !< Moved copy.
   type(cut_cells_3d_t)                 :: k0, k1  !< Cuts at the step's ends.
   type(swept_cells_3d_t)               :: sw      !< Swept geometry.
   type(piece_map_3d_t)                 :: map     !< Piece map.
   real(R8P),    allocatable            :: g(:), m0(:,:), m1(:,:) !< Grid (all axes), masses per field [np, NF].
   integer(I4P), allocatable            :: grp(:,:) !< Group of each track per grouping [nt, 3].
   real(R8P)                            :: h, cen(3), tot0(NF), bout(NF), lost(NF), rmass(NF), drift, efree, qfree, eest
   real(R8P)                            :: dmax, leak(2, 3), cel, floor, fmulti
   integer(I4P)                         :: k, st, f, nst(-1:7), nrem(NF), nmulti, nuntr, nskip, ndead_nt, nsw, nbadcut
   integer(I8P)                         :: t0, t1, rate

   call system_clock(t0, rate)
   allocate(p(0:ns))
   call poses3(body, m, path, delta, ns, p, g, cen, dmax)
   h = g(2) - g(1) ; cel = h**3
   s = body ; call s%move(p(0)) ; call cut_cells_3d(s, g, g, g, k0)
   allocate(m0(k0%np, NF))
   do k = 1_I4P, k0%np
      do f = 1_I4P, NF
         m0(k, f) = k0%piece_volume(k) * exact_q3(f, p(0), p(0), cen, k0%piece_centroid(:, k))
      enddo
   enddo
   do f = 1_I4P, NF
      tot0(f) = sum(m0(:, f))
   enddo
   bout = 0._R8P ; lost = 0._R8P ; rmass = 0._R8P ; nrem = 0_I4P ; drift = 0._R8P ; efree = 0._R8P ; qfree = 0._R8P
   eest = 0._R8P ; floor = 0._R8P ; fmulti = 0._R8P ; nsw = 0_I4P ; nst = 0_I4P ; nmulti = 0_I4P ; nuntr = 0_I4P
   nskip = 0_I4P ; ndead_nt = 0_I4P ; nbadcut = 0_I4P
   do k = 1_I4P, ns
      call swept_cells_3d(body, p(k - 1), p(k), g, g, g, sw, exact=path /= PATH_GAUSS, status=st, map=map)
      s = body ; call s%move(p(k)) ; call cut_cells_3d(s, g, g, g, k1)
      if (st /= 0_I4P .or. map%np0 /= k0%np .or. map%np1 /= k1%np) then
         print '(A,I0)', '  map and cuts disagree at step ', k ; good = .false. ; return
      endif
      nst = nst + map%n_status ; nmulti = nmulti + map%multi_track_cells ; nuntr = nuntr + map%untrusted_cells
      nskip = nskip + map%skipped_nodes ; nbadcut = nbadcut + count(k1%pieces < 0_I4P)
      if (path == PATH_GAUSS) eest = max(eest, sw%error_estimate)
      floor = max(floor, maxval(abs(map%residual)) / cel)
      do f = 1_I4P, map%nt
         if (map%tracks(map%track_cell(1, f), map%track_cell(2, f), map%track_cell(3, f)) > 1_I4P) &
            fmulti = max(fmulti, abs(map%residual(f)) / cel)
      enddo
      call groupings3(map, k0, k1, grp)
      allocate(m1(k1%np, NF))
      do f = 1_I4P, NF
         call advance3(map, k0, k1, g, f, p(k), p(0), cen, grp(:, GSEL(f)), m0(:, f), m1(:, f), bout(f), lost(f), &
                       nrem(f), rmass(f), nsw)
      enddo
      drift = max(drift, abs(sum(m1(:, 1)) + bout(1) - tot0(1)) / tot0(1))
      efree = max(efree, maxval(abs(m1(:, 2) - k1%piece_volume)) / cel)
      qfree = max(qfree, maxval(abs(m1(:, 2) / k1%piece_volume - 1._R8P), mask=k1%piece_volume >= 0.1_R8P * cel))
      call move_alloc(m1, m0) ; k0 = k1
      ndead_nt = ndead_nt + map%n_status(TRACK_DEAD) + map%n_status(TRACK_TRANSIENT)
   enddo
   call leakage3(k0, p(ns), cen, m0, leak)
   call system_clock(t1)
   print '(A,A,1X,I0,A,A,A,F4.1,A,F5.2,A,I0,A,F6.1,A)', '- ', trim(label), m, '^3, ', trim(PATHS(path)), ', delta ', &
         delta, ' cells (measured ', dmax / h, '), ', ns, ' steps, ', real(t1 - t0, R8P) / rate, ' s'
   print '(A,ES9.2,A,ES9.2,A)', '  1. mass drift (max over steps, relative) ', drift, ' (bound ', 1e-12_R8P * ns, ')'
   print '(A,ES9.2,A,ES9.2)', '  2. uniform flow: max |M - V| / |cell| ', efree, ', max |q - 1| (pieces >= |cell| / 10) ', &
         qfree
   print '(A,ES9.2,A,ES9.2,A,ES9.2)', '     largest per-track |residual| / |cell| of a step ', floor, &
         ' (tracks of multi-track cells ', fmulti, '), nq-vs-2nq face estimate ', eest
   print '(A,8I8)', '  3. tracks continuing/fresh/dead/split/merge/split-merge/transient/untrusted ', &
         nst([TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, TRACK_SPLIT, TRACK_MERGE, TRACK_SPLIT_MERGE, TRACK_TRANSIENT, &
              TRACK_UNTRUSTED])
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', '     multi-track cells ', nmulti, ', untrusted cells ', nuntr, ', skipped nodes ', &
         nskip, ', untrusted end cuts (pieces -1) ', nbadcut, ', Gauss-Seidel sweeps (largest) ', nsw
   print '(A,I0,A,I0,A,ES9.2,A,ES9.2,A,ES9.2)', '     dead + transient tracks ', ndead_nt, ', remainders handed on: ', &
         nrem(1), ' (mass field, total |mass| ', rmass(1), '; uniform field ', rmass(2), '), lost ', maxval(abs(lost))
   print '(A,3(ES10.3,A,ES10.3,A))', '  4. leak a below / b above: tracks ', leak(1, 1), ' / ', leak(2, 1), &
         ', overlap ', leak(1, 2), ' / ', leak(2, 2), ', lumped ', leak(1, 3), ' / ', leak(2, 3), ''
   good = drift <= 1e-12_R8P * ns .and. all(lost <= 1e-14_R8P * abs(tot0))
   if (path /= PATH_GAUSS) then
      good = good .and. efree <= 2._R8P * floor
   else
      good = good .and. efree <= 2._R8P * floor .and. fmulti <= 1e-11_R8P
   endif
   print '(A,L1)', '  ok ', good
   endfunction run3

   subroutine poses3(body, m, path, delta, ns, p, g, cen, dmax)
   !< P_k, k = 0..ns, one period: translation A (sin wk, sin 2wk / 2, sin 3wk / 3) and, on the screw paths, a
   !< rotation by T sin wk about the axis (0.2, 1, 0.3) through the body's box centre; A and T give a largest vertex
   !< displacement per step of about delta h (35 % from the translation, 65 % from the rotation on the screw paths; from
   !< the measured per-step displacements of the unit motions, so few steps per period reach it too). The
   !< grid g (every axis): m cells over the cube of half-side 1.05 (r + |A|) about the centre (r: the body's largest
   !< distance from it; grown if the centre is off the diagonal), which holds the whole path. dmax: the measured
   !< largest displacement.
   type(surface_stl_object),  intent(in)  :: body        !< Body.
   integer(I4P),              intent(in)  :: m           !< Cells per axis.
   integer(I4P),              intent(in)  :: path        !< PATH_*.
   real(R8P),                 intent(in)  :: delta       !< Displacement per step [cells].
   integer(I4P),              intent(in)  :: ns          !< Steps.
   type(rigid_pose_t),        intent(out) :: p(0:ns)     !< Poses.
   real(R8P), allocatable,    intent(out) :: g(:)        !< Grid.
   real(R8P),                 intent(out) :: cen(3)      !< Body box centre.
   real(R8P),                 intent(out) :: dmax        !< Largest vertex displacement per step.
   real(R8P), allocatable                 :: v(:,:)      !< Vertices [3, nv].
   type(vector_R8P)                       :: b0, b1      !< Box.
   real(R8P)                              :: wk, r, ca, ct, at, th, half, a(3), b(3), ax(3), dt1, dr1 !< Work.
   integer(I4P)                           :: k, q, nv    !< Counters.

   b0 = body%get_bmin() ; b1 = body%get_bmax()
   cen = 0.5_R8P * [b0%x + b1%x, b0%y + b1%y, b0%z + b1%z]
   call body_vertices(body, v) ; nv = size(v, 2, kind=I4P)
   r = 0._R8P
   do q = 1_I4P, nv
      r = max(r, sqrt(sum((v(:, q) - cen)**2)))
   enddo
   wk = 2._R8P * PI / ns
   ca = 1._R8P ; ct = 0._R8P ! shares of the displacement: translation, rotation
   if (path /= PATH_TRANSLATION) then
      ca = 0.35_R8P ; ct = 0.65_R8P
   endif
   ax = [0.2_R8P, 1._R8P, 0.3_R8P]
   ! the largest per-step displacement of the unit motions (translation amplitude 1; rotation amplitude 1 rad)
   dt1 = 0._R8P ; dr1 = 0._R8P
   do k = 1_I4P, ns
      dt1 = max(dt1, sqrt(sum((tpath(wk, k) - tpath(wk, k - 1_I4P))**2)))
      do q = 1_I4P, nv
         a = matmul(rot_unit(ax, wk, k), v(:, q) - cen) ; b = matmul(rot_unit(ax, wk, k - 1_I4P), v(:, q) - cen)
         dr1 = max(dr1, sqrt(sum((a - b)**2)))
      enddo
   enddo
   ! at = ca delta h / dt1, h = 2 half / m, half = 1.05 (r + at): solve for at; then th = ct delta h / dr1
   at = (2.1_R8P * ca * delta / (m * dt1)) * r / (1._R8P - 2.1_R8P * ca * delta / (m * dt1))
   half = 1.05_R8P * (r + at) + maxval(abs(cen - cen(1))) ! the same edges on every axis: cover an off-diagonal centre
   g = [(cen(1) - half + 2._R8P * half * k / m, k = 0, m)]
   th = ct * delta * (2._R8P * half / m) / dr1
   do k = 0_I4P, ns
      p(k) = pose_from_axis_angle(ax, th * sin(wk * k))
      p(k)%t = cen - matmul(p(k)%r, cen) + at * tpath(wk, k)
   enddo
   dmax = 0._R8P
   do k = 1_I4P, ns
      do q = 1_I4P, nv
         a = matmul(p(k)%r, v(:, q)) + p(k)%t ; b = matmul(p(k - 1)%r, v(:, q)) + p(k - 1)%t
         dmax = max(dmax, sqrt(sum((a - b)**2)))
      enddo
   enddo
   endsubroutine poses3

   pure function tpath(wk, k) result(t)
   !< The unit translation path at step k: (sin wk, sin 2wk / 2, sin 3wk / 3).
   real(R8P),    intent(in) :: wk !< Angular step.
   integer(I4P), intent(in) :: k  !< Step.
   real(R8P)                :: t(3) !< Translation.

   t = [sin(wk * k), 0.5_R8P * sin(2._R8P * wk * k), sin(3._R8P * wk * k) / 3._R8P]
   endfunction tpath

   function rot_unit(ax, wk, k) result(r)
   !< The rotation of the unit rotation path at step k: angle sin wk about ax.
   real(R8P),    intent(in) :: ax(3), wk !< Axis, angular step.
   integer(I4P), intent(in) :: k         !< Step.
   real(R8P)                :: r(3, 3)   !< Rotation.
   type(rigid_pose_t)       :: p         !< Pose.

   p = pose_from_axis_angle(ax, sin(wk * k)) ; r = p%r
   endfunction rot_unit

   subroutine body_vertices(body, v)
   !< The facets' vertices of the body (with repetitions) [3, 3 nf].
   type(surface_stl_object), intent(in)  :: body   !< Body.
   real(R8P), allocatable,   intent(out) :: v(:,:) !< Vertices.
   integer(I4P)                          :: q, nf  !< Work.

   nf = body%get_facets_number()
   allocate(v(3, 3 * nf))
   do q = 1_I4P, nf
      associate (fc => body%facet_at(q))
         v(:, 3 * q - 2) = [fc%vertex(1)%x, fc%vertex(1)%y, fc%vertex(1)%z]
         v(:, 3 * q - 1) = [fc%vertex(2)%x, fc%vertex(2)%y, fc%vertex(2)%z]
         v(:, 3 * q) = [fc%vertex(3)%x, fc%vertex(3)%y, fc%vertex(3)%z]
      end associate
   enddo
   endsubroutine body_vertices

   function exact_q3(fid, pk, p0, cen, w) result(qv)
   !< The exact solution of field fid at the world point w at pose pk (the 3-D twin of `exact_q`): q0 =
   !< 1 + sin(2 pi x) cos(2 pi y) cos(2 pi z) / 2 at the point's position at p0 relative to the body centre; 1; a = 1
   !< above the plane z = cen(3) of the body frame, b = 1 below.
   integer(I4P),       intent(in) :: fid    !< Field.
   type(rigid_pose_t), intent(in) :: pk, p0 !< Pose, initial pose.
   real(R8P),          intent(in) :: cen(3) !< Body box centre (body frame).
   real(R8P),          intent(in) :: w(3)   !< World point.
   real(R8P)                      :: qv     !< Value.
   real(R8P)                      :: b(3), w0(3) !< Body point, its position at p0.

   b = matmul(transpose(pk%r), w - pk%t)
   select case (fid)
   case (1_I4P)
      w0 = matmul(p0%r, b) + p0%t - cen
      qv = 1._R8P + 0.5_R8P * sin(2._R8P * PI * w0(1)) * cos(2._R8P * PI * w0(2)) * cos(2._R8P * PI * w0(3))
   case (2_I4P)
      qv = 1._R8P
   case (3_I4P, 5_I4P, 7_I4P)
      qv = merge(1._R8P, 0._R8P, b(3) >= cen(3))
   case default
      qv = merge(0._R8P, 1._R8P, b(3) >= cen(3))
   endselect
   endfunction exact_q3

   subroutine groupings3(map, k0, k1, grp)
   !< `groupings` in 3-D. Overlap: an old piece of one track and a new piece of another in the same cell both carry
   !< fluid on one of the cell's faces (`piece_face`; coarser than the 2-D edge intervals: it may join pieces that only
   !< share a face, which makes the overlap control lump more, never less).
   type(piece_map_3d_t),      intent(in)  :: map      !< Map.
   type(cut_cells_3d_t),      intent(in)  :: k0, k1   !< Cuts.
   integer(I4P), allocatable, intent(out) :: grp(:,:) !< Groups [nt, 3].
   integer(I4P), allocatable :: par(:)                !< Union-find over tracks.
   integer(I4P)              :: i, j, l, p, q, ra, rb, t !< Work.

   allocate(grp(map%nt, 3), par(map%nt))
   do t = 1_I4P, map%nt
      par(t) = t ; grp(t, M_TRACKS) = t
      grp(t, M_LUMPED) = map%first_track(map%track_cell(1, t), map%track_cell(2, t), map%track_cell(3, t))
   enddo
   do l = 1_I4P, map%nz
      do j = 1_I4P, map%ny
         do i = 1_I4P, map%nx
            if (k0%pieces(i, j, l) <= 0_I4P .or. k1%pieces(i, j, l) <= 0_I4P) cycle
            do p = k0%first_piece(i, j, l), k0%first_piece(i, j, l) + k0%pieces(i, j, l) - 1_I4P
               do q = k1%first_piece(i, j, l), k1%first_piece(i, j, l) + k1%pieces(i, j, l) - 1_I4P
                  if (.not. any(k0%piece_face(:, p) > 0._R8P .and. k1%piece_face(:, q) > 0._R8P)) cycle
                  ra = root(par, map%old_track(p)) ; rb = root(par, map%new_track(q))
                  if (ra /= rb) par(max(ra, rb)) = min(ra, rb)
               enddo
            enddo
         enddo
      enddo
   enddo
   do t = 1_I4P, map%nt
      grp(t, M_OVERLAP) = root(par, t)
   enddo
   endsubroutine groupings3

   subroutine advance3(map, k0, k1, g, fid, pk, p0, cen, grp, m0, m1, bout, lost, nrem, rmass, nsweep)
   !< `advance` in 3-D: one step of one field with one grouping (see the program notes), six faces per cell.
   type(piece_map_3d_t), intent(in)    :: map    !< Map.
   type(cut_cells_3d_t), intent(in)    :: k0, k1 !< Cuts.
   real(R8P),            intent(in)    :: g(:)   !< Grid (every axis).
   integer(I4P),         intent(in)    :: fid    !< Field.
   type(rigid_pose_t),   intent(in)    :: pk, p0 !< Pose at the step's end, initial pose.
   real(R8P),            intent(in)    :: cen(3) !< Body box centre.
   integer(I4P),         intent(in)    :: grp(:) !< Group of each track.
   real(R8P),            intent(in)    :: m0(:)  !< Old masses.
   real(R8P),            intent(out)   :: m1(:)  !< New masses.
   real(R8P),            intent(inout) :: bout   !< Domain-boundary outflow, accumulated.
   real(R8P),            intent(inout) :: lost   !< Mass with nowhere to go, accumulated.
   integer(I4P),         intent(inout) :: nrem   !< Remainders handed on, accumulated.
   real(R8P),            intent(inout) :: rmass  !< Their total |mass|, accumulated.
   integer(I4P),         intent(inout) :: nsweep !< Largest number of sweeps.
   integer(I4P), parameter   :: DI(6) = [-1, 1, 0, 0, 0, 0], DJ(6) = [0, 0, -1, 1, 0, 0], DK(6) = [0, 0, 0, 0, -1, 1]
   real(R8P),    allocatable :: q(:,:,:), mc(:,:,:), v1(:,:,:), so(:,:,:,:), sa(:,:,:,:), ob(:,:,:,:), qb(:,:,:,:)
   real(R8P),    allocatable :: fo6(:,:,:,:), mg(:), vg(:), ag(:,:) !< Outward face fluxes per cell [6, n, n, n].
   integer(I4P), allocatable :: ng(:)             !< New pieces per group.
   integer(I4P)              :: n, i, j, l, t, p, f, gg, sw, ii, jj, ll, ia, ib, is, ja, jb, js, la, lb, ls !< Work.
   real(R8P)                 :: o, sh, a, b, dq, qn, w(3) !< Work.

   n = map%nx
   allocate(q(n, n, n), mc(n, n, n), v1(n, n, n), so(6, n, n, n), sa(6, n, n, n), ob(6, n, n, n), qb(6, n, n, n), &
            fo6(6, n, n, n), mg(map%nt), vg(map%nt), ng(map%nt), ag(6, map%nt))
   mc = 0._R8P ; v1 = 0._R8P
   do p = 1_I4P, k0%np
      mc(k0%piece_cell(1, p), k0%piece_cell(2, p), k0%piece_cell(3, p)) = &
         mc(k0%piece_cell(1, p), k0%piece_cell(2, p), k0%piece_cell(3, p)) + m0(p)
   enddo
   do p = 1_I4P, k1%np
      v1(k1%piece_cell(1, p), k1%piece_cell(2, p), k1%piece_cell(3, p)) = &
         v1(k1%piece_cell(1, p), k1%piece_cell(2, p), k1%piece_cell(3, p)) + k1%piece_volume(p)
   enddo
   so = 0._R8P ; sa = 0._R8P
   do t = 1_I4P, map%nt
      i = map%track_cell(1, t) ; j = map%track_cell(2, t) ; l = map%track_cell(3, t)
      so(:, i, j, l) = so(:, i, j, l) + map%outflow(:, t) ; sa(:, i, j, l) = sa(:, i, j, l) + map%faces(:, t)
   enddo
   ! one outflow per face, from the cell on its low side (its own on the domain's low boundary), seen outward from
   ! each of its cells; the exact values at the middles of the domain-boundary faces
   qb = 0._R8P
   do l = 1_I4P, n
      do j = 1_I4P, n
         do i = 1_I4P, n
            ob(:, i, j, l) = so(:, i, j, l)
            if (i > 1_I4P) ob(1, i, j, l) = -so(2, i - 1, j, l)
            if (j > 1_I4P) ob(3, i, j, l) = -so(4, i, j - 1, l)
            if (l > 1_I4P) ob(5, i, j, l) = -so(6, i, j, l - 1)
            do f = 1_I4P, 6_I4P
               ii = i + DI(f) ; jj = j + DJ(f) ; ll = l + DK(f)
               if (ii >= 1 .and. ii <= n .and. jj >= 1 .and. jj <= n .and. ll >= 1 .and. ll <= n) cycle
               w = 0.5_R8P * ([g(i), g(j), g(l)] + [g(i + 1), g(j + 1), g(l + 1)])
               select case (f)
               case (1_I4P) ; w(1) = g(i)
               case (2_I4P) ; w(1) = g(i + 1)
               case (3_I4P) ; w(2) = g(j)
               case (4_I4P) ; w(2) = g(j + 1)
               case (5_I4P) ; w(3) = g(l)
               case default ; w(3) = g(l + 1)
               endselect
               qb(f, i, j, l) = exact_q3(fid, pk, p0, cen, w)
            enddo
         enddo
      enddo
   enddo
   ! the end-of-step cell values (backward Euler), Gauss-Seidel sweeps in alternating directions
   q = 0._R8P
   where (v1 > 0._R8P) q = mc / v1
   do sw = 1_I4P, 4000_I4P
      dq = 0._R8P
      ia = 1_I4P ; ib = n ; is = 1_I4P ; ja = 1_I4P ; jb = n ; js = 1_I4P ; la = 1_I4P ; lb = n ; ls = 1_I4P
      if (mod(sw, 2_I4P) == 0_I4P) then
         ia = n ; ib = 1_I4P ; is = -1_I4P
      endif
      if (mod((sw - 1_I4P) / 2_I4P, 2_I4P) == 1_I4P) then
         ja = n ; jb = 1_I4P ; js = -1_I4P
      endif
      if (mod((sw - 1_I4P) / 4_I4P, 2_I4P) == 1_I4P) then
         la = n ; lb = 1_I4P ; ls = -1_I4P
      endif
      do l = la, lb, ls
         do j = ja, jb, js
            do i = ia, ib, is
               a = v1(i, j, l) ; b = mc(i, j, l)
               do f = 1_I4P, 6_I4P
                  o = ob(f, i, j, l)
                  if (o > 0._R8P) then
                     a = a + o
                  elseif (o < 0._R8P) then
                     ii = i + DI(f) ; jj = j + DJ(f) ; ll = l + DK(f)
                     if (ii < 1 .or. ii > n .or. jj < 1 .or. jj > n .or. ll < 1 .or. ll > n) then
                        b = b - o * qb(f, i, j, l)
                     else
                        b = b - o * q(ii, jj, ll)
                     endif
                  endif
               enddo
               if (a > 0._R8P) then
                  qn = b / a ; dq = max(dq, abs(qn - q(i, j, l)) / max(1._R8P, abs(qn))) ; q(i, j, l) = qn
               endif
            enddo
         enddo
      enddo
      if (dq <= 1e-15_R8P) exit
   enddo
   nsweep = max(nsweep, sw)
   ! outward face fluxes per cell: the face's outflow (from its low-side cell) times the upwind value, the same number
   ! seen from both cells
   do l = 1_I4P, n
      do j = 1_I4P, n
         do i = 1_I4P, n
            do f = 1_I4P, 6_I4P
               o = ob(f, i, j, l)
               ii = i + DI(f) ; jj = j + DJ(f) ; ll = l + DK(f)
               if (o > 0._R8P) then
                  fo6(f, i, j, l) = o * q(i, j, l)
               elseif (ii < 1 .or. ii > n .or. jj < 1 .or. jj > n .or. ll < 1 .or. ll > n) then
                  fo6(f, i, j, l) = o * qb(f, i, j, l)
               else
                  fo6(f, i, j, l) = o * q(ii, jj, ll)
               endif
               if (ii < 1 .or. ii > n .or. jj < 1 .or. jj > n .or. ll < 1 .or. ll > n) bout = bout + fo6(f, i, j, l)
            enddo
            if (map%tracks(i, j, l) == 0_I4P) lost = lost + sum(abs(fo6(:, i, j, l)))
         enddo
      enddo
   enddo
   mg = 0._R8P ; vg = 0._R8P ; ng = 0_I4P ; ag = 0._R8P
   do p = 1_I4P, k0%np
      gg = grp(map%old_track(p)) ; mg(gg) = mg(gg) + m0(p)
   enddo
   do t = 1_I4P, map%nt
      i = map%track_cell(1, t) ; j = map%track_cell(2, t) ; l = map%track_cell(3, t) ; gg = grp(t)
      do f = 1_I4P, 6_I4P
         if (ob(f, i, j, l) /= 0._R8P) then ! q times the track's own outflow, the sides' difference by aperture
            sh = map%outflow(f, t) * (fo6(f, i, j, l) / ob(f, i, j, l))
            if (sa(f, i, j, l) > 0._R8P) then
               sh = sh + (fo6(f, i, j, l) - so(f, i, j, l) * (fo6(f, i, j, l) / ob(f, i, j, l))) * &
                         (map%faces(f, t) / sa(f, i, j, l))
            else
               sh = sh + (fo6(f, i, j, l) - so(f, i, j, l) * (fo6(f, i, j, l) / ob(f, i, j, l))) / map%tracks(i, j, l)
            endif
         elseif (sa(f, i, j, l) > 0._R8P) then
            sh = fo6(f, i, j, l) * (map%faces(f, t) / sa(f, i, j, l))
         else
            sh = fo6(f, i, j, l) / map%tracks(i, j, l)
         endif
         mg(gg) = mg(gg) - sh
      enddo
      ag(:, gg) = ag(:, gg) + map%faces(:, t)
   enddo
   do p = 1_I4P, k1%np
      gg = grp(map%new_track(p)) ; vg(gg) = vg(gg) + k1%piece_volume(p) ; ng(gg) = ng(gg) + 1_I4P
   enddo
   m1 = 0._R8P
   do p = 1_I4P, k1%np
      gg = grp(map%new_track(p))
      if (vg(gg) > 0._R8P) then
         m1(p) = mg(gg) * (k1%piece_volume(p) / vg(gg))
      else
         m1(p) = mg(gg) / ng(gg)
      endif
   enddo
   do t = 1_I4P, map%nt
      if (grp(t) /= t .or. ng(t) > 0_I4P) cycle
      if (mg(t) == 0._R8P) cycle
      p = receiver3(map, k1, t, ag(:, t))
      if (p == 0_I4P) then
         lost = lost + abs(mg(t))
      else
         m1(p) = m1(p) + mg(t) ; nrem = nrem + 1_I4P ; rmass = rmass + abs(mg(t))
      endif
   enddo
   endsubroutine advance3

   function receiver3(map, k1, t, ap) result(r)
   !< `receiver` in 3-D: across the faces of the cell in decreasing order of the group's aperture, the neighbour's new
   !< piece with the largest fluid fraction on that face; else the largest new piece of the 3^3, then 5^3 block.
   type(piece_map_3d_t), intent(in) :: map   !< Map.
   type(cut_cells_3d_t), intent(in) :: k1    !< New cut.
   integer(I4P),         intent(in) :: t     !< Group (its track id).
   real(R8P),            intent(in) :: ap(6) !< The group's apertures.
   integer(I4P)                     :: r     !< Receiver piece.
   integer(I4P), parameter          :: DI(6) = [-1, 1, 0, 0, 0, 0], DJ(6) = [0, 0, -1, 1, 0, 0], DK(6) = [0, 0, 0, 0, -1, 1]
   integer(I4P), parameter          :: OPP(6) = [2, 1, 4, 3, 6, 5] !< The neighbour's face.
   real(R8P)                        :: a(6), best !< Work.
   integer(I4P)                     :: i, j, l, f, k, ii, jj, ll, q, rad, n !< Work.

   r = 0_I4P ; n = map%nx
   i = map%track_cell(1, t) ; j = map%track_cell(2, t) ; l = map%track_cell(3, t) ; a = ap
   do k = 1_I4P, 6_I4P
      f = maxloc(a, dim=1) ; if (a(f) <= 0._R8P) exit
      a(f) = -1._R8P
      ii = i + DI(f) ; jj = j + DJ(f) ; ll = l + DK(f)
      if (ii < 1 .or. ii > n .or. jj < 1 .or. jj > n .or. ll < 1 .or. ll > n) cycle
      if (k1%pieces(ii, jj, ll) <= 0_I4P) cycle
      best = 0._R8P
      do q = k1%first_piece(ii, jj, ll), k1%first_piece(ii, jj, ll) + k1%pieces(ii, jj, ll) - 1_I4P
         if (k1%piece_face(OPP(f), q) > best) then
            best = k1%piece_face(OPP(f), q) ; r = q
         endif
      enddo
      if (r > 0_I4P) return
   enddo
   do rad = 1_I4P, 2_I4P
      best = 0._R8P
      do ll = max(1_I4P, l - rad), min(n, l + rad)
         do jj = max(1_I4P, j - rad), min(n, j + rad)
            do ii = max(1_I4P, i - rad), min(n, i + rad)
               if (k1%pieces(ii, jj, ll) <= 0_I4P) cycle
               do q = k1%first_piece(ii, jj, ll), k1%first_piece(ii, jj, ll) + k1%pieces(ii, jj, ll) - 1_I4P
                  if (k1%piece_volume(q) > best) then
                     best = k1%piece_volume(q) ; r = q
                  endif
               enddo
            enddo
         enddo
      enddo
      if (r > 0_I4P) return
   enddo
   endfunction receiver3

   subroutine leakage3(k, p, cen, m, leak)
   !< `leakage` in 3-D: the mass of a in pieces whose body-frame centroid is below the plane z = cen(3), of b above.
   type(cut_cells_3d_t), intent(in)  :: k          !< Final cut.
   type(rigid_pose_t),   intent(in)  :: p          !< Final pose.
   real(R8P),            intent(in)  :: cen(3)     !< Body box centre.
   real(R8P),            intent(in)  :: m(:,:)     !< Masses.
   real(R8P),            intent(out) :: leak(2, 3) !< (a below, b above) x (tracks, overlap, lumped).
   real(R8P)                         :: cb(3)      !< Body-frame centroid.
   integer(I4P)                      :: q, gg      !< Work.

   leak = 0._R8P
   do q = 1_I4P, k%np
      cb = matmul(transpose(p%r), k%piece_centroid(:, q) - p%t)
      do gg = 1_I4P, 3_I4P
         if (cb(3) < cen(3)) then
            leak(1, gg) = leak(1, gg) + m(q, 1 + 2 * gg)
         else
            leak(2, gg) = leak(2, gg) + m(q, 2 + 2 * gg)
         endif
      enddo
   enddo
   endsubroutine leakage3

   subroutine bodies3(sb)
   !< The marching-cubes sphere (24^3 samples) and the sign bench's cone centred on its box, as `fossil_val_device_chain`.
   type(surface_stl_object), intent(inout) :: sb(2) !< Bodies.
   type(body_t), allocatable               :: bb(:) !< Sign-bench bodies.
   type(vector_R8P)                        :: b0, b1 !< Box.

   call mc_body(BODY_SPHERE, 24_I4P, sb(1))
   call bench_bodies(bb) ; call to_surface(bb(7)%soup, sb(2))
   b0 = sb(2)%get_bmin() ; b1 = sb(2)%get_bmax()
   call sb(2)%translate(delta=-0.5_R8P * (b0 + b1)) ; call sb(2)%analyze
   endsubroutine bodies3

   subroutine bodies(nv, cb)
   !< The ellipse (centre (0.5, 0.5), semi-axes 0.4, 0.25) and the NACA 0012 (chord 0.7 from (0.15, 0.5)) of nv vertices.
   integer(I4P),       intent(in)    :: nv    !< Vertices.
   type(curve_object), intent(inout) :: cb(2) !< Bodies.
   real(R8P), allocatable            :: u(:), v(:) !< Points.
   real(R8P)                         :: xi    !< Abscissa.
   integer(I4P)                      :: k, m, j !< Counters.

   allocate(u(nv), v(nv))
   do k = 1_I4P, nv
      u(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / nv) ; v(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / nv)
   enddo
   call cb(1)%adopt_loops(u=u, v=v, loop_sizes=[nv])
   m = nv / 2_I4P ; j = 0_I4P ! upper side from the trailing edge, then the lower side (counter-clockwise)
   do k = m, 0_I4P, -1_I4P
      xi = 0.5_R8P * (1._R8P - cos(PI * k / m)) ; j = j + 1_I4P ; u(j) = xi ; v(j) = half(xi)
   enddo
   do k = 1_I4P, m - 1_I4P
      xi = 0.5_R8P * (1._R8P - cos(PI * k / m)) ; j = j + 1_I4P ; u(j) = xi ; v(j) = -half(xi)
   enddo
   u(1:j) = 0.15_R8P + 0.7_R8P * u(1:j) ; v(1:j) = 0.5_R8P + 0.7_R8P * v(1:j)
   call cb(2)%adopt_loops(u=u(1:j), v=v(1:j), loop_sizes=[j])
   endsubroutine bodies

   pure function half(xc) result(yc)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xc !< Chord position.
   real(R8P)             :: yc !< Half thickness.

   yc = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1036_R8P * xc**4)
   endfunction half
endprogram fossil_val_moving_conservation
