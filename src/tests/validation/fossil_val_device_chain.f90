!< FOSSIL validation: the per-step chain of a moving-body cut-cell solver, CPU against device (issue #100 step 6b).
program fossil_val_device_chain
!< FOSSIL validation: the per-step chain of a moving-body cut-cell solver, CPU against device (issue #100 step 6b).
!<
!< A body moves rigidly through the absolute poses P_0 = identity, P_1, ..., P_nstep (body frame). Time step k:
!< 1. move: the body-frame copy moved to P_k (`s = body ; call s%move(P_k)`): no drift over the steps, and the same
!<    operations as the end-pose copies of `swept_cells_*`, so the swept alpha1 is the cut alpha bit for bit;
!< 2. tag: `box_classify_many(_device)`, one box per cell (the default fast-winding-number sign);
!< 3. cut: the cut cells and their piece links: 3-D `cut_cells_3d(..., links=)` / `cut_cells_3d_device(..., links=)`;
!<    2-D `curve%cut_cells` + `cut_links_2d` / `curve%cut_cells_device(resident=.true.)` + `cut_links_2d_device`;
!< 4. swept: `swept_cells_3d(_device)` / `swept_cells_2d(_device)` of the body between P_{k-1} and P_k. Path "Gauss": a
!<    rotation (2 degrees per step about an axis through the block centre) plus a translation per step, Gauss-Legendre
!<    nq = 4 (`exact=.false.`, the solver path; it also computes the 8-node error estimate). Path "translation": a pure
!<    translation per step, the exact translation path (3-D at the smallest size only: tens of seconds per host call).
!< 5. map (issue #105 step 4, path "Gauss" only): the same swept call with `map=`, the piece map of the step
!<    (`swept_cells_3d(_device)(map=)` / `swept_cells_2d(_device)(map=)`).
!< Bodies: the step-0 bodies of issue #100, the marching-cubes sphere and the cone of the sign bench (3-D); the
!< 1024-gon ellipse and the 1022-vertex NACA 0012 of the curve bench (2-D).
!<
!< Check (the first pass over the steps, which is also the warm-up of the device paths): at every step the device tags,
!< the device cut and links copied back and resident (then `update_host`) and the device swept geometry equal the CPU's
!< bit for bit (`fossil_cut_cells_compare`, `fossil_swept_3d_compare`, `fossil_swept_compare`), and on the Gauss path
!< the device piece map equals the CPU's (`fossil_piece_map_compare`) and the swept geometry of the call with `map=`
!< equals that of the call without (CPU and device): `device chain vs serial mismatches` must be 0. Also on the CPU:
!< the swept alpha0 / alpha1 equal the cut alpha of the previous / current step bit for bit (`swept vs cut alpha
!< mismatches`, 0).
!<
!< Timing: per step and per stage (move, tag, cut + links, swept, map) and the whole step, best and median of nrep runs
!< of nstep steps, in three modes: the CPU; the device with every output copied back; the device with the cut and the
!< links resident (the tags, the swept geometry and the map have no resident form: copied back). The map is timed on
!< the Gauss path at the smallest size (every size with `full`): both swept calls, without and with `map=`, run in the
!< step, and the map column is the best (median) of the call with `map=` less the best (median) of the call without (on
!< the device the call with `map=` also cuts the end poses in full, 3-D, or on the host, 2-D: part of the map's cost).
!< The step is move + tag + cut + links + the swept call with `map=` where the map is timed, else the swept call. Each
!< device mode starts with one untimed step (the pageable-copy warm-up of issue #100 step 6a). The device modes map the
!< moved surface (curve) with `enter_device` inside the move stage. Each table prints the load average (/proc/loadavg):
!< the machine is shared.
!<
!< Sizes (3 repetitions everywhere). Built with OpenACC (`fobis build --mode tests-nvidia-chain --varset local_nvf`), the
!< default run (the device gate; 10 min measured on a shared RTX 4070 and CPU, about 7 min estimated idle): 3-D 32^3 and
!< 64^3 (3 steps), the exact translation at 32^3 (2 steps, the cone), 2-D 256^2 (5 steps). `fossil_val_device_chain
!< full` adds the step-0 sizes 128^3 and 1024^2 (3 steps, Gauss only) and the sphere's exact translation (about 35 min
!< before the map stage; with it, longer: not measured). Without OpenACC (gfortran, the validation tier) the device
!< routines run their loops on the host: 3-D 8^3 (the translation 6^3, the cone), 2-D 48^2, 2 steps (seconds).
!< Usage: fossil_val_device_chain [full]
use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_links_t, &
                                     cut_cells_3d, cut_cells_3d_device, cut_links_2d, cut_links_2d_device, rigid_pose_t, &
                                     rigid_pose_2d_t, pose_identity, pose_from_axis_angle, pose_2d, swept_cells_2d_t, &
                                     swept_cells_3d_t, swept_cells_2d, swept_cells_2d_device, swept_cells_3d, swept_cells_3d_device
use fossil,                   only : piece_map_2d_t, piece_map_3d_t
use fossil,                   only : BOX_CUT
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use fossil_cut_cells_compare, only : cut_cells_2d_mismatches, cut_cells_3d_mismatches, cut_links_mismatches
use fossil_swept_3d_compare,  only : swept_cells_3d_mismatches
use fossil_swept_compare,     only : swept_cells_2d_mismatches
use fossil_piece_map_compare, only : piece_map_2d_mismatches, piece_map_3d_mismatches
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P
use, intrinsic :: iso_fortran_env, only : output_unit

implicit none

#ifdef _OPENACC
logical,      parameter   :: ON_DEVICE = .true.                    !< Built with OpenACC.
integer(I4P), parameter   :: SIZE3T = 32_I4P                       !< 3-D cells per axis of the exact translation.
integer(I4P), parameter   :: NSTEP3 = 3_I4P, NSTEP2 = 5_I4P        !< Time steps.
#else
logical,      parameter   :: ON_DEVICE = .false.                   !< Built without OpenACC: host twins.
integer(I4P), parameter   :: SIZE3T = 6_I4P                        !< 3-D cells per axis of the exact translation.
integer(I4P), parameter   :: NSTEP3 = 2_I4P, NSTEP2 = 2_I4P        !< Time steps.
#endif
integer(I4P), parameter   :: NREP = 3_I4P                          !< Timed repetitions.
integer(I4P), parameter   :: PATH_GAUSS = 1_I4P, PATH_TRANSLATION = 2_I4P !< Motions.
real(R8P),    parameter   :: PI = acos(-1._R8P)                    !< Pi.
character(12), parameter  :: MODES(3) = [character(12) :: 'CPU', 'copied back', 'resident'] !< Timing modes.
type(surface_stl_object)  :: sb(2)                                 !< 3-D bodies (body frame).
type(curve_object)        :: cb(2)                                 !< 2-D bodies (body frame).
type(body_t), allocatable :: bb(:)                                 !< Sign-bench bodies.
character(12)             :: names3(2), names2(2)                  !< Body names.
character(32)             :: arg                                   !< Argument.
real(R8P),    allocatable :: uu(:), vv(:)                          !< Work.
integer(I4P), allocatable :: sizes3(:), sizes2(:)                  !< Cells per axis (3-D, 2-D).
integer(I4P)              :: bad, abad, b, g, k                    !< Work.
logical                   :: full                                  !< The full step-0 table.

full = .false.
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg) ; full = trim(arg) == 'full'
endif
if (.not. ON_DEVICE) then
   sizes3 = [8_I4P] ; sizes2 = [48_I4P]
elseif (full) then
   sizes3 = [32_I4P, 64_I4P, 128_I4P] ; sizes2 = [256_I4P, 1024_I4P]
else
   sizes3 = [32_I4P, 64_I4P] ; sizes2 = [256_I4P]
endif
print '(A,L1,A,L1)', 'device chain (issue #100 step 6b, map #105); built with OpenACC: ', ON_DEVICE, ', full: ', full
call mc_body(BODY_SPHERE, 24_I4P, sb(1))
call bench_bodies(bb) ; call to_surface_centred(bb(7), sb(2))
names3 = [character(12) :: 'mc sphere', 'cone']
allocate(uu(1024), vv(1024))
do k = 1, 1024
   uu(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 1024) ; vv(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 1024)
enddo
call cb(1)%adopt_loops(u=uu, v=vv, loop_sizes=[1024_I4P])
deallocate(uu, vv) ; allocate(uu(1022), vv(1022))
call naca0012(511_I4P, uu, vv)
uu = 0.15_R8P + 0.7_R8P * uu ; vv = 0.5_R8P + 0.7_R8P * vv
call cb(2)%adopt_loops(u=uu, v=vv, loop_sizes=[1022_I4P])
names2 = [character(12) :: 'ellipse', 'NACA 0012']
bad = 0_I4P ; abad = 0_I4P
print '(A)', 'per-step times [s] (per stage: best / median of the repetitions); speed-up = CPU best / mode best (total)'
do g = 1, size(sizes3)
   do b = 1, 2
      call chain_3d(names3(b), sb(b), sizes3(g), PATH_GAUSS, NSTEP3, bad, abad)
   enddo
enddo
do b = merge(1, 2, full), 2 ! the sphere's exact translation only in full: 0.4 s per step at 6^3 (host), 6 s at 32^3
   call chain_3d(names3(b), sb(b), SIZE3T, PATH_TRANSLATION, min(NSTEP3, 2_I4P), bad, abad)
enddo
do g = 1, size(sizes2)
   do b = 1, 2
      ! 1024^2 (full): 3 steps, Gauss only (the 2-D swept faces are host work, about 5 s per call on either path)
      call chain_2d(names2(b), cb(b), sizes2(g), PATH_GAUSS, merge(3_I4P, NSTEP2, sizes2(g) > 256_I4P), bad, abad)
      if (sizes2(g) <= 256_I4P) call chain_2d(names2(b), cb(b), sizes2(g), PATH_TRANSLATION, NSTEP2, bad, abad)
   enddo
enddo
print '(A,I0)', 'swept vs cut alpha mismatches: ', abad
print '(A,I0)', 'device chain vs serial mismatches: ', bad
print '(A,L1)', 'Are all tests passed? ', bad == 0_I4P .and. abad == 0_I4P
if (bad /= 0_I4P .or. abad /= 0_I4P) error stop 1

contains
   subroutine chain_3d(label, body, m, path, nstep, bad, abad)
   !< The 3-D chain of `body` on m^3 cells over its box grown by 10 %: the check pass, then the timing table.
   character(*),             intent(in)    :: label    !< Body name.
   type(surface_stl_object), intent(in)    :: body     !< Body (body frame).
   integer(I4P),             intent(in)    :: m        !< Cells per axis.
   integer(I4P),             intent(in)    :: path     !< PATH_GAUSS or PATH_TRANSLATION.
   integer(I4P),             intent(in)    :: nstep    !< Time steps.
   integer(I4P),             intent(inout) :: bad      !< Device vs serial mismatches.
   integer(I4P),             intent(inout) :: abad     !< Swept vs cut alpha mismatches.
   type(surface_stl_object)                :: s        !< Moving surface.
   type(rigid_pose_t),       allocatable   :: p(:)     !< Poses P_0..P_nstep.
   type(cut_cells_3d_t)                    :: ch, cd, cr !< Cut cells: CPU, copied back, resident.
   type(cut_links_t)                       :: lh, ld, lr !< Links: CPU, copied back, resident.
   type(swept_cells_3d_t)                  :: wh, wd, wm !< Swept geometry: CPU, device, with the map.
   type(piece_map_3d_t)                    :: mh, md   !< Piece maps: CPU, device.
   real(R8P),                allocatable   :: gx(:), gy(:), gz(:), lo(:,:), hi(:,:), aprev(:,:,:) !< Grid, boxes, alpha.
   integer(I4P),             allocatable   :: th(:), td(:) !< Tags: CPU, device.
   real(R8P)                               :: t(6, NREP, 3) !< Per-step times [stage, repetition, mode].
   integer(I4P)                            :: k, r, v, nb, st, st2, b0, a0 !< Work.
   logical                                 :: tm       !< Time the map.

   tm = path == PATH_GAUSS .and. (full .or. m == sizes3(1))
   call box_grid(body, m, gx, gy, gz)
   call cell_boxes_3d(gx, gy, gz, lo, hi)
   nb = size(lo, 2, kind=I4P) ; allocate(th(nb), td(nb))
   call poses_3d(body, gx(2) - gx(1), path, nstep, p)
   b0 = bad ; a0 = abad
   do k = 1, nstep ! check
      s = body ; call s%move(p(k))
      call s%box_classify_many(lo=lo, hi=hi, codes=th, status=st)
      call s%enter_device
      call s%box_classify_many_device(lo=lo, hi=hi, codes=td, status=st2)
      bad = bad + count(th /= td, kind=I4P) + abs(st - st2) + abs(st)
      call cut_cells_3d(s, gx, gy, gz, ch, status=st, links=lh)
      call cut_cells_3d_device(s, gx, gy, gz, cd, status=st2, links=ld)
      bad = bad + abs(st - st2) + abs(st) + cut_cells_3d_mismatches(ch, cd) + cut_links_mismatches(lh, ld)
      call cut_cells_3d_device(s, gx, gy, gz, cr, status=st2, links=lr, resident=.true.)
      call cr%update_host ; call lr%update_host
      bad = bad + abs(st - st2) + cut_cells_3d_mismatches(ch, cr) + cut_links_mismatches(lh, lr)
      call s%exit_device
      call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, .false., wh, st)
      call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, .true., wd, st2)
      bad = bad + abs(st - st2) + abs(st) + swept_cells_3d_mismatches(wh, wd)
      if (path == PATH_GAUSS) then
         call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, .false., wm, st, mh)
         bad = bad + abs(st) + swept_cells_3d_mismatches(wh, wm)
         call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, .true., wm, st2, md)
         bad = bad + abs(st2) + swept_cells_3d_mismatches(wh, wm) + piece_map_3d_mismatches(mh, md)
      endif
      abad = abad + count(wh%alpha1 /= ch%alpha, kind=I4P)
      if (k > 1) abad = abad + count(wh%alpha0 /= aprev, kind=I4P)
      aprev = ch%alpha
   enddo
   call cr%exit_device ; call lr%exit_device
   print '(A,A,1X,I0,A,A,A,I0,A,I0,A,I0,A,I0,A,A,I0)', '3-D ', trim(label), m, '^3, ', trim(path_name(path)), &
         ': cells ', nb, ', cut ', count(th == BOX_CUT), ', pieces ', ch%np, ', links ', links_number(lh), &
         trim(map_note(path == PATH_GAUSS, mh%nt, mh%evaluated_cells)), ', mismatches ', bad - b0 + abad - a0
   do v = 1, 3
      if (v > 1) call step_3d(body, p, 1_I4P, v, gx, gy, gz, lo, hi, path, tm, s, th, cd, ld, wd, md, t(:, 1, v)) ! warm-up
      call cd%exit_device ; call ld%exit_device
      do r = 1, NREP
         do k = 1, nstep
            call step_3d(body, p, k, v, gx, gy, gz, lo, hi, path, tm, s, th, ch, lh, wh, mh, t(:, r, v))
         enddo
         call s%exit_device
      enddo
      call ch%exit_device ; call lh%exit_device
   enddo
   call print_table(t, nstep, tm)
   endsubroutine chain_3d

   subroutine step_3d(body, p, k, v, gx, gy, gz, lo, hi, path, tm, s, tags, cut, links, sw, map, t)
   !< One timed 3-D step k in mode v (1 CPU, 2 device copied back, 3 device resident); adds the stage times to t
   !< (move, tag, cut + links, swept, the swept call with `map=` if tm, total; reset at k = 1).
   type(surface_stl_object), intent(in)    :: body     !< Body (body frame).
   type(rigid_pose_t),       intent(in)    :: p(0:)    !< Poses.
   integer(I4P),             intent(in)    :: k, v     !< Step, mode.
   real(R8P),                intent(in)    :: gx(:), gy(:), gz(:), lo(:,:), hi(:,:) !< Grid, cell boxes.
   integer(I4P),             intent(in)    :: path     !< Motion.
   logical,                  intent(in)    :: tm       !< Time the swept call with `map=`.
   type(surface_stl_object), intent(inout) :: s        !< Moving surface.
   integer(I4P),             intent(inout) :: tags(:)  !< Tags.
   type(cut_cells_3d_t),     intent(inout) :: cut      !< Cut cells.
   type(cut_links_t),        intent(inout) :: links    !< Links.
   type(swept_cells_3d_t),   intent(inout) :: sw       !< Swept geometry.
   type(piece_map_3d_t),     intent(inout) :: map      !< Piece map.
   real(R8P),                intent(inout) :: t(6)     !< Stage times (summed over the steps).
   integer(I8P)                            :: c(0:5)   !< Clock.

   if (k == 1_I4P) t = 0._R8P
   call tick(c(0))
   if (v > 1) call s%exit_device
   s = body ; call s%move(p(k))
   if (v > 1) call s%enter_device
   call tick(c(1))
   if (v == 1) then
      call s%box_classify_many(lo=lo, hi=hi, codes=tags)
   else
      call s%box_classify_many_device(lo=lo, hi=hi, codes=tags)
   endif
   call tick(c(2))
   select case (v)
   case (1)
      call cut_cells_3d(s, gx, gy, gz, cut, links=links)
   case (2)
      call cut_cells_3d_device(s, gx, gy, gz, cut, links=links)
   case default
      call cut_cells_3d_device(s, gx, gy, gz, cut, links=links, resident=.true.)
   endselect
   call tick(c(3))
   call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, v > 1, sw)
   call tick(c(4))
   if (tm) call swept_3d(body, p(k - 1), p(k), gx, gy, gz, path, v > 1, sw, map=map)
   call tick(c(5))
   call add_times(c, tm, t)
   endsubroutine step_3d

   subroutine swept_3d(body, p0, p1, gx, gy, gz, path, dev, sw, status, map)
   !< The swept geometry of the chain: Gauss nq = 4 (`exact=.false.`) on PATH_GAUSS, the exact translation else; the
   !< piece map of the step if `map` is present.
   type(surface_stl_object), intent(in)            :: body   !< Body (body frame).
   type(rigid_pose_t),       intent(in)            :: p0, p1 !< End poses.
   real(R8P),                intent(in)            :: gx(:), gy(:), gz(:) !< Grid.
   integer(I4P),             intent(in)            :: path   !< Motion.
   logical,                  intent(in)            :: dev    !< On the device.
   type(swept_cells_3d_t),   intent(inout)         :: sw     !< Swept geometry.
   integer(I4P),             intent(out), optional :: status !< Status.
   type(piece_map_3d_t),     intent(inout), optional :: map  !< Piece map.

   if (path == PATH_GAUSS) then
      if (dev) then
         call swept_cells_3d_device(body, p0, p1, gx, gy, gz, sw, nq=4_I4P, exact=.false., status=status, map=map)
      else
         call swept_cells_3d(body, p0, p1, gx, gy, gz, sw, nq=4_I4P, exact=.false., status=status, map=map)
      endif
   else
      if (dev) then
         call swept_cells_3d_device(body, p0, p1, gx, gy, gz, sw, status=status, map=map)
      else
         call swept_cells_3d(body, p0, p1, gx, gy, gz, sw, status=status, map=map)
      endif
   endif
   endsubroutine swept_3d

   subroutine chain_2d(label, body, n, path, nstep, bad, abad)
   !< The 2-D chain of `body` on n^2 cells over [0.0031, 0.9969]^2 grown by 2 %: the check pass, then the timing table.
   character(*),       intent(in)    :: label    !< Body name.
   type(curve_object), intent(in)    :: body     !< Body (body frame).
   integer(I4P),       intent(in)    :: n        !< Cells per axis.
   integer(I4P),       intent(in)    :: path     !< PATH_GAUSS or PATH_TRANSLATION.
   integer(I4P),       intent(in)    :: nstep    !< Time steps.
   integer(I4P),       intent(inout) :: bad      !< Device vs serial mismatches.
   integer(I4P),       intent(inout) :: abad     !< Swept vs cut alpha mismatches.
   type(curve_object)                :: c        !< Moving curve.
   type(rigid_pose_2d_t), allocatable :: p(:)    !< Poses P_0..P_nstep.
   type(cut_cells_2d_t)              :: ch, cr   !< Cut cells: CPU, device (resident, or copied back after it).
   type(cut_links_t)                 :: lh, ld, lr !< Links: CPU, copied back, resident.
   type(swept_cells_2d_t)            :: wh, wd, wm !< Swept geometry: CPU, device, with the map.
   type(piece_map_2d_t)              :: mh, md   !< Piece maps: CPU, device.
   real(R8P),          allocatable   :: x(:), y(:), lo(:,:), hi(:,:), aprev(:,:) !< Grid, boxes, alpha.
   integer(I4P),       allocatable   :: th(:), td(:) !< Tags: CPU, device.
   real(R8P)                         :: t(6, NREP, 3) !< Per-step times [stage, repetition, mode].
   integer(I4P)                      :: k, r, v, nb, st, st2, b0, a0 !< Work.
   logical                           :: tm       !< Time the map.

   tm = path == PATH_GAUSS .and. (full .or. n == sizes2(1))
   call grid_edges([0.0031_R8P, 0.0031_R8P], [0.9969_R8P, 0.9969_R8P], n, x, y)
   call cell_boxes_2d(x, y, lo, hi)
   nb = size(lo, 2, kind=I4P) ; allocate(th(nb), td(nb))
   call poses_2d(x(2) - x(1), path, nstep, p)
   b0 = bad ; a0 = abad
   do k = 1, nstep ! check
      c = body ; call c%move(p(k))
      call c%box_classify_many(lo=lo, hi=hi, codes=th, status=st)
      call c%enter_device
      call c%box_classify_many_device(lo=lo, hi=hi, codes=td, status=st2)
      bad = bad + count(th /= td, kind=I4P) + abs(st - st2) + abs(st)
      call c%cut_cells(x=x, y=y, cut=ch, status=st) ; call cut_links_2d(ch, lh)
      call c%cut_cells_device(x=x, y=y, cut=cr, status=st2, resident=.true.)
      call cut_links_2d_device(cr, ld) ; call cut_links_2d_device(cr, lr, resident=.true.)
      call cr%update_host ; call lr%update_host
      bad = bad + abs(st - st2) + abs(st) + cut_cells_2d_mismatches(ch, cr) + cut_links_mismatches(lh, ld) + &
            cut_links_mismatches(lh, lr)
      call c%exit_device
      call swept_2d(body, p(k - 1), p(k), x, y, path, .false., wh, st)
      call swept_2d(body, p(k - 1), p(k), x, y, path, .true., wd, st2)
      bad = bad + abs(st - st2) + abs(st) + swept_cells_2d_mismatches(wh, wd)
      if (path == PATH_GAUSS) then
         call swept_2d(body, p(k - 1), p(k), x, y, path, .false., wm, st, mh)
         bad = bad + abs(st) + swept_cells_2d_mismatches(wh, wm)
         call swept_2d(body, p(k - 1), p(k), x, y, path, .true., wm, st2, md)
         bad = bad + abs(st2) + swept_cells_2d_mismatches(wh, wm) + piece_map_2d_mismatches(mh, md)
      endif
      abad = abad + count(wh%alpha1 /= ch%alpha, kind=I4P)
      if (k > 1) abad = abad + count(wh%alpha0 /= aprev, kind=I4P)
      aprev = ch%alpha
   enddo
   call cr%exit_device ; call lr%exit_device
   print '(A,A,1X,I0,A,A,A,I0,A,I0,A,I0,A,I0,A,A,I0)', '2-D ', trim(label), n, '^2, ', trim(path_name(path)), &
         ': cells ', nb, ', cut ', count(th == BOX_CUT), ', pieces ', ch%np, ', links ', links_number(lh), &
         trim(map_note(path == PATH_GAUSS, mh%nt, mh%evaluated_cells)), ', mismatches ', bad - b0 + abad - a0
   do v = 1, 3
      if (v > 1) call step_2d(body, p, 1_I4P, v, x, y, lo, hi, path, tm, c, th, cr, ld, wd, md, t(:, 1, v)) ! warm-up
      call ld%exit_device
      do r = 1, NREP
         do k = 1, nstep
            call step_2d(body, p, k, v, x, y, lo, hi, path, tm, c, th, cr, lh, wh, mh, t(:, r, v))
         enddo
         call c%exit_device
      enddo
      call cr%exit_device ; call lh%exit_device
   enddo
   call print_table(t, nstep, tm)
   endsubroutine chain_2d

   subroutine step_2d(body, p, k, v, x, y, lo, hi, path, tm, c, tags, cut, links, sw, map, t)
   !< One timed 2-D step k in mode v (1 CPU, 2 device copied back: the cut resident, the links copied back, then the cut
   !< by `update_host` and released, the twin of the 3-D copied-back path; 3 device resident); adds the stage times to t.
   !< The swept call with `map=` if tm.
   type(curve_object),     intent(in)    :: body     !< Body (body frame).
   type(rigid_pose_2d_t),  intent(in)    :: p(0:)    !< Poses.
   integer(I4P),           intent(in)    :: k, v     !< Step, mode.
   real(R8P),              intent(in)    :: x(:), y(:), lo(:,:), hi(:,:) !< Grid, cell boxes.
   integer(I4P),           intent(in)    :: path     !< Motion.
   logical,                intent(in)    :: tm       !< Time the swept call with `map=`.
   type(curve_object),     intent(inout) :: c        !< Moving curve.
   integer(I4P),           intent(inout) :: tags(:)  !< Tags.
   type(cut_cells_2d_t),   intent(inout) :: cut      !< Cut cells.
   type(cut_links_t),      intent(inout) :: links    !< Links.
   type(swept_cells_2d_t), intent(inout) :: sw       !< Swept geometry.
   type(piece_map_2d_t),   intent(inout) :: map      !< Piece map.
   real(R8P),              intent(inout) :: t(6)     !< Stage times (summed over the steps).
   integer(I8P)                          :: q(0:5)   !< Clock.

   if (k == 1_I4P) t = 0._R8P
   call tick(q(0))
   if (v > 1) call c%exit_device
   c = body ; call c%move(p(k))
   if (v > 1) call c%enter_device
   call tick(q(1))
   if (v == 1) then
      call c%box_classify_many(lo=lo, hi=hi, codes=tags)
   else
      call c%box_classify_many_device(lo=lo, hi=hi, codes=tags)
   endif
   call tick(q(2))
   select case (v)
   case (1)
      call c%cut_cells(x=x, y=y, cut=cut) ; call cut_links_2d(cut, links)
   case (2)
      call c%cut_cells_device(x=x, y=y, cut=cut, resident=.true.) ; call cut_links_2d_device(cut, links)
      call cut%update_host ; call cut%exit_device
   case default
      call c%cut_cells_device(x=x, y=y, cut=cut, resident=.true.) ; call cut_links_2d_device(cut, links, resident=.true.)
   endselect
   call tick(q(3))
   call swept_2d(body, p(k - 1), p(k), x, y, path, v > 1, sw)
   call tick(q(4))
   if (tm) call swept_2d(body, p(k - 1), p(k), x, y, path, v > 1, sw, map=map)
   call tick(q(5))
   call add_times(q, tm, t)
   endsubroutine step_2d

   subroutine swept_2d(body, p0, p1, x, y, path, dev, sw, status, map)
   !< The swept geometry of the chain: Gauss nq = 4 (`exact=.false.`) on PATH_GAUSS, the exact translation else; the
   !< piece map of the step if `map` is present.
   type(curve_object),     intent(in)            :: body   !< Body (body frame).
   type(rigid_pose_2d_t),  intent(in)            :: p0, p1 !< End poses.
   real(R8P),              intent(in)            :: x(:), y(:) !< Grid.
   integer(I4P),           intent(in)            :: path   !< Motion.
   logical,                intent(in)            :: dev    !< On the device.
   type(swept_cells_2d_t), intent(inout)         :: sw     !< Swept geometry.
   integer(I4P),           intent(out), optional :: status !< Status.
   type(piece_map_2d_t),   intent(inout), optional :: map  !< Piece map.

   if (path == PATH_GAUSS) then
      if (dev) then
         call swept_cells_2d_device(body, p0, p1, x, y, sw, nq=4_I4P, exact=.false., status=status, map=map)
      else
         call swept_cells_2d(body, p0, p1, x, y, sw, nq=4_I4P, exact=.false., status=status, map=map)
      endif
   else
      if (dev) then
         call swept_cells_2d_device(body, p0, p1, x, y, sw, status=status, map=map)
      else
         call swept_cells_2d(body, p0, p1, x, y, sw, status=status, map=map)
      endif
   endif
   endsubroutine swept_2d

   subroutine poses_3d(body, h, path, nstep, p)
   !< P_k, k = 0..nstep: PATH_GAUSS rotates by 2k degrees about the axis (0.2, 1, 0.3) through the box centre and
   !< translates by k (0.3, 0.2, -0.1) h; PATH_TRANSLATION translates by k (0.7, 0.4, -0.3) h (h: cell size).
   type(surface_stl_object),        intent(in)  :: body   !< Body.
   real(R8P),                       intent(in)  :: h      !< Cell size.
   integer(I4P),                    intent(in)  :: path   !< Motion.
   integer(I4P),                    intent(in)  :: nstep  !< Steps.
   type(rigid_pose_t), allocatable, intent(out) :: p(:)   !< Poses [0:nstep].
   type(vector_R8P)                             :: b0, b1 !< Box.
   real(R8P)                                    :: cen(3) !< Box centre.
   integer(I4P)                                 :: k      !< Step.

   b0 = body%get_bmin() ; b1 = body%get_bmax()
   cen = 0.5_R8P * [b0%x + b1%x, b0%y + b1%y, b0%z + b1%z]
   allocate(p(0:nstep))
   do k = 0, nstep
      if (path == PATH_GAUSS) then
         p(k) = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], k * 2._R8P * PI / 180._R8P)
         p(k)%t = cen - matmul(p(k)%r, cen) + k * [0.3_R8P, 0.2_R8P, -0.1_R8P] * h
      else
         p(k) = pose_identity()
         p(k)%t = k * [0.7_R8P, 0.4_R8P, -0.3_R8P] * h
      endif
   enddo
   endsubroutine poses_3d

   subroutine poses_2d(h, path, nstep, p)
   !< P_k, k = 0..nstep: PATH_GAUSS rotates by 2k degrees about the block centre (0.5, 0.5) and translates by
   !< k (0.3, 0.2) h; PATH_TRANSLATION translates by k (0.7, 0.4) h (h: cell size).
   real(R8P),                          intent(in)  :: h     !< Cell size.
   integer(I4P),                       intent(in)  :: path  !< Motion.
   integer(I4P),                       intent(in)  :: nstep !< Steps.
   type(rigid_pose_2d_t), allocatable, intent(out) :: p(:)  !< Poses [0:nstep].
   real(R8P), parameter                            :: CEN(2) = [0.5_R8P, 0.5_R8P] !< Block centre.
   integer(I4P)                                    :: k     !< Step.

   allocate(p(0:nstep))
   do k = 0, nstep
      if (path == PATH_GAUSS) then
         p(k) = pose_2d(k * 2._R8P * PI / 180._R8P)
         p(k)%t = CEN - matmul(p(k)%r, CEN) + k * [0.3_R8P, 0.2_R8P] * h
      else
         p(k) = pose_2d(0._R8P, k * [0.7_R8P, 0.4_R8P] * h)
      endif
   enddo
   endsubroutine poses_2d

   subroutine tick(c)
   !< The clock after the device queue drained.
   integer(I8P), intent(out) :: c !< Clock count.

   !$acc wait
   call system_clock(c)
   endsubroutine tick

   subroutine add_times(c, with_map, t)
   !< Add the stage times between the clock counts c(0..5) to t: move, tag, cut + links, swept (c(3)..c(4)), the swept
   !< call with `map=` (c(4)..c(5); 0 without it) and the step (move + tag + cut + links + the swept call with `map=`, or
   !< the swept call without it).
   integer(I8P), intent(in)    :: c(0:5)   !< Clock counts.
   logical,      intent(in)    :: with_map !< The swept call with `map=` was timed.
   real(R8P),    intent(inout) :: t(6)     !< Stage times, total.
   integer(I8P)                :: rate     !< Clock rate.
   real(R8P)                   :: d(5)     !< Intervals.
   integer(I4P)                :: i        !< Counter.

   call system_clock(count_rate=rate)
   do i = 1, 5
      d(i) = real(c(i) - c(i - 1), R8P) / rate
   enddo
   t(1:4) = t(1:4) + d(1:4)
   if (with_map) then
      t(5) = t(5) + d(5)
      t(6) = t(6) + d(1) + d(2) + d(3) + d(5)
   else
      t(6) = t(6) + d(1) + d(2) + d(3) + d(4)
   endif
   endsubroutine add_times

   subroutine print_table(t, nstep, with_map)
   !< The timing table of one case: per mode, per stage the best / median per-step time over the repetitions, the
   !< total and its speed-up against the CPU; the load average first. The map column is the best (median) of the swept
   !< call with `map=` less the best (median) of the call without; it reads `-` without the map.
   real(R8P),    intent(in) :: t(6, NREP, 3) !< Times summed over the steps [stage, repetition, mode].
   integer(I4P), intent(in) :: nstep         !< Steps.
   logical,      intent(in) :: with_map      !< The map was timed.
   real(R8P)                :: w(NREP)       !< Work.
   real(R8P)                :: best(6, 3), med(6, 3) !< Statistics per step.
   character(20)            :: col(6)        !< Columns.
   integer(I4P)             :: v, i          !< Counters.

   do v = 1, 3
      do i = 1, 6
         w = t(i, :, v) / nstep
         best(i, v) = minval(w) ; med(i, v) = median(w)
      enddo
      best(5, v) = best(5, v) - best(4, v) ; med(5, v) = med(5, v) - med(4, v)
   enddo
   call print_load
   print '(3X,A12,6(2X,A20),2X,A)', 'mode', 'move', 'tag', 'cut + links', 'swept', 'map', 'step', 'speed-up'
   do v = 1, 3
      do i = 1, 6
         write(col(i), '(ES9.2,A,ES9.2)') best(i, v), ' /', med(i, v)
      enddo
      if (.not. with_map) col(5) = '                   -'
      print '(3X,A12,6(2X,A20),2X,F7.2,A)', MODES(v), col, best(6, 1) / best(6, v), 'x'
   enddo
   flush(output_unit)
   endsubroutine print_table

   pure function map_note(with_map, nt, ne) result(note)
   !< `, map tracks nt, evaluated cells ne` with the map, else empty.
   logical,      intent(in) :: with_map !< The map was built.
   integer(I4P), intent(in) :: nt, ne   !< Tracks, evaluated cells.
   character(64)            :: note     !< Note.

   note = ''
   if (with_map) write(note, '(A,I0,A,I0)') ', map tracks ', nt, ', evaluated cells ', ne
   endfunction map_note

   subroutine print_load
   !< The load average of the machine (/proc/loadavg), if readable.
   character(64) :: line !< Its line.
   integer       :: u, ios !< Unit, status.

   open(newunit=u, file='/proc/loadavg', action='read', status='old', iostat=ios)
   if (ios /= 0) return
   read(u, '(A)', iostat=ios) line
   close(u)
   if (ios == 0) print '(3X,A,A)', 'load average: ', trim(line)
   endsubroutine print_load

   function median(a) result(m)
   !< The median of a (the lower middle value for an even size).
   real(R8P), intent(in) :: a(:) !< Values.
   real(R8P)             :: m    !< Median.
   real(R8P)             :: b(size(a)), z !< Sorted copy, work.
   integer               :: i, j !< Counters.

   b = a
   do i = 2, size(b)
      z = b(i) ; j = i - 1
      do while (j >= 1)
         if (b(j) <= z) exit
         b(j + 1) = b(j) ; j = j - 1
      enddo
      b(j + 1) = z
   enddo
   m = b((size(b) + 1) / 2)
   endfunction median

   pure function links_number(l) result(n)
   !< The links of l (0 if none).
   type(cut_links_t), intent(in) :: l !< Links.
   integer(I4P)                  :: n !< Their number.

   n = 0_I4P
   if (allocated(l%piece)) n = size(l%piece, kind=I4P)
   endfunction links_number

   pure function path_name(path) result(nm)
   !< The name of a motion.
   integer(I4P), intent(in) :: path !< Motion.
   character(24)            :: nm   !< Name.

   nm = merge('Gauss nq=4 (rotation)   ', 'exact translation       ', path == PATH_GAUSS)
   endfunction path_name

   subroutine box_grid(s, m, gx, gy, gz)
   !< m^3 uniform cells over the box of s grown by 10 % (the grid of the distance bench).
   type(surface_stl_object), intent(in)  :: s              !< Body.
   integer(I4P),             intent(in)  :: m              !< Cells per axis.
   real(R8P), allocatable,   intent(out) :: gx(:), gy(:), gz(:) !< Cell edges.
   type(vector_R8P)                      :: b0, b1         !< Box.
   integer(I4P)                          :: q              !< Counter.

   b0 = s%get_bmin() ; b1 = s%get_bmax()
   gx = [(b0%x - 0.1_R8P * (b1%x - b0%x) + 1.2_R8P * (b1%x - b0%x) * q / m, q = 0, m)]
   gy = [(b0%y - 0.1_R8P * (b1%y - b0%y) + 1.2_R8P * (b1%y - b0%y) * q / m, q = 0, m)]
   gz = [(b0%z - 0.1_R8P * (b1%z - b0%z) + 1.2_R8P * (b1%z - b0%z) * q / m, q = 0, m)]
   endsubroutine box_grid

   subroutine cell_boxes_3d(gx, gy, gz, lo, hi)
   !< One box per cell, cell (i, j, k) at i + (j - 1) nx + (k - 1) nx ny.
   real(R8P),              intent(in)  :: gx(:), gy(:), gz(:) !< Cell edges.
   real(R8P), allocatable, intent(out) :: lo(:,:), hi(:,:)    !< Boxes [3, nx ny nz].
   integer(I4P)                        :: i, j, k, nx, ny, nz, c !< Work.

   nx = size(gx, kind=I4P) - 1_I4P ; ny = size(gy, kind=I4P) - 1_I4P ; nz = size(gz, kind=I4P) - 1_I4P
   allocate(lo(3, nx * ny * nz), hi(3, nx * ny * nz))
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + (j - 1) * nx + (k - 1) * nx * ny
            lo(:, c) = [gx(i), gy(j), gz(k)] ; hi(:, c) = [gx(i + 1), gy(j + 1), gz(k + 1)]
         enddo
      enddo
   enddo
   endsubroutine cell_boxes_3d

   subroutine cell_boxes_2d(x, y, lo, hi)
   !< One box per cell of the grid x, y, cell (i, j) at i + (j - 1) nx.
   real(R8P),              intent(in)  :: x(:), y(:)       !< Cell edges.
   real(R8P), allocatable, intent(out) :: lo(:,:), hi(:,:) !< Boxes [2, nx ny].
   integer(I4P)                        :: i, j, nx, ny     !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P
   allocate(lo(2, nx * ny), hi(2, nx * ny))
   do j = 1, ny
      do i = 1, nx
         lo(:, i + (j - 1) * nx) = [x(i), y(j)] ; hi(:, i + (j - 1) * nx) = [x(i + 1), y(j + 1)]
      enddo
   enddo
   endsubroutine cell_boxes_2d

   subroutine grid_edges(lo, hi, n, x, y)
   !< n x n cells over the box [lo, hi] grown by 2 % (the grid of the curve bench).
   real(R8P),              intent(in)  :: lo(2), hi(2) !< Box.
   integer(I4P),           intent(in)  :: n            !< Cells per axis.
   real(R8P), allocatable, intent(out) :: x(:), y(:)   !< Cell edges.
   real(R8P)                           :: a(2), b(2)   !< Grown box.
   integer(I4P)                        :: k            !< Counter.

   a = lo - 0.02_R8P * (hi - lo) ; b = hi + 0.02_R8P * (hi - lo)
   allocate(x(n + 1), y(n + 1))
   do k = 1, n + 1
      x(k) = a(1) + (b(1) - a(1)) * (k - 1) / real(n, R8P) ; y(k) = a(2) + (b(2) - a(2)) * (k - 1) / real(n, R8P)
   enddo
   endsubroutine grid_edges

   subroutine to_surface_centred(b, s)
   !< The sign-bench body b as a surface, centred on its box (the perf bench's cone).
   type(body_t),             intent(in)    :: b      !< Body.
   type(surface_stl_object), intent(inout) :: s      !< Surface.
   type(vector_R8P)                        :: b0, b1 !< Box.

   call to_surface(b%soup, s)
   b0 = s%get_bmin() ; b1 = s%get_bmax()
   call s%translate(delta=-0.5_R8P * (b0 + b1)) ; call s%analyze
   endsubroutine to_surface_centred

   subroutine naca0012(m, x, y)
   !< NACA 0012 with a closed trailing edge, 2m vertices, cosine spacing, counter-clockwise.
   integer(I4P), intent(in)  :: m          !< Points per side.
   real(R8P),    intent(out) :: x(:), y(:) !< Vertices [2m].
   real(R8P)                 :: xi         !< Abscissa.
   integer(I4P)              :: kk, j      !< Counters.

   j = 0
   do kk = m, 0, -1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = half_thickness(xi)
   enddo
   do kk = 1, m - 1
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1 ; x(j) = xi ; y(j) = -half_thickness(xi)
   enddo
   endsubroutine naca0012

   pure function half_thickness(x) result(y)
   !< NACA 4-digit half thickness, 12 %, closed trailing edge.
   real(R8P), intent(in) :: x !< Abscissa.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction half_thickness
endprogram fossil_val_device_chain
