!< FOSSIL test: the piece map of a time step in 3-D, the tracks of every cell (issue #102 step 2).
program fossil_test_piece_map_3d
!< FOSSIL test: the piece map of a time step in 3-D, the tracks of every cell (issue #102 step 2).
!<
!< Hand-built cases on the 3 x 3 x 3 block [-1, 2]^3, each a surface and a motion whose tracks in the middle cell
!< [0, 1]^3 are known (the end-pose overlap, which links an old and a new piece iff they overlap, is given for
!< contrast). The bodies are boxes, tilted by a small rotation so that no face of them lies on a grid plane:
!< - a thin plate crossing the whole cell within the step (translation; the Gauss path; rotation about a far pivot,
!<   exact and Gauss): one dead and one fresh track (overlap: continuing);
!< - a plate 2e-11 thick crossing the cell (issue #104): each corner covered for 1e-11 of the step, between two events
!<   closer than the 1e-10 under which they used to act as one (continuing), far apart against their round-off: the
!<   state is evaluated (`narrow_cells`): dead + fresh, exact and Gauss;
!< - a narrow slot between two slabs moving by more than its width inside the cell: continuing (overlap: dead + fresh);
!< - a plate entering from the top and crossing the cell: split; leaving: merge;
!< - a block uncovering the cell: fresh (cleared); covering it: dead (covered);
!< - a hollow cube inside the cell, its cavity a fluid piece bounded by facets only, translating and rotating (exact
!<   and Gauss): two continuing tracks;
!< - a fluid groove (V notch) whose apex line enters the cell through a face: continuing (its facets, its only
!<   features, persist through the events);
!< - a slab sweeping over the cell corner (1, 1, 1) and stopping past it: the main piece continues and the corner piece
!<   is fresh (overlap: a split);
!< - tangential contacts under a rotation (`contact_cells`): the corner (1, 1, 1) grazing a facet at the turning point
!<   of the rotation, and the apex of a fluid pit grazing the face x = 1 from inside the cell. The cell's piece shrinks
!<   towards a point and grows back. Exactly at the contact (delta = 0) the event is not certified by a sign change: the
!<   contact is within the round-off of its function, so whether the piece vanishes (dead + fresh) or not (continuing)
!<   is undecidable in double precision: untrusted (it read continuing before). The corner 1e-9 inside the facet's
!<   solid at the turning point: dead + fresh; 1e-9 short of it: continuing.
!< Every case also passes the map checks (`check_map`): every piece at pose0 and at pose1 in exactly one track of its
!< cell, numbered as the cut cells of a copy moved to the pose (the track volumes equal to their pieces' volumes bit for
!< bit), per-cell sums over the tracks equal to the swept geometry (faces, wall, residual, swept = V1 - V0, outflows =
!< residual - swept), residual = V1 - V0 + sum of the outflows per track, every per-track residual at round-off on the
!< exact paths (`fossil_piece_map_3d_checks`). The cone on a grid and the refusals: `fossil_test_piece_map_3d_bodies`.
use fossil,                   only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, &
                                     swept_cells_3d, piece_map_3d_t, cut_cells_3d_t, cut_cells_3d, STATUS_OK, &
                                     TRACK_CONTINUING, TRACK_FRESH, TRACK_DEAD, TRACK_SPLIT, TRACK_MERGE, &
                                     TRACK_UNTRUSTED
use fossil_piece_map_3d_checks, only : check_map
use fossil_sign_bench_bodies, only : soup_t, body_t, make_prism, to_surface
use penf,                     only : I4P, R8P

implicit none

type(surface_stl_object)        :: s
type(swept_cells_3d_t)          :: sw
type(piece_map_3d_t)            :: map
type(soup_t)                    :: sp
type(body_t)                    :: gb
type(rigid_pose_t)              :: p0, p1
real(R8P)                       :: x(4), n(3), c(3), ax(3)
real(R8P), parameter            :: TA(3) = [0.3_R8P, 0.5_R8P, 0.8_R8P], TANG = 0.02_R8P !< Tilt of the box bodies.
logical                         :: ok(20)
integer(I4P)                    :: nput !< Facets written by `put`.

ok = .true.
x = [-1._R8P, 0._R8P, 1._R8P, 2._R8P]
! 1-4. thin plate crossing the cell within the step: dead + fresh
call boxes(s, reshape([-3._R8P, -3._R8P, -0.6_R8P, 4._R8P, 4._R8P, -0.5_R8P], [3, 2, 1]))
ok(1) = case_3d('plate crossing, translation', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), .true., &
                [TRACK_DEAD, TRACK_FRESH])
ok(2) = case_3d('plate crossing, translation, Gauss', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), &
                .false., [TRACK_DEAD, TRACK_FRESH])
! 19-20. a plate 2e-11 thick crossing the cell: each corner covered for 1e-11 of the step (issue #104): dead + fresh
call boxes(s, reshape([-3._R8P, -3._R8P, -0.55_R8P, 4._R8P, 4._R8P, -0.55_R8P + 2e-11_R8P], [3, 2, 1]))
ok(19) = case_3d('plate 2e-11 thick crossing, translation', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), &
                 .true., [TRACK_DEAD, TRACK_FRESH]) .and. map%narrow_cells >= 1_I4P
ok(20) = case_3d('plate 2e-11 thick crossing, translation, Gauss', tilt([0._R8P, 0._R8P, 0._R8P]), &
                 tilt([0._R8P, 0._R8P, 2._R8P]), .false., [TRACK_DEAD, TRACK_FRESH]) .and. map%narrow_cells >= 1_I4P
print '(A,I0)', '  narrow cells (state between events 1e-11 apart evaluated): ', map%narrow_cells
call boxes(s, reshape([-0.05_R8P, 3._R8P, -3._R8P, 0.05_R8P, 8._R8P, 4._R8P], [3, 2, 1]))
ax = [0.02_R8P, 0.01_R8P, 1._R8P]
ok(3) = case_3d('plate crossing, rotation', pose_from_axis_angle(ax, 0.2_R8P, [0.5_R8P, -5._R8P, 0._R8P]), &
                pose_from_axis_angle(ax, -0.2_R8P, [0.5_R8P, -5._R8P, 0._R8P]), .true., [TRACK_DEAD, TRACK_FRESH])
ok(4) = case_3d('plate crossing, rotation, Gauss', pose_from_axis_angle(ax, 0.2_R8P, [0.5_R8P, -5._R8P, 0._R8P]), &
                pose_from_axis_angle(ax, -0.2_R8P, [0.5_R8P, -5._R8P, 0._R8P]), .false., [TRACK_DEAD, TRACK_FRESH])
! 5. a slot of width 0.05 between two slabs moving by 0.2 inside the cell: continuing
call boxes(s, reshape([-3._R8P, -3._R8P, -3._R8P, 4._R8P, 0.40_R8P, 4._R8P, -3._R8P, 0.45_R8P, -3._R8P, 4._R8P, 4._R8P, &
                       4._R8P], [3, 2, 2]))
ok(5) = case_3d('narrow slot moving past its width', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0.6_R8P, 0.2_R8P, -0.5_R8P]), .true., &
                [TRACK_CONTINUING])
! 6-7. a plate entering from the top and crossing the cell: split; the reverse motion: merge
call boxes(s, reshape([-3._R8P, 0.45_R8P, 1.2_R8P, 4._R8P, 0.55_R8P, 5._R8P], [3, 2, 1]))
ok(6) = case_3d('plate entering: split', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, -1.6_R8P]), .true., &
                [TRACK_SPLIT])
ok(7) = case_3d('plate leaving: merge', tilt([0._R8P, 0._R8P, -1.6_R8P]), tilt([0._R8P, 0._R8P, 0._R8P]), .true., &
                [TRACK_MERGE])
! 8-9. a block uncovering the cell: cleared (fresh); covering it: covered (dead)
call boxes(s, reshape([-3._R8P, -3._R8P, -3._R8P, 1.2_R8P, 4._R8P, 4._R8P], [3, 2, 1]))
ok(8) = case_3d('block uncovering: cleared', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([-1.5_R8P, 0._R8P, 0._R8P]), .true., &
                [TRACK_FRESH])
ok(9) = case_3d('block covering: covered', tilt([-1.5_R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 0._R8P]), .true., &
                [TRACK_DEAD])
ok(9) = ok(9) .and. map%tracks(2, 2, 2) == 1_I4P
! 10-12. a hollow cube inside the cell, its cavity a fluid piece: two continuing tracks
call boxes(s, reshape([0.2_R8P, 0.2_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, 0.8_R8P, 0.35_R8P, 0.35_R8P, 0.35_R8P, 0.65_R8P, &
                       0.65_R8P, 0.65_R8P], [3, 2, 2]), inner=2_I4P)
ok(10) = case_3d('cavity, translation', tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0.05_R8P, 0.03_R8P, 0.02_R8P]), .true., &
                 [TRACK_CONTINUING, TRACK_CONTINUING])
c = [0.5_R8P, 0.5_R8P, 0.5_R8P]
p1 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.3_R8P) ; p1%t = c - matmul(p1%r, c) + [0.04_R8P, -0.03_R8P, 0.02_R8P]
ok(11) = case_3d('cavity, rotation', tilt([0._R8P, 0._R8P, 0._R8P]), p1, .true., [TRACK_CONTINUING, TRACK_CONTINUING])
ok(12) = case_3d('cavity, rotation, Gauss', tilt([0._R8P, 0._R8P, 0._R8P]), p1, .false., &
                 [TRACK_CONTINUING, TRACK_CONTINUING])
! 13. a slab normal to (1, 1, 1) sweeping over the corner (1, 1, 1) and stopping past it: continuing + fresh
call boxes(s, reshape([-0.025_R8P, -5._R8P, -5._R8P, 0.025_R8P, 5._R8P, 5._R8P], [3, 2, 1]))
n = [1._R8P, 1._R8P, 1._R8P] / sqrt(3._R8P)
p0 = pose_from_axis_angle([0._R8P, -1._R8P, 1._R8P], acos(n(1)), [1._R8P, 1._R8P, 1._R8P] + 0.3_R8P * n)
p1 = p0 ; p1%t = [1._R8P, 1._R8P, 1._R8P] - 0.15_R8P * n
ok(13) = case_3d('slab over a corner', p0, p1, .true., [TRACK_CONTINUING, TRACK_FRESH])
! 14. a fluid groove (a V notch along z) whose apex line enters the cell through its -x face: continuing (the
! groove's facets are its only features and persist through the events)
call make_prism(reshape([-2._R8P, -2._R8P, 3._R8P, -2._R8P, 3._R8P, 0.19_R8P, -0.1_R8P, 0.5_R8P, 3._R8P, 0.81_R8P, 3._R8P, &
                         3._R8P, -2._R8P, 3._R8P], [2, 7]), 7._R8P, 4, 'groove', gb)
call to_surface(gb%soup, s)
ok(14) = case_3d('groove apex entering', tilt([0._R8P, 0._R8P, -3._R8P]), tilt([0.3_R8P, 0._R8P, -3._R8P]), .true., &
                 [TRACK_CONTINUING])

! 15-17. the corner (1, 1, 1) grazing a facet at the turning point of the rotation (delta: the gap left at the contact)
call grazing_corner(0._R8P)
ok(15) = case_3d('corner grazing, delta 0', graze_pose(-1._R8P), graze_pose(1._R8P), .true., [TRACK_UNTRUSTED]) .and. &
         map%contact_cells >= 1_I4P
call grazing_corner(-1e-9_R8P)
ok(16) = case_3d('corner grazing, delta -1e-9', graze_pose(-1._R8P), graze_pose(1._R8P), .true., [TRACK_DEAD, TRACK_FRESH])
call grazing_corner(1e-9_R8P)
ok(17) = case_3d('corner grazing, delta 1e-9', graze_pose(-1._R8P), graze_pose(1._R8P), .true., [TRACK_CONTINUING])
! 18. the apex of a fluid pit grazing the face x = 1 from inside the cell
call pit(s)
p0 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], -0.15_R8P, [-1._R8P, 0.5_R8P, 0.5_R8P])
p1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.15_R8P, [-1._R8P, 0.5_R8P, 0.5_R8P])
ok(18) = case_3d('pit apex grazing a face', p0, p1, .true., [TRACK_UNTRUSTED]) .and. map%contact_cells >= 1_I4P


print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine grazing_corner(delta)
   !< A big block whose face, rotating about a pivot beyond the corner (1, 1, 1) on the diagonal (`graze_pose`), comes
   !< closest to the corner at the middle of the step, delta short of it (the corner on the face for delta = 0): the
   !< cell's only fluid is the corner tetrahedron.
   real(R8P), intent(in) :: delta !< Gap left at the contact (< 0: the corner inside the solid).

   call boxes(s, reshape([-22._R8P - delta, -10._R8P, -12._R8P, -2._R8P - delta, 13._R8P, 10._R8P], [3, 2, 1]))
   endsubroutine grazing_corner

   function graze_pose(sgn) result(p)
   !< The end poses of `grazing_corner`: the body origin on the pivot (1, 1, 1) + 2 (1, 1, 1) / sqrt(3), its +x axis
   !< along the diagonal turned by sgn 0.3 about (1, -1, 0).
   real(R8P), intent(in) :: sgn  !< -1 (pose0) or 1 (pose1).
   type(rigid_pose_t)    :: p    !< Pose.
   type(rigid_pose_t)    :: r, a !< Work.
   real(R8P)             :: d(3) !< Diagonal.

   d = [1._R8P, 1._R8P, 1._R8P] / sqrt(3._R8P)
   r = pose_from_axis_angle([0._R8P, -1._R8P, 1._R8P], acos(d(1)))
   a = pose_from_axis_angle([1._R8P, -1._R8P, 0._R8P], sgn * 0.3_R8P)
   p%r = matmul(a%r, r%r) ; p%t = [1._R8P, 1._R8P, 1._R8P] + 2._R8P * d
   endfunction graze_pose

   subroutine pit(surf)
   !< A block with a square pyramidal pit (fluid) on its +x face, the apex at (2, 0, 0): the block [-8, 3] x [-8, 8]^2,
   !< the pit's opening [-0.4, 0.4]^2 on x = 3.
   type(surface_stl_object), intent(inout) :: surf    !< Surface.
   real(R8P)                               :: q(3, 8) !< The +x face's outer and inner rings.

   if (allocated(sp%v)) deallocate(sp%v)
   allocate(sp%v(3, 3, 22)) ; nput = 0_I4P
   call put4([-8._R8P, -8._R8P, -8._R8P], [-8._R8P, -8._R8P, 8._R8P], [-8._R8P, 8._R8P, 8._R8P], [-8._R8P, 8._R8P, -8._R8P])
   call put4([-8._R8P, -8._R8P, -8._R8P], [3._R8P, -8._R8P, -8._R8P], [3._R8P, -8._R8P, 8._R8P], [-8._R8P, -8._R8P, 8._R8P])
   call put4([-8._R8P, 8._R8P, -8._R8P], [-8._R8P, 8._R8P, 8._R8P], [3._R8P, 8._R8P, 8._R8P], [3._R8P, 8._R8P, -8._R8P])
   call put4([-8._R8P, -8._R8P, -8._R8P], [-8._R8P, 8._R8P, -8._R8P], [3._R8P, 8._R8P, -8._R8P], [3._R8P, -8._R8P, -8._R8P])
   call put4([-8._R8P, -8._R8P, 8._R8P], [3._R8P, -8._R8P, 8._R8P], [3._R8P, 8._R8P, 8._R8P], [-8._R8P, 8._R8P, 8._R8P])
   q(:, 1) = [3._R8P, -8._R8P, -8._R8P] ; q(:, 2) = [3._R8P, 8._R8P, -8._R8P] ; q(:, 3) = [3._R8P, 8._R8P, 8._R8P]
   q(:, 4) = [3._R8P, -8._R8P, 8._R8P] ; q(:, 5) = [3._R8P, -0.4_R8P, -0.4_R8P] ; q(:, 6) = [3._R8P, 0.4_R8P, -0.4_R8P]
   q(:, 7) = [3._R8P, 0.4_R8P, 0.4_R8P] ; q(:, 8) = [3._R8P, -0.4_R8P, 0.4_R8P]
   call put4(q(:, 1), q(:, 2), q(:, 6), q(:, 5)) ; call put4(q(:, 2), q(:, 3), q(:, 7), q(:, 6))
   call put4(q(:, 3), q(:, 4), q(:, 8), q(:, 7)) ; call put4(q(:, 4), q(:, 1), q(:, 5), q(:, 8))
   call put(q(:, 5), q(:, 6), [2._R8P, 0._R8P, 0._R8P]) ; call put(q(:, 6), q(:, 7), [2._R8P, 0._R8P, 0._R8P])
   call put(q(:, 7), q(:, 8), [2._R8P, 0._R8P, 0._R8P]) ; call put(q(:, 8), q(:, 5), [2._R8P, 0._R8P, 0._R8P])
   call to_surface(sp, surf)
   endsubroutine pit

   subroutine put4(a1, a2, a3, a4)
   !< The quad a1 a2 a3 a4 (outward by the right-hand rule) as two facets of sp.
   real(R8P), intent(in) :: a1(3), a2(3), a3(3), a4(3) !< Vertices.

   call put(a1, a2, a3) ; call put(a1, a3, a4)
   endsubroutine put4

   subroutine put(a1, a2, a3)
   !< One facet of sp.
   real(R8P), intent(in) :: a1(3), a2(3), a3(3) !< Vertices.

   nput = nput + 1_I4P ; sp%v(:, 1, nput) = a1 ; sp%v(:, 2, nput) = a2 ; sp%v(:, 3, nput) = a3
   endsubroutine put

   function tilt(t) result(p)
   !< The tilt of the box bodies (a small rotation about the origin) with translation t.
   real(R8P), intent(in) :: t(3) !< Translation.
   type(rigid_pose_t)    :: p    !< Pose.

   p = pose_from_axis_angle(TA, TANG, t)
   endfunction tilt

   subroutine boxes(surf, b, inner)
   !< A surface of the boxes b(:, 1, k) .. b(:, 2, k), outward; box `inner` (if given) inward (a cavity).
   type(surface_stl_object), intent(inout)        :: surf   !< Surface.
   real(R8P),                intent(in)           :: b(:,:,:) !< Boxes [3, 2, nb] (or [3, 2] through the caller's reshape).
   integer(I4P),             intent(in), optional :: inner  !< Box with inward faces.
   integer(I4P), parameter :: QD(4, 6) = reshape([0, 4, 6, 2, 1, 3, 7, 5, 0, 1, 5, 4, 2, 6, 7, 3, 0, 2, 3, 1, 4, 5, 7, 6], &
                                                 [4, 6]) !< Face quads (corner bits x + 2 y + 4 z), outward.
   real(R8P)    :: v(3, 0:7) !< Corners.
   integer(I4P) :: k, f, q, m !< Work.

   if (allocated(sp%v)) deallocate(sp%v)
   allocate(sp%v(3, 3, 12 * size(b, 3)))
   m = 0
   do k = 1, size(b, 3)
      do q = 0, 7
         v(:, q) = [b(1, 1 + mod(q, 2), k), b(2, 1 + mod(q / 2, 2), k), b(3, 1 + q / 4, k)]
      enddo
      do f = 1, 6
         if (present(inner)) then
            if (inner == k) then
               m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(3, f)), v(:, QD(2, f))], [3, 3])
               m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(4, f)), v(:, QD(3, f))], [3, 3])
               cycle
            endif
         endif
         m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(2, f)), v(:, QD(3, f))], [3, 3])
         m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(3, f)), v(:, QD(4, f))], [3, 3])
      enddo
   enddo
   call to_surface(sp, surf)
   endsubroutine boxes

   logical function case_3d(name, q0, q1, exact, expected) result(pass)
   !< One hand-built case: the statuses of the tracks of the middle cell, in track order, and the map checks.
   character(*),       intent(in) :: name        !< Case.
   type(rigid_pose_t), intent(in) :: q0, q1      !< End poses.
   logical,            intent(in) :: exact       !< Exact paths.
   integer(I4P),       intent(in) :: expected(:) !< Statuses of the middle cell's tracks.
   integer(I4P), allocatable      :: got(:)      !< Statuses.
   character(64)                  :: ov          !< End-pose overlap verdict.
   integer(I4P)                   :: st_, t0     !< Work.
   logical                        :: checks      !< Map checks.

   call swept_cells_3d(s, q0, q1, x, x, x, sw, exact=exact, status=st_, map=map)
   t0 = map%first_track(2, 2, 2)
   got = map%status(t0:t0 + map%tracks(2, 2, 2) - 1_I4P)
   call overlap(q0, q1, ov)
   checks = check_map(s, sw, map, q0, q1, x, exact)
   pass = st_ == STATUS_OK .and. size(got) == size(expected) .and. checks
   if (pass) pass = all(got == expected)
   print '(A,A,A,*(I3))', '  ', name, ': middle cell tracks', got
   print '(A,A,A,I0,A,I0,A,L1)', '    (end-pose overlap: ', trim(ov), ')  untrusted cells ', map%untrusted_cells, &
         '  evaluated cells ', map%evaluated_cells, '  ok ', pass
   endfunction case_3d


   subroutine overlap(q0, q1, verdict)
   !< The end-pose overlap of the middle cell: one old and one new piece overlap (alpha of the intersection > 0)?
   !< Approximated by both ends having fluid with the same piece counts (the overlap's answer for these cases).
   type(rigid_pose_t), intent(in)  :: q0, q1  !< End poses.
   character(*),       intent(out) :: verdict !< Verdict.
   type(surface_stl_object)        :: c0, c1  !< Copies.
   type(cut_cells_3d_t)            :: k0, k1  !< Cuts.

   c0 = s ; call c0%move(q0) ; call cut_cells_3d(c0, x, x, x, k0)
   c1 = s ; call c1%move(q1) ; call cut_cells_3d(c1, x, x, x, k1)
   write(verdict, '(A,I0,A,I0,A)') 'pieces ', k0%pieces(2, 2, 2), ' -> ', k1%pieces(2, 2, 2), ' at the ends'
   endsubroutine overlap

endprogram fossil_test_piece_map_3d
