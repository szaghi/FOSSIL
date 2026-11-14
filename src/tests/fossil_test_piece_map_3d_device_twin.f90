!< FOSSIL test: the 3-D piece map of `swept_cells_3d_device` equals that of `swept_cells_3d` bit for bit (issue #105 step 3).
program fossil_test_piece_map_3d_device_twin
!< FOSSIL test: the 3-D piece map of `swept_cells_3d_device` equals that of `swept_cells_3d` bit for bit (issue #105 step 3).
!<
!< Built without OpenACC (this CI build) the device loops of the map (`piece_map_3d_build(device=)`: the evaluations'
!< sizes, both evaluation passes and the fill, the tracks per cell, the nodes of every measure round with their retries
!< and sums, the round's test per pair, the assembly) run on the host: the test catches the driver bugs of the device
!< twin (the flat arrays handed to the launch sets, the padded slices, the retry lists, the launch sets themselves, the
!< per-cell sums), not device codegen or residency (the distance bench does those on the GPU, `piece map 3-D device vs
!< serial mismatches`). Every map compared by its bits (`piece_map_3d_mismatches`), the swept geometry and the status
!< too:
!< - the cone of the sign bench on 3^3 cells: exact translation, exact rotation (doubling rounds, posed evaluations);
!< - on the 3 x 3 x 3 block [-1, 2]^3 (the hand cases of `fossil_test_piece_map_3d`): a thin plate crossing the middle
!<   cell (dead + fresh: multi-track cells measured per event interval), translating and by Gauss; a hollow cube in the
!<   cell, its cavity a second fluid piece (end evaluations with piece volumes), rotating, exact and Gauss + refine;
!<   the plate's translation and both cavity cases again in launch sets of 4 KB of workspace (`chunk_bytes`: many sets
!<   per pass and per round);
!< - the refusal of a map with node poses: the same statuses;
!< - negative control: the plate's end pose one ulp off gives mismatches.
!< (The retries with the exact corner states, of an evaluation and of a node, are not reached here, nor by the host
!< map's tests.)
use fossil,                   only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, pose_interpolate, swept_cells_3d_t, &
                                     swept_cells_3d, swept_cells_3d_device, piece_map_3d_t, gauss_legendre_01
use fossil_piece_map_compare, only : piece_map_3d_mismatches
use fossil_swept_3d_compare,  only : swept_cells_3d_mismatches
use fossil_sign_bench_bodies, only : soup_t, body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P

implicit none

real(R8P), parameter            :: TA(3) = [0.3_R8P, 0.5_R8P, 0.8_R8P], TANG = 0.02_R8P !< Tilt of the box bodies.
type(surface_stl_object)        :: s           !< Body.
type(body_t), allocatable       :: bb(:)       !< Bench bodies.
type(soup_t)                    :: sp          !< Box soup.
type(rigid_pose_t)              :: p0, p1, pr  !< Poses.
type(rigid_pose_t), allocatable :: pq(:)       !< Node poses.
type(piece_map_3d_t)            :: a, b        !< Maps (negative control, refusal).
type(swept_cells_3d_t)          :: sa, sb      !< Swept geometry (negative control, refusal).
real(R8P),          allocatable :: xg(:), sq(:), wq(:) !< Grid, Gauss rule.
real(R8P)                       :: x3(4), c(3) !< The 3 x 3 x 3 block [-1, 2]^3, a centre.
integer(I4P)                    :: i, k, m, total !< Work.
logical                         :: ok          !< All cases.

ok = .true. ; total = 0_I4P
! the cone on 3^3 cells
call bench_bodies(bb) ; call to_surface(bb(7)%soup, s)
xg = [(-0.3_R8P + 0.6_R8P * i, i = 0, 3)]
p0 = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.4_R8P, [0.31_R8P, 0.27_R8P, 0.12_R8P])
p1 = p0 ; p1%t = p0%t + [0.37_R8P, 0.21_R8P, 0.29_R8P]
pr = pose_from_axis_angle([0.2_R8P, -0.3_R8P, 1._R8P], 0.05_R8P) ; pr%r = matmul(pr%r, p0%r)
pr%t = p0%t + [0.1_R8P, 0.05_R8P, 0.08_R8P]
call check('cone, translation', twin(xg, p0, p1, .true.))
call check('cone, rotation', twin(xg, p0, pr, .true.))
! hand cases on [-1, 2]^3
x3 = [-1._R8P, 0._R8P, 1._R8P, 2._R8P]
call boxes(reshape([-3._R8P, -3._R8P, -0.6_R8P, 4._R8P, 4._R8P, -0.5_R8P], [3, 2, 1]))
call check('plate crossing, translation', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), .true.))
call check('plate crossing, Gauss', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), .false.))
call check('plate crossing, translation, 4 KB sets', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), tilt([0._R8P, 0._R8P, 2._R8P]), &
                                                          .true., chunk=4096_I8P))
call boxes(reshape([0.2_R8P, 0.2_R8P, 0.2_R8P, 0.8_R8P, 0.8_R8P, 0.8_R8P, 0.35_R8P, 0.35_R8P, 0.35_R8P, 0.65_R8P, 0.65_R8P, &
                    0.65_R8P], [3, 2, 2]), inner=2_I4P)
c = [0.5_R8P, 0.5_R8P, 0.5_R8P]
p1 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.3_R8P) ; p1%t = c - matmul(p1%r, c) + [0.04_R8P, -0.03_R8P, 0.02_R8P]
call check('cavity, rotation', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), p1, .true.))
call check('cavity, Gauss + refine', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), p1, .false., refine=.true.))
call check('cavity, rotation, 4 KB sets', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), p1, .true., chunk=4096_I8P))
call check('cavity, Gauss + refine, 4 KB sets', twin(x3, tilt([0._R8P, 0._R8P, 0._R8P]), p1, .false., refine=.true., &
                                                     chunk=4096_I8P))
! refusal: a map with node poses
call gauss_legendre_01(2_I4P, sq, wq)
allocate(pq(2))
do i = 1, 2
   pq(i) = pose_interpolate(p0, p1, sq(i))
enddo
call swept_cells_3d(s, p0, p1, x3, x3, x3, sa, nq=2_I4P, poses=pq, status=k, map=a)
call swept_cells_3d_device(s, p0, p1, x3, x3, x3, sb, nq=2_I4P, poses=pq, status=m, map=b)
call check('refusal: map with node poses', abs(k - m) + merge(0_I4P, 1_I4P, k /= 0_I4P))
! negative control
call boxes(reshape([-3._R8P, -3._R8P, -0.6_R8P, 4._R8P, 4._R8P, -0.5_R8P], [3, 2, 1]))
p0 = tilt([0._R8P, 0._R8P, 0._R8P]) ; p1 = tilt([0._R8P, 0._R8P, 2._R8P])
call swept_cells_3d(s, p0, p1, x3, x3, x3, sa, map=a)
p1%t(3) = nearest(p1%t(3), 1._R8P)
call swept_cells_3d_device(s, p0, p1, x3, x3, x3, sb, map=b)
m = piece_map_3d_mismatches(a, b)
print '(A,I0)', 'negative control (one ulp): mismatches ', m
ok = ok .and. m > 0_I4P
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   function twin(g, q0, q1, exact, refine, chunk) result(n)
   !< Mismatches of the device twin against the host call (same arguments, nq = 4, the grid g^3): map, swept geometry,
   !< status.
   real(R8P),          intent(in)           :: g(:)   !< Cell edges (every direction).
   type(rigid_pose_t), intent(in)           :: q0, q1 !< End poses.
   logical,            intent(in)           :: exact  !< Exact paths.
   logical,            intent(in), optional :: refine !< Gauss path refined per event interval.
   integer(I8P),       intent(in), optional :: chunk  !< Workspace bytes per launch set.
   integer(I4P)                             :: n      !< Mismatches.
   type(piece_map_3d_t)                     :: mh, md !< Maps.
   type(swept_cells_3d_t)                   :: wh, wd !< Swept geometry.
   integer(I4P)                             :: sh, sd !< Statuses.

   call swept_cells_3d(s, q0, q1, g, g, g, wh, nq=4_I4P, exact=exact, status=sh, map=mh, refine=refine)
   call swept_cells_3d_device(s, q0, q1, g, g, g, wd, nq=4_I4P, exact=exact, status=sd, map=md, refine=refine, &
                              chunk_bytes=chunk)
   n = piece_map_3d_mismatches(mh, md, verbose=.true.) + swept_cells_3d_mismatches(wh, wd) + abs(sh - sd) + abs(sh)
   if (sh == 0_I4P .and. mh%nt > 0_I4P) print '(A,I0,A,I0,A,I0,A,I0)', '   tracks ', mh%nt, ', evaluated cells ', &
      mh%evaluated_cells, ', multi-track cells ', mh%multi_track_cells, ', evaluations ', mh%evaluations
   endfunction twin

   subroutine check(label, n)
   !< One case.
   character(*), intent(in) :: label !< Case.
   integer(I4P), intent(in) :: n     !< Mismatches.

   print '(A,I0)', label//': mismatches ', n
   total = total + n ; ok = ok .and. n == 0_I4P
   endsubroutine check

   function tilt(t) result(p)
   !< The tilt of the box bodies (a small rotation about the origin) with translation t.
   real(R8P), intent(in) :: t(3) !< Translation.
   type(rigid_pose_t)    :: p    !< Pose.

   p = pose_from_axis_angle(TA, TANG, t)
   endfunction tilt

   subroutine boxes(bx, inner)
   !< The surface s of the boxes bx(:, 1, k) .. bx(:, 2, k), outward; box `inner` (if given) inward (a cavity).
   real(R8P),    intent(in)           :: bx(:,:,:) !< Boxes [3, 2, nb].
   integer(I4P), intent(in), optional :: inner     !< Box with inward faces.
   integer(I4P), parameter :: QD(4, 6) = reshape([0, 4, 6, 2, 1, 3, 7, 5, 0, 1, 5, 4, 2, 6, 7, 3, 0, 2, 3, 1, 4, 5, 7, 6], &
                                                 [4, 6]) !< Face quads (corner bits x + 2 y + 4 z), outward.
   real(R8P)    :: v(3, 0:7)   !< Corners.
   integer(I4P) :: k, f, q, m  !< Work.
   logical      :: inw         !< Inward box.

   if (allocated(sp%v)) deallocate(sp%v)
   allocate(sp%v(3, 3, 12 * size(bx, 3)))
   m = 0
   do k = 1, size(bx, 3)
      do q = 0, 7
         v(:, q) = [bx(1, 1 + mod(q, 2), k), bx(2, 1 + mod(q / 2, 2), k), bx(3, 1 + q / 4, k)]
      enddo
      inw = .false. ; if (present(inner)) inw = inner == k
      do f = 1, 6
         if (inw) then
            m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(3, f)), v(:, QD(2, f))], [3, 3])
            m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(4, f)), v(:, QD(3, f))], [3, 3])
         else
            m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(2, f)), v(:, QD(3, f))], [3, 3])
            m = m + 1 ; sp%v(:, :, m) = reshape([v(:, QD(1, f)), v(:, QD(3, f)), v(:, QD(4, f))], [3, 3])
         endif
      enddo
   enddo
   call to_surface(sp, s)
   endsubroutine boxes
endprogram fossil_test_piece_map_3d_device_twin
