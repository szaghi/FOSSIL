!< FOSSIL validation: the 3-D piece map of a time step on bodies (issue #102 step 2).
program fossil_val_piece_map_3d
!< FOSSIL validation: the 3-D piece map of a time step on bodies (issue #102 step 2).
!<
!< The marching-cubes sphere (24^3 samples) and the cone of the sign bench (32 faces, centred) on n^3 cells over a box
!< around them, one step of the moving-body motions of the step-0 measurement: a translation along (0.7, 0.4, -0.3),
!< and a rotation about (0.2, 1, 0.3) through the box centre plus a translation along (0.3, 0.2, -0.1), scaled so that
!< the largest vertex displacement is 2 cells. Checks (`fossil_piece_map_3d_checks`) at 16^3 and 24^3: translation
!< (exact), rotation (exact, 16^3; Gauss, nq = 4): numbering as explicit cuts of moved copies, per-cell sums equal to
!< the swept geometry, the residual identity, per-track residuals at round-off on the exact paths; the untrusted
!< cells and the track statuses are counted. Then the cost: the time of `swept_cells_3d` with the map against without,
!< at 32^3 (and 64^3 with the argument `full`), translation and Gauss rotation.
use fossil,                     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, &
                                       swept_cells_3d, piece_map_3d_t, STATUS_OK
use fossil_boolean_bodies,      only : BODY_SPHERE, mc_body
use fossil_piece_map_3d_checks, only : check_map
use fossil_sign_bench_bodies,   only : body_t, bench_bodies, to_surface
use penf,                       only : I4P, I8P, R8P
use vecfor,                     only : vector_R8P

implicit none

integer(I4P), parameter   :: PATH_TRANSLATION = 1_I4P, PATH_ROTATION = 2_I4P, PATH_GAUSS = 3_I4P !< Motions.
character(20), parameter  :: PNAME(3) = [character(20) :: 'translation', 'rotation (exact)', 'rotation (Gauss)'] !< Names.
type(surface_stl_object)  :: sb(2)       !< Bodies (body frame).
type(body_t), allocatable :: bb(:)       !< Sign-bench bodies.
character(12)             :: names(2)    !< Body names.
real(R8P)                 :: half(2)     !< Half side of the box around each body.
character(32)             :: arg         !< Argument.
integer(I4P), allocatable :: sizes(:)    !< Timed sizes.
integer(I4P)              :: b, g, m     !< Work.
logical                   :: ok, full    !< Checks, the 64^3 timing.
type(vector_R8P)          :: b0, b1      !< Box.

full = .false.
if (command_argument_count() >= 1) then
   call get_command_argument(1, arg) ; full = trim(arg) == 'full'
endif
call mc_body(BODY_SPHERE, 24_I4P, sb(1))
call bench_bodies(bb) ; call to_surface(bb(7)%soup, sb(2))
b0 = sb(2)%get_bmin() ; b1 = sb(2)%get_bmax()
call sb(2)%translate(delta=-0.5_R8P * (b0 + b1)) ; call sb(2)%analyze
names = [character(12) :: 'mc sphere', 'cone']
half = [0.75_R8P, 1.1_R8P]
ok = .true.
print '(A)', 'checks (per-cell sums, numbering, residuals): body, cells, motion'
do b = 1, 2
   do m = 16, 24, 8
      ok = check_case(b, m, PATH_TRANSLATION) .and. ok
      if (m == 16) ok = check_case(b, m, PATH_ROTATION) .and. ok
      ok = check_case(b, m, PATH_GAUSS) .and. ok
   enddo
enddo
sizes = [32_I4P] ; if (full) sizes = [32_I4P, 64_I4P]
print '(A)', 'cost: time of swept_cells_3d without / with the map [s] (wall, -O2)'
do g = 1, size(sizes)
   do b = 1, 2
      call time_case(b, sizes(g), PATH_TRANSLATION)
      call time_case(b, sizes(g), PATH_GAUSS)
   enddo
enddo
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine motion(b, m, path, x, p0, p1)
   !< The grid of body b on m^3 cells and the end poses of the step (largest vertex displacement 2 cells).
   integer(I4P),           intent(in)  :: b, m, path !< Body, cells, motion.
   real(R8P), allocatable, intent(out) :: x(:)       !< Cell edges (the same on the three axes).
   type(rigid_pose_t),     intent(out) :: p0, p1     !< End poses.
   type(rigid_pose_t)                  :: st         !< Step rotation.
   real(R8P)                           :: h, d(3), r !< Work.
   integer(I4P)                        :: i          !< Work.

   h = 2._R8P * half(b) / m
   x = [(-half(b) + h * i, i = 0, m)]
   p0 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.1_R8P)
   if (path == PATH_TRANSLATION) then
      d = [0.7_R8P, 0.4_R8P, -0.3_R8P] ; d = 2._R8P * h * d / norm2(d)
      p1 = p0 ; p1%t = p0%t + d
   else
      r = merge(0.5_R8P, 0.95_R8P, b == 1) ! largest vertex radius
      d = [0.3_R8P, 0.2_R8P, -0.1_R8P] ; d = 0.5_R8P * h * d / norm2(d)
      st = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 1.5_R8P * h / r)
      p1%r = matmul(st%r, p0%r) ; p1%t = matmul(st%r, p0%t) + d
   endif
   endsubroutine motion

   logical function check_case(b, m, path) result(pass)
   !< The map checks of body b on m^3 cells for one motion.
   integer(I4P), intent(in) :: b, m, path !< Body, cells, motion.
   type(swept_cells_3d_t)   :: sw         !< Swept geometry.
   type(piece_map_3d_t)     :: map        !< Piece map.
   type(rigid_pose_t)       :: p0, p1     !< End poses.
   real(R8P), allocatable   :: x(:)       !< Cell edges.
   integer(I4P)             :: st         !< Status.

   call motion(b, m, path, x, p0, p1)
   call swept_cells_3d(sb(b), p0, p1, x, x, x, sw, exact=path /= PATH_GAUSS, status=st, map=map)
   print '(A,A,A,I0,A,A)', '  ', trim(names(b)), ' ', m, '^3 ', trim(PNAME(path))
   pass = st == STATUS_OK .and. check_map(sb(b), sw, map, p0, p1, x, path /= PATH_GAUSS)
   print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,8(I0,1X),A,L1,A,L1)', '    tracks ', map%nt, '  multi-track cells ', &
         map%multi_track_cells, '  evaluated cells ', map%evaluated_cells, '  evaluations ', map%evaluations, &
         '  untrusted cells ', map%untrusted_cells, '  statuses (-1, 1..7) ', map%n_status(-1), map%n_status(1:7), &
         '  exact ', map%exact, '  ok ', pass
   if (path == PATH_GAUSS) print '(A,ES9.2,A,ES9.2)', '    Gauss: max |per-track residual| ', maxval(abs(map%residual)), &
                                 '  sw error estimate ', sw%error_estimate
   endfunction check_case

   subroutine time_case(b, m, path)
   !< The time of `swept_cells_3d` of body b on m^3 cells without and with the map.
   integer(I4P), intent(in) :: b, m, path !< Body, cells, motion.
   type(swept_cells_3d_t)   :: sw         !< Swept geometry.
   type(piece_map_3d_t)     :: map        !< Piece map.
   type(rigid_pose_t)       :: p0, p1     !< End poses.
   real(R8P), allocatable   :: x(:)       !< Cell edges.
   integer(I8P)             :: c0, c1, c2, rate !< Clock.
   real(R8P)                :: t0, t1     !< Times.

   call motion(b, m, path, x, p0, p1)
   call system_clock(c0, rate)
   call swept_cells_3d(sb(b), p0, p1, x, x, x, sw, exact=path /= PATH_GAUSS)
   call system_clock(c1)
   call swept_cells_3d(sb(b), p0, p1, x, x, x, sw, exact=path /= PATH_GAUSS, map=map)
   call system_clock(c2)
   t0 = real(c1 - c0, R8P) / rate ; t1 = real(c2 - c1, R8P) / rate
   print '(A,A,A,I0,A,A,A,F9.3,A,F9.3,A,F6.2,A,I0,A,I0,A,I0,A,I0)', '  ', trim(names(b)), ' ', m, '^3 ', PNAME(path), &
         ': without ', t0, '  with ', t1, '  map / swept ', (t1 - t0) / t0, '  evaluated cells ', map%evaluated_cells, &
         '  evaluations ', map%evaluations, '  multi-track ', map%multi_track_cells, '  untrusted ', map%untrusted_cells
   endsubroutine time_case
endprogram fossil_val_piece_map_3d
