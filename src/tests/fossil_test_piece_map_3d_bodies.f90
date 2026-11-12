!< FOSSIL test: the piece map of a time step in 3-D on a body (issue #102 step 2).
program fossil_test_piece_map_3d_bodies
!< FOSSIL test: the piece map of a time step in 3-D on a body (issue #102 step 2).
!<
!< The cone of the sign bench (32 faces) on 4^3 cells, translating (exact), then rotating (exact and Gauss): the map
!< checks of `fossil_piece_map_3d_checks` on every cell (numbering as an explicit `cut_cells_3d` of moved copies,
!< per-cell sums over the tracks equal to the swept geometry, the residual identity, per-track residuals at round-off
!< on the exact paths) and no untrusted cell. Then the refusal of `map` with explicit node poses. Last, the cone
!< centred and unrotated (its base parallel to the z planes) in the first step of `fossil_val_moving_conservation`'s
!< translation at 2 cells per step on 24^3, on the 3 x 3 x 3 cells around its cell (15, 12, 7): the base reaches the
!< cell's plane z = z_7 with corners and edges at one instant computed apart by round-off (intervals narrower than
!< EVENT_GAP), right after the evaluation of the covered cell (no piece). Before the follow-up of issue #102 the
!< quadrature nodes there were told from the others by the evaluations' first nodes, equal after an evaluation
!< without pieces: not skipped, they failed, and the cell (dead + fresh) was untrusted (2 of the 27 cells here; 28 cells
!< over the toy's 20 steps).
!< The hand-built cases are in `fossil_test_piece_map_3d`.
use fossil,                     only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, &
                                       swept_cells_3d, piece_map_3d_t, gauss_legendre_01, pose_interpolate, STATUS_OK, &
                                       STATUS_INVALID_INPUT, TRACK_DEAD, TRACK_FRESH
use fossil_piece_map_3d_checks, only : check_map
use fossil_sign_bench_bodies,   only : body_t, bench_bodies, to_surface
use penf,                       only : I4P, R8P
use vecfor,                     only : vector_R8P

implicit none

type(surface_stl_object)        :: s
type(swept_cells_3d_t)          :: sw
type(piece_map_3d_t)            :: map
type(body_t), allocatable       :: bb(:)
type(rigid_pose_t)              :: p0, p1
type(rigid_pose_t), allocatable :: pq(:)
real(R8P), allocatable          :: sq(:), wq(:), xg(:)
real(R8P)                       :: gx(4), gy(4), gz(4)
type(vector_R8P)                :: b0, b1
integer(I4P)                    :: i, st
logical                         :: ok(5)

! 1-3. the cone (32 faces) on 4^3 cells, translating, then rotating: the map checks on every cell
call bench_bodies(bb) ; call to_surface(bb(7)%soup, s)
xg = [(-0.3_R8P + 0.45_R8P * i, i = 0, 4)]
p0 = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.4_R8P, [0.31_R8P, 0.27_R8P, 0.12_R8P])
p1 = p0 ; p1%t = p0%t + [0.37_R8P, 0.21_R8P, 0.29_R8P]
ok(1) = case_grid('cone, translation', xg, p0, p1, .true.)
p1 = pose_from_axis_angle([0.2_R8P, -0.3_R8P, 1._R8P], 0.5_R8P) ; p1%r = matmul(p1%r, p0%r)
p1%t = p0%t + [0.2_R8P, 0.1_R8P, 0.15_R8P]
ok(2) = case_grid('cone, rotation', xg, p0, p1, .true.)
ok(3) = case_grid('cone, rotation, Gauss', xg, p0, p1, .false.)
! 4. map with explicit node poses: refused
call gauss_legendre_01(4_I4P, sq, wq)
allocate(pq(4))
do i = 1, 4
   pq(i) = pose_interpolate(p0, p1, sq(i))
enddo
call swept_cells_3d(s, p0, p1, xg, xg, xg, sw, nq=4_I4P, poses=pq, status=st, map=map)
ok(4) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'map with explicit node poses refused: ', ok(4)

! 5. the unrotated cone: events one instant apart right after a covered cell (`fossil_val_moving_conservation`)
call bench_bodies(bb) ; call to_surface(bb(7)%soup, s)
b0 = s%get_bmin() ; b1 = s%get_bmax()
call s%translate(delta=-0.5_R8P * (b0 + b1)) ; call s%analyze
p0 = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0._R8P)
p1 = p0 ; p1%t = [1.7304484372741907e-01_R8P, 1.6457542623823845e-01_R8P, 1.5101242748543242e-01_R8P]
gx = [1.4127485809335272e-01_R8P, 2.8254971618670544e-01_R8P, 4.2382457428005793e-01_R8P, 5.6509943237341087e-01_R8P]
gy = [-2.8254971618670544e-01_R8P, -1.4127485809335272e-01_R8P, 2.2204460492503131e-16_R8P, 1.4127485809335272e-01_R8P]
gz = [-9.8892400665346880e-01_R8P, -8.4764914856011597e-01_R8P, -7.0637429046676337e-01_R8P, -5.6509943237341065e-01_R8P]
call swept_cells_3d(s, p0, p1, gx, gy, gz, sw, status=st, map=map)
i = map%first_track(2, 2, 2)
ok(5) = st == STATUS_OK .and. map%untrusted_cells == 0_I4P .and. map%tracks(2, 2, 2) == 2_I4P
if (ok(5)) ok(5) = map%status(i) == TRACK_DEAD .and. map%status(i + 1) == TRACK_FRESH
print '(A,I0,A,I0,A,L1)', 'unrotated cone, base on a z plane after a covered cell: untrusted cells ', map%untrusted_cells, &
      '  skipped nodes ', map%skipped_nodes, '  middle cell dead + fresh: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   logical function case_grid(name, g, q0, q1, exact) result(pass)
   !< The map checks on every cell of the block g^3, and no untrusted cell.
   character(*),       intent(in) :: name   !< Case.
   real(R8P),          intent(in) :: g(:)   !< Cell edges.
   type(rigid_pose_t), intent(in) :: q0, q1 !< End poses.
   logical,            intent(in) :: exact  !< Exact paths.
   integer(I4P)                   :: st_    !< Status.

   call swept_cells_3d(s, q0, q1, g, g, g, sw, exact=exact, status=st_, map=map)
   pass = st_ == STATUS_OK .and. check_map(s, sw, map, q0, q1, g, exact) .and. map%untrusted_cells == 0_I4P
   print '(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A,L1)', '  ', name, ': tracks ', map%nt, '  multi-track cells ', &
         map%multi_track_cells, '  evaluated cells ', map%evaluated_cells, '  evaluations ', map%evaluations, &
         '  untrusted ', map%untrusted_cells, '  ok ', pass
   endfunction case_grid
endprogram fossil_test_piece_map_3d_bodies
