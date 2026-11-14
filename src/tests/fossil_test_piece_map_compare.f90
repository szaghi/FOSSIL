!< FOSSIL test: the piece map comparator (issue #105 step 1a).
program fossil_test_piece_map_compare
!< FOSSIL test: the piece map comparator (issue #105 step 1a).
!<
!< `piece_map_2d_mismatches` / `piece_map_3d_mismatches` (`fossil_piece_map_compare`) on small maps: two maps of the
!< same call give 0, in 2-D (an ellipse of 64 segments on 12^2 cells, exact rotation and Gauss with refine) and in 3-D
!< (the marching-cubes sphere of 8^3 samples on 6^3 cells, exact translation and Gauss rotation). The negative
!< controls, each on a copy of a map: one real element one ulp off (volume0, outflow, faces) gives exactly 1; an
!< integer element changed by one (status, old_piece) gives 1; a counter (evaluations, n_status) gives 1; the exact
!< flag gives 1; an array deallocated in one map gives its size; a -0 for a +0 gives 1 (bits, not values).
use fossil,                   only : curve_object, rigid_pose_2d_t, pose_2d, swept_cells_2d_t, swept_cells_2d, &
                                     piece_map_2d_t, surface_stl_object, rigid_pose_t, pose_from_axis_angle, &
                                     swept_cells_3d_t, swept_cells_3d, piece_map_3d_t
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_piece_map_compare, only : piece_map_2d_mismatches, piece_map_3d_mismatches
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter     :: PI = acos(-1._R8P) !< Pi.
type(curve_object)       :: c                  !< 2-D body.
type(surface_stl_object) :: s                  !< 3-D body.
type(rigid_pose_2d_t)    :: q0, q1             !< 2-D poses.
type(rigid_pose_t)       :: p0, p1             !< 3-D poses.
type(swept_cells_2d_t)   :: w2                 !< 2-D swept geometry.
type(swept_cells_3d_t)   :: w3                 !< 3-D swept geometry.
type(piece_map_2d_t)     :: a2, b2, e2         !< 2-D maps (e2: perturbed copy).
type(piece_map_3d_t)     :: a3, b3, e3         !< 3-D maps.
real(R8P), allocatable   :: u(:), v(:), x(:)   !< Ellipse, cell edges.
integer(I4P)             :: k, m               !< Work.
logical                  :: ok                 !< All checks.

ok = .true.
! 2-D: an ellipse rotating about its centre plus a small translation
allocate(u(64), v(64))
do k = 1_I4P, 64_I4P
   u(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 64) ; v(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 64)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
x = [(-0.02_R8P + 1.04_R8P * k / 12._R8P, k = 0, 12)]
q0 = pose_2d(0.1_R8P) ; q0%t = [0.5_R8P, 0.5_R8P] - matmul(q0%r, [0.5_R8P, 0.5_R8P])
q1 = pose_2d(0.25_R8P) ; q1%t = [0.5_R8P, 0.5_R8P] - matmul(q1%r, [0.5_R8P, 0.5_R8P]) + [0.02_R8P, 0.013_R8P]
call swept_cells_2d(c, q0, q1, x, x, w2, exact=.true., map=a2)
call swept_cells_2d(c, q0, q1, x, x, w2, exact=.true., map=b2)
call expect('2-D exact rotation, same call twice', piece_map_2d_mismatches(a2, b2, verbose=.true.), 0_I4P)
if (a2%nt < 1_I4P .or. a2%evaluated_cells < 1_I4P) call fail('2-D map has no evaluated cell')
call swept_cells_2d(c, q0, q1, x, x, w2, exact=.false., refine=.true., map=b2)
call swept_cells_2d(c, q0, q1, x, x, w2, exact=.false., refine=.true., map=e2)
call expect('2-D Gauss + refine, same call twice', piece_map_2d_mismatches(b2, e2, verbose=.true.), 0_I4P)
e2 = a2 ; e2%volume0(1) = nearest(e2%volume0(1), 1._R8P)
call expect('2-D volume0 one ulp off', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; m = size(e2%outflow, 2, kind=I4P) ; e2%outflow(4, m) = nearest(e2%outflow(4, m), -1._R8P)
call expect('2-D outflow one ulp off', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; e2%status(1) = e2%status(1) + 1_I4P
call expect('2-D status changed', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; e2%evaluations = e2%evaluations + 1_I4P
call expect('2-D evaluations changed', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; e2%n_status(-1) = e2%n_status(-1) + 1_I4P
call expect('2-D n_status(-1) changed', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; e2%exact = .not. e2%exact
call expect('2-D exact flag flipped', piece_map_2d_mismatches(a2, e2), 1_I4P)
e2 = a2 ; deallocate(e2%new_track)
call expect('2-D new_track deallocated', piece_map_2d_mismatches(a2, e2), max(1_I4P, size(a2%new_track, kind=I4P)))
e2 = a2 ; e2%residual(1) = 0._R8P ; b2 = e2 ; e2%residual(1) = -0._R8P
call expect('2-D -0 against +0', piece_map_2d_mismatches(b2, e2), 1_I4P)
! 3-D: the marching-cubes sphere, translating (exact) and rotating (Gauss)
call mc_body(BODY_SPHERE, 8_I4P, s)
x = [(-0.8_R8P + 1.6_R8P * k / 6._R8P, k = 0, 6)]
p0 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.1_R8P)
p1 = p0 ; p1%t = p0%t + [0.21_R8P, 0.12_R8P, -0.09_R8P]
call swept_cells_3d(s, p0, p1, x, x, x, w3, exact=.true., map=a3)
call swept_cells_3d(s, p0, p1, x, x, x, w3, exact=.true., map=b3)
call expect('3-D exact translation, same call twice', piece_map_3d_mismatches(a3, b3, verbose=.true.), 0_I4P)
if (a3%nt < 1_I4P .or. a3%evaluated_cells < 1_I4P) call fail('3-D map has no evaluated cell')
p1 = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0.3_R8P) ; p1%r = matmul(p1%r, p0%r)
p1%t = [0.05_R8P, 0.03_R8P, -0.02_R8P]
call swept_cells_3d(s, p0, p1, x, x, x, w3, exact=.false., map=b3)
call swept_cells_3d(s, p0, p1, x, x, x, w3, exact=.false., map=e3)
call expect('3-D Gauss rotation, same call twice', piece_map_3d_mismatches(b3, e3, verbose=.true.), 0_I4P)
e3 = a3 ; m = size(e3%faces, 2, kind=I4P) ; e3%faces(6, m) = nearest(e3%faces(6, m), 1._R8P)
call expect('3-D faces one ulp off', piece_map_3d_mismatches(a3, e3), 1_I4P)
e3 = a3 ; e3%old_piece(1) = e3%old_piece(1) + 1_I4P
call expect('3-D old_piece changed', piece_map_3d_mismatches(a3, e3), 1_I4P)
e3 = a3 ; e3%narrow_cells = e3%narrow_cells + 1_I4P
call expect('3-D narrow_cells changed', piece_map_3d_mismatches(a3, e3), 1_I4P)
e3 = a3 ; deallocate(e3%tracks)
call expect('3-D tracks deallocated', piece_map_3d_mismatches(a3, e3), max(1_I4P, size(a3%tracks, kind=I4P)))
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) stop 1
contains
   subroutine expect(label, got, want)
   !< Check a mismatch count.
   character(*), intent(in) :: label !< Case.
   integer(I4P), intent(in) :: got   !< Mismatches counted.
   integer(I4P), intent(in) :: want  !< Expected.

   print '(A,A,A,I0,A,I0)', '  ', label, ': mismatches ', got, ', expected ', want
   if (got /= want) call fail(label)
   endsubroutine expect

   subroutine fail(label)
   !< Record a failure.
   character(*), intent(in) :: label !< Case.

   print '(A,A)', '  FAILED: ', label
   ok = .false.
   endsubroutine fail
endprogram fossil_test_piece_map_compare
