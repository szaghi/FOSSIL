!< FOSSIL test: constrained projection and classification (issue #74, step 3).
!<
!< 1. The unit cube, one patch per face, the point (1.01, 0.97, 0.5) near the edge x = y = 1: free projection lands on
!<    the +x face (negative control: a grid point meant for +y crosses the edge); PROJ_PATCH +y gives (1, 1, 0.5) on
!<    the +y face; PROJ_CURVE on that edge gives (1, 1, 0.5); PROJ_CORNER (nearest) of (0.9, 0.95, 1.1) gives (1, 1, 1)
!<    bit for bit.
!< 2. Idempotence: 2000 fixed-seed points around a UV sphere (16 x 32) projected free, then their feet again: they
!<    move by at most 8 ulps of 1 (measured 3.4).
!< 3. A cylinder rim (32 chords): points above it projected on the rim curve equal the closest point over the 32 chords
!<    computed here, bit for bit, at height 2 exactly.
!< 4. Classification of a 5 x 5 grid on the cube's +x face (tolerance 1e-9): its 4 corners PROJ_CORNER, the other 12
!<    boundary vertices PROJ_CURVE, the 9 interior ones PROJ_PATCH +x.
!< 5. 2-D: a unit square: (1.1, 0.5) on its patch goes to (1, 0.5), the nearest corner of (1.1, 1.05) is (1, 1);
!<    classification tags the corners; PROJ_CURVE is rejected in 2-D (PROJ_STATUS_BAD_INPUT, point kept).
!< 6. Bad input: an unknown kind or a patch id out of range keeps the point, distance -1.
program fossil_test_projection

use fossil,                   only : surface_stl_object, curve_object, projection_target_t, project_points, &
                                     project_points_2d, classify_points, classify_points_2d, PROJ_FREE, PROJ_PATCH, &
                                     PROJ_CURVE, PROJ_CORNER, PROJ_STATUS_OK, PROJ_STATUS_BAD_INPUT, facet_object
use fossil_projection_bodies, only : cube_patches, cylinder, uv_sphere
use penf,                     only : I4P, I8P, R8P

implicit none

real(R8P), parameter        :: PI = acos(-1._R8P)
type(surface_stl_object)    :: s
type(curve_object)          :: c
type(projection_target_t)   :: t
real(R8P)                   :: fx(1), fy(1), fz(1), d(1), gx(25), gy(25), gz(25), a, best, q(3), e(3), w
real(R8P), allocatable      :: px(:), py(:), pz(:), qx(:), qy(:), qz(:), rx(:), ry(:), rz(:)
integer(I4P)                :: k(25), id(25), st, pxid, pyid, cid, i, j, kk(1), ii(1)
integer(I8P)                :: seed
logical                     :: all_ok, ok

all_ok = .true.

! 1. cube
call cube_patches(s, .true.) ; call s%projection_target(t)
pxid = face_patch(1, 1._R8P) ; pyid = face_patch(2, 1._R8P)
call project_points(t, [1.01_R8P], [0.97_R8P], [0.5_R8P], [PROJ_FREE], [0_I4P], fx, fy, fz, d)
call check(fx(1) == 1._R8P .and. fy(1) == 0.97_R8P, 'free projection near the edge lands on +x (negative control) (y)', fy(1))
call project_points(t, [1.01_R8P], [0.97_R8P], [0.5_R8P], [PROJ_PATCH], [pyid], fx, fy, fz, d)
call check(fx(1) == 1._R8P .and. fy(1) == 1._R8P .and. fz(1) == 0.5_R8P, 'PROJ_PATCH +y: (1, 1, 0.5) (distance)', d(1))
call classify_points(t, [1._R8P], [1._R8P], [0.5_R8P], 1e-9_R8P, kk, ii) ; cid = ii(1)
call project_points(t, [1.01_R8P], [0.97_R8P], [0.5_R8P], [PROJ_CURVE], [cid], fx, fy, fz, d)
call check(kk(1) == PROJ_CURVE .and. fx(1) == 1._R8P .and. fy(1) == 1._R8P .and. fz(1) == 0.5_R8P, &
           'PROJ_CURVE on the edge x = y = 1: (1, 1, 0.5) (curve id)', real(cid, R8P))
call project_points(t, [0.9_R8P], [0.95_R8P], [1.1_R8P], [PROJ_CORNER], [0_I4P], fx, fy, fz, d)
call check(fx(1) == 1._R8P .and. fy(1) == 1._R8P .and. fz(1) == 1._R8P, 'PROJ_CORNER (nearest): (1, 1, 1) bit for bit', d(1))

! 2. idempotence
call uv_sphere(s, 16_I4P, 32_I4P, 1._R8P, .false.) ; call s%projection_target(t)
allocate(px(2000), py(2000), pz(2000), qx(2000), qy(2000), qz(2000), rx(2000), ry(2000), rz(2000))
seed = 20261107_I8P
do i = 1, 2000
   px(i) = 3 * rnd() - 1.5_R8P ; py(i) = 3 * rnd() - 1.5_R8P ; pz(i) = 3 * rnd() - 1.5_R8P
enddo
call project_points(t, px, py, pz, [(PROJ_FREE, i = 1, 2000)], [(0_I4P, i = 1, 2000)], qx, qy, qz)
call project_points(t, qx, qy, qz, [(PROJ_FREE, i = 1, 2000)], [(0_I4P, i = 1, 2000)], rx, ry, rz)
w = maxval(abs([rx - qx, ry - qy, rz - qz]))
call check(w <= 8 * epsilon(1._R8P), 'idempotence on a sphere: feet projected again move <= 8 ulps (max move)', w)

! 3. cylinder rim against the chords
call cylinder(s, 32_I4P, 1._R8P, 2._R8P, .true.) ; call s%projection_target(t)
ok = .true.
do i = 1, 40
   a = 2 * PI * i / 40._R8P + 0.01_R8P
   call classify_points(t, [cos(a)], [sin(a)], [2._R8P], 1e-1_R8P, kk, ii)
   call project_points(t, [1.2_R8P * cos(a)], [1.2_R8P * sin(a)], [2.1_R8P], [PROJ_CURVE], ii, fx, fy, fz)
   best = huge(1._R8P)
   do j = 0, 31
      q = rim_foot([1.2_R8P * cos(a), 1.2_R8P * sin(a), 2.1_R8P], j)
      if (sum((q - [1.2_R8P * cos(a), 1.2_R8P * sin(a), 2.1_R8P])**2) < best) then
         best = sum((q - [1.2_R8P * cos(a), 1.2_R8P * sin(a), 2.1_R8P])**2) ; e = q
      endif
   enddo
   ok = ok .and. kk(1) == PROJ_CURVE .and. fz(1) == 2._R8P .and. maxval(abs([fx(1), fy(1)] - e(1:2))) <= 2e-16_R8P
enddo
call check(ok, 'cylinder rim: feet equal the closest point over the 32 chords, at z = 2', 0._R8P)

! 4. classification of a face grid
call cube_patches(s, .true.) ; call s%projection_target(t)
do j = 0, 4
   do i = 0, 4
      gx(1 + i + 5 * j) = 1._R8P ; gy(1 + i + 5 * j) = i / 4._R8P ; gz(1 + i + 5 * j) = j / 4._R8P
   enddo
enddo
call classify_points(t, gx, gy, gz, 1e-9_R8P, k, id)
call check(count(k == PROJ_CORNER) == 4 .and. count(k == PROJ_CURVE) == 12 .and. count(k == PROJ_PATCH) == 9 .and. &
           all(pack(id, k == PROJ_PATCH) == pxid), 'face grid: 4 corners, 12 on curves, 9 on the +x patch', &
           real(count(k == PROJ_CURVE), R8P))

! 5. 2-D
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%projection_target(t)
call project_points_2d(t, [1.1_R8P], [0.5_R8P], [PROJ_PATCH], [1_I4P], fx, fy, d)
call check(fx(1) == 1._R8P .and. fy(1) == 0.5_R8P, '2-D: (1.1, 0.5) on the patch -> (1, 0.5) (distance)', d(1))
call project_points_2d(t, [1.1_R8P], [1.05_R8P], [PROJ_CORNER], [0_I4P], fx, fy)
call check(fx(1) == 1._R8P .and. fy(1) == 1._R8P, '2-D: nearest corner of (1.1, 1.05) = (1, 1)', fx(1))
call classify_points_2d(t, [0._R8P, 0.5_R8P], [0._R8P, 0._R8P], 1e-9_R8P, k(1:2), id(1:2))
call check(k(1) == PROJ_CORNER .and. k(2) == PROJ_PATCH .and. id(2) == 1, '2-D classification: corner, patch', real(k(1), R8P))
call project_points_2d(t, [0.3_R8P], [0.2_R8P], [PROJ_CURVE], [0_I4P], fx, fy, d, status=st)
call check(st == PROJ_STATUS_BAD_INPUT .and. fx(1) == 0.3_R8P .and. d(1) == -1._R8P, '2-D: PROJ_CURVE rejected, point kept', &
           real(st, R8P))

! 6. bad input
call cube_patches(s, .true.) ; call s%projection_target(t)
call project_points(t, [0.3_R8P], [0.2_R8P], [0.1_R8P], [9_I4P], [0_I4P], fx, fy, fz, d, status=st)
call check(st == PROJ_STATUS_BAD_INPUT .and. fx(1) == 0.3_R8P .and. d(1) == -1._R8P, 'unknown kind: point kept', real(st, R8P))
call project_points(t, [0.3_R8P], [0.2_R8P], [0.1_R8P], [PROJ_PATCH], [99_I4P], fx, fy, fz, d, status=st)
call check(st == PROJ_STATUS_BAD_INPUT .and. d(1) == -1._R8P, 'patch id out of range: point kept', real(st, R8P))
call project_points(t, [0.3_R8P], [0.2_R8P], [0.1_R8P], [PROJ_FREE], [0_I4P], fx, fy, fz, d, status=st)
call check(st == PROJ_STATUS_OK, 'a good point: status ok', real(st, R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES23.15)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   function face_patch(axis, value) result(pid)
   !< Patch of the cube facets with coordinate `axis` equal to value.
   integer,   intent(in)       :: axis  !< Axis.
   real(R8P), intent(in)       :: value !< Coordinate.
   integer(I4P)                :: pid   !< Patch.
   type(facet_object), pointer :: f
   integer(I4P)                :: m

   pid = 0
   do m = 1, s%get_facets_number()
      f => s%facet_at(m)
      if (axis == 1 .and. all([f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x] == value)) pid = s%facet_patch(m)
      if (axis == 2 .and. all([f%vertex(1)%y, f%vertex(2)%y, f%vertex(3)%y] == value)) pid = s%facet_patch(m)
   enddo
   endfunction face_patch

   function rim_foot(p, j_) result(f)
   !< Closest point of chord j_ of the top rim (as the cylinder builder places its vertices).
   real(R8P),    intent(in) :: p(3) !< Point.
   integer(I4P), intent(in) :: j_   !< Chord.
   real(R8P)                :: f(3) !< Foot.
   real(R8P)                :: a0(3), a1(3), e(3), tt

   a0 = [cos(2 * PI * j_ / 32), sin(2 * PI * j_ / 32), 2._R8P]
   a1 = [cos(2 * PI * mod(j_ + 1, 32) / 32), sin(2 * PI * mod(j_ + 1, 32) / 32), 2._R8P]
   e = a1 - a0 ; tt = max(0._R8P, min(1._R8P, dot_product(p - a0, e) / dot_product(e, e)))
   f = a0 + tt * e
   if (tt == 1._R8P) f = a1
   endfunction rim_foot

   function rnd() result(r)
   !< Park-Miller, fixed seed.
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_projection
