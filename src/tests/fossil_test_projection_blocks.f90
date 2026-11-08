!< FOSSIL test: projection of block faces onto bodies (issue #74, step 4).
!<
!< 1. A block over the cube face z = 0.5 (x, y in [-0.5, 0.5], 9 x 9) projected onto a UV sphere of radius 1 (16 x 32,
!<    no feature): every vertex on the sphere (projected again it moves by round-off only), radius within the chord
!<    sag of 1, the face's interior vertices moved outward by up to 1 - 0.5 sqrt(3)... (reported).
!< 2. A squircle face (9 x 9, the disc map of the square) of radius 1.005 at z = 1.99 under the top cap of a cylinder
!<    (radius 1, height 2, 32 sides, patches side / bottom / top), tolerance 0.02: its 32 boundary vertices go to the rim
!<    (z = 2, on a rim chord), its interior to the cap (z = 2, radius < 1); negative control: free projection puts the
!<    boundary vertices on the side, below the rim.
!< 3. 2-D: a ring of 32 vertices at radius 1.005 sqrt(2) around the square [-1, 1]^2 (the side j = 1 of an O-grid),
!<    tolerance 0.02: the 4 vertices at 45 degrees go to the corners bit for bit, the others onto the sides.
program fossil_test_projection_blocks

use fossil,                   only : surface_stl_object, curve_object, projection_target_t, curvilinear_block_t, &
                                     curvilinear_block_2d_t, project_block_face, project_block_side, project_points, &
                                     PROJ_FREE, PROJ_CORNER, PROJ_CURVE, PROJ_PATCH, PROJ_STATUS_OK
use fossil_projection_bodies, only : cylinder, uv_sphere
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter         :: PI = acos(-1._R8P)
type(surface_stl_object)     :: s
type(curve_object)           :: c
type(projection_target_t)    :: t
type(curvilinear_block_t)    :: b, bf
type(curvilinear_block_2d_t) :: b2
real(R8P)                    :: x(9, 9, 2), y(9, 9, 2), z(9, 9, 2), u(33, 2), v(33, 2), r, w, uu, vv, sag
real(R8P), allocatable       :: moved(:,:), m1(:), fx(:), fy(:), fz(:)
integer(I4P), allocatable    :: kinds(:,:), k1(:)
integer(I4P)                 :: i, j, st, nb
logical                      :: all_ok, ok

all_ok = .true.

! 1. a face onto a sphere
call uv_sphere(s, 16_I4P, 32_I4P, 1._R8P, .false.) ; call s%projection_target(t)
do j = 1, 9
   do i = 1, 9
      x(i, j, :) = -0.5_R8P + (i - 1) / 8._R8P ; y(i, j, :) = -0.5_R8P + (j - 1) / 8._R8P ; z(i, j, :) = [0.5_R8P, 0.25_R8P]
   enddo
enddo
call b%initialize(x, y, z)
call project_block_face(b, 5_I4P, t, 1e-6_R8P, moved=moved, status=st)
allocate(fx(81), fy(81), fz(81))
call project_points(t, pack(b%x(:, :, 1), .true.), pack(b%y(:, :, 1), .true.), pack(b%z(:, :, 1), .true.), &
                    [(PROJ_FREE, i = 1, 81)], [(0_I4P, i = 1, 81)], fx, fy, fz)
w = maxval(abs([fx - pack(b%x(:, :, 1), .true.), fy - pack(b%y(:, :, 1), .true.), fz - pack(b%z(:, :, 1), .true.)]))
sag = 1._R8P - cos(PI / 16._R8P)
r = maxval(abs(sqrt(b%x(:, :, 1)**2 + b%y(:, :, 1)**2 + b%z(:, :, 1)**2) - 1._R8P))
call check(st == PROJ_STATUS_OK .and. w <= 8 * epsilon(1._R8P) .and. r <= sag, &
           'face onto a sphere: on the surface (re-projection moves <= 8 ulps), radius within the chord sag', r)
print '(A,ES12.4)', '     largest move ', maxval(moved)

! 2. a squircle face under a cylinder cap
call cylinder(s, 32_I4P, 1._R8P, 2._R8P, .true.) ; call s%projection_target(t)
do j = 1, 9
   do i = 1, 9
      uu = -1._R8P + (i - 1) / 4._R8P ; vv = -1._R8P + (j - 1) / 4._R8P
      x(i, j, :) = 1.005_R8P * uu * sqrt(1._R8P - vv**2 / 2) ; y(i, j, :) = 1.005_R8P * vv * sqrt(1._R8P - uu**2 / 2)
      z(i, j, :) = [1.99_R8P, 1.5_R8P]
   enddo
enddo
call b%initialize(x, y, z) ; bf = b
call project_block_face(b, 5_I4P, t, 0.02_R8P, moved=moved, kinds=kinds, status=st)
ok = st == PROJ_STATUS_OK ; nb = 0
do j = 1, 9
   do i = 1, 9
      r = hypot(b%x(i, j, 1), b%y(i, j, 1))
      ok = ok .and. b%z(i, j, 1) == 2._R8P
      if (i == 1 .or. i == 9 .or. j == 1 .or. j == 9) then
         nb = nb + 1
         ok = ok .and. kinds(i, j) == PROJ_CURVE .and. r >= cos(PI / 32) - 1e-15_R8P .and. r <= 1._R8P + 1e-15_R8P
      else
         ok = ok .and. kinds(i, j) == PROJ_PATCH .and. r < 1._R8P
      endif
   enddo
enddo
call check(ok .and. nb == 32, 'squircle face: 32 boundary vertices on the rim, interior on the cap, all at z = 2', &
           real(nb, R8P))
call project_block_face(bf, 5_I4P, t, 0._R8P, kinds=kinds) ! tolerance 0: every vertex PROJ_PATCH of its nearest facet
call check(any(bf%z(1, :, 1) < 2._R8P) .and. any(abs(hypot(bf%x(1, :, 1), bf%y(1, :, 1)) - 1._R8P) < 0.05_R8P), &
           'negative control: unconstrained, boundary vertices land on the side below the rim (min z)', minval(bf%z(:, :, 1)))

! 3. 2-D ring around a square
call c%adopt_loops(u=[-1._R8P, 1._R8P, 1._R8P, -1._R8P], v=[-1._R8P, -1._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%projection_target(t)
do i = 1, 33
   u(i, 1) = 1.005_R8P * sqrt(2._R8P) * cos(PI / 4 + 2 * PI * (i - 1) / 32) ; u(i, 2) = 2 * u(i, 1)
   v(i, 1) = 1.005_R8P * sqrt(2._R8P) * sin(PI / 4 + 2 * PI * (i - 1) / 32) ; v(i, 2) = 2 * v(i, 1)
enddo
call b2%initialize(u, v)
call project_block_side(b2, 3_I4P, t, 0.02_R8P, moved=m1, kinds=k1, status=st)
ok = st == PROJ_STATUS_OK .and. count(k1 == PROJ_CORNER) == 5 ! the first and last vertex coincide
do i = 1, 33
   ok = ok .and. abs(max(abs(b2%u(i, 1)), abs(b2%v(i, 1))) - 1._R8P) <= 1e-15_R8P
   if (mod(i - 1, 8) == 0) ok = ok .and. abs(b2%u(i, 1)) == 1._R8P .and. abs(b2%v(i, 1)) == 1._R8P
enddo
call check(ok, '2-D ring: the 45-degree vertices on the corners bit for bit, the others on the sides', &
           real(count(k1 == PROJ_CORNER), R8P))

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
endprogram fossil_test_projection_blocks
