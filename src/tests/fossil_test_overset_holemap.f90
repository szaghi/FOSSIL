!< FOSSIL test: Cartesian hole map and object X-rays against the exact hole flags (issue #72, step 4).
program fossil_test_overset_holemap
!< FOSSIL test: Cartesian hole map and object X-rays against the exact hole flags (issue #72, step 4).
!<
!< 1. 2-D cylinder in a background block: the hole map (crossing parity) equals the exact flags (hole_flags_2d);
!< 2. 2-D thin trailing edge (a lens 0.008 thick at most, the map elements 0.06 wide): the hole map with the parity
!<    still equals the exact flags; PEGASUS 5's line of sight leaves vertices unresolved and misses holes;
!< 3. 2-D X-rays: a vertex differs from the exact flags only within half the ray spacing of the wall;
!< 4. 3-D sphere (closed shell wall) in a background block: hole map = exact flags; X-rays within half the raster
!<    diagonal of the wall.

use fossil,         only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, &
                           surface_stl_object, curve_object, wall_cutter, wall_cutter_2d
use fossil_overset, only : hole_flags, hole_flags_2d, hole_map, hole_map_2d, xray_flags, xray_flags_2d, &
                           block_logical_t, block_real_t
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(surface_stl_object)                  :: s
type(block_logical_t),        allocatable :: hx(:), hm(:), hp(:), hr(:)
type(block_real_t),           allocatable :: dist(:)
real(R8P)                                 :: u(41), v(41), h, bb(2)
integer(I4P)                              :: i, nu, nbad
logical                                   :: all_ok

all_ok = .true.

! 1. cylinder
allocate(q(2))
call o_grid(q(1)) ; call background_2d(q(2), -2.03_R8P, 2.03_R8P, -2.03_R8P, 2.03_R8P, 41_I4P, 41_I4P)
call wall_cutter_2d(q(1:1), c)
call hole_flags_2d(q, c, hx)
call hole_map_2d(q, c, hm, [24_I4P, 24_I4P])
call check(all(hm(2)%v .eqv. hx(2)%v) .and. all(hm(1)%v .eqv. hx(1)%v), '2-D cylinder: hole map = exact flags', &
           real(count(hm(2)%v), R8P))

! 2. thin trailing edge
do i = 1, 20 ! upper side from the trailing edge to the leading edge, then the lower side back
   u(i) = 1._R8P - (i - 1) / 19._R8P ; v(i) = 0.004_R8P * sin(PI * u(i)) + 0.0001_R8P
enddo
do i = 21, 40
   u(i) = (i - 21) / 19._R8P ; v(i) = -0.004_R8P * sin(PI * u(i)) - 0.0001_R8P
enddo
call c%adopt_loops(u=u(1:40), v=v(1:40), loop_sizes=[40_I4P])
deallocate(q) ; allocate(q(1))
call background_2d(q(1), -0.1_R8P, 1.1_R8P, -0.0503_R8P, 0.0503_R8P, 241_I4P, 41_I4P)
call hole_flags_2d(q, c, hx)
call hole_map_2d(q, c, hm, [20_I4P, 20_I4P])
call hole_map_2d(q, c, hp, [20_I4P, 20_I4P], pegasus=.true., unresolved=nu)
call check(all(hm(1)%v .eqv. hx(1)%v) .and. count(hx(1)%v) > 0, 'thin trailing edge: hole map with parity = exact', &
           real(count(hx(1)%v), R8P))
call check(nu > 0 .and. count(hp(1)%v) < count(hx(1)%v), 'thin trailing edge: PEGASUS 5 line of sight misses holes', &
           real(count(hx(1)%v) - count(hp(1)%v), R8P))

! 3. X-rays, 2-D
deallocate(q) ; allocate(q(2))
call o_grid(q(1)) ; call background_2d(q(2), -2.03_R8P, 2.03_R8P, -2.03_R8P, 2.03_R8P, 41_I4P, 41_I4P)
call wall_cutter_2d(q(1:1), c)
call hole_flags_2d(q, c, hx, distances=dist)
call xray_flags_2d(q, c, hr, 13_I4P)
bb = c%get_bmax() - c%get_bmin() ; h = bb(1) / 12
nbad = count((hr(2)%v .neqv. hx(2)%v) .and. abs(dist(2)%v) > 0.5_R8P * h)
call check(nbad == 0, '2-D X-rays: differences only within half the ray spacing of the wall', &
           real(count(hr(2)%v .neqv. hx(2)%v), R8P))

! 4. 3-D sphere
allocate(b(2))
call shell(b(1)) ; call background_3d(b(2))
call wall_cutter(b(1:1), s)
call hole_flags(b, s, hx, distances=dist)
call hole_map(b, s, hm, [12_I4P, 12_I4P, 12_I4P])
call check(all(hm(2)%v .eqv. hx(2)%v) .and. all(hm(1)%v .eqv. hx(1)%v), '3-D sphere: hole map = exact flags', &
           real(count(hm(2)%v), R8P))
call xray_flags(b, s, hr, [13_I4P, 13_I4P])
h = sqrt(2._R8P) * 2 / 12
nbad = count((hr(2)%v .neqv. hx(2)%v) .and. abs(dist(2)%v) > 0.5_R8P * h)
call check(nbad == 0, '3-D X-rays: differences only within half the raster diagonal of the wall', &
           real(count(hr(2)%v .neqv. hx(2)%v), R8P))

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label, val)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.
   real(R8P),    intent(in) :: val   !< Measured value.

   print '(A,L1,2A,ES11.3)', '  ', cond, '  ', label, val
   all_ok = all_ok .and. cond
   endsubroutine check

   subroutine o_grid(blk)
   !< O-grid around the unit cylinder, radius 1 to 1.5, wall on face 3.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: uu(49, 7), vv(49, 7), r !< Vertices.
   integer(I4P)                                :: i, j !< Counters.
   type(topo_record_t)                         :: rec  !< Wall.

   do j = 1, 7
      r = 1._R8P + 0.5_R8P * (j - 1) / 6
      do i = 1, 49
         uu(i, j) = r * cos(2 * PI * (i - 1) / 48) ; vv(i, j) = r * sin(2 * PI * (i - 1) / 48)
      enddo
   enddo
   uu(49,:) = uu(1,:) ; vv(49,:) = vv(1,:)
   call blk%initialize(uu, vv)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [49, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine o_grid

   subroutine background_2d(blk, x0, x1, y0, y1, ni, nj)
   !< Cartesian block.
   type(curvilinear_block_2d_t), intent(inout) :: blk              !< Block.
   real(R8P),                    intent(in)    :: x0, x1, y0, y1   !< Box.
   integer(I4P),                 intent(in)    :: ni, nj           !< Vertices.
   real(R8P)                                   :: uu(ni, nj), vv(ni, nj) !< Vertices.
   integer(I4P)                                :: i, j             !< Counters.

   do j = 1, nj
      do i = 1, ni
         uu(i, j) = x0 + (x1 - x0) * (i - 1) / (ni - 1) ; vv(i, j) = y0 + (y1 - y0) * (j - 1) / (nj - 1)
      enddo
   enddo
   call blk%initialize(uu, vv)
   endsubroutine background_2d

   subroutine shell(blk)
   !< Spherical shell 1 <= r <= 1.5 (azimuth seam copied, poles collapsed), wall on face 5.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(25, 13, 3), y(25, 13, 3), z(25, 13, 3), r, th, ph !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.
   type(topo_record_t)                      :: rec     !< Wall.

   do k = 1, 3
      r = 1._R8P + 0.25_R8P * (k - 1)
      do j = 1, 13
         ph = PI * (j - 1) / 12
         do i = 1, 25
            th = 2 * PI * (i - 1) / 24
            x(i, j, k) = r * sin(ph) * cos(th) ; y(i, j, k) = r * sin(ph) * sin(th) ; z(i, j, k) = r * cos(ph)
         enddo
         if (j == 1 .or. j == 13) then
            x(:, j, k) = 0._R8P ; y(:, j, k) = 0._R8P ; z(:, j, k) = merge(r, -r, j == 1)
         endif
      enddo
   enddo
   x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:) ; z(25,:,:) = z(1,:,:)
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [25, 13, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine shell

   subroutine background_3d(blk)
   !< Cartesian block [-1.53, 1.53]^3, 17^3 vertices.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(17, 17, 17), y(17, 17, 17), z(17, 17, 17) !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 17
      do j = 1, 17
         do i = 1, 17
            x(i, j, k) = -1.53_R8P + 0.19125_R8P * (i - 1) ; y(i, j, k) = -1.53_R8P + 0.19125_R8P * (j - 1)
            z(i, j, k) = -1.53_R8P + 0.19125_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   endsubroutine background_3d
endprogram fossil_test_overset_holemap
