!< FOSSIL test: exact hole flags of block vertices against wall cutters (issue #72, step 3).
program fossil_test_overset_holes
!< FOSSIL test: exact hole flags of block vertices against wall cutters (issue #72, step 3).
!<
!< 1. 2-D: the wall of an O-grid around the unit cylinder cuts a Cartesian background block: the holes are exactly the
!<    vertices inside the cutter polygon (exact winding sum); the O-grid itself, cut too, keeps every vertex (its wall
!<    vertices are on the cutter: field points); the device path gives the same flags;
!< 2. offset 0.2: the holes grow to the vertices within 0.2 of the wall, and only those;
!< 3. 3-D: the closed wall of a spherical shell (seam, collapsed poles) cuts a background block: the holes equal the
!<    exact inside test (signed by the solid angle); the shell keeps every vertex; device = CPU.

use fossil,                   only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, &
                                     surface_stl_object, curve_object, SIGN_SOLID_ANGLE, wall_cutter, wall_cutter_2d
use fossil_overset,           only : hole_flags, hole_flags_2d, block_logical_t, block_real_t
use vecfor,                   only : vector_R8P
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(surface_stl_object)                  :: s
type(block_logical_t),        allocatable :: h(:), hd(:), ho(:)
type(block_real_t),           allocatable :: dist(:)
type(vector_R8P),             allocatable :: pts(:)
real(R8P),                    allocatable :: dref(:)
integer(I4P)                              :: i, j, k, n, nbad
logical                                   :: all_ok

all_ok = .true.

! 1-2. 2-D
allocate(q(2))
call o_grid(q(1)) ; call background_2d(q(2))
call wall_cutter_2d(q(1:1), c)
call hole_flags_2d(q, c, h)
nbad = 0 ; n = 0
do j = 1, q(2)%nj
   do i = 1, q(2)%ni
      n = n + 1
      if (h(2)%v(n) .neqv. (c%winding_number([q(2)%u(i, j), q(2)%v(i, j)], use_brute_force=.true.) > 0.5_R8P)) &
         nbad = nbad + 1
   enddo
enddo
call check(nbad == 0 .and. count(h(2)%v) > 0 .and. .not. any(h(1)%v), &
           '2-D: background holes = inside the cutter; the O-grid keeps its vertices', real(count(h(2)%v), R8P))
call hole_flags_2d(q, c, hd, device=.true.)
call check(all(hd(2)%v .eqv. h(2)%v) .and. all(hd(1)%v .eqv. h(1)%v), '2-D device flags = CPU', 0._R8P)
call hole_flags_2d(q, c, ho, offset=0.2_R8P, distances=dist)
nbad = count(ho(2)%v .neqv. (dist(2)%v < 0.2_R8P)) + count(h(2)%v .and. .not. ho(2)%v)
call check(nbad == 0 .and. count(ho(2)%v) > count(h(2)%v) .and. all(dist(2)%v(pack([(n, n = 1, size(ho(2)%v))], &
           ho(2)%v .and. .not. h(2)%v)) >= 0._R8P), '2-D offset 0.2: holes grow to the near vertices only', &
           real(count(ho(2)%v), R8P))

! 3. 3-D
allocate(b(2))
call shell(b(1)) ; call background_3d(b(2))
call wall_cutter(b(1:1), s)
call hole_flags(b, s, h)
allocate(pts(size(h(2)%v)), dref(size(h(2)%v))) ; n = 0
do k = 1, b(2)%nk
   do j = 1, b(2)%nj
      do i = 1, b(2)%ni
         n = n + 1 ; pts(n) = vector_R8P(b(2)%x(i, j, k), b(2)%y(i, j, k), b(2)%z(i, j, k))
      enddo
   enddo
enddo
call s%distance_many(points=pts, distances=dref, is_signed=.true., sign_algorithm=SIGN_SOLID_ANGLE, is_square_root=.true.)
call check(all(h(2)%v .eqv. (dref < 0._R8P)) .and. count(h(2)%v) > 0 .and. .not. any(h(1)%v), &
           '3-D: background holes = exact inside test; the shell keeps its vertices', real(count(h(2)%v), R8P))
call hole_flags(b, s, hd, device=.true.)
call check(all(hd(2)%v .eqv. h(2)%v), '3-D device flags = CPU', 0._R8P)

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
   real(R8P)                                   :: u(49, 7), v(49, 7), r !< Vertices.
   integer(I4P)                                :: i, j !< Counters.
   type(topo_record_t)                         :: rec  !< Wall.

   do j = 1, 7
      r = 1._R8P + 0.5_R8P * (j - 1) / 6
      do i = 1, 49
         u(i, j) = r * cos(2 * PI * (i - 1) / 48) ; v(i, j) = r * sin(2 * PI * (i - 1) / 48)
      enddo
   enddo
   u(49,:) = u(1,:) ; v(49,:) = v(1,:)
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [49, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine o_grid

   subroutine background_2d(blk)
   !< Cartesian block [-2.03, 2.03]^2, 41 x 41 vertices.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: u(41, 41), v(41, 41) !< Vertices.
   integer(I4P)                                :: i, j !< Counters.

   do j = 1, 41
      do i = 1, 41
         u(i, j) = -2.03_R8P + 0.1015_R8P * (i - 1) ; v(i, j) = -2.03_R8P + 0.1015_R8P * (j - 1)
      enddo
   enddo
   call blk%initialize(u, v)
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
endprogram fossil_test_overset_holes
