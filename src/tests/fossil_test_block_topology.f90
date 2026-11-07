!< FOSSIL test: curvilinear blocks and the detection of their topology records (issue #71, step 2).
program fossil_test_block_topology
!< FOSSIL test: curvilinear blocks and the detection of their topology records (issue #71, step 2).
!<
!< Analytic grids, where the records are known:
!< 1. a 2-D O-grid around a cylinder, seam copied exactly: one record on each of faces 1 and 2, onto each other;
!<    the same grid with the seam from cos/sin(2 pi) (round-off): nothing exact, both records with tolerance 1e-12;
!< 2. a 2-D C-grid around an ellipse with a wake cut: two records on face 3 with the index reversal, nothing else;
!<    the same C-grid extruded in 3-D: the same two records, spanning k;
!< 3. a 3-D polar grid with a collapsed axis: face 3 singular along i, faces 1 and 2 abut over the whole face;
!< 4. a sheared Cartesian 3-D block: no record;
!< 5. a periodic 2-D channel (dyadic coordinates: x + L - L = x), translation given: faces 1 and 2 periodic, +-L;
!< 6. two 2-D blocks sharing a face, the second with reversed j: one record each, onto each other;
!< 7. every record round-trips: mapped vertices coincide (after the translation) and the donor maps them back;
!< 8. bad input is rejected.

use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, detect_topology, &
                                     detect_topology_2d, record_map, TOPO_ABUT, TOPO_PERIODIC, TOPO_SINGULAR, &
                                     BLOCK_STATUS_OK, BLOCK_STATUS_BAD_INPUT
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: b2(:)
type(curvilinear_block_t),    allocatable :: b3(:)
real(R8P),                    allocatable :: u(:,:), v(:,:), x(:,:,:), y(:,:,:), z(:,:,:)
integer(I4P)                              :: i, j, k, ni, nj, nk, iw, st
logical                                   :: ok, all_ok

all_ok = .true.

! 1. O-grid, exact seam
ni = 33 ; nj = 9
call o_grid(ni, nj, .true., u, v)
allocate(b2(1)) ; call b2(1)%initialize(u, v, status=st) ; call check(st == BLOCK_STATUS_OK, 'O-grid initialize')
call detect_topology_2d(b2)
ok = b2(1)%nrec == 2 .and. b2(1)%rec(1)%face == 1 .and. b2(1)%rec(2)%face == 2
if (ok) ok = b2(1)%rec(1)%donor_face == 2 .and. b2(1)%rec(2)%donor_face == 1 .and. &
             all(b2(1)%rec(1)%transform == [1, 2, 3]) .and. b2(1)%rec(1)%lo(2) == 1 .and. b2(1)%rec(1)%hi(2) == nj
call check(ok, 'O-grid seam: faces 1 and 2 abut')
call round_trip_2d(b2, 0._R8P)
! the seam by round-off
call o_grid(ni, nj, .false., u, v)
call b2(1)%initialize(u, v) ; call detect_topology_2d(b2)
call check(b2(1)%nrec == 0, 'O-grid seam with round-off: nothing exact')
call b2(1)%initialize(u, v) ; call detect_topology_2d(b2, tolerance=1e-12_R8P)
call check(b2(1)%nrec == 2, 'O-grid seam with round-off: found with tolerance')
call round_trip_2d(b2, 1e-12_R8P)

! 2. C-grid
ni = 61 ; nj = 7 ; iw = 13
call c_grid(ni, nj, iw, u, v)
call b2(1)%initialize(u, v) ; call detect_topology_2d(b2)
ok = b2(1)%nrec == 2
if (ok) ok = all(b2(1)%rec(:)%face == 3) .and. all(b2(1)%rec(:)%donor_face == 3) .and. &
             b2(1)%rec(1)%lo(1) == 1 .and. b2(1)%rec(1)%hi(1) == iw .and. &
             b2(1)%rec(2)%lo(1) == ni - iw + 1 .and. b2(1)%rec(2)%hi(1) == ni .and. &
             b2(1)%rec(1)%transform(1) == -1 .and. b2(1)%rec(1)%transform(2) == -2
call check(ok, 'C-grid: wake cut, two records with the index reversal')
call round_trip_2d(b2, 0._R8P)
nk = 4
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   x(:,:,k) = u ; y(:,:,k) = v ; z(:,:,k) = 0.25_R8P * (k - 1)
enddo
allocate(b3(1)) ; call b3(1)%initialize(x, y, z) ; call detect_topology(b3)
ok = b3(1)%nrec == 2
if (ok) ok = all(b3(1)%rec(:)%face == 3) .and. all(b3(1)%rec(:)%lo(3) == 1) .and. all(b3(1)%rec(:)%hi(3) == nk) .and. &
             b3(1)%rec(1)%transform(3) == 3
call check(ok, 'C-grid extruded: the same two records, spanning k')
call round_trip_3d(b3)
deallocate(x, y, z)

! 3. polar grid with a collapsed axis
ni = 17 ; nj = 6 ; nk = 5
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         x(i, j, k) = 0.2_R8P * (j - 1) * cos(2 * PI * (i - 1) / (ni - 1))
         y(i, j, k) = 0.2_R8P * (j - 1) * sin(2 * PI * (i - 1) / (ni - 1))
         z(i, j, k) = 0.3_R8P * (k - 1)
      enddo
   enddo
enddo
x(ni,:,:) = x(1,:,:) ; y(ni,:,:) = y(1,:,:)
call b3(1)%initialize(x, y, z) ; call detect_topology(b3)
ok = b3(1)%nrec == 3
if (ok) ok = count(b3(1)%rec(:)%kind == TOPO_SINGULAR) == 1 .and. count(b3(1)%rec(:)%kind == TOPO_ABUT) == 2
if (ok) then
   do i = 1, 3
      if (b3(1)%rec(i)%kind == TOPO_SINGULAR) ok = ok .and. b3(1)%rec(i)%face == 3 .and. b3(1)%rec(i)%collapsed == 1
      if (b3(1)%rec(i)%kind == TOPO_ABUT) ok = ok .and. b3(1)%rec(i)%face <= 2 .and. b3(1)%rec(i)%lo(2) == 1 .and. &
                                               b3(1)%rec(i)%hi(2) == nj .and. b3(1)%rec(i)%hi(3) == nk
   enddo
endif
call check(ok, 'polar grid: singular axis, seam over the whole face')
call round_trip_3d(b3)
deallocate(x, y, z)

! 4. sheared Cartesian block
ni = 6 ; nj = 5 ; nk = 4
allocate(x(ni, nj, nk), y(ni, nj, nk), z(ni, nj, nk))
do k = 1, nk
   do j = 1, nj
      do i = 1, ni
         x(i, j, k) = 0.1_R8P * i + 0.03_R8P * j ; y(i, j, k) = 0.1_R8P * j + 0.02_R8P * k ; z(i, j, k) = 0.1_R8P * k
      enddo
   enddo
enddo
call b3(1)%initialize(x, y, z) ; call detect_topology(b3)
call check(b3(1)%nrec == 0, 'sheared Cartesian block: no record')
deallocate(x, y, z)

! 5. periodic channel
ni = 9 ; nj = 5
deallocate(u, v) ; allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   do i = 1, ni
      u(i, j) = 0.25_R8P * (i - 1) + 0.0625_R8P * mod(j, 3) ; v(i, j) = 0.2_R8P * (j - 1)
   enddo
enddo
call b2(1)%initialize(u, v) ; call detect_topology_2d(b2)
call check(b2(1)%nrec == 0, 'channel: nothing without the translation')
call b2(1)%initialize(u, v) ; call detect_topology_2d(b2, translations=reshape([2._R8P, 0._R8P], [2, 1]))
ok = b2(1)%nrec == 2
if (ok) ok = all(b2(1)%rec(:)%kind == TOPO_PERIODIC) .and. b2(1)%rec(1)%face == 1 .and. &
             b2(1)%rec(1)%translation(1) == 2._R8P .and. b2(1)%rec(2)%translation(1) == -2._R8P
call check(ok, 'channel: periodic faces 1 and 2, translations +-L')
call round_trip_2d(b2, 0._R8P)

! 6. two blocks sharing a face, the second with reversed j
deallocate(b2) ; allocate(b2(2))
ni = 5 ; nj = 4
deallocate(u, v) ; allocate(u(ni, nj), v(ni, nj))
do j = 1, nj
   do i = 1, ni
      u(i, j) = 0.1_R8P * (i - 1) ; v(i, j) = 0.1_R8P * (j - 1) + 0.16_R8P * ((i - 1) / 4._R8P) ** 2
   enddo
enddo
call b2(1)%initialize(u, v, name='left')
do j = 1, nj
   do i = 1, ni
      u(i, j) = 0.4_R8P + 0.1_R8P * (i - 1) ; v(i, j) = 0.1_R8P * (nj - j) + 0.16_R8P
   enddo
enddo
call b2(2)%initialize(u, v, name='right')
call detect_topology_2d(b2)
ok = b2(1)%nrec == 1 .and. b2(2)%nrec == 1
if (ok) ok = b2(1)%rec(1)%face == 2 .and. b2(1)%rec(1)%donor_block == 2 .and. b2(1)%rec(1)%donor_face == 1 .and. &
             b2(1)%rec(1)%transform(1) == 1 .and. b2(1)%rec(1)%transform(2) == -2 .and. b2(2)%rec(1)%donor_block == 1
call check(ok, 'two blocks: one interface record each, j reversed')
call round_trip_2d(b2, 0._R8P)

! 8. bad input
call b2(1)%initialize(u(1:1,:), v(1:1,:), status=st)
call check(st == BLOCK_STATUS_BAD_INPUT, 'bad input: one vertex line')
u(2, 2) = huge(1._R8P) ; u(2, 2) = u(2, 2) * 2
call b2(1)%initialize(u, v, status=st)
call check(st == BLOCK_STATUS_BAD_INPUT, 'bad input: infinite coordinate')

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   subroutine check(cond, label)
   !< Record a check.
   logical,      intent(in) :: cond  !< Condition.
   character(*), intent(in) :: label !< Label.

   print '(A,L1,2A)', '  ', cond, '  ', label
   all_ok = all_ok .and. cond
   endsubroutine check

   subroutine o_grid(ni, nj, exact_seam, u, v)
   !< O-grid around the unit cylinder, radius 1 to 2.
   integer(I4P),           intent(in)  :: ni, nj     !< Vertices.
   logical,                intent(in)  :: exact_seam !< Copy the seam.
   real(R8P), allocatable, intent(out) :: u(:,:), v(:,:) !< Vertices.
   integer(I4P)                        :: i, j       !< Counters.
   real(R8P)                           :: r, t       !< Polar coordinates.

   allocate(u(ni, nj), v(ni, nj))
   do j = 1, nj
      r = 1._R8P + (j - 1._R8P) / (nj - 1)
      do i = 1, ni
         t = 2 * PI * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   if (exact_seam) then
      u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   endif
   endsubroutine o_grid

   subroutine c_grid(ni, nj, iw, u, v)
   !< C-grid around an ellipse (chord 1, thickness 0.2) with a wake to x = 3: iw vertices per wake side.
   integer(I4P),           intent(in)  :: ni, nj, iw !< Vertices, wake vertices.
   real(R8P), allocatable, intent(out) :: u(:,:), v(:,:) !< Vertices.
   real(R8P)                           :: p(2), nrm(2), phi, r !< Work.
   integer(I4P)                        :: i, j       !< Counters.

   allocate(u(ni, nj), v(ni, nj))
   do i = 1, ni
      if (i <= iw) then
         p = [3._R8P - 2._R8P * (i - 1) / (iw - 1), 0._R8P] ; nrm = [0._R8P, -1._R8P]
      elseif (i >= ni - iw + 1) then
         p = [3._R8P - 2._R8P * (ni - i) / (iw - 1), 0._R8P] ; nrm = [0._R8P, 1._R8P]
      else
         phi = 2 * PI * (i - iw) / (ni - 2 * iw + 1)
         p = [0.5_R8P + 0.5_R8P * cos(phi), -0.1_R8P * sin(phi)]
         nrm = [0.1_R8P * cos(phi), -0.5_R8P * sin(phi)] ; nrm = nrm / norm2(nrm)
      endif
      do j = 1, nj
         r = 0.15_R8P * (j - 1)
         u(i, j) = p(1) + r * nrm(1) ; v(i, j) = p(2) + r * nrm(2)
      enddo
   enddo
   endsubroutine c_grid

   subroutine round_trip_2d(b, tol)
   !< Every record maps its vertices onto coincident donor vertices, and a donor record maps them back.
   type(curvilinear_block_2d_t), intent(in) :: b(:) !< Blocks.
   real(R8P),                    intent(in) :: tol  !< Tolerance.
   integer(I4P)                             :: bb, r, i, j, q(3), s, w(3) !< Work.
   logical                                  :: ok, back !< Flags.

   ok = .true.
   do bb = 1, size(b)
      do r = 1, b(bb)%nrec
         associate(rec => b(bb)%rec(r))
            if (rec%kind /= TOPO_ABUT .and. rec%kind /= TOPO_PERIODIC) cycle
            do j = rec%lo(2), rec%hi(2)
               do i = rec%lo(1), rec%hi(1)
                  q = record_map(rec, [i, j, 1])
                  associate(d => b(rec%donor_block))
                     ok = ok .and. abs(b(bb)%u(i, j) + rec%translation(1) - d%u(q(1), q(2))) <= tol .and. &
                                   abs(b(bb)%v(i, j) + rec%translation(2) - d%v(q(1), q(2))) <= tol
                     back = .false.
                     do s = 1, d%nrec
                        if (d%rec(s)%kind /= rec%kind .or. d%rec(s)%donor_block /= bb) cycle
                        if (any(q < d%rec(s)%lo) .or. any(q > d%rec(s)%hi)) cycle
                        w = record_map(d%rec(s), q) ; back = back .or. all(w == [i, j, 1])
                     enddo
                     ok = ok .and. back
                  endassociate
               enddo
            enddo
         endassociate
      enddo
   enddo
   call check(ok, '  records round-trip')
   endsubroutine round_trip_2d

   subroutine round_trip_3d(b)
   !< Every record maps its vertices onto coincident donor vertices, and a donor record maps them back.
   type(curvilinear_block_t), intent(in) :: b(:) !< Blocks.
   integer(I4P)                          :: bb, r, i, j, k, q(3), s, w(3) !< Work.
   logical                               :: ok, back !< Flags.

   ok = .true.
   do bb = 1, size(b)
      do r = 1, b(bb)%nrec
         associate(rec => b(bb)%rec(r))
            if (rec%kind /= TOPO_ABUT .and. rec%kind /= TOPO_PERIODIC) cycle
            do k = rec%lo(3), rec%hi(3)
               do j = rec%lo(2), rec%hi(2)
                  do i = rec%lo(1), rec%hi(1)
                     q = record_map(rec, [i, j, k])
                     associate(d => b(rec%donor_block))
                        ok = ok .and. b(bb)%x(i, j, k) == d%x(q(1), q(2), q(3)) .and. &
                                      b(bb)%y(i, j, k) == d%y(q(1), q(2), q(3)) .and. &
                                      b(bb)%z(i, j, k) == d%z(q(1), q(2), q(3))
                        back = .false.
                        do s = 1, d%nrec
                           if (d%rec(s)%kind /= rec%kind .or. d%rec(s)%donor_block /= bb) cycle
                           if (any(q < d%rec(s)%lo) .or. any(q > d%rec(s)%hi)) cycle
                           w = record_map(d%rec(s), q) ; back = back .or. all(w == [i, j, k])
                        enddo
                        ok = ok .and. back
                     endassociate
                  enddo
               enddo
            enddo
         endassociate
      enddo
   enddo
   call check(ok, '  records round-trip')
   endsubroutine round_trip_3d
endprogram fossil_test_block_topology
