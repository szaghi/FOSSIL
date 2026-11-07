!< FOSSIL test: overset classification, stencil validity and the fixed point (issue #72, step 6).
program fossil_test_overset_classify
!< FOSSIL test: overset classification, stencil validity and the fixed point (issue #72, step 6).
!<
!< 1. explicit: an O-grid around a cylinder (wall inside, overlap boundary outside) in a Cartesian background cut by
!<    the cylinder wall, q = 2: every vertex of the O-grid's outer ring receives from the background, the ring of
!<    background vertices next to the hole receives from the O-grid, nothing else does, no orphan;
!< 2. the classification does not depend on the order of the blocks (the blocks swapped give the same classes);
!< 3. two Cartesian blocks overlapping by 0.7 cell (vertex columns not coincident), both with an overlap boundary
!<    inside the other: with explicit interpolation their boundary receptors are orphans (each stencil holds the
!<    other's receptors); with implicit interpolation they are valid, and the circular pairs are counted (one per row);
!< 4. implicit hole cutting (no cutter): every receptor has a smaller donor cell or is mandatory, no hole;
!< 5. 3-D: a cylindrical O-grid in a Cartesian background: its outer ring receives, no orphan.

use fossil,         only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, &
                           wall_cutter_2d
use fossil_overset, only : donor_set_t, donor_candidates, donor_candidates_2d, overset_mask_t, classify, classify_2d, &
                           hole_flags_2d, block_logical_t, SELECT_EXPLICIT, SELECT_IMPLICIT
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:), qs(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(donor_set_t),            allocatable :: dn(:), ds(:)
type(overset_mask_t),         allocatable :: m(:), ms(:)
type(block_logical_t),        allocatable :: h(:), hs(:)
integer(I4P)                              :: it, circ, i, j, nring
logical                                   :: ok, all_ok

all_ok = .true.

! 1. explicit
allocate(q(2))
call o_grid(q(1)) ; call background_2d(q(2), -3._R8P, 3._R8P, -3._R8P, 3._R8P, 41_I4P, 41_I4P)
call q(1)%build_search ; call q(2)%build_search
call wall_cutter_2d(q(1:1), c)
call hole_flags_2d(q, c, h)
call donor_candidates_2d(q, dn)
call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P, holes=h, iterations=it)
nring = 0
do j = 1, 41
   do i = 1, 41
      if (h(2)%v(i + 41 * (j - 1))) cycle
      if (near_hole(i, j)) nring = nring + 1
   enddo
enddo
ok = all(m(1)%iota(33 * 12 + 1:33 * 13) == -2) .and. count(m(1)%iota == -2) == 33 .and. &
     count(m(2)%iota == -1) == nring .and. all(m(2)%iota(pack([(i, i = 1, 41 * 41)], h(2)%v)) == 0) .and. &
     all(m(1)%orphan == 0) .and. all(m(2)%orphan == 0)
call check(ok, 'explicit: the O-grid outer ring and the hole ring receive, nothing else', real(nring, R8P))

! 2. block order
allocate(qs(2)) ; qs(1) = q(2) ; qs(2) = q(1) ; allocate(hs(2)) ; hs(1) = h(2) ; hs(2) = h(1)
call donor_candidates_2d(qs, ds)
call classify_2d(qs, ds, ms, SELECT_EXPLICIT, 2_I4P, holes=hs)
ok = all(ms(1)%iota == remap(m(2)%iota)) .and. all(ms(2)%iota == remap(m(1)%iota))
call check(ok, 'the classes do not depend on the block order', real(it, R8P))

! 3. one-cell overlap
deallocate(q, dn, m) ; allocate(q(2))
call background_2d(q(1), 0._R8P, 1._R8P, 0._R8P, 1._R8P, 11_I4P, 11_I4P) ; call add_overlap(q(1), 2_I4P, 11_I4P)
call background_2d(q(2), 0.93_R8P, 1.93_R8P, 0._R8P, 1._R8P, 11_I4P, 11_I4P) ; call add_overlap(q(2), 1_I4P, 11_I4P)
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn)
call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P)
ok = count(m(1)%orphan > 0) > 0 .and. count(m(2)%orphan > 0) > 0
call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P, explicit_interpolation=.false., circular=circ)
ok = ok .and. all(m(1)%orphan == 0) .and. all(m(2)%orphan == 0) .and. circ >= 11
call check(ok, 'one-cell overlap: orphans with explicit interpolation, circular pairs with implicit', real(circ, R8P))

! 4. implicit hole cutting
deallocate(q, dn, m) ; allocate(q(2))
call o_grid(q(1)) ; call background_2d(q(2), -3._R8P, 3._R8P, -3._R8P, 3._R8P, 41_I4P, 41_I4P)
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn)
call classify_2d(q, dn, m, SELECT_IMPLICIT, 2_I4P)
ok = .true.
do i = 1, 2
   do j = 1, size(m(i)%iota)
      if (m(i)%iota(j) >= 0) cycle
      ok = ok .and. (dn(i)%ratio(m(i)%donor(j)) < 1._R8P .or. (i == 1 .and. j > 33 * 12))
   enddo
enddo
ok = ok .and. count(m(1)%iota < 0) > 0 .and. count(m(2)%iota < 0) > 0 .and. all(m(1)%iota /= 0) .and. &
     all(m(2)%iota /= 0)
call check(ok, 'implicit: receptors have smaller donors or are mandatory; no hole', &
           real(count(m(1)%iota < 0) + count(m(2)%iota < 0), R8P))

! 5. 3-D
allocate(b(2))
call cylinder_3d(b(1)) ; call background_3d(b(2))
call b(1)%build_search ; call b(2)%build_search
call donor_candidates(b, dn)
call classify(b, dn, m, SELECT_EXPLICIT, 2_I4P)
call check(count(m(1)%iota == -2) == 25 * 4 .and. all(m(1)%orphan == 0), '3-D: the O-grid outer ring receives', &
           real(count(m(1)%iota == -2), R8P))

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

   logical function near_hole(i, j) result(yes)
   !< A background vertex shares a cell with a hole (Chebyshev ring 1).
   integer(I4P), intent(in) :: i, j   !< Vertex.
   integer(I4P)             :: di, dj !< Offsets.

   yes = .false.
   do dj = -1, 1
      do di = -1, 1
         if (i + di < 1 .or. j + dj < 1 .or. i + di > 41 .or. j + dj > 41) cycle
         yes = yes .or. h(2)%v(i + di + 41 * (j + dj - 1))
      enddo
   enddo
   endfunction near_hole

   elemental integer(I4P) function remap(k) result(r)
   !< The class with blocks 1 and 2 swapped.
   integer(I4P), intent(in) :: k !< Class.

   r = k
   if (k == -1) r = -2
   if (k == -2) r = -1
   endfunction remap

   subroutine o_grid(blk)
   !< O-grid 1 <= r <= 2 around the unit cylinder, 33 x 13; wall on face 3, overlap boundary on face 4.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: u(33, 13), v(33, 13), r !< Vertices.
   integer(I4P)                                :: i, j !< Counters.
   type(topo_record_t)                         :: rec  !< Records.

   do j = 1, 13
      r = 1._R8P + (j - 1) / 12._R8P
      do i = 1, 33
         u(i, j) = r * cos(2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(2 * PI * (i - 1) / 32)
      enddo
   enddo
   u(33,:) = u(1,:) ; v(33,:) = v(1,:)
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine o_grid

   subroutine background_2d(blk, x0, x1, y0, y1, ni, nj)
   !< Cartesian block.
   type(curvilinear_block_2d_t), intent(inout) :: blk            !< Block.
   real(R8P),                    intent(in)    :: x0, x1, y0, y1 !< Box.
   integer(I4P),                 intent(in)    :: ni, nj         !< Vertices.
   real(R8P)                                   :: u(ni, nj), v(ni, nj) !< Vertices.
   integer(I4P)                                :: i, j           !< Counters.

   do j = 1, nj
      do i = 1, ni
         u(i, j) = x0 + (x1 - x0) * (i - 1) / (ni - 1) ; v(i, j) = y0 + (y1 - y0) * (j - 1) / (nj - 1)
      enddo
   enddo
   call blk%initialize(u, v)
   endsubroutine background_2d

   subroutine add_overlap(blk, face, n)
   !< Name an i face of the block an overlap boundary.
   type(curvilinear_block_2d_t), intent(inout) :: blk  !< Block.
   integer(I4P),                 intent(in)    :: face, n !< Face (1 or 2), vertices along it.
   type(topo_record_t)                         :: rec  !< Record.

   rec%kind = TOPO_BOUNDARY ; rec%face = face ; rec%name = 'overlap'
   rec%lo = [merge(1_I4P, blk%ni, face == 1), 1_I4P, 1_I4P] ; rec%hi = [rec%lo(1), n, 1_I4P]
   call blk%add_record(rec)
   endsubroutine add_overlap

   subroutine cylinder_3d(blk)
   !< Cylindrical O-grid 1 <= r <= 1.5, 0 <= z <= 1 (seam copied), overlap boundary on face 4.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(25, 5, 4), y(25, 5, 4), z(25, 5, 4), r !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.
   type(topo_record_t)                      :: rec     !< Record.

   do k = 1, 4
      do j = 1, 5
         r = 1._R8P + 0.125_R8P * (j - 1)
         do i = 1, 25
            x(i, j, k) = r * cos(2 * PI * (i - 1) / 24) ; y(i, j, k) = r * sin(2 * PI * (i - 1) / 24)
            z(i, j, k) = (k - 1) / 3._R8P
         enddo
      enddo
   enddo
   x(25,:,:) = x(1,:,:) ; y(25,:,:) = y(1,:,:)
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 4 ; rec%lo = [1, 5, 1] ; rec%hi = [25, 5, 4] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine cylinder_3d

   subroutine background_3d(blk)
   !< Cartesian block around the cylinder.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(25, 25, 8), y(25, 25, 8), z(25, 25, 8) !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 8
      do j = 1, 25
         do i = 1, 25
            x(i, j, k) = -1.83_R8P + 0.1525_R8P * (i - 1) ; y(i, j, k) = -1.83_R8P + 0.1525_R8P * (j - 1)
            z(i, j, k) = -0.33_R8P + 0.25_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   endsubroutine background_3d
endprogram fossil_test_overset_classify
