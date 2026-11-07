!< FOSSIL test: wall projection for the donor search of viscous walls (issue #72, step 8).
program fossil_test_overset_projection
!< FOSSIL test: wall projection for the donor search of viscous walls (issue #72, step 8).
!<
!< Two arcs of an annulus overlap on a curved wall, block 1 fine (41 vertices over 0.9 pi), block 2 coarse (9 over
!< 0.4 pi), each wall discretized by its own chords (Rogers et al. 2003):
!< 1. concave wall (the outer circle, fluid inside): without projection the fine wall vertices between the coarse
!<    vertices lie outside the coarse block (no candidate); with the projection every one finds a donor;
!< 2. convex wall (the inner circle, fluid outside): without projection the fine wall vertices find donors off the
!<    coarse wall (eta > 0, about the chord sag over the cell height); with the projection on it (eta below 1e-7: the
!<    1e-9 nudge into the fluid over the first cell height);
!< 3. the shift decays along the wall-normal line to zero at max_distance;
!< 4. 3-D: the convex case extruded in z gives eta below 1e-7 on the wall too.

use fossil,         only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY
use fossil_overset, only : donor_set_t, donor_candidates, donor_candidates_2d, wall_projection, wall_projection_2d, &
                           block_real_t
use fossil_block_locate, only : LOC_INSIDE
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(donor_set_t),            allocatable :: dn(:)
type(block_real_t)                        :: sh(2, 2)
type(block_real_t)                        :: sh3(2, 2)
integer(I4P)                              :: i, nmiss, nfound
real(R8P)                                 :: e0, e1
logical                                   :: all_ok

all_ok = .true.
allocate(q(2))

! 1. concave
call arc(q(1), 0._R8P, 0.9_R8P * PI, 41_I4P, 4_I4P) ; call arc(q(2), 0.62_R8P * PI, 1.02_R8P * PI, 9_I4P, 4_I4P)
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn)
nmiss = count_wall(q(1), dn(1), 4_I4P, .false.)
call wall_projection_2d(q, 1_I4P, 2_I4P, sh(1, 2)%v, max_distance=0.3_R8P)
call donor_candidates_2d(q, dn, shifts=sh)
nfound = count_wall(q(1), dn(1), 4_I4P, .true.)
call check(nmiss > 0 .and. nfound == overlap_wall(q(1)), &
           'concave wall: orphan wall vertices without projection, all found with it', real(nmiss, R8P))

! 2-3. convex
call arc(q(1), 0._R8P, 0.9_R8P * PI, 41_I4P, 3_I4P) ; call arc(q(2), 0.62_R8P * PI, 1.02_R8P * PI, 9_I4P, 3_I4P)
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn) ; e0 = wall_eta(q(1), dn(1))
deallocate(sh(1, 2)%v) ; call wall_projection_2d(q, 1_I4P, 2_I4P, sh(1, 2)%v, max_distance=0.3_R8P)
call donor_candidates_2d(q, dn, shifts=sh) ; e1 = wall_eta(q(1), dn(1))
call check(e0 > 1e-3_R8P .and. e1 <= 1e-7_R8P, 'convex wall: donors on the coarse wall with the projection (eta)', e1)
call check(all(sh(1, 2)%v(2 * (30 + 41 * 6) - 1:2 * (30 + 41 * 6)) == 0._R8P) .and. &
           norm2(sh(1, 2)%v(2 * (30 + 41) - 1:2 * (30 + 41))) < norm2(sh(1, 2)%v(2 * 30 - 1:2 * 30)), &
           'the shift decays along the wall-normal line', norm2(sh(1, 2)%v(2 * 30 - 1:2 * 30)))

! 4. 3-D
allocate(b(2))
call arc3(b(1), 0._R8P, 0.9_R8P * PI, 41_I4P) ; call arc3(b(2), 0.62_R8P * PI, 1.02_R8P * PI, 9_I4P)
call b(1)%build_search ; call b(2)%build_search
call wall_projection(b, 1_I4P, 2_I4P, sh3(1, 2)%v, max_distance=0.3_R8P)
call donor_candidates(b, dn, shifts=sh3)
e1 = 0._R8P
do i = 1, size(dn(1)%start) - 1
   if (dn(1)%start(i + 1) == dn(1)%start(i)) cycle
   if (mod((i - 1) / 41, 9) /= 0) cycle ! the wall layer j = 1
   e1 = max(e1, abs(dn(1)%xi(2, dn(1)%start(i))))
enddo
call check(e1 <= 1e-7_R8P, '3-D convex wall: eta = 0 with the projection', e1)

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

   subroutine arc(blk, t0, t1, ni, wall_face)
   !< Annulus 1 <= r <= 2 over [t0, t1], 9 radial vertices clustered at the wall; the wall on face 3 or 4.
   type(curvilinear_block_2d_t), intent(inout) :: blk       !< Block.
   real(R8P),                    intent(in)    :: t0, t1    !< Arc.
   integer(I4P),                 intent(in)    :: ni, wall_face !< Vertices, wall face.
   real(R8P)                                   :: u(ni, 9), v(ni, 9), r, t, s !< Vertices.
   integer(I4P)                                :: i, j      !< Counters.
   type(topo_record_t)                         :: rec       !< Wall.

   do j = 1, 9
      s = ((j - 1) / 8._R8P) ** 2 ; if (wall_face == 4) s = 1 - ((9 - j) / 8._R8P) ** 2
      r = 1._R8P + s
      do i = 1, ni
         t = t0 + (t1 - t0) * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = wall_face ; rec%name = 'wall'
   rec%lo = [1, merge(1, 9, wall_face == 3), 1] ; rec%hi = [ni, rec%lo(2), 1]
   call blk%add_record(rec)
   endsubroutine arc

   subroutine arc3(blk, t0, t1, ni)
   !< The convex 2-D arc extruded over 0 <= z <= 1 (3 layers), wall on face 3.
   type(curvilinear_block_t), intent(inout) :: blk    !< Block.
   real(R8P),                 intent(in)    :: t0, t1 !< Arc.
   integer(I4P),              intent(in)    :: ni     !< Vertices.
   real(R8P)                                :: x(ni, 9, 3), y(ni, 9, 3), z(ni, 9, 3), r, t !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.
   type(topo_record_t)                      :: rec    !< Wall.

   do k = 1, 3
      do j = 1, 9
         r = 1._R8P + ((j - 1) / 8._R8P) ** 2
         do i = 1, ni
            t = t0 + (t1 - t0) * (i - 1) / (ni - 1)
            x(i, j, k) = r * cos(t) ; y(i, j, k) = r * sin(t) ; z(i, j, k) = 0.5_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%name = 'wall' ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 3]
   call blk%add_record(rec)
   endsubroutine arc3

   integer(I4P) function overlap_wall(blk) result(n)
   !< Wall vertices of the fine arc strictly inside the coarse arc's angular range.
   type(curvilinear_block_2d_t), intent(in) :: blk !< Fine block.
   integer(I4P)                             :: i   !< Counter.

   n = 0
   do i = 1, 41
      if (0.9_R8P * PI * (i - 1) / 40 > 0.62_R8P * PI + 1e-9_R8P) n = n + 1
   enddo
   endfunction overlap_wall

   integer(I4P) function count_wall(blk, d, face, found) result(n)
   !< Wall vertices (face 3 or 4) of the fine block in the coarse block's range with (found) or without a candidate.
   type(curvilinear_block_2d_t), intent(in) :: blk   !< Fine block.
   type(donor_set_t),            intent(in) :: d     !< Its candidates.
   integer(I4P),                 intent(in) :: face  !< Wall face.
   logical,                      intent(in) :: found !< Count found (else missing).
   integer(I4P)                             :: i, v  !< Counters.

   n = 0
   do i = 1, 41
      if (0.9_R8P * PI * (i - 1) / 40 <= 0.62_R8P * PI + 1e-9_R8P) cycle
      v = i + 41 * merge(0, 8, face == 3)
      if (found .eqv. (d%start(v + 1) > d%start(v))) then
         if (found) then
            if (d%status(d%start(v)) == LOC_INSIDE) n = n + 1
         else
            n = n + 1
         endif
      endif
   enddo
   endfunction count_wall

   real(R8P) function wall_eta(blk, d) result(e)
   !< Largest local eta of the donors of the fine wall vertices in the coarse block's range.
   type(curvilinear_block_2d_t), intent(in) :: blk !< Fine block.
   type(donor_set_t),            intent(in) :: d   !< Candidates.
   integer(I4P)                             :: i   !< Counter.

   e = 0._R8P
   do i = 1, 41
      if (0.9_R8P * PI * (i - 1) / 40 <= 0.62_R8P * PI + 1e-9_R8P) cycle
      if (d%start(i + 1) == d%start(i)) cycle
      e = max(e, abs(d%xi(2, d%start(i))))
   enddo
   endfunction wall_eta
endprogram fossil_test_overset_projection
