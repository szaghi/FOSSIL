!< FOSSIL test: fringe levels, orphan causes and overlap width (issue #72, step 7).
program fossil_test_overset_fringes
!< FOSSIL test: fringe levels, orphan causes and overlap width (issue #72, step 7).
!<
!< 1. fringe levels (SUGGAR node rules) for a cylinder O-grid in a background cut by the cylinder: INNER 1 = the
!<    vertices sharing a cell with a hole, INNER 2 the next ring; OUTER 1 the O-grid's overlap boundary (33 vertices),
!<    OUTER 2 the ring inside it (33);
!< 2. each orphan cause produced on purpose: no candidate (an overlap boundary outside the background),
!<    a stencil hitting a hole (holes given next to the overlap boundary), receptors in the stencil (two blocks
!<    overlapping by 0.7 cell, explicit interpolation), a failed stencil (width above STENCIL_MAX), a wall vertex
!<    outside every donor (a concave wall discretized at two resolutions);
!< 3. overlap width of concentric annuli (1 <= r <= 2 and 1.5 <= r <= 3): d = 0.5.

use fossil,         only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d
use fossil_overset, only : donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, hole_flags_2d, &
                           block_logical_t, fringe_levels_2d, overlap_width_2d, SELECT_EXPLICIT, ORPHAN_NO_CANDIDATE, &
                           ORPHAN_STENCIL_HOLE, ORPHAN_RECEPTOR_DONOR, ORPHAN_BOUNDARY, ORPHAN_WALL, FRINGE_INNER_1, &
                           FRINGE_INNER_2, FRINGE_OUTER_1, FRINGE_OUTER_2
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curve_object)                        :: c
type(donor_set_t),            allocatable :: dn(:)
type(overset_mask_t),         allocatable :: m(:)
type(block_logical_t),        allocatable :: h(:)
real(R8P),                    allocatable :: d(:,:)
integer(I4P)                              :: i, j, n1, n2
logical                                   :: ok, all_ok

all_ok = .true.

! 1. fringe levels
allocate(q(2))
call annulus(q(1), 1._R8P, 2._R8P, 0._R8P, 2 * PI, 33_I4P, 13_I4P, 'wall', 'overlap')
call box_2d(q(2), -3._R8P, 3._R8P, -3._R8P, 3._R8P, 41_I4P, 41_I4P)
call q(1)%build_search ; call q(2)%build_search
call wall_cutter_2d(q(1:1), c) ; call hole_flags_2d(q, c, h)
call donor_candidates_2d(q, dn) ; call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P, holes=h)
call fringe_levels_2d(q, m)
n1 = 0 ; n2 = 0
do j = 1, 41
   do i = 1, 41
      if (h(2)%v(i + 41 * (j - 1))) cycle
      if (ring_of(h(2)%v, i, j)) then
         n1 = n1 + 1
      endif
   enddo
enddo
ok = count(m(2)%fringe == FRINGE_INNER_1) == n1 .and. count(m(2)%fringe == FRINGE_INNER_2) > 0 .and. &
     count(m(1)%fringe == FRINGE_OUTER_1) == 33 .and. count(m(1)%fringe == FRINGE_OUTER_2) == 33 .and. &
     all(m(2)%fringe(pack([(i, i = 1, 41 * 41)], m(2)%iota < 0)) == FRINGE_INNER_1)
call check(ok, 'fringe levels: INNER 1 = the hole ring, INNER 2 next; OUTER 1 and 2 on the O-grid', real(n1, R8P))

! 2a. no candidate: the background does not reach the O-grid's outer ring
call box_2d(q(2), -1.8_R8P, 1.8_R8P, -1.8_R8P, 1.8_R8P, 41_I4P, 41_I4P) ; call q(2)%build_search
call donor_candidates_2d(q, dn) ; call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P)
ok = count(m(1)%orphan == ORPHAN_NO_CANDIDATE) > 0
! 2b. holes next to the overlap boundary
call box_2d(q(2), -3._R8P, 3._R8P, -3._R8P, 3._R8P, 41_I4P, 41_I4P) ; call q(2)%build_search
call donor_candidates_2d(q, dn)
do j = 1, 41
   do i = 1, 41
      h(2)%v(i + 41 * (j - 1)) = norm2([q(2)%u(i, j) - 2._R8P, q(2)%v(i, j)]) < 0.3_R8P
   enddo
enddo
h(1)%v = .false.
call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P, holes=h)
ok = ok .and. count(m(1)%orphan == ORPHAN_STENCIL_HOLE) > 0
! 2c. failed stencils
call classify_2d(q, dn, m, SELECT_EXPLICIT, 9_I4P)
ok = ok .and. count(m(1)%orphan == ORPHAN_BOUNDARY) == 33
call check(ok, 'orphans: no candidate, stencil on a hole, failed stencil', 0._R8P)
! 2d. receptors in the stencil
deallocate(q) ; allocate(q(2))
call box_2d(q(1), 0._R8P, 1._R8P, 0._R8P, 1._R8P, 11_I4P, 11_I4P) ; call name_face(q(1), 2_I4P, 'overlap')
call box_2d(q(2), 0.93_R8P, 1.93_R8P, 0._R8P, 1._R8P, 11_I4P, 11_I4P) ; call name_face(q(2), 1_I4P, 'overlap')
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn) ; call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P)
call check(count(m(1)%orphan == ORPHAN_RECEPTOR_DONOR) == 11 .and. count(m(2)%orphan == ORPHAN_RECEPTOR_DONOR) == 11, &
           'orphans: receptors in the stencil (explicit interpolation)', 0._R8P)
! 2e. concave wall at two resolutions
call annulus(q(1), 1._R8P, 2._R8P, 0._R8P, 0.9_R8P * PI, 41_I4P, 9_I4P, '', 'wall')
call name_face(q(1), 2_I4P, 'overlap')
call annulus(q(2), 1._R8P, 2._R8P, 0.62_R8P * PI, 1.02_R8P * PI, 9_I4P, 9_I4P, '', 'wall')
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn) ; call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P)
call check(m(1)%orphan(41 * 9) == ORPHAN_WALL, 'orphans: a concave-wall vertex outside the coarser donor', 0._R8P)

! 3. overlap width
call annulus(q(1), 1._R8P, 2._R8P, 0._R8P, 2 * PI, 33_I4P, 9_I4P, '', 'overlap')
call annulus(q(2), 1.5_R8P, 3._R8P, 0._R8P, 2 * PI, 33_I4P, 9_I4P, 'overlap', '')
call q(1)%build_search ; call q(2)%build_search
call donor_candidates_2d(q, dn) ; call classify_2d(q, dn, m, SELECT_EXPLICIT, 2_I4P)
call overlap_width_2d(q, m, d)
call check(abs(d(1, 2) - 0.5_R8P) <= 1e-14_R8P .and. d(1, 2) == d(2, 1), 'overlap width of concentric annuli: 0.5', &
           d(1, 2))

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

   logical function ring_of(hole, i, j) result(yes)
   !< A vertex of the 41 x 41 background shares a cell with a hole.
   logical,      intent(in) :: hole(:) !< Holes.
   integer(I4P), intent(in) :: i, j    !< Vertex.
   integer(I4P)             :: di, dj  !< Offsets.

   yes = .false.
   do dj = -1, 1
      do di = -1, 1
         if (i + di < 1 .or. j + dj < 1 .or. i + di > 41 .or. j + dj > 41) cycle
         yes = yes .or. hole(i + di + 41 * (j + dj - 1))
      enddo
   enddo
   endfunction ring_of

   subroutine annulus(blk, r0, r1, t0, t1, ni, nj, inner, outer)
   !< Annulus r0 <= r <= r1 over [t0, t1] (seam copied on a full circle); faces 3 and 4 named inner and outer.
   type(curvilinear_block_2d_t), intent(inout) :: blk          !< Block.
   real(R8P),                    intent(in)    :: r0, r1, t0, t1 !< Extent.
   integer(I4P),                 intent(in)    :: ni, nj       !< Vertices.
   character(*),                 intent(in)    :: inner, outer !< Face names ('' none).
   real(R8P)                                   :: u(ni, nj), v(ni, nj), r, t !< Vertices.
   integer(I4P)                                :: i, j         !< Counters.
   type(topo_record_t)                         :: rec          !< Record.

   do j = 1, nj
      r = r0 + (r1 - r0) * (j - 1) / (nj - 1)
      do i = 1, ni
         t = t0 + (t1 - t0) * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   if (t1 - t0 >= 2 * PI) then
      u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   endif
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY
   if (len_trim(inner) > 0) then
      rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 1] ; rec%name = inner ; call blk%add_record(rec)
   endif
   if (len_trim(outer) > 0) then
      rec%face = 4 ; rec%lo = [1, nj, 1] ; rec%hi = [ni, nj, 1] ; rec%name = outer ; call blk%add_record(rec)
   endif
   endsubroutine annulus

   subroutine box_2d(blk, x0, x1, y0, y1, ni, nj)
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
   endsubroutine box_2d

   subroutine name_face(blk, face, name)
   !< Name a whole i face (1 or 2) of the block.
   type(curvilinear_block_2d_t), intent(inout) :: blk  !< Block.
   integer(I4P),                 intent(in)    :: face !< Face.
   character(*),                 intent(in)    :: name !< Name.
   type(topo_record_t)                         :: rec  !< Record.

   rec%kind = TOPO_BOUNDARY ; rec%face = face ; rec%name = name
   rec%lo = [merge(1_I4P, blk%ni, face == 1), 1_I4P, 1_I4P] ; rec%hi = [rec%lo(1), blk%nj, 1_I4P]
   call blk%add_record(rec)
   endsubroutine name_face
endprogram fossil_test_overset_fringes
