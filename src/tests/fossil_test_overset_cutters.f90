!< FOSSIL test: wall cutters of curvilinear blocks (issue #72, step 2).
program fossil_test_overset_cutters
!< FOSSIL test: wall cutters of curvilinear blocks (issue #72, step 2).
!<
!< 2-D winding numbers by the exact sum.
!< 1. a 2-D O-grid around the unit cylinder, wall on face 3: one closed loop, solid on the left (winding 1 at the
!<    centre, 0 outside; signed distance negative inside);
!< 2. two 2-D blocks whose wall arcs overlap by 0.4 pi (fine and coarse): shared(1, 2); trimmed, the winding number at
!<    the centre is 1 within 0.02 (a residual strip one coarse element wide); untrimmed (tolerance 0) it counts the
!<    overlap twice (1.4);
!< 3. a Yin-Yang pair of spherical shells (Yang coarser): shared; trimmed, winding 1 within 0.02 at the centre and 0
!<    within 0.02 outside; untrimmed above 1.05; the gap between bilinear faces and triangles is reported (> 0, small);
!< 4. a background block without walls shares nothing; a block set without walls gives status 1.

use fossil_overset,           only : wall_cutter, wall_cutter_2d
use fossil_curvilinear_block, only : curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY
use fossil_surface_stl_object, only : surface_stl_object
use fossil_curve_object,      only : curve_object
use fossil,                   only : SIGN_FAST_WINDING_NUMBER
use vecfor,                   only : vector_R8P
use penf,                     only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P)
type(curvilinear_block_2d_t), allocatable :: q(:)
type(curvilinear_block_t),    allocatable :: b(:)
type(curve_object)                        :: c
type(surface_stl_object)                  :: s
logical,      allocatable                 :: sh(:,:)
integer(I4P), allocatable                 :: own(:)
real(R8P)                                 :: w, gap, d(1)
integer(I4P)                              :: st
logical                                   :: all_ok

all_ok = .true.

! 1. O-grid
allocate(q(1))
call arc_block(q(1), 0._R8P, 2 * PI, 65_I4P)
call wall_cutter_2d(q, c, owner=own)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call c%distance_many(reshape([0.5_R8P, 0._R8P], [2, 1]), d, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
call check(abs(w - 1) < 1e-12_R8P .and. abs(c%winding_number([2._R8P, 0._R8P], use_brute_force=.true.)) < 1e-12_R8P &
           .and. d(1) < 0, &
           'O-grid wall: closed loop, solid on the left', w)

! 2. overlapping arcs
deallocate(q) ; allocate(q(2))
call arc_block(q(1), -0.2_R8P * PI, 1.4_R8P * PI, 49_I4P)
call arc_block(q(2), 0.8_R8P * PI, 1.4_R8P * PI, 23_I4P)
call wall_cutter_2d(q, c, shared=sh)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call check(sh(1, 2) .and. abs(w - 1) < 0.02_R8P, 'overlapping arcs: shared wall, trimmed winding 1', w)
call wall_cutter_2d(q, c, tolerance=0._R8P)
w = c%winding_number([0._R8P, 0._R8P], use_brute_force=.true.)
call check(w > 1.15_R8P, 'overlapping arcs untrimmed: overlap counted twice', w)

! 3. Yin-Yang
allocate(b(3))
call yin_block(b(1), 25_I4P, 9_I4P, .false.)
call yin_block(b(2), 17_I4P, 7_I4P, .true.)
call box_block(b(3))
call wall_cutter(b, s, gap=gap, shared=sh, status=st)
w = s%winding_number(vector_R8P(0._R8P, 0._R8P, 0._R8P))
call check(st == 0 .and. sh(1, 2) .and. .not. sh(1, 3) .and. .not. sh(2, 3) .and. abs(w - 1) < 0.02_R8P .and. &
           abs(s%winding_number(vector_R8P(3._R8P, 0.1_R8P, 0.2_R8P))) < 0.02_R8P, &
           'Yin-Yang: shared wall, trimmed winding 1 inside, 0 outside; box shares nothing', w)
call check(gap > 0._R8P .and. gap < 0.05_R8P, 'Yin-Yang: bilinear-to-triangle gap reported', gap)
call wall_cutter(b, s, tolerance=0._R8P)
w = s%winding_number(vector_R8P(0._R8P, 0._R8P, 0._R8P))
call check(w > 1.05_R8P, 'Yin-Yang untrimmed: overlap counted twice', w)

! 4. no walls
call wall_cutter(b(3:3), s, status=st)
call check(st == 1, 'no wall faces: status 1', 0._R8P)

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

   subroutine arc_block(blk, t0, span, ni)
   !< 2-D block over the arc [t0, t0 + span] of the annulus 1 <= r <= 1.5, wall on face 3; a full circle copies the seam.
   type(curvilinear_block_2d_t), intent(inout) :: blk    !< Block.
   real(R8P),                    intent(in)    :: t0, span !< Arc.
   integer(I4P),                 intent(in)    :: ni     !< Vertices along the arc.
   real(R8P)                                   :: u(ni, 5), v(ni, 5), t, r !< Vertices.
   integer(I4P)                                :: i, j   !< Counters.
   type(topo_record_t)                         :: rec    !< Wall record.

   do j = 1, 5
      r = 1._R8P + 0.125_R8P * (j - 1)
      do i = 1, ni
         t = t0 + span * (i - 1) / (ni - 1) ; u(i, j) = r * cos(t) ; v(i, j) = r * sin(t)
      enddo
   enddo
   if (span >= 2 * PI) then
      u(ni,:) = u(1,:) ; v(ni,:) = v(1,:)
   endif
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine arc_block

   subroutine yin_block(blk, ni, nj, yang)
   !< Spherical shell 1 <= r <= 2 over |lon| <= 3 pi / 4 + 0.15, |lat| <= pi / 4 + 0.15; Yang = (-x, z, y) of Yin.
   type(curvilinear_block_t), intent(inout) :: blk    !< Block.
   integer(I4P),              intent(in)    :: ni, nj !< Vertices in longitude, latitude.
   logical,                   intent(in)    :: yang   !< Yang copy.
   real(R8P)                                :: x(ni, nj, 3), y(ni, nj, 3), z(ni, nj, 3), lo, la, r, p(3) !< Data.
   integer(I4P)                             :: i, j, k !< Counters.
   type(topo_record_t)                      :: rec    !< Wall record.

   do k = 1, 3
      r = 1._R8P + 0.5_R8P * (k - 1)
      do j = 1, nj
         la = (PI / 4 + 0.15_R8P) * (2._R8P * (j - 1) / (nj - 1) - 1)
         do i = 1, ni
            lo = (3 * PI / 4 + 0.15_R8P) * (2._R8P * (i - 1) / (ni - 1) - 1)
            p = r * [cos(la) * cos(lo), cos(la) * sin(lo), sin(la)]
            if (yang) p = [-p(1), p(3), p(2)]
            x(i, j, k) = p(1) ; y(i, j, k) = p(2) ; z(i, j, k) = p(3)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   rec%kind = TOPO_BOUNDARY ; rec%face = 5 ; rec%lo = [1, 1, 1] ; rec%hi = [ni, nj, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   endsubroutine yin_block

   subroutine box_block(blk)
   !< A Cartesian background block around the sphere, no wall.
   type(curvilinear_block_t), intent(inout) :: blk !< Block.
   real(R8P)                                :: x(5, 5, 5), y(5, 5, 5), z(5, 5, 5) !< Vertices.
   integer(I4P)                             :: i, j, k !< Counters.

   do k = 1, 5
      do j = 1, 5
         do i = 1, 5
            x(i, j, k) = -3 + 1.5_R8P * (i - 1) ; y(i, j, k) = -3 + 1.5_R8P * (j - 1) ; z(i, j, k) = -3 + 1.5_R8P * (k - 1)
         enddo
      enddo
   enddo
   call blk%initialize(x, y, z)
   endsubroutine box_block
endprogram fossil_test_overset_cutters
