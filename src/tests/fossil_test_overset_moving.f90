!< FOSSIL test: overset assembly of a moving block (issue #72, step 9).
program fossil_test_overset_moving
!< FOSSIL test: overset assembly of a moving block (issue #72, step 9).
!<
!< A cylinder O-grid kept in its body frame translates through a Cartesian background (h = 0.15) with a pose (the
!< cutter cuts the background only: posed wall vertices land off the wall by round-off):
!< 1. the posed assembly equals the one of the O-grid moved by hand (holes, candidate cells; xi to 1e-12: the two
!<    roundings differ);
!< 2. the incremental search (walks from the previous step's donors) equals the search from scratch, bit for bit;
!< 3. a step of 0.6 h: exposed points exist (old holes or receptors now field) and no hole becomes a field point
!<    directly (Meakin & Suhs); a step of 3 h: holes become field points, flagged.

use fossil,         only : curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, curve_object, wall_cutter_2d, &
                           rigid_pose_2d_t, pose_2d
use fossil_overset, only : donor_set_t, donor_candidates_2d, overset_mask_t, classify_2d, hole_flags_2d, &
                           block_logical_t, exposed_points, SELECT_EXPLICIT
use penf,           only : I4P, R8P

implicit none

real(R8P), parameter                      :: PI = acos(-1._R8P), H = 0.15_R8P
type(curvilinear_block_2d_t), allocatable :: q(:), qm(:)
type(curve_object)                        :: c, cm
type(donor_set_t),            allocatable :: d0(:), d1(:), di(:), dm(:)
type(overset_mask_t),         allocatable :: m0(:), m1(:), m2(:)
type(block_logical_t),        allocatable :: h0(:), h1(:), hm(:), ex(:)
type(rigid_pose_2d_t)                     :: p(2)
integer(I4P)                              :: nviol
logical                                   :: ok, all_ok

all_ok = .true.
allocate(q(2), qm(2))
call o_grid(q(1), 0._R8P) ; call background_2d(q(2))
call q(1)%build_search ; call q(2)%build_search
call wall_cutter_2d(q(1:1), c)

! 1. posed = moved by hand
p(1) = pose_2d(0._R8P, [0.31_R8P, -0.17_R8P]) ; p(2) = pose_2d(0._R8P)
call hole_flags_2d(q, c, h0, cuts=[.false., .true.], poses=p, cutter_pose=p(1))
call donor_candidates_2d(q, d0, poses=p)
call o_grid(qm(1), 0._R8P, [0.31_R8P, -0.17_R8P]) ; qm(2) = q(2)
call qm(1)%build_search ; call qm(2)%build_search
call wall_cutter_2d(qm(1:1), cm)
call hole_flags_2d(qm, cm, hm, cuts=[.false., .true.]) ; call donor_candidates_2d(qm, dm)
ok = all(h0(2)%v .eqv. hm(2)%v) .and. all(h0(1)%v .eqv. hm(1)%v) .and. all(d0(2)%start == dm(2)%start)
if (ok) ok = all(d0(2)%cell == dm(2)%cell) .and. all(d0(2)%status == dm(2)%status) .and. &
             maxval(abs(d0(2)%xi - dm(2)%xi)) <= 1e-12_R8P .and. all(d0(1)%start == dm(1)%start)
call check(ok, 'posed assembly = the O-grid moved by hand (same cells, xi to 1e-12)', real(count(h0(2)%v), R8P))

! 2-3. a step of 0.6 h
call classify_2d(q, d0, m0, SELECT_EXPLICIT, 2_I4P, holes=h0)
p(1) = pose_2d(0._R8P, [0.31_R8P + 0.6_R8P * H, -0.17_R8P])
call hole_flags_2d(q, c, h1, cuts=[.false., .true.], poses=p, cutter_pose=p(1))
call donor_candidates_2d(q, d1, poses=p)
call donor_candidates_2d(q, di, poses=p, previous=d0)
call check(same(d1(1), di(1)) .and. same(d1(2), di(2)), 'incremental search = search from scratch, bit for bit', 0._R8P)
call classify_2d(q, d1, m1, SELECT_EXPLICIT, 2_I4P, holes=h1)
call exposed_points(m0, m1, ex, nviol)
call check(count(ex(2)%v) > 0 .and. nviol == 0, 'step of 0.6 h: exposed points, no hole becomes a field point', &
           real(count(ex(2)%v), R8P))
p(1) = pose_2d(0._R8P, [0.31_R8P + 3.6_R8P * H, -0.17_R8P])
call hole_flags_2d(q, c, h1, cuts=[.false., .true.], poses=p, cutter_pose=p(1))
call donor_candidates_2d(q, d1, poses=p)
call classify_2d(q, d1, m2, SELECT_EXPLICIT, 2_I4P, holes=h1)
call exposed_points(m0, m2, ex, nviol)
call check(nviol > 0, 'step of 3 h: holes become field points, flagged', real(nviol, R8P))

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

   logical function same(a, b) result(yes)
   !< Two candidate sets equal bit for bit.
   type(donor_set_t), intent(in) :: a, b !< Sets.

   yes = all(a%start == b%start)
   if (yes) yes = all(a%block == b%block) .and. all(a%cell == b%cell) .and. all(a%xi == b%xi) .and. &
                  all(a%status == b%status)
   endfunction same

   subroutine o_grid(blk, t0, shift)
   !< O-grid 1 <= r <= 2 around the unit cylinder (33 x 13), optionally shifted; wall face 3, overlap face 4.
   type(curvilinear_block_2d_t), intent(inout)        :: blk      !< Block.
   real(R8P),                    intent(in)           :: t0       !< Start angle.
   real(R8P),                    intent(in), optional :: shift(2) !< Shift.
   real(R8P)                                          :: u(33, 13), v(33, 13), r !< Vertices.
   integer(I4P)                                       :: i, j     !< Counters.
   type(topo_record_t)                                :: rec      !< Records.

   do j = 1, 13
      r = 1._R8P + (j - 1) / 12._R8P
      do i = 1, 33
         u(i, j) = r * cos(t0 + 2 * PI * (i - 1) / 32) ; v(i, j) = r * sin(t0 + 2 * PI * (i - 1) / 32)
      enddo
   enddo
   u(33,:) = u(1,:) ; v(33,:) = v(1,:)
   if (present(shift)) then
      u = u + shift(1) ; v = v + shift(2)
   endif
   call blk%initialize(u, v)
   rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [33, 1, 1] ; rec%name = 'wall'
   call blk%add_record(rec)
   rec%face = 4 ; rec%lo = [1, 13, 1] ; rec%hi = [33, 13, 1] ; rec%name = 'overlap'
   call blk%add_record(rec)
   endsubroutine o_grid

   subroutine background_2d(blk)
   !< Cartesian block [-3, 3]^2, h = 0.15.
   type(curvilinear_block_2d_t), intent(inout) :: blk !< Block.
   real(R8P)                                   :: u(41, 41), v(41, 41) !< Vertices.
   integer(I4P)                                :: i, j !< Counters.

   do j = 1, 41
      do i = 1, 41
         u(i, j) = -3._R8P + H * (i - 1) ; v(i, j) = -3._R8P + H * (j - 1)
      enddo
   enddo
   call blk%initialize(u, v)
   endsubroutine background_2d
endprogram fossil_test_overset_moving
