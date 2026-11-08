!< FOSSIL test: wall-spacing and quality report of curvilinear blocks (issue #74, step 5).
!<
!< 1. 2-D O-grid around the unit circle (65 x 21), radii 1 + h1 (q^k - 1) / (q - 1) with h1 = 1e-3, q = 1.2, wall side
!<    j = 1: y1 along the grid line and to the wall cutter (wall_cutter_2d) = h1, growth ratio q (1e-12), with
!<    u_tau = 0.05 and nu = 1e-5: y1+ = 5 everywhere, all above the target 1; the report file is written.
!< 2. A channel block (wall y = 0, 9 x 9, the same spacing): cell-centred heights are h1 / 2 (the vertex ones h1).
!< 3. The flat-plate estimate (u_inf = 10, nu = 1.5e-5, x from the leading edge): u_tau = u_inf sqrt(0.0288 Re_x^-0.2)
!<    and the report flags it as an estimate.
!< 4. A 2-D block with a folded cell: one cell not certified valid in the report.
!< 5. 3-D: the O-grid extruded along z (3 layers) with its wall face j = 1: y1 = h1 both ways, growth q.
program fossil_test_block_report

use fossil, only : surface_stl_object, curve_object, curvilinear_block_t, curvilinear_block_2d_t, topo_record_t, &
                   TOPO_BOUNDARY, wall_cutter, wall_cutter_2d, block_report_t, block_report, block_report_2d, &
                   save_block_report, utau_flat_plate
use penf,   only : I4P, R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P), H1 = 1e-3_R8P, Q = 1.2_R8P
type(curvilinear_block_2d_t)    :: b2(1)
type(curvilinear_block_t)       :: b3(1)
type(curve_object)              :: c
type(surface_stl_object)        :: s
type(block_report_t)            :: rep
type(topo_record_t)             :: rec
real(R8P)                       :: u(65, 21), v(65, 21), x(65, 21, 3), y(65, 21, 3), z(65, 21, 3), r, ut
real(R8P)                       :: cu(9, 9), cv(9, 9)
integer(I4P)                    :: i, j, st
integer                         :: uu
character(80)                   :: line
logical                         :: all_ok

all_ok = .true.

! 1. 2-D O-grid
do j = 1, 21
   r = 1._R8P + H1 * (Q**(j - 1) - 1._R8P) / (Q - 1._R8P)
   do i = 1, 65
      u(i, j) = r * cos(2 * PI * mod(i - 1, 64) / 64) ; v(i, j) = r * sin(2 * PI * mod(i - 1, 64) / 64)
   enddo
enddo
call b2(1)%initialize(u, v)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [65, 1, 1] ; rec%name = 'wall'
call b2(1)%add_record(rec)
call wall_cutter_2d(b2, c)
call block_report_2d(b2(1), c, 1e-5_R8P, rep, u_tau=0.05_R8P, status=st)
call check(st == 0 .and. rep%nw == 65 .and. maxval(abs(rep%y1_line - H1)) <= 1e-15_R8P .and. &
           maxval(abs(rep%y1_wall - H1)) <= 1e-15_R8P, '2-D O-grid: y1 = h1 along the line and to the wall (max error)', &
           maxval(abs(rep%y1_wall - H1)))
call check(abs(rep%growth_max - Q) <= 1e-12_R8P .and. abs(rep%y1p_min - 5._R8P) <= 1e-10_R8P .and. &
           abs(rep%y1p_max - 5._R8P) <= 1e-10_R8P .and. rep%above == 1._R8P .and. rep%utau_given, &
           '2-D O-grid: growth q, y1+ = 5, all above the target (growth)', rep%growth_max)
call save_block_report(rep, 'fossil_test_block_report.txt', st)
open(newunit=uu, file='fossil_test_block_report.txt', status='old', iostat=st)
read(uu, '(A)') line
close(uu, status='delete')
call check(st == 0 .and. line(1:22) == '# FOSSIL block report ', 'report file written (and removed)', real(st, R8P))

! 2. channel, cell centred
do j = 1, 9
   do i = 1, 9
      cu(i, j) = 0.1_R8P * (i - 1) ; cv(i, j) = H1 * (Q**(j - 1) - 1._R8P) / (Q - 1._R8P)
   enddo
enddo
call b2(1)%initialize(cu, cv)
rec%hi = [9, 1, 1] ; call b2(1)%add_record(rec)
call wall_cutter_2d(b2, c)
call block_report_2d(b2(1), c, 1e-5_R8P, rep, u_tau=0.05_R8P, cell_centred=.true.)
call check(rep%nw == 8 .and. maxval(abs(rep%y1_line - H1 / 2)) <= 1e-15_R8P .and. &
           maxval(abs(rep%y1_wall - H1 / 2)) <= 1e-15_R8P, 'channel, cell centred: y1 = h1 / 2 (wall cells)', &
           real(rep%nw, R8P))

! 3. flat-plate estimate
ut = 10._R8P * sqrt(0.0288_R8P * (10._R8P * 0.45_R8P / 1.5e-5_R8P)**(-0.2_R8P))
call check(abs(utau_flat_plate(10._R8P, 1.5e-5_R8P, 0.45_R8P) - ut) <= 1e-14_R8P * ut, &
           'flat-plate u_tau = u_inf sqrt(0.0288 Re_x^-0.2) (u_tau)', ut)
call block_report_2d(b2(1), c, 1.5e-5_R8P, rep, u_inf=10._R8P, x_le=-0.05_R8P)
call check(rep%utau_estimated .and. .not. rep%utau_given .and. rep%y1p_max > 0._R8P, &
           'the report flags the estimate (max y1+)', rep%y1p_max)

! 4. a folded cell
cu(5, 5) = 0.75_R8P ; cv(5, 5) = 2 * cv(5, 6)
call b2(1)%initialize(cu, cv) ; call b2(1)%add_record(rec)
call block_report_2d(b2(1), c, 1e-5_R8P, rep, u_tau=0.05_R8P)
call check(rep%not_valid >= 1 .and. rep%minj_min <= 0._R8P, 'folded cell: not certified valid (count)', &
           real(rep%not_valid, R8P))

! 5. 3-D O-grid
do j = 1, 21
   do i = 1, 65
      x(i, j, :) = u(i, j) ; y(i, j, :) = v(i, j) ; z(i, j, :) = [0._R8P, 0.1_R8P, 0.2_R8P]
   enddo
enddo
call b3(1)%initialize(x, y, z)
rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [65, 1, 3] ; call b3(1)%add_record(rec)
call wall_cutter(b3, s)
call block_report(b3(1), s, 1e-5_R8P, rep, u_tau=0.05_R8P, status=st)
call check(st == 0 .and. rep%nw == 195 .and. maxval(abs(rep%y1_line - H1)) <= 1e-15_R8P .and. &
           maxval(abs(rep%y1_wall - H1)) <= 1e-12_R8P .and. abs(rep%growth_max - Q) <= 1e-12_R8P, &
           '3-D O-grid: y1 = h1 both ways, growth q (max wall error)', maxval(abs(rep%y1_wall - H1)))
print '(A,ES12.4,A,I0)', '     3-D skewness max ', rep%skew_max, ', cells not certified valid ', rep%not_valid

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
endprogram fossil_test_block_report
