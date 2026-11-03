!< FOSSIL test: curve_object core (issue #36 step 2).
program fossil_test_curve_core
!< FOSSIL test: curve_object core (issue #36 step 2).
!<
!< - Square [0,2]^2: area 4, perimeter 8, centroid (1,1), valid.
!< - 64-gon inscribed in the unit circle: area (n/2) sin(2 pi/n), perimeter 2 n sin(pi/n).
!< - Annulus (square [0,4]^2 with a clockwise hole [1,3]^2): area 12, centroid (2,2), valid,
!<   hole depth 1; winding numbers 1 in the ring, 0 in the hole and outside.
!< - Island in a hole (0..6 / 1..5 / 2..4): valid, depths 0, 1, 2, area 24.
!< - Defects: the hole oriented counter-clockwise (WRONG_ORIENTATION on loop 2), a 2-vertex loop,
!<   a repeated vertex, a collinear (zero-area) loop.
!< - I/O: the 64-gon saved and reloaded bit for bit; a hand-written file with comments, blank
!<   lines between loops and an explicitly repeated closing vertex loads as the annulus.
!< - Status: mismatched arrays -> STATUS_INVALID_INPUT; missing file -> STATUS_FILE_NOT_FOUND.

use fossil, only : curve_object, CURVE_VALID, CURVE_TOO_FEW_VERTICES, CURVE_ZERO_LENGTH_SEGMENT, CURVE_ZERO_AREA_LOOP, &
                   CURVE_WRONG_ORIENTATION, STATUS_OK, STATUS_INVALID_INPUT, STATUS_FILE_NOT_FOUND
use penf,   only : I4P, R8P

implicit none

character(*), parameter :: F1 = 'fossil_test_curve_core_1.dat', F2 = 'fossil_test_curve_core_2.dat'
real(R8P),    parameter :: PI = acos(-1._R8P)
type(curve_object)      :: c, r
real(R8P)               :: u(64), v(64)
integer(I4P)            :: i, n, st, l, unit
logical                 :: ok(8)

! square
call c%adopt_loops(u=[0._R8P, 2._R8P, 2._R8P, 0._R8P], v=[0._R8P, 0._R8P, 2._R8P, 2._R8P], loop_sizes=[4_I4P], status=st)
ok(1) = st == STATUS_OK .and. c%get_area() == 4._R8P .and. c%get_perimeter() == 8._R8P .and. &
        all(c%get_centroid() == [1._R8P, 1._R8P]) .and. c%check_validity() == CURVE_VALID .and. c%get_normal_axis() == 3
print '(A,L1)', 'square: ', ok(1)

! 64-gon
n = 64
do i = 1, n
   u(i) = cos(2._R8P * PI * (i - 1) / n) ; v(i) = sin(2._R8P * PI * (i - 1) / n)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[n], normal_axis=1_I4P)
ok(2) = abs(c%get_area() - 0.5_R8P * n * sin(2._R8P * PI / n)) < 1.e-14_R8P .and. &
        abs(c%get_perimeter() - 2._R8P * n * sin(PI / n)) < 1.e-14_R8P .and. all(abs(c%get_centroid()) < 1.e-15_R8P) .and. &
        c%check_validity() == CURVE_VALID .and. c%get_normal_axis() == 1
print '(A,L1,ES12.3)', '64-gon: ', ok(2), c%get_area() - 0.5_R8P * n * sin(2._R8P * PI / n)

! I/O round trip of the 64-gon
call c%save_into_file(file_name=F1, status=st)
call r%load_from_file(file_name=F1, normal_axis=1_I4P, status=st)
ok(3) = st == STATUS_OK .and. r%get_vertices_number() == n .and. r%get_loops_number() == 1
if (ok(3)) then
   do i = 1, n
      ok(3) = ok(3) .and. all(r%get_vertex(i) == c%get_vertex(i))
   enddo
endif
print '(A,L1)', 'save/load bit for bit: ', ok(3)

! annulus
call annulus(hole_ccw=.false., cur=c)
ok(4) = c%get_area() == 12._R8P .and. all(c%get_centroid() == [2._R8P, 2._R8P]) .and. c%check_validity() == CURVE_VALID .and. &
        c%nesting_depth(1_I4P) == 0 .and. c%nesting_depth(2_I4P) == 1 .and. c%get_loop_area(2_I4P) == -4._R8P .and. &
        c%winding_number_exact([0.5_R8P, 0.5_R8P]) == 1 .and. c%winding_number_exact([2._R8P, 2._R8P]) == 0 .and. &
        c%winding_number_exact([5._R8P, 5._R8P]) == 0 .and. c%winding_number_exact([3.5_R8P, 2._R8P]) == 1
print '(A,L1)', 'annulus: ', ok(4)

! hand-written file: comments, blank lines, explicit closing vertex
open(newunit=unit, file=F2, status='replace', action='write')
write(unit, '(A)') '# annulus, hand written', '0 0', '4 0', '4 4  # corner', '0 4', '0 0', '', '', &
                   '# hole, clockwise', '1 1', '1 3', '3 3', '3 1'
close(unit)
call r%load_from_file(file_name=F2, status=st)
ok(5) = st == STATUS_OK .and. r%get_loops_number() == 2 .and. r%get_vertices_number() == 8 .and. r%get_area() == 12._R8P .and. &
        r%check_validity() == CURVE_VALID
print '(A,L1)', 'hand-written file: ', ok(5)

! island in a hole
call c%adopt_loops(u=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P,  2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   v=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P,  2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
ok(6) = c%check_validity() == CURVE_VALID .and. c%nesting_depth(3_I4P) == 2 .and. c%get_area() == 24._R8P .and. &
        c%winding_number_exact([3._R8P, 3._R8P]) == 1 .and. c%winding_number_exact([1.5_R8P, 3._R8P]) == 0
print '(A,L1)', 'island in a hole: ', ok(6)

! defects
ok(7) = .true.
call annulus(hole_ccw=.true., cur=c)
ok(7) = ok(7) .and. c%check_validity(loop=l) == CURVE_WRONG_ORIENTATION .and. l == 2
call c%adopt_loops(u=[0._R8P, 1._R8P], v=[0._R8P, 0._R8P], loop_sizes=[2_I4P])
ok(7) = ok(7) .and. c%check_validity() == CURVE_TOO_FEW_VERTICES
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 0._R8P, 1._R8P], loop_sizes=[4_I4P])
ok(7) = ok(7) .and. c%check_validity() == CURVE_ZERO_LENGTH_SEGMENT
call c%adopt_loops(u=[0._R8P, 1._R8P, 2._R8P], v=[0._R8P, 1._R8P, 2._R8P], loop_sizes=[3_I4P])
ok(7) = ok(7) .and. c%check_validity() == CURVE_ZERO_AREA_LOOP
print '(A,L1)', 'defects detected: ', ok(7)

! status
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P], v=[0._R8P, 0._R8P], loop_sizes=[3_I4P], status=st)
ok(8) = st == STATUS_INVALID_INPUT .and. c%get_vertices_number() == 0
call c%load_from_file(file_name='fossil_test_curve_core_missing.dat', status=st)
ok(8) = ok(8) .and. st == STATUS_FILE_NOT_FOUND
print '(A,L1)', 'status codes: ', ok(8)

open(newunit=unit, file=F1, status='old') ; close(unit, status='delete')
open(newunit=unit, file=F2, status='old') ; close(unit, status='delete')
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine annulus(hole_ccw, cur)
   !< Square [0,4]^2 with the hole [1,3]^2, clockwise unless hole_ccw.
   logical,            intent(in)    :: hole_ccw !< Wrong hole orientation.
   type(curve_object), intent(inout) :: cur      !< Curve.

   if (hole_ccw) then
      call cur%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], &
                           v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], loop_sizes=[4_I4P, 4_I4P])
   else
      call cur%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], &
                           v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
   endif
   endsubroutine annulus
endprogram fossil_test_curve_core
