!< FOSSIL test: curve repair (issue #36 step 6).
program fossil_test_curve_repair
!< FOSSIL test: curve repair (issue #36 step 6).
!<
!< - An open square whose end stops 1e-3 short of its start: closed with tolerance 1e-2, valid.
!< - Two open halves of a square, the second starting 1e-4 off the first's end: joined and closed
!<   (2 gaps), valid, area 1 to 1e-3. Same with the second half reversed: joined by its end.
!< - Duplicated vertices (an exact repeat and one 1e-12 away): merged with tolerance 1e-9, 4
!<   vertices left.
!< - An annulus with a counter-clockwise hole: the hole reversed, valid.
!< - Nested islands all wrongly oriented (outer clockwise, hole counter-clockwise, island
!<   clockwise): 3 loops reversed, valid, depths 0, 1, 2.
!< - A 2-vertex loop and a collinear (zero-area) loop: removed.
!< - A gap wider than the tolerance stays open (reported by loops_open).

use fossil, only : curve_object, CURVE_VALID
use penf,   only : I4P, R8P

implicit none

type(curve_object) :: c
integer(I4P)       :: gaps, merged, removed, reversed, nopen
logical            :: ok(7)

! open square with a small gap
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 0.001_R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
call c%repair(tolerance=1.e-2_R8P, gaps_closed=gaps, loops_open=nopen)
ok(1) = gaps == 1 .and. nopen == 0 .and. c%check_validity() == CURVE_VALID .and. c%is_loop_closed(1_I4P)
print '(A,L1)', 'open square closed: ', ok(1)

! two halves, then the second reversed
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P,  1._R8P, 0._R8P, 0._R8P], &
                   v=[0._R8P, 0._R8P, 1._R8P,  1.0001_R8P, 1._R8P, 0.0001_R8P], loop_sizes=[3_I4P, 3_I4P], &
                   closed=[.false., .false.])
call c%repair(tolerance=1.e-3_R8P, gaps_closed=gaps, loops_open=nopen)
ok(2) = gaps == 2 .and. nopen == 0 .and. c%get_loops_number() == 1 .and. c%check_validity() == CURVE_VALID .and. &
        abs(c%get_area() - 1._R8P) < 1.e-3_R8P
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P,  0._R8P, 0._R8P, 1._R8P], &
                   v=[0._R8P, 0._R8P, 1._R8P,  0.0001_R8P, 1._R8P, 1.0001_R8P], loop_sizes=[3_I4P, 3_I4P], &
                   closed=[.false., .false.])
call c%repair(tolerance=1.e-3_R8P, gaps_closed=gaps, loops_open=nopen)
ok(2) = ok(2) .and. gaps == 2 .and. nopen == 0 .and. c%get_loops_number() == 1 .and. c%check_validity() == CURVE_VALID
print '(A,L1)', 'two halves joined (one reversed): ', ok(2)

! duplicated vertices
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P], &
                   v=[0._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P, 1._R8P + 1.e-12_R8P], loop_sizes=[6_I4P])
call c%repair(tolerance=1.e-9_R8P, vertices_merged=merged)
ok(3) = merged == 2 .and. c%get_vertices_number() == 4 .and. c%check_validity() == CURVE_VALID
print '(A,L1)', 'duplicated vertices merged: ', ok(3)

! reversed hole
call c%adopt_loops(u=[0._R8P, 4._R8P, 4._R8P, 0._R8P, 1._R8P, 3._R8P, 3._R8P, 1._R8P], &
                   v=[0._R8P, 0._R8P, 4._R8P, 4._R8P, 1._R8P, 1._R8P, 3._R8P, 3._R8P], loop_sizes=[4_I4P, 4_I4P])
call c%repair(tolerance=0._R8P, loops_reversed=reversed)
ok(4) = reversed == 1 .and. c%check_validity() == CURVE_VALID .and. c%get_area() == 12._R8P
print '(A,L1)', 'reversed hole: ', ok(4)

! nested islands, all wrong
call c%adopt_loops(u=[0._R8P, 0._R8P, 6._R8P, 6._R8P,  1._R8P, 5._R8P, 5._R8P, 1._R8P,  2._R8P, 2._R8P, 4._R8P, 4._R8P], &
                   v=[0._R8P, 6._R8P, 6._R8P, 0._R8P,  1._R8P, 1._R8P, 5._R8P, 5._R8P,  2._R8P, 4._R8P, 4._R8P, 2._R8P], &
                   loop_sizes=[4_I4P, 4_I4P, 4_I4P])
call c%repair(tolerance=0._R8P, loops_reversed=reversed)
ok(5) = reversed == 3 .and. c%check_validity() == CURVE_VALID .and. c%get_area() == 24._R8P .and. &
        c%nesting_depth(3_I4P) == 2
print '(A,L1)', 'nested islands oriented: ', ok(5)

! degenerate loops
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P,  5._R8P, 6._R8P,  7._R8P, 8._R8P, 9._R8P], &
                   v=[0._R8P, 0._R8P, 1._R8P, 1._R8P,  5._R8P, 5._R8P,  7._R8P, 8._R8P, 9._R8P], loop_sizes=[4_I4P, 2_I4P, 3_I4P])
call c%repair(tolerance=0._R8P, loops_removed=removed)
ok(6) = removed == 2 .and. c%get_loops_number() == 1 .and. c%check_validity() == CURVE_VALID
print '(A,L1)', 'degenerate loops removed: ', ok(6)

! gap wider than the tolerance
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 0.5_R8P], &
                   loop_sizes=[5_I4P], closed=[.false.])
call c%repair(tolerance=1.e-2_R8P, gaps_closed=gaps, loops_open=nopen)
ok(7) = gaps == 0 .and. nopen == 1 .and. .not. c%is_loop_closed(1_I4P)
print '(A,L1)', 'wide gap left open: ', ok(7)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_curve_repair
