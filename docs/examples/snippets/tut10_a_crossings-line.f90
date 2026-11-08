use fossil,     only : curve_object, STATUS_OK
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(curve_object)        :: section
real(R8P)                 :: x(49), solid
real(R8P),    allocatable :: t(:)
integer(I4P), allocatable :: segments(:)
logical,      allocatable :: entering(:)
integer(I4P)              :: i, j, best, n, nbest, status

call section%load_from_file(file_name='dragon-section.dat', status=status)       ! chapter 7
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]          ! the cell edges of a 48 x 48 block, h = 0.025
best = 1 ; nbest = 0
do j = 1, 49                                          ! the grid line x = x(j), along y, with the most crossings
   call section%line_crossings(axis=2_I4P, a=x(j), t0=x(1), t1=x(49), t=t, segments=segments, entering=entering, n=n)
   if (n > nbest) then
      best = j ; nbest = n
   endif
enddo
call section%line_crossings(axis=2_I4P, a=x(best), t0=x(1), t1=x(49), t=t, segments=segments, entering=entering, n=n)
print '(A,I0,A,F7.4,A,I0,A)', 'grid line ', best, ' (x = ', x(best), '): ', n, ' crossings'
do i = 1, n
   print '(A,F9.5,A,I4,A,A)', '  y = ', t(i), '  segment ', segments(i), '  ', merge('enters', 'leaves', entering(i))
enddo
solid = sum(t(2:n:2) - t(1:n-1:2))                    ! in +y: enter, leave, enter, leave, ...
print '(A,F9.6)', 'solid length of the line: ', solid
