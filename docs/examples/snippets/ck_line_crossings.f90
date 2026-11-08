program ck_line_crossings
use fossil, only : curve_object
use penf,   only : I4P, R8P
implicit none
type(curve_object)        :: diamond
real(R8P),    allocatable :: t(:)
integer(I4P), allocatable :: segments(:)
logical,      allocatable :: entering(:)
real(R8P)                 :: a
integer(I4P)              :: n, k, j

call diamond%adopt_loops(u=[0.5_R8P, 1._R8P, 0.5_R8P, 0._R8P], v=[0._R8P, 0.5_R8P, 1._R8P, 0.5_R8P], &
                         loop_sizes=[4_I4P])              ! a diamond, counter-clockwise: solid inside
do k = 0, 4                                               ! the lines v = 0, 0.25, ..., 1 along u
   a = 0.25_R8P * k
   call diamond%line_crossings(axis=1_I4P, a=a, t0=-1._R8P, t1=2._R8P, t=t, segments=segments, entering=entering, n=n)
   print '(A,F5.2,A,I0,A)', 'v = ', a, ': ', n, ' crossings'
   do j = 1, n
      print '(A,F7.4,A,I0,A,L1)', '   u = ', t(j), ', segment ', segments(j), ', entering ', entering(j)
   enddo
enddo
endprogram ck_line_crossings
