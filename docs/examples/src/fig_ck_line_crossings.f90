!run figck-line-crossings fig_ck_line_crossings
!fig ck-line-crossings
!< Figure of the cookbook recipe "Find where a grid line crosses a curve": the diamond of ck_line_crossings and the five
!< grid lines v = 0, 0.25, ..., 1; each crossing line_crossings returns, teal where the line enters the solid and ember
!< where it leaves. v = 0 touches the bottom vertex and gets two crossings there (enter, leave); v = 1 touches the top
!< vertex and gets none.
program fig_ck_line_crossings
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(curve_object)            :: diamond
real(R8P),    allocatable     :: t(:)
integer(I4P), allocatable     :: segments(:)
logical,      allocatable     :: entering(:)
type(vector_R8P)              :: a(4), b(4)
type(vector_R8P), allocatable :: pin(:), pout(:)
real(R8P)                     :: v
integer(I4P)                  :: n, k, j, u
character(30)                 :: text

call diamond%adopt_loops(u=[0.5_R8P, 1._R8P, 0.5_R8P, 0._R8P], v=[0._R8P, 0.5_R8P, 1._R8P, 0.5_R8P], &
                         loop_sizes=[4_I4P])
a = [vector_R8P(0.5_R8P, 0._R8P, 0._R8P), vector_R8P(1._R8P, 0.5_R8P, 0._R8P), vector_R8P(0.5_R8P, 1._R8P, 0._R8P), &
     vector_R8P(0._R8P, 0.5_R8P, 0._R8P)]
b = [a(2), a(3), a(4), a(1)]
call fig_open(u, 'ck-line-crossings', 'line_crossings: entering the solid (teal), leaving it (ember)', &
              -90._R8P, 90._R8P, height=420_I4P)
call fig_segments(u, 'body', a, b)
allocate(pin(0), pout(0))
do k = 0, 4
   v = 0.25_R8P * k
   call diamond%line_crossings(axis=1_I4P, a=v, t0=-1._R8P, t1=2._R8P, t=t, segments=segments, entering=entering, n=n)
   call fig_segment(u, 'box', vector_R8P(-0.65_R8P, v, 0._R8P), vector_R8P(1.25_R8P, v, 0._R8P))
   write(text, '(A,F4.2,A,I0,A)') 'v = ', v, ': ', n, ' crossings'
   call fig_label(u, 'ghost', vector_R8P(-0.62_R8P, v + 0.02_R8P, 0._R8P), trim(text))
   do j = 1, n
      if (entering(j)) then ; pin = [pin, vector_R8P(t(j), v, 0._R8P)] ; else ; pout = [pout, vector_R8P(t(j), v, 0._R8P)]
      endif
   enddo
enddo
call fig_points(u, 'outside', pout)
call fig_points(u, 'inside', pin)
call fig_close(u)
endprogram fig_ck_line_crossings
