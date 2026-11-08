use fossil, only : curve_object, CURVE_VALID, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter :: PI = acos(-1._R8P)
type(curve_object)   :: c
real(R8P)            :: u(16), v(16), d, foot(2), p(2, 3)
integer(I4P)         :: k, l, defect, seg, vtx, gaps, reversed, still_open

u(1:4) = [0._R8P, 2._R8P, 2._R8P, 0._R8P] ; v(1:4) = [0._R8P, 0._R8P, 1._R8P, 1._R8P]  ! outer plate, counter-clockwise
do k = 0, 7                                                                          ! hole, an octagon: counter-clockwise,
   u(5 + k) = 0.6_R8P + 0.25_R8P * cos(PI * k / 4) ; v(5 + k) = 0.5_R8P + 0.25_R8P * sin(PI * k / 4) ! the wrong way
enddo
u(13:16) = [3._R8P, 4._R8P, 3._R8P, 3._R8P] ; v(13:16) = [0._R8P, 0._R8P, 1._R8P, 1.e-7_R8P] ! a triangle that misses
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, 8_I4P, 4_I4P], closed=[.true., .true., .false.]) ! closing: open chain
