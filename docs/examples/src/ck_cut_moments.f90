!run ck_cut_moments ck_cut_moments
program ck_cut_moments
use fossil, only : curve_object, cut_cells_2d_t
use penf,   only : I4P, R8P
implicit none
type(curve_object)   :: disc
type(cut_cells_2d_t) :: cut
real(R8P)            :: u(48), v(48), x(17), sm(3), c(2), area_err, xx, xx_piece
integer(I4P)         :: i, p, q

u = [(0.5_R8P + 0.3_R8P * cos(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]   ! a disc, 48 segments
v = [(0.5_R8P + 0.3_R8P * sin(2._R8P * acos(-1._R8P) * i / 48), i = 0, 47)]
call disc%adopt_loops(u=u, v=v, loop_sizes=[48_I4P])
x = [(0.0625_R8P * i, i = 0, 16)]                                               ! 16 x 16 cells, h = 0.0625
call disc%cut_cells(x=x, y=x, cut=cut, moments_order=2_I4P, quadrature=.true.)
area_err = 0._R8P ; xx = 0._R8P
do p = 1, cut%np                                                                ! every fluid piece
   area_err = max(area_err, abs(sum(cut%quad_weight(cut%quad_start(p):cut%quad_start(p + 1) - 1)) - cut%piece_area(p)))
   do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
      xx = xx + cut%quad_weight(q) * cut%quad_point(1, q) ** 2                  ! int x^2 over the fluid
   enddo
enddo
call disc%mass_properties(second_moments=sm, about=[0._R8P, 0._R8P])            ! int x^2 over the disc
print '(A,I0,A,ES9.2)', 'pieces ', cut%np, ', max |sum of weights - area| ', area_err
print '(A,F16.13,A,F16.13)', 'int x^2 over the fluid: quadrature ', xx, ', exact 1/3 - disc ', 1._R8P / 3 - sm(1)
p = cut%first_piece(4, 8)                                                       ! a cut cell
c = [0.5_R8P * (x(4) + x(5)), 0.5_R8P * (x(8) + x(9))]
xx_piece = 0._R8P
do q = cut%quad_start(p), cut%quad_start(p + 1) - 1
   xx_piece = xx_piece + cut%quad_weight(q) * ((cut%quad_point(1, q) - c(1)) / 0.0625_R8P) ** 2
enddo
print '(A,F8.5,A,ES20.13,A,ES20.13)', 'cell (4, 8): alpha ', cut%alpha(4, 8), ', moment X^2 ', cut%piece_moments(4, p), &
                                      ', by quadrature ', xx_piece
endprogram ck_cut_moments
