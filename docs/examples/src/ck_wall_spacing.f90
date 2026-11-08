!run ck_wall_spacing ck_wall_spacing
program ck_wall_spacing
use fossil, only : curve_object, curvilinear_block_2d_t, topo_record_t, TOPO_BOUNDARY, hyperbolic_grid_2d, &
                   geometric_spacing, wall_cutter_2d, block_report_t, block_report_2d
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P), NU = 1.5e-5_R8P, UTAU = 0.4_R8P
type(curvilinear_block_2d_t) :: grid(1)
type(curve_object)           :: wall
type(topo_record_t)          :: rec
type(block_report_t)         :: report
real(R8P)                    :: u(129), v(129), xc, yt, h1
integer(I4P)                 :: i, invalid, status

do i = 1, 129 ! NACA0012, chord 1: upper side from the trailing edge, then the lower side (counter-clockwise)
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / 129._R8P))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= 65)
enddo
h1 = 1._R8P * NU / UTAU                                                   ! y1+ = 1 with this u_tau
call hyperbolic_grid_2d(u, v, geometric_spacing(h1, 1.2_R8P, 40_I4P), grid(1), invalid=invalid, status=status)
print '(A,I0,A,I0,A,I0)', 'O-grid ', grid(1)%ni, ' x ', grid(1)%nj, ', cells not certified valid ', invalid
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)
call wall_cutter_2d(grid, wall)
call block_report_2d(grid(1), wall, NU, report, u_tau=UTAU, target_y1p=1._R8P)
print '(A,3F8.4)', 'y1+ min mean max ', report%y1p_min, report%y1p_mean, report%y1p_max
print '(A,F7.4,A,F7.4)', 'growth ratio max ', report%growth_max, ', skewness max ', report%skew_max
endprogram ck_wall_spacing
