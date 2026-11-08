program ck_block_report_3d
use fossil, only : surface_stl_object, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, wall_cutter, &
                   geometric_spacing, block_report_t, block_report, utau_flat_plate
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter      :: NU = 1.5e-5_R8P, U_INF = 10._R8P   ! air at 10 m/s
type(curvilinear_block_t) :: plate(1)
type(topo_record_t)       :: rec
type(surface_stl_object)  :: cutter
type(block_report_t)      :: report
real(R8P)                 :: x(41, 31, 5), y(41, 31, 5), z(41, 31, 5), h(30), h1
integer(I4P)              :: i, j, k, status

h1 = NU / utau_flat_plate(U_INF, NU, 0.5_R8P)               ! y1+ = 1 at mid-plate, by the flat-plate estimate
h = geometric_spacing(h1, 1.2_R8P, 30_I4P)
do k = 1, 5 ; do j = 1, 31 ; do i = 1, 41                   ! a block over the plate 0.05 <= x <= 1.05, 0 <= z <= 0.1
   x(i, j, k) = 0.05_R8P + 0.025_R8P * (i - 1) ; y(i, j, k) = sum(h(1:j - 1)) ; z(i, j, k) = 0.025_R8P * (k - 1)
enddo ; enddo ; enddo
call plate(1)%initialize(x, y, z)
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [41, 1, 5] ; rec%name = 'wall'
call plate(1)%add_record(rec)                               ! j = 1: the wall
call wall_cutter(plate, cutter)
call block_report(plate(1), cutter, NU, report, u_inf=U_INF, x_le=0._R8P, target_y1p=1._R8P, status=status)
if (status /= 0) error stop 'block_report failed'
print '(A,ES10.3,A,I0,A,L1)', 'first height ', h1, ', wall points ', report%nw, ', u_tau estimated ', &
                              report%utau_estimated
print '(A,3F8.4)', 'y1+ min mean max ', report%y1p_min, report%y1p_mean, report%y1p_max
print '(A,F7.4)', 'fraction of wall points above y1+ = 1: ', report%above
print '(A,F7.4,A,I0)', 'growth ratio max ', report%growth_max, ', cells not certified valid ', report%not_valid
endprogram ck_block_report_3d
