!< Figure of the cookbook recipe "Check the first-cell height of a 3-D wall block": y1+ at the 41 wall points of the
!< mid-span line k = 3 of ck_block_report_3d, drawn as bars along the plate (height 0.3 per unit of y1+), with the
!< target y1+ = 1 as a line: above it near the leading edge, below it downstream, where u_tau falls.
program fig_ck_block_report_3d
use fossil, only : surface_stl_object, curvilinear_block_t, topo_record_t, TOPO_BOUNDARY, wall_cutter, &
                   geometric_spacing, block_report_t, block_report, utau_flat_plate
use fossil_fig, only : fig_open, fig_segments, fig_segment, fig_points, fig_label, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter      :: NU = 1.5e-5_R8P, U_INF = 10._R8P   ! air at 10 m/s
type(curvilinear_block_t) :: plate(1)
type(topo_record_t)       :: rec
type(surface_stl_object)  :: cutter
type(block_report_t)      :: report
real(R8P)                 :: x(41, 31, 5), y(41, 31, 5), z(41, 31, 5), h(30), h1
integer(I4P)              :: i, j, k, status, u
type(vector_R8P), allocatable :: a(:), b(:)
character(90)             :: caption

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
allocate(a(0), b(0))
do i = 1, report%nw
   if (report%point(3, i) /= 3) cycle
   a = [a, vector_R8P(x(report%point(1, i), 1, 1), 0._R8P, 0._R8P)]
   b = [b, vector_R8P(x(report%point(1, i), 1, 1), 0.3_R8P * report%y1p(i), 0._R8P)]
enddo
write(caption, '(A,F6.4,A,F6.4,A)') 'y1+ along the plate, from ', report%y1p_max, ' to ', report%y1p_min, &
                                    '; the line: target 1'
call fig_open(u, 'ck-block-report-3d', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
call fig_segments(u, 'probe', a, b)
call fig_points(u, 'foot', b)
call fig_segment(u, 'hole', vector_R8P(0.03_R8P, 0.3_R8P, 0._R8P), vector_R8P(1.07_R8P, 0.3_R8P, 0._R8P))
call fig_segment(u, 'body', vector_R8P(0.03_R8P, 0._R8P, 0._R8P), vector_R8P(1.07_R8P, 0._R8P, 0._R8P))
call fig_label(u, 'ghost', vector_R8P(0.05_R8P, -0.04_R8P, 0._R8P), 'x = 0.05')
call fig_label(u, 'ghost', vector_R8P(0.95_R8P, -0.04_R8P, 0._R8P), 'x = 1.05')
call fig_close(u)
endprogram fig_ck_block_report_3d
