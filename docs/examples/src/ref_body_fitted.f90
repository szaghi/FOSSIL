!run ref_body_fitted ref_body_fitted
program ref_body_fitted
use fossil,              only : surface_stl_object, projection_target_t, classify_points, project_points, &
                                PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER, curve_object, curvilinear_block_2d_t, &
                                topo_record_t, TOPO_BOUNDARY, hyperbolic_grid_2d, geometric_spacing, wall_cutter_2d, &
                                block_report_t, block_report_2d, save_block_report, utau_flat_plate
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
implicit none
character(6), parameter      :: KIND_NAME(4) = ['free  ', 'patch ', 'curve ', 'corner']
real(R8P), parameter         :: PI = acos(-1._R8P), NU = 1.5e-5_R8P, U_INF = 10._R8P
type(surface_stl_object)     :: wing
type(facet_object), pointer  :: f
type(projection_target_t)    :: target
type(curvilinear_block_2d_t) :: grid(1)
type(curve_object)           :: wall
type(topo_record_t)          :: rec
type(block_report_t)         :: report
integer(I4P), allocatable    :: flap(:)
real(R8P)                    :: px(4), py(4), pz(4), fx(4), fy(4), fz(4), d(4), gx(4), gy(4), gz(4), e(4)
real(R8P)                    :: u(129), v(129), xc, yt
integer(I4P)                 :: kinds(4), ids(4), i, invalid, status
character(160)               :: line

!region target
call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
allocate(flap(0))
do i = 1, wing%get_facets_number() ! the rear facets become the patch 'flap'
   f => wing%facet_at(i)
   if (min(f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x) >= 0.45_R8P) flap = [flap, i]
enddo
call wing%set_patch(flap, 'flap')
call wing%projection_target(target) ! sharp edges above 30 degrees, patch boundaries, corners
print '(A,I0,A,I0,A,I0)', 'target: patches ', target%np, ', constraint curves ', target%nc, ', corners ', target%nk
!endregion target

!region project
px = [0.745_R8P, 0.74_R8P, 0.57_R8P, 0.2_R8P]  ! near the trailing-edge corner at the root, near the trailing edge,
py = [0.004_R8P, 0.004_R8P, 0.03_R8P, 0.07_R8P] ! near the upper patch boundary, above the upper surface
pz = [0.01_R8P, 0.5_R8P, 0.4_R8P, 0.5_R8P]
call classify_points(target, px, py, pz, tolerance=0.02_R8P, kinds=kinds, ids=ids)
call project_points(target, px, py, pz, kinds, ids, fx, fy, fz, distances=d)
call project_points(target, px, py, pz, [(PROJ_FREE, i = 1, 4)], [(0_I4P, i = 1, 4)], gx, gy, gz, distances=e)
do i = 1, 4
   print '(A,I0,3A,I0,A,3F9.5,A,F8.5)', 'point ', i, ': ', trim(KIND_NAME(kinds(i) + 1)), ' ', ids(i), ', foot', fx(i), fy(i), &
                                       fz(i), ', moved ', d(i)
   print '(A,3F9.5,A,F8.5)', '         free foot          ', gx(i), gy(i), gz(i), ', moved ', e(i)
enddo
! a point above the main patch, constrained to the flap
call project_points(target, [0.45_R8P], [0.05_R8P], [0.5_R8P], [PROJ_PATCH], [2_I4P], fx(1:1), fy(1:1), fz(1:1), &
                    distances=d(1:1))
print '(A,3F9.5,A,F8.5)', '(0.45, 0.05, 0.5) on the flap ', fx(1), fy(1), fz(1), ', moved ', d(1)
!endregion project

!region report
do i = 1, 129 ! a NACA 0012 section, chord 1, counter-clockwise from the trailing edge
   xc = 0.5_R8P * (1._R8P + cos(2 * PI * (i - 1) / 129._R8P))
   yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
   u(i) = xc ; v(i) = merge(yt, -yt, i <= 65)
enddo
call hyperbolic_grid_2d(u, v, geometric_spacing(2.e-5_R8P, 1.15_R8P, 50_I4P), grid(1), invalid=invalid, status=status)
print '(A,I0,A,I0,A,I0,A,I0)', 'O-grid ', grid(1)%ni, ' x ', grid(1)%nj, ', status ', status, ', not valid ', invalid
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)
call wall_cutter_2d(grid, wall)
print '(A,F8.5)', 'flat-plate u_tau at mid-chord ', utau_flat_plate(U_INF, NU, 0.5_R8P)
call block_report_2d(grid(1), wall, NU, report, u_inf=U_INF, x_le=0._R8P, target_y1p=1._R8P, cell_centred=.true.)
print '(A,I0,A,L1,A,L1)', 'wall cells ', report%nw, ', cell centred ', report%cell_centred, ', u_tau estimated ', &
                          report%utau_estimated
print '(A,3F8.4,A,F6.3)', 'y1+ min mean max ', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', fraction above 1 ', report%above
i = maxloc(abs(report%y1_line - report%y1_wall), 1)
print '(A,ES10.3,A,ES10.3,A,I0)', 'y1 along the line ', report%y1_line(i), ', to the wall ', report%y1_wall(i), &
                                  ' at wall cell ', i
print '(A,F7.4,A,F7.4,A,I0)', 'growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                              ', cells not valid ', report%not_valid
call save_block_report(report, 'naca_report.txt')
open(newunit=i, file='naca_report.txt', action='read')
read(i, '(A)') line ; print '(A)', trim(line)
read(i, '(A)') line ; print '(A)', trim(line)
read(i, '(A)') line ; print '(A)', trim(line)
close(i)
!endregion report
endprogram ref_body_fitted
