!run tut15 tut15_body_fitted
!< Tutorial 15: body-fitted grids on the NACA0012 wing: a constrained projection of a block face onto the wing, then a
!< hyperbolic O-grid around its section with the first height set from a y1+ target, and the wall-spacing report.
program tut15_body_fitted
!region target
use fossil, only : surface_stl_object, curve_object, projection_target_t, curvilinear_block_t, curvilinear_block_2d_t, &
                   topo_record_t, TOPO_BOUNDARY, project_points, classify_points, project_block_face, slice_surface, &
                   PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER, hyperbolic_grid_2d, geometric_spacing, &
                   wall_cutter_2d, block_report_t, block_report_2d, save_block_report, utau_flat_plate
use penf,   only : I4P, R8P

implicit none
real(R8P),    parameter      :: PI = acos(-1._R8P)
type(surface_stl_object)     :: wing
type(projection_target_t)    :: target
type(curvilinear_block_t)    :: block
type(curve_object)           :: section, wall
type(curvilinear_block_2d_t) :: grid(1)
type(topo_record_t)          :: rec
type(block_report_t)         :: report
real(R8P)                    :: x(41, 21, 2), y(41, 21, 2), z(41, 21, 2), xc, yt
real(R8P)                    :: px(3), py(3), pz(3), fx(3), fy(3), fz(3)
real(R8P),    allocatable    :: moved(:,:), spacing(:), u(:), v(:)
integer(I4P), allocatable    :: kinds(:,:)
integer(I4P)                 :: pk(3), pid(3), i, j, invalid, status, unit, first, last, m
real(R8P)                    :: nu, u_inf, u_tau, h1, a(2), b(2)
character(200)               :: line

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
call wing%projection_target(target)                                          ! patches, constraint curves, corners
print '(A,I0,A,I0,A,I0,A,I0)', 'wing: ', wing%get_facets_number(), ' facets, patches ', target%np, &
                               ', constraint curves ', target%nc, ', corners ', target%nk
!endregion target
!region points
px = 0.75_R8P ; py = 0.00626_R8P ; pz = [0.1_R8P, 0.5_R8P, 0.9_R8P]  ! 0.005 above the upper trailing-edge line
call classify_points(target, px, py, pz, tolerance=0.02_R8P, kinds=pk, ids=pid)
print '(A,3I2)', 'classified kinds (2 = PROJ_CURVE): ', pk
call project_points(target, px, py, pz, [PROJ_FREE, PROJ_FREE, PROJ_FREE], [0_I4P, 0_I4P, 0_I4P], fx, fy, fz)
print '(A,3(" (",F7.4,",",F8.5,")"))', 'free feet x, y:        ', (fx(i), fy(i), i = 1, 3)
call project_points(target, px, py, pz, pk, pid, fx, fy, fz)
print '(A,3(" (",F7.4,",",F8.5,")"))', 'constrained feet x, y: ', (fx(i), fy(i), i = 1, 3)
!endregion points
!region face
do j = 1, 21                                 ! a 41 x 21 grid 0.005 above the upper side, cosine chordwise
   do i = 1, 41
      xc = 0.5_R8P * (1._R8P - cos(PI * (i - 1) / 40._R8P))
      yt = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.126_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1015_R8P * xc**4)
      x(i, j, :) = xc - 0.25_R8P ; y(i, j, :) = [yt + 0.005_R8P, yt + 0.1_R8P] ; z(i, j, :) = (j - 1) / 20._R8P
   enddo
enddo
call block%initialize(x, y, z)
call project_block_face(block, face=5_I4P, target=target, tolerance=0.01_R8P, moved=moved, kinds=kinds, status=status)
print '(A,I0,A,I0,A,I0,A,F7.4)', 'face: corners ', count(kinds == PROJ_CORNER), ', on curves ', &
      count(kinds == PROJ_CURVE), ', on patches ', count(kinds == PROJ_PATCH), ', largest move ', maxval(moved)
print '(A,L1,A,L1)', 'trailing-edge row on x = 0.75: ', all(block%x(41, :, 1) == 0.75_R8P), &
                     ', root and tip rows on z = 0, 1: ', all(block%z(:, 1, 1) == 0._R8P) .and. all(block%z(:, 21, 1) == 1._R8P)
!endregion face
!region ogrid
call slice_surface(wing, normal_axis=3_I4P, level=0.5_R8P, curve=section)   ! the mid-span section, (u, v) = (x, y)
call section%get_loop_range(1_I4P, first, last)
allocate(u(2 * (last - first + 1)), v(2 * (last - first + 1)))
do i = first, last                                                      ! each segment split in 2: the same polygon
   a = section%get_vertex(i) ; b = section%get_vertex(merge(first, i + 1_I4P, i == last))
   do j = 0, 1
      u(2 * (i - first) + j + 1) = a(1) + (b(1) - a(1)) * j / 2._R8P
      v(2 * (i - first) + j + 1) = a(2) + (b(2) - a(2)) * j / 2._R8P
   enddo
enddo
print '(A,I0,A,I0,A)', 'section: ', last - first + 1, ' vertices, ', size(u), ' wall points'
nu = 1.e-6_R8P ; u_inf = 1._R8P                                         ! chord 1, Re = u_inf c / nu = 1e6
u_tau = utau_flat_plate(u_inf, nu, 0.5_R8P)                             ! flat-plate estimate at mid-chord
h1 = 1._R8P * nu / u_tau                                                ! y1+ = 1 by design: y1 = y1+ nu / u_tau
spacing = geometric_spacing(h1, 1.2_R8P, 40_I4P)                        ! 40 layers, ratio 1.2
call hyperbolic_grid_2d(u, v, spacing, grid(1), invalid=invalid, status=status)
print '(A,ES10.3,A,ES10.3,A,F7.3)', 'u_tau ', u_tau, ', h1 ', h1, ', outer height ', sum(spacing)
print '(A,I0,A,I0,A,I0,A,I0)', 'O-grid ', grid(1)%ni, ' x ', grid(1)%nj, ', status ', status, &
                               ', cells not certified valid ', invalid
!endregion ogrid
!region report
rec%kind = TOPO_BOUNDARY ; rec%face = 3 ; rec%lo = [1, 1, 1] ; rec%hi = [grid(1)%ni, 1, 1] ; rec%name = 'wall'
call grid(1)%add_record(rec)                                            ! side j = 1 is the wall
call wall_cutter_2d(grid, wall)                                         ! the wall as a curve
call block_report_2d(grid(1), wall, nu, report, u_tau=u_tau, target_y1p=1._R8P)
print '(A,3F8.4,A,F6.3)', 'u_tau given:      y1+ min mean max', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
call block_report_2d(grid(1), wall, nu, report, u_inf=u_inf, x_le=-0.25_R8P, target_y1p=1._R8P)
print '(A,3F8.4,A,F6.3)', 'flat-plate u_tau: y1+ min mean max', report%y1p_min, report%y1p_mean, report%y1p_max, &
                          ', above target ', report%above
m = maxloc(report%y1p, dim=1)
print '(A,F7.4,A,I0,A,2F9.5,A)', '  max ', report%y1p(m), ' at wall point ', report%point(1, m), ', (x, y) = (', &
                                 grid(1)%u(report%point(1, m), 1), grid(1)%v(report%point(1, m), 1), ')'
print '(A,F7.4,A,F7.4,A,F7.4,A,I0)', 'growth max ', report%growth_max, ', skewness max ', report%skew_max, &
                                     ', orthogonality min ', report%orth_min, ', not valid ', report%not_valid
call save_block_report(report, 'tut15-report.txt')
open(newunit=unit, file='tut15-report.txt', action='read')
do i = 1, 7
   read(unit, '(A)') line
   if (i /= 5) print '(A)', trim(line)                                  ! line 5, the cell quality, is long
enddo
close(unit)
!endregion report
endprogram tut15_body_fitted
