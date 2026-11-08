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
