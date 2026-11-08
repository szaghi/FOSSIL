use fossil,     only : surface_stl_object, curve_object, projection_target_t, project_points, classify_points, &
                       slice_surface, PROJ_FREE, PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter   :: N = 6
type(surface_stl_object)  :: wing
type(projection_target_t) :: target
type(curve_object)        :: section
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                 :: px(N), py(N), pz(N), ffx(N), ffy(N), ffz(N), fx(N), fy(N), fz(N)
integer(I4P)              :: kinds(N), ids(N), free(N), i, k, u
character(8)              :: num
character(16), parameter  :: NAME(0:3) = [character(16) :: 'PROJ_FREE', 'PROJ_PATCH', 'PROJ_CURVE', 'PROJ_CORNER']

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.) ! chord x in [-0.25, 0.75], span z in [0, 1]
call wing%projection_target(target)                                          ! built once: patches, curves, corners
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'wing: ', wing%get_facets_number(), ' facets; target: patches ', target%np, &
      ', constraint curves ', target%nc, ' (', size(target%seg, 2), ' segments), corners ', target%nk
