!< Figures of the reference page "Marching cubes": the unit sphere extract_isosurface draws from its analytic signed
!< distance on the page's 64^3 grid over [-2, 2]^3, and a corner of the unit cube resample_via_distance_field
!< re-tessellates at resolution 64 (zoomed: the whole has ~40k facets, only those within 0.75 of the corner are drawn).
program fig_ref_marching_cubes
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, SDF_STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_triangle, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 64
type(surface_stl_object)        :: sphere, cube_in, cube_resampled
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: corner
real(R8P),          allocatable :: values(:, :, :)
real(R8P)                       :: dx
character(80)                   :: caption
integer(I4P)                    :: u, i, j, k, status

allocate(values(NG, NG, NG))                                                       ! the page's sphere
dx = 4._R8P / real(NG - 1, R8P)
do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   values(i, j, k) = norm2(-2._R8P + dx * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=values, bmin=vector_R8P(-2._R8P, -2._R8P, -2._R8P), &
                        bmax=vector_R8P(2._R8P, 2._R8P, 2._R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
write(caption, '(A,I0,A,F5.3,A)') 'extract_isosurface, 64^3: ', sphere%get_facets_number(), ' facets, volume ', &
                                  sphere%get_volume(), ' (4/3 pi = 4.189)'
call fig_open(u, 'ref-mc-sphere', trim(caption), -60._R8P, 25._R8P)
call fig_surface(u, 'body', sphere)
call fig_close(u)

call cube_in%load_from_file(file_name='cube-01.stl', guess_format=.true.)          ! the page's roundtrip
call cube_in%resample_via_distance_field(resolution=64_I4P, surface_out=cube_resampled, status=status)
if (status /= SDF_STATUS_OK) error stop 'resample failed'
write(caption, '(A,I0,A,F7.5,A)') 'resample_via_distance_field(64): ', cube_resampled%get_facets_number(), &
                                  ' facets, volume ', cube_resampled%get_volume(), ', a corner'
call fig_open(u, 'ref-mc-cube-resampled', trim(caption), -60._R8P, 25._R8P)
corner = vector_R8P(1._R8P, 0._R8P, 1._R8P)
call fig_zoom(u, corner, 0.2_R8P)
do i = 1, cube_resampled%get_facets_number()                                       ! the facets near the corner
   f => cube_resampled%facet_at(i)
   if (maxval(abs([f%centroid%x - corner%x, f%centroid%y - corner%y, f%centroid%z - corner%z])) < 0.75_R8P) &
      call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_close(u)
endprogram fig_ref_marching_cubes
