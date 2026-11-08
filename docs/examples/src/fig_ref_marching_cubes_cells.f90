!run figref-marching-cubes-cells fig_ref_marching_cubes_cells
!fig ref-marching-cubes-resample
!fig ref-marching-cubes-cells
!< Figures of the reference page "Marching cubes". "What it does": a rounded box (marching cubes on 40^3 samples)
!< decimated to 1/16 of its facets, whose triangles are long slivers, re-tessellated by
!< resample_via_distance_field(resolution = 32): the facet sizes are now set by the grid. "Pipeline": the field of a
!< sphere of radius 2.3 about (3.4, -0.4, 3.4) sampled on 4^3 corners (3^3 cells of side 1): the corners inside
!< (value < 0, teal) and outside (ember), the triangles extract_isosurface emits, cell by cell, and their vertices
!< (gold), each interpolated on a cell edge between an inside and an outside corner.
program fig_ref_marching_cubes_cells
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 40
type(surface_stl_object)        :: box, resampled, cells
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P), allocatable   :: a(:), b(:), corners(:), feet(:)
logical,          allocatable   :: inside(:)
type(vector_R8P)                :: shift
real(R8P)                       :: phi(NG, NG, NG), v(4, 4, 4), p(3), q(3)
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, m, status

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   p = -1.1_R8P + 2.2_R8P / (NG - 1) * ([i, j, k] - 1)
   q = abs(p) - [0.8_R8P, 0.45_R8P, 0.3_R8P]
   phi(i, j, k) = norm2(max(q, 0._R8P)) + min(maxval(q), 0._R8P) - 0.15_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.1_R8P, -1.1_R8P, -1.1_R8P), &
                        bmax=vector_R8P(1.1_R8P, 1.1_R8P, 1.1_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call box%adopt_facets(facets=facets)
call box%sanitize
call box%decimate(target_facets=box%get_facets_number() / 16, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
call box%resample_via_distance_field(resolution=32_I4P, surface_out=resampled, status=status)
if (status /= MC_STATUS_OK) error stop 'resample failed'
write(caption, '(A,I0,A,F6.4,A,I0,A,L1,A,F6.4)') 'input ', box%get_facets_number(), ' facets, volume ', &
      box%get_volume(), '; resampled at 32: ', resampled%get_facets_number(), ', watertight ', &
      resampled%is_watertight(), ', volume ', resampled%get_volume()
print '(A)', trim(caption)
shift = vector_R8P(2.1_R8P, 0._R8P, 0._R8P)
call resampled%translate(delta=shift)
call fig_open(u, 'ref-marching-cubes-resample', trim(caption), -90._R8P, 35._R8P, 720_I4P, 300_I4P)
call fig_surface(u, 'body', box)
call fig_surface(u, 'body', resampled)
call fig_label(u, 'accent', vector_R8P(-0.9_R8P, -0.6_R8P, -0.6_R8P), 'decimated input')
call fig_label(u, 'accent', vector_R8P(-0.9_R8P, -0.6_R8P, -0.6_R8P) + shift, 'resample_via_distance_field(32)')
call fig_close(u)

allocate(corners(64), inside(64))
m = 0
do k = 1, 4 ; do j = 1, 4 ; do i = 1, 4
   v(i, j, k) = norm2([i - 1._R8P, j - 1._R8P, k - 1._R8P] - [3.4_R8P, -0.4_R8P, 3.4_R8P]) - 2.3_R8P
   m = m + 1 ; corners(m) = vector_R8P(i - 1._R8P, j - 1._R8P, k - 1._R8P) ; inside(m) = v(i, j, k) < 0._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=v, bmin=vector_R8P(0._R8P, 0._R8P, 0._R8P), bmax=vector_R8P(3._R8P, 3._R8P, 3._R8P), &
                        iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
allocate(feet(3 * size(facets)))
do i = 1, size(facets)
   feet(3 * i - 2 : 3 * i) = facets(i)%vertex
enddo
write(caption, '(A,I0,A,I0,A,I0,A)') '27 cells: ', count(inside), ' corners inside (teal), ', count(.not. inside), &
      ' outside (ember), ', size(facets), ' triangles'
print '(A)', trim(caption)
call cells%adopt_facets(facets=facets)
allocate(a(0), b(0))
do k = 0, 3 ; do j = 0, 3                                              ! the grid lines
   a = [a, vector_R8P(0._R8P, real(j, R8P), real(k, R8P)), vector_R8P(real(j, R8P), 0._R8P, real(k, R8P)), &
        vector_R8P(real(j, R8P), real(k, R8P), 0._R8P)]
   b = [b, vector_R8P(3._R8P, real(j, R8P), real(k, R8P)), vector_R8P(real(j, R8P), 3._R8P, real(k, R8P)), &
        vector_R8P(real(j, R8P), real(k, R8P), 3._R8P)]
enddo ; enddo
call fig_open(u, 'ref-marching-cubes-cells', trim(caption), -55._R8P, 28._R8P, 640_I4P, 480_I4P)
call fig_segments(u, 'box', a, b)
call fig_surface(u, 'body', cells)
call fig_points(u, 'inside', pack(corners, inside))
call fig_points(u, 'outside', pack(corners, .not. inside))
call fig_points(u, 'foot', feet)
call fig_close(u)
endprogram fig_ref_marching_cubes_cells
