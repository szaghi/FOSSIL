!run figck-hole-map fig_ck_hole_map
!fig ck-hole-map
!< Figure of the cookbook recipe "Cut holes with a hole map or X-rays": the background of ck_hole_map (97 x 25
!< vertices, h = 0.0125) around the NACA 0012 section, seen from above and zoomed on the nose. Hole vertices of the
!< hole map in teal (equal to the signed-distance holes), field vertices faint; the vertices the 33 X-rays get wrong
!< in ember, the section in yellow.
program fig_ck_hole_map
use fossil,     only : surface_stl_object, curve_object, slice_surface, curvilinear_block_2d_t, block_logical_t, &
                       hole_flags_2d, hole_map_2d, xray_flags_2d
use fossil_fig, only : fig_open, fig_zoom, fig_points, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)           :: wing
type(curve_object)                 :: cutter
type(curvilinear_block_2d_t)       :: background(1)
type(block_logical_t), allocatable :: exact(:), map(:), xray(:)
type(vector_R8P), allocatable      :: pts(:), a(:), b(:)
real(R8P)                          :: x(97, 25), y(97, 25), uv(2)
integer(I4P)                       :: i, j, n, un
character(100)                     :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, cutter)
do j = 1, 25 ; do i = 1, 97
   x(i, j) = -0.35_R8P + 0.0125_R8P * (i - 1) ; y(i, j) = -0.15_R8P + 0.0125_R8P * (j - 1)
enddo ; enddo
call background(1)%initialize(x, y)
call hole_flags_2d(background, cutter, exact)
call hole_map_2d(background, cutter, map, [16_I4P, 16_I4P])
call xray_flags_2d(background, cutter, xray, 33_I4P)
allocate(pts(97 * 25))
n = 0
do j = 1, 25 ; do i = 1, 97
   n = n + 1 ; pts(n) = vector_R8P(x(i, j), y(i, j), 0._R8P)
enddo ; enddo
n = cutter%get_vertices_number() ; allocate(a(n), b(n))
do i = 1, n
   uv = cutter%get_vertex(i) ; a(i) = vector_R8P(uv(1), uv(2), 1e-3_R8P)
   uv = cutter%get_vertex(cutter%get_next(i)) ; b(i) = vector_R8P(uv(1), uv(2), 1e-3_R8P)
enddo
write(caption, '(A,I0,A,I0,A)') 'hole map: ', count(map(1)%v), ' holes (teal); X-rays: ', &
                                count(xray(1)%v .neqv. exact(1)%v), ' vertices wrong (ember)'
call fig_open(un, 'ck-hole-map', trim(caption), -90._R8P, 90._R8P, height=360_I4P)
call fig_zoom(un, vector_R8P(-0.12_R8P, 0._R8P, 0._R8P), 0.13_R8P)
call fig_points(un, 'ghost', pack(pts, .not. map(1)%v))
call fig_points(un, 'inside', pack(pts, map(1)%v))
call fig_segments(un, 'accent', a, b)
call fig_points(un, 'hole', pack(pts, xray(1)%v .neqv. exact(1)%v))
call fig_close(un)
endprogram fig_ck_hole_map
