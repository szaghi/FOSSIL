!< Figures of the cookbook recipe "Close the holes of a surface": the cube loaded without its top face, with its
!< boundary loop, and the same box after sanitize(fill_holes, max_hole_size=1), the patch in ember (the facets
!< sanitize appends after the 40 it kept). Seen from 50 degrees up: lower, the painter sort lets the walls cover a
!< corner of the patch.
program fig_ck_watertight
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: box
integer(I4P)             :: u, n0
character(80)            :: caption

call box%load_from_file(file_name='cube-01.stl', guess_format=.true.,          &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), &
                        clip_max=vector_R8P( 2._R8P,  2._R8P, 0.9_R8P))
n0 = box%get_facets_number()
write(caption, '(A,I0,A)') 'the open box: ', n0, ' facets, its boundary loop (is_watertight F)'
call fig_open(u, 'ck-watertight-open', trim(caption), -55._R8P, 50._R8P)
call fig_surface(u, 'body', box)
call fig_boundary(u, 'hole', box)
call fig_close(u)
call box%sanitize(fill_holes=.true., max_hole_size=1._R8P)
write(caption, '(A,I0,A)') 'sanitize(fill_holes): ', box%get_facets_number() - n0, ' facets close it (watertight T)'
call fig_open(u, 'ck-watertight-closed', trim(caption), -55._R8P, 50._R8P)
call fig_surface(u, 'body', box, last=n0)
call fig_surface(u, 'patch', box, first=n0 + 1_I4P)
call fig_close(u)
endprogram fig_ck_watertight
