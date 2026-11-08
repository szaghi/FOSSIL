!run tut04-closest tut04_c_closest
!fig tut04-closest
!< Tutorial 4, third run: the closest facet, feature, point and normal of one point, from compute_distance. The figure
!< zooms on the closest point: the direction to the point, the vertex pseudo-normal and the normal of
!< facet 1289.
program tut04_c_closest
use fossil_fig,          only : fig_open, fig_zoom, fig_surface, fig_segment, fig_points, fig_close
!region rich
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)    :: surface
type(facet_object), pointer :: f
type(vector_R8P)            :: p, closest, normal
real(R8P)                   :: d
integer(I4P)                :: facet, edge, vertex, i, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
p = 0._R8P * ex_R8P + 0.6_R8P * ey_R8P + 0._R8P * ez_R8P
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=facet, edge_index=edge, vertex_index=vertex, closest_point=closest, normal=normal)
print '(A,F10.6,A,I0,A,I0,A,I0)', 'd ', d, '  facet ', facet, '  edge ', edge, '  vertex ', vertex
print '(A,3F10.6)', 'closest point ', closest%x, closest%y, closest%z
print '(A,3F10.6)', 'normal        ', normal%x, normal%y, normal%z
f => surface%facet_at(facet)
print '(A,3F10.6)', 'facet normal  ', f%normal%x, f%normal%y, f%normal%z
!endregion rich
call fig_open(u, 'tut04-closest', 'closest vertex: direction to the point (gold), pseudo-normal (teal), facet normal (grey)', &
              45._R8P, -25._R8P)
call fig_zoom(u, closest, 0.06_R8P)
call fig_surface(u, 'body', surface)
call fig_segment(u, 'accent', closest, closest + 0.05_R8P * (p - closest) / d)
call fig_segment(u, 'inside', closest, closest + 0.05_R8P * normal)
call fig_segment(u, 'box', f%centroid, f%centroid + 0.03_R8P * f%normal)
call fig_points(u, 'foot', [closest])
call fig_close(u)
endprogram tut04_c_closest
