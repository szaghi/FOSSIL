!run figck-closest_point fig_ck_closest_point
!fig ck-closest_point
!< Figure of the cookbook recipe "Find the closest point and its normal": the cube, the query point of the recipe,
!< the closest point on the edge between the faces x = 1 and z = 1 with the distance segment (teal), the facet
!< compute_distance returns (ember) and the normal there, the edge pseudo-normal (0.5 long).
program fig_ck_closest_point
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: cube
type(vector_R8P)         :: point, foot, normal
real(R8P)                :: d
integer(I4P)             :: u, facet
character(80)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
point = vector_R8P(1.5_R8P, 0.5_R8P, 2._R8P)
call cube%compute_distance(point=point, distance=d, is_signed=.true., is_square_root=.true., facet_index=facet, &
                           closest_point=foot, normal=normal)
write(caption, '(A,F7.5,A,I0,A)') 'compute_distance: d = ', d, ', facet ', facet, ', the normal at the closest point'
call fig_open(u, 'ck-closest_point', trim(caption), -60._R8P, 15._R8P)
call fig_surface(u, 'body', cube)
call fig_surface(u, 'patch', cube, first=facet, last=facet)
call fig_segment(u, 'probe', point, foot)
call fig_segment(u, 'accent', foot, foot + 0.5_R8P * normal)
call fig_points(u, 'probe', [point])
call fig_points(u, 'foot', [foot])
call fig_label(u, 'probe', point, 'point')
call fig_label(u, 'accent', foot + 0.5_R8P * normal, 'normal')
call fig_close(u)
endprogram fig_ck_closest_point
