!run figck-query-point fig_ck_query_point
!fig ck-query-point
!< Figure of the cookbook recipe "Tell whether the closest point is on a face, an edge or a vertex": the unit cube,
!< the three query points of ck_query_point joined to their feet, each foot labelled with its feature and carrying
!< the pseudo-normal query_point returns (along the segment).
program fig_ck_query_point
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_points, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
character(6), parameter  :: KIND_NAME(3) = [character(6) :: 'face', 'edge', 'vertex']
type(surface_stl_object) :: cube
type(vector_R8P)         :: p(3), foot(3), normal
real(R8P)                :: d
integer(I4P)             :: i, feature, u

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
p = [vector_R8P(0.3_R8P, 0.15_R8P, 1.3_R8P), vector_R8P(1.3_R8P, 0.3_R8P, 1.3_R8P), &
     vector_R8P(1.3_R8P, 1.3_R8P, 1.3_R8P)]
call fig_open(u, 'ck-query-point', 'query_point: three points, their feet and the feature of each foot', &
              -35._R8P, 28._R8P)
call fig_surface(u, 'body', cube)
do i = 1, 3
   call cube%query_point(point=p(i), distance=d, closest_point=foot(i), feature=feature, normal=normal)
   call fig_segment(u, 'probe', p(i), foot(i))
   call fig_label(u, 'foot', foot(i) + vector_R8P(0.04_R8P, -0.06_R8P, -0.08_R8P), trim(KIND_NAME(feature)))
enddo
call fig_points(u, 'probe', p)
call fig_points(u, 'foot', foot)
call fig_close(u)
endprogram fig_ck_query_point
