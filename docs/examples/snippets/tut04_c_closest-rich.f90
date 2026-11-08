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
