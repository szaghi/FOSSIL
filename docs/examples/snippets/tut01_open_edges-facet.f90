use fossil,                    only : surface_stl_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_facet_object,       only : facet_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
type(surface_stl_object)          :: surface
type(facet_object), pointer       :: f
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:)
integer(I4P), allocatable         :: edges(:,:), kinds(:)
integer(I4P)                      :: i, v, open_edges, u
character(96)                     :: caption

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
f => surface%facet_at(1_I4P)                                            ! a pointer into the surface, read-only use
do v = 1, 3
   print '(A,I0,A,3F9.4)', 'facet 1 vertex ', v, '   ', f%vertex(v)%x, f%vertex(v)%y, f%vertex(v)%z
enddo
print '(A,3F9.4)', 'facet 1 normal     ', f%normal%x, f%normal%y, f%normal%z
print '(A,3I6)',   'facet 1 neighbours ', f%fcon_edge                    ! across edges 1-2, 2-3, 3-1; 0: none
