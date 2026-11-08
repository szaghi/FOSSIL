use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface, back
type(facet_object), pointer :: f
integer(I4P), allocatable   :: ids(:)
real(R8P),    allocatable   :: xc(:)
integer(I4P)                :: i, p, n, u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
n = surface%get_facets_number()
allocate(xc(n))
do i = 1, n
   f => surface%facet_at(i) ; xc(i) = f%centroid%x
enddo
ids = [(i, i = 1, n)]
call surface%set_patch(facets=ids, name='middle')                        ! every facet, then two subsets
call surface%set_patch(facets=pack(ids, xc < -4._R8P), name='front')
call surface%set_patch(facets=pack(ids, xc >  4._R8P), name='rear')
call surface%save_into_file(file_name='dragon-patches.stl', is_ascii=.true.) ! one solid block per patch
