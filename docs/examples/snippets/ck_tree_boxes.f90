!< Cookbook: read the surface's own acceleration tree (issue #89): the leaves of the BVH, their facets, and the leaf
!< boxes that touch a region (a block of an AMR grid, a debugging probe).
program ck_tree_boxes
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)      :: dragon
type(vector_R8P), allocatable :: bmin(:), bmax(:)
integer(I4P),     allocatable :: depth(:), start(:), ids(:)
logical,          allocatable :: leaf(:)
type(vector_R8P)              :: lo, hi
integer(I4P)                  :: i, touching, facets

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%tree_boxes(bmin=bmin, bmax=bmax, depth=depth, leaf=leaf, facet_start=start, facet_ids=ids)
print '(A,I0,A,I0,A,I0)', 'boxes ', size(bmin), ', leaves ', count(leaf), ', deepest ', maxval(depth)
print '(A,I0,A,I0,A,I0)', 'facets in the leaves ', size(ids), ' of ', dragon%get_facets_number(), &
                          ', largest leaf ', maxval(start(2:) - start(:size(start) - 1))
lo = vector_R8P(-2._R8P, -1._R8P, 8._R8P) ; hi = vector_R8P(2._R8P, 1._R8P, 12._R8P)   ! a region of the domain
touching = 0 ; facets = 0
do i = 1, size(bmin)
   if (.not. leaf(i)) cycle
   if (any([bmax(i)%x < lo%x, bmax(i)%y < lo%y, bmax(i)%z < lo%z, bmin(i)%x > hi%x, bmin(i)%y > hi%y, bmin(i)%z > hi%z])) &
      cycle
   touching = touching + 1 ; facets = facets + start(i + 1) - start(i)
enddo
print '(A,I0,A,I0,A)', 'leaf boxes touching the region ', touching, ' (', facets, ' facets to test)'
endprogram ck_tree_boxes
