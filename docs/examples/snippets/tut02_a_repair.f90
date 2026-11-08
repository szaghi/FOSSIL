!< Tutorial 2, first run: repair the dragon with sanitize, check what the validity predicates say, save it.
program tut02_a_repair
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: facets_before

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
facets_before = surface%get_facets_number()
call surface%sanitize(fill_holes=.true.)
print '(A,I0,A,I0)', 'facets      ', facets_before, ' -> ', surface%get_facets_number()
print '(A,I0)',      'degenerate  ', surface%get_degenerate_facets_removed()
print '(A,I0)',      'duplicate   ', surface%get_duplicate_facets_removed()
print '(A,I0)',      'non-manif.  ', surface%get_non_manifold_edges_number()
print '(A,F12.4)',   'volume      ', surface%get_volume()
print '(A,L1,A,L1,A,L1)', 'watertight  ', surface%is_watertight(), '   manifold ', surface%is_manifold(), &
                          '   volume ', surface%is_volume()
call surface%save_into_file(file_name='dragon-clean.stl')
print '(A)', 'saved dragon-clean.stl'
endprogram tut02_a_repair
