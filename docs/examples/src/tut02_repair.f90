!run tut02 tut02_repair
!< Tutorial 2: repair the dragon with sanitize, check it is a closed volume.
program tut02_repair
!region load
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P

implicit none
type(surface_stl_object) :: surface, original
integer(I4P)             :: facets_before
real(R8P)                :: lower, upper

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
facets_before = surface%get_facets_number()
!endregion load
!region sanitize
call surface%sanitize(fill_holes=.true.)
!endregion sanitize
!region check
print '(A,I0,A,I0)', 'facets      ', facets_before, ' -> ', surface%get_facets_number()
print '(A,I0)',      'degenerate  ', surface%get_degenerate_facets_removed()
print '(A,I0)',      'duplicate   ', surface%get_duplicate_facets_removed()
print '(A,I0)',      'non-manif.  ', surface%get_non_manifold_edges_number()
print '(A,F12.4)',   'volume      ', surface%get_volume()
print '(A,L1,A,L1,A,L1)', 'watertight  ', surface%is_watertight(), '   manifold ', surface%is_manifold(), &
                          '   volume ', surface%is_volume()
call surface%save_into_file(file_name='dragon-clean.stl')
print '(A)', 'saved dragon-clean.stl'
!endregion check
!region deviation
call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
call surface%deviation(original, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'repaired -> original ', lower, ' <= h <= ', upper
call original%deviation(surface, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'original -> repaired ', lower, ' <= h <= ', upper
!endregion deviation
endprogram tut02_repair
