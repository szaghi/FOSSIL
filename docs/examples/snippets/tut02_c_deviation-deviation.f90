use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: repaired, original
type(facet_object), pointer :: f
real(R8P), allocatable      :: d(:)
real(R8P)                   :: lower, upper
integer(I4P)                :: i, u
character(96)               :: caption

call repaired%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
call repaired%deviation(original, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'repaired -> original ', lower, ' <= h <= ', upper
call original%deviation(repaired, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'original -> repaired ', lower, ' <= h <= ', upper
