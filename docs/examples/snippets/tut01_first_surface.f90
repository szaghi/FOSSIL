!< Tutorial 1: load an STL surface, read what analyze computed, save it.
program tut01_first_surface
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.) ! analyze runs inside
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
print '(A,I0)',       'facets      ', surface%get_facets_number()
print '(A,3F10.4)',   'box min     ', bmin%x, bmin%y, bmin%z
print '(A,3F10.4)',   'box max     ', bmax%x, bmax%y, bmax%z
print '(A,F12.4)',    'volume      ', surface%get_volume()
print '(A,L1,A,L1)',  'watertight  ', surface%is_watertight(), '   manifold ', surface%is_manifold()
call surface%save_into_file(file_name='dragon-ascii.stl', is_ascii=.true.)
print '(A)', 'saved dragon-ascii.stl'
endprogram tut01_first_surface
