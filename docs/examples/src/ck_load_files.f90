!run ck_load_files ck_load_files
program ck_load_files
use fossil, only : surface_stl_object
use penf,   only : I4P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: scene
type(vector_R8P)         :: bmin, bmax
integer(I4P)             :: p, i, n

call scene%load_from_files(file_names=[character(18) :: 'cube-01.stl', 'naca0012-ascii.stl']) ! binary + ASCII
print '(A,I0,A,I0,A)', 'one surface: ', scene%get_facets_number(), ' facets, ', scene%get_patches_number(), ' patches'
do p = 1, scene%get_patches_number()
   n = count([(scene%facet_patch(i) == p, i = 1, scene%get_facets_number())])
   print '(A,I0,A,I0,A)', 'patch ', p, ': ', n, ' facets'
enddo
bmin = scene%get_bmin() ; bmax = scene%get_bmax()
print '(A,3F8.3,A,3F8.3)', 'box ', bmin%x, bmin%y, bmin%z, ' to', bmax%x, bmax%y, bmax%z
endprogram ck_load_files
