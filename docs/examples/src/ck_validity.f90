!run ck_validity ck_validity
program ck_validity
use fossil, only : surface_stl_object
implicit none
type(surface_stl_object) :: surface
character(len=18)        :: files(3) = [character(len=18) :: 'cube-01.stl', 'naca0012-ascii.stl', 'dragon.stl']
integer                  :: f

print '(A)', 'file                watertight manifold volume'
do f = 1, size(files)
   call surface%load_from_file(file_name=trim(files(f)), guess_format=.true.)
   print '(A,2X,L10,L9,L7)', files(f), surface%is_watertight(), surface%is_manifold(), surface%is_volume()
enddo
endprogram ck_validity
