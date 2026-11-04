!run ck_save_formats ck_save_formats
program ck_save_formats
use fossil, only : surface_stl_object
use penf,   only : I8P
implicit none
type(surface_stl_object) :: surface
integer(I8P)             :: bytes

call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call surface%save_into_file(file_name='cube-ascii.stl', is_ascii=.true.)
call surface%save_into_file(file_name='cube-binary.stl', is_ascii=.false.)
inquire(file='cube-ascii.stl', size=bytes) ; print '(A,I0,A)', 'cube-ascii.stl  ', bytes, ' bytes'
inquire(file='cube-binary.stl', size=bytes) ; print '(A,I0,A)', 'cube-binary.stl ', bytes, ' bytes'
endprogram ck_save_formats
