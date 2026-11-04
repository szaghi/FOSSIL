program ck_flip
use fossil, only : surface_stl_object
use penf,   only : I4P
implicit none
type(surface_stl_object) :: cube
integer(I4P)             :: flipped

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%reverse_normals ; call cube%analyze ! an inside-out cube
print '(A,F9.5,A,L1)', 'inside-out: volume ', cube%get_volume(), ', is_volume ', cube%is_volume()
call cube%sanitize_normals(flipped_count=flipped) ; call cube%analyze
print '(A,F9.5,A,L1,A,I0,A)', 'fixed:      volume ', cube%get_volume(), ', is_volume ', cube%is_volume(), &
      ' (', flipped, ' facets flipped)'
endprogram ck_flip
