!run ck_watertight ck_watertight
program ck_watertight
use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: box

! an open box: the facets of the top face (z = 1) are not loaded
call box%load_from_file(file_name='cube-01.stl', guess_format=.true.,          &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), &
                        clip_max=vector_R8P( 2._R8P,  2._R8P, 0.9_R8P))
print '(A,I0,A,L1)', 'open box:  ', box%get_facets_number(), ' facets, watertight ', box%is_watertight()
call box%sanitize(fill_holes=.true., max_hole_size=1._R8P)
print '(A,I0,A,L1,A,F8.5)', 'sanitized: ', box%get_facets_number(), ' facets, watertight ', box%is_watertight(), &
      ', volume ', box%get_volume()
endprogram ck_watertight
