!run ck_clip ck_clip
program ck_clip
use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: dragon, tail

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
print '(A,I0)', 'dragon facets ', dragon%get_facets_number()
! keep the facets inside the box (a facet is kept if it touches it), the others go to the remainder
call dragon%clip(bmin=vector_R8P(-12._R8P, -6._R8P, 4._R8P), bmax=vector_R8P(0._R8P, 6._R8P, 21._R8P), &
                 remainder=tail)
print '(A,I0,A,I0)', 'clipped ', dragon%get_facets_number(), ', remainder ', tail%get_facets_number()
print '(A,L1,A,L1)', 'watertight: clipped ', dragon%is_watertight(), ', remainder ', tail%is_watertight()
endprogram ck_clip
