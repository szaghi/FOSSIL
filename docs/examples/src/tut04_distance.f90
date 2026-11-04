!run tut04 tut04_distance
!< Tutorial 4: signed distance of a few points, inside test, closest point and normal.
program tut04_distance
!region distance
use fossil, only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p, closest, normal
real(R8P)                :: d_pn, d_wn, d
integer(I4P)             :: k, facet, edge, vertex

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
print '(A)', '    y      d (pseudo-n.)  d (winding n.) inside'
do k = 0, 4
   p = 0._R8P * ex_R8P + (-0.6_R8P + 0.3_R8P * k) * ey_R8P + 0._R8P * ez_R8P
   d_pn = surface%distance(point=p, is_signed=.true., is_square_root=.true.)
   d_wn = surface%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   print '(F7.3,2F15.6,L5)', p%y, d_pn, d_wn, surface%is_point_inside(point=p, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo
!endregion distance
!region rich
p = 0._R8P * ex_R8P + 0.6_R8P * ey_R8P + 0._R8P * ez_R8P
call surface%compute_distance(point=p, distance=d, is_signed=.true., is_square_root=.true., &
                              facet_index=facet, edge_index=edge, vertex_index=vertex, closest_point=closest, normal=normal)
print '(A,F10.6,A,I0,A,I0,A,I0)', 'd ', d, '  facet ', facet, '  edge ', edge, '  vertex ', vertex
print '(A,3F10.6)', 'closest point ', closest%x, closest%y, closest%z
print '(A,3F10.6)', 'normal        ', normal%x, normal%y, normal%z
!endregion rich
endprogram tut04_distance
