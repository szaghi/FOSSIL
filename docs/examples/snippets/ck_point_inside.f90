program ck_point_inside
use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                   SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: cube
type(vector_R8P)         :: p, q
integer(I4P)             :: a
integer(I4P), parameter  :: algorithms(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                             SIGN_FAST_WINDING_NUMBER]
character(len=14)        :: names(4) = [character(len=14) :: 'pseudo-normal', 'ray', 'solid angle', 'winding number']

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
p = vector_R8P(0.3_R8P, 0.6_R8P, 0.9_R8P) ! inside
q = vector_R8P(1.2_R8P, 0.5_R8P, 0.5_R8P) ! outside
print '(A)', 'inside?         p  q'
do a = 1, 4
   print '(A,2L3)', names(a), cube%is_point_inside(point=p, sign_algorithm=algorithms(a)), &
                              cube%is_point_inside(point=q, sign_algorithm=algorithms(a))
enddo
endprogram ck_point_inside
