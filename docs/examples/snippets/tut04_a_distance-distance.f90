use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: p(0:4), foot(0:4)
real(R8P)                :: d_pn, d_wn
logical                  :: inside(0:4)
integer(I4P)             :: k, u
character(16)            :: label

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
print '(A)', '    y      d (pseudo-n.)  d (winding n.) inside'
do k = 0, 4
   p(k) = 0._R8P * ex_R8P + (-0.6_R8P + 0.3_R8P * k) * ey_R8P + 0._R8P * ez_R8P
   d_pn = surface%distance(point=p(k), is_signed=.true., is_square_root=.true.)
   d_wn = surface%distance(point=p(k), is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   inside(k) = surface%is_point_inside(point=p(k), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   print '(F7.3,2F15.6,L5)', p(k)%y, d_pn, d_wn, inside(k)
enddo
