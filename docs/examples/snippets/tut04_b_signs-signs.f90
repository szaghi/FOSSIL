use fossil,     only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_RAY_INTERSECTIONS, &
                       SIGN_SOLID_ANGLE
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter  :: ALG(4) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_RAY_INTERSECTIONS, &
                                      SIGN_SOLID_ANGLE]
type(surface_stl_object) :: surface
type(vector_R8P)         :: p(3), foot, normal
real(R8P)                :: d
integer(I4P)             :: k, a, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
p(1) = 0._R8P * ex_R8P
p(2) = 0.6_R8P * ey_R8P
p(3) = 0.36_R8P * ex_R8P - 0.2_R8P * ey_R8P + 0.2_R8P * ez_R8P
print '(A)', '        point          pseudo-n.  winding n.  rays  solid angle'
do k = 1, 3
   print '(3F7.3,4L10)', p(k)%x, p(k)%y, p(k)%z, &
                         (surface%is_point_inside(point=p(k), sign_algorithm=ALG(a)), a = 1, 4)
enddo
