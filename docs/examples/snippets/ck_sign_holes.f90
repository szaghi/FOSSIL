program ck_sign_holes
use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: box
type(vector_R8P)         :: p
logical                  :: truth
integer(I4P)             :: i, j, k, wrong_pn, wrong_fwn

! a cube without its top face
call box%load_from_file(file_name='cube-01.stl', guess_format=.true.,          &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), &
                        clip_max=vector_R8P( 2._R8P,  2._R8P, 0.9_R8P))
wrong_pn = 0 ; wrong_fwn = 0
do k = 1, 20 ; do j = 1, 20 ; do i = 1, 20 ! 8000 points in [-0.25, 1.25]^3
   p = vector_R8P(-0.25_R8P + 1.5_R8P * (i - 0.5_R8P) / 20, -0.25_R8P + 1.5_R8P * (j - 0.5_R8P) / 20, &
                  -0.25_R8P + 1.5_R8P * (k - 0.5_R8P) / 20)
   truth = all([p%x, p%y, p%z] > 0._R8P .and. [p%x, p%y, p%z] < 1._R8P)
   if (box%is_point_inside(point=p, sign_algorithm=SIGN_PSEUDO_NORMAL) .neqv. truth) wrong_pn = wrong_pn + 1
   if (box%is_point_inside(point=p, sign_algorithm=SIGN_FAST_WINDING_NUMBER) .neqv. truth) wrong_fwn = wrong_fwn + 1
enddo ; enddo ; enddo
print '(A,I0,A)', 'pseudo-normal:  ', wrong_pn, ' of 8000 points misclassified'
print '(A,I0,A)', 'winding number: ', wrong_fwn, ' of 8000 points misclassified'
endprogram ck_sign_holes
