!< Figure of the cookbook recipe "Get a reliable sign on a mesh with holes": the cube without its top face (ghost,
!< its boundary loop in ember) and the points of the recipe's 20^3 grid that SIGN_PSEUDO_NORMAL misclassifies,
!< coloured by the side it puts them on (teal inside, ember outside). SIGN_FAST_WINDING_NUMBER misclassifies none.
program fig_ck_sign_holes
use fossil,     only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_box, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: box
type(vector_R8P)         :: p(8000)
logical                  :: truth(8000), pn(8000), fwn(8000)
integer(I4P)             :: u, i, j, k, n
character(80)            :: caption

call box%load_from_file(file_name='cube-01.stl', guess_format=.true.,          &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), &
                        clip_max=vector_R8P( 2._R8P,  2._R8P, 0.9_R8P))
n = 0
do k = 1, 20 ; do j = 1, 20 ; do i = 1, 20
   n = n + 1
   p(n) = vector_R8P(-0.25_R8P + 1.5_R8P * (i - 0.5_R8P) / 20, -0.25_R8P + 1.5_R8P * (j - 0.5_R8P) / 20, &
                     -0.25_R8P + 1.5_R8P * (k - 0.5_R8P) / 20)
   truth(n) = all([p(n)%x, p(n)%y, p(n)%z] > 0._R8P .and. [p(n)%x, p(n)%y, p(n)%z] < 1._R8P)
   pn(n)  = box%is_point_inside(point=p(n), sign_algorithm=SIGN_PSEUDO_NORMAL)
   fwn(n) = box%is_point_inside(point=p(n), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo ; enddo ; enddo
write(caption, '(A,I0,A,I0)') 'misclassified of 8000: pseudo-normal ', count(pn .neqv. truth), &
                              ', winding number ', count(fwn .neqv. truth)
call fig_open(u, 'ck-sign_holes', trim(caption), -60._R8P, 30._R8P)
call fig_surface(u, 'ghost', box)
call fig_boundary(u, 'hole', box)
call fig_points(u, 'inside', pack(p, pn .and. .not. truth))
call fig_points(u, 'outside', pack(p, .not. pn .and. truth))
call fig_close(u)
endprogram fig_ck_sign_holes
