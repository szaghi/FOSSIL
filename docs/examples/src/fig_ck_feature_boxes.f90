!run figck-feature-boxes fig_ck_feature_boxes
!fig ck-feature-boxes
!< Figure of the cookbook recipe "Refine the cells along sharp edges": the unit cube (translucent) and the cells of
!< side 0.1 that box_feature_many finds touching one of its 12 sharp edges, outlined.
program fig_ck_feature_boxes
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter  :: N = 13
type(surface_stl_object) :: cube
real(R8P)                :: lo(3, N**3), hi(3, N**3)
logical                  :: touched(N**3)
integer(I4P)             :: i, j, k, c, u
character(80)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%build_features(angle=30._R8P)
c = 0
do k = 1, N ; do j = 1, N ; do i = 1, N
   c = c + 1 ; lo(:, c) = -0.13_R8P + 0.1_R8P * ([i, j, k] - 1) ; hi(:, c) = lo(:, c) + 0.1_R8P
enddo ; enddo ; enddo
call cube%box_feature_many(lo=lo, hi=hi, touched=touched)
write(caption, '(A,I0,A,I0,A)') 'box_feature_many: ', count(touched), ' of ', N**3, ' cells touch a sharp edge'
call fig_open(u, 'ck-feature-boxes', trim(caption), -55._R8P, 25._R8P)
call fig_surface(u, 'ghost', cube)
do c = 1, N**3
   if (touched(c)) call fig_box(u, 'accent', vector_R8P(lo(1, c), lo(2, c), lo(3, c)), &
                                vector_R8P(hi(1, c), hi(2, c), hi(3, c)))
enddo
call fig_close(u)
endprogram fig_ck_feature_boxes
