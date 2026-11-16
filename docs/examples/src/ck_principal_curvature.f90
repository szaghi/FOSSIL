!run ck_principal_curvature ck_principal_curvature
program ck_principal_curvature
use fossil, only : surface_stl_object, facet_object, extract_isosurface, normal_curvature
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: body
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(61, 51, 41), d(1), k(2, 1), dirs(3, 2, 1), k2(2, 1), dirs2(3, 2, 1)
integer(I4P)                    :: i, j, l

! an ellipsoid with semi-axes 2, 1.5, 1 by marching cubes on [-2.4, 2.4] x [-2, 2] x [-1.6, 1.6], step 0.08
do l = 1, 41 ; do j = 1, 51 ; do i = 1, 61
   phi(i, j, l) = norm2([(-2.4_R8P + 0.08_R8P * (i - 1)) / 2, (-2._R8P + 0.08_R8P * (j - 1)) / 1.5_R8P, &
                         -1.6_R8P + 0.08_R8P * (l - 1)]) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-2.4_R8P, -2._R8P, -1.6_R8P), &
                        bmax=vector_R8P(2.4_R8P, 2._R8P, 1.6_R8P), iso=0._R8P, surface=facets)
call body%adopt_facets(facets=facets)
call body%build_curvature                       ! per-vertex principal curvatures and directions, once
call body%query_many(points=[vector_R8P(2.5_R8P, 0._R8P, 0._R8P)], distances=d, principal_curvatures=k, &
                     principal_directions=dirs)
print '(A,2F8.4,A)', 'k1, k2 at the tip  ', k(:, 1), '  (exact 2, 0.8889)'
print '(A,3F8.4)',   'd1 ', dirs(:, 1, 1)
print '(A,3F8.4)',   'd2 ', dirs(:, 2, 1)
print '(A,F8.4,A)',  'normal curvature along (0, 1, 1) ', normal_curvature(k(:, 1), dirs(:, :, 1), [0._R8P, 1._R8P, 1._R8P]), &
                     '  (exact 1.4444)'
call body%enter_device                          ! after build_curvature: the shape operator is mapped too
call body%query_many_device(points=[vector_R8P(2.5_R8P, 0._R8P, 0._R8P)], distances=d, principal_curvatures=k2, &
                            principal_directions=dirs2)
print '(A,L1)', 'device = host bit for bit: ', all(k2 == k) .and. all(dirs2 == dirs)
call body%exit_device
endprogram ck_principal_curvature
