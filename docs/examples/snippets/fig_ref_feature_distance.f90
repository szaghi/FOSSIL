!< Figure of the reference page "Feature edges", section "Distance to the nearest feature": the unit cube (its 12 edges
!< are the features at 30 degrees), build_features, then query_many(feature_distances=) on a 48 x 48 grid of points in
!< the plane z = 0.5 over [-0.5, 1.5]^2, through the cube. Each grid cell is coloured by the distance of its centre to
!< the nearest feature edge (the 12 edges drawn in ember), from 0 to 0.6 and above: the field grows from the four
!< vertical edges, where the plane crosses them, towards the plane's corners and the cube's axis.
program fig_ref_feature_distance
use fossil,     only : surface_stl_object, STATUS_OK
use fossil_fig, only : fig_open, fig_box, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter  :: N = 48
real(R8P),    parameter  :: Z = 0.5_R8P, H = 2._R8P / N
type(surface_stl_object) :: cube
type(vector_R8P)         :: points(N * N)
real(R8P)                :: d(N * N), fd(N * N), x, y
character(96)            :: caption
integer(I4P)             :: u, i, j, k, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%build_features(status=status)
if (status /= STATUS_OK) error stop 'build_features failed'
do j = 1, N ; do i = 1, N
   points(i + (j - 1) * N) = vector_R8P(-0.5_R8P + (i - 0.5_R8P) * H, -0.5_R8P + (j - 0.5_R8P) * H, Z)
enddo ; enddo
call cube%query_many(points=points, distances=d, feature_distances=fd, status=status)
if (status /= STATUS_OK) error stop 'query_many failed'
write(caption, '(A,F4.2,A,F5.3,A,F5.3)') 'distance to the nearest feature edge on z = ', Z, ': ', minval(fd), &
                                         ' to ', maxval(fd)
call fig_open(u, 'ref-feature-distance', trim(caption), -60._R8P, 40._R8P)
call fig_box(u, 'hole', vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(1._R8P, 1._R8P, 1._R8P))   ! the features
do k = 1, N * N
   x = points(k)%x ; y = points(k)%y
   call fig_quad(u, 'frac', vector_R8P(x - H / 2, y - H / 2, Z), vector_R8P(x + H / 2, y - H / 2, Z), &
                 vector_R8P(x + H / 2, y + H / 2, Z), vector_R8P(x - H / 2, y + H / 2, Z), value=min(1._R8P, fd(k) / 0.6_R8P))
enddo
call fig_close(u)
endprogram fig_ref_feature_distance
