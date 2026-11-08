!run figck-markers fig_ck_markers
!fig ck-markers
!< Figure of the cookbook recipe "Place Lagrangian markers for an IB method": the 288 markers of ck_markers on cube-01
!< (dx = 0.1, sharp edges preserved) on the three faces turned to the camera, each drawn with its unit normal scaled to
!< 0.08, with the cube's edges.
program fig_ck_markers
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_box, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)      :: cube
real(R8P), allocatable        :: centroids(:,:), areas(:), normals(:,:)
type(vector_R8P), allocatable :: a(:), b(:)
integer(I4P)                  :: i, u
character(100)                :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%surface_markers(dx=0.1_R8P, centroids=centroids, areas=areas, normals=normals, preserve_features=.true.)
allocate(a(0), b(0))
do i = 1, size(areas)                                              ! the markers on the faces turned to the camera
   if (dot_product(normals(:, i), [0.433_R8P, -0.75_R8P, 0.5_R8P]) <= 0._R8P) cycle
   a = [a, vector_R8P(centroids(1, i) + 0.01_R8P * normals(1, i), centroids(2, i) + 0.01_R8P * normals(2, i), &
                      centroids(3, i) + 0.01_R8P * normals(3, i))]  ! lifted off the face
   b = [b, vector_R8P(centroids(1, i) + 0.08_R8P * normals(1, i), centroids(2, i) + 0.08_R8P * normals(2, i), &
                      centroids(3, i) + 0.08_R8P * normals(3, i))]
enddo
write(caption, '(A,I0,A)') 'surface_markers(dx = 0.1): ', size(areas), ' markers (three faces shown)'
call fig_open(u, 'ck-markers', trim(caption), -60._R8P, 30._R8P)
call fig_box(u, 'box', cube%get_bmin(), cube%get_bmax())
call fig_segments(u, 'accent', a, b)
call fig_points(u, 'probe', a)
call fig_close(u)
endprogram fig_ck_markers
