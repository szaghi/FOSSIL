program ck_vertex_curvature
use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, CURV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: k_defect(:), h_cotan(:), h_jet(:), k_jet(:)
real(R8P)                       :: phi(31, 31, 31)
integer(I4P)                    :: i, j, k, status, s(3)

! the unit sphere by marching cubes, h = 0.1; 30 grid nodes lie on the sphere, (1, 0, 0), (0.6, 0, -0.8), ...
do k = 1, 31 ; do j = 1, 31 ; do i = 1, 31
   phi(i, j, k) = norm2(-1.5_R8P + 0.1_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
call sphere%gaussian_curvature(K=k_defect, status=s(1))    ! angle defect
call sphere%mean_curvature(H=h_cotan, status=s(2))         ! cotangent Laplacian
call sphere%jet_curvature(H=h_jet, K=k_jet, status=s(3))   ! cubic fit over the 2-ring
if (any(s /= CURV_STATUS_OK)) error stop 'curvature failed'
print '(A,I0,A)', 'vertices ', size(h_jet), ', exact H = K = 1'
print '(A,F8.4,A,F8.4)', 'angle defect K: min ', minval(k_defect), ', max ', maxval(k_defect)
print '(A,F8.4,A,F8.4)', 'cotangent H:    min ', minval(h_cotan), ', max ', maxval(h_cotan)
print '(A,F8.4,A,F8.4)', 'jet K:          min ', minval(k_jet), ', max ', maxval(k_jet)
print '(A,F8.4,A,F8.4)', 'jet H:          min ', minval(h_jet), ', max ', maxval(h_jet)
endprogram ck_vertex_curvature
