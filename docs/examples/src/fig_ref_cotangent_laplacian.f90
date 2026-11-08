!run figref-cotangent-laplacian fig_ref_cotangent_laplacian
!fig ref-laplacian-sphere-Lx
!< Figure of the reference page "Cotangent Laplacian": the remeshed unit sphere of the curvature page (marching cubes
!< on 31^3 samples of [-1.5, 1.5]^3, no slivers since issue #87, isotropic_remesh), cotangent_laplacian(L, M), and the
!< action of the operator on the x-coordinate field: (M^-1 L x)_i = (L x)_i / M_ii, which on the unit sphere tends to
!< the Laplace-Beltrami of x, -2 x. Each facet is coloured by the mean of its three vertex values, divided by 2
!< (teal -1 at +x, ember +1 at -x); the caption gives the RMS and maximum deviation from -2 x over the vertices.
program fig_ref_cotangent_laplacian
use fossil,              only : surface_stl_object, csr_matrix_t, extract_isosurface, MC_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P),       parameter   :: M = 31
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(csr_matrix_t)              :: lap, mass
real(R8P),          allocatable :: x(:), lx(:)
real(R8P)                       :: phi(M, M, M), h
character(96)                   :: caption
integer(I4P)                    :: u, i, j, k, n, status

h = 3._R8P / real(M - 1, R8P)
do k = 1, M ; do j = 1, M ; do i = 1, M
   phi(i, j, k) = norm2(-1.5_R8P + h * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
call sphere%isotropic_remesh(status=status)
if (.not. sphere%is_watertight()) error stop 'sphere not watertight'
call sphere%cotangent_laplacian(L=lap, M=mass, status=status)
n = lap%get_nrows()
allocate(x(n), lx(n))
do i = 1, sphere%get_facets_number()                                         ! x of every pool vertex
   f => sphere%facet_at(i)
   do k = 1, 3
      x(f%vertex_id(k)) = f%vertex(k)%x
   enddo
enddo
call lap%multiply_vector(x=x, y=lx)
do i = 1, n
   lx(i) = lx(i) / mass%get_value(i, i)                                     ! M^-1 L x
enddo
write(caption, '(A,F5.3,A,F4.2,A)') 'M^-1 L x (exact -2x): RMS error ', sqrt(sum((lx + 2._R8P * x)**2) / n), &
                                    ', max ', maxval(abs(lx + 2._R8P * x)), '; colour value / 2'
call fig_open(u, 'ref-laplacian-sphere-Lx', trim(caption), -90._R8P, 15._R8P)    ! +x to the right
do i = 1, sphere%get_facets_number()
   f => sphere%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), value=sum(lx(f%vertex_id)) / 6._R8P)
enddo
call fig_close(u)
endprogram fig_ref_cotangent_laplacian
