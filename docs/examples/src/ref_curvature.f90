!run ref_curvature ref_curvature
program ref_curvature
use fossil, only : surface_stl_object, facet_object, csr_matrix_t, extract_isosurface, MC_STATUS_OK, CURV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P), RT = 0.45_R8P   ! tube radius of the torus
type(surface_stl_object) :: sphere, torus
type(csr_matrix_t)       :: L, M
real(R8P), allocatable   :: K(:), H(:), Hjet(:), Kjet(:), area(:), k1(:), k2(:), exact(:)
integer(I4P)             :: status
!region sphere
call implicit_body(sphere, torus=.false.)                ! unit sphere by marching cubes, then sanitize
call sphere%gaussian_curvature(K=K, status=status)       ! angle defect / barycentric area
if (status /= CURV_STATUS_OK) error stop 'gaussian_curvature failed'
call sphere%mean_curvature(H=H)                          ! |M^-1 L x| / 2, signed by the pseudo-normal
call sphere%jet_curvature(H=Hjet, K=Kjet)                ! cubic fit over the 2-ring
print '(A,I0,A,L1)', 'sphere: vertices ', size(K), ', watertight ', sphere%is_watertight()
print '(A)', '                   mean     2%    98%   (percentiles over the vertices)'
print '(A,3F7.3)', '  K angle defect', sum(K) / size(K), pct(K, 2), pct(K, 98)
print '(A,3F7.3)', '  H cotangent   ', sum(H) / size(H), pct(H, 2), pct(H, 98)
print '(A,3F7.3)', '  H jet         ', sum(Hjet) / size(Hjet), pct(Hjet, 2), pct(Hjet, 98)
!endregion sphere
!region bonnet
call vertex_area(sphere, area)                           ! M_ii of cotangent_laplacian: the area each vertex owns
print '(A,ES8.1)', 'sphere: sum K_i A_i / (2 pi) - 2 ', sum(K * area) / (2 * PI) - 2._R8P ! Euler characteristic 2
call implicit_body(torus, torus=.true.)                  ! a torus, R = 1, r = 0.45: genus 1
call torus%gaussian_curvature(K=K)
call vertex_area(torus, area)
print '(A,ES8.1)', 'torus:  sum K_i A_i / (2 pi)     ', sum(K * area) / (2 * PI)          ! Euler characteristic 0
!endregion bonnet
!region principal
call torus%jet_curvature(H=Hjet, K=Kjet)
k1 = Hjet + sqrt(max(0._R8P, Hjet**2 - Kjet))            ! principal curvatures from H and K, k1 >= k2
k2 = Hjet - sqrt(max(0._R8P, Hjet**2 - Kjet))
call torus_exact(torus, exact)                           ! (rho - R) / (r rho) at every vertex, rho = |(x, y)|
k1 = k1 - 1._R8P / RT ; k2 = k2 - exact                 ! errors
print '(A,2F7.4)', 'torus:  k1 - 1/r,   rms and max ', sqrt(sum(k1**2) / size(k1)), maxval(abs(k1))
print '(A,2F7.4)', '        k2 - exact, rms and max ', sqrt(sum(k2**2) / size(k2)), maxval(abs(k2))
!endregion principal
contains
   subroutine implicit_body(s, torus)
   !< Marching cubes of a unit sphere or of a torus (R = 1, r = 0.45, axis z) on 41^3 samples, then sanitize.
   type(surface_stl_object), intent(inout) :: s
   logical,                  intent(in)    :: torus
   type(facet_object), allocatable         :: facets(:)
   real(R8P)                               :: phi(41, 41, 41), p(3)
   integer(I4P)                            :: i, j, k
   do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41
      p = -1.5_R8P + 0.075_R8P * ([i, j, k] - 1)
      if (torus) then
         phi(i, j, k) = sqrt((norm2(p(1:2)) - 1._R8P)**2 + p(3)**2) - RT
      else
         phi(i, j, k) = norm2(p) - 1._R8P
      endif
   enddo ; enddo ; enddo
   call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                           bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
   if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
   call s%adopt_facets(facets=facets)
   call s%sanitize(fill_holes=.true.)                    ! drops the slivers marching cubes leaves
   endsubroutine implicit_body

   subroutine vertex_area(s, a)
   !< Barycentric area of every pool vertex, the diagonal of the mass matrix.
   type(surface_stl_object), intent(in) :: s
   real(R8P), allocatable, intent(out)  :: a(:)
   integer(I4P)                         :: n
   call s%cotangent_laplacian(L=L, M=M)
   allocate(a(M%get_nrows()))
   do n = 1, size(a)
      a(n) = M%get_value(n, n)
   enddo
   endsubroutine vertex_area

   subroutine torus_exact(s, kx)
   !< The smaller principal curvature of the torus R = 1, r = 0.45 at every pool vertex: (rho - R) / (r rho).
   type(surface_stl_object), intent(in) :: s
   real(R8P), allocatable, intent(out)  :: kx(:)
   type(facet_object), pointer          :: f
   real(R8P)                            :: rho
   integer(I4P)                         :: n, v
   allocate(kx(size(Kjet)))
   do n = 1, s%get_facets_number()
      f => s%facet_at(n)
      do v = 1, 3
         rho = sqrt(f%vertex(v)%x**2 + f%vertex(v)%y**2)
         kx(f%vertex_id(v)) = (rho - 1._R8P) / (RT * rho)
      enddo
   enddo
   endsubroutine torus_exact

   function pct(v, q) result(x)
   !< The q-th percentile of v (by sorting a copy).
   real(R8P),    intent(in) :: v(:)
   integer(I4P), intent(in) :: q
   real(R8P)                :: x, w(size(v)), t
   integer(I4P)             :: a, b
   w = v
   do a = 2, size(w)                                     ! insertion sort
      t = w(a) ; b = a - 1
      do while (b >= 1)
         if (w(b) <= t) exit
         w(b + 1) = w(b) ; b = b - 1
      enddo
      w(b + 1) = t
   enddo
   x = w(max(1, nint(q * size(w) / 100._R8P)))
   endfunction pct
endprogram ref_curvature
