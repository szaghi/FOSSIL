!run figref-smoothing-steps fig_ref_smoothing_steps
!fig ref-smoothing-stability
!fig ref-smoothing-steps
!< Figures of the reference page "Mesh smoothing" on the page's example surface, the marching-cubes unit sphere on 41^3
!< samples (sanitized with fill_holes). "Why uniform weights": the explicit step bound of the cotangent operator,
!< M_ii / L_ii per vertex (cotangent_laplacian: barycentric mass over the diagonal), each facet coloured by its
!< smallest vertex value on a log scale (sepia = the smallest): the global explicit step is set by the few vertices of
!< tiny mass, while the uniform Laplacian is stable for any lambda <= 1 everywhere. "Pipeline": the same sphere with
!< the page's noise (sigma 0.035), zoomed: from each vertex, its lambda step (ember: smooth(EXPLICIT, lambda = 0.5,
!< iterations = 1)) and then the mu step back out (teal: the first Taubin pair, smooth(TAUBIN, iterations = 1)).
program fig_ref_smoothing_steps
use fossil,     only : surface_stl_object, facet_object, csr_matrix_t, extract_isosurface, MC_STATUS_OK, &
                       SMOOTH_METHOD_EXPLICIT, SMOOTH_METHOD_TAUBIN, SMOOTH_STATUS_OK, LAPL_STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_triangle, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)        :: sphere, noisy, lam, pair
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f, g, h
type(csr_matrix_t)              :: L, M
type(vector_R8P), allocatable   :: a(:), b(:), c(:)
type(vector_R8P)                :: centre
real(R8P),        allocatable   :: tau(:)
real(R8P)                       :: phi(41, 41, 41), lo, hi, med
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status

do k = 1, 41 ; do j = 1, 41 ; do i = 1, 41
   phi(i, j, k) = norm2(-1.5_R8P + 0.075_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
call sphere%sanitize(fill_holes=.true.)

call sphere%cotangent_laplacian(L=L, M=M, status=status)
if (status /= LAPL_STATUS_OK) error stop 'cotangent_laplacian failed'
allocate(tau(M%n_rows))
do i = 1, M%n_rows                                                   ! M_ii / |L_ii|
   tau(i) = M%val(M%row_ptr(i))
   do j = L%row_ptr(i), L%row_ptr(i + 1) - 1
      if (L%col_idx(j) == i) tau(i) = tau(i) / abs(L%val(j))
   enddo
enddo
lo = log10(minval(tau)) ; hi = log10(maxval(tau)) ; med = median(tau)
write(caption, '(A,ES7.1,A,ES7.1,A,ES7.1,A)') 'step bound M_ii/L_ii: min ', minval(tau), ', median ', med, &
      ', max ', maxval(tau), ' (sepia: smallest)'
print '(A)', trim(caption)
call fig_open(u, 'ref-smoothing-stability', trim(caption), -60._R8P, 25._R8P)
do i = 1, sphere%get_facets_number()
   f => sphere%facet_at(i)
   call fig_triangle(u, 'frac', f%vertex(1), f%vertex(2), f%vertex(3), &
                     1._R8P - (log10(minval(tau(f%vertex_id))) - lo) / (hi - lo))
enddo
call fig_close(u)

allocate(facets(sphere%get_facets_number()))                         ! the page's noise
do i = 1, size(facets)
   f => sphere%facet_at(i)
   facets(i) = f
   do j = 1, 3
      facets(i)%vertex(j) = facets(i)%vertex(j) + 0.035_R8P * noise(facets(i)%vertex(j))
   enddo
   call facets(i)%compute_normal
enddo
call noisy%adopt_facets(facets=facets)
lam = noisy
call lam%smooth(method=SMOOTH_METHOD_EXPLICIT, lambda=0.5_R8P, iterations=1_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'explicit smoothing failed'
pair = noisy
call pair%smooth(method=SMOOTH_METHOD_TAUBIN, iterations=1_I4P, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'Taubin smoothing failed'
centre = vector_R8P(0.5_R8P, -0.6_R8P, 0.62_R8P)
allocate(a(0), b(0), c(0))
do i = 1, noisy%get_facets_number()
   f => noisy%facet_at(i) ; g => lam%facet_at(i) ; h => pair%facet_at(i)
   do j = 1, 3
      if (norm2([f%vertex(j)%x - centre%x, f%vertex(j)%y - centre%y, f%vertex(j)%z - centre%z]) > 0.3_R8P) cycle
      a = [a, f%vertex(j)] ; b = [b, g%vertex(j)] ; c = [c, h%vertex(j)]
   enddo
enddo
write(caption, '(A)') 'one Taubin pair: the lambda = 0.5 step (ember), then the mu = -0.53 step (teal)'
call fig_open(u, 'ref-smoothing-steps', trim(caption), -60._R8P, 25._R8P)
call fig_zoom(u, centre, 0.3_R8P)
do i = 1, noisy%get_facets_number()                                  ! the facets about the centre only
   f => noisy%facet_at(i)
   if (norm2([f%centroid%x - centre%x, f%centroid%y - centre%y, f%centroid%z - centre%z]) < 0.35_R8P) &
      call fig_triangle(u, 'ghost', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_segments(u, 'outside', a, b)
call fig_segments(u, 'probe', b, c)
call fig_points(u, 'accent', a)
call fig_close(u)
print '(A,I0,A)', 'steps drawn from ', size(a), ' vertex copies'
contains
   function median(x) result(m)
   !< Lower median of x (selection by counting, small arrays).
   real(R8P), intent(in) :: x(:) !< Values.
   real(R8P)             :: m    !< Median.
   integer(I4P)          :: p

   do p = 1, size(x)
      if (count(x < x(p)) <= (size(x) - 1) / 2 .and. count(x <= x(p)) > (size(x) - 1) / 2) then
         m = x(p) ; return
      endif
   enddo
   m = x(1)
   endfunction median

   function noise(p) result(gd)
   !< Three standard Gaussian deviates (Box-Muller), a deterministic function of the point `p` (as ref_smoothing).
   type(vector_R8P), intent(in) :: p  !< Point.
   type(vector_R8P)             :: gd !< Deviates.
   real(R8P)                    :: r(6)
   integer(I4P)                 :: s

   do s = 1, 6
      r(s) = hash(p%x * (12.9898_R8P + s) + p%y * (78.233_R8P + 3 * s) + p%z * (37.719_R8P + 7 * s))
   enddo
   r(1:3) = max(r(1:3), 1.e-12_R8P)
   gd%x = sqrt(-2._R8P * log(r(1))) * cos(2._R8P * acos(-1._R8P) * r(4))
   gd%y = sqrt(-2._R8P * log(r(2))) * cos(2._R8P * acos(-1._R8P) * r(5))
   gd%z = sqrt(-2._R8P * log(r(3))) * cos(2._R8P * acos(-1._R8P) * r(6))
   endfunction noise

   pure function hash(x) result(hh)
   !< A uniform deviate in [0, 1) from a real (fract(sin(x) * 43758.5453)).
   real(R8P), intent(in) :: x  !< Seed.
   real(R8P)             :: hh !< Deviate.

   hh = sin(x) * 43758.5453_R8P
   hh = hh - floor(hh)
   endfunction hash
endprogram fig_ref_smoothing_steps
