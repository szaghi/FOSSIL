!< Figures of the reference page "Mesh smoothing": the marching-cubes unit sphere (radius 1 on a 31^3 grid, 30 nodes on
!< the sphere: no slivers since issue #87, so no repair) with deterministic Gaussian noise added to every vertex (sigma
!< 1 % of the bounding-box diagonal, the same displacement for every copy of a vertex: it is a function of the vertex
!< coordinates), then the same surface after smooth at its defaults (5 Taubin pairs, lambda 0.5, mu -0.53). The
!< captions give the volume and the rms of |v| - 1 over the vertices.
program fig_ref_smoothing
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, SMOOTH_METHOD_TAUBIN, &
                       SMOOTH_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 31
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
real(R8P)                       :: phi(NG, NG, NG), sigma
type(vector_R8P)                :: d
integer(I4P)                    :: i, j, k, status

do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   phi(i, j, k) = norm2(-1.5_R8P + 3._R8P / (NG - 1) * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
d = sphere%get_bmax() - sphere%get_bmin()
sigma = 0.01_R8P * d%normL2()
allocate(facets(sphere%get_facets_number()))                                       ! the noisy copy
do i = 1, size(facets)
   f => sphere%facet_at(i)
   facets(i) = f
   do j = 1, 3
      facets(i)%vertex(j) = facets(i)%vertex(j) + sigma * noise(facets(i)%vertex(j))
   enddo
   call facets(i)%compute_normal
enddo
call sphere%adopt_facets(facets=facets)
call scene('ref-smoothing-sphere-noisy', 'noise sigma 1 % of the box diagonal')
call sphere%smooth(method=SMOOTH_METHOD_TAUBIN, status=status)
if (status /= SMOOTH_STATUS_OK) error stop 'smooth failed'
call scene('ref-smoothing-sphere-taubin', 'smooth, 5 Taubin pairs')
contains
   function noise(p) result(g)
   !< Three standard Gaussian deviates (Box-Muller), a deterministic function of the point `p`.
   type(vector_R8P), intent(in) :: p !< Point.
   type(vector_R8P)             :: g !< Deviates.
   real(R8P)                    :: r(6)
   integer(I4P)                 :: s

   do s = 1, 6
      r(s) = hash(p%x * (12.9898_R8P + s) + p%y * (78.233_R8P + 3 * s) + p%z * (37.719_R8P + 7 * s))
   enddo
   r(1:3) = max(r(1:3), 1.e-12_R8P)
   g%x = sqrt(-2._R8P * log(r(1))) * cos(2._R8P * acos(-1._R8P) * r(4))
   g%y = sqrt(-2._R8P * log(r(2))) * cos(2._R8P * acos(-1._R8P) * r(5))
   g%z = sqrt(-2._R8P * log(r(3))) * cos(2._R8P * acos(-1._R8P) * r(6))
   endfunction noise

   pure function hash(x) result(h)
   !< A uniform deviate in [0, 1) from a real (the classic fract(sin(x) * 43758.5453)).
   real(R8P), intent(in) :: x !< Seed.
   real(R8P)             :: h !< Deviate.

   h = sin(x) * 43758.5453_R8P
   h = h - floor(h)
   endfunction hash

   subroutine scene(name, title)
   !< The sphere, its volume and its rms radial error in the caption.
   character(*), intent(in) :: name, title !< Figure name, caption head.
   real(R8P)                :: e
   character(100)           :: caption
   integer(I4P)             :: u, i, j

   e = 0._R8P
   do i = 1, sphere%get_facets_number()
      f => sphere%facet_at(i)
      do j = 1, 3
         e = e + (f%vertex(j)%normL2() - 1._R8P)**2
      enddo
   enddo
   e = sqrt(e / (3 * sphere%get_facets_number()))
   write(caption, '(A,A,F5.3,A,F6.4)') title, ': volume ', sphere%get_volume(), ', rms(|v| - 1) ', e
   call fig_open(u, name, trim(caption), -60._R8P, 25._R8P)
   call fig_surface(u, 'body', sphere)
   call fig_close(u)
   endsubroutine scene
endprogram fig_ref_smoothing
