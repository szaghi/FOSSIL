!run figref-curvature fig_ref_curvature
!fig ref-curvature-sphere-gaussian
!fig ref-curvature-sphere-mean
!fig ref-curvature-sphere-jet
!< Figures of the reference page "Per-vertex curvature": a unit sphere by marching cubes (the |x| - 1 field on 31^3
!< samples of [-1.5, 1.5]^3; the raw extraction has slits, so sanitize(fill_holes) repairs it), isotropically remeshed
!< (default target length and iterations), then gaussian_curvature, mean_curvature and, for comparison, the H of
!< jet_curvature (same colour scale as mean_curvature's). Each facet is coloured by the
!< mean of its three vertex values minus the exact 1 (teal below, ember above, saturating at +-SCALE); the caption
!< gives the median and the 2nd-98th percentile range of the per-vertex values.
program fig_ref_curvature
use fossil,              only : surface_stl_object, extract_isosurface, MC_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P),       parameter   :: M = 31
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:)
real(R8P),          allocatable :: gk(:), hm(:), hj(:), kj(:)
real(R8P)                       :: phi(M, M, M), h
integer(I4P)                    :: i, j, k, status

h = 3._R8P / real(M - 1, R8P)
do k = 1, M ; do j = 1, M ; do i = 1, M
   phi(i, j, k) = norm2(-1.5_R8P + h * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=facets)
call sphere%sanitize(fill_holes=.true.)                                      ! the raw extraction has slits
call sphere%isotropic_remesh(status=status)
if (.not. sphere%is_watertight()) error stop 'sphere not watertight'
call sphere%gaussian_curvature(K=gk, status=status)
call sphere%mean_curvature(H=hm, status=status)
call draw('ref-curvature-sphere-gaussian', 'K', gk, 1._R8P)
call sphere%jet_curvature(H=hj, K=kj, status=status)
call draw('ref-curvature-sphere-mean',     'H', hm, 0.5_R8P)
call draw('ref-curvature-sphere-jet',      'jet H', hj, 0.5_R8P)
contains
   subroutine draw(name, label, v, scale)
   !< The sphere coloured by (facet mean of v) - 1, saturating at +-scale.
   character(*), intent(in) :: name  !< Figure name.
   character(*), intent(in) :: label !< Field name.
   real(R8P),    intent(in) :: v(:)  !< Per-vertex values.
   real(R8P),    intent(in) :: scale !< Colour saturation.
   type(facet_object), pointer :: f
   real(R8P)                   :: s(size(v))
   character(96)               :: caption
   integer(I4P)                :: u, i, n

   s = v ; call sort(s) ; n = size(s, kind=I4P)
   write(caption, '(2A,F4.2,A,F5.2,A,F4.2,3A,F3.1)') label, ': median ', s(n / 2), ', 2-98 % ', s(max(1_I4P, n / 50)), &
      ' to ', s(n - n / 50), ' (exact 1); colour ', label, ' - 1 to +-', scale
   call fig_open(u, name, trim(caption), -60._R8P, 25._R8P)
   do i = 1, sphere%get_facets_number()
      f => sphere%facet_at(i)
      call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), &
                        value=(sum(v(f%vertex_id)) / 3._R8P - 1._R8P) / scale)
   enddo
   call fig_close(u)
   endsubroutine draw

   pure subroutine sort(a)
   !< Insertion sort (ascending).
   real(R8P), intent(inout) :: a(:) !< Values.
   real(R8P)                :: t
   integer(I4P)             :: i, j

   do i = 2, size(a, kind=I4P)
      t = a(i) ; j = i - 1
      do while (j >= 1)
         if (a(j) <= t) exit
         a(j + 1) = a(j) ; j = j - 1
      enddo
      a(j + 1) = t
   enddo
   endsubroutine sort
endprogram fig_ref_curvature
