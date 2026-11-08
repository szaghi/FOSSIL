!run figref-hole-filling fig_ref_hole_filling
!fig ref-hole-filling-plain
!fig ref-hole-filling-refined
!fig ref-hole-filling-faired
!< Figures of the reference page "Hole filling": a unit sphere (32 x 16 longitude-latitude facets) with the cap
!< z > 0.7 removed, filled three ways: the minimum-weight patch alone (flat: the loop is a planar ring), refined to the
!< density of the loop (refine), refined and faired to the thin plate (fair). The patch in ember. Printed per variant:
!< the patch facets, the holes filled, the largest distance of a patch centroid from the sphere, the volume, watertightness.
program fig_ref_hole_filling
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: NLON = 32, NLAT = 16
real(R8P),    parameter         :: PI = acos(-1._R8P)
type(surface_stl_object)        :: sphere, holed, filled
type(facet_object), allocatable :: facets(:), kept(:)
type(facet_object), pointer     :: f
real(R8P)                       :: th0, th1, ph0, ph1
integer(I4P)                    :: i, j, m, n

allocate(facets(2 * NLON * (NLAT - 1)))                                 ! no facet degenerates at the poles
n = 0
do j = 0, NLAT - 1
   th0 = PI * j / NLAT ; th1 = PI * (j + 1) / NLAT                      ! colatitude
   do i = 0, NLON - 1
      ph0 = 2 * PI * i / NLON ; ph1 = 2 * PI * mod(i + 1, NLON) / NLON   ! the seam closes exactly
      if (j > 0) call tri(at(th0, ph0), at(th1, ph0), at(th0, ph1))     ! outward: counter-clockwise seen from outside
      if (j < NLAT - 1) call tri(at(th0, ph1), at(th1, ph0), at(th1, ph1))
   enddo
enddo
call sphere%adopt_facets(facets=facets)
print '(A,I0,A,L1,A,F8.5)', 'sphere: ', sphere%get_facets_number(), ' facets, watertight ', sphere%is_watertight(), &
                            ', volume ', sphere%get_volume()
allocate(kept(n))
m = 0
do i = 1, sphere%get_facets_number()
   f => sphere%facet_at(i)
   if (f%centroid%z > 0.7_R8P) cycle                                     ! the cap
   m = m + 1 ; kept(m) = f
enddo
allocate(facets, source=kept(1:m))                                    ! adopt_facets took the array
call holed%adopt_facets(facets=facets)
call variant('ref-hole-filling-plain', 'fill_holes: minimum-weight patch (ember), flat on a planar loop', .false., .false.)
call variant('ref-hole-filling-refined', 'fill_holes(refine=.true.): patch refined to the density of its loop', &
             .true., .false.)
call variant('ref-hole-filling-faired', 'fill_holes(fair=.true.): refined, then faired to the thin plate', .true., .true.)
contains
   function at(th, ph) result(p)
   !< Point of the unit sphere at colatitude th, longitude ph.
   real(R8P), intent(in) :: th, ph
   real(R8P)             :: p(3)

   p = [sin(th) * cos(ph), sin(th) * sin(ph), cos(th)]
   if (th >= PI) p = [0._R8P, 0._R8P, -1._R8P]                           ! the south pole, exactly
   endfunction at

   subroutine tri(a, b, c)
   !< Append the facet a, b, c.
   real(R8P), intent(in) :: a(3), b(3), c(3)

   n = n + 1
   facets(n)%vertex(1) = vector_R8P(a(1), a(2), a(3))
   facets(n)%vertex(2) = vector_R8P(b(1), b(2), b(3))
   facets(n)%vertex(3) = vector_R8P(c(1), c(2), c(3))
   call facets(n)%compute_metrix
   endsubroutine tri

   subroutine variant(name, caption, refine, fair)
   !< Fill the holed sphere one way, print its measures, draw it.
   character(*), intent(in) :: name, caption
   logical,      intent(in) :: refine, fair
   real(R8P)                :: dmax
   integer(I4P)             :: u, k, holes

   filled = holed
   call filled%fill_holes(max_hole_size=1._R8P, holes_filled=holes, refine=refine, fair=fair)
   dmax = 0._R8P
   do k = m + 1, filled%get_facets_number()
      f => filled%facet_at(k)
      dmax = max(dmax, abs(1._R8P - f%centroid%normL2()))
   enddo
   print '(A,L1,A,L1,A,I0,A,I0,A,F6.3,A,F8.5,A,L1)', 'refine ', refine, ', fair ', fair, ': holes ', holes, &
         ', patch facets ', filled%get_facets_number() - m, ', max distance from the sphere ', dmax, ', volume ', &
         filled%get_volume(), ', watertight ', filled%is_watertight()
   call fig_open(u, name, caption, -60._R8P, 18._R8P)
   call fig_surface(u, 'body', filled, first=1_I4P, last=m)
   call fig_surface(u, 'patch', filled, first=m + 1_I4P)
   call fig_close(u)
   endsubroutine variant
endprogram fig_ref_hole_filling
