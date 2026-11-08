!run figck-thin-walls fig_ck_thin_walls
!fig ck-thin-walls
!< Figure of the cookbook recipe "Find walls too thin for my grid": the NACA 0012 wing, the facets thinner than the
!< cell size h = chord / 32 (medial thickness, thin_regions with factor 1) in ember: they gather at the trailing edge.
program fig_ck_thin_walls
use fossil,     only : surface_stl_object, facet_object, THICKNESS_MEDIAL
use fossil_fig, only : fig_open, fig_triangle, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)    :: wing
type(facet_object), pointer :: f
type(vector_R8P)            :: bmin, bmax
real(R8P),    allocatable   :: thickness(:), gap(:), area(:)
integer(I4P), allocatable   :: region(:)
real(R8P)                   :: h
integer(I4P)                :: u, i
character(96)               :: caption

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL)
bmin = wing%get_bmin() ; bmax = wing%get_bmax()
h = (bmax%x - bmin%x) / 32._R8P
call wing%thin_regions(thickness=thickness, gap=gap, h=h, factor=1._R8P, region=region, region_area=area)
write(caption, '(A,ES9.2,A,I0,A,I0,A)') 'h = chord / 32 = ', h, ': ', count(region > 0), ' facets thinner (ember), ', &
                                        size(area), ' regions'
call fig_open(u, 'ck-thin-walls', trim(caption), -60._R8P, 35._R8P)
do i = 1, wing%get_facets_number()
   f => wing%facet_at(i)
   if (region(i) > 0) then
      call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
   else
      call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
   endif
enddo
call fig_close(u)
endprogram fig_ck_thin_walls
