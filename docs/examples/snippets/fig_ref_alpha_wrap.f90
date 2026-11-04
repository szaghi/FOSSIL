!< Figures of the reference page "Alpha wrapping": the unit cube minus its first facet (the page's broken input; the
!< hole's boundary in ember, the camera looking at it), then holed%alpha_wrap(alpha=0.1, offset=0.033) from the same
!< camera, cut open by the plane y = the hole's centroid (the wrap facets on the near side are not drawn): the hole is
!< wider than alpha, so the flood fill gets inside and the wrap is a closed hollow shell lining both sides of the input
!< (volume ~0.35, not the cube's 1), entering through the hole.
program fig_ref_alpha_wrap
use fossil,              only : surface_stl_object, AWRAP_STATUS_OK, AWRAP_STATUS_NOT_CONVERGED
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_boundary, fig_triangle, fig_close
use penf,                only : I4P, R8P

implicit none
real(R8P), parameter                    :: DEG = 57.29577951308232_R8P
type(surface_stl_object)                :: source, holed, wrapped
type(facet_object), pointer             :: facets(:), f
type(facet_object), allocatable, target :: holed_facets(:)
character(96)                           :: caption
real(R8P)                               :: azimuth, elevation, cut
integer(I4P)                            :: u, i, status

call source%load_from_file(file_name='cube-01.stl', guess_format=.true.)
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1))                                 ! drop the first facet, as the page does
do i = 1, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1)
enddo
azimuth = atan2(facets(1)%normal%y, facets(1)%normal%x) * DEG - 30._R8P     ! look at the hole, a little askew
elevation = asin(max(-1._R8P, min(1._R8P, facets(1)%normal%z))) * DEG * 0.5_R8P + 20._R8P
call holed%adopt_facets(facets=holed_facets)
write(caption, '(A,I0,A,L1)') 'cube-01.stl minus facet 1: ', holed%get_facets_number(), ' facets, watertight ', &
                              holed%is_watertight()
call fig_open(u, 'ref-alpha-wrap-input', trim(caption), azimuth, elevation)
call fig_surface(u, 'body', holed)
call fig_boundary(u, 'hole', holed)
call fig_close(u)
call holed%alpha_wrap(alpha=0.1_R8P, offset=0.033_R8P, wrapped=wrapped, status=status)
if (status /= AWRAP_STATUS_OK .and. status /= AWRAP_STATUS_NOT_CONVERGED) error stop 'alpha_wrap failed'
write(caption, '(A,I0,A,L1,A,F6.4)') 'alpha_wrap(0.1, 0.033), cut open: ', wrapped%get_facets_number(), &
                                     ' facets, watertight ', wrapped%is_watertight(), ', volume ', wrapped%get_volume()
call fig_open(u, 'ref-alpha-wrap-output', trim(caption), azimuth, elevation)
cut = facets(1)%centroid%y                                                 ! cutaway: the wrap in y <= the hole's y
do i = 1, wrapped%get_facets_number()
   f => wrapped%facet_at(i)
   if (f%centroid%y <= cut) call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_close(u)
endprogram fig_ref_alpha_wrap
