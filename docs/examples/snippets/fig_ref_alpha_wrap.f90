!< Figures of the reference page "Alpha wrapping": the unit cube minus its first facet (the page's broken input; the
!< hole's boundary in ember, the camera looking at it); then holed%alpha_wrap(alpha=0.15, offset=0.033), cut open by
!< the plane y = the hole's centroid (the wrap facets on the near side are not drawn): the hole is wider than the ball,
!< so the wrap is one closed shell lining both sides of every wall, two offset layers kept apart by the barrier;
!< then holed%alpha_wrap(alpha=0.6), cut the same way: the ball does not fit through the hole, the wrap caps it.
program fig_ref_alpha_wrap
use fossil,              only : surface_stl_object, AWRAP_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_boundary, fig_triangle, fig_close
use penf,                only : I4P, R8P

implicit none
real(R8P), parameter                    :: DEG = 57.29577951308232_R8P
type(surface_stl_object)                :: source, holed, wrapped
type(facet_object), pointer             :: facets(:)
type(facet_object), allocatable, target :: holed_facets(:)
character(96)                           :: caption
real(R8P)                               :: azimuth, elevation, cut
integer(I4P)                            :: u, i

call source%load_from_file(file_name='cube-01.stl', guess_format=.true.)
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1))                                 ! drop the first facet, as the page does
do i = 1, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1)
enddo
azimuth = atan2(facets(1)%normal%y, facets(1)%normal%x) * DEG - 30._R8P     ! look at the hole, a little askew
elevation = asin(max(-1._R8P, min(1._R8P, facets(1)%normal%z))) * DEG * 0.5_R8P + 20._R8P
cut = facets(1)%centroid%y                                                 ! cutaway: the wrap in y <= the hole's y
call holed%adopt_facets(facets=holed_facets)
write(caption, '(A,I0,A,L1)') 'cube-01.stl minus facet 1: ', holed%get_facets_number(), ' facets, watertight ', &
                              holed%is_watertight()
call fig_open(u, 'ref-alpha-wrap-input', trim(caption), azimuth, elevation)
call fig_surface(u, 'body', holed)
call fig_boundary(u, 'hole', holed)
call fig_close(u)
call cutaway('ref-alpha-wrap-output', 0.15_R8P, 0.033_R8P)
call cutaway('ref-alpha-wrap-capped', 0.6_R8P, 0.02_R8P)
contains
   subroutine cutaway(name, alpha, offset)
   !< The wrap of the holed cube, cut open by the plane y = cut (the facets beyond it are not drawn).
   character(*), intent(in)    :: name          !< Figure name.
   real(R8P),    intent(in)    :: alpha, offset !< Wrap parameters.
   type(facet_object), pointer :: f
   integer(I4P)                :: status, j

   call holed%alpha_wrap(alpha=alpha, offset=offset, wrapped=wrapped, status=status)
   if (status /= AWRAP_STATUS_OK) error stop 'alpha_wrap failed'
   write(caption, '(A,F4.2,A,F5.3,A,I0,A,L1,A,F6.4)') 'alpha_wrap(', alpha, ', ', offset, '), cut open: ', &
         wrapped%get_facets_number(), ' facets, watertight ', wrapped%is_watertight(), ', volume ', wrapped%get_volume()
   call fig_open(u, name, trim(caption), azimuth, elevation)
   do j = 1, wrapped%get_facets_number()
      f => wrapped%facet_at(j)
      if (f%centroid%y <= cut) call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
   enddo
   call fig_close(u)
   endsubroutine cutaway
endprogram fig_ref_alpha_wrap
