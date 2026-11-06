!< Figures of the reference page "Thickness and gaps", under-resolution report: a plate 2 x 1 x 0.5 with a fin 0.05 thick
!< and 0.7 high (boolean union of two boxes), medial thickness, and thin_regions (factor 1) for a coarse cell h = 0.1 and
!< a fine one h = 0.04, below the fin thickness: the fin's sides drop out, its cap and end faces stay (a ball tangent off
!< their centre line meets the side walls, so the local feature size there is below 0.05). The facets of region 1 in
!< ember, of the other regions in gold, the rest plain.
program fig_ref_thin_regions
use fossil,              only : surface_stl_object, STATUS_OK, THICKNESS_MEDIAL
use fossil_boolean,      only : BOOL_UNION, BOOL_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)  :: part, fin
real(R8P), allocatable    :: th(:)
integer(I4P)              :: status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P])
call box(fin, [0.9_R8P, 0._R8P, 0.3_R8P], [0.95_R8P, 1._R8P, 1._R8P])
call part%boolean(fin, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%thickness(thickness=th, method=THICKNESS_MEDIAL, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call draw('ref-thin-regions-coarse', 0.1_R8P, '')
call draw('ref-thin-regions-fine', 0.04_R8P, '')
contains
   subroutine box(s, lo, hi)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)

   call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call s%resize(factor=(hi(1) - lo(1)) * ex_R8P + (hi(2) - lo(2)) * ey_R8P + (hi(3) - lo(3)) * ez_R8P)
   call s%translate(delta=lo(1) * ex_R8P + lo(2) * ey_R8P + lo(3) * ez_R8P)
   call s%analyze
   endsubroutine box

   subroutine draw(name, h, note)
   !< The part with the thin regions for the cell size h.
   character(*), intent(in)    :: name, note
   real(R8P),    intent(in)    :: h
   type(facet_object), pointer :: fp
   real(R8P),    allocatable   :: area(:)
   integer(I4P), allocatable   :: region(:)
   character(96)               :: caption
   integer(I4P)                :: u, i

   call part%thin_regions(thickness=th, h=h, factor=1._R8P, region=region, region_area=area)
   write(caption, '(A,F4.2,2A,I0,A,I0,A)') 'h = ', h, note, ': ', count(region > 0), ' facets thinner than h, ', &
                                           size(area), ' regions'
   print '(A)', trim(caption)
   call fig_open(u, name, trim(caption), -60._R8P, 25._R8P)
   do i = 1, part%get_facets_number()
      fp => part%facet_at(i)
      if (region(i) == 1) then
         call fig_triangle(u, 'hole', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      elseif (region(i) > 1) then
         call fig_triangle(u, 'accent', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      else
         call fig_triangle(u, 'body', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      endif
   enddo
   call fig_close(u)
   endsubroutine draw
endprogram fig_ref_thin_regions
