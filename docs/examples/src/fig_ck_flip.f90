!run figck-flip fig_ck_flip
!fig ck-flip-before
!fig ck-flip-after
!< Figures of the cookbook recipe "Fix an inside-out surface": the facet normals of the inside-out cube (they point
!< in, ember) and after sanitize_normals (they point out), as segments from the facet centroids with a dot at the
!< tip; the cube is a ghost so the normals inside show.
program fig_ck_flip
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_box, fig_segments, fig_points, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
type(surface_stl_object) :: cube
integer(I4P)             :: u, flipped
character(80)            :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%reverse_normals ; call cube%analyze
write(caption, '(A,F6.2,A)') 'inside-out: the normals point in, volume ', cube%get_volume()
call scene('ck-flip-before', trim(caption), 'outside')
call cube%sanitize_normals(flipped_count=flipped) ; call cube%analyze
write(caption, '(A,I0,A,F5.2)') 'sanitize_normals: ', flipped, ' facets flipped, volume ', cube%get_volume()
call scene('ck-flip-after', trim(caption), 'accent')
contains
   subroutine scene(name, caption, style)
   !< The cube (ghost, with its bounding box as outline) and its facet normals, 0.2 long.
   character(*), intent(in)    :: name, caption !< Figure name and caption.
   character(*), intent(in)    :: style         !< Style of the normals.
   type(facet_object), pointer :: f
   type(vector_R8P)            :: a(cube%get_facets_number()), b(cube%get_facets_number())
   integer(I4P)                :: i

   do i = 1, cube%get_facets_number()
      f => cube%facet_at(i)
      a(i) = f%centroid ; b(i) = f%centroid + 0.2_R8P * f%normal
   enddo
   call fig_open(u, name, caption, -55._R8P, 30._R8P)
   call fig_surface(u, 'ghost', cube)
   call fig_box(u, 'box', cube%get_bmin(), cube%get_bmax())
   call fig_segments(u, style, a, b)
   call fig_points(u, style, b)
   call fig_close(u)
   endsubroutine scene
endprogram fig_ck_flip
