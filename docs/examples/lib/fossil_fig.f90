!< Scene writer of the documentation figures (issue #47): the records of scripts/mesh2svg.py.
module fossil_fig
!< Scene writer of the documentation figures (issue #47): the records of scripts/mesh2svg.py.
!<
!< A figure program computes its geometry with FOSSIL and writes it here; `!fig NAME` in the program makes
!< scripts/docs_examples.sh draw NAME.fig into docs/examples/images/NAME.svg. Styles (roles of the palette): body,
!< ghost (the body before an operation, translucent), patch (what an operation added), hole (boundary edges), probe,
!< foot, inside, outside, box, accent, cell; field (a value in [-1, 1]) and frac (a value in [0, 1]) colour a mesh.
!< Coordinates are written with 8 significant digits: the drawing, not the data, is the product.
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
private
public :: fig_open, fig_zoom, fig_surface, fig_boundary, fig_box, fig_segment, fig_segments, fig_points, fig_label
public :: fig_triangle, fig_quad, fig_close

character(*), parameter :: F3 = '(3ES16.8)', F6 = '(6ES16.8)', F9 = '(9ES16.8)', F10 = '(10ES16.8)'

contains
   subroutine fig_open(u, name, caption, azimuth, elevation, width, height)
   !< Open NAME.fig with its canvas, camera (degrees) and caption.
   integer(I4P), intent(out)          :: u         !< Unit of the scene.
   character(*), intent(in)           :: name      !< Figure name (the file is NAME.fig).
   character(*), intent(in)           :: caption   !< Caption under the drawing.
   real(R8P),    intent(in)           :: azimuth   !< View azimuth about +z from +x.
   real(R8P),    intent(in)           :: elevation !< View elevation above the xy plane (90: from above).
   integer(I4P), intent(in), optional :: width     !< Canvas width in px (default 640).
   integer(I4P), intent(in), optional :: height    !< Canvas height in px (default 420).
   integer(I4P)                       :: w, h

   w = 640_I4P ; if (present(width)) w = width
   h = 420_I4P ; if (present(height)) h = height
   open(newunit=u, file=name//'.fig', status='replace', action='write')
   write(u, '(A,I0,A,I0)') 'size ', w, ' ', h
   write(u, '(A,F8.2,F8.2)') 'camera ', azimuth, elevation
   write(u, '(A)') 'caption '//caption
   endsubroutine fig_open

   subroutine fig_zoom(u, centre, radius)
   !< Frame the sphere about `centre` instead of the whole scene.
   integer(I4P),     intent(in) :: u      !< Unit of the scene.
   type(vector_R8P), intent(in) :: centre !< Centre of the view.
   real(R8P),        intent(in) :: radius !< Radius shown.

   write(u, '(A,4ES16.8)') 'zoom ', centre%x, centre%y, centre%z, radius
   endsubroutine fig_zoom

   subroutine fig_surface(u, style, surface, first, last)
   !< The facets first..last (default all) of a surface.
   integer(I4P),             intent(in)           :: u       !< Unit of the scene.
   character(*),             intent(in)           :: style   !< Style.
   type(surface_stl_object), intent(in)           :: surface !< Surface.
   integer(I4P),             intent(in), optional :: first   !< First facet.
   integer(I4P),             intent(in), optional :: last    !< Last facet.
   type(facet_object), pointer                    :: f
   integer(I4P)                                   :: i, i0, i1

   i0 = 1_I4P ; if (present(first)) i0 = first
   i1 = surface%get_facets_number() ; if (present(last)) i1 = last
   write(u, '(A,A,A,I0)') 'mesh ', style, ' ', max(0_I4P, i1 - i0 + 1_I4P)
   do i = i0, i1
      f => surface%facet_at(i)
      write(u, F9) f%vertex(1)%x, f%vertex(1)%y, f%vertex(1)%z, f%vertex(2)%x, f%vertex(2)%y, f%vertex(2)%z, &
                   f%vertex(3)%x, f%vertex(3)%y, f%vertex(3)%z
   enddo
   endsubroutine fig_surface

   subroutine fig_boundary(u, style, surface)
   !< The boundary edges of a surface: facet edges without a neighbour (`fcon_edge == 0`).
   integer(I4P),             intent(in) :: u       !< Unit of the scene.
   character(*),             intent(in) :: style   !< Style.
   type(surface_stl_object), intent(in) :: surface !< Surface.
   type(facet_object), pointer          :: f
   integer(I4P)                         :: i, e, n

   n = 0
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i) ; n = n + count(f%fcon_edge == 0_I4P)
   enddo
   write(u, '(A,A,A,I0)') 'segs ', style, ' ', n
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i)
      do e = 1, 3
         if (f%fcon_edge(e) /= 0_I4P) cycle
         write(u, F6) f%vertex(e)%x, f%vertex(e)%y, f%vertex(e)%z, &
                      f%vertex(mod(e, 3) + 1)%x, f%vertex(mod(e, 3) + 1)%y, f%vertex(mod(e, 3) + 1)%z
      enddo
   enddo
   endsubroutine fig_boundary

   subroutine fig_box(u, style, bmin, bmax)
   !< The 12 edges of an axis-aligned box.
   integer(I4P),     intent(in) :: u          !< Unit of the scene.
   character(*),     intent(in) :: style      !< Style.
   type(vector_R8P), intent(in) :: bmin, bmax !< Box corners.
   real(R8P)                    :: b(2, 3)
   integer(I4P)                 :: j, k

   b(1, :) = [bmin%x, bmin%y, bmin%z] ; b(2, :) = [bmax%x, bmax%y, bmax%z]
   write(u, '(A,A,A)') 'segs ', style, ' 12'
   do k = 1, 2
      do j = 1, 2
         write(u, F6) b(1, 1), b(j, 2), b(k, 3), b(2, 1), b(j, 2), b(k, 3)
         write(u, F6) b(j, 1), b(1, 2), b(k, 3), b(j, 1), b(2, 2), b(k, 3)
         write(u, F6) b(j, 1), b(k, 2), b(1, 3), b(j, 1), b(k, 2), b(2, 3)
      enddo
   enddo
   endsubroutine fig_box

   subroutine fig_segment(u, style, a, b)
   !< One segment.
   integer(I4P),     intent(in) :: u     !< Unit of the scene.
   character(*),     intent(in) :: style !< Style.
   type(vector_R8P), intent(in) :: a, b  !< Ends.

   write(u, '(A,A,A)') 'segs ', style, ' 1'
   write(u, F6) a%x, a%y, a%z, b%x, b%y, b%z
   endsubroutine fig_segment

   subroutine fig_segments(u, style, a, b)
   !< Segments a(i)-b(i).
   integer(I4P),     intent(in) :: u          !< Unit of the scene.
   character(*),     intent(in) :: style      !< Style.
   type(vector_R8P), intent(in) :: a(:), b(:) !< Ends.
   integer(I4P)                 :: i

   write(u, '(A,A,A,I0)') 'segs ', style, ' ', size(a)
   do i = 1, size(a)
      write(u, F6) a(i)%x, a(i)%y, a(i)%z, b(i)%x, b(i)%y, b(i)%z
   enddo
   endsubroutine fig_segments

   subroutine fig_points(u, style, p)
   !< Points.
   integer(I4P),     intent(in) :: u     !< Unit of the scene.
   character(*),     intent(in) :: style !< Style.
   type(vector_R8P), intent(in) :: p(:)  !< Points.
   integer(I4P)                 :: i

   write(u, '(A,A,A,I0)') 'pts ', style, ' ', size(p)
   do i = 1, size(p)
      write(u, F3) p(i)%x, p(i)%y, p(i)%z
   enddo
   endsubroutine fig_points

   subroutine fig_label(u, style, p, text)
   !< A label anchored at a point.
   integer(I4P),     intent(in) :: u     !< Unit of the scene.
   character(*),     intent(in) :: style !< Style (its colour).
   type(vector_R8P), intent(in) :: p     !< Anchor.
   character(*),     intent(in) :: text  !< Text.

   write(u, '(A,A,A,3ES16.8,A,A)') 'label ', style, ' ', p%x, p%y, p%z, ' ', text
   endsubroutine fig_label

   subroutine fig_triangle(u, style, a, b, c, value)
   !< One triangle, with an optional value (styles field, frac).
   integer(I4P),     intent(in)           :: u       !< Unit of the scene.
   character(*),     intent(in)           :: style   !< Style.
   type(vector_R8P), intent(in)           :: a, b, c !< Vertices.
   real(R8P),        intent(in), optional :: value   !< Value.

   write(u, '(A,A,A)') 'mesh ', style, ' 1'
   if (present(value)) then
      write(u, F10) a%x, a%y, a%z, b%x, b%y, b%z, c%x, c%y, c%z, value
   else
      write(u, F9) a%x, a%y, a%z, b%x, b%y, b%z, c%x, c%y, c%z
   endif
   endsubroutine fig_triangle

   subroutine fig_quad(u, style, a, b, c, d, value)
   !< A planar quadrilateral a-b-c-d (two triangles, no seam), with an optional value (styles field, frac).
   integer(I4P),     intent(in)           :: u          !< Unit of the scene.
   character(*),     intent(in)           :: style      !< Style.
   type(vector_R8P), intent(in)           :: a, b, c, d !< Vertices, in order around the quad.
   real(R8P),        intent(in), optional :: value      !< Value.

   write(u, '(A,A,A)') 'mesh ', style, ' 2'
   if (present(value)) then
      write(u, F10) a%x, a%y, a%z, b%x, b%y, b%z, c%x, c%y, c%z, value
      write(u, F10) a%x, a%y, a%z, c%x, c%y, c%z, d%x, d%y, d%z, value
   else
      write(u, F9) a%x, a%y, a%z, b%x, b%y, b%z, c%x, c%y, c%z
      write(u, F9) a%x, a%y, a%z, c%x, c%y, c%z, d%x, d%y, d%z
   endif
   endsubroutine fig_quad

   subroutine fig_close(u)
   !< Close the scene.
   integer(I4P), intent(in) :: u !< Unit of the scene.

   close(u)
   endsubroutine fig_close
endmodule fossil_fig
