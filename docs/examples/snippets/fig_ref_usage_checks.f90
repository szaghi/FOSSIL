!< Figures of the usage reference: "Validity predicates", the dragon as loaded, then after sanitize, with its boundary
!< edges (red, one facet) and non-manifold edges (gold, 3+ facets) from feature_edges and the three predicates;
!< "Connect nearby vertices", the
!< cube of the page whose first facet is shifted by 1e-4, zoomed on one corner of that facet, before and after
!< connect_nearby_vertices (two frames: the crack, then the vertices merged at their centroid); "Full workflow
!< example", the dragon after the page's workflow: sanitized, its centroid moved to the origin, its longest side scaled
!< to 1.
program fig_ref_usage_checks
use fossil,                    only : surface_stl_object, facet_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_fig,                only : fig_open, fig_zoom, fig_surface, fig_boundary, fig_segments, fig_box, fig_points, &
                                      fig_label, fig_anim, fig_frame, fig_close
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object)          :: dragon, cube
type(facet_object), allocatable   :: facets(:)
type(facet_object), pointer       :: f
type(vertex_pool_object), pointer :: pool
type(vector_R8P), allocatable     :: a(:), b(:), p(:)
type(vector_R8P)                  :: corner, bmin, bmax, c
integer(I4P), allocatable         :: edges(:,:), kinds(:)
character(120)                    :: caption
integer(I4P)                      :: u, i, j

! validity: what the predicates see on the dragon as loaded, then after sanitize
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call fig_open(u, 'anim-ref-usage-validity', '', -60._R8P, 22._R8P)
call fig_anim(u, 3._R8P)
call edges_frame('as loaded')
call dragon%sanitize
call edges_frame('after sanitize')
call fig_close(u)

! connect_nearby_vertices: the crack of the page's cube, at one corner of the shifted facet
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number()))
do i = 1, size(facets)
   f => cube%facet_at(i) ; facets(i) = f
enddo
corner = facets(1)%vertex(2)                                      ! (0, 0.5, 0): on the cube's edge x = z = 0
do i = 1, 3
   facets(1)%vertex(i) = facets(1)%vertex(i) + 1.e-4_R8P * ex_R8P
enddo
call cube%adopt_facets(facets=facets)
call fig_open(u, 'anim-ref-usage-connect', '', 200._R8P, -25._R8P)
call fig_zoom(u, corner, 4.e-4_R8P)
call fig_anim(u, 2.5_R8P)
call scene('before: facet 1 shifted by 1e-4 along x, its edges disconnected (red)')
call fig_label(u, 'accent', corner, 'neighbours')
call fig_label(u, 'accent', corner + 1.e-4_R8P * ex_R8P, 'facet 1')
call cube%connect_nearby_vertices
call cube%analyze
call scene('after connect_nearby_vertices: the group merged at its centroid')
call fig_close(u)

! the full workflow
c = dragon%get_centroid() ; bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
call dragon%translate(delta=-c)
call dragon%resize(factor=(1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])) * &
                          (ex_R8P + ey_R8P + ez_R8P))
call dragon%analyze
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax() ; c = dragon%get_centroid()
write(caption, '(A,F6.4,A,L1)') 'normalized: longest side ', maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z]), &
                                ', centroid within 1e-12 of the origin ', max(abs(c%x), abs(c%y), abs(c%z)) < 1.e-12_R8P
print '(A)', trim(caption)
call fig_open(u, 'ref-usage-workflow', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', dragon)
call fig_box(u, 'box', bmin, bmax)
call fig_segments(u, 'probe', [vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(0._R8P, 0._R8P, 0._R8P), &
                  vector_R8P(0._R8P, 0._R8P, 0._R8P)], [0.7_R8P * ex_R8P, 0.45_R8P * ey_R8P, 0.6_R8P * ez_R8P])
call fig_label(u, 'probe', 0.7_R8P * ex_R8P, 'x')
call fig_label(u, 'probe', 0.45_R8P * ey_R8P, 'y')
call fig_label(u, 'probe', 0.6_R8P * ez_R8P, 'z')
call fig_points(u, 'accent', [c])
call fig_label(u, 'accent', c, 'centroid = origin')
write(caption, '(A,3F8.4)') 'bmin', bmin%x, bmin%y, bmin%z
call fig_label(u, 'box', bmin, trim(caption))
call fig_close(u)
contains
   subroutine edges_frame(state)
   !< One frame of the validity figure: the dragon, its boundary edges (red) and non-manifold edges (gold).
   character(*), intent(in) :: state !< State of the dragon.

   call dragon%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
   pool => dragon%get_vertex_pool()
   if (allocated(a)) deallocate(a, b)
   allocate(a(size(kinds)), b(size(kinds)))
   do i = 1, size(kinds)
      a(i) = pool%coord(edges(1, i)) ; b(i) = pool%coord(edges(2, i))
   enddo
   write(caption, '(A,A,I0,A,I0,A)') state, ': ', count(kinds == FEATURE_EDGE_BOUNDARY), &
         ' boundary (red), ', count(kinds == FEATURE_EDGE_NONMANIFOLD), ' non-manifold (gold) edges'
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call fig_surface(u, 'ghost', dragon)
   if (count(kinds == FEATURE_EDGE_BOUNDARY) > 0) &
      call fig_segments(u, 'hole', pack(a, kinds == FEATURE_EDGE_BOUNDARY), pack(b, kinds == FEATURE_EDGE_BOUNDARY))
   call fig_segments(u, 'accent', pack(a, kinds == FEATURE_EDGE_NONMANIFOLD), pack(b, kinds == FEATURE_EDGE_NONMANIFOLD))
   endsubroutine edges_frame

   subroutine scene(text)
   !< One frame of the crack: the facets, the boundary edges, the vertices near the corner.
   character(*), intent(in) :: text !< Caption.

   print '(A)', text
   call fig_frame(u, caption=text)
   call fig_surface(u, 'body', cube)
   call fig_boundary(u, 'hole', cube)
   if (allocated(p)) deallocate(p)
   allocate(p(0))
   do i = 1, cube%get_facets_number()
      f => cube%facet_at(i)
      do j = 1, 3
         if (norm2([f%vertex(j)%x - corner%x, f%vertex(j)%y - corner%y, f%vertex(j)%z - corner%z]) < 1.e-3_R8P) &
            p = [p, f%vertex(j)]
      enddo
   enddo
   call fig_points(u, 'accent', p)
   endsubroutine scene
endprogram fig_ref_usage_checks
