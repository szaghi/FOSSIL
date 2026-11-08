!run anim-tut02-stages anim_tut02_stages
!fig anim-tut02-stages
!< Animation of tutorial 2: the passes of sanitize(fill_holes=.true.) run one by one on dragon.stl, one frame each,
!< with what each pass found or changed over a translucent copy of the dragon (decimated to 900 facets, drawing only):
!< the facets on non-manifold edges (gold dots at their centroids), the boundary edges the removal of duplicates opens
!< (ember), the facets the orientation pass flips (teal dots), the patch facets fill_holes adds (ember dots). The counts
!< are those of the full surface.
program anim_tut02_stages
use fossil,                    only : surface_stl_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_facet_object,       only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_fig,                only : fig_open, fig_anim, fig_frame, fig_surface, fig_points, fig_segments, fig_close
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
type(surface_stl_object)      :: surface, light
type(facet_object), pointer   :: f
type(vector_R8P), allocatable :: normal(:)
integer(I4P), allocatable     :: kinds(:)
integer(I4P)                  :: u, i, n0, flipped, filled, status
character(100)                :: caption

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
light = surface
call light%decimate(target_facets=900_I4P, status=status)                   ! the drawing of the body
call fig_open(u, 'anim-tut02-stages', '', -60._R8P, 22._R8P)
call fig_anim(u, 2.2_R8P)
call fig_surface(u, 'ghost', light)                                           ! in every frame

call edge_kinds
write(caption, '(A,I0,A,I0,A)') 'as loaded: ', surface%get_facets_number(), ' facets, ', &
                                count(kinds == FEATURE_EDGE_NONMANIFOLD), ' non-manifold edges, the facets on them'
call frame_start ; call nonmanifold_facets

call surface%remove_degenerate_facets
call surface%connect_nearby_vertices
call surface%analyze
call surface%remove_duplicate_facets
call surface%analyze
call edge_kinds
write(caption, '(A,I0,A,I0,A)') 'remove_duplicate_facets: ', surface%get_duplicate_facets_removed(), ' removed, ', &
                                count(kinds == FEATURE_EDGE_BOUNDARY), ' boundary edges open'
call frame_start ; call boundary_edges

allocate(normal(surface%get_facets_number()))
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i) ; normal(i) = f%normal
enddo
call surface%sanitize_normals(flipped_count=flipped)
call surface%analyze
write(caption, '(A,I0,A)') 'sanitize_normals: ', flipped, ' facets flipped to a consistent outward winding'
call frame_start
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   if ((f%normal .dot. normal(i)) < 0._R8P) call fig_points(u, 'inside', [f%centroid])
enddo

n0 = surface%get_facets_number()
call surface%fill_holes(holes_filled=filled)
call surface%sanitize_normals(flipped_count=flipped)
call surface%analyze
call edge_kinds
write(caption, '(A,I0,A,I0,A,I0,A)') 'fill_holes: ', filled, ' holes closed by ', surface%get_facets_number() - n0, &
                                     ' facets, ', count(kinds == FEATURE_EDGE_BOUNDARY), ' boundary edges left'
call frame_start
do i = n0 + 1, surface%get_facets_number()
   f => surface%facet_at(i) ; call fig_points(u, 'outside', [f%centroid])
enddo

write(caption, '(A,I0,A)') 'left: ', count(kinds == FEATURE_EDGE_NONMANIFOLD), &
                           ' non-manifold edges, sanitize cannot remove them'
call frame_start ; call nonmanifold_facets
call fig_close(u)
contains
   subroutine edge_kinds
   !< Kinds of the boundary and non-manifold edges of the surface now.
   integer(I4P), allocatable :: edges(:,:)

   call surface%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
   endsubroutine edge_kinds

   subroutine frame_start
   !< A frame with the current caption, printed too.
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   endsubroutine frame_start

   subroutine nonmanifold_facets
   !< The facets with an edge that has no single neighbour (all on non-manifold edges here).
   do i = 1, surface%get_facets_number()
      f => surface%facet_at(i)
      if (any(f%fcon_edge == 0_I4P)) call fig_points(u, 'foot', [f%centroid])
   enddo
   endsubroutine nonmanifold_facets

   subroutine boundary_edges
   !< The boundary edges, from feature_edges.
   type(vertex_pool_object), pointer :: pool
   type(vector_R8P), allocatable     :: a(:), b(:)
   integer(I4P), allocatable         :: edges(:,:)
   integer(I4P)                      :: k

   call surface%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
   pool => surface%get_vertex_pool()
   a = [(pool%coord(edges(1, k)), k = 1, size(kinds))]
   b = [(pool%coord(edges(2, k)), k = 1, size(kinds))]
   call fig_segments(u, 'hole', pack(a, kinds == FEATURE_EDGE_BOUNDARY), pack(b, kinds == FEATURE_EDGE_BOUNDARY))
   endsubroutine boundary_edges
endprogram anim_tut02_stages
