!< Tutorial 2, second run: what sanitize left. The repaired dragon (tut02_a_repair writes it) has no boundary edge;
!< its facet edges without a single neighbour all lie on non-manifold edges. It is still made of many edge-connected
!< shells: the figure draws the largest one translucent and the others solid.
program tut02_b_left
use fossil_fig,                only : fig_open, fig_triangle, fig_points, fig_close
use fossil,                    only : surface_stl_object, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD
use fossil_facet_object,       only : facet_object
use penf,                      only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface
type(facet_object), pointer :: f
integer(I4P), allocatable   :: edges(:,:), kinds(:), size_of(:)
integer(I4P)                :: i, n_open, u, main
character(96)               :: caption

call surface%load_from_file(file_name='dragon-clean.stl', guess_format=.true.)
call surface%feature_edges(edges=edges, kinds=kinds, angle=180._R8P)
print '(A,I0,A,I0)', 'boundary edges ', count(kinds == FEATURE_EDGE_BOUNDARY), &
                     ', non-manifold edges ', count(kinds == FEATURE_EDGE_NONMANIFOLD)
n_open = 0
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i) ; n_open = n_open + count(f%fcon_edge == 0_I4P)
enddo
print '(A,I0)', 'facet edges without a single neighbour ', n_open
print '(A,I0,A,I0)', 'shells ', surface%get_shells_number(), ', bodies ', surface%get_bodies_number()
allocate(size_of(surface%get_shells_number())) ; size_of = 0                ! facets per shell
do i = 1, surface%get_facets_number()
   size_of(surface%facet_shell(i)) = size_of(surface%facet_shell(i)) + 1
enddo
main = maxloc(size_of, dim=1)
write(caption, '(A,I0,A,I0,A,I0,A)') 'largest shell: ', size_of(main), ' facets; ', size(size_of) - 1, &
                                     ' small shells: ', sum(size_of) - size_of(main), ' facets (ember)'
print '(A)', trim(caption)
call fig_open(u, 'tut02-left', trim(caption), -60._R8P, 22._R8P)
do i = 1, surface%get_facets_number()
   f => surface%facet_at(i)
   call fig_triangle(u, merge('ghost', 'patch', surface%facet_shell(i) == main), f%vertex(1), f%vertex(2), f%vertex(3))
   if (surface%facet_shell(i) /= main) call fig_points(u, 'outside', [f%centroid])     ! small: mark them
enddo
call fig_close(u)
endprogram tut02_b_left
