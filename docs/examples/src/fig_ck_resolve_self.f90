!run figck-resolve-self fig_ck_resolve_self
!fig ck-resolve-self-folded
!fig ck-resolve-self-resolved
!< Figures of the cookbook recipe "Resolve the self-intersections of one surface": the unit cube with its corner pushed
!< through its own bottom face (translucent, the crossing segments in ember), and what resolve_self_intersections
!< returns (issue #78): the region the surface winds around once, the pocket below the bottom face gone.
program fig_ck_resolve_self
use fossil,     only : surface_stl_object, facet_object, intersection_pair_t
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)               :: cube
type(facet_object), allocatable        :: facets(:)
type(facet_object), pointer            :: f
type(intersection_pair_t), allocatable :: pairs(:)
integer(I4P)                           :: u, i, j
character(90)                          :: caption

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number()))
do i = 1, size(facets)
   f => cube%facet_at(i) ; facets(i) = f
   do j = 1, 3
      if (facets(i)%vertex(j)%x == 1._R8P .and. facets(i)%vertex(j)%y == 1._R8P .and. facets(i)%vertex(j)%z == 1._R8P) &
         facets(i)%vertex(j) = vector_R8P(0.3_R8P, 0.35_R8P, -0.3_R8P)
   enddo
   call facets(i)%compute_normal
enddo
call cube%adopt_facets(facets=facets)
call cube%find_self_intersections(pairs=pairs)
write(caption, '(A,I0,A)') 'corner pushed through the bottom: ', size(pairs), ' crossing pairs (ember)'
call fig_open(u, 'ck-resolve-self-folded', trim(caption), -60._R8P, 20._R8P)
call fig_surface(u, 'ghost', cube)
call fig_segments(u, 'hole', pairs%p, pairs%q)
call fig_close(u)
call cube%resolve_self_intersections
write(caption, '(A,I0,A,F6.4,A)') 'resolve_self_intersections: ', cube%get_facets_number(), ' facets, volume ', &
                                  cube%get_volume(), ', watertight'
call fig_open(u, 'ck-resolve-self-resolved', trim(caption), -60._R8P, 20._R8P)
call fig_surface(u, 'body', cube)
call fig_close(u)
endprogram fig_ck_resolve_self
