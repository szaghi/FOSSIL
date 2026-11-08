!run figref-self-intersection-fold fig_ref_self_intersection_fold
!fig ref-self-intersection-fold
!fig ref-self-intersection-resolve
!< Figures of the reference page "Self-intersection" on one closed shell that passes through itself: the unit cube with
!< its corner (1, 1, 1) pushed to (0.3, 0.35, -0.3), through its own bottom face. "What it does": the surface drawn
!< translucent, the facets of the reported pairs in ember, the crossing segments in red; it is watertight and one
!< shell. "Pipeline (resolution)": the same surface after resolve_self_intersections (issue #78), with the pairs
!< find_self_intersections still reports on it: none.
program fig_ref_self_intersection_fold
use fossil,     only : surface_stl_object, facet_object, intersection_pair_t, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_triangle, fig_segment, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object)               :: cube
type(facet_object), allocatable        :: facets(:)
type(facet_object), pointer            :: f
type(intersection_pair_t), allocatable :: pairs(:)
type(vector_R8P), parameter            :: CORNER = vector_R8P(1._R8P, 1._R8P, 1._R8P)
character(120)                         :: caption
integer(I4P)                           :: u, i, j, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number()))
do i = 1, size(facets)
   f => cube%facet_at(i) ; facets(i) = f
   do j = 1, 3
      if (abs(facets(i)%vertex(j)%x - 1._R8P) + abs(facets(i)%vertex(j)%y - 1._R8P) + &
          abs(facets(i)%vertex(j)%z - 1._R8P) < 1.e-12_R8P) facets(i)%vertex(j) = vector_R8P(0.3_R8P, 0.35_R8P, -0.3_R8P)
   enddo
   call facets(i)%compute_normal
enddo
call cube%adopt_facets(facets=facets)
call cube%find_self_intersections(pairs=pairs, status=status)
write(caption, '(A,L1,A,I0,A,I0,A)') 'corner pushed through the bottom: watertight ', cube%is_watertight(), &
      ', shells ', cube%get_shells_number(), ', ', size(pairs), ' crossing pairs'
print '(A)', trim(caption)
call scene('ref-self-intersection-fold')
call cube%resolve_self_intersections(status=status)
print '(A,L1)', 'resolve_self_intersections: BOOL_STATUS_OK ', status == BOOL_STATUS_OK
call cube%find_self_intersections(pairs=pairs)
write(caption, '(A,I0,A,L1,A,I0,A,F6.4)') 'resolved: ', cube%get_facets_number(), ' facets, watertight ', &
      cube%is_watertight(), ', ', size(pairs), ' crossing pairs, volume ', cube%get_volume()
print '(A)', trim(caption)
call scene('ref-self-intersection-resolve')
contains
   subroutine scene(name)
   !< The surface translucent, the facets of the pairs in ember, the crossing segments in red.
   character(*), intent(in) :: name !< Figure name.
   type(facet_object), pointer :: g
   logical, allocatable        :: hit(:)
   integer(I4P)                :: m

   allocate(hit(cube%get_facets_number())) ; hit = .false.
   do m = 1, size(pairs)
      hit(pairs(m)%a) = .true. ; hit(pairs(m)%b) = .true.
   enddo
   call fig_open(u, name, trim(caption), -60._R8P, 20._R8P)
   call fig_surface(u, 'ghost', cube)
   do m = 1, cube%get_facets_number()
      if (.not. hit(m)) cycle
      g => cube%facet_at(m)
      call fig_triangle(u, 'patch', g%vertex(1), g%vertex(2), g%vertex(3))
   enddo
   do m = 1, size(pairs)
      call fig_segment(u, 'hole', pairs(m)%p, pairs(m)%q)
   enddo
   call fig_close(u)
   endsubroutine scene
endprogram fig_ref_self_intersection_fold
