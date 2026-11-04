!run figck-load_clipped fig_ck_load_clipped
!fig ck-load_clipped
!< Figure of the cookbook recipe "Load only a region of a large file": the clip box, the facets the clipped load keeps
!< (whole, so they reach out of the box) and the rest of the dragon as a ghost (the facets of the full load that the
!< clipped one does not have, matched by their vertices).
program fig_ck_load_clipped
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_surface, fig_box, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
type(surface_stl_object)        :: dragon, head
type(facet_object), pointer     :: f, g
type(facet_object), allocatable :: rest(:)
type(vector_R8P)                :: clip_min, clip_max
logical,            allocatable :: kept(:)
integer(I4P)                    :: u, i, j, k
character(80)                   :: caption

clip_min = vector_R8P(-11._R8P, -5._R8P, 15._R8P) ; clip_max = vector_R8P(10._R8P, 5._R8P, 20._R8P)
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call head%load_from_file(file_name='dragon.stl', guess_format=.true., clip_min=clip_min, clip_max=clip_max)
allocate(kept(dragon%get_facets_number())) ; kept = .false.
do i = 1, dragon%get_facets_number()
   f => dragon%facet_at(i)
   do j = 1, head%get_facets_number()
      g => head%facet_at(j)
      if (all([(f%vertex(k)%x == g%vertex(k)%x .and. f%vertex(k)%y == g%vertex(k)%y .and. &
                f%vertex(k)%z == g%vertex(k)%z, k=1, 3)])) then
         kept(i) = .true. ; exit
      endif
   enddo
enddo
allocate(rest(count(.not.kept))) ; j = 0
do i = 1, dragon%get_facets_number()
   if (kept(i)) cycle
   j = j + 1 ; f => dragon%facet_at(i) ; rest(j) = f
enddo
call dragon%adopt_facets(facets=rest)
write(caption, '(A,I0,A,I0,A)') 'load_from_file(clip_min, clip_max): ', head%get_facets_number(), ' of ', &
                                count(kept) + size(rest), ' facets, kept whole'
call fig_open(u, 'ck-load_clipped', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_surface(u, 'body', head)
call fig_box(u, 'accent', clip_min, clip_max)
call fig_close(u)
endprogram fig_ck_load_clipped
