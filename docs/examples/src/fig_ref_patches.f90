!run figref-patches fig_ref_patches
!fig ref-patches
!fig anim-patches
!< Figures of the reference page "Named patches" (issue #66): a box with its patches (z- inlet in teal, z+ outlet in
!< ember, the sides wall), and an animation: a ring of query points around the box, rising from below to above it,
!< each point coloured by the patch of its closest facet (query_many(patches=)).
program fig_ref_patches
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_points, fig_label, fig_anim, fig_frame, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
character(6), parameter     :: STYLE(3) = ['body  ', 'probe ', 'patch ']
type(surface_stl_object)    :: duct
type(facet_object), pointer :: f
type(vector_R8P)            :: ring(48)
real(R8P)                   :: d(48), z
integer(I4P)                :: patch(48), u, i, k, p, n(3)
integer(I4P), allocatable   :: inlet(:), outlet(:)
character(96)               :: caption

call duct%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(inlet(0), outlet(0))
do i = 1, duct%get_facets_number()
   f => duct%facet_at(i)
   if (f%normal%z < -0.5_R8P) inlet = [inlet, i]
   if (f%normal%z > 0.5_R8P) outlet = [outlet, i]
enddo
call duct%set_patch(facets=[(i, i = 1, duct%get_facets_number())], name='wall')
call duct%set_patch(facets=inlet, name='inlet')
call duct%set_patch(facets=outlet, name='outlet')

call fig_open(u, 'ref-patches', 'wall, inlet (teal, z-) and outlet (ember, z+): three patches of one body', &
              -55._R8P, 25._R8P)
call draw_body
call fig_label(u, 'probe', vector_R8P(0.5_R8P, 0.5_R8P, -0.15_R8P), 'inlet')
call fig_label(u, 'patch', vector_R8P(0.5_R8P, 0.5_R8P, 1.15_R8P), 'outlet')
call fig_label(u, 'accent', vector_R8P(1.15_R8P, 0.5_R8P, 0.5_R8P), 'wall')
call fig_close(u)

call fig_open(u, 'anim-patches', '', -55._R8P, 20._R8P)
call fig_anim(u, 0.4_R8P)
call draw_body
do k = 0, 11
   z = -0.4_R8P + 1.8_R8P * real(k, R8P) / 11._R8P
   do i = 1, 48
      ring(i) = vector_R8P(0.5_R8P + 0.9_R8P * cos(0.1309_R8P * i), 0.5_R8P + 0.9_R8P * sin(0.1309_R8P * i), z)
   enddo
   call duct%query_many(points=ring, distances=d, patches=patch)
   do p = 1, 3
      n(p) = count(patch_is(p))
   enddo
   write(caption, '(A,F5.2,A,I0,A,I0,A,I0)') 'z = ', z, ': closest patch wall ', n(1), ', inlet ', n(2), ', outlet ', n(3)
   call fig_frame(u, trim(caption))
   do p = 1, 3
      call fig_points(u, merge('accent', STYLE(p), p == 1), pack(ring, patch_is(p)))
   enddo
enddo
call fig_close(u)
print '(A)', 'ref-patches, anim-patches written'
contains
   subroutine draw_body
   integer(I4P) :: j, q

   do j = 1, duct%get_facets_number()
      f => duct%facet_at(j)
      q = 1
      if (duct%patch_name(f%patch) == 'inlet') q = 2
      if (duct%patch_name(f%patch) == 'outlet') q = 3
      call fig_triangle(u, trim(STYLE(q)), f%vertex(1), f%vertex(2), f%vertex(3))
   enddo
   endsubroutine draw_body

   function patch_is(q) result(m)
   integer(I4P), intent(in) :: q
   logical                  :: m(48)
   character(6), parameter  :: NM(3) = ['wall  ', 'inlet ', 'outlet']
   integer(I4P)             :: j

   do j = 1, 48
      m(j) = duct%patch_name(patch(j)) == trim(NM(q))
   enddo
   endfunction patch_is
endprogram fig_ref_patches
