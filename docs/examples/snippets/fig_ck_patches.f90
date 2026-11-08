!< Figure of the cookbook recipe "Wall, inlet and outlet from one STL": the box of the recipe with its three patches,
!< the wall translucent, the inlet (z = 0) teal, the outlet (z = 1) ember, and the recipe's three probes joined to
!< their closest points, each labelled with the patch query_many returns.
program fig_ck_patches
use fossil,     only : surface_stl_object, facet_object
use fossil_fig, only : fig_open, fig_triangle, fig_points, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)    :: duct
type(facet_object), pointer :: f
integer(I4P), allocatable   :: inlet(:), outlet(:)
type(vector_R8P)            :: probes(3)
real(R8P)                   :: d(3), foot(3, 3)
integer(I4P)                :: i, u, patch(3)

call duct%load_from_file(file_name='cube-01.stl', guess_format=.true.)        ! the recipe's patches
allocate(inlet(0), outlet(0))
do i = 1, duct%get_facets_number()
   f => duct%facet_at(i)
   if (f%normal%z < -0.5_R8P) inlet = [inlet, i]
   if (f%normal%z > 0.5_R8P) outlet = [outlet, i]
enddo
call duct%set_patch(facets=[(i, i = 1, duct%get_facets_number())], name='wall')
call duct%set_patch(facets=inlet, name='inlet')
call duct%set_patch(facets=outlet, name='outlet')
probes = [vector_R8P(0.5_R8P, 0.5_R8P, -0.2_R8P), vector_R8P(0.5_R8P, 0.5_R8P, 1.3_R8P), &
          vector_R8P(1.4_R8P, 0.5_R8P, 0.5_R8P)]
call duct%query_many(points=probes, distances=d, closest_points=foot, patches=patch)

call fig_open(u, 'ck-patches', 'patches wall, inlet (teal), outlet (ember); each probe takes its closest one', &
              -60._R8P, 8._R8P)
do i = 1, duct%get_facets_number()
   f => duct%facet_at(i)
   select case(trim(duct%patch_name(duct%facet_patch(i))))
   case('inlet')
      call fig_triangle(u, 'inside', f%vertex(1), f%vertex(2), f%vertex(3))
   case('outlet')
      call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
   case default
      call fig_triangle(u, 'ghost', f%vertex(1), f%vertex(2), f%vertex(3))
   endselect
enddo
call fig_points(u, 'probe', probes)
do i = 1, 3
   call fig_segment(u, 'accent', probes(i), vector_R8P(foot(1, i), foot(2, i), foot(3, i)))
   call fig_label(u, 'accent', probes(i), duct%patch_name(patch(i)))
enddo
call fig_close(u)
endprogram fig_ck_patches
