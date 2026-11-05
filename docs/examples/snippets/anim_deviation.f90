!< Figure and animation of the reference page "Surface deviation". The dragon is decimated again and again to half of
!< its facets (quadric edge collapses, down to about 200 facets); each level is coloured by its deviation from the
!< original: each facet by the distance of its centroid to the original (distance_many), over the largest of the level
!< (cream 0, ember the largest), the certified bounds lower <= h <= upper from deviation(original, lower=, upper=) in
!< the caption. The vertices of a quadric decimation stay close to the surface: the deviation peaks inside the
!< facets that bridge the concave parts. The still ref-deviation-dragon is the level of a quarter of the facets.
program anim_deviation
use fossil,              only : surface_stl_object, facet_object, DEC_STATUS_OK, DEV_STATUS_OK
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_surface, fig_triangle, fig_close
use vecfor,              only : vector_R8P
use penf,                only : I4P, R8P

implicit none
integer(I4P), parameter  :: NLEV = 5
type(surface_stl_object) :: original, coarse(NLEV)
real(R8P),   allocatable :: dc(:)
type(vector_R8P), allocatable :: centres(:)
real(R8P)                :: lower(NLEV), upper(NLEV)
character(96)            :: caption
integer(I4P)             :: u, k, status, triangles

call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
coarse(1) = original
do k = 1, NLEV                                                    ! the sweep and its bounds
   if (k > 1) coarse(k) = coarse(k - 1)
   call coarse(k)%decimate(target_facets=coarse(k)%get_facets_number() / 2, status=status)
   if (status /= DEC_STATUS_OK) error stop 'decimate failed'
   call coarse(k)%deviation(original, lower=lower(k), upper=upper(k), triangles=triangles, status=status)
   if (status /= DEV_STATUS_OK) error stop 'deviation did not converge'
   print '(I6,A,ES10.3,A,ES10.3,A,I0,A)', coarse(k)%get_facets_number(), ' facets: ', lower(k), ' <= h <= ', upper(k), &
                                         '  (', triangles, ' triangles)'
enddo

call fig_open(u, 'anim-deviation', '', -60._R8P, 22._R8P)
call fig_anim(u, 1.2_R8P)
call fig_frame(u, caption='dragon.stl as loaded: 6588 facets, h = 0')
call fig_surface(u, 'body', original)
do k = 1, NLEV
   write(caption, '(I0,A,ES9.2,A,ES9.2)') coarse(k)%get_facets_number(), ' facets: ', lower(k), ' <= h <= ', upper(k)
   if (k == NLEV) then
      call fig_frame(u, caption=trim(caption), seconds=2.5_R8P)
   else
      call fig_frame(u, caption=trim(caption))
   endif
   call coloured(u, coarse(k))
enddo
call fig_close(u)

write(caption, '(A,I0,A,ES9.2,A,ES9.2)') 'deviation from the original, ', coarse(2)%get_facets_number(), ' facets: ', &
                                         lower(2), ' <= h <= ', upper(2)
call fig_open(u, 'ref-deviation-dragon', trim(caption), -60._R8P, 22._R8P)
call coloured(u, coarse(2))
call fig_close(u)

contains
   subroutine coloured(unit, s)
   !< The facets of s coloured by the distance of their centroids to the original, over the largest of them.
   integer(I4P),             intent(in) :: unit
   type(surface_stl_object), intent(in) :: s
   real(R8P)                            :: scale
   type(facet_object),       pointer    :: f
   integer(I4P)                         :: i

   if (allocated(centres)) deallocate(centres, dc)
   allocate(centres(s%get_facets_number()), dc(s%get_facets_number()))
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      centres(i) = f%centroid
   enddo
   call original%distance_many(points=centres, distances=dc, is_signed=.false., is_square_root=.true.)
   scale = max(maxval(dc), tiny(1._R8P))
   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      call fig_triangle(unit, 'field', f%vertex(1), f%vertex(2), f%vertex(3), value=min(1._R8P, dc(i) / scale))
   enddo
   endsubroutine coloured
endprogram anim_deviation
