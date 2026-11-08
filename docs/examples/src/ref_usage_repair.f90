!run ref-usage-sanitize ref_usage_repair sanitize
!run ref-usage-connect ref_usage_repair connect
!< Usage reference, repair: sanitize on the dragon as it comes, connect_nearby_vertices alone on a cube whose first facet
!< was shifted by 1e-4 (a crack along its three edges).
program ref_usage_repair
use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P

implicit none
type(surface_stl_object)        :: surface
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
integer(I4P)                    :: i
character(16)                   :: what

call get_command_argument(1, what)
select case (trim(what))
case ('sanitize')
!region sanitize
call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
!endregion sanitize
case ('connect')
   call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)
!region crack
   allocate(facets(surface%get_facets_number()))
   do i = 1, size(facets)
      f => surface%facet_at(i) ; facets(i) = f
   enddo
   do i = 1, 3                                                  ! shift the first facet: its edges disconnect
      facets(1)%vertex(i) = facets(1)%vertex(i) + 1.e-4_R8P * ex_R8P
   enddo
   call surface%adopt_facets(facets=facets)                     ! adopt_facets runs analyze
!endregion crack
!region connect
print '(A)', surface%statistics(prefix='before ')   ! shows disconnected edge counts
call surface%connect_nearby_vertices
call surface%analyze
print '(A)', surface%statistics(prefix='after  ')   ! disconnected edges should be reduced
!endregion connect
endselect
endprogram ref_usage_repair
