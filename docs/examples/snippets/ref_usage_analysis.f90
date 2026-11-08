!< Usage reference, analysis: the statistics of the cube fixture after sanitize, the validity predicates of the raw
!< dragon and of the repaired cube.
program ref_usage_analysis
use fossil, only : surface_stl_object

implicit none
type(surface_stl_object) :: surface
character(16)            :: what

call get_command_argument(1, what)
select case (trim(what))
case ('statistics')
call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
case ('validity')
   call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
   print '(A)', 'dragon.stl as loaded:'
   call report(surface)
   call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call surface%sanitize
   print '(A)', 'cube-01.stl after sanitize:'
   call report(surface)
endselect
contains
   subroutine report(surface)
   !< The validity predicates of a surface.
   type(surface_stl_object), intent(in) :: surface !< Surface.

if (.not. surface%is_watertight()) print '(A)', 'surface has boundary or non-manifold edges'
if (.not. surface%is_manifold())   print '(A)', 'surface has non-manifold edges (3+ incidences)'
if (.not. surface%is_volume()) then
   print '(A)', 'volume is not physically meaningful: do not trust get_volume()/get_centroid()'
endif
   print '(A,L1,A,L1,A,L1)', '  is_watertight ', surface%is_watertight(), '  is_manifold ', surface%is_manifold(), &
                             '  is_volume ', surface%is_volume()
   endsubroutine report
endprogram ref_usage_analysis
