!run ref-usage-io ref_usage_io
!< Usage reference, I/O: load (format guessed or explicit, clipped on the fly), reject non-finite input, save.
program ref_usage_io
!region load
use fossil, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object) :: surface, wing, copy, head
type(vector_R8P)         :: bmin, bmax
integer(I4P)             :: stat, u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)  ! ASCII or binary, from the file size
print '(A,I0)', 'dragon.stl             facets ', surface%get_facets_number()
!endregion load
!region save
call surface%save_into_file(file_name='dragon-out-ascii.stl')                    ! ASCII (the default)
call surface%save_into_file(file_name='dragon-out-binary.stl', is_ascii=.false.) ! binary
!endregion save
!region explicit
call wing%load_from_file(file_name='naca0012-ascii.stl', is_ascii=.true.)     ! explicitly ASCII
call copy%load_from_file(file_name='dragon-out-binary.stl', is_ascii=.false.)  ! explicitly binary
print '(A,I0)', 'naca0012-ascii.stl     facets ', wing%get_facets_number()
print '(A,I0)', 'dragon-out-binary.stl  facets ', copy%get_facets_number()
!endregion explicit
!region clipped
bmin%x = -15._R8P ; bmin%y = -5._R8P ; bmin%z =  0._R8P
bmax%x =   0._R8P ; bmax%y =  5._R8P ; bmax%z = 20._R8P
call head%load_from_file(file_name='dragon.stl', guess_format=.true., clip_min=bmin, clip_max=bmax)
print '(A,I0)', 'clipped load           facets ', head%get_facets_number()
bmin = head%get_bmin() ; bmax = head%get_bmax()
print '(A,F9.4,A,F9.4,A)', 'clipped load           x in [', bmin%x, ',', bmax%x, ']'
!endregion clipped
!region reject
open(newunit=u, file='maybe_bad.stl', action='write')  ! one facet with a NaN coordinate
write(u, '(A)') 'solid bad', ' facet normal 0 0 1', '  outer loop', '   vertex 0 0 0', '   vertex 1 0 0', &
                '   vertex NaN 1 0', '  endloop', ' endfacet', 'endsolid bad'
close(u)
call surface%load_from_file(file_name='maybe_bad.stl', guess_format=.true., status=stat)
if (stat == STATUS_INVALID_INPUT) then
   write(*, '(A)') 'STL contains non-finite coordinates; aborting.'
   stop
endif
!endregion reject
endprogram ref_usage_io
