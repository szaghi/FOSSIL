!run figck-deviation fig_ck_deviation
!fig ck-deviation
!< Figure of the cookbook recipe "Measure how far a repair moved the surface": the dragon after sanitize(fill_holes),
!< the facets it kept in stone and the patches it added in ember, with the boundary edges of the original: the
!< repaired surface leaves the original only across its holes, where the patches span them.
program fig_ck_deviation
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_close
use penf,       only : I4P, R8P
implicit none
type(surface_stl_object) :: original, repaired
integer(I4P)             :: u, kept
character(96)            :: caption

call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
repaired = original
call repaired%remove_duplicate_facets
kept = repaired%get_facets_number()                               ! sanitize appends the patches after these
repaired = original
call repaired%sanitize(fill_holes=.true.)
write(caption, '(A,I0,A)') 'sanitize(fill_holes): ', repaired%get_facets_number() - kept, &
                           ' patch facets (ember) over the holes of the original'
call fig_open(u, 'ck-deviation', trim(caption), -60._R8P, 22._R8P)
call fig_surface(u, 'body', repaired, last=kept)
call fig_surface(u, 'patch', repaired, first=kept + 1_I4P)
call fig_boundary(u, 'hole', original)
call fig_close(u)
endprogram fig_ck_deviation
