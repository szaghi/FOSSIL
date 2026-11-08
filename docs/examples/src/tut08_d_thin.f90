!run tut08-thin tut08_d_thin
!fig tut08-thin
!< Tutorial 8, fourth run: walls and gaps thinner than a grid spacing, measured on a closed wrap of the dragon before
!< any grid. The figure draws the wrap with the facets of the thin regions in ember.
program tut08_d_thin
use fossil_fig,          only : fig_open, fig_triangle, fig_close
!region thin
use fossil,              only : surface_stl_object, THICKNESS_MEDIAL
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface, wrap
type(facet_object), pointer :: f
real(R8P), allocatable      :: thickness(:), gap(:), area(:)
integer(I4P), allocatable   :: region(:)
real(R8P)                   :: h
integer(I4P)                :: i, u
character(96)               :: caption

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
h = 1.2_R8P / 40._R8P                                                     ! the spacing of the 41^3 grid
call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap)   ! one closed surface (alpha wrapping)
call wrap%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL)
call wrap%thin_regions(thickness=thickness, gap=gap, h=h, factor=1._R8P, region=region, region_area=area)
print '(A,ES9.2,A,ES9.2,A,F6.3)', 'wrap: thinnest wall ', minval(thickness), ', narrowest gap ', minval(gap), ', h ', h
print '(A,I0,A,I0,A,I0,A,F7.4)', 'under h: ', count(region > 0), ' of ', size(region), ' facets in ', size(area), &
                                 ' regions, the largest of area ', area(1)
!endregion thin
write(caption, '(A,I0,A,I0,A)') 'the wrap: ', count(region > 0), ' of ', size(region), &
                                ' facets thinner than h = 0.03 (ember)'
call fig_open(u, 'tut08-thin', trim(caption), -60._R8P, 22._R8P)
do i = 1, wrap%get_facets_number()
   f => wrap%facet_at(i)
   call fig_triangle(u, merge('patch', 'body ', region(i) > 0), f%vertex(1), f%vertex(2), f%vertex(3))
enddo
call fig_close(u)
endprogram tut08_d_thin
