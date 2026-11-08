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
