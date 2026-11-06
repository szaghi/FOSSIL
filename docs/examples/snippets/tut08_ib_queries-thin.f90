call surface%alpha_wrap(alpha=0.03_R8P, offset=0.003_R8P, wrapped=wrap)   ! one closed surface (chapter 11)
call wrap%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL)
call wrap%thin_regions(thickness=thickness, gap=gap, h=h, factor=1._R8P, region=region, region_area=area)
print '(A,ES9.2,A,ES9.2,A,F6.3)', 'wrap: thinnest wall ', minval(thickness), ', narrowest gap ', minval(gap), ', h ', h
print '(A,I0,A,I0,A,I0,A,F7.4)', 'under h: ', count(region > 0), ' of ', size(region), ' facets in ', size(area), &
                                 ' regions, the largest of area ', area(1)
