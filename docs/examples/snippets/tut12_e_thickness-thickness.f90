call mesh%thickness(thickness=thickness, method=THICKNESS_RAY, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call mesh%thin_regions(thickness=thickness, h=0.1_R8P, factor=1._R8P, region=region, region_area=area)
print '(A,I0,A,I0,A,I0)', 'thinner than h = 0.1: ', count(region > 0), ' facets in ', size(area), &
                          ' regions; invalid ', invalid
do i = 1, size(area)
   print '(A,I0,A,I0,A,F7.4,A,F7.4)', '  region ', i, ': ', count(region == i), ' facets, area ', area(i), &
                                      ', thinnest ', minval(thickness, mask=region == i)
enddo
