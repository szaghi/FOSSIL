call source%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube, 8 facets per face
facets => source%facets_ref()
allocate(holed_facets(size(facets) - 1))                                  ! drop facet 1: a hole on the face x = 0
do i = 1, size(holed_facets, kind=I4P)
   holed_facets(i) = facets(i + 1)
enddo
call holed%adopt_facets(facets=holed_facets)
print '(A,I0,A,L1,A,L1)', 'input:   ', holed%get_facets_number(), ' facets, watertight ', holed%is_watertight(), &
                          ', volume ', holed%is_volume()
