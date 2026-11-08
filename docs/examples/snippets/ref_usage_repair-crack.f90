allocate(facets(surface%get_facets_number()))
do i = 1, size(facets)
   f => surface%facet_at(i) ; facets(i) = f
enddo
do i = 1, 3                                                  ! shift the first facet: its edges disconnect
   facets(1)%vertex(i) = facets(1)%vertex(i) + 1.e-4_R8P * ex_R8P
enddo
call surface%adopt_facets(facets=facets)                     ! adopt_facets runs analyze
