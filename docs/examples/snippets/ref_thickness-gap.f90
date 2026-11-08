call upper%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call upper%resize(factor=ex_R8P + ey_R8P + 0.1_R8P * ez_R8P)
call upper%translate(delta=vector_R8P(0._R8P, 0._R8P, 0.3_R8P))           ! a second plate 0.2 above
pair = plate
call pair%merge_solids(other=upper)
deallocate(flat) ; allocate(flat(pair%get_facets_number()))
do i = 1, size(flat)
   f => pair%facet_at(i)
   flat(i) = (f%normal%z > 0.5_R8P .and. f%centroid%z < 0.2_R8P) .or. (f%normal%z < -0.5_R8P .and. f%centroid%z > 0.2_R8P)
enddo                                                                       ! the two facing faces
call pair%thickness(gap=gap, method=THICKNESS_MEDIAL)
call show('gap between the plates, medial', gap)
