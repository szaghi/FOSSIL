call plate%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call plate%resize(factor=ex_R8P + ey_R8P + 0.1_R8P * ez_R8P) ; call plate%analyze   ! 1 x 1 x 0.1
allocate(flat(plate%get_facets_number()))
do i = 1, size(flat)
   f => plate%facet_at(i)
   flat(i) = abs(f%normal%z) > 0.5_R8P                                     ! the two large faces
enddo
call plate%thickness(thickness=thickness, gap=gap, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call show('ray, 31 rays in 60 degrees', thickness)
call plate%thickness(thickness=thickness, rays=1_I4P)
call show('ray, axial ray only       ', thickness)
call plate%thickness(thickness=thickness, method=THICKNESS_MEDIAL)
call show('medial                    ', thickness)
print '(A,L1,A,I0)', 'gap huge on every facet: ', all(gap == huge(1._R8P)), ', invalid ', invalid
