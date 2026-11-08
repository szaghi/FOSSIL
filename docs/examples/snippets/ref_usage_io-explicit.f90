call wing%load_from_file(file_name='naca0012-ascii.stl', is_ascii=.true.)     ! explicitly ASCII
call copy%load_from_file(file_name='dragon-out-binary.stl', is_ascii=.false.)  ! explicitly binary
print '(A,I0)', 'naca0012-ascii.stl     facets ', wing%get_facets_number()
print '(A,I0)', 'dragon-out-binary.stl  facets ', copy%get_facets_number()
