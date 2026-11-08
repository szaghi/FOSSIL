call surface%save_into_file(file_name='dragon_part_1.stl')
call remainder%save_into_file(file_name='dragon_part_2.stl')
call surface_1%load_from_file(file_name='dragon_part_1.stl', guess_format=.true.)
call surface_2%load_from_file(file_name='dragon_part_2.stl', guess_format=.true.)
call surface_1%merge_solids(other=surface_2)              ! merge_solids re-runs analyze
call surface_1%save_into_file(file_name='dragon_merged.stl')
print '(A,I0,A)', 'merge          ', surface_1%get_facets_number(), ' facets'
