call body%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call other%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call other%translate(delta=vector_R8P(0.3_R8P, 0.4_R8P, 0.45_R8P))
call body%merge_solids(other=other)                         ! one surface, two shells passing through each other
call report('merged:   ', body)
