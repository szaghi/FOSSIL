call surface%save_into_file(file_name='dragon-out-ascii.stl')                    ! ASCII (the default)
call surface%save_into_file(file_name='dragon-out-binary.stl', is_ascii=.false.) ! binary
