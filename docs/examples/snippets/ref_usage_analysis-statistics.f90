call surface%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call surface%sanitize
print '(A)', surface%statistics()
