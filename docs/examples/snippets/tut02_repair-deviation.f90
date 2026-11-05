call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
call surface%deviation(original, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'repaired -> original ', lower, ' <= h <= ', upper
call original%deviation(surface, lower=lower, upper=upper)
print '(A,ES10.3,A,ES10.3)', 'original -> repaired ', lower, ' <= h <= ', upper
