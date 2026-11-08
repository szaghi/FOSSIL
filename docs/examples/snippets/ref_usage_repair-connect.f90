print '(A)', surface%statistics(prefix='before ')   ! shows disconnected edge counts
call surface%connect_nearby_vertices
call surface%analyze
print '(A)', surface%statistics(prefix='after  ')   ! disconnected edges should be reduced
