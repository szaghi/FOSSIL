x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]          ! cell edges of one 32^3 block over the domain
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags)
print '(A)', 'tag      -5     -4     -3     -2     -1      0      1      2      3      4      5'
print '(A,11I7)', 'cells', (count(tags == i), i = -5, 5)
