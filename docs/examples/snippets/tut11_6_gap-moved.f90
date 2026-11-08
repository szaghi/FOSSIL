moved = surface
call moved%move(p1)
call moved%distance_to(box, distance=gap1)
print '(A,L1)', 'pose_self = moved copy, bit for bit: ', gap == gap1
