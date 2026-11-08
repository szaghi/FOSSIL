call dragon%tree_boxes(bmin=bmin, bmax=bmax, depth=depth, leaf=leaf)                ! the whole SAH BVH
print '(A,I0,A,I0,A,I0)', 'BVH boxes ', size(bmin), ', leaves ', count(leaf), ', depth ', maxval(depth)
call dragon_octree%tree_boxes(bmin=bmin, bmax=bmax, depth=depth, leaf=leaf)
print '(A,I0,A,I0,A,I0)', 'octree boxes ', size(bmin), ', leaves ', count(leaf), ', depth ', maxval(depth)
