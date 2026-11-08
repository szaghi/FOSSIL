zc = 0.5_R8P * (x(kc) + x(kc + 1))                      ! the mid-plane of layer 16
call slice_surface(surface, 3_I4P, zc, section, status=status)
call section%tag_block(x=x, y=x, band=4, tags=tags2, status=status)   ! the same 32 x 32 cells, in the plane
print '(A,11I7)', '2-D   ', (count(tags2 == i), i = -5, 5)
print '(A,I0,A,I0)', 'layer 16: cut cells in 3-D ', count(tags(:, :, kc) == 0), ', on its mid-plane section ', &
                     count(tags2 == 0)
print '(A,I0)', 'cells cut in 2-D but not in the layer: ', count(tags2 == 0 .and. tags(:, :, kc) /= 0)
