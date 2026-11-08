! a block of 48 x 24 cells around the airfoil, stretched away from y = 0 by 1.1
x = [(-0.35_R8P + 1.2_R8P * i / 48._R8P, i = 0, 48)]
allocate(y(25)) ; y(13) = 0._R8P
do j = 1, 12
   y(13 + j) = y(12 + j) + 0.006_R8P * 1.1_R8P**(j - 1)
   y(13 - j) = -y(13 + j)
enddo
allocate(tags(48, 24), ftags(48, 24))
call section%tag_block(x=x, y=y, band=3, tags=tags, feature_band=2, feature_tags=ftags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
do i = -4, 4
   print '(A,I3,A,I4)', 'tag', i, ': cells', count(tags == i)
enddo
print '(A,I0,A,I0)', 'cells holding a corner ', count(ftags == 0), ', within 2 cells of one ', count(ftags <= 2)
