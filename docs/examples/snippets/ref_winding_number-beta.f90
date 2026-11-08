call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
err = 0._R8P
do k = 1, 5 ; do j = 1, 5 ; do i = 1, 5                       ! 125 points over the dragon's box
   q(1) = bmin + ((i - 0.5_R8P) / 5) * (bmax%x - bmin%x) * vector_R8P(1._R8P, 0._R8P, 0._R8P) + &
                 ((j - 0.5_R8P) / 5) * (bmax%y - bmin%y) * vector_R8P(0._R8P, 1._R8P, 0._R8P) + &
                 ((k - 0.5_R8P) / 5) * (bmax%z - bmin%z) * vector_R8P(0._R8P, 0._R8P, 1._R8P)
   err(1) = max(err(1), abs(dragon%winding_number(point=q(1)) - dragon%winding_number(point=q(1), beta=-1._R8P)))
   err(2) = max(err(2), abs(dragon%winding_number(point=q(1), beta=4._R8P) - &
                            dragon%winding_number(point=q(1), beta=-1._R8P)))
enddo ; enddo ; enddo
print '(A,ES8.2,A,ES8.2)', 'dragon, max |fast - exact|: beta 2 ', err(1), ', beta 4 ', err(2)
