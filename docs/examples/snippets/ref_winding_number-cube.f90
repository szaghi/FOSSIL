call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
q = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(0.5_R8P, 0.5_R8P, 0._R8P), vector_R8P(0.5_R8P, 0._R8P, 0._R8P), &
     vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(2._R8P, 2._R8P, 2._R8P)]
do i = 1, 5
   print '(A,A,F12.9)', WHERE(i), ' w = ', cube%winding_number(point=q(i))
enddo
