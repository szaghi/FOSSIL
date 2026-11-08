allocate(ones(n), Lones(n), x(n), Lx(n))
ones = 1._R8P
call L%multiply_vector(x=ones, y=Lones)                  ! constants are in the kernel
print '(A,ES9.2)', 'max |L 1|           ', maxval(abs(Lones))
mass = 0._R8P
do i = 1, n
   mass = mass + M%get_value(i, i)
enddo
print '(A,F12.9,A,F12.9)', 'sum of M_ii      ', mass, ', area ', cube%get_area()
do i = 1, cube%get_facets_number()                       ! the x coordinate of every pool vertex
   f => cube%facet_at(i)
   do k = 1, 3
      x(f%vertex_id(k)) = f%vertex(k)%x
   enddo
enddo
call L%multiply_vector(x=x, y=Lx)
print '(A,F12.9,A)', 'x^T L x          ', dot_product(x, Lx), ', exact -4 (area of the faces not normal to x)'
