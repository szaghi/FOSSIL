program ref_cotangent_laplacian
use fossil, only : surface_stl_object, facet_object, csr_matrix_t, LAPL_STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)    :: cube
type(facet_object), pointer :: f
type(csr_matrix_t)          :: L, M
real(R8P), allocatable      :: x(:), Lx(:), ones(:), Lones(:)
real(R8P)                   :: mass
integer(I4P)                :: n, i, k, status
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%cotangent_laplacian(L=L, M=M, status=status)
if (status /= LAPL_STATUS_OK) error stop 'cotangent_laplacian failed'
n = L%get_nrows()
print '(A,I0,A,I0,A,L1)', 'vertices ', n, ', nonzeros of L ', L%get_nnz(), ', symmetric ', L%is_symmetric(tol=0._R8P)
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
endprogram ref_cotangent_laplacian
