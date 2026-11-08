program ck_laplacian
use fossil, only : surface_stl_object, csr_matrix_t
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube
type(csr_matrix_t)       :: L, M
real(R8P), allocatable   :: x(:), y(:)
integer(I4P)             :: status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%cotangent_laplacian(L=L, M=M, status=status)  ! over the unique vertices
print '(A,I0,A,I0,A,I0)', 'status ', status, ', vertices ', L%get_nrows(), ', nonzeros of L ', L%get_nnz()
allocate(x(L%get_nrows()), y(L%get_nrows()))
x = 1._R8P
call L%multiply_vector(x, y)                              ! constants are in the kernel: L 1 = 0
print '(A,ES10.2)', 'max |L 1| ', maxval(abs(y))
print '(A,F10.6,A,F10.6)', 'sum of the vertex masses ', sum(M%val), ', surface area ', cube%get_area()
endprogram ck_laplacian
