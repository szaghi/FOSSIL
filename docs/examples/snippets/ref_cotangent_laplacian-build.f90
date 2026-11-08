call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%cotangent_laplacian(L=L, M=M, status=status)
if (status /= LAPL_STATUS_OK) error stop 'cotangent_laplacian failed'
n = L%get_nrows()
print '(A,I0,A,I0,A,L1)', 'vertices ', n, ', nonzeros of L ', L%get_nnz(), ', symmetric ', L%is_symmetric(tol=0._R8P)
