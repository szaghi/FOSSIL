!< FOSSIL xfail: resize with both factor and x/y/z must error stop (audit #14 D8).
program fossil_xfail_resize_ambiguous_args
use fossil,  only : surface_stl_object
use vecfor,  only : vector_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: f
f = vector_R8P(2._8, 2._8, 2._8)
call surface%load_from_file(file_name='fossil_test_resize-factor.stl', is_ascii=.false.)
call surface%resize(x=2._8, factor=f)   ! must error stop
print '(A)', 'ERROR: resize did not stop on ambiguous arguments.'
end program fossil_xfail_resize_ambiguous_args
