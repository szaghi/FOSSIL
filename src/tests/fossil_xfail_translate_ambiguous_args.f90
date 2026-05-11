!< FOSSIL xfail: translate with both delta and x/y/z must error stop (audit #14 D8).
program fossil_xfail_translate_ambiguous_args
use fossil,  only : surface_stl_object
use vecfor,  only : vector_R8P
implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: d
d = vector_R8P(1._8, 0._8, 0._8)
call surface%load_from_file(file_name='fossil_test_translate-delta.stl', is_ascii=.false.)
call surface%translate(x=1._8, delta=d)   ! must error stop
print '(A)', 'ERROR: translate did not stop on ambiguous arguments.'
end program fossil_xfail_translate_ambiguous_args
