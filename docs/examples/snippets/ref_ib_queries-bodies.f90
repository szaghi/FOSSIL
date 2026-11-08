use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
type(surface_stl_object)   :: a, b
type(vector_R8P)           :: ghost(3), gap(2), wall(3)
real(R8P)                  :: d(3), probe(3, 3), img(3, 3), second(3, 2), dg(2), dw(3)
real(R8P), allocatable     :: mc(:,:), ma(:), mn(:,:)
real(R8P)                  :: fraction
integer(I4P)               :: bod(3), bw(3), i
logical                    :: vis(3), amb(2)

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)  ! [0, 1]^3
b = a
call b%translate(x=1.2_R8P)                                           ! [1.2, 2.2] x [0, 1]^2
call a%merge_solids(b)
call a%analyze
print '(A,I0,A,I0,A,I0)', 'facets ', a%get_facets_number(), ', shells ', a%get_shells_number(), &
                          ', bodies ', a%get_bodies_number()
