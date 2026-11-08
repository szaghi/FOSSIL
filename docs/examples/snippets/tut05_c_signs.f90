!< Tutorial 5, third run: the same 16^3 batch with the two sign algorithms that have a device path.
program tut05_c_signs
use fossil, only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16                      ! points per direction
type(surface_stl_object)      :: surface
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_pn(:), d_wn(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p

allocate(points(n**3), d_pn(n**3), d_wn(n**3))
h = 1.2_R8P / (n - 1)                                        ! the grid of the first run
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
call surface%distance_many(points=points, distances=d_pn, is_signed=.true., is_square_root=.true.)  ! SIGN_PSEUDO_NORMAL
call surface%distance_many(points=points, distances=d_wn, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I5,A,F10.6)', 'pseudo-normal:  inside ', count(d_pn < 0._R8P), '  min d ', minval(d_pn)
print '(A,I5,A,F10.6)', 'winding number: inside ', count(d_wn < 0._R8P), '  min d ', minval(d_wn)
print '(A,I0,A,ES10.3)', 'signs differ at ', count(d_pn * d_wn < 0._R8P), ' points, max ||d_pn| - |d_wn|| ', &
                        maxval(abs(abs(d_pn) - abs(d_wn)))
endprogram tut05_c_signs
