!< Tutorial 16: the distance batch on the device (OpenACC), checked against the host batch.
program tut16_gpu
use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                    ! points per direction
type(surface_stl_object)      :: surface
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_host(:), d_device(:)
character(14),   parameter    :: names(2) = ['pseudo-normal ', 'winding number']
integer(I4P)                  :: i, j, k, p, s, algorithm(2), status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)  ! SAH BVH, the device tree
allocate(points(n**3), d_host(n**3), d_device(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + 1.2_R8P * (i - 1) / (n - 1)) * ex_R8P + (-0.6_R8P + 1.2_R8P * (j - 1) / (n - 1)) * ey_R8P + &
               (-0.6_R8P + 1.2_R8P * (k - 1) / (n - 1)) * ez_R8P
enddo ; enddo ; enddo

call surface%enter_device                                   ! copy the BVH arrays to the device, once
algorithm = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER]
do s = 1, 2
   call surface%distance_many_device(points=points, distances=d_device, is_signed=.true., is_square_root=.true., &
                                     sign_algorithm=algorithm(s), status=status)
   if (status /= STATUS_OK) error stop 'distance_many_device failed'
   call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                              sign_algorithm=algorithm(s))
   print '(A,A,I0,A,I0)', names(s), ': device vs host mismatches ', &
                             count(d_device /= d_host), ' of ', size(points)
enddo
call surface%exit_device                                    ! release the device copies
endprogram tut16_gpu
