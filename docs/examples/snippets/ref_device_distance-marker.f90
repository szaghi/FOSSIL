use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, SIGN_SOLID_ANGLE, &
                   AABB_TREE_OCTREE, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
integer(I4P), parameter       :: NC = 24                     ! cells per direction
type(surface_stl_object)      :: body, octree
type(vector_R8P)              :: bmin, bmax, h
type(vector_R8P), allocatable :: centre(:)
real(R8P),        allocatable :: d(:), d_host(:), d_pn(:)
real(R8P)                     :: half_diag
integer(I4P)                  :: i, j, k, c, status

call body%load_from_file(file_name='dragon.stl', guess_format=.true.)   ! open: holes, so the winding-number sign
bmin = body%get_bmin() ; bmax = body%get_bmax() ; h = bmax - bmin
bmin = bmin - 0.1_R8P * h ; bmax = bmax + 0.1_R8P * h                    ! the block: the box grown by 10 %
h = (bmax - bmin) / real(NC, R8P)
half_diag = 0.5_R8P * sqrt(h%x**2 + h%y**2 + h%z**2)
allocate(centre(NC**3), d(NC**3), d_host(NC**3), d_pn(NC**3))
c = 0
do k = 1, NC ; do j = 1, NC ; do i = 1, NC
   c = c + 1
   centre(c) = vector_R8P(bmin%x + (i - 0.5_R8P) * h%x, bmin%y + (j - 0.5_R8P) * h%y, bmin%z + (k - 0.5_R8P) * h%z)
enddo ; enddo ; enddo

call body%enter_device                                       ! the tree arrays to the device, once
call body%distance_many_device(points=centre, distances=d, is_signed=.true., is_square_root=.true., &
                               sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'distance_many_device failed'
call body%distance_many_device(points=centre, distances=d_pn, is_signed=.true., is_square_root=.true., &
                               sign_algorithm=SIGN_PSEUDO_NORMAL, status=status)
call body%exit_device                                        ! before the surface changes or dies
call body%distance_many(points=centre, distances=d_host, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells ', NC**3, ': inside ', count(d < -half_diag), ', cut ', &
                               count(abs(d) <= half_diag), ', outside ', count(d > half_diag)
print '(A,I0)', 'device vs host, winding number: mismatches ', count(d /= d_host)
print '(A,I0,A,I0,A)', 'pseudo-normal vs winding number: ', count((d_pn < 0._R8P) .neqv. (d < 0._R8P)), &
                       ' signs differ, ', count(abs(d_pn) /= abs(d)), ' magnitudes differ'
