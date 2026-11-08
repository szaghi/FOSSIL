!run ck_device ck_device
program ck_device
use fossil, only : surface_stl_object, SIGN_FAST_WINDING_NUMBER, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
integer(I4P), parameter  :: n = 20
type(surface_stl_object) :: dragon
type(vector_R8P)         :: points(n**3), lo, hi
real(R8P)                :: d_dev(n**3), d_host(n**3), foot_dev(3, n**3), foot_host(3, n**3)
integer(I4P)             :: id_dev(n**3), id_host(n**3), i, j, k, status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
lo = dragon%get_bmin() ; hi = dragon%get_bmax()
points = [(((vector_R8P(lo%x + (hi%x - lo%x) * i / (n - 1), lo%y + (hi%y - lo%y) * j / (n - 1), &
                        lo%z + (hi%z - lo%z) * k / (n - 1)), i = 0, n - 1), j = 0, n - 1), k = 0, n - 1)]
call dragon%enter_device                                       ! copy the BVH arrays to the device, once
call dragon%query_many_device(points=points, distances=d_dev, closest_points=foot_dev, facets=id_dev, &
                              sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
if (status /= STATUS_OK) error stop 'query_many_device failed'
call dragon%exit_device                                        ! release the device copies
call dragon%query_many(points=points, distances=d_host, closest_points=foot_host, facets=id_host, &
                       sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I0,A,I0)', 'points ', n**3, ', inside ', count(d_dev < 0._R8P)
print '(A,I0,A,I0)', 'device vs host mismatches: distance or facet ', count(d_dev /= d_host .or. id_dev /= id_host), &
                     ', foot coordinates ', count(foot_dev /= foot_host)
endprogram ck_device
