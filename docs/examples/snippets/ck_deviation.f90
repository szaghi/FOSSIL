program ck_deviation
use fossil, only : surface_stl_object, DEV_STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: original, repaired
real(R8P)                :: lower, mean, rms
integer(I4P)             :: status

call original%load_from_file(file_name='dragon.stl', guess_format=.true.)
repaired = original
call repaired%sanitize(fill_holes=.true.)                       ! duplicates removed, 54 holes patched
call repaired%deviation(original, lower=lower, mean=mean, rms=rms, status=status)
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
print '(A,ES10.3,A,ES10.3,A,ES10.3)', 'repaired -> original: h >= ', lower, ', mean ', mean, ', rms ', rms
call original%deviation(repaired, lower=lower, mean=mean, rms=rms)
print '(A,ES10.3,A,ES10.3,A,ES10.3)', 'original -> repaired: h >= ', lower, ', mean ', mean, ', rms ', rms
endprogram ck_deviation
