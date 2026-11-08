!run ck_decimate_tolerance ck_decimate_tolerance
program ck_decimate_tolerance
use fossil, only : surface_stl_object, DEC_STATUS_NO_PROGRESS, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: wing, coarse
type(vector_R8P)         :: diagonal
real(R8P)                :: tol, h
integer(I4P)             :: status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
diagonal = wing%get_bmax() - wing%get_bmin()
tol = 1.e-3_R8P * diagonal%normL2()                                    ! 0.1 % of the bounding-box diagonal
coarse = wing
call coarse%decimate(target_facets=1_I4P, max_error=tol, status=status) ! as few facets as the tolerance allows
if (status /= DEC_STATUS_NO_PROGRESS) error stop 'decimate failed'
call coarse%hausdorff(wing, upper=h, tolerance=tol / 4, status=status)  ! check it: max_error is no certificate
if (status /= DEV_STATUS_OK) error stop 'hausdorff failed'
print '(A,I0,A,I0,A)', 'decimated: ', wing%get_facets_number(), ' -> ', coarse%get_facets_number(), ' facets'
print '(A,ES9.2,A,ES9.2)', 'max_error ', tol, ', certified Hausdorff distance <= ', h
endprogram ck_decimate_tolerance
