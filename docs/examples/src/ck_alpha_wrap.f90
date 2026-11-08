!run ck_alpha_wrap ck_alpha_wrap
program ck_alpha_wrap
use fossil, only : surface_stl_object, AWRAP_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: dragon, wrap
type(vector_R8P)         :: lo, hi
real(R8P)                :: diagonal
integer(I4P)             :: status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)     ! holes, duplicates, non-manifold edges
lo = dragon%get_bmin() ; hi = dragon%get_bmax()
diagonal = norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
call dragon%alpha_wrap(alpha=diagonal / 20._R8P, offset=diagonal / 200._R8P, wrapped=wrap, status=status)
if (status /= AWRAP_STATUS_OK) error stop 'alpha_wrap failed'
print '(A,I0,A,L1,A,L1)', 'input: ', dragon%get_facets_number(), ' facets, watertight ', dragon%is_watertight(), &
                          ', manifold ', dragon%is_manifold()
print '(A,I0,A,L1,A,L1,A,F8.2)', 'wrap:  ', wrap%get_facets_number(), ' facets, watertight ', wrap%is_watertight(), &
                                 ', manifold ', wrap%is_manifold(), ', volume ', wrap%get_volume()
endprogram ck_alpha_wrap
