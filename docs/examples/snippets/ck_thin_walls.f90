program ck_thin_walls
use fossil, only : surface_stl_object, THICKNESS_MEDIAL, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)  :: wing
type(vector_R8P)          :: bmin, bmax
real(R8P), allocatable    :: thickness(:), gap(:), area(:)
integer(I4P), allocatable :: region(:)
real(R8P)                 :: h
integer(I4P)              :: level, invalid, status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%thickness(thickness=thickness, gap=gap, method=THICKNESS_MEDIAL, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
bmin = wing%get_bmin() ; bmax = wing%get_bmax()
print '(A,I0,A,L1,A,ES9.2,A,I0)', 'facets ', wing%get_facets_number(), ', watertight ', wing%is_watertight(), &
                                  ', thinnest wall ', minval(thickness), ', invalid ', invalid
do level = 0, 4                                              ! h = chord / 8, / 16, ... / 128
   h = (bmax%x - bmin%x) / 2._R8P ** (level + 3)
   call wing%thin_regions(thickness=thickness, gap=gap, h=h, factor=1._R8P, region=region, region_area=area)
   print '(A,ES9.2,A,I5,A,I3,A,ES9.2)', 'h ', h, ': facets under h ', count(region > 0), ', regions ', size(area), &
                                        ', area ', sum(area)
enddo
endprogram ck_thin_walls
