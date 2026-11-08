!< Usage reference, full workflow: load, repair, centre the centroid at the origin, rescale the longest side to 1, save.
program ref_usage_workflow
use fossil, only : surface_stl_object
use penf,   only : R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
type(surface_stl_object) :: surface
type(vector_R8P)         :: bmin, bmax, c

! 1. Load (analyze runs internally).
call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)

! 2. Run the full repair pipeline.
call surface%sanitize

! 3. Centre at the origin and rescale to fit in a unit box.
c = surface%get_centroid() ; bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-c)
call surface%resize(factor=(1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])) * &
                           (ex_R8P + ey_R8P + ez_R8P))
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax() ; c = surface%get_centroid()
print '(A,3F9.4)', 'box min      ', bmin%x, bmin%y, bmin%z
print '(A,3F9.4)', 'box max      ', bmax%x, bmax%y, bmax%z
print '(A,L1)',    'centroid within 1e-12 of the origin ', max(abs(c%x), abs(c%y), abs(c%z)) < 1.e-12_R8P

! 4. Save.
call surface%save_into_file(file_name='dragon-normalized.stl')
endprogram ref_usage_workflow
