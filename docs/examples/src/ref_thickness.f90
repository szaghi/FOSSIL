!run ref_thickness ref_thickness
program ref_thickness
use fossil, only : surface_stl_object, facet_object, THICKNESS_RAY, THICKNESS_MEDIAL, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object)    :: plate, upper, pair
type(facet_object), pointer :: f
real(R8P), allocatable      :: thickness(:), gap(:)
logical, allocatable        :: flat(:)
integer(I4P)                :: invalid, status, i
!region plate
call plate%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call plate%resize(factor=ex_R8P + ey_R8P + 0.1_R8P * ez_R8P) ; call plate%analyze   ! 1 x 1 x 0.1
allocate(flat(plate%get_facets_number()))
do i = 1, size(flat)
   f => plate%facet_at(i)
   flat(i) = abs(f%normal%z) > 0.5_R8P                                     ! the two large faces
enddo
call plate%thickness(thickness=thickness, gap=gap, invalid=invalid, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
call show('ray, 31 rays in 60 degrees', thickness)
call plate%thickness(thickness=thickness, rays=1_I4P)
call show('ray, axial ray only       ', thickness)
call plate%thickness(thickness=thickness, method=THICKNESS_MEDIAL)
call show('medial                    ', thickness)
print '(A,L1,A,I0)', 'gap huge on every facet: ', all(gap == huge(1._R8P)), ', invalid ', invalid
!endregion plate
!region gap
call upper%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call upper%resize(factor=ex_R8P + ey_R8P + 0.1_R8P * ez_R8P)
call upper%translate(delta=vector_R8P(0._R8P, 0._R8P, 0.3_R8P))           ! a second plate 0.2 above
pair = plate
call pair%merge_solids(other=upper)
deallocate(flat) ; allocate(flat(pair%get_facets_number()))
do i = 1, size(flat)
   f => pair%facet_at(i)
   flat(i) = (f%normal%z > 0.5_R8P .and. f%centroid%z < 0.2_R8P) .or. (f%normal%z < -0.5_R8P .and. f%centroid%z > 0.2_R8P)
enddo                                                                       ! the two facing faces
call pair%thickness(gap=gap, method=THICKNESS_MEDIAL)
call show('gap between the plates, medial', gap)
!endregion gap
contains
   subroutine show(label, v)
   character(*), intent(in) :: label
   real(R8P),    intent(in) :: v(:)
   print '(A,A,2F9.5,A,I0,A)', label, ': min, max ', minval(v, mask=flat), maxval(v, mask=flat), ' over ', &
         count(flat), ' facets'
   endsubroutine show
endprogram ref_thickness
