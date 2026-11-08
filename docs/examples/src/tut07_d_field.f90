!run tut07-field tut07_d_field
!< Tutorial 7, fourth run: a 2-D signed-distance field of the section, the 2-D twin of chapter 6.
program tut07_d_field
!region field
use fossil, only : curve_object, distance_field_object, STATUS_OK
use penf,   only : I4P, R8P

implicit none
type(curve_object)          :: section
type(distance_field_object) :: field
real(R8P)                   :: x(81)
integer(I4P)                :: i, status

call section%load_from_file(file_name='dragon-section.dat')
x = [(-0.6_R8P + 0.015_R8P * i, i = 0, 80)]
call field%initialize(x=x, y=x)                                   ! no z: a 2-D field
call field%compute(curve=section, status=status)
if (status /= STATUS_OK) error stop '2-D field failed'
print '(A,I0,A,I0)', '2-D field: band cells ', field%get_band_cells(), ', inside ', count(field%phi < 0._R8P)
print '(A,F8.5,A,F8.5)', 'inside area ', count(field%phi < 0._R8P) * 0.015_R8P**2, ', curve area ', section%get_area()
!endregion field
endprogram tut07_d_field
