program ck_mass_properties
use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : ex_R8P, ey_R8P, ez_R8P
implicit none
type(surface_stl_object) :: block
real(R8P)                :: mass, center(3), inertia(3,3), moments(3), axes(3,3)
integer(I4P)             :: status, i

call block%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! the unit cube
call block%resize(factor=0.1_R8P * ex_R8P + 0.2_R8P * ey_R8P + 0.3_R8P * ez_R8P)  ! a 0.1 x 0.2 x 0.3 m block
call block%analyze
call block%mass_properties(mass=mass, center=center, inertia=inertia, principal_moments=moments, &
                           principal_axes=axes, density=7850._R8P, status=status)           ! steel, kg/m^3
if (status /= STATUS_OK) error stop 'not a closed body'
print '(A,F10.4,A)',  'mass              ', mass, ' kg'
print '(A,3F10.4,A)', 'centre of mass    ', center, ' m'
print '(A)',          'inertia tensor about the centre, kg m^2:'
do i = 1, 3
   print '(3F12.6)', inertia(i,:)
enddo
print '(A,3F10.6)',   'principal moments ', moments
call block%mass_properties(inertia=inertia, density=7850._R8P, about=[0._R8P, 0._R8P, 0._R8P])
print '(A,3F10.6)',   'about the corner: diagonal ', (inertia(i,i), i = 1, 3)
endprogram ck_mass_properties
