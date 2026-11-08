program ck_curve_corners
use fossil, only : surface_stl_object, curve_object, slice_surface, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: wing
type(curve_object)        :: section
integer(I4P), allocatable :: corners(:)
real(R8P)                 :: uv(2)
integer(I4P)              :: i, status

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call slice_surface(wing, 3_I4P, 0.5_R8P, section)                    ! the airfoil section
call section%feature_corners(corners=corners, angle=30._R8P, status=status)  ! turning by more than 30 degrees
if (status /= STATUS_OK) error stop 'feature_corners failed'
print '(A,I0,A,I0)', 'vertices ', section%get_vertices_number(), ', corners ', size(corners)
do i = 1, size(corners)
   uv = section%get_vertex(corners(i))
   print '(A,I0,A,2F9.5)', 'corner at vertex ', corners(i), ': ', uv
enddo
endprogram ck_curve_corners
