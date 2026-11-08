use fossil, only : surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P

implicit none
integer(I4P), parameter  :: nc = 32                     ! cells per direction of the tagged block
type(surface_stl_object) :: surface
real(R8P),   allocatable :: thickness(:), gap(:)
real(R8P)                :: x(nc + 1)
integer(I4P)             :: tags(nc, nc, nc), ftags(nc, nc, nc), ttags(nc, nc, nc), i, j, k, status, u

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 1.2_R8P * i / nc, i = 0, nc)]
call surface%build_features(angle=60._R8P, status=status) ! creases above 60 degrees, boundary, non-manifold
if (status /= STATUS_OK) error stop 'build_features failed'
call surface%thickness(thickness=thickness, gap=gap, rays=9) ! per facet: shape diameter and gap, 9 rays
call surface%build_thin(thickness=thickness, gap=gap, status=status)
if (status /= STATUS_OK) error stop 'build_thin failed'
call surface%tag_block(x=x, y=x, z=x, band=4, tags=tags, feature_band=2, feature_tags=ftags, &
                       thin_factor=1._R8P, thin_band=2, thin_tags=ttags, status=status)
if (status /= STATUS_OK) error stop 'tag_block failed'
print '(A,I0,A,I0)', 'facets ', surface%get_facets_number(), ', thinner than a cell side (0.0375): ', &
                     count(min(thickness, gap) < x(2) - x(1))
print '(A)', 'cells      at 0      1      2  beyond'
print '(A,4I7)', 'cut     ', count(tags == 0), count(abs(tags) == 1), count(abs(tags) == 2), count(abs(tags) > 2)
print '(A,4I7)', 'feature ', (count(ftags == i), i = 0, 3)
print '(A,4I7)', 'thin    ', (count(ttags == i), i = 0, 3)
print '(A,I0,A,I0,A,I0)', 'cut cells: on a feature edge ', count(tags == 0 .and. ftags == 0), &
      ', on a thin wall only ', count(tags == 0 .and. ftags /= 0 .and. ttags == 0), &
      ', neither ', count(tags == 0 .and. ftags /= 0 .and. ttags /= 0)
