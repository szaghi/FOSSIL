program ck_eikonal
use fossil, only : distance_field_object, surface_stl_object, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)    :: cube
type(distance_field_object) :: field
real(R8P)                   :: x(41), q(41, 41, 41)
integer(I4P)                :: i, flagged, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
x = [(-0.5_R8P + 0.05_R8P * i, i = 0, 40)]
call field%initialize(x=x, y=x, z=x)
call field%compute(surface=cube, status=status)
call field%gradient_diagnostic(q=q, flagged=flagged, status=status)   ! q = |1 - |grad phi||, flagged: q > 0.05
if (status /= STATUS_OK) error stop 'diagnostic failed'
print '(A,I0,A,I0)', 'flagged nodes ', flagged, ' of ', size(q)
print '(A,I0,A,I0)', '  inside ', count(q > 0.05_R8P .and. field%phi < 0._R8P), ', outside ', &
                     count(q > 0.05_R8P .and. field%phi > 0._R8P)
print '(A,I0)', '  in the band ', count(q > 0.05_R8P .and. field%in_band)
print '(A,F7.4,A,F7.4)', 'max q inside ', maxval(q, mask=field%phi < 0._R8P), ', outside ', &
                         maxval(q, mask=field%phi > 0._R8P)
endprogram ck_eikonal
