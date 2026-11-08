program ck_project_points
use fossil, only : surface_stl_object, projection_target_t, classify_points, project_points, PROJ_STATUS_OK
use penf,   only : I4P, R8P
implicit none
character(6), parameter   :: KIND_NAME(0:3) = [character(6) :: 'free', 'patch', 'curve', 'corner'] ! PROJ_*
type(surface_stl_object)  :: cube
type(projection_target_t) :: target
real(R8P)                 :: x(5), y(5), z(5), fx(5), fy(5), fz(5), d(5)
integer(I4P)              :: kinds(5), ids(5), i, status

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%projection_target(target)                           ! patches, sharp-edge curves, corners
x = [0.40_R8P, 0.92_R8P, 0.92_R8P, 1.08_R8P, 0.93_R8P]        ! points a grid generator left near the top face
y = [0.50_R8P, 0.50_R8P, 0.20_R8P, 0.65_R8P, 1.07_R8P]
z = [1.08_R8P, 1.07_R8P, 0.92_R8P, 1.04_R8P, 0.95_R8P]
call classify_points(target, x, y, z, tolerance=0.12_R8P, kinds=kinds, ids=ids)  ! corner, else curve, else patch
call project_points(target, x, y, z, kinds, ids, fx, fy, fz, distances=d, status=status)
if (status /= PROJ_STATUS_OK) error stop 'project_points failed'
do i = 1, 5
   print '(A,I0,A,3F6.2,3A,3F6.2,A,F7.4)', 'point ', i, ' (', x(i), y(i), z(i), ') -> ', KIND_NAME(kinds(i)), &
                                          ' (', fx(i), fy(i), fz(i), '), moved ', d(i)
enddo
endprogram ck_project_points
