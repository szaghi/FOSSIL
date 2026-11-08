use fossil,     only : surface_stl_object, projection_target_t, curvilinear_block_t, project_block_face, &
                       PROJ_PATCH, PROJ_CURVE, PROJ_CORNER
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_label, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P

implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
real(R8P), parameter      :: OFFSET(2) = [0.005_R8P, 0.0005_R8P], TOLERANCE(2) = [0.01_R8P, 0.0014_R8P]
type(surface_stl_object)  :: wing
type(projection_target_t) :: target
type(curvilinear_block_t) :: block(2), block0(2)
real(R8P),    allocatable :: moved(:,:)
integer(I4P), allocatable :: kinds(:,:,:), k1(:,:)
integer(I4P)              :: t, i, status, u
character(8)              :: num

call wing%load_from_file(file_name='naca0012-ascii.stl', guess_format=.true.)
call wing%projection_target(target)
allocate(kinds(41, 21, 2))
block(1) = face_block(0._R8P)
print '(A,2F8.5)', 'last two chordwise steps at the trailing edge: ', block(1)%x(40, 1, 1) - block(1)%x(39, 1, 1), &
                                                                   block(1)%x(41, 1, 1) - block(1)%x(40, 1, 1)
do t = 1, 2
   block(t) = face_block(OFFSET(t))                ! layers k = 1, 2: face 5 is k = 1, the one near the wing
   block0(t) = block(t)
   call project_block_face(block(t), face=5_I4P, target=target, tolerance=TOLERANCE(t), moved=moved, &
                           kinds=k1, status=status)
   kinds(:, :, t) = k1
   print '(A,F7.4,A,F7.4,A,I0)', 'face ', OFFSET(t), ' off the wing, tolerance ', TOLERANCE(t), ': status ', status
   print '(A,I0,A,I0,A,I0,A,F7.4)', '  on corners ', count(kinds(:, :, t) == PROJ_CORNER), ', on curves ', &
         count(kinds(:, :, t) == PROJ_CURVE), ', on the patch ', count(kinds(:, :, t) == PROJ_PATCH), &
         ', largest move ', maxval(moved)
   print '(A,L1,A,L1,A,I0)', '  trailing-edge row on x = 0.75: ', all(block(t)%x(41, :, 1) == 0.75_R8P), &
         ', root and tip rows on z = 0, 1: ', all(block(t)%z(:, 1, 1) == 0._R8P) .and. all(block(t)%z(:, 21, 1) == 1._R8P), &
         ', face edges of zero length: ', zero_edges(block(t))
enddo
