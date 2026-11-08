program ck_project_face
use fossil, only : surface_stl_object, curvilinear_block_t, projection_target_t, project_block_face, &
                   PROJ_CORNER, PROJ_CURVE, PROJ_PATCH
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: cube
type(projection_target_t) :: target
type(curvilinear_block_t) :: block
real(R8P)                 :: x(9, 9, 2), y(9, 9, 2), z(9, 9, 2)
real(R8P), allocatable    :: moved(:,:)
integer(I4P), allocatable :: kinds(:,:)
integer(I4P)              :: i, j

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)   ! [0, 1]^3, sharp edges
call cube%projection_target(target)                                      ! patches, edges, corners
! a 9 x 9 block face near the face x = 1, a little off it and a little wider than it
do j = 1, 9
   do i = 1, 9
      x(i, j, :) = [1.02_R8P, 0.8_R8P] ; y(i, j, :) = -0.01_R8P + 1.02_R8P * (i - 1) / 8
      z(i, j, :) = -0.01_R8P + 1.02_R8P * (j - 1) / 8 + 0.005_R8P * sin(3._R8P * i)
   enddo
enddo
call block%initialize(x, y, z)
call project_block_face(block, face=5_I4P, target=target, tolerance=0.05_R8P, moved=moved, kinds=kinds)
print '(A,I0,A,I0,A,I0)', 'corners ', count(kinds == PROJ_CORNER), ', on edges ', count(kinds == PROJ_CURVE), &
                          ', on the face ', count(kinds == PROJ_PATCH)
print '(A,3F8.4)', 'face corner (1, 1) now at ', block%x(1, 1, 1), block%y(1, 1, 1), block%z(1, 1, 1)
print '(A,3F8.4)', 'edge vertex (5, 1) now at ', block%x(5, 1, 1), block%y(5, 1, 1), block%z(5, 1, 1)
print '(A,F8.4)',  'largest move ', maxval(moved)
endprogram ck_project_face
