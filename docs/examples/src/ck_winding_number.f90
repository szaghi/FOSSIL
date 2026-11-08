!run ck_winding_number ck_winding_number
program ck_winding_number
use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: box
type(vector_R8P)         :: p
real(R8P)                :: z(4) = [0.5_R8P, 0.9_R8P, 1.1_R8P, 2._R8P]
integer(I4P)             :: i

! the unit cube without its top face (z = 1): an open box
call box%load_from_file(file_name='cube-01.stl', guess_format=.true., &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), clip_max=vector_R8P(2._R8P, 2._R8P, 0.9_R8P))
print '(A)', '   point (0.5, 0.5, z)   w fast      w exact     inside'
do i = 1, 4
   p = vector_R8P(0.5_R8P, 0.5_R8P, z(i))
   print '(A,F4.1,A,2F12.8,L6)', '   z = ', z(i), '         ', box%winding_number(point=p), &
                                 box%winding_number(point=p, beta=-1._R8P), box%winding_number(point=p) > 0.5_R8P
enddo
endprogram ck_winding_number
