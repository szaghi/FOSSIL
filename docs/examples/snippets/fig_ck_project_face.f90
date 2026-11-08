!< Figure of the recipe "Put a block face on a body with sharp edges": the unit cube (ghost), the 9 x 9 face grid before the
!< projection (ghost lines, off the face and wider than it) and after it (accent lines): its corners on the cube's
!< corners, its edge lines on the cube's edges, its interior on the face x = 1.
program fig_ck_project_face
use fossil,     only : surface_stl_object, curvilinear_block_t, projection_target_t, project_block_face
use fossil_fig, only : fig_open, fig_surface, fig_segments, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
type(surface_stl_object)  :: cube
type(projection_target_t) :: target
type(curvilinear_block_t) :: block
real(R8P)                 :: x(9, 9, 2), y(9, 9, 2), z(9, 9, 2)
integer(I4P)              :: i, j, u

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%projection_target(target)
do j = 1, 9
   do i = 1, 9
      x(i, j, :) = [1.02_R8P, 0.8_R8P] ; y(i, j, :) = -0.01_R8P + 1.02_R8P * (i - 1) / 8
      z(i, j, :) = -0.01_R8P + 1.02_R8P * (j - 1) / 8 + 0.005_R8P * sin(3._R8P * i)
   enddo
enddo
call block%initialize(x, y, z)
call fig_open(u, 'ck-project-face', 'project_block_face: corners, edges, face', 25._R8P, 20._R8P)
call fig_surface(u, 'ghost', cube)
call grid_lines(u, 'ghost', block)
call project_block_face(block, face=5_I4P, target=target, tolerance=0.05_R8P)
call grid_lines(u, 'accent', block)
call fig_close(u)
contains
   subroutine grid_lines(un, style, b)
   integer(I4P),              intent(in) :: un
   character(*),              intent(in) :: style
   type(curvilinear_block_t), intent(in) :: b
   type(vector_R8P)                      :: p(144), q(144)
   integer(I4P)                          :: k, m, n

   n = 0
   do m = 1, 9
      do k = 1, 8
         n = n + 1 ; p(n) = vector_R8P(b%x(k, m, 1), b%y(k, m, 1), b%z(k, m, 1))
         q(n) = vector_R8P(b%x(k + 1, m, 1), b%y(k + 1, m, 1), b%z(k + 1, m, 1))
         n = n + 1 ; p(n) = vector_R8P(b%x(m, k, 1), b%y(m, k, 1), b%z(m, k, 1))
         q(n) = vector_R8P(b%x(m, k + 1, 1), b%y(m, k + 1, 1), b%z(m, k + 1, 1))
      enddo
   enddo
   call fig_segments(un, style, p(1:n), q(1:n))
   endsubroutine grid_lines
endprogram fig_ck_project_face
