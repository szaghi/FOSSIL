!< Figure of the cookbook recipe "Swept volume of a moving body on a 3-D grid": the box of ck_swept_cells_3d at its
!< start pose (translucent) and end pose (stone), seen from above, over the cell layer 0.415 <= z <= 0.498 coloured by
!< its swept volume (ember: fluid gained where the box left, teal: fluid lost where it arrived).
program fig_ck_swept_cells_3d
use fossil, only : surface_stl_object, rigid_pose_t, pose_from_axis_angle, swept_cells_3d_t, swept_cells_3d, STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_quad, fig_close
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
type(surface_stl_object) :: box
type(swept_cells_3d_t)   :: sw
real(R8P)                :: x(13)
integer(I4P)             :: i, j, status, u
type(rigid_pose_t)       :: p0, p1
type(surface_stl_object) :: b0, b1
real(R8P)                :: s
character(90)            :: caption

call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call box%translate(delta=vector_R8P(-0.5_R8P, -0.5_R8P, -0.5_R8P))
call box%resize(factor=vector_R8P(0.3_R8P, 0.1_R8P, 0.2_R8P)) ; call box%analyze   ! 0.3 x 0.1 x 0.2 about its origin
x = [(0.083_R8P * i, i = 0, 12)]                              ! 12^3 cells on [0, 0.996]^3
p0 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0._R8P, [0.4_R8P, 0.45_R8P, 0.47_R8P])
p1 = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], PI / 9, [0.55_R8P, 0.5_R8P, 0.52_R8P])
call swept_cells_3d(box, p0, p1, x, x, x, sw, status=status)
if (status /= STATUS_OK) error stop 'swept_cells_3d failed'
b0 = box ; call b0%move(p0) ; b1 = box ; call b1%move(p1)
s = maxval(abs(sw%swept(:, :, 6)))
write(caption, '(A,F9.6)') 'swept volume, layer 0.415 < z < 0.498, colour to +-', s
call fig_open(u, 'ck-swept-cells-3d', trim(caption), -90._R8P, 90._R8P, height=440_I4P)
call fig_zoom(u, vector_R8P(0.48_R8P, 0.48_R8P, 0._R8P), 0.28_R8P)
do j = 1, 12 ; do i = 1, 12
   call fig_quad(u, 'field', vector_R8P(x(i), x(j), 0._R8P), vector_R8P(x(i + 1), x(j), 0._R8P), &
                 vector_R8P(x(i + 1), x(j + 1), 0._R8P), vector_R8P(x(i), x(j + 1), 0._R8P), value=sw%swept(i, j, 6) / s)
enddo ; enddo
call fig_surface(u, 'ghost', b0)
call fig_surface(u, 'body', b1)
call fig_close(u)
endprogram fig_ck_swept_cells_3d
