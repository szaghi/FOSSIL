!run fig03 fig_tut03
!fig tut03-placed
!fig tut03-clip
!fig tut03-merge
!< Figures of tutorial 3: the dragon repaired as in chapter 2, then the operations of chapter 3: the body centred,
!< turned a quarter about z and resized into the domain cube (x and y edges of the box marked), the
!< clip and its remainder (ghost), the half and its mirror image. It saves fig-dragon-placed.stl, the placed body of the figures of
!< chapters 5 and 6 (they run before the tutorial programs, which write dragon-placed.stl).
program fig_tut03
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface, half, rest, twin
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: scale
integer(I4P)             :: u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! chapter 2
call surface%sanitize(fill_holes=.true.)
call surface%save_into_file(file_name='fig-dragon-clean.stl')
call surface%load_from_file(file_name='fig-dragon-clean.stl', guess_format=.true.) ! chapter 3, place
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%save_into_file(file_name='fig-dragon-placed.stl')
call fig_open(u, 'tut03-placed', 'rotate a quarter turn about z, resize: in the domain [-0.6, 0.6]^3', -60._R8P, 50._R8P)
call fig_surface(u, 'body', surface)
call fig_box(u, 'box', -0.6_R8P * (ex_R8P + ey_R8P + ez_R8P), 0.6_R8P * (ex_R8P + ey_R8P + ez_R8P))
call fig_axes(u, -0.6_R8P * (ex_R8P + ey_R8P + ez_R8P), 0.6_R8P * (ex_R8P + ey_R8P + ez_R8P))
call fig_close(u)
half = surface                                                                     ! clip
call half%clip(bmin=bmin, bmax=bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P, remainder=rest)
call fig_open(u, 'tut03-clip', 'clip: the facets overlapping the box y <= 0, the remainder (ghost)', -60._R8P, 22._R8P)
call fig_surface(u, 'body', half)
call fig_surface(u, 'ghost', rest)
call fig_box(u, 'box', bmin, bmax%x * ex_R8P + 0._R8P * ey_R8P + bmax%z * ez_R8P)
call fig_close(u)
twin = half                                                                        ! mirror and merge
call twin%mirror(normal=ey_R8P)
call twin%analyze
call fig_open(u, 'tut03-merge', 'mirror(normal=ey_R8P) and merge_solids: the half and its mirror image', &
              -60._R8P, 22._R8P)
call fig_surface(u, 'body', half)
call fig_surface(u, 'patch', twin)
call fig_close(u)
contains
   subroutine fig_axes(u, bmin, bmax)
   !< The x and y edges of a box from its corner bmin, labelled (after the box: drawn over its edges).
   integer(I4P),     intent(in) :: u          !< Unit of the scene.
   type(vector_R8P), intent(in) :: bmin, bmax !< Box corners.
   type(vector_R8P)             :: a, b       !< Edge ends.

   a = bmin ; a%x = bmax%x
   b = bmin ; b%y = bmax%y
   call fig_segment(u, 'accent', bmin, a)
   call fig_segment(u, 'accent', bmin, b)
   call fig_label(u, 'accent', a, 'x')
   call fig_label(u, 'accent', b, 'y')
   endsubroutine fig_axes
endprogram fig_tut03
