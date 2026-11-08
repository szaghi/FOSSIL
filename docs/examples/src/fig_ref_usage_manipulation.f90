!run figref-usage-manipulation fig_ref_usage_manipulation
!fig ref-usage-translate
!fig ref-usage-rotate
!fig ref-usage-mirror
!fig ref-usage-resize
!fig ref-usage-clip
!fig ref-usage-merge
!< Figures of the usage reference, manipulation: each operation of ref_usage_manipulation on the dragon as loaded
!< (ghost), the result solid; the origin, the pivot of rotate, mirror and resize, marked.
program fig_ref_usage_manipulation
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_segment, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: pi = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: dragon, surface, remainder
type(vector_R8P)         :: o, bmin, bmax
integer(I4P)             :: u

o = 0._R8P * ex_R8P
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)

surface = dragon                                                            ! translate
call surface%translate(x=4._R8P, y=3.5_R8P, z=0.5_R8P)
call surface%analyze
call fig_open(u, 'ref-usage-translate', 'translate by (4, 3.5, 0.5): the dragon as loaded (ghost) and moved', &
              -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_surface(u, 'body', surface)
call fig_segment(u, 'accent', dragon%get_centroid(), surface%get_centroid())
call fig_close(u)

surface = dragon                                                            ! rotate
call surface%rotate(axis=ex_R8P, angle=pi / 2._R8P)
call surface%analyze
call fig_open(u, 'ref-usage-rotate', 'rotate(axis=ex_R8P, angle=pi/2): a quarter turn about the x axis', &
              -35._R8P, 18._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_surface(u, 'body', surface)
call fig_segment(u, 'accent', -12._R8P * ex_R8P, 12._R8P * ex_R8P)
call fig_points(u, 'accent', [o])
call fig_label(u, 'accent', 12._R8P * ex_R8P, 'x axis')
call fig_label(u, 'accent', o, 'origin')
call fig_close(u)

surface = dragon                                                            ! mirror
call surface%mirror(normal=ex_R8P + ey_R8P)
call surface%analyze
call fig_open(u, 'ref-usage-mirror', 'mirror(normal=ex_R8P + ey_R8P), seen from above: the plane x + y = 0', &
              -90._R8P, 90._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_surface(u, 'body', surface)
call fig_segment(u, 'accent', vector_R8P(-11._R8P, 11._R8P, 20._R8P), vector_R8P(11._R8P, -11._R8P, 20._R8P))
call fig_label(u, 'accent', vector_R8P(11._R8P, -11._R8P, 20._R8P), 'x + y = 0')
call fig_close(u)

surface = dragon                                                            ! resize
call surface%resize(factor=2._R8P * (ex_R8P + ey_R8P + ez_R8P))
call surface%analyze
call fig_open(u, 'ref-usage-resize', 'resize(factor=2): about the origin, so the dragon also moves away from it', &
              -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', dragon)
call fig_surface(u, 'body', surface)
call fig_points(u, 'accent', [o])
call fig_label(u, 'accent', o, 'origin')
call fig_close(u)

surface = dragon                                                            ! clip
bmin = vector_R8P(-15._R8P, -5._R8P, 0._R8P) ; bmax = vector_R8P(0._R8P, 5._R8P, 20._R8P)
call surface%clip(bmin=bmin, bmax=bmax, remainder=remainder)
call fig_open(u, 'ref-usage-clip', 'clip: the facets overlapping the box (whole), the remainder (ghost)', &
              -60._R8P, 22._R8P)
call fig_surface(u, 'ghost', remainder)
call fig_surface(u, 'body', surface)
call fig_box(u, 'box', bmin, bmax)
call fig_close(u)

call fig_open(u, 'ref-usage-merge', 'merge_solids: dragon_part_1 and dragon_part_2 in one surface', -60._R8P, 22._R8P)
call fig_surface(u, 'body', surface)
call fig_surface(u, 'patch', remainder)
call fig_close(u)
endprogram fig_ref_usage_manipulation
