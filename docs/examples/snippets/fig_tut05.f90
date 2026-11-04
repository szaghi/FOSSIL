!< Figures of tutorial 5: the points of the 16^3 grid that distance_many puts inside the placed dragon (fig_tut03
!< writes it), with the pseudo-normal sign and with the fast winding number. The points where the two signs differ
!< are in ember on the side that puts them inside.
program fig_tut05
use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER
use fossil_fig, only : fig_open, fig_surface, fig_box, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 16
type(surface_stl_object)      :: bvh
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_pn(:), d_wn(:)
real(R8P)                     :: h
integer(I4P)                  :: i, j, k, p, u
character(80)                 :: caption

allocate(points(n**3), d_pn(n**3), d_wn(n**3))
h = 1.2_R8P / (n - 1)
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n
   p = p + 1
   points(p) = (-0.6_R8P + (i - 1) * h) * ex_R8P + (-0.6_R8P + (j - 1) * h) * ey_R8P + (-0.6_R8P + (k - 1) * h) * ez_R8P
enddo ; enddo ; enddo
call bvh%load_from_file(file_name='fig-dragon-placed.stl', guess_format=.true.)
call bvh%distance_many(points=points, distances=d_pn, is_signed=.true., is_square_root=.true.)
call bvh%distance_many(points=points, distances=d_wn, is_signed=.true., is_square_root=.true., &
                       sign_algorithm=SIGN_FAST_WINDING_NUMBER)
write(caption, '(A,I0,A,I0,A)') 'SIGN_PSEUDO_NORMAL: ', count(d_pn < 0._R8P), ' points inside, ', &
                                count(d_pn < 0._R8P .and. d_wn >= 0._R8P), ' not by the winding number'
call scene('tut05-pseudo-normal', trim(caption), d_pn, d_wn)
write(caption, '(A,I0,A,I0,A)') 'SIGN_FAST_WINDING_NUMBER: ', count(d_wn < 0._R8P), ' points inside, ', &
                                count(d_wn < 0._R8P .and. d_pn >= 0._R8P), ' not by the pseudo-normal'
call scene('tut05-winding-number', trim(caption), d_wn, d_pn)
contains
   subroutine scene(name, caption, d, other)
   !< The dragon (ghost), the domain and the points inside by `d`: teal if `other` agrees, ember if not.
   character(*), intent(in) :: name, caption !< Figure name and caption.
   real(R8P),    intent(in) :: d(:), other(:) !< Distances, by this sign and by the other.

   call fig_open(u, name, caption, -60._R8P, 22._R8P)
   call fig_surface(u, 'ghost', bvh)
   call fig_box(u, 'box', -0.6_R8P * (ex_R8P + ey_R8P + ez_R8P), 0.6_R8P * (ex_R8P + ey_R8P + ez_R8P))
   call fig_points(u, 'inside', pack(points, d < 0._R8P .and. other < 0._R8P))
   call fig_points(u, 'outside', pack(points, d < 0._R8P .and. other >= 0._R8P))
   call fig_close(u)
   endsubroutine scene
endprogram fig_tut05
