!run figref-winding-number fig_ref_winding_number
!fig ref-winding-cube
!< Figure of the reference page "Winding number": the unit cube without its top (cube-01.stl clipped at z = 0.9, as
!< the cookbook's sign-with-holes recipe does; its open rim in ember) and the winding number on the slice y = 0.5,
!< one cell of a 48 x 48 grid over [-0.25, 1.25]^2 in x-z coloured by w at its centre (style field: ember towards +1,
!< cream 0, teal towards -1; w < 0 just outside the walls below the rim, where the missing top is seen from its inner
!< side).
program fig_ref_winding_number
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_boundary, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter  :: N = 48
real(R8P),    parameter  :: LO = -0.25_R8P, H = 1.5_R8P / N, Y = 0.5_R8P
type(surface_stl_object) :: box
real(R8P)                :: x0, z0, w, wc, wmin, wmax
character(80)            :: caption
integer(I4P)             :: u, i, k

call box%load_from_file(file_name='cube-01.stl', guess_format=.true.,          &
                        clip_min=vector_R8P(-1._R8P, -1._R8P, -1._R8P), &
                        clip_max=vector_R8P( 2._R8P,  2._R8P, 0.9_R8P))
wc = box%winding_number(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.45_R8P))
write(caption, '(A,F5.3,A)') 'winding_number on y = 0.5 (teal -1, ember +1): w = ', wc, ' at the centre'
call fig_open(u, 'ref-winding-cube', trim(caption), -75._R8P, 25._R8P)
call fig_boundary(u, 'hole', box)
wmin = huge(1._R8P) ; wmax = -huge(1._R8P)
do k = 1, N ; do i = 1, N
   x0 = LO + (i - 1) * H ; z0 = LO + (k - 1) * H
   w = box%winding_number(point=vector_R8P(x0 + 0.5_R8P * H, Y, z0 + 0.5_R8P * H))
   wmin = min(wmin, w) ; wmax = max(wmax, w)
   call fig_quad(u, 'field', vector_R8P(x0, Y, z0), vector_R8P(x0 + H, Y, z0), vector_R8P(x0 + H, Y, z0 + H), &
                 vector_R8P(x0, Y, z0 + H), value=w)
enddo ; enddo
call fig_close(u)
print '(A,2F9.5)', 'w range on the slice: ', wmin, wmax
endprogram fig_ref_winding_number
