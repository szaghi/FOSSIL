!run figref-winding-number-field fig_ref_winding_number_field
!fig ref-winding-dragon
!fig ref-winding-error
!< Figures of the reference page "Generalized winding number", What it does and Pipeline, on dragon.stl as it is (open:
!< its boundary edges, the rims of its holes, in red):
!< - the winding number w on the plane through the middle of the dragon's box, normal to y, 48 x 32 cells coloured by w
!<   at their centres (cream 0, ember 1, teal below 0): 1 inside the body, 0 outside; the small holes change neither.
!< - on the same plane, the error of the fast evaluation (beta = 2, the default) against the exact sum (beta <= 0),
!<   |w_fast - w_exact| on a log scale from 1e-6 (cream) to 1e-3 (sepia): the far-field expansions err at every point,
!<   inside and outside alike, and stay below 1e-3.
program fig_ref_winding_number_field
use fossil,     only : surface_stl_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter  :: NX = 48, NZ = 32
type(surface_stl_object) :: dragon
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: w(NX, NZ), err(NX, NZ), hx, hz, y, x0, z0, we
integer(I4P)             :: i, k, u
character(100)           :: caption

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
hx = (bmax%x - bmin%x) * 1.1_R8P / NX ; hz = (bmax%z - bmin%z) * 1.1_R8P / NZ
y = 0.5_R8P * (bmin%y + bmax%y)
do k = 1, NZ
   do i = 1, NX
      x0 = 0.5_R8P * (bmin%x + bmax%x) - 0.5_R8P * NX * hx + (i - 0.5_R8P) * hx
      z0 = 0.5_R8P * (bmin%z + bmax%z) - 0.5_R8P * NZ * hz + (k - 0.5_R8P) * hz
      w(i, k) = dragon%winding_number(point=vector_R8P(x0, y, z0))
      we = dragon%winding_number(point=vector_R8P(x0, y, z0), beta=-1._R8P)
      err(i, k) = abs(w(i, k) - we)
   enddo
enddo
print '(A,L1,A,I0)', 'dragon: watertight ', dragon%is_watertight(), ', facets ', dragon%get_facets_number()
print '(A,F7.3,A,F7.3,A,I0,A,I0)', 'w from ', minval(w), ' to ', maxval(w), '; cells with 0.1 < w < 0.9: ', &
                                    count(w > 0.1_R8P .and. w < 0.9_R8P), ' of ', NX * NZ
print '(A,ES9.2,A,I0,A)', 'max |w_fast - w_exact| ', maxval(err), ', ', count(err > 1.e-3_R8P), ' cells above 1e-3'
write(caption, '(A,F6.3,A,F6.3,A,I0,A,I0,A)') 'open dragon: w from ', minval(w), ' to ', maxval(w), '; ', &
                                             count(w > 0.1_R8P .and. w < 0.9_R8P), ' of ', NX * NZ, &
                                             ' cells between 0.1 and 0.9'
call draw('ref-winding-dragon', w, 'field')
write(caption, '(A,ES8.1,A)') '|w_fast - w_exact|, beta = 2: max ', maxval(err), '; log scale 1e-6 to 1e-3'
call draw('ref-winding-error', (log10(max(err, 1.e-6_R8P)) + 6._R8P) / 3._R8P, 'frac')
contains
   subroutine draw(name, v, style)
   !< The plane's cells coloured by v (clipped to the style's range), the dragon translucent, its holes' rims.
   character(*), intent(in) :: name, style
   real(R8P),    intent(in) :: v(:,:)
   real(R8P)                :: lo

   lo = merge(-1._R8P, 0._R8P, style == 'field')
   call fig_open(u, name, trim(caption), -90._R8P, 0._R8P)
   call fig_surface(u, 'ghost', dragon)
   call fig_boundary(u, 'hole', dragon)
   do k = 1, NZ
      do i = 1, NX
         x0 = 0.5_R8P * (bmin%x + bmax%x) - 0.5_R8P * NX * hx + (i - 1) * hx
         z0 = 0.5_R8P * (bmin%z + bmax%z) - 0.5_R8P * NZ * hz + (k - 1) * hz
         call fig_quad(u, style, vector_R8P(x0, y, z0), vector_R8P(x0 + hx, y, z0), vector_R8P(x0 + hx, y, z0 + hz), &
                       vector_R8P(x0, y, z0 + hz), value=max(lo, min(1._R8P, v(i, k))))
      enddo
   enddo
   call fig_close(u)
   endsubroutine draw
endprogram fig_ref_winding_number_field
