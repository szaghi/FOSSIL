!run tut06-sweep tut06_c_sweep
!fig anim-tut06-sweep
!< Tutorial 6, third run: the three stages of compute one by one. After compute_band only the band holds distances;
!< compute_signs gives every other node its sign; sweep fills in the magnitudes, here stopped after 1, 4 and 8
!< directional sweeps (each from the same signed field) and run to convergence. The animation shows phi on the grid
!< plane x = 0 at each stage (colour sign(phi) sqrt(|phi| / 0.6), as the field figure; a node without a magnitude at
!< full colour), the uniform core only; the captions count the nodes still without a magnitude and give the largest
!< difference from the converged field over the others.
program tut06_c_sweep
use fossil_fig, only : fig_open, fig_zoom, fig_anim, fig_frame, fig_quad, fig_close
!region stages
use fossil,     only : surface_stl_object, distance_field_object, SIGN_FAST_WINDING_NUMBER
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter     :: STOPS(3) = [1, 4, 8]
real(R8P),    parameter     :: BIG = 1.e300_R8P             ! above it: no magnitude yet
type(surface_stl_object)    :: surface
type(distance_field_object) :: field, banded, signed, stopped(3)
real(R8P), allocatable      :: x(:), lo(:), hi(:)
real(R8P)                   :: h, xs
integer(I4P)                :: i, m, n, s, sweeps, u
character(96)               :: caption(3)

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)
x = [(-0.6_R8P + 0.04_R8P * i, i = 0, 30)]                    ! the grid of the first run
h = 0.04_R8P
do m = 1, 6
   h = 1.2_R8P * h
   x = [x(1) - h, x, x(size(x)) + h]
enddo
call field%initialize(x=x, y=x, z=x)
call field%compute_band(surface=surface, sign_algorithm=SIGN_FAST_WINDING_NUMBER)   ! exact, |d| <= band
print '(A,I0,A)', 'compute_band:  ', field%get_band_cells(), ' band cells'
banded = field
call field%compute_signs(surface=surface)                         ! flood fill: a sign for every other node
print '(A,I0,A)', 'compute_signs: ', count(field%phi < 0._R8P), ' nodes inside'
signed = field
call field%sweep(sweeps=sweeps)                                   ! fast sweeping to convergence
print '(A,I0,A)', 'sweep:         ', sweeps, ' directional sweeps to converge'
do s = 1, size(STOPS)                                             ! the signed field swept a few times only
   stopped(s) = signed
   call stopped(s)%sweep(max_sweeps=STOPS(s))
   write(caption(s), '(A,I2,A,I0,A,ES8.1)') 'stopped after ', STOPS(s), ' sweeps: ', &
                                            count(abs(stopped(s)%phi) > BIG), ' nodes unknown, error ', &
                                            maxval(abs(stopped(s)%phi - field%phi), abs(stopped(s)%phi) <= BIG)
   print '(A)', trim(caption(s))
enddo
!endregion stages
n = size(x) ; i = minloc(abs(x), dim=1) ; xs = x(i)                ! the plane x = 0, its nodes as dual cells
lo = [x(1), (0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1)]
hi = [(0.5_R8P * (x(m) + x(m + 1)), m = 1, n - 1), x(n)]
call fig_open(u, 'anim-tut06-sweep', '', 0._R8P, 0._R8P, height=500_I4P)
call fig_zoom(u, p(0._R8P, 0._R8P), 0.62_R8P)
call fig_anim(u, 1.2_R8P)
call draw(banded, 'compute_band: exact distances in the band only', band_only=.true.)
call draw(signed, 'compute_signs: every node signed, magnitudes still unknown')
do s = 1, size(STOPS)
   call draw(stopped(s), trim(caption(s)))
enddo
write(caption(1), '(A,I0,A)') 'sweep: converged after ', sweeps, ' directional sweeps'
call draw(field, trim(caption(1)))
call fig_close(u)
contains
   subroutine draw(f, caption, band_only)
   !< A frame: phi of `f` on the plane x = 0, each node as its dual cell; a node without a magnitude yet takes the
   !< colour of its sign at full scale.
   type(distance_field_object), intent(in)           :: f         !< Field.
   character(*),                intent(in)           :: caption   !< Caption.
   logical,                     intent(in), optional :: band_only !< Draw the band nodes only.
   integer(I4P)                                      :: j, k
   real(R8P)                                         :: v

   call fig_frame(u, caption=caption)
   do k = 1, n ; do j = 1, n
      if (max(abs(x(j)), abs(x(k))) > 0.6_R8P + 1.e-9_R8P) cycle          ! the uniform core only
      if (present(band_only)) then
         if (.not. f%in_band(i, j, k)) cycle
      endif
      v = f%phi(i, j, k)
      call fig_quad(u, 'field', p(lo(j), lo(k)), p(hi(j), lo(k)), p(hi(j), hi(k)), p(lo(j), hi(k)), &
                    sign(sqrt(min(1._R8P, abs(v) / 0.6_R8P)), v))
   enddo ; enddo
   endsubroutine draw

   function p(y, z) result(q)
   !< A point of the plane x = xs.
   real(R8P), intent(in) :: y, z !< In-plane coordinates.
   type(vector_R8P)      :: q    !< Point.

   q%x = xs ; q%y = y ; q%z = z
   endfunction p
endprogram tut06_c_sweep
