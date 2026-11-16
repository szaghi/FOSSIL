!< Figure of the reference page "Surface loads", What accuracy to expect: the 2-D half of the validation
!< fossil_val_surface_loads. A circular cylinder of radius R = 0.5 in the potential flow of speed U = 1 with circulation
!< Gamma = -0.8 (2 pi U R): drag 0 (d'Alembert), lift -rho U Gamma (Kutta-Joukowski). Its wall is a polygon of N = 32 ..
!< 1024 segments, vertices on the circle, moved along it towards a point (no symmetry). The pressure is
!< p = p_inf + rho (U^2 - |u|^2) / 2 with u the exact velocity on the circle at the radial projection of each point. Two
!< paths: the load samples of order 2 with integrate_loads (FD), and the wall segments of the cut cells on n^2 cells,
!< n = N / 2, with cut_loads_2d (FV), p at their midpoints. The figure plots the four errors against N.
program fig_ref_surface_loads_accuracy
use fossil,     only : curve_object, load_samples_t, load_samples_2d, integrate_loads, cut_cells_2d_t, cut_walls_2d_t, &
                       cut_walls_2d, cut_loads_2d
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P),    parameter :: PI = acos(-1._R8P), R = 0.5_R8P, U = 1._R8P, RHO = 1._R8P, P_INF = 0.7_R8P
real(R8P),    parameter :: GAM = -0.8_R8P * 2._R8P * PI * U * R, LIFT = -RHO * U * GAM, ZC(2) = [0.13_R8P, -0.07_R8P]
integer(I4P), parameter :: NL = 6, ROW(4) = [1, 3, 2, 4]
type(curve_object)      :: wall
type(load_samples_t)    :: ls
type(cut_cells_2d_t)    :: cut
type(cut_walls_2d_t)    :: walls
real(R8P), allocatable  :: p(:), x(:), y(:)
real(R8P)               :: e(NL, 4), nn(NL), f(2), m, h
type(vector_R8P)        :: pts(NL)
character(16)           :: text
integer(I4P)            :: k, i, n, c, uf
character(*), parameter :: STYLE(4) = [character(7) :: 'probe', 'accent', 'outside', 'cell']
character(*), parameter :: LABEL(4) = [character(24) :: 'FD drag', 'FD lift', 'FV drag', 'FV lift']

do k = 1, NL
   n = 32 * 2**(k - 1) ; nn(k) = n
   call polygon(n, wall)
   call load_samples_2d(wall, samples=ls, order=2_I4P)                 ! FD: the load samples of order 2
   p = [(pressure(ls%points(:, i)), i = 1, ls%n)]                      ! p of the true wall at each sample
   call integrate_loads(ls, f, m, pressure=p)
   e(k, 1) = abs(f(1)) / (RHO * U**2 * R) ; e(k, 2) = abs(f(2) - LIFT) / abs(LIFT)
   h = 1.6_R8P / (n / 2)                                               ! FV: n / 2 cells per axis
   x = [(ZC(1) - 0.8_R8P + 0.31_R8P * h + i * h, i = 0, n / 2)]
   y = [(ZC(2) - 0.8_R8P + 0.17_R8P * h + i * h, i = 0, n / 2)]
   call wall%cut_cells(x, y, cut, boundaries=.true.)
   call cut_walls_2d(cut, walls)                                       ! the wall segments of the cut pieces
   p = [(pressure(walls%centroid(:, i)), i = 1, walls%n)]              ! p of the true wall at each midpoint
   call cut_loads_2d(walls, f, m, pressure=p)
   e(k, 3) = abs(f(1)) / (RHO * U**2 * R) ; e(k, 4) = abs(f(2) - LIFT) / abs(LIFT)
enddo
print '(A)', 'cylinder R = 0.5, U = 1, Gamma = -0.8 (2 pi U R): errors |F_x| / (rho U^2 R) and |F_y - L| / |L|'
print '(A)', '     N   FD drag    FD lift    FV drag    FV lift   (FD: order-2 samples; FV: N / 2 cells per axis)'
do k = 1, NL
   print '(I6,4ES11.2)', nint(nn(k)), e(k, :)
enddo
print '(A6,4F11.2)', 'order', (log(e(NL - 1, c) / e(NL, c)) / log(2._R8P), c = 1, 4)

call fig_open(uf, 'ref-surface-loads-accuracy', 'drag and lift errors of the cylinder against N', -90._R8P, 90._R8P)
call fig_segment(uf, 'box', pt(nn(1), 1.e-6_R8P), pt(nn(NL), 1.e-6_R8P))
call fig_segment(uf, 'box', pt(nn(1), 1.e-6_R8P), pt(nn(1), 1.e-2_R8P))
do i = -6, -2                                                          ! y ticks: 1e-6 .. 1e-2
   write(text, '(A,I0)') '1e', i
   call fig_label(uf, 'box', pt(nn(1), 10._R8P**i) - vector_R8P(1.1_R8P, 0._R8P, 0._R8P), trim(text))
   call fig_segment(uf, 'box', pt(nn(1), 10._R8P**i) - vector_R8P(0.2_R8P, 0._R8P, 0._R8P), pt(nn(1), 10._R8P**i))
enddo
do k = 1, NL                                                           ! x ticks: N
   write(text, '(I0)') nint(nn(k))
   call fig_label(uf, 'box', pt(nn(k), 10._R8P**(-6.4_R8P)) - vector_R8P(0.3_R8P, 0._R8P, 0._R8P), trim(text))
enddo
call fig_label(uf, 'box', pt(nn(NL), 10._R8P**(-6.4_R8P)) + vector_R8P(0.9_R8P, 0._R8P, 0._R8P), 'N')
call fig_segment(uf, 'box', pt(nn(2), 2.e-3_R8P), pt(nn(5), 2.e-3_R8P / 64))  ! slope -2
call fig_label(uf, 'box', pt(nn(4), 2.e-3_R8P / 16) + vector_R8P(0.4_R8P, 0.5_R8P, 0._R8P), 'N^-2')
do c = 1, 4
   do k = 1, NL
      pts(k) = pt(nn(k), e(k, c))
   enddo
   call fig_segments(uf, trim(STYLE(c)), pts(1:NL - 1), pts(2:NL)) ; call fig_points(uf, trim(STYLE(c)), pts)
   i = ROW(c)                                                          ! legend, in the order of the curves
   call fig_segment(uf, trim(STYLE(c)), vector_R8P(5.6_R8P, -2.9_R8P - 0.45_R8P * i, 0._R8P), &
                    vector_R8P(6.3_R8P, -2.9_R8P - 0.45_R8P * i, 0._R8P))
   call fig_label(uf, trim(STYLE(c)), vector_R8P(6.5_R8P, -2.9_R8P - 0.45_R8P * i, 0._R8P), trim(LABEL(c)))
enddo
call fig_close(uf)
contains
   real(R8P) function pressure(xp) result(pp)
   !< p_inf + rho (U^2 - |u|^2) / 2, u = conj(dw/dz) on the circle at the radial projection of xp.
   real(R8P), intent(in) :: xp(2) !< Point.
   complex(R8P)          :: z, dw
   real(R8P)             :: d(2)

   d = xp - ZC ; d = R * d / norm2(d)                                  ! the nearest point of the true wall
   z = cmplx(d(1), d(2), R8P)
   dw = U * (1._R8P - R**2 / z**2) - cmplx(0._R8P, 1._R8P, R8P) * GAM / (2._R8P * PI * z)
   pp = P_INF + 0.5_R8P * RHO * (U**2 - abs(dw)**2)
   endfunction pressure

   subroutine polygon(nseg, c)
   !< The circle about ZC as nseg segments, vertices on it, a regular polygon projected radially from 0.3 R off the centre.
   integer(I4P),       intent(in)    :: nseg !< Segments.
   type(curve_object), intent(inout) :: c    !< Curve.
   real(R8P)    :: uu(nseg), vv(nseg), w(2), t
   integer(I4P) :: j

   do j = 1, nseg
      t = 0.1_R8P + 2._R8P * PI * (j - 1) / nseg
      w = [cos(t), sin(t)] + 0.3_R8P * [cos(0.7_R8P), sin(0.7_R8P)] ; w = w / norm2(w)
      uu(j) = ZC(1) + R * w(1) ; vv(j) = ZC(2) + R * w(2)
   enddo
   call c%adopt_loops(u=uu, v=vv, loop_sizes=[nseg])
   endsubroutine polygon

   function pt(nk, ek) result(v)
   !< A plot point: x = 6 log10(N / 32), y = 1.5 log10(e).
   real(R8P), intent(in) :: nk, ek !< Segments, error.
   type(vector_R8P)      :: v

   v = vector_R8P(6._R8P * log10(nk / 32._R8P), 1.5_R8P * log10(ek), 0._R8P)
   endfunction pt
endprogram fig_ref_surface_loads_accuracy
