!< Figure of the reference page "Cut-cell geometry", section "What it does": the 2-D outputs on a solid ellipse
!< (semi-axes 0.36 and 0.2, turned by 0.4 rad, a 128-gon) on 8 x 6 cells of 0.125. The cells are shaded by their
!< fluid fraction (cream to sepia); the crossings of every grid line (line_crossings) are gold dots; the fluid
!< intervals of the edges the curve crosses (the partial apertures, from piece_iv_*) are red; the wetted normal
!< integral of every cut cell is a teal stroke from its wetted centroid, drawn at 0.8 times its length.
program fig_ref_cut_cells_outputs
use fossil,     only : curve_object, cut_cells_2d_t
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NX = 8, NY = 6, NS = 128
real(R8P),    parameter       :: H = 0.125_R8P, PI = acos(-1._R8P), A = 0.4_R8P, LIFT = 1.e-3_R8P
type(curve_object)            :: ellipse
type(cut_cells_2d_t)          :: cut
real(R8P)                     :: x(NX + 1), y(NY + 1), u(NS), v(NS), ph, wl, res, xc
real(R8P),        allocatable :: t(:)
integer(I4P),     allocatable :: seg(:)
logical,          allocatable :: ent(:)
type(vector_R8P), allocatable :: s(:), e(:), xs(:)
integer(I4P)                  :: i, j, k, m, n, w, id, nc

x = [(H * i, i = 0, NX)] ; y = [(H * j, j = 0, NY)]
do k = 1, NS
   ph = 2 * PI * (k - 1) / NS
   u(k) = 0.5_R8P + 0.36_R8P * cos(ph) * cos(A) - 0.2_R8P * sin(ph) * sin(A)
   v(k) = 0.38_R8P + 0.36_R8P * cos(ph) * sin(A) + 0.2_R8P * sin(ph) * cos(A)
enddo
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[NS])
call ellipse%cut_cells(x=x, y=y, cut=cut)
allocate(xs(0)) ; nc = 0
do j = 1, NY + 1                                              ! crossings of the u lines, then of the v lines
   call ellipse%line_crossings(axis=1_I4P, a=y(j), t0=x(1), t1=x(NX + 1), t=t, segments=seg, entering=ent, n=n)
   xs = [xs, (q(t(k), y(j), 3 * LIFT), k = 1, n)] ; nc = nc + n
enddo
do i = 1, NX + 1
   call ellipse%line_crossings(axis=2_I4P, a=x(i), t0=y(1), t1=y(NY + 1), t=t, segments=seg, entering=ent, n=n)
   xs = [xs, (q(x(i), t(k), 3 * LIFT), k = 1, n)] ; nc = nc + n
enddo
wl = 0._R8P ; res = 0._R8P
do j = 1, NY
   do i = 1, NX
      wl = wl + cut%wet_length(i, j)
      res = max(res, abs((cut%aperture(2, i, j) - cut%aperture(1, i, j)) * H - cut%wet_normal(1, i, j)), &
                     abs((cut%aperture(4, i, j) - cut%aperture(3, i, j)) * H - cut%wet_normal(2, i, j)))
   enddo
enddo
print '(A,I0,A,I0,A,I0)', 'cells ', NX * NY, ': cut ', count(cut%alpha > 0._R8P .and. cut%alpha < 1._R8P), &
                          ', grid-line crossings ', nc
print '(A,F8.5,A,F8.5)', 'solid area ', sum(1._R8P - cut%alpha) * H**2, ', polygon area ', &
                         0.5_R8P * sum(u * cshift(v, 1) - cshift(u, 1) * v)
print '(A,F8.5,A,ES9.2)', 'wetted length ', wl, ', largest divergence residual ', res

call fig_open(w, 'ref-cut-cells-outputs', 'fluid fraction (shading), crossings (gold), partial apertures (red), '// &
              'wetted normals (teal)', -90._R8P, 90._R8P, width=760_I4P, height=520_I4P)
do j = 1, NY
   do i = 1, NX
      call fig_quad(w, 'frac', q(x(i), y(j), 0._R8P), q(x(i + 1), y(j), 0._R8P), q(x(i + 1), y(j + 1), 0._R8P), &
                    q(x(i), y(j + 1), 0._R8P), value=max(0._R8P, cut%alpha(i, j)))
   enddo
enddo
call fig_segments(w, 'box', [([(q(x(i), y(j), LIFT), i = 1, NX)], j = 1, NY + 1)], &
                  [([(q(x(i + 1), y(j), LIFT), i = 1, NX)], j = 1, NY + 1)])
call fig_segments(w, 'box', [([(q(x(i), y(j), LIFT), j = 1, NY)], i = 1, NX + 1)], &
                  [([(q(x(i), y(j + 1), LIFT), j = 1, NY)], i = 1, NX + 1)])
call fig_segments(w, 'body', [(q(u(k), v(k), 2 * LIFT), k = 1, NS)], &
                  [(q(u(mod(k, NS) + 1), v(mod(k, NS) + 1), 2 * LIFT), k = 1, NS)])
allocate(s(size(cut%piece_iv_edge)), e(size(cut%piece_iv_edge))) ; m = 0
do k = 1, cut%np                                              ! the fluid intervals of the cut pieces' edges
   if (cut%piece_area(k) >= H**2) cycle
   do n = cut%piece_iv_start(k), cut%piece_iv_start(k + 1) - 1
      if (cut%piece_iv_hi(n) - cut%piece_iv_lo(n) >= H) cycle ! whole fluid edges left grey
      id = cut%piece_iv_edge(n) ; m = m + 1
      if (id <= (NX + 1) * NY) then
         i = mod(id - 1, NX + 1) + 1
         s(m) = q(x(i), cut%piece_iv_lo(n), 2 * LIFT) ; e(m) = q(x(i), cut%piece_iv_hi(n), 2 * LIFT)
      else
         id = id - (NX + 1) * NY ; j = (id - 1) / NX + 1
         s(m) = q(cut%piece_iv_lo(n), y(j), 2 * LIFT) ; e(m) = q(cut%piece_iv_hi(n), y(j), 2 * LIFT)
      endif
   enddo
enddo
call fig_segments(w, 'hole', s(1:m), e(1:m))
m = 0
do j = 1, NY
   do i = 1, NX
      if (cut%wet_length(i, j) <= 0._R8P) cycle
      m = m + 1
      s(m) = q(cut%wet_centroid(1, i, j), cut%wet_centroid(2, i, j), 3 * LIFT)
      xc = 0.8_R8P
      e(m) = q(cut%wet_centroid(1, i, j) + xc * cut%wet_normal(1, i, j), &
               cut%wet_centroid(2, i, j) + xc * cut%wet_normal(2, i, j), 3 * LIFT)
   enddo
enddo
call fig_segments(w, 'probe', s(1:m), e(1:m))
call fig_points(w, 'probe', s(1:m))
call fig_points(w, 'foot', xs)
call fig_close(w)
contains
   function q(a_, b_, dz) result(r)
   !< A point of the plane z = 0, lifted by dz towards the viewer.
   real(R8P), intent(in) :: a_, b_, dz !< Coordinates, lift.
   type(vector_R8P)      :: r          !< Point.

   r = vector_R8P(a_, b_, dz)
   endfunction q
endprogram fig_ref_cut_cells_outputs
