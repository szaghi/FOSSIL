!< Figure of the reference page "Surface loads", section "Moving fields between the wall and the grid": an ellipse of 24
!< segments with the vertex field T = cos 2 theta, drawn as a band along the wall coloured by the values
!< `transfer_field_2d` gives at six points per segment (linear between the two vertices); the ghost points of the
!< closure table on a grid of spacing 0.1 (grey) joined to their feet (gold), where `transfer_to_closure` reads T. Then
!< the way back: T at the load samples (order 2) and `samples_to_vertices`, a linear field exactly, cos 2 theta smoothed.
program fig_ref_surface_loads_transfer
use fossil,     only : curve_object, ib_closure_t, fd_closure_2d, transfer_field_2d, transfer_to_closure, &
                       load_samples_t, load_samples_2d, transfer_to_samples, samples_to_vertices
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter       :: NV = 24, NP = 6
real(R8P),    parameter       :: PI = acos(-1._R8P), BAND = 0.06_R8P, LIFT = 1.e-3_R8P
type(curve_object)            :: el
type(ib_closure_t)            :: cl
type(load_samples_t)          :: ls
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P),        allocatable :: tv(:,:), pts(:,:), val(:,:), tf(:,:), x(:), y(:), ts(:,:), back(:,:), lv(:,:)
integer(I4P),     allocatable :: seg(:)
real(R8P)                     :: u(NV), v(NV), p(2), q(2), d(2), n(2), s0, s1
integer(I4P)                  :: i, e, k, un

do i = 1, NV
   u(i) = cos(2 * PI * (i - 1) / NV) ; v(i) = 0.5_R8P * sin(2 * PI * (i - 1) / NV)
enddo
call el%adopt_loops(u=u, v=v, loop_sizes=[NV])
allocate(tv(1, NV), lv(1, NV))
do i = 1, NV
   tv(1, i) = cos(4 * PI * (i - 1) / NV) ; lv(1, i) = u(i) + 2 * v(i)       ! T = cos 2 theta, and a linear field
enddo
allocate(pts(2, NV * NP), seg(NV * NP), val(1, NV * NP))                    ! six points per segment
do e = 1, NV
   p = [u(e), v(e)] ; q = [u(mod(e, NV) + 1), v(mod(e, NV) + 1)]
   do k = 1, NP
      pts(:, (e - 1) * NP + k) = p + (k - 0.5_R8P) / NP * (q - p) ; seg((e - 1) * NP + k) = e
   enddo
enddo
call transfer_field_2d(el, seg, pts, val, vertex_field=tv)
allocate(x(26), y(16))
x = [(-1.25_R8P + 0.1_R8P * i, i = 0, 25)] ; y = [(-0.75_R8P + 0.1_R8P * i, i = 0, 15)]
call fd_closure_2d(el, x, y, [2_I4P, 2_I4P], cl)                           ! the ghost points of the grid and their feet
allocate(tf(1, cl%n))
call transfer_to_closure(el, cl, tf, vertex_field=tv)

call fig_open(un, 'ref-surface-loads-transfer', &
              'T = cos 2 theta on the wall, read at the feet (gold) of the ghosts (grey)', &
              -90._R8P, 90._R8P)
do e = 1, NV
   p = [u(e), v(e)] ; q = [u(mod(e, NV) + 1), v(mod(e, NV) + 1)] ; d = q - p
   n = [d(2), -d(1)] / sqrt(sum(d * d))                                       ! outward
   do k = 1, NP
      s0 = (k - 1._R8P) / NP ; s1 = real(k, R8P) / NP
      call fig_quad(un, 'field', vec(p + s0 * d), vec(p + s1 * d), vec(p + s1 * d + BAND * n), vec(p + s0 * d + BAND * n), &
                    val(1, (e - 1) * NP + k))
   enddo
enddo
allocate(a(NV), b(NV))
do e = 1, NV
   a(e) = vec([u(e), v(e)], LIFT) ; b(e) = vec([u(mod(e, NV) + 1), v(mod(e, NV) + 1)], LIFT)
enddo
call fig_segments(un, 'body', a, b)
deallocate(a, b) ; allocate(a(cl%n), b(cl%n))
do i = 1, cl%n
   a(i) = vec(cl%points(:, i), 2 * LIFT) ; b(i) = vec(cl%foot(:, i), 2 * LIFT)
enddo
call fig_segments(un, 'probe', a, b)
call fig_points(un, 'ghost', a)
call fig_points(un, 'foot', b)
call fig_close(un)
print '(A,I0,A,I0,A,F8.5,A,F8.5)', 'ellipse of ', NV, ' segments: ', cl%n, ' ghost points; T at their feet in ', &
      minval(tf), ' .. ', maxval(tf)
call load_samples_2d(el, samples=ls, order=2_I4P)                          ! the way back, through the load samples
allocate(ts(1, ls%n), back(1, NV))
call transfer_to_samples(el, ls, ts, vertex_field=lv)
call samples_to_vertices(el, ls, ts, back)
print '(A,L1)', 'linear field, vertices -> samples -> vertices, equal to 1e-14: ', maxval(abs(back - lv)) < 1.e-14_R8P
call transfer_to_samples(el, ls, ts, vertex_field=tv)
call samples_to_vertices(el, ls, ts, back)
print '(A,L1)', 'cos 2 theta from the vertices (linear on each segment), the same round trip, unchanged to 1e-14: ', &
      maxval(abs(back - tv)) < 1.e-14_R8P
do i = 1, ls%n                                                             ! cos 2 theta itself at the samples
   ts(1, i) = cos(2 * atan2(ls%points(2, i) / 0.5_R8P, ls%points(1, i)))
enddo
call samples_to_vertices(el, ls, ts, back)
print '(A,F8.5,A)', 'cos 2 theta at the samples, back to the vertices: largest error ', maxval(abs(back - tv)), &
      ' (a smoothing projection)'

contains
   function vec(w, z) result(pt)
   !< A point of the drawing plane.
   real(R8P), intent(in)           :: w(2) !< Coordinates.
   real(R8P), intent(in), optional :: z    !< Lift.
   type(vector_R8P)                :: pt   !< Point.

   pt = vector_R8P(w(1), w(2), 0._R8P)
   if (present(z)) pt%z = z
   endfunction vec
endprogram fig_ref_surface_loads_transfer
