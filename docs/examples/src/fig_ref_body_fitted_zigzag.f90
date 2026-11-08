!run figref-body-fitted-zigzag fig_ref_body_fitted_zigzag
!fig ref-body-fitted-zigzag
!< Figure of the reference page "Body-fitted grid support", section Why projection needs constraints: a grid line of 13
!< points running along the top front edge of cube-01.stl, just outside it, above the top face and in front of the
!< front face in turn (teal). Projected to the closest point (PROJ_FREE) the feet alternate between the two faces
!< and the line zig-zags across the edge (ember); projected onto the constraint curve (PROJ_CURVE, id 0 = any curve)
!< they lie on the edge (gold).
program fig_ref_body_fitted_zigzag
use fossil,     only : surface_stl_object, projection_target_t, project_points, PROJ_FREE, PROJ_CURVE
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
integer(I4P), parameter   :: N = 13
type(surface_stl_object)  :: cube
type(projection_target_t) :: target
type(vector_R8P)          :: lo, hi, a(N), f(N), g(N), m(N - 1)
real(R8P)                 :: x(N), y(N), z(N), fx(N), fy(N), fz(N), gx(N), gy(N), gz(N), df(N), dg(N), l, e
integer(I4P)              :: i, w, onside

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call cube%projection_target(target)
lo = cube%get_bmin() ; hi = cube%get_bmax() ; l = hi%x - lo%x
print '(A,3F7.3,A,3F7.3,A,I0,A,I0)', 'cube ', lo%x, lo%y, lo%z, ' to ', hi%x, hi%y, hi%z, ', constraint curves ', &
                                     target%nc, ', corners ', target%nk
e = 0.08_R8P * l
do i = 1, N   ! outside, 0.03 L above the top face and 0.08 L behind the edge, or 0.03 L in front and 0.08 L below it
   x(i) = lo%x + l * (0.1_R8P + 0.8_R8P * (i - 1) / (N - 1))
   if (mod(i, 2) == 1) then
      y(i) = lo%y + e ; z(i) = hi%z + 0.03_R8P * l
   else
      y(i) = lo%y - 0.03_R8P * l ; z(i) = hi%z - e
   endif
enddo
call project_points(target, x, y, z, [(PROJ_FREE, i = 1, N)], [(0_I4P, i = 1, N)], fx, fy, fz, distances=df)
call project_points(target, x, y, z, [(PROJ_CURVE, i = 1, N)], [(0_I4P, i = 1, N)], gx, gy, gz, distances=dg)
onside = count(abs(fy - lo%y) < 1.e-12_R8P * l)
print '(A,I0,A,I0,A,F7.4,A,F7.4)', 'free feet on the front face ', onside, ', on the top face ', N - onside, &
                                   ', moved ', maxval(df), ', off the edge by ', maxval(sqrt((fy - lo%y)**2 + (hi%z - fz)**2))
print '(A,L1,A,F7.4)', 'curve feet on the edge ', all(abs(gy - lo%y) + abs(gz - hi%z) < 1.e-12_R8P * l), &
                       ', moved at most ', maxval(dg)
a = [(vector_R8P(x(i), y(i), z(i)), i = 1, N)]
f = [(vector_R8P(fx(i), fy(i), fz(i)), i = 1, N)]
g = [(vector_R8P(gx(i), gy(i), gz(i)), i = 1, N)]
call fig_open(w, 'ref-body-fitted-zigzag', 'grid points near an edge (teal): free feet zig-zag (ember), '// &
              'PROJ_CURVE feet on the edge (gold)', -62._R8P, 28._R8P, width=760_I4P, height=420_I4P)
call fig_zoom(w, vector_R8P(lo%x + 0.5_R8P * l, lo%y, hi%z), 0.36_R8P * l)
call fig_surface(w, 'body', cube)
m = [(vector_R8P(0.5_R8P * (fx(i) + fx(i + 1)), lo%y, hi%z), i = 1, N - 1)]   ! drawn over the edge, on the faces
call fig_segments(w, 'hole', [f(1:N - 1), m], [m, f(2:N)])
call fig_segments(w, 'foot', g(1:N - 1), g(2:N))
call fig_points(w, 'probe', a)
call fig_points(w, 'hole', f)
call fig_points(w, 'foot', g)
call fig_close(w)
endprogram fig_ref_body_fitted_zigzag
