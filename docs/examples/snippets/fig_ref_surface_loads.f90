!< Figures of the reference page "Surface loads", the load samples (`load_samples_2d`, `load_samples_3d`):
!< - samples: an ellipse of 40 segments, the three Gauss points of every segment and its outward normal scaled by the
!<   weight;
!< - rules: three triangles, the points of the rules of order 1, 2 and 3;
!< - probes: two squares 0.25 apart, the probes 0.1 and 0.5 out of every sample, the blocked ones in ember, and the
!<   solver's cells enclosing the probes of the top side;
!< - sphere: a marching-cubes sphere, the samples of order 2 of its lower half with their normals;
!< - force: the ellipse in the pressure p = 1 - v (a fluid at rest under gravity -v): the contributions -p n w of the
!<   samples, which add up to the buoyancy (0, area).
program fig_ref_surface_loads
use fossil,     only : surface_stl_object, curve_object, facet_object, extract_isosurface, load_samples_t, &
                       load_samples_2d, load_samples_3d
use fossil_fig, only : fig_open, fig_zoom, fig_segments, fig_points, fig_surface, fig_quad, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter            :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)              :: ellipse, squares
type(surface_stl_object)        :: tris, sphere
type(load_samples_t)            :: ls
type(facet_object), allocatable :: facets(:)
type(vector_R8P),   allocatable :: a(:), b(:), p1(:), p2(:)
real(R8P)                       :: u(40), v(40), f(17, 17, 17), x(17), y(9), force(2), area
integer(I4P)                    :: un, i, j, k, s, o, ii, jj, nblk

u = [(0.4_R8P * cos(2 * PI * (i - 1) / 40), i = 1, 40)] ; v = [(0.2_R8P * sin(2 * PI * (i - 1) / 40), i = 1, 40)]
call ellipse%adopt_loops(u=u, v=v, loop_sizes=[40_I4P])

! ---- samples: three Gauss points per segment, normals scaled by the weights
call load_samples_2d(ellipse, samples=ls, order=3_I4P)
call fig_open(un, 'ref-surface-loads-samples', 'segments (gold), three Gauss points each (teal), normal x weight (ember)', &
              -90._R8P, 90._R8P)
call curve_segments(u, v, a, b) ; call fig_segments(un, 'accent', a, b)
deallocate(a, b) ; allocate(a(ls%n), b(ls%n), p1(ls%n))
do s = 1, ls%n
   a(s) = vec(ls%points(1, s), ls%points(2, s), LIFT)
   b(s) = vec(ls%points(1, s) + 3 * ls%weight(s) * ls%normal(1, s), ls%points(2, s) + 3 * ls%weight(s) * ls%normal(2, s), &
              LIFT)
   p1(s) = vec(ls%points(1, s), ls%points(2, s), 2 * LIFT)
enddo
call fig_segments(un, 'patch', a, b) ; call fig_points(un, 'probe', p1)
call fig_close(un)
print '(A,I0,A,I0,A,F10.7,A,F10.7)', 'ellipse: segments ', ls%ne, ', samples ', ls%n, ', sum of weights ', sum(ls%weight), &
      ', perimeter ', ellipse%get_perimeter()

! ---- rules: three triangles, the points of order 1, 2, 3
allocate(facets(3))
do k = 1, 3
   call tri(facets(k), [1.2_R8P * (k - 1), 0._R8P], [1.2_R8P * (k - 1) + 1._R8P, 0._R8P], [1.2_R8P * (k - 1) + 0.35_R8P, 0.9_R8P])
enddo
call tris%adopt_facets(facets=facets)
call fig_open(un, 'ref-surface-loads-rules', 'the rules of order 1, 2 and 3 on a facet: 1, 3 and 6 points (teal)', -90._R8P, &
              90._R8P, height=260_I4P)
call fig_surface(un, 'body', tris)
deallocate(p1) ; allocate(p1(0))
do o = 1, 3
   call load_samples_3d(tris, samples=ls, order=o)
   do s = 1, ls%n
      if (ls%element(s) /= o) cycle
      p1 = [p1, vector_R8P(ls%points(1, s), ls%points(2, s), 2 * LIFT)]
      print '(A,I0,A,2F9.5,A,F9.6)', 'order ', o, ': point', ls%points(1:2, s) - [1.2_R8P * (o - 1), 0._R8P], &
            '  weight / area ', ls%weight(s) / 0.45_R8P
   enddo
enddo
call fig_points(un, 'probe', p1)
call fig_close(un)

! ---- probes: two squares, probes 0.1 and 0.5 out, the blocked ones in ember, the cells of the top side's probes
call squares%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P, 1.25_R8P, 2.25_R8P, 2.25_R8P, 1.25_R8P], &
                         v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
x = [(-0.5_R8P + 0.1875_R8P * i, i = 0, 16)] ; y = [(-0.5_R8P + 0.25_R8P * j, j = 0, 8)]
call load_samples_2d(squares, x, y, ls, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P])
call fig_open(un, 'ref-surface-loads-probes', 'samples (gold), probes 0.1 and 0.5 out (teal), blocked (ember), cells', &
              -90._R8P, 90._R8P)
do s = 1, ls%n
   if (ls%element(s) /= 3) cycle                                ! the top side of the first square
   do k = 1, ls%m
      ii = ls%probe_cell(1, k, s) ; jj = ls%probe_cell(2, k, s)
      call fig_quad(un, 'frac', vec(x(ii), y(jj), -LIFT), vec(x(ii + 1), y(jj), -LIFT), vec(x(ii + 1), y(jj + 1), -LIFT), &
                    vec(x(ii), y(jj + 1), -LIFT), 0.6_R8P)
   enddo
enddo
deallocate(a, b, p1) ; allocate(a(0), b(0), p1(0), p2(0))
do j = 1, size(y) ; do i = 1, size(x)
   p1 = [p1, vec(x(i), y(j))]
enddo ; enddo
call fig_points(un, 'box', p1)
call fig_segments(un, 'accent', [vec(0._R8P, 0._R8P), vec(1._R8P, 0._R8P), vec(1._R8P, 1._R8P), vec(0._R8P, 1._R8P), &
                                 vec(1.25_R8P, 0._R8P), vec(2.25_R8P, 0._R8P), vec(2.25_R8P, 1._R8P), vec(1.25_R8P, 1._R8P)], &
                  [vec(1._R8P, 0._R8P), vec(1._R8P, 1._R8P), vec(0._R8P, 1._R8P), vec(0._R8P, 0._R8P), &
                   vec(2.25_R8P, 0._R8P), vec(2.25_R8P, 1._R8P), vec(1.25_R8P, 1._R8P), vec(1.25_R8P, 0._R8P)])
deallocate(p1) ; allocate(p1(0))
nblk = 0
do s = 1, ls%n
   p1 = [p1, vec(ls%points(1, s), ls%points(2, s), 2 * LIFT)]
   do k = ls%m, 1, -1
      if (ls%visible(k, s)) then
         a = [a, vec(ls%points(1, s), ls%points(2, s), LIFT)] ; b = [b, vec(ls%probe(1, k, s), ls%probe(2, k, s), LIFT)]
      else
         nblk = nblk + 1
         p2 = [p2, vec(ls%probe(1, k, s), ls%probe(2, k, s), 2 * LIFT)]
      endif
   enddo
enddo
call fig_segments(un, 'probe', a, b) ; call fig_points(un, 'patch', p2) ; call fig_points(un, 'foot', p1)
call fig_close(un)
print '(A,I0,A,I0,A,I0)', 'squares: samples ', ls%n, ', probes ', ls%n * ls%m, ', blocked ', nblk

! ---- sphere: samples of order 2 and normals, lower half
do k = 1, 17 ; do j = 1, 17 ; do i = 1, 17
   f(i, j, k) = norm2([i - 9._R8P, j - 9._R8P, k - 9._R8P] / 8._R8P) - 0.5_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=f, bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(1._R8P, 1._R8P, 1._R8P), &
                        iso=0._R8P, surface=facets)
call sphere%adopt_facets(facets=facets)
call load_samples_3d(sphere, samples=ls, order=2_I4P, lengths=[0.08_R8P])
call fig_open(un, 'ref-surface-loads-sphere', 'samples of order 2 (gold) and their facet normals (ember), lower half', &
              -60._R8P, 25._R8P)
call fig_surface(un, 'ghost', sphere)
deallocate(a, b, p1) ; allocate(a(0), b(0), p1(0))
do s = 1, ls%n
   if (ls%points(3, s) > 0._R8P) cycle
   p1 = [p1, vector_R8P(ls%points(1, s), ls%points(2, s), ls%points(3, s))]
   a = [a, vector_R8P(ls%points(1, s), ls%points(2, s), ls%points(3, s))]
   b = [b, vector_R8P(ls%probe(1, 1, s), ls%probe(2, 1, s), ls%probe(3, 1, s))]
enddo
call fig_segments(un, 'patch', a, b) ; call fig_points(un, 'foot', p1)
call fig_close(un)
print '(A,I0,A,I0,A,F9.6,A,F9.6)', 'sphere: facets ', ls%ne, ', samples ', ls%n, ', area ', sum(ls%weight), &
      ', 4 pi r^2 ', PI

! ---- force: -p n w in p = 1 - v, summed: the buoyancy
call load_samples_2d(ellipse, samples=ls, order=1_I4P)
call fig_open(un, 'ref-surface-loads-force', 'the ellipse in p = 1 - v: the contributions -p n w (ember) add up to (0, area)', &
              -90._R8P, 90._R8P)
call curve_segments(u, v, a, b) ; call fig_segments(un, 'accent', a, b)
deallocate(a, b) ; allocate(a(ls%n), b(ls%n))
force = 0._R8P
do s = 1, ls%n
   force = force - (1._R8P - ls%points(2, s)) * ls%normal(:, s) * ls%weight(s)
   a(s) = vec(ls%points(1, s), ls%points(2, s), LIFT)
   b(s) = vec(ls%points(1, s) - 4 * (1._R8P - ls%points(2, s)) * ls%normal(1, s) * ls%weight(s), &
              ls%points(2, s) - 4 * (1._R8P - ls%points(2, s)) * ls%normal(2, s) * ls%weight(s), LIFT)
enddo
call fig_segments(un, 'patch', a, b)
call fig_close(un)
area = ellipse%get_area()
print '(A,2F11.7,A,F10.7,A,L1)', 'force:', anint(force * 1.e7_R8P) / 1.e7_R8P + 0._R8P, '  area', area, &
      '  equal to 1e-14: ', &
      abs(force(1)) < 1.e-14_R8P .and. abs(force(2) - area) < 1.e-14_R8P

contains
   function vec(x, y, z) result(p)
   !< A point of the drawing plane.
   real(R8P), intent(in)           :: x, y !< Coordinates.
   real(R8P), intent(in), optional :: z    !< Lift.
   type(vector_R8P)                :: p    !< Point.

   p = vector_R8P(x, y, 0._R8P)
   if (present(z)) p%z = z
   endfunction vec

   subroutine curve_segments(u, v, a, b)
   !< The segments of a closed polyline.
   real(R8P),                     intent(in)  :: u(:), v(:) !< Vertices.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Segment ends.
   integer(I4P) :: i, n                                     !< Counters.

   n = size(u) ; allocate(a(n), b(n))
   do i = 1, n
      a(i) = vec(u(i), v(i)) ; b(i) = vec(u(mod(i, n) + 1), v(mod(i, n) + 1))
   enddo
   endsubroutine curve_segments

   subroutine tri(fo, p, q, r)
   !< A facet in the plane z = 0, counter-clockwise.
   type(facet_object), intent(inout) :: fo         !< Facet.
   real(R8P),          intent(in)    :: p(2), q(2), r(2) !< Vertices.

   fo%vertex(1) = vector_R8P(p(1), p(2), 0._R8P) ; fo%vertex(2) = vector_R8P(q(1), q(2), 0._R8P)
   fo%vertex(3) = vector_R8P(r(1), r(2), 0._R8P)
   call fo%compute_metrix
   endsubroutine tri
endprogram fig_ref_surface_loads
