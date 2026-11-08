!< Figure of the reference page "Immersed-boundary queries", section "What it does": one `query_many` batch on a curve
!< of two bodies, a unit square (sharp corners) and a 24-gon of radius 0.4. Every grid point within 0.3 of a wall is
!< joined to its foot point; teal where the foot lies inside a segment (the normal is the segment normal), ember where
!< the foot is a vertex (the normal is the vertex bisector): at the square's corners the points fan out around one
!< foot, the region where a normal from a finite-difference gradient of phi is least accurate.
program fig_ref_ib_queries_rich
use fossil,     only : curve_object
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), H = 0.1_R8P, BAND = 0.3_R8P
integer(I4P), parameter       :: NI = 33, NJ = 19, NC = 24
type(curve_object)            :: c
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: u(4 + NC), v(4 + NC), p(2, NI * NJ), d(NI * NJ), foot(2, NI * NJ), nrm(2, NI * NJ), err
integer(I4P)                  :: vtx(NI * NJ), bod(NI * NJ), i, j, k, n, uu
logical                       :: near(NI * NJ)
character(120)                :: caption

u(1:4) = [0._R8P, 1._R8P, 1._R8P, 0._R8P] ; v(1:4) = [0._R8P, 0._R8P, 1._R8P, 1._R8P]
do k = 0, NC - 1
   u(5 + k) = 2._R8P + 0.4_R8P * cos(2 * PI * k / NC) ; v(5 + k) = 0.5_R8P + 0.4_R8P * sin(2 * PI * k / NC)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[4_I4P, NC])
n = 0
do j = 1, NJ
   do i = 1, NI
      n = n + 1 ; p(:, n) = [-0.45_R8P + H * (i - 1), -0.45_R8P + H * (j - 1)]
   enddo
enddo
call c%query_many(points=p, distances=d, closest_points=foot, vertices=vtx, normals=nrm, bodies=bod)
near = abs(d) < BAND
err = 0._R8P
do k = 1, n
   if (near(k) .and. vtx(k) == 0) err = max(err, norm2(p(:, k) - foot(:, k) - d(k) * nrm(:, k)))
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'bodies ', c%get_bodies_number(), ', points within 0.3: ', count(near), &
                               ', body 1: ', count(near .and. bod == 1), ', body 2: ', count(near .and. bod == 2)
print '(A,I0,A,I0,A,I0)', 'feet inside a segment ', count(near .and. vtx == 0), ', at a vertex ', &
                          count(near .and. vtx /= 0), ', of them at the square corners ', &
                          count(near .and. vtx /= 0 .and. bod == 1)
print '(A,ES9.2)', 'feet inside a segment: max |p - foot - d n| ', err

write(caption, '(A,I0,A,I0,A,I0,A)') 'query_many: ', count(near), ' points within 0.3; foot inside a segment (teal) ', &
      count(near .and. vtx == 0), ', at a vertex (ember) ', count(near .and. vtx /= 0)
call fig_open(uu, 'ref-ib-queries-rich', trim(caption), -90._R8P, 90._R8P, width=820_I4P, height=420_I4P)
allocate(a(4 + NC), b(4 + NC))
do k = 1, 4
   a(k) = pt(u(k), v(k), 0._R8P) ; b(k) = pt(u(mod(k, 4) + 1), v(mod(k, 4) + 1), 0._R8P)
enddo
do k = 1, NC
   a(4 + k) = pt(u(4 + k), v(4 + k), 0._R8P) ; b(4 + k) = pt(u(4 + mod(k, NC) + 1), v(4 + mod(k, NC) + 1), 0._R8P)
enddo
call fig_segments(uu, 'accent', a, b)
deallocate(a, b)
a = [(pt(p(1, k), p(2, k), 1.e-3_R8P), k = 1, n)] ; b = [(pt(foot(1, k), foot(2, k), 1.e-3_R8P), k = 1, n)]
call fig_segments(uu, 'box', pack(a, near), pack(b, near))
call fig_points(uu, 'probe', pack(a, near .and. vtx == 0))
call fig_points(uu, 'outside', pack(a, near .and. vtx /= 0))
call fig_label(uu, 'accent', pt(0.5_R8P, 0.5_R8P, 2.e-3_R8P), 'body 1')
call fig_label(uu, 'accent', pt(2._R8P, 0.5_R8P, 2.e-3_R8P), 'body 2')
call fig_close(uu)
contains
   pure function pt(x, y, z) result(q)
   !< A point.
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: q

   q%x = x ; q%y = y ; q%z = z
   endfunction pt
endprogram fig_ref_ib_queries_rich
