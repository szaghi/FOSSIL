!run figref-deviation-2d fig_ref_deviation_2d
!fig ref-deviation-2d
!< Figure of the reference page "Surface deviation", 2-D section: a NACA 0012 section (200 points, cosine spacing) and a
!< coarse copy of it (12 points); curve%deviation of the fine section from the coarse one, its certified bounds in the
!< caption, and the fine section's vertex farthest from the coarse one (largest per_vertex) joined to its closest point.
program fig_ref_deviation_2d
use fossil,     only : curve_object, DEV_STATUS_OK
use fossil_fig, only : fig_open, fig_segment, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter       :: PI = acos(-1._R8P), LIFT = 1.e-3_R8P
type(curve_object)            :: fine, coarse
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P),        allocatable :: uf(:), vf(:), uc(:), vc(:), pv(:)
real(R8P)                     :: lower, upper, foot(2), dist, x(2)
character(96)                 :: caption, text
integer(I4P)                  :: u, i, k, status

call naca(100_I4P, uf, vf)
call naca(6_I4P, uc, vc)
call fine%adopt_loops(u=uf, v=vf, loop_sizes=[size(uf, kind=I4P)])
call coarse%adopt_loops(u=uc, v=vc, loop_sizes=[size(uc, kind=I4P)])
call fine%deviation(coarse, lower=lower, upper=upper, per_vertex=pv, tolerance=1.e-6_R8P, status=status)
if (status /= DEV_STATUS_OK) error stop 'deviation failed'
k = maxloc(pv, dim=1)
x = fine%get_vertex(k)
call coarse%compute_distance(x, distance=dist, closest_point=foot)
print '(A,ES12.5,A,ES12.5)', 'h(fine, coarse): ', lower, ' <= h <= ', upper
print '(A,I0,A,2F9.5,A,ES12.5)', 'farthest vertex ', k, ' at', x, ', distance ', abs(dist)

write(caption, '(A,ES9.3,A,ES9.3)') 'NACA 0012, 200 points against 12: ', lower, ' <= h <= ', upper
call fig_open(u, 'ref-deviation-2d', trim(caption), -90._R8P, 90._R8P, height=300_I4P)
call polyline(uc, vc, 'accent', 0._R8P)
call polyline(uf, vf, 'body', LIFT)
call fig_segment(u, 'patch', p(x(1), x(2), 2 * LIFT), p(foot(1), foot(2), 2 * LIFT))
call fig_points(u, 'patch', [p(x(1), x(2), 2 * LIFT)])
write(text, '(A,ES9.2)') 'max d = ', abs(dist)
call fig_label(u, 'patch', p(x(1), x(2) + 0.03_R8P, 2 * LIFT), trim(text))
call fig_label(u, 'accent', p(0.45_R8P, -0.09_R8P, LIFT), 'coarse (12 points)')
call fig_label(u, 'body', p(0.45_R8P, 0.08_R8P, LIFT), 'fine (200 points)')
call fig_close(u)
contains
   subroutine naca(np, uu, vv)
   !< NACA 0012 of unit chord, np points per side (cosine spacing, closed trailing edge), counter-clockwise.
   integer(I4P),           intent(in)  :: np
   real(R8P), allocatable, intent(out) :: uu(:), vv(:)
   real(R8P)                           :: xs(np + 1), ts(np + 1)
   integer(I4P)                        :: j, n

   xs = [(0.5_R8P * (1._R8P - cos(PI * j / np)), j = 0, np)]
   ts = 0.6_R8P * (0.2969_R8P * sqrt(xs) - 0.1260_R8P * xs - 0.3516_R8P * xs ** 2 + 0.2843_R8P * xs ** 3 - &
                   0.1036_R8P * xs ** 4)
   allocate(uu(2 * np), vv(2 * np))
   n = 0
   do j = np + 1, 2, -1
      n = n + 1 ; uu(n) = xs(j) ; vv(n) = ts(j)
   enddo
   do j = 1, np
      n = n + 1 ; uu(n) = xs(j) ; vv(n) = -ts(j)
   enddo
   endsubroutine naca

   subroutine polyline(uu, vv, style, dz)
   !< A closed polyline.
   real(R8P),    intent(in) :: uu(:), vv(:), dz
   character(*), intent(in) :: style
   integer(I4P)             :: j

   if (allocated(a)) deallocate(a, b)
   allocate(a(size(uu)), b(size(uu)))
   do j = 1, size(uu)
      a(j) = p(uu(j), vv(j), dz)
   enddo
   b = cshift(a, 1)
   call fig_segments(u, style, a, b)
   endsubroutine polyline

   function p(pu, pv_, dz) result(r)
   !< A point of the plane, lifted by dz towards the viewer.
   real(R8P), intent(in) :: pu, pv_, dz
   type(vector_R8P)      :: r

   r = pu * ex_R8P + pv_ * ey_R8P + dz * ez_R8P
   endfunction p
endprogram fig_ref_deviation_2d
