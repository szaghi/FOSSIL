!< Figure of the reference page "Immersed-boundary queries", section "Distance to the nearest feature edge": a closed
!< cylinder (radius 0.5, height 1, 48 sectors) whose only feature edges at 30 degrees are the two rims (the side seams
!< turn by 7.5 degrees). `build_features`, then one `query_many` on a 48 x 48 grid of the plane y = 0 with
!< `feature_distances=`. Left half: the wall distance |phi|; right half: the distance to the nearest feature edge, on
!< the same scale (cream 0, sepia 0.5 and above). |phi| vanishes along the whole side and the caps; the feature distance
!< only at the rims, so it marks where the wall stops being smooth.
program fig_ref_ib_queries_feature
use fossil,              only : surface_stl_object, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_quad, fig_segments, fig_points, fig_label, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P), R = 0.5_R8P, H = 2._R8P / 48, SCALE = 0.5_R8P
integer(I4P), parameter         :: NS = 48, N = 48
type(surface_stl_object)        :: cyl
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: p(N * N), a(2 * NS + 4), b(2 * NS + 4)
integer(I4P), allocatable       :: edges(:,:)
real(R8P)                       :: d(N * N), fd(N * N), x, z, t0, t1
integer(I4P)                    :: i, j, k, m, u, status
character(120)                  :: caption

allocate(f(4 * NS)) ; m = 0
do k = 0, NS - 1
   t0 = 2 * PI * k / NS ; t1 = 2 * PI * mod(k + 1, NS) / NS
   m = m + 1 ; call tri(f(m), [R * cos(t0), R * sin(t0), 0._R8P], [R * cos(t1), R * sin(t1), 0._R8P], &
                        [R * cos(t1), R * sin(t1), 1._R8P])
   m = m + 1 ; call tri(f(m), [R * cos(t0), R * sin(t0), 0._R8P], [R * cos(t1), R * sin(t1), 1._R8P], &
                        [R * cos(t0), R * sin(t0), 1._R8P])
   m = m + 1 ; call tri(f(m), [0._R8P, 0._R8P, 1._R8P], [R * cos(t0), R * sin(t0), 1._R8P], [R * cos(t1), R * sin(t1), 1._R8P])
   m = m + 1 ; call tri(f(m), [0._R8P, 0._R8P, 0._R8P], [R * cos(t1), R * sin(t1), 0._R8P], [R * cos(t0), R * sin(t0), 0._R8P])
enddo
call cyl%adopt_facets(facets=f) ; call cyl%analyze
call cyl%feature_edges(edges=edges)
call cyl%build_features(status=status)
if (status /= STATUS_OK .or. .not. cyl%is_watertight()) error stop 'cylinder: build_features failed'
print '(A,I0,A,I0,A,F8.5)', 'cylinder: ', cyl%get_facets_number(), ' facets, feature edges ', size(edges, 2), &
                            ', volume ', cyl%get_volume()
do j = 1, N
   do i = 1, N
      p(i + N * (j - 1)) = vector_R8P(-1._R8P + H * (i - 0.5_R8P), 0._R8P, -0.5_R8P + H * (j - 0.5_R8P))
   enddo
enddo
call cyl%query_many(points=p, distances=d, feature_distances=fd, status=status)
if (status /= STATUS_OK) error stop 'query_many failed'
print '(A,F6.4,A,F6.4)', 'plane y = 0: min |phi| ', minval(abs(d)), ', min feature distance ', minval(fd)
print '(A,F6.4,A,F6.4)', 'at the side (x = 0.479, z = 0.479): |phi| ', abs(d(36 + N * 23)), &
                         ', feature distance ', fd(36 + N * 23)
write(caption, '(A,I0,A)') 'y = 0: |phi| (left), distance to the ', size(edges, 2), ' rim edges (right); 0 cream, 0.5 sepia'
call fig_open(u, 'ref-ib-queries-feature', trim(caption), -90._R8P, 0._R8P, width=640_I4P, height=560_I4P)
do k = 1, N * N
   x = p(k)%x ; z = p(k)%z
   call fig_quad(u, 'frac', vector_R8P(x - H / 2, 0._R8P, z - H / 2), vector_R8P(x + H / 2, 0._R8P, z - H / 2), &
                 vector_R8P(x + H / 2, 0._R8P, z + H / 2), vector_R8P(x - H / 2, 0._R8P, z + H / 2), &
                 value=min(1._R8P, merge(abs(d(k)), fd(k), x < 0._R8P) / SCALE))
enddo
a(1:4) = [vector_R8P(-R, -0.05_R8P, 0._R8P), vector_R8P(R, -0.05_R8P, 0._R8P), vector_R8P(-R, -0.05_R8P, 0._R8P), &
          vector_R8P(R, -0.05_R8P, 1._R8P)]
b(1:4) = [vector_R8P(-R, -0.05_R8P, 1._R8P), vector_R8P(R, -0.05_R8P, 1._R8P), vector_R8P(R, -0.05_R8P, 0._R8P), &
          vector_R8P(-R, -0.05_R8P, 1._R8P)]
call fig_segments(u, 'accent', a(1:4), b(1:4))  ! the section of the cylinder by the plane
call fig_points(u, 'hole', a(1:4))               ! where the plane crosses the rims
call fig_label(u, 'accent', vector_R8P(-0.5_R8P, -0.06_R8P, 1.38_R8P), '|phi|')
call fig_label(u, 'accent', vector_R8P(0.5_R8P, -0.06_R8P, 1.38_R8P), 'feature distance')
call fig_close(u)
contains
   subroutine tri(ft, p1, p2, p3)
   !< A facet p1, p2, p3 (counter-clockwise seen from outside).
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call ft%compute_metrix
   endsubroutine tri
endprogram fig_ref_ib_queries_feature
