!< Animation of the reference page "Thickness and gaps": the shrinking ball of the medial method at three points of a
!< NACA 0012 wing (60 points per side, extruded), seen along the span. Each ball starts at half the axial chord, then
!< shrinks to the ball through p and the closest surface point q to its centre, tangent at p, until q lies on it; the
!< frames show the iterations (ember circle, gold q). The program replays the iteration of `thickness` with the same
!< steps (`intersect_ray_first(exclude=)`, `query_many`, `ball_radius`) and checks that its final diameters are those of
!< `thickness(method=THICKNESS_MEDIAL)`.
program anim_shrinking_ball
use fossil,              only : surface_stl_object, STATUS_OK, THICKNESS_MEDIAL
use fossil_facet_object, only : facet_object
use fossil_ray_query,    only : ray_hit_t
use fossil_thickness,    only : ball_radius
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_segments, fig_points, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
real(R8P),    parameter         :: PI = acos(-1._R8P)
integer(I4P), parameter         :: NP = 60, NS = 3, KMAX = 12
real(R8P),    parameter         :: XS(NS) = [0.3_R8P, 0.7_R8P, 0.93_R8P], SIDE(NS) = [1._R8P, -1._R8P, 1._R8P]
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
type(ray_hit_t)                 :: hit
type(vector_R8P)                :: c(1), dd
type(vector_R8P), allocatable   :: a(:), b(:)
real(R8P), allocatable          :: th(:), u(:), v(:)
real(R8P)                       :: r(0:KMAX, NS), q(3, 0:KMAX, NS), p(3, NS), n(3, NS), dq(1), qq(3, 1), rn, w(3), best
integer(I4P)                    :: id(NS), last(NS), uu, i, k, j, status, nv
logical                         :: has, ok
character(96)                   :: caption

call wing(f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, method=THICKNESS_MEDIAL, status=status)
if (status /= STATUS_OK) error stop 'thickness failed'
do j = 1, NS                                           ! the side facet nearest to (XS(j), +-y_t, 0.5)
   best = huge(1._R8P)
   do i = 1, s%get_facets_number()
      fp => s%facet_at(i)
      if (abs(fp%normal%z) > 0.5_R8P .or. fp%normal%y * SIDE(j) < 0._R8P) cycle
      if (abs(fp%centroid%x - XS(j)) + abs(fp%centroid%z - 0.5_R8P) < best) then
         best = abs(fp%centroid%x - XS(j)) + abs(fp%centroid%z - 0.5_R8P) ; id(j) = i
      endif
   enddo
   fp => s%facet_at(id(j))
   p(:, j) = [fp%centroid%x, fp%centroid%y, fp%centroid%z] ; n(:, j) = -[fp%normal%x, fp%normal%y, fp%normal%z]
enddo
ok = .true.
do j = 1, NS                                           ! replay the iteration
   dd = vector_R8P(n(1, j), n(2, j), n(3, j)) ; fp => s%facet_at(id(j))
   call s%intersect_ray_first(ray_origin=fp%centroid, ray_direction=dd, hit=hit, has_hit=has, exclude=id(j))
   r(0, j) = 0.5_R8P * hit%t ; w = [hit%point%x, hit%point%y, hit%point%z]
   q(:, 0, j) = merge(0._R8P, w, abs(w) < 1e-30_R8P) ; last(j) = 0
   do k = 1, KMAX
      c(1) = vector_R8P(p(1, j) + r(k - 1, j) * n(1, j), p(2, j) + r(k - 1, j) * n(2, j), p(3, j) + r(k - 1, j) * n(3, j))
      call s%query_many(points=c, distances=dq, closest_points=qq, is_signed=.false.)
      r(k, j) = r(k - 1, j) ; q(:, k, j) = merge(0._R8P, qq(:, 1), abs(qq(:, 1)) < 1e-30_R8P)   ! no denormals in the scene
      if (dq(1) >= r(k - 1, j) * (1._R8P - 1e-12_R8P)) exit
      rn = ball_radius(p(:, j), n(:, j), qq(:, 1))
      if (rn >= r(k - 1, j)) exit
      w = qq(:, 1) - p(:, j) - rn * n(:, j)
      if (-dot_product(n(:, j), w) > cos(PI / 6._R8P) * rn) exit
      r(k, j) = rn ; last(j) = k
   enddo
   r(last(j) + 1:, j) = r(last(j), j)
   do i = min(k, KMAX) + 1, KMAX                        ! the frames after convergence repeat the last contact
      q(:, i, j) = q(:, i - 1, j)
   enddo
   ok = ok .and. 2._R8P * r(last(j), j) == th(id(j))
   print '(A,F5.2,A,I0,A,F8.5,A,F8.5)', 'x = ', XS(j), ': ', last(j), ' shrinks, diameter ', 2._R8P * r(last(j), j), &
                                        ', thickness() ', th(id(j))
enddo
if (.not. ok) error stop 'replay differs from thickness()'
call fig_open(uu, 'anim-shrinking-ball', '', -90._R8P, 90._R8P, width=900_I4P, height=360_I4P)
call fig_anim(uu, 1.2_R8P)
do k = 0, maxval(last)
   write(caption, '(A,I0,A,3F7.4)') 'iteration ', k, ': diameters ', 2._R8P * r(k, :)
   call fig_frame(uu, caption=trim(caption), seconds=merge(3._R8P, 1.2_R8P, k == maxval(last)))
   call fig_segments(uu, 'body', a, b)
   do j = 1, NS
      call circle(p(:, j) + r(k, j) * n(:, j), r(k, j))
      call fig_points(uu, 'accent', [vector_R8P(q(1, k, j), q(2, k, j), 0._R8P)])
      call fig_points(uu, 'hole', [vector_R8P(p(1, j), p(2, j), 0._R8P)])
   enddo
enddo
call fig_close(uu)
contains
   subroutine wing(ff)
   !< NACA 0012 of unit chord, NP points per side (cosine spacing, closed trailing edge), extruded over z in [0, 1]; the
   !< section outline kept in (a, b) for drawing.
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: xs(0:NP), ys(0:NP)
   integer(I4P)                                 :: ii, m, i2

   do ii = 0, NP
      xs(ii) = 0.5_R8P * (1._R8P - cos(PI * ii / NP))
      ys(ii) = 0.6_R8P * (0.2969_R8P * sqrt(xs(ii)) - 0.1260_R8P * xs(ii) - 0.3516_R8P * xs(ii) ** 2 + &
                          0.2843_R8P * xs(ii) ** 3 - 0.1036_R8P * xs(ii) ** 4)
   enddo
   ys(NP) = 0._R8P
   allocate(u(2 * NP), v(2 * NP)) ; nv = 0
   do ii = NP, 1, -1
      nv = nv + 1 ; u(nv) = xs(ii) ; v(nv) = ys(ii)
   enddo
   do ii = 0, NP - 1
      nv = nv + 1 ; u(nv) = xs(ii) ; v(nv) = -ys(ii)
   enddo
   allocate(ff(4 * nv), a(nv), b(nv))
   m = 0
   do ii = 1, nv
      i2 = mod(ii, nv) + 1
      a(ii) = vector_R8P(u(ii), v(ii), 0._R8P) ; b(ii) = vector_R8P(u(i2), v(i2), 0._R8P)
      m = m + 1 ; call tri(ff(m), [u(ii), v(ii), 0._R8P], [u(i2), v(i2), 0._R8P], [u(i2), v(i2), 1._R8P])
      m = m + 1 ; call tri(ff(m), [u(ii), v(ii), 0._R8P], [u(i2), v(i2), 1._R8P], [u(ii), v(ii), 1._R8P])
      m = m + 1 ; call tri(ff(m), [0.5_R8P, 0._R8P, 0._R8P], [u(ii), v(ii), 0._R8P], [u(i2), v(i2), 0._R8P])
      m = m + 1 ; call tri(ff(m), [0.5_R8P, 0._R8P, 1._R8P], [u(ii), v(ii), 1._R8P], [u(i2), v(i2), 1._R8P])
   enddo
   endsubroutine wing

   subroutine tri(ft, p1, p2, p3)
   !< A facet oriented away from the wing's axis (0.5, 0, 0.5) (the section is convex).
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), nn(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   nn = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(nn, (p1 + p2 + p3) / 3._R8P - [0.5_R8P, 0._R8P, 0.5_R8P]) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri

   subroutine circle(cc, rr)
   !< The circle of centre cc and radius rr in the section plane, 72 segments.
   real(R8P), intent(in)  :: cc(3), rr
   type(vector_R8P)       :: ca(72), cb(72)
   integer(I4P)           :: l

   do l = 1, 72
      ca(l) = vector_R8P(cc(1) + rr * cos(2 * PI * (l - 1) / 72), cc(2) + rr * sin(2 * PI * (l - 1) / 72), 0._R8P)
      cb(l) = vector_R8P(cc(1) + rr * cos(2 * PI * l / 72), cc(2) + rr * sin(2 * PI * l / 72), 0._R8P)
   enddo
   call fig_segments(uu, 'hole', ca, cb)
   endsubroutine circle
endprogram anim_shrinking_ball
