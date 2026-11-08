!< Figures of the reference page "Body-body distance", What it does, Method and Accuracy:
!< - the dragon and a cube beside its head: distance_to's witness points joined by a segment, the two facets that
!<   realize the distance (facet_self, facet_other) in ember;
!< - the pair distance of two disjoint triangles: the 6 vertex-triangle and 9 edge-edge candidates (closest_on_triangle,
!<   closest_segment_segment) drawn faint, the least of them, equal to tri_tri_distance, in gold;
!< - the error of the segment-segment distance against a quadruple-precision reference, in units of eps L, against the
!<   angle between the segments (1e-2 .. 1e-15; the largest of 40 fixed-seed pairs per decade): the library's kernel
!<   (teal) stays under 1, the textbook solve of Ericson 5.1.9 (ember, written here) is off by up to 1e7 where the
!<   angle is small but a e - b^2 still nonzero.
program fig_ref_body_distance_method
use, intrinsic :: iso_fortran_env, only : real128
use fossil,               only : surface_stl_object
use fossil_body_distance, only : tri_tri_distance, closest_on_triangle, closest_segment_segment
use fossil_fig,           only : fig_open, fig_zoom, fig_surface, fig_segment, fig_segments, fig_points, fig_label, &
                                 fig_triangle, fig_close
use penf,                 only : I4P, I8P, R8P
use vecfor,               only : vector_R8P

implicit none
integer,      parameter  :: QP = real128
integer(I4P), parameter  :: K0 = 2, K1 = 15, NPAIR = 40
type(surface_stl_object) :: dragon, cube
type(vector_R8P)         :: bmin, bmax, ca(15), cb(15), pts(K1 - K0 + 1), ptt(K1 - K0 + 1)
real(R8P)                :: d, pa(3), pb(3), ta(3, 3), tb(3, 3), d2, x(3), y(3), cand(15), size_d
real(R8P)                :: p1(3), q1(3), p2(3), q2(3), c1(3), c2(3), u3(3), w3(3), e3(3), ang, h, l, dref
real(R8P)                :: elib(K0:K1), etxt(K0:K1)
integer(I4P)             :: fa, fb, i, j, k, m, n, u
integer(I8P)             :: seed
logical                  :: meet, hit
character(96)            :: caption
character(16)            :: text

! ---- what it does: the dragon and a cube beside its head
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
size_d = bmax%x - bmin%x
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)          ! unit cube, side 0.12 of the dragon
call cube%resize(x=0.12_R8P * size_d, y=0.12_R8P * size_d, z=0.12_R8P * size_d)
call cube%translate(delta=vector_R8P(bmax%x + 0.06_R8P * size_d, 0.5_R8P * (bmin%y + bmax%y), &
                                     bmin%z + 0.55_R8P * (bmax%z - bmin%z)))
call cube%analyze
call dragon%distance_to(cube, distance=d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb, intersecting=hit)
print '(A,F8.4,A,I0,A,I0,A,L1)', 'distance_to: ', d, ', facet_self ', fa, ', facet_other ', fb, ', intersecting ', hit
write(caption, '(A,F6.3,A,I0,A,I0)') 'distance_to = ', d, ' between dragon facet ', fa, ' and cube facet ', fb
call fig_open(u, 'ref-body-distance-witness', trim(caption), -60._R8P, 15._R8P)
call fig_zoom(u, vector_R8P(0.5_R8P * (pa(1) + pb(1)), 0.5_R8P * (pa(2) + pb(2)), 0.5_R8P * (pa(3) + pb(3))), &
              0.16_R8P * size_d)
call fig_surface(u, 'body', dragon)
call fig_surface(u, 'ghost', cube)
call fig_surface(u, 'patch', dragon, first=fa, last=fa)
call fig_surface(u, 'patch', cube, first=fb, last=fb)
call fig_segment(u, 'accent', vector_R8P(pa(1), pa(2), pa(3)), vector_R8P(pb(1), pb(2), pb(3)))
call fig_points(u, 'foot', [vector_R8P(pa(1), pa(2), pa(3)), vector_R8P(pb(1), pb(2), pb(3))])
call fig_close(u)

! ---- method: the 15 candidates of a pair of disjoint triangles
ta(:, 1) = [0._R8P, 0._R8P, 0._R8P] ; ta(:, 2) = [1._R8P, 0._R8P, 0._R8P] ; ta(:, 3) = [0.3_R8P, 0.9_R8P, 0._R8P]
tb(:, 1) = [0.6_R8P, 0.2_R8P, 0.35_R8P] ; tb(:, 2) = [1.3_R8P, 0.8_R8P, 0.6_R8P] ; tb(:, 3) = [0.4_R8P, 1.1_R8P, 0.9_R8P]
call tri_tri_distance(ta, tb, d2, pa, pb, meet)
n = 0
do i = 1, 3                                                       ! vertex-triangle: 3 + 3
   n = n + 1 ; x = tb(:, i) ; y = closest_on_triangle(x, ta) ; call keep(y, x)
   n = n + 1 ; x = ta(:, i) ; y = closest_on_triangle(x, tb) ; call keep(x, y)
enddo
do i = 1, 3                                                       ! edge-edge: 3 x 3
   do j = 1, 3
      n = n + 1
      call closest_segment_segment(ta(:, i), ta(:, mod(i, 3) + 1), tb(:, j), tb(:, mod(j, 3) + 1), x, y)
      call keep(x, y)
   enddo
enddo
m = minloc(cand, dim=1)
print '(A,F8.5,A,F8.5,A,I0)', 'least of the 15 candidates ', cand(m), ', tri_tri_distance ', sqrt(d2), &
                              ', candidate ', m
write(caption, '(A,F6.4,A,F6.4)') '15 candidates; the least, ', cand(m), ' = tri_tri_distance ', sqrt(d2)
call fig_open(u, 'ref-body-distance-candidates', trim(caption), -110._R8P, 25._R8P)
call fig_triangle(u, 'body', v(ta(:, 1)), v(ta(:, 2)), v(ta(:, 3)))
call fig_triangle(u, 'probe', v(tb(:, 1)), v(tb(:, 2)), v(tb(:, 3)))
call fig_segments(u, 'box', pack(ca, [(i /= m, i = 1, 15)]), pack(cb, [(i /= m, i = 1, 15)]))
call fig_segment(u, 'accent', ca(m), cb(m))
call fig_points(u, 'foot', [ca(m), cb(m)])
call fig_close(u)

! ---- accuracy: segment-segment error against the angle
seed = 1069_I8P ; elib = 0._R8P ; etxt = 0._R8P
do k = K0, K1
   do i = 1, NPAIR
      u3 = unit3([rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P])
      w3 = unit3(cross3(u3, [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]))
      ang = 10._R8P**(-k - rnd()) ; h = 10._R8P**(-3._R8P - 9._R8P * rnd())
      p1 = [rnd(), rnd(), rnd()] ; q1 = p1 + (0.2_R8P + rnd()) * u3
      e3 = unit3(u3 + ang * w3)
      p2 = p1 + h * unit3(cross3(u3, w3)) + (rnd() - 0.5_R8P) * 1.5_R8P * u3
      q2 = p2 + (0.2_R8P + rnd()) * e3
      dref = real(segseg_ref(p1, q1, p2, q2), R8P)
      l = max(maxval(abs([p1, q1, p2, q2])), 1._R8P)
      call closest_segment_segment(p1, q1, p2, q2, c1, c2)
      elib(k) = max(elib(k), abs(norm2(c1 - c2) - dref) / (epsilon(1._R8P) * l))
      call textbook(p1, q1, p2, q2, c1, c2)
      etxt(k) = max(etxt(k), abs(norm2(c1 - c2) - dref) / (epsilon(1._R8P) * l))
   enddo
   print '(A,I2,A,F8.2,A,ES9.2)', 'angle 1e-', k, ': closest_segment_segment ', elib(k), ' eps L, textbook ', etxt(k)
enddo
do k = K0, K1                                                     ! x = decade, y = log10 of the error (floor 0.1)
   pts(k - K0 + 1) = vector_R8P(real(k, R8P), log10(max(elib(k), 0.1_R8P)), 0._R8P)
   ptt(k - K0 + 1) = vector_R8P(real(k, R8P), log10(max(etxt(k), 0.1_R8P)), 0._R8P)
enddo
write(caption, '(A,F4.2,A,ES8.1,A)') 'max error / (eps L): library ', maxval(elib), ' (teal), textbook ', &
                                     maxval(etxt), ' (ember)'
call fig_open(u, 'ref-body-distance-accuracy', trim(caption), -90._R8P, 90._R8P)
call fig_segment(u, 'box', vector_R8P(real(K0, R8P), -1._R8P, 0._R8P), vector_R8P(real(K1, R8P), -1._R8P, 0._R8P))
call fig_segment(u, 'box', vector_R8P(real(K0, R8P), -1._R8P, 0._R8P), vector_R8P(real(K0, R8P), 8._R8P, 0._R8P))
do j = -1, 7, 2                                                   ! y ticks: 1e-1 .. 1e7
   write(text, '(A,I0)') '1e', j
   call fig_label(u, 'box', vector_R8P(real(K0, R8P) - 0.9_R8P, real(j, R8P), 0._R8P), trim(text))
   call fig_segment(u, 'box', vector_R8P(real(K0, R8P) - 0.15_R8P, real(j, R8P), 0._R8P), &
                    vector_R8P(real(K0, R8P), real(j, R8P), 0._R8P))
enddo
do k = K0, K1, 3                                                  ! x ticks: angle 1e-2 .. 1e-14
   write(text, '(A,I0)') '1e-', k
   call fig_label(u, 'box', vector_R8P(real(k, R8P) - 0.4_R8P, -1.6_R8P, 0._R8P), trim(text))
enddo
call fig_label(u, 'box', vector_R8P(real(K1, R8P) - 0.5_R8P, -0.8_R8P, 0._R8P), 'angle')
call fig_label(u, 'box', vector_R8P(real(K0, R8P) + 0.2_R8P, 7.8_R8P, 0._R8P), '|d - d_ref| / (eps L)')
call fig_segment(u, 'accent', vector_R8P(real(K0, R8P), log10(16._R8P), 0._R8P), &
                 vector_R8P(real(K1, R8P), log10(16._R8P), 0._R8P))
call fig_label(u, 'accent', vector_R8P(real(K1, R8P) - 3._R8P, log10(16._R8P), 0._R8P), 'test bound 16 eps L')
call fig_segments(u, 'probe', pts(1:K1 - K0), pts(2:K1 - K0 + 1)) ; call fig_points(u, 'probe', pts)
call fig_segments(u, 'outside', ptt(1:K1 - K0), ptt(2:K1 - K0 + 1)) ; call fig_points(u, 'outside', ptt)
call fig_close(u)
contains
   subroutine keep(a, b)
   !< Candidate n: the pair a (on the first triangle), b (on the second).
   real(R8P), intent(in) :: a(3), b(3)

   ca(n) = v(a) ; cb(n) = v(b) ; cand(n) = norm2(a - b)
   endsubroutine keep

   function v(a) result(w)
   !< An array as a vector.
   real(R8P), intent(in) :: a(3)
   type(vector_R8P)      :: w

   w = vector_R8P(a(1), a(2), a(3))
   endfunction v

   subroutine textbook(a1, b1, a2, b2, s1, s2)
   !< Ericson's ClosestPtSegmentSegment (Real-Time Collision Detection, 5.1.9), the solve used before issue #69.
   real(R8P), intent(in)  :: a1(3), b1(3), a2(3), b2(3)
   real(R8P), intent(out) :: s1(3), s2(3)
   real(R8P)              :: d1(3), d2(3), r(3), a, b, c, e, f, den, s, t

   d1 = b1 - a1 ; d2 = b2 - a2 ; r = a1 - a2
   a = dot_product(d1, d1) ; e = dot_product(d2, d2) ; f = dot_product(d2, r)
   if (e <= 0._R8P) then
      t = 0._R8P ; s = min(1._R8P, max(0._R8P, -dot_product(d1, r) / a))
   else
      c = dot_product(d1, r) ; b = dot_product(d1, d2) ; den = a * e - b * b
      s = 0._R8P
      if (den /= 0._R8P) s = min(1._R8P, max(0._R8P, (b * f - c * e) / den))
      t = (b * s + f) / e
      if (t < 0._R8P) then
         t = 0._R8P ; s = min(1._R8P, max(0._R8P, -c / a))
      elseif (t > 1._R8P) then
         t = 1._R8P ; s = min(1._R8P, max(0._R8P, (b - c) / a))
      endif
   endif
   s1 = a1 + s * d1 ; s2 = a2 + t * d2
   endsubroutine textbook

   function segseg_ref(a1, b1, a2, b2) result(dd)
   !< Distance between segments a1-b1 and a2-b2 in quadruple precision: the four endpoint-segment distances, and the
   !< line-line distance when the closest points of the lines fall inside both segments.
   real(R8P), intent(in) :: a1(3), b1(3), a2(3), b2(3)
   real(QP)              :: dd, x1(3), y1(3), x2(3), y2(3), d1(3), d2(3), r(3), nq(3), nn, s, tt

   x1 = real(a1, QP) ; y1 = real(b1, QP) ; x2 = real(a2, QP) ; y2 = real(b2, QP)
   dd = min(ps(x1, x2, y2), ps(y1, x2, y2), ps(x2, x1, y1), ps(y2, x1, y1))
   d1 = y1 - x1 ; d2 = y2 - x2 ; r = x2 - x1
   nq = crq(d1, d2) ; nn = dot_product(nq, nq)
   if (nn > 0._QP) then
      s = dot_product(crq(r, d2), nq) / nn ; tt = dot_product(crq(r, d1), nq) / nn
      if (s > 0._QP .and. s < 1._QP .and. tt > 0._QP .and. tt < 1._QP) dd = min(dd, abs(dot_product(r, nq)) / sqrt(nn))
   endif
   endfunction segseg_ref

   function ps(xq, a, b) result(dd)
   !< Distance from xq to the segment a-b, quadruple precision.
   real(QP), intent(in) :: xq(3), a(3), b(3)
   real(QP)             :: dd, ab(3), s

   ab = b - a ; s = 0._QP
   if (dot_product(ab, ab) > 0._QP) s = min(1._QP, max(0._QP, dot_product(xq - a, ab) / dot_product(ab, ab)))
   dd = sqrt(dot_product(xq - a - s * ab, xq - a - s * ab))
   endfunction ps

   pure function crq(a, b) result(c)
   real(QP), intent(in) :: a(3), b(3)
   real(QP)             :: c(3)

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction crq

   pure function cross3(a, b) result(c)
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: c(3)

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   pure function unit3(a) result(b)
   real(R8P), intent(in) :: a(3)
   real(R8P)             :: b(3)

   b = a / norm2(a)
   endfunction unit3

   function rnd() result(r)
   !< Park-Miller minimal standard generator: fixed seed, the same numbers on every run.
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fig_ref_body_distance_method
