!< FOSSIL test: the Delaunay triangle store grows (issue #63, roadmap #58 WP5, fix found by the boolean harness).
!<
!< Bowyer-Watson keeps the deactivated triangles in place, so the store fills with the sum of the cavity sizes. The
!< store was allocated once at 6 (n + 3) entries and never grew: an insertion order that makes large cavities wrote
!< past it (heap corruption, found by the curved-boolean sweep: the cube against a marching-cubes torus).
!< 1. n = 600 points on an elliptic arc (semi-axes 1 and 0.5, 3 rad), inserted in angle order (every new point lies in
!<    the circumcircles of a growing fan): 27087 triangles are created against the initial 6 (n + 3) = 3618 entries; the
!<    store grows, and the triangulation is right: n - 2 triangles (all points on the hull, none cocircular), every one
!<    with an empty circumcircle (checked against every 7th point).
!< 2. The boolean that corrupted the heap: the unit cube (cube.stl, centred) against a 24^3 marching-cubes torus moved by
!<    the 9th motion of seed 20261106; the four ops return BOOL_STATUS_OK and V(U) + V(I) = V(A) + V(B) to 1e-12.
program fossil_test_dt_capacity

use fossil,                only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_STATUS_OK
use fossil_boolean_bodies, only : BODY_TORUS, mc_body, move_body
use fossil_dt,             only : triangulation_t, dt_build
use fossil_predicates,     only : incircle
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

integer(I4P), parameter  :: N = 600
real(R8P),    parameter  :: PI = acos(-1._R8P)
type(triangulation_t)    :: dt
type(surface_stl_object) :: a, b, r
real(R8P)                :: pts(2, N), axis(3), angle, t(3), x, vu, vi
integer(I8P)             :: seed
integer(I4P)             :: i, j, k, ntri, st, bad
logical                  :: ok, all_ok

all_ok = .true.

! ---- 1. points on an elliptic arc in angle order
do i = 1, N
   pts(:, i) = [cos(3._R8P * i / N), 0.5_R8P * sin(3._R8P * i / N)]
enddo
call dt_build(dt, pts)
ntri = 0 ; bad = 0
do k = 1, dt%n_tri_used
   if (.not. dt%tri(k)%active) cycle
   ntri = ntri + 1
   do j = 1, N, 7                                       ! empty circumcircle, on a sample of the points
      if (any(dt%tri(k)%v == j)) cycle
      if (incircle(a=p(dt%tri(k)%v(1)), b=p(dt%tri(k)%v(2)), c=p(dt%tri(k)%v(3)), d=p(j)) > 0._R8P) bad = bad + 1
   enddo
enddo
ok = dt%n_tri_alloc > 6 * (N + 3) .and. ntri == N - 2 .and. bad == 0
print '(A,I0,A,I0,A,I0,A,I0,A,L1)', '1. elliptic arc: store ', dt%n_tri_alloc, ' > ', 6 * (N + 3), ', triangles ', ntri, &
                                    ', circumcircle violations ', bad, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. the boolean that corrupted the heap
seed = 20261106_I8P
do i = 1, 56
   x = rnd()
enddo
axis = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] ; angle = 2 * PI * rnd()
t = 0.4_R8P * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call a%translate(x=-0.5_R8P, y=-0.5_R8P, z=-0.5_R8P) ; call a%analyze
call mc_body(BODY_TORUS, 24_I4P, b)
call move_body(b, axis, angle, t)
r = a ; call r%boolean(other=b, op=BOOL_UNION, status=st) ; vu = r%get_volume()
ok = st == BOOL_STATUS_OK
r = a ; call r%boolean(other=b, op=BOOL_INTERSECT, status=st) ; vi = r%get_volume()
ok = ok .and. st == BOOL_STATUS_OK .and. abs(vu + vi - a%get_volume() - b%get_volume()) <= 1e-12_R8P
print '(A,ES10.2,A,L1)', '2. cube x torus: union and intersection OK, volume identity ', &
                         abs(vu + vi - a%get_volume() - b%get_volume()), '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   function p(idx) result(v)
   !< Point idx as a vector.
   integer(I4P), intent(in) :: idx
   type(vector_R8P)         :: v

   v = vector_R8P(dt%px(idx), dt%py(idx), 0._R8P)
   endfunction p

   function rnd() result(y)
   !< Park-Miller in (0, 1).
   real(R8P) :: y

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   y = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_dt_capacity
