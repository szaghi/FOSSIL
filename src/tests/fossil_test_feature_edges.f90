!< FOSSIL test: feature edges, corners and polylines (issue #61, roadmap #58 WP3, step 1).
!<
!< 1. cube-01: every feature edge is sharp and lies on one of the cube's 12 edges, their total length is 12, the corners
!<    are the 8 cube corners, the polylines are 12 (one per cube edge, corner to corner).
!< 2. A cylinder (regular 32-gon extruded, caps fanned): the two cap rims are sharp (90 degrees), the side edges are not
!<    (11.25 degrees between neighbours), no corners: two closed polylines of 32 edges.
!< 3. cube-01 without one facet: 3 boundary edges (the hole), and with include_boundary=.false. none.
!< 4. Three triangles sharing one edge: that edge is non-manifold, the other 6 are boundary.
!< 5. dragon.stl: the counts are unchanged by a translation and decrease (not strictly) as the angle grows.
program fossil_test_feature_edges

use fossil, only : surface_stl_object, facet_object, FEATURE_EDGE_SHARP, FEATURE_EDGE_BOUNDARY, &
                   FEATURE_EDGE_NONMANIFOLD, STATUS_OK
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
type(surface_stl_object)        :: s, s2
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp(:)
type(vertex_pool_object), pointer :: pool
type(vector_R8P)                :: a, b
integer(I4P), allocatable       :: edges(:,:), kinds(:), corners(:), ps(:), pv(:)
real(R8P)                       :: total, angles(5)
integer(I4P)                    :: status, i, n(5), n0
logical                         :: ok, all_ok, on_edge

all_ok = .true.

! ---- 1. cube-01
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call s%feature_edges(edges=edges, kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
pool => s%get_vertex_pool()
total = 0._R8P ; on_edge = .true.
do i = 1, size(kinds)
   a = pool%coord(edges(1, i)) ; b = pool%coord(edges(2, i))
   a = b - a ; total = total + a%normL2() ; a = pool%coord(edges(1, i))
   on_edge = on_edge .and. count([on_face(a%x) .and. on_face(b%x), on_face(a%y) .and. on_face(b%y), &
                                  on_face(a%z) .and. on_face(b%z)]) >= 2
enddo
ok = status == STATUS_OK .and. all(kinds == FEATURE_EDGE_SHARP) .and. abs(total - 12._R8P) <= 1e-12_R8P .and. on_edge .and. &
     size(corners) == 8 .and. size(ps) - 1 == 12
do i = 1, size(corners)
   a = pool%coord(corners(i))
   ok = ok .and. on_face(a%x) .and. on_face(a%y) .and. on_face(a%z)
enddo
print '(A,I0,A,F8.4,A,I0,A,I0,A,L1)', '1. cube-01: feature edges ', size(kinds), ', length ', total, ', corners ', &
                                      size(corners), ', polylines ', size(ps) - 1, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 2. cylinder
call make_cylinder(s, 32_I4P)
call s%feature_edges(edges=edges, kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
ok = status == STATUS_OK .and. size(kinds) == 64 .and. all(kinds == FEATURE_EDGE_SHARP) .and. size(corners) == 0 .and. &
     size(ps) - 1 == 2
if (ok) ok = ps(2) - ps(1) == 33 .and. ps(3) - ps(2) == 33 .and. pv(ps(1)) == pv(ps(2) - 1)
print '(A,I0,A,I0,A,I0,A,L1)', '2. cylinder: feature edges ', size(kinds), ', corners ', size(corners), ', polylines ', &
                               size(ps) - 1, '  ok ', ok
all_ok = all_ok .and. ok

! ---- 3. open cube
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
fp => s%facets_ref()
allocate(f(size(fp) - 1))
do i = 1, size(f)
   f(i) = fp(i + 1) ; f(i)%vertex_id = 0_I4P
enddo
call s2%adopt_facets(facets=f)
call s2%feature_edges(kinds=kinds, status=status)
ok = status == STATUS_OK .and. count(kinds == FEATURE_EDGE_BOUNDARY) == 3
call s2%feature_edges(kinds=kinds, include_boundary=.false., status=status)
ok = ok .and. count(kinds == FEATURE_EDGE_BOUNDARY) == 0
print '(A,L1)', '3. open cube: 3 boundary edges, none when excluded: ', ok
all_ok = all_ok .and. ok

! ---- 4. three triangles on one edge
if (allocated(f)) deallocate(f)   ! adopt_facets moved it
allocate(f(3))
do i = 1, 3
   f(i)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; f(i)%vertex(2) = vector_R8P(1._R8P, 0._R8P, 0._R8P)
   f(i)%vertex(3) = vector_R8P(0.5_R8P, cos(2 * PI * i / 3), sin(2 * PI * i / 3))
   call f(i)%compute_metrix
enddo
call s%destroy
call s%adopt_facets(facets=f)
call s%feature_edges(kinds=kinds, status=status)
ok = status == STATUS_OK .and. count(kinds == FEATURE_EDGE_NONMANIFOLD) == 1 .and. count(kinds == FEATURE_EDGE_BOUNDARY) == 6
print '(A,L1)', '4. three sheets on one edge: 1 non-manifold, 6 boundary: ', ok
all_ok = all_ok .and. ok

! ---- 5. dragon: translation invariance, monotone in the angle
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call s%sanitize
call s%feature_edges(kinds=kinds, status=status)
n0 = size(kinds)
s2 = s
call s2%translate(delta=1._R8P * ex_R8P + 2._R8P * ey_R8P + 3._R8P * ez_R8P)
call s2%analyze
call s2%feature_edges(kinds=kinds, status=status)
ok = status == STATUS_OK .and. size(kinds) == n0
angles = [10._R8P, 20._R8P, 30._R8P, 45._R8P, 60._R8P]
do i = 1, 5
   call s%feature_edges(kinds=kinds, angle=angles(i))
   n(i) = size(kinds)
enddo
ok = ok .and. all(n(2:5) <= n(1:4)) .and. n(3) == n0
print '(A,I0,A,5I6,A,L1)', '5. dragon: ', n0, ' feature edges; for 10..60 degrees', n, '  ok ', ok
all_ok = all_ok .and. ok

print '(A,L1)', 'Are all tests passed? ', all_ok

contains
   pure function on_face(x) result(yes)
   !< x on a face of the unit cube [0, 1] (to 1e-12).
   real(R8P), intent(in) :: x
   logical               :: yes

   yes = abs(x) <= 1e-12_R8P .or. abs(x - 1._R8P) <= 1e-12_R8P
   endfunction on_face

   subroutine make_cylinder(surf, m)
   !< A regular m-gon of unit radius extruded over z in [0, 1], caps fanned from vertex 1, outward.
   type(surface_stl_object), intent(inout) :: surf
   integer(I4P),             intent(in)    :: m
   type(facet_object), allocatable         :: ff(:)
   real(R8P)                               :: u(m), v(m)
   integer(I4P)                            :: k, j, c

   u = [(cos(2 * PI * (k - 1) / m), k = 1, m)] ; v = [(sin(2 * PI * (k - 1) / m), k = 1, m)]
   allocate(ff(2 * (m - 2) + 2 * m))
   c = 0
   do k = 2, m - 1
      c = c + 1 ; call tri(ff(c), [u(1), v(1), 1._R8P], [u(k), v(k), 1._R8P], [u(k+1), v(k+1), 1._R8P])
      c = c + 1 ; call tri(ff(c), [u(1), v(1), 0._R8P], [u(k+1), v(k+1), 0._R8P], [u(k), v(k), 0._R8P])
   enddo
   do k = 1, m
      j = mod(k, m) + 1
      c = c + 1 ; call tri(ff(c), [u(k), v(k), 0._R8P], [u(j), v(j), 0._R8P], [u(j), v(j), 1._R8P])
      c = c + 1 ; call tri(ff(c), [u(k), v(k), 0._R8P], [u(j), v(j), 1._R8P], [u(k), v(k), 1._R8P])
   enddo
   call surf%destroy
   call surf%adopt_facets(facets=ff)
   endsubroutine make_cylinder

   subroutine tri(ft, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call ft%compute_metrix
   endsubroutine tri
endprogram fossil_test_feature_edges
