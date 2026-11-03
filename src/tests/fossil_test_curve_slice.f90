!< FOSSIL test: slicing a surface by an axis-aligned plane (issue #36 step 5).
program fossil_test_curve_slice
!< FOSSIL test: slicing a surface by an axis-aligned plane (issue #36 step 5).
!<
!< cube-02 (unit cube, 4 x 4 grid squares per face):
!< - z = 0.4 (generic): one closed loop, valid, area 1 to 1e-15; with collinear vertices merged,
!<   exactly the four corners of the unit square, area exactly 1;
!< - z = 0.5 and z = 0.25: the plane passes through grid vertices and along facet edges; merged,
!<   the same four corners;
!< - z = 1 (the top face lies in the plane): the outline of the top face, area 1;
!<   z = 0 (the bottom face lies in the plane, its vertices count as above): empty;
!< - x = 0.5 and y = 0.75: unit squares in the (y, z) and (z, x) planes, valid.
!< Icosphere (levels 2, 3, 4) at z = 0.3: valid, no self-intersection, area converging to
!< pi (1 - 0.09) by at least 3x per level.
!< cube-02 without one facet crossing z = 0.4: one open chain, and the fast winding number still
!< puts the section's centre inside.
!< Icosphere (level 3) with every third facet removed: many open chains (an open chain adds its end
!< point: the vertex buffer must hold up to twice the segments), every vertex kept.

use fossil, only : surface_stl_object, facet_object, curve_object, slice_surface, CURVE_VALID, SIGN_FAST_WINDING_NUMBER
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
type(surface_stl_object)        :: cube, holed, sphere
type(curve_object)              :: c
type(facet_object), allocatable :: f(:)
type(facet_object), pointer     :: fp
integer(I4P), allocatable       :: pairs(:,:)
real(R8P)                       :: err(3)
integer(I4P)                    :: i, n, nopen, lev
logical                         :: ok(7)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)

! generic level
call slice_surface(cube, 3_I4P, 0.4_R8P, c)
ok(1) = c%get_loops_number() == 1 .and. c%check_validity() == CURVE_VALID .and. abs(c%get_area() - 1._R8P) < 1.e-15_R8P
call slice_surface(cube, 3_I4P, 0.4_R8P, c, merge_collinear=.true.)
ok(1) = ok(1) .and. unit_square(c)
print '(A,L1)', 'cube z = 0.4: ', ok(1)

! through vertices and along edges
call slice_surface(cube, 3_I4P, 0.5_R8P, c, merge_collinear=.true.)
ok(2) = unit_square(c)
call slice_surface(cube, 3_I4P, 0.25_R8P, c, merge_collinear=.true.)
ok(2) = ok(2) .and. unit_square(c)
print '(A,L1)', 'cube z = 0.5, 0.25 (vertices, edges in the plane): ', ok(2)

! faces in the plane
call slice_surface(cube, 3_I4P, 1._R8P, c, merge_collinear=.true.)
ok(3) = unit_square(c)
call slice_surface(cube, 3_I4P, 0._R8P, c)
ok(3) = ok(3) .and. c%get_loops_number() == 0
print '(A,L1)', 'cube z = 1 (top face), z = 0 (bottom face, empty): ', ok(3)

! other axes
call slice_surface(cube, 1_I4P, 0.5_R8P, c, merge_collinear=.true.)
ok(4) = unit_square(c) .and. c%get_normal_axis() == 1
call slice_surface(cube, 2_I4P, 0.75_R8P, c, merge_collinear=.true.)
ok(4) = ok(4) .and. unit_square(c) .and. c%get_normal_axis() == 2
print '(A,L1)', 'cube x = 0.5, y = 0.75: ', ok(4)

! icosphere convergence
ok(5) = .true.
do lev = 2, 4
   call icosphere(lev, f)
   call sphere%adopt_facets(facets=f)
   call slice_surface(sphere, 3_I4P, 0.3_R8P, c)
   call c%self_intersections(pairs, n)
   err(lev - 1) = abs(c%get_area() - PI * (1._R8P - 0.09_R8P))
   ok(5) = ok(5) .and. c%get_loops_number() == 1 .and. c%check_validity() == CURVE_VALID .and. n == 0
enddo
ok(5) = ok(5) .and. err(2) < err(1) / 3._R8P .and. err(3) < err(2) / 3._R8P
print '(A,L1,3ES11.3)', 'icosphere z = 0.3, area errors (levels 2, 3, 4): ', ok(5), err

! open section
allocate(f(cube%get_facets_number())) ; n = 0 ; nopen = 0
do i = 1, cube%get_facets_number()
   fp => cube%facet_at(i)
   if (nopen == 0 .and. fp%centroid%x > 1._R8P - 1.e-9_R8P .and. abs(fp%centroid%z - 0.4_R8P) < 0.1_R8P) then
      nopen = 1 ; cycle ! drop one facet of the face x = 1 crossing z = 0.4
   endif
   n = n + 1 ; f(n) = fp ; f(n)%vertex_id = 0_I4P
enddo
f = f(1:n)
call holed%adopt_facets(facets=f)
call slice_surface(holed, 3_I4P, 0.4_R8P, c, open_loops=nopen)
ok(6) = nopen == 1 .and. c%get_loops_number() == 1 .and. .not. c%is_loop_closed(1_I4P) .and. &
        c%is_point_inside([0.5_R8P, 0.5_R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER) .and. &
        .not. c%is_point_inside([1.5_R8P, 0.5_R8P], sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,L1)', 'open section (one facet removed): ', ok(6)

! many open chains
call icosphere(3_I4P, f)
f = pack(f, [(mod(i, 3) /= 0, i = 1, size(f))])
call holed%adopt_facets(facets=f)
call slice_surface(holed, 3_I4P, 0.3_R8P, c, open_loops=nopen)
ok(7) = nopen > 1 .and. nopen == c%get_loops_number()
print '(A,L1,A,I0)', 'many open chains: ', ok(7), '  chains ', nopen

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function unit_square(cur) result(yes)
   !< True if cur is exactly the counter-clockwise unit square (4 vertices, area 1).
   type(curve_object), intent(in) :: cur !< Curve.
   logical                        :: yes !< Result.
   real(R8P)                      :: corner(2, 4) !< Corners.
   integer(I4P)                   :: k, j !< Counters.
   logical                        :: found !< Corner found.

   corner = reshape([0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P], [2, 4])
   yes = cur%get_loops_number() == 1 .and. cur%get_vertices_number() == 4 .and. cur%get_area() == 1._R8P .and. &
         cur%check_validity() == CURVE_VALID
   if (.not. yes) return
   do k = 1, 4
      found = .false.
      do j = 1, 4
         found = found .or. all(cur%get_vertex(j) == corner(:, k))
      enddo
      yes = yes .and. found
   enddo
   endfunction unit_square

   subroutine icosphere(levels, fo)
   !< Unit icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere.
   integer(I4P),                    intent(in)  :: levels !< Subdivisions.
   type(facet_object), allocatable, intent(out) :: fo(:)  !< Facets.
   type(vector_R8P), allocatable :: tri(:,:), nxt(:,:)
   type(vector_R8P)              :: v(12), m12, m23, m31
   integer(I4P), parameter       :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                        2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                        4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                        5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   real(R8P)                     :: g
   integer(I4P)                  :: k, l, m

   g = (1._R8P + sqrt(5._R8P)) / 2._R8P
   v = [vec(-1._R8P, g, 0._R8P), vec(1._R8P, g, 0._R8P), vec(-1._R8P, -g, 0._R8P), vec(1._R8P, -g, 0._R8P), &
        vec(0._R8P, -1._R8P, g), vec(0._R8P, 1._R8P, g), vec(0._R8P, -1._R8P, -g), vec(0._R8P, 1._R8P, -g), &
        vec(g, 0._R8P, -1._R8P), vec(g, 0._R8P, 1._R8P), vec(-g, 0._R8P, -1._R8P), vec(-g, 0._R8P, 1._R8P)]
   do k = 1, 12
      v(k) = unit(v(k))
   enddo
   allocate(tri(3, 20))
   do k = 1, 20
      tri(:, k) = v(IDX(:, k))
   enddo
   do l = 1, levels
      m = size(tri, dim=2)
      allocate(nxt(3, 4 * m))
      do k = 1, m
         m12 = unit(tri(1, k) + tri(2, k)) ; m23 = unit(tri(2, k) + tri(3, k)) ; m31 = unit(tri(3, k) + tri(1, k))
         nxt(:, 4 * k - 3) = [tri(1, k), m12, m31]
         nxt(:, 4 * k - 2) = [m12, tri(2, k), m23]
         nxt(:, 4 * k - 1) = [m31, m23, tri(3, k)]
         nxt(:, 4 * k)     = [m12, m23, m31]
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   m = size(tri, dim=2)
   allocate(fo(m))
   do k = 1, m
      fo(k)%vertex = tri(:, k)
      call fo(k)%compute_metrix
   enddo
   endsubroutine icosphere

   pure function vec(x, y, z) result(r)
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: r

   r%x = x ; r%y = y ; r%z = z
   endfunction vec

   pure function unit(a) result(r)
   type(vector_R8P), intent(in) :: a
   type(vector_R8P)             :: r
   real(R8P)                    :: nn

   nn = sqrt(a%x**2 + a%y**2 + a%z**2)
   r%x = a%x / nn ; r%y = a%y / nn ; r%z = a%z / nn
   endfunction unit
endprogram fossil_test_curve_slice
