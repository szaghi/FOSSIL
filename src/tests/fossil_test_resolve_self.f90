!< FOSSIL test: resolve_self_intersections selects the self-union (issue #78).
!<
!< 1. Two unit cubes overlapping by a corner, merged into one surface (`merge_solids`): the result equals the union
!<    of the two shells as separate surfaces (`boolean(BOOL_UNION)`): same volume and area to 1e-12, watertight, one
!<    shell, no self-intersection left.
!< 2. The folded cube of docs/guide/advanced/self-intersection.md: the unit cube with its corner (1, 1, 1) pushed to
!<    (0.3, 0.35, -0.3), through its own bottom face. One closed shell, 10 crossing pairs. The region it winds
!<    around once is the cube minus the dent, cut at z = 0; the pocket below the bottom face is wound -1 and goes.
!<    Checked: status OK, watertight, one shell, no crossing pair, no vertex below z = 0, and the volume equals the
!<    input's signed volume plus the pocket's volume (the cone of apex (0.3, 0.35, -0.3) over its section at z = 0,
!<    computed from the vertices around the corner).
!< 3. A clean cube is returned unchanged (same facets, bit for bit).
!< 4. The crossing cubes of 1 with a facet of the first cube's face x = 0 removed (outside the second cube): not
!<    closed, so BOOL_STATUS_NOT_CLOSED and the surface unchanged.
!< 5. dragon.stl (1191 pairs): it is open, so BOOL_STATUS_NOT_CLOSED and the surface unchanged. Before issue #91 it
!<    returned BOOL_STATUS_CDT_FAILED: the cuts of facet 36 cross each other (three facets meet at a point inside it)
!<    and the arrangement had no triple-point construction.
!< 6. Three unit cubes crossing pairwise (issue #91), translated by (0.3, 0.3, 0.3) and (0.6, -0.15, 0.45), merged into
!<    one surface: the faces x = 1 of the first, y = 0.3 of the second and z = 0.45 of the third meet at
!<    (1, 0.3, 0.45), strictly inside a facet of each (a three-plane point, IPK_TPI). The result equals the boolean
!<    chain (A u B) u C: same volume and area to 1e-12, watertight, manifold, one shell, no self-intersection left.
!<    Base 204cd62: BOOL_STATUS_CDT_FAILED.
!< Negative control (base 9335a01): case 1 returns status 0, 224 facets, volume 2 (the arrangement unselected), not
!< watertight; case 2 returns 88 facets, not watertight, 20 crossing pairs.
program fossil_test_resolve_self

use fossil, only : surface_stl_object, facet_object, intersection_pair_t, BOOL_UNION, BOOL_STATUS_OK, &
                   BOOL_STATUS_NOT_CLOSED
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)               :: a, b, c, soup, resolved, cube, holed, dragon
type(facet_object), allocatable        :: facets(:)
type(facet_object), pointer            :: f, g
type(intersection_pair_t), allocatable :: pairs(:)
type(vector_R8P), parameter            :: CORNER = vector_R8P(1._R8P, 1._R8P, 1._R8P)
type(vector_R8P), parameter            :: APEX = vector_R8P(0.3_R8P, 0.35_R8P, -0.3_R8P)
real(R8P)                              :: v_in, pocket, zmin
integer(I4P)                           :: status, i, j, n, cut
logical                                :: ok(6)

! 1. two crossing cubes in one surface
call a%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
b = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)) ; call b%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%analyze
resolved = soup
call resolved%resolve_self_intersections(status=status)
call resolved%find_self_intersections(pairs=pairs)
call a%boolean(other=b, op=BOOL_UNION)
print '(A,I0,A,I0,A,F19.16,A,L1,A,I0,A,I0)', '1. status ', status, ', facets ', resolved%get_facets_number(), &
   ', volume ', resolved%get_volume(), ', watertight ', resolved%is_watertight(), ', shells ', &
   resolved%get_shells_number(), ', pairs ', size(pairs)
print '(A,I0,A,F19.16)', '   union of the shells: facets ', a%get_facets_number(), ', volume ', a%get_volume()
ok(1) = status == BOOL_STATUS_OK .and. resolved%is_watertight() .and. resolved%get_shells_number() == 1 .and. &
        size(pairs) == 0 .and. abs(resolved%get_volume() - a%get_volume()) < 1.e-12_R8P .and. &
        abs(area(resolved) - area(a)) < 1.e-12_R8P
print '(A,L1)', '1. crossing cubes = union of the shells: ', ok(1)

! 2. the folded cube
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
allocate(facets(cube%get_facets_number()))
pocket = 0._R8P
do i = 1, size(facets)
   f => cube%facet_at(i) ; facets(i) = f
   do j = 1, 3
      if (is_corner(facets(i)%vertex(j))) then ! the cone over the section of facet i's two other edges from the corner
         pocket = pocket + cone(facets(i)%vertex(mod(j, 3) + 1), facets(i)%vertex(mod(j + 1, 3) + 1))
         facets(i)%vertex(j) = APEX
      endif
   enddo
   call facets(i)%compute_normal
enddo
call cube%adopt_facets(facets=facets)
v_in = cube%get_volume()
call cube%find_self_intersections(pairs=pairs)
print '(A,I0,A,F19.16,A,F19.16)', '2. input: crossing pairs ', size(pairs), ', signed volume ', v_in, ', pocket ', pocket
call cube%resolve_self_intersections(status=status)
call cube%find_self_intersections(pairs=pairs)
zmin = huge(1._R8P)
do i = 1, cube%get_facets_number()
   g => cube%facet_at(i)
   do j = 1, 3
      zmin = min(zmin, g%vertex(j)%z)
   enddo
enddo
print '(A,I0,A,I0,A,F19.16,A,L1,A,I0,A,I0,A,ES10.3)', '2. status ', status, ', facets ', cube%get_facets_number(), &
   ', volume ', cube%get_volume(), ', watertight ', cube%is_watertight(), ', shells ', cube%get_shells_number(), &
   ', pairs ', size(pairs), ', z min ', zmin
ok(2) = status == BOOL_STATUS_OK .and. cube%is_watertight() .and. cube%get_shells_number() == 1 .and. &
        size(pairs) == 0 .and. zmin >= 0._R8P .and. abs(cube%get_volume() - (v_in + pocket)) < 1.e-12_R8P
print '(A,L1)', '2. folded cube = the region wound once: ', ok(2)

! 3. a clean cube is unchanged
call a%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
b = a
call b%resolve_self_intersections(status=status)
ok(3) = status == BOOL_STATUS_OK .and. b%get_facets_number() == a%get_facets_number()
if (ok(3)) then
   do i = 1, a%get_facets_number()
      f => a%facet_at(i) ; g => b%facet_at(i)
      do j = 1, 3
         ok(3) = ok(3) .and. f%vertex(j)%x == g%vertex(j)%x .and. f%vertex(j)%y == g%vertex(j)%y .and. &
                 f%vertex(j)%z == g%vertex(j)%z
      enddo
   enddo
endif
print '(A,L1)', '3. clean cube unchanged: ', ok(3)

! 4. not closed: the crossing cubes with a facet of the first cube's face x = 0 removed
allocate(facets(soup%get_facets_number() - 1)) ! adopt_facets took the earlier array
cut = 0 ; n = 0
do i = 1, soup%get_facets_number()
   f => soup%facet_at(i)
   if (cut == 0 .and. all([f%vertex(1)%x, f%vertex(2)%x, f%vertex(3)%x] == 0._R8P)) then
      cut = i ; cycle
   endif
   n = n + 1 ; facets(n) = f
enddo
call holed%adopt_facets(facets=facets)
n = holed%get_facets_number() ; v_in = holed%get_volume()
call holed%resolve_self_intersections(status=status)
print '(A,I0,A,I0,A,I0)', '4. removed facet ', cut, ', status ', status, ', facets ', holed%get_facets_number()
ok(4) = status == BOOL_STATUS_NOT_CLOSED .and. holed%get_facets_number() == n .and. holed%get_volume() == v_in
print '(A,L1)', '4. not closed: NOT_CLOSED, unchanged: ', ok(4)

! 5. the dragon: its cuts cross at triple points (issue #91), it is open: NOT_CLOSED, unchanged
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
n = dragon%get_facets_number() ; v_in = dragon%get_volume()
call dragon%resolve_self_intersections(status=status)
print '(A,I0,A,I0)', '5. dragon: status ', status, ', facets ', dragon%get_facets_number()
ok(5) = status == BOOL_STATUS_NOT_CLOSED .and. dragon%get_facets_number() == n .and. dragon%get_volume() == v_in
print '(A,L1)', '5. dragon: NOT_CLOSED, unchanged: ', ok(5)

! 6. three cubes crossing pairwise, a triple point inside a facet of each (issue #91)
call a%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
b = a ; c = a
call b%translate(delta=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)) ; call b%analyze
call c%translate(delta=vector_R8P(0.6_R8P, -0.15_R8P, 0.45_R8P)) ; call c%analyze
soup = a
call soup%merge_solids(other=b) ; call soup%merge_solids(other=c) ; call soup%analyze
resolved = soup
call resolved%resolve_self_intersections(status=status)
call resolved%find_self_intersections(pairs=pairs)
call a%boolean(other=b, op=BOOL_UNION) ; call a%boolean(other=c, op=BOOL_UNION)
print '(A,I0,A,I0,A,F19.16,A,L1,L1,A,I0,A,I0)', '6. status ', status, ', facets ', resolved%get_facets_number(), &
   ', volume ', resolved%get_volume(), ', watertight, manifold ', resolved%is_watertight(), resolved%is_manifold(), &
   ', shells ', resolved%get_shells_number(), ', pairs ', size(pairs)
print '(A,I0,A,F19.16,A,ES10.3,A,ES10.3)', '   boolean chain: facets ', a%get_facets_number(), ', volume ', &
   a%get_volume(), ', |dV| ', abs(resolved%get_volume() - a%get_volume()), ', |dA| ', abs(area(resolved) - area(a))
ok(6) = status == BOOL_STATUS_OK .and. resolved%is_watertight() .and. resolved%is_manifold() .and. &
        resolved%get_shells_number() == 1 .and. size(pairs) == 0 .and. &
        abs(resolved%get_volume() - a%get_volume()) < 1.e-12_R8P .and. abs(area(resolved) - area(a)) < 1.e-12_R8P
print '(A,L1)', '6. three crossing cubes = the boolean chain: ', ok(6)
print '(A,L1)', 'Are all tests passed? ', all(ok)

contains
   function area(s) result(ar)
   !< Total area of the facets.
   type(surface_stl_object), intent(in) :: s
   real(R8P)                            :: ar
   type(facet_object), pointer          :: h
   integer(I4P)                         :: k

   ar = 0._R8P
   do k = 1, s%get_facets_number()
      h => s%facet_at(k) ; ar = ar + h%area()
   enddo
   endfunction area

   pure function is_corner(v) result(yes)
   !< True for the moved corner (1, 1, 1).
   type(vector_R8P), intent(in) :: v
   logical                      :: yes

   yes = v%x == CORNER%x .and. v%y == CORNER%y .and. v%z == CORNER%z
   endfunction is_corner

   pure function cone(p, q) result(vol)
   !< Volume of the tetrahedron APEX, p', q' with p', q' the points where the segments APEX-p and APEX-q cross z = 0
   !< (p, q at z >= 0): the part below z = 0 of the cone of facet (APEX, p, q), with the sign of the facet's turn as
   !< seen from below. The sum over the facets around the corner is the volume of the pocket below the bottom face.
   type(vector_R8P), intent(in) :: p, q
   real(R8P)                    :: vol
   real(R8P)                    :: tp, tq, px, py, qx, qy

   tp = -APEX%z / (p%z - APEX%z) ; tq = -APEX%z / (q%z - APEX%z)
   px = APEX%x + tp * (p%x - APEX%x) - APEX%x ; py = APEX%y + tp * (p%y - APEX%y) - APEX%y
   qx = APEX%x + tq * (q%x - APEX%x) - APEX%x ; qy = APEX%y + tq * (q%y - APEX%y) - APEX%y
   vol = (px * qy - py * qx) / 2._R8P * (-APEX%z) / 3._R8P
   endfunction cone
endprogram fossil_test_resolve_self
