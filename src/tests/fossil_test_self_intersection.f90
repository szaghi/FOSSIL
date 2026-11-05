!< FOSSIL test: self-intersection detection (issue #18 §1.2).
!<
!< Five invariants:
!<   1. Hand-crafted intersecting pair — build two facets that cross along a
!<      known segment, wrap them as a 2-facet surface, and assert the detector
!<      finds exactly one pair with segment endpoints matching the analytic
!<      intersection within tolerance.
!<   2. Clean meshes return zero intersections — cube.stl (12 facets, exact),
!<      cube-05.stl and a faceted cone (`make_apex_body(256, 15)`, 510 facets,
!<      curved, with sliver base-fan triangles) must report zero pairs. This is the no-false-positive guarantee against the
!<      adjacency filter (every edge in a closed mesh is shared by two
!<      facets; if the filter were too narrow, every shared edge would
!<      register as a hit).
!<   3. Dirty mesh is detected non-empty — the dragon.stl fixture in this
!<      repo has known real self-intersections (~759 pairs at the time of
!<      writing). Asserting > 0 catches silent regressions where the broad
!<      phase or narrow phase falsely prunes everything. The exact count is
!<      not asserted because it depends on EPS-sensitive behavior at the
!<      narrow phase boundary; only its non-zero-ness is load-bearing.
!<
!<   4. Facets sharing one vertex (issue #46): they were skipped before the narrow phase, so crossings away from
!<      the shared vertex went unseen. Two triangles sharing the origin and crossing along (0,0,0)-(0.5,0.5,0) are
!<      one pair with that segment; two sharing it and touching only there are none. Invariant 2 covers the other
!<      side: the fans around the vertices of closed meshes (the cone apex: 256 facets) still report none. Not
!<      covered: two sheets crossing along a shared edge (a non-manifold edge, e.g. two cube-01 merged with the
!<      offset (0.5, 0.5, 0.5)): pairs sharing an edge are still skipped.
!<
!<   5. One shared vertex, decided exactly (issue #55): pair_classify's on-plane snap (1e-10) reported crossings on
!<      near-coplanar pairs that meet at the shared vertex only. Four such pairs recorded from the alpha wrap of a
!<      marching-cubes sphere (issue #50; checked in exact rational arithmetic: they touch at the vertex only) report
!<      none; a pair sharing a segment along an edge reports one; a true crossing tilted by 1e-11 rad (below the snap,
!<      above the 1e-12 rad coplanar skip) reports one.
!<
!< Note on dragon.stl: a future PR may either replace this fixture with a
!< cleaned version or use it as ground truth for §1.2 resolution. Either way,
!< invariant 3 documents the current state explicitly so it doesn't surprise.

program fossil_test_self_intersection

use fossil, only : surface_stl_object, intersection_pair_t
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
use fossil_sign_bench_bodies, only : body_t, make_apex_body, to_surface

implicit none

type(surface_stl_object) :: pair_surface, cube, dragon, cube05, cone
type(body_t) :: cb
type(intersection_pair_t), allocatable :: pairs(:)
type(facet_object), allocatable :: tri(:)
type(vector_R8P) :: p_expected_lo, p_expected_hi
type(vector_R8P) :: p_lo, p_hi
real(R8P) :: err_lo, err_hi
integer(I4P) :: status
integer(I4P) :: n_dragon
logical :: are_tests_passed(5)

real(R8P), parameter :: TOL_SEGMENT = 1.0e-12_R8P  ! analytic intersection — full FP precision

are_tests_passed = .false.

! ----- 1. Hand-crafted intersecting pair.
!         Triangle A: in z=0 plane, vertices (0,0,0), (2,0,0), (1,2,0).
!         Triangle B: in y=0.5 plane, vertices (0.5,0.5,-1), (1.5,0.5,-1), (1,0.5,1).
!         They cross along the segment y=0.5, z=0, x in [0.75, 1.25].
allocate(tri(2))
tri(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(3) = vector_R8P(1._R8P, 2._R8P, 0._R8P)
tri(1)%id        = 1_I4P
tri(2)%vertex(1) = vector_R8P(0.5_R8P, 0.5_R8P, -1._R8P)
tri(2)%vertex(2) = vector_R8P(1.5_R8P, 0.5_R8P, -1._R8P)
tri(2)%vertex(3) = vector_R8P(1.0_R8P, 0.5_R8P,  1._R8P)
tri(2)%id        = 2_I4P
call pair_surface%adopt_facets(facets=tri)

call pair_surface%find_self_intersections(pairs=pairs, status=status)

! Expected segment endpoints (order-agnostic).
p_expected_lo = vector_R8P(0.75_R8P, 0.5_R8P, 0._R8P)
p_expected_hi = vector_R8P(1.25_R8P, 0.5_R8P, 0._R8P)

if (status == 0_I4P .and. size(pairs) == 1) then
   ! Order the returned segment so its low-x endpoint comes first.
   if (pairs(1)%p%x <= pairs(1)%q%x) then
      p_lo = pairs(1)%p ; p_hi = pairs(1)%q
   else
      p_lo = pairs(1)%q ; p_hi = pairs(1)%p
   endif
   err_lo = sqrt((p_lo%x - p_expected_lo%x)**2 + (p_lo%y - p_expected_lo%y)**2 + (p_lo%z - p_expected_lo%z)**2)
   err_hi = sqrt((p_hi%x - p_expected_hi%x)**2 + (p_hi%y - p_expected_hi%y)**2 + (p_hi%z - p_expected_hi%z)**2)
   are_tests_passed(1) = pairs(1)%a == 1 .and. pairs(1)%b == 2 .and. &
                         err_lo <= TOL_SEGMENT .and. err_hi <= TOL_SEGMENT
endif

! ----- 2. Clean meshes return zero intersections.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube%find_self_intersections(pairs=pairs, status=status)
if (.not. (status == 0_I4P .and. size(pairs) == 0)) goto 100

call cube05%load_from_file(file_name='src/tests/cube-05.stl', guess_format=.true.)
call cube05%find_self_intersections(pairs=pairs, status=status)
print '(A,I0)', 'cube-05: nfound=', size(pairs)
if (.not. (status == 0_I4P .and. size(pairs) == 0)) goto 100

call make_apex_body(256, 15._R8P, 'cone 15 (256)', cb)
call to_surface(cb%soup, cone)
call cone%find_self_intersections(pairs=pairs, status=status)
print '(A,I0)', 'cone 15 (256): nfound=', size(pairs)
if (.not. (status == 0_I4P .and. size(pairs) == 0)) goto 100

are_tests_passed(2) = .true.

100 continue

! ----- 3. Dirty mesh (dragon) reports > 0 real intersections.
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call dragon%find_self_intersections(pairs=pairs, status=status)
n_dragon = size(pairs)
are_tests_passed(3) = (status == 0_I4P .and. n_dragon > 0_I4P)

print '(A,I0)', 'dragon: nfound=', n_dragon

! ----- 4. Facets sharing one vertex.
!         A: z=0, (0,0,0), (2,0,0), (0,2,0). B: plane x=y through the shared origin, crossing A along
!         (0,0,0)-(0.5,0.5,0). C: shares the origin, on the other side, touching A there only.
if (allocated(tri)) deallocate(tri)
allocate(tri(2))
tri(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(3) = vector_R8P(0._R8P, 2._R8P, 0._R8P)
tri(1)%id        = 1_I4P
tri(2)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(2)%vertex(2) = vector_R8P(0.5_R8P, 0.5_R8P, -1._R8P)
tri(2)%vertex(3) = vector_R8P(0.5_R8P, 0.5_R8P, 1._R8P)
tri(2)%id        = 2_I4P
call pair_surface%adopt_facets(facets=tri)
call pair_surface%find_self_intersections(pairs=pairs, status=status)
are_tests_passed(4) = status == 0_I4P .and. size(pairs) == 1
if (are_tests_passed(4)) are_tests_passed(4) = &
   min(norm2([pairs(1)%p%x, pairs(1)%p%y, pairs(1)%p%z]), norm2([pairs(1)%q%x, pairs(1)%q%y, pairs(1)%q%z])) <= TOL_SEGMENT &
   .and. abs(max(pairs(1)%p%x, pairs(1)%q%x) - 0.5_R8P) <= TOL_SEGMENT .and. &
   abs(max(pairs(1)%p%y, pairs(1)%q%y) - 0.5_R8P) <= TOL_SEGMENT .and. &
   max(abs(pairs(1)%p%z), abs(pairs(1)%q%z)) <= TOL_SEGMENT
print '(A,I0)', 'shared vertex, crossing: nfound=', size(pairs)
allocate(tri(2)) ! adopt_facets took the previous array
tri(1)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(2) = vector_R8P(2._R8P, 0._R8P, 0._R8P)
tri(1)%vertex(3) = vector_R8P(0._R8P, 2._R8P, 0._R8P)
tri(1)%id        = 1_I4P
tri(2)%vertex(1) = vector_R8P(0._R8P, 0._R8P, 0._R8P)
tri(2)%vertex(2) = vector_R8P(-1._R8P, 0._R8P, 1._R8P)
tri(2)%vertex(3) = vector_R8P(0._R8P, -1._R8P, 1._R8P)
tri(2)%id        = 2_I4P
call pair_surface%adopt_facets(facets=tri)
call pair_surface%find_self_intersections(pairs=pairs, status=status)
print '(A,I0)', 'shared vertex, touching: nfound=', size(pairs)
are_tests_passed(4) = are_tests_passed(4) .and. status == 0_I4P .and. size(pairs) == 0

! ----- 5. One shared vertex, exact (issue #55).
are_tests_passed(5) = count_pairs([-5.77306988902962059E-01_R8P, -5.53817895460635534E-01_R8P, 6.46182104539364532E-01_R8P, &
                                   -5.53817895460635534E-01_R8P, -5.78817895460635556E-01_R8P, 6.44968696221600912E-01_R8P, &
                                   -5.53817895460635534E-01_R8P, -5.77306988902962059E-01_R8P, 6.46182104539364532E-01_R8P], &
                                  [-5.52306988902962037E-01_R8P, -5.78817895460635556E-01_R8P, 6.46182104539364532E-01_R8P, &
                                   -5.46177596370816087E-01_R8P, -5.53817895460635534E-01_R8P, 6.71182104539364444E-01_R8P, &
                                   -5.53817895460635534E-01_R8P, -5.77306988902962059E-01_R8P, 6.46182104539364532E-01_R8P], &
                                  'recorded, edge on the other plane') == 0
are_tests_passed(5) = count_pairs([-2.88178954606355120E-02_R8P, -5.28817895460635512E-01_R8P, 8.82476194912600542E-01_R8P, &
                                   -2.88178954606355120E-02_R8P, -5.48116400805985893E-01_R8P, 8.71182104539364621E-01_R8P, &
                                   -3.81789546063560081E-03_R8P, -5.28817895460635512E-01_R8P, 8.82476194912600542E-01_R8P], &
                                  [-2.88178954606355120E-02_R8P, -5.28817895460635512E-01_R8P, 8.82476194912600542E-01_R8P, &
                                   -5.38178954606354232E-02_R8P, -5.05398251952148159E-01_R8P, 8.96182104539364532E-01_R8P, &
                                   -5.38178954606354232E-02_R8P, -5.28817895460635512E-01_R8P, 8.82476194912600542E-01_R8P], &
                                  'recorded, edge on the other plane') == 0 .and. are_tests_passed(5)
are_tests_passed(5) = count_pairs([-6.78817895460635534E-01_R8P, 4.21182104539364444E-01_R8P, -6.46405453701377919E-01_R8P, &
                                   -7.03817895460635556E-01_R8P, 4.08608085662430343E-01_R8P, -6.28817895460635490E-01_R8P, &
                                   -6.96405453701377963E-01_R8P, 4.21182104539364444E-01_R8P, -6.28817895460635490E-01_R8P], &
                                  [-7.03817895460635556E-01_R8P, 4.21182104539364444E-01_R8P, -6.21405453701377786E-01_R8P, &
                                   -7.03817895460635556E-01_R8P, 4.46182104539364577E-01_R8P, -6.06667839126415420E-01_R8P, &
                                   -6.96405453701377963E-01_R8P, 4.21182104539364444E-01_R8P, -6.28817895460635490E-01_R8P], &
                                  'recorded, straddling') == 0 .and. are_tests_passed(5)
are_tests_passed(5) = count_pairs([-6.03817895460635579E-01_R8P, 4.21182104539364444E-01_R8P, -7.21405453700546428E-01_R8P, &
                                   -6.28817895460635490E-01_R8P, 4.08608085662430343E-01_R8P, -7.03817895460635556E-01_R8P, &
                                   -6.21405453701377786E-01_R8P, 4.21182104539364444E-01_R8P, -7.03817895460635556E-01_R8P], &
                                  [-6.28817895460635490E-01_R8P, 4.21182104539364444E-01_R8P, -6.96405453701377963E-01_R8P, &
                                   -6.28817895460635490E-01_R8P, 4.46182104539364577E-01_R8P, -6.81667839126415487E-01_R8P, &
                                   -6.21405453701377786E-01_R8P, 4.21182104539364444E-01_R8P, -7.03817895460635556E-01_R8P], &
                                  'recorded, straddling') == 0 .and. are_tests_passed(5)
!         a segment shared along an edge: B's edge (0,0,0)-(1,0,0) lies on A's edge (0,0,0)-(2,0,0)
are_tests_passed(5) = count_pairs([0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P], &
                                  [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], &
                                  'segment along an edge') == 1 .and. are_tests_passed(5)
!         a true crossing tilted by 1e-11 rad: B's free vertices 1e-11 above and below A's plane z = 0
are_tests_passed(5) = count_pairs([0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P], &
                                  [0._R8P, 0._R8P, 0._R8P, 1._R8P, 0.5_R8P, 1.e-11_R8P, 0.5_R8P, 1._R8P, -1.e-11_R8P], &
                                  'crossing tilted by 1e-11 rad') == 1 .and. are_tests_passed(5)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   function count_pairs(va, vb, label) result(n)
   !< Self-intersections of the two-facet surface (va, vb): 9 coordinates each.
   real(R8P),    intent(in) :: va(9), vb(9) !< Vertices of the two facets.
   character(*), intent(in) :: label        !< Case, printed.
   integer(I4P)             :: n
   type(facet_object), allocatable :: two(:)
   type(surface_stl_object)        :: s2
   type(intersection_pair_t), allocatable :: found(:)
   integer(I4P)                    :: v

   allocate(two(2))
   do v = 1, 3
      two(1)%vertex(v) = vector_R8P(va(3 * v - 2), va(3 * v - 1), va(3 * v))
      two(2)%vertex(v) = vector_R8P(vb(3 * v - 2), vb(3 * v - 1), vb(3 * v))
   enddo
   two(1)%id = 1_I4P ; two(2)%id = 2_I4P
   call s2%adopt_facets(facets=two)
   call s2%find_self_intersections(pairs=found)
   n = size(found)
   print '(A,A,A,I0)', 'shared vertex, ', label, ': nfound=', n
   endfunction count_pairs
endprogram fossil_test_self_intersection
