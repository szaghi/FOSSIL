!< FOSSIL test: SAT triangle-AABB overlap predicate.
!<
!< Verifies the 13-axis test in `fossil_utils%triangle_overlaps_aabb` on hand-built
!< geometry covering every failure mode of the old "all-three-vertices-inside" rule:
!<
!<   1. Triangle entirely inside the box                              -> keep.
!<   2. Triangle entirely outside on the +x side                      -> drop.
!<   3. Triangle with 1 vertex inside, 2 outside                      -> keep.
!<   4. Triangle straddling a face, all 3 vertices outside the box    -> keep.
!<   5. Triangle in the plane x+y+z=2.9 whose body straddles the
!<      (+x,+y,+z) corner of the box. All three vertices are outside
!<      the box on every face axis, yet the triangle clips the corner.
!<      The old all-vertices-inside test wrongly drops this           -> keep.
!<   6. Coplanar triangle touching the +x face from outside           -> keep.
!<   7. Triangle separated from the box on the y axis                 -> drop.
!<   8. Triangle separated from the box on the z axis                 -> drop.
!<   9. Load-bearing edge-cross case. Triangle in plane z=0.9 whose
!<      vertices are A=(3,0,0.9), B=(0,3,0.9), C=(4,4,0.9). Triangle
!<      bbox overlaps the box on every face axis, plane z=0.9 cuts
!<      the box, but the triangle body lies entirely on the side
!<      x+y>3 while the box cut z=0.9 reaches at most x+y=2.
!<      Categories 1 & 2 cannot reject; only an edge-cross axis can  -> drop.

program fossil_test_clip_sat

use fossil_utils, only : triangle_overlaps_aabb
use penf,         only : R8P
use vecfor,       only : vector_R8P

implicit none

type(vector_R8P) :: bmin, bmax
type(vector_R8P) :: v1, v2, v3
logical          :: are_tests_passed(9)

are_tests_passed = .false.

! Unit cube centred on the origin.
bmin = vector_R8P(-1._R8P, -1._R8P, -1._R8P)
bmax = vector_R8P( 1._R8P,  1._R8P,  1._R8P)

! 1. Fully inside.
v1 = vector_R8P(-0.5_R8P, -0.5_R8P, 0._R8P)
v2 = vector_R8P( 0.5_R8P, -0.5_R8P, 0._R8P)
v3 = vector_R8P( 0._R8P,   0.5_R8P, 0._R8P)
are_tests_passed(1) = triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 2. Fully outside, on +x side. Box face-normal axis separates.
v1 = vector_R8P(2._R8P, -0.5_R8P, -0.5_R8P)
v2 = vector_R8P(3._R8P,  0.5_R8P,  0._R8P)
v3 = vector_R8P(2.5_R8P, 0._R8P,  0.5_R8P)
are_tests_passed(2) = .not. triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 3. One vertex inside, two outside. Old predicate wrongly dropped this.
v1 = vector_R8P( 0._R8P,  0._R8P, 0._R8P)
v2 = vector_R8P( 3._R8P,  0._R8P, 0._R8P)
v3 = vector_R8P( 0._R8P,  3._R8P, 0._R8P)
are_tests_passed(3) = triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 4. All 3 vertices outside the box yet the body straddles a face.
!    Triangle in the z = 0 plane that wraps around the +x face.
v1 = vector_R8P( 2._R8P, -2._R8P, 0._R8P)
v2 = vector_R8P( 2._R8P,  2._R8P, 0._R8P)
v3 = vector_R8P(-2._R8P,  0._R8P, 0._R8P)
are_tests_passed(4) = triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 5. Triangle in plane x+y+z=2.9 whose vertices are A=(2.9,0,0), B=(0,2.9,0),
!    C=(0,0,2.9). Each vertex has two coordinates inside [-1,1] and one outside,
!    so the OLD all-three-vertices-inside test dropped this. The plane just barely
!    crosses the +x+y+z corner (offset 2.9 vs box radius 3 along (1,1,1)), and the
!    triangle body contains the cut region.
v1 = vector_R8P(2.9_R8P, 0._R8P, 0._R8P)
v2 = vector_R8P(0._R8P, 2.9_R8P, 0._R8P)
v3 = vector_R8P(0._R8P, 0._R8P, 2.9_R8P)
are_tests_passed(5) = triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 6. Coplanar triangle resting on the +x face of the box (closed-set contact).
v1 = vector_R8P(1._R8P, -0.5_R8P, -0.5_R8P)
v2 = vector_R8P(1._R8P,  0.5_R8P, -0.5_R8P)
v3 = vector_R8P(1._R8P,  0._R8P,   0.5_R8P)
are_tests_passed(6) = triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 7. Fully outside on +y side.
v1 = vector_R8P(-0.5_R8P, 2._R8P, -0.5_R8P)
v2 = vector_R8P( 0.5_R8P, 2._R8P,  0._R8P)
v3 = vector_R8P( 0._R8P,  3._R8P,  0.5_R8P)
are_tests_passed(7) = .not. triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 8. Fully outside on -z side.
v1 = vector_R8P(-0.5_R8P, -0.5_R8P, -2._R8P)
v2 = vector_R8P( 0.5_R8P,  0.5_R8P, -2._R8P)
v3 = vector_R8P( 0._R8P,   0._R8P,  -3._R8P)
are_tests_passed(8) = .not. triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

! 9. Load-bearing edge-cross case. Plane z=0.9 cuts the box. Triangle bbox
!    overlaps the box on every face axis. But the triangle body lies on the
!    x+y>3 side of plane z=0.9, while the box cut at z=0.9 reaches at most
!    x+y=2. Categories 1 & 2 alone cannot reject this; only one of the 9
!    edge-cross axes finds the separator. If the predicate omits category 3,
!    it falsely reports overlap here.
v1 = vector_R8P(3._R8P, 0._R8P, 0.9_R8P)
v2 = vector_R8P(0._R8P, 3._R8P, 0.9_R8P)
v3 = vector_R8P(4._R8P, 4._R8P, 0.9_R8P)
are_tests_passed(9) = .not. triangle_overlaps_aabb(bmin, bmax, v1, v2, v3)

print '(A,9L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_clip_sat
