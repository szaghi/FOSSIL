!< FOSSIL test: exact triangle-triangle classification and canonical points (issue #33, step 2).
program fossil_test_tritri_exact
!< FOSSIL test: exact triangle-triangle classification and canonical points (issue #33, step 2).
!<
!< Hand-built pairs, each with its expected PAIR_* label and, where one exists, its common
!< segment:
!<  - separated (parallel planes; and planes crossing outside both triangles);
!<  - proper crossing: segment equal to the analytic one;
!<  - vertex on the other plane, rest on one side: touching, the point is that vertex bit for bit;
!<  - vertex on the other plane, rest straddling: crossing, one endpoint is that vertex bit for bit;
!<  - edge lying in the other plane through its interior: crossing along that edge;
!<  - coplanar.
!< Canonical points: the point where an edge of B crosses A, computed from the two B triangles
!< that share that edge (with their vertices in different orders), is the same bit for bit; the
!< edge-plane and edge-edge formulas do not depend on the order of their inputs.

use fossil_tritri_exact, only : pair_classify, edge_plane_point, edge_edge_point, &
                                PAIR_SEPARATED, PAIR_CROSSING, PAIR_TOUCHING, PAIR_COPLANAR
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

type(vector_R8P) :: a(3), b(3), b2(3), p, q, p2, q2, x1, x2
integer(I4P)     :: kind, kind2
logical          :: ok

ok = .true.
! A: unit right triangle in z = 0
a = [v(0._R8P, 0._R8P, 0._R8P), v(1._R8P, 0._R8P, 0._R8P), v(0._R8P, 1._R8P, 0._R8P)]

! separated: parallel plane z = 1
b = [v(0._R8P, 0._R8P, 1._R8P), v(1._R8P, 0._R8P, 1._R8P), v(0._R8P, 1._R8P, 1._R8P)]
call check('separated, parallel', PAIR_SEPARATED)
! separated: B crosses z = 0 but far from A
b = [v(5._R8P, 5._R8P, -1._R8P), v(6._R8P, 5._R8P, 1._R8P), v(5._R8P, 6._R8P, 1._R8P)]
call check('separated, planes cross outside', PAIR_SEPARATED)
! proper crossing: vertical triangle in x = 0.25 through A
b = [v(0.25_R8P, -1._R8P, -1._R8P), v(0.25_R8P, 2._R8P, -1._R8P), v(0.25_R8P, 0.5_R8P, 1._R8P)]
call check('proper crossing', PAIR_CROSSING)
ok = ok .and. on_segment_ends(p, q, v(0.25_R8P, 0._R8P, 0._R8P), v(0.25_R8P, 0.75_R8P, 0._R8P))
! vertex on plane, rest above: touching at (0.2, 0.2, 0)
b = [v(0.2_R8P, 0.2_R8P, 0._R8P), v(1._R8P, 0.2_R8P, 1._R8P), v(0.2_R8P, 1._R8P, 1._R8P)]
call check('vertex on plane, touching', PAIR_TOUCHING)
ok = ok .and. same(p, b(1))
! vertex on plane, rest straddling: crossing from that vertex
b = [v(0.2_R8P, 0.2_R8P, 0._R8P), v(0.2_R8P, 0.6_R8P, 1._R8P), v(0.2_R8P, 0.6_R8P, -1._R8P)]
call check('vertex on plane, crossing', PAIR_CROSSING)
ok = ok .and. (same(p, b(1)) .or. same(q, b(1)))
! edge in the plane through A's interior: B's edge (0.1,0.1,0)-(0.5,0.1,0), apex above
b = [v(0.1_R8P, 0.1_R8P, 0._R8P), v(0.5_R8P, 0.1_R8P, 0._R8P), v(0.3_R8P, 0.1_R8P, 1._R8P)]
call check('edge in plane', PAIR_CROSSING)
ok = ok .and. on_segment_ends(p, q, b(1), b(2))
! coplanar
b = [v(0.2_R8P, 0.2_R8P, 0._R8P), v(2._R8P, 0.2_R8P, 0._R8P), v(0.2_R8P, 2._R8P, 0._R8P)]
call check('coplanar', PAIR_COPLANAR)

! canonical points: edge E = (0.3,0.2,-1)-(0.3,0.2,1) of two B triangles crossing A
b  = [v(0.3_R8P, 0.2_R8P, -1._R8P), v(0.3_R8P, 0.2_R8P, 1._R8P), v(0.9_R8P, 0.9_R8P, 0.3_R8P)]
b2 = [v(0.3_R8P, 0.2_R8P, 1._R8P), v(-0.4_R8P, -0.3_R8P, 0.1_R8P), v(0.3_R8P, 0.2_R8P, -1._R8P)]
call pair_classify(ta=a, tb=b, kind=kind, p=p, q=q)
call pair_classify(ta=a, tb=b2, kind=kind2, p=p2, q=q2)
ok = ok .and. kind == PAIR_CROSSING .and. kind2 == PAIR_CROSSING
x1 = edge_plane_point(b(1), b(2), a)
ok = ok .and. (same(p, x1) .or. same(q, x1)) .and. (same(p2, x1) .or. same(q2, x1))
print '(A,L1)', 'shared edge crossing, same bits from both B triangles: ', &
      (same(p, x1) .or. same(q, x1)) .and. (same(p2, x1) .or. same(q2, x1))
! order independence of the formulas
x1 = edge_plane_point(b(1), b(2), a)
x2 = edge_plane_point(b(2), b(1), [a(3), a(1), a(2)])
ok = ok .and. same(x1, x2)
x1 = edge_edge_point(v(0._R8P, 0._R8P, 0._R8P), v(1._R8P, 0.3_R8P, 0._R8P), v(0.2_R8P, 1._R8P, 0._R8P), &
                     v(0.7_R8P, -1._R8P, 0._R8P))
x2 = edge_edge_point(v(0.7_R8P, -1._R8P, 0._R8P), v(0.2_R8P, 1._R8P, 0._R8P), v(1._R8P, 0.3_R8P, 0._R8P), &
                     v(0._R8P, 0._R8P, 0._R8P))
ok = ok .and. same(x1, x2)
print '(A,L1)', 'edge-plane and edge-edge points independent of input order: ', same(x1, x2)

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   pure function v(x, y, z) result(r)
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: r

   r%x = x ; r%y = y ; r%z = z
   endfunction v

   pure function same(a1, a2) result(s)
   type(vector_R8P), intent(in) :: a1, a2
   logical                      :: s

   s = a1%x == a2%x .and. a1%y == a2%y .and. a1%z == a2%z
   endfunction same

   pure function on_segment_ends(p1, q1, e1, e2) result(s)
   !< (p1, q1) equals (e1, e2) in either order, to 1e-14.
   type(vector_R8P), intent(in) :: p1, q1, e1, e2
   logical                      :: s

   s = (near(p1, e1) .and. near(q1, e2)) .or. (near(p1, e2) .and. near(q1, e1))
   endfunction on_segment_ends

   pure function near(a1, a2) result(s)
   type(vector_R8P), intent(in) :: a1, a2
   logical                      :: s

   s = abs(a1%x - a2%x) < 1.e-14_R8P .and. abs(a1%y - a2%y) < 1.e-14_R8P .and. abs(a1%z - a2%z) < 1.e-14_R8P
   endfunction near

   subroutine check(label, expected)
   character(*), intent(in) :: label
   integer(I4P), intent(in) :: expected

   call pair_classify(ta=a, tb=b, kind=kind, p=p, q=q)
   print '(A,A,I0,A,I0)', label, ': kind ', kind, '  expected ', expected
   ok = ok .and. kind == expected
   endsubroutine check
endprogram fossil_test_tritri_exact
