!< FOSSIL test: find_self_intersections with exactly degenerate facets (issue #98, plan #99 step 3).
module fossil_test_selfint_degenerate_lib
!< Helpers of `fossil_test_selfint_degenerate`: surfaces from coordinate lists, the cube with a T-junction sliver.
use fossil, only : surface_stl_object, facet_object, intersection_pair_t
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: tri, pnt, selfint, cube_t_junction, expect

contains
   pure function tri(p1, p2, p3) result(t)
   !< The triangle (p1, p2, p3) as 9 coordinates.
   real(R8P), intent(in) :: p1(3), p2(3), p3(3) !< Vertices.
   real(R8P)             :: t(9)                !< Coordinates.

   t = [p1, p2, p3]
   endfunction tri

   pure function pnt(p) result(t)
   !< The point-degenerate triangle (p, p, p).
   real(R8P), intent(in) :: p(3) !< Point.
   real(R8P)             :: t(9) !< Coordinates.

   t = [p, p, p]
   endfunction pnt

   subroutine selfint(tris, touching, found, watertight)
   !< Self-intersections of the surface made of the columns of `tris` (9 coordinates each, in this order).
   real(R8P),                              intent(in)            :: tris(:, :)  !< Facets, one per column.
   logical,                                intent(in)            :: touching    !< include_touching.
   type(intersection_pair_t), allocatable, intent(out)           :: found(:)    !< Pairs.
   logical,                                intent(out), optional :: watertight  !< The surface is watertight.
   type(facet_object),        allocatable                        :: f(:)        !< Facets.
   type(surface_stl_object)                                      :: s           !< Surface.
   integer(I4P)                                                  :: i, v        !< Counters.

   allocate(f(size(tris, dim=2)))
   do i = 1, size(tris, dim=2)
      do v = 1, 3
         f(i)%vertex(v) = vector_R8P(tris(3 * v - 2, i), tris(3 * v - 1, i), tris(3 * v, i))
      enddo
      f(i)%id = i
   enddo
   call s%adopt_facets(facets=f)
   call s%find_self_intersections(pairs=found, include_touching=touching)
   if (present(watertight)) watertight = s%is_watertight()
   endsubroutine selfint

   pure function cube_t_junction(sliver_first) result(tris)
   !< The unit cube, its top face (z = 1) split along the diagonal (0,0,1)-(1,1,1): one half a single facet, the other
   !< half two facets meeting at the diagonal's midpoint M (a T-junction), closed by the collinear sliver
   !< ((0,0,1), (1,1,1), M). Watertight, outward oriented. The sliver is the first facet or the last.
   logical,   intent(in) :: sliver_first   !< Sliver first (else last).
   real(R8P)             :: tris(9, 14)    !< Facets.
   real(R8P)             :: q(3, 4, 5)     !< The five other faces, counter-clockwise from outside.
   real(R8P)             :: a(3), b(3), c(3), d(3), m(3) !< Top face corners, midpoint.
   integer(I4P)          :: k, n           !< Counters.

   q(:, :, 1) = reshape([0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P], &
                        [3, 4])
   q(:, :, 2) = reshape([0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P], &
                        [3, 4])
   q(:, :, 3) = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P], &
                        [3, 4])
   q(:, :, 4) = reshape([1._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P], &
                        [3, 4])
   q(:, :, 5) = reshape([0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P], &
                        [3, 4])
   a = [0._R8P, 0._R8P, 1._R8P] ; b = [1._R8P, 0._R8P, 1._R8P] ; c = [1._R8P, 1._R8P, 1._R8P]
   d = [0._R8P, 1._R8P, 1._R8P] ; m = [0.5_R8P, 0.5_R8P, 1._R8P]
   n = 0_I4P
   if (sliver_first) then
      n = n + 1_I4P ; tris(:, n) = tri(a, c, m)
   endif
   do k = 1, 5
      n = n + 1_I4P ; tris(:, n) = tri(q(:, 1, k), q(:, 2, k), q(:, 3, k))
      n = n + 1_I4P ; tris(:, n) = tri(q(:, 1, k), q(:, 3, k), q(:, 4, k))
   enddo
   n = n + 1_I4P ; tris(:, n) = tri(a, b, c)
   n = n + 1_I4P ; tris(:, n) = tri(a, m, d)
   n = n + 1_I4P ; tris(:, n) = tri(m, c, d)
   if (.not. sliver_first) then
      n = n + 1_I4P ; tris(:, n) = tri(a, c, m)
   endif
   endfunction cube_t_junction

   subroutine expect(c, ta, tb, kind_plain, kind_touching, ok, x)
   !< Case c: the two-facet surface (ta, tb), in both orders, without and with include_touching: one pair of the given
   !< kind (0: none), at the point x (p = q) if given.
   integer(I4P), intent(in)        :: c                         !< Case.
   real(R8P), intent(in)           :: ta(9), tb(9)              !< Facets.
   integer(I4P), intent(in)        :: kind_plain, kind_touching !< Expected kinds (0: no pair).
   logical,   intent(inout)        :: ok(:)                     !< Case flags.
   real(R8P), intent(in), optional :: x(3)                      !< Expected contact point.
   type(intersection_pair_t), allocatable :: got(:)
   integer(I4P)                    :: order, t, kind
   logical                         :: good

   good = .true.
   do order = 1, 2
      do t = 1, 2
         if (order == 1) then
            call selfint(reshape([ta, tb], [9, 2]), t == 2, got)
         else
            call selfint(reshape([tb, ta], [9, 2]), t == 2, got)
         endif
         kind = kind_plain ; if (t == 2) kind = kind_touching
         if (kind == 0) then
            good = good .and. size(got) == 0
         else
            good = good .and. size(got) == 1
            if (size(got) == 1) then
               good = good .and. got(1)%kind == kind
               if (present(x)) good = good .and. all([got(1)%p%x, got(1)%p%y, got(1)%p%z] == x) .and. &
                                     all([got(1)%q%x, got(1)%q%y, got(1)%q%z] == x)
            endif
         endif
         if (.not. good) then
            print '(A,I0,A,I0,A,L1,A,I0)', '   case ', c, ' order ', order, ' touching ', t == 2, ': pairs ', size(got)
            if (size(got) > 0) print '(A,I0)', '      kind ', got(1)%kind
            ok(c) = .false. ; return
         endif
      enddo
   enddo
   print '(A,I0,A,L1)', '   case ', c, ': ', good
   endsubroutine expect
endmodule fossil_test_selfint_degenerate_lib

!< FOSSIL test: find_self_intersections with exactly degenerate facets (issue #98, plan #99 step 3).
!<
!< A degenerate facet (collinear vertices: the segment between its extreme vertices; coincident ones: a point) has no
!< plane: `pair_cut` cannot classify it. Expected, decided exactly: no contact, no report; a contact at a point interior
!< to the segment and to the other facet (the segment crossing it) SELFINT_CROSSING; a contact only at a boundary point
!< or an end of the segment (or a point facet) SELFINT_TOUCHING, with include_touching only; facets sharing vertices:
!< only the contact beyond the shared vertices (and beyond the shared edge) counts. Every case in both facet orders.
!< F = (0,0,0), (1,0,0), (0,1,0).
!< 1. The #96 probe (the face x = 0 and the segment (0,0,1.5)-(1,1,2), 0.5 apart): none (boxes apart: the broad phase
!<    prunes it, before and after).
!< 2. A segment crossing F's plane at (0.8, 0.8, 0), outside F, its box overlapping F's: none, even with touching.
!< 3. A segment through F's interior, transversal: CROSSING at (0.25, 0.25, 0).
!< 4. A segment in F's plane across F's interior: CROSSING.
!< 5. A segment ending on F's interior: TOUCHING at (0.25, 0.25, 0), with touching only.
!< 6. A segment across F's edge y = 0: TOUCHING with touching only. 7. In F's plane along that edge: the same.
!< 8. Two segments crossing at interior points: CROSSING; collinear overlapping: CROSSING; an end of one inside the
!<    other: TOUCHING with touching only.
!< 9. A point facet in F's interior: TOUCHING with touching only; beside F in its plane: none.
!< 10. Sharing F's vertex (0,0,0): a segment into F's interior, also with the vertex in its middle: CROSSING; along F's
!<    edge: TOUCHING with touching only; outward in the plane and off the plane: none. Two segments sharing an end,
!<    collinear and overlapping: CROSSING; not collinear: none.
!< 11. The unit cube with a T-junction closed by a collinear sliver (watertight): no report; with touching only the two
!<    coplanar T-junction pairs of regular facets, none with the sliver.
!< 12. A balanced fan of 4 facets around one edge, one of them degenerate (no half-plane): no SELFINT_SHEETS.
!< Negative control (base 1147c33 + step 1): every case but 1 fails (2: a false TOUCHING; 3-10: missing or wrong
!< reports, in one order or both; 11: a false TOUCHING of the sliver with a side facet; 12: two false SHEETS).
program fossil_test_selfint_degenerate

use fossil,                             only : intersection_pair_t, SELFINT_CROSSING, SELFINT_TOUCHING
use fossil_test_selfint_degenerate_lib, only : tri, pnt, selfint, cube_t_junction, expect
use penf,                               only : I4P, R8P

implicit none

real(R8P), parameter                   :: O(3) = [0._R8P, 0._R8P, 0._R8P], X1(3) = [1._R8P, 0._R8P, 0._R8P]
real(R8P), parameter                   :: Y1(3) = [0._R8P, 1._R8P, 0._R8P]
real(R8P)                              :: f(9), sa(9)
type(intersection_pair_t), allocatable :: found(:)
logical                                :: ok(12), wt, first
integer(I4P)                           :: k

ok = .true.
f = tri(O, X1, Y1)
! 1. the #96 probe
call expect(1, tri([0._R8P, 1._R8P, 1._R8P], O, [0._R8P, 0._R8P, 1._R8P]), &
            tri([0._R8P, 0._R8P, 1.5_R8P], [0.5_R8P, 0.5_R8P, 1.75_R8P], [1._R8P, 1._R8P, 2._R8P]), 0, 0, ok)
! 2. crossing F's plane outside F, boxes overlapping
call expect(2, f, tri([0.8_R8P, 0.8_R8P, -0.5_R8P], [0.8_R8P, 0.8_R8P, 0.5_R8P], [0.8_R8P, 0.8_R8P, 0.25_R8P]), 0, 0, ok)
! 3. transversal through F's interior
call expect(3, f, tri([0.25_R8P, 0.25_R8P, -0.5_R8P], [0.25_R8P, 0.25_R8P, 0.5_R8P], [0.25_R8P, 0.25_R8P, 0.25_R8P]), &
            SELFINT_CROSSING, SELFINT_CROSSING, ok, [0.25_R8P, 0.25_R8P, 0._R8P])
! 4. in F's plane, across its interior
call expect(4, f, tri([-0.5_R8P, 0.25_R8P, 0._R8P], [0.75_R8P, 0.25_R8P, 0._R8P], [0.125_R8P, 0.25_R8P, 0._R8P]), &
            SELFINT_CROSSING, SELFINT_CROSSING, ok)
! 5. an end on F's interior
call expect(5, f, tri([0.25_R8P, 0.25_R8P, 0._R8P], [0.25_R8P, 0.25_R8P, 1._R8P], [0.25_R8P, 0.25_R8P, 0.5_R8P]), &
            0, SELFINT_TOUCHING, ok, [0.25_R8P, 0.25_R8P, 0._R8P])
! 6. across F's edge y = 0; 7. along it, in the plane
call expect(6, f, tri([0.5_R8P, 0._R8P, -0.5_R8P], [0.5_R8P, 0._R8P, 0.5_R8P], [0.5_R8P, 0._R8P, 0.25_R8P]), &
            0, SELFINT_TOUCHING, ok, [0.5_R8P, 0._R8P, 0._R8P])
call expect(7, f, tri([0.25_R8P, 0._R8P, 0._R8P], [0.75_R8P, 0._R8P, 0._R8P], [0.5_R8P, 0._R8P, 0._R8P]), &
            0, SELFINT_TOUCHING, ok)
! 8. two segments: crossing, collinear overlapping, an end inside the other
sa = tri(O, [1._R8P, 1._R8P, 0._R8P], [0.25_R8P, 0.25_R8P, 0._R8P])
call expect(8, sa, tri(X1, Y1, [0.25_R8P, 0.75_R8P, 0._R8P]), SELFINT_CROSSING, SELFINT_CROSSING, ok, &
            [0.5_R8P, 0.5_R8P, 0._R8P])
call expect(8, sa, tri([0.5_R8P, 0.5_R8P, 0._R8P], [2._R8P, 2._R8P, 0._R8P], [1.5_R8P, 1.5_R8P, 0._R8P]), &
            SELFINT_CROSSING, SELFINT_CROSSING, ok)
call expect(8, sa, tri([0.5_R8P, 0.5_R8P, 0._R8P], [0.5_R8P, 1.5_R8P, 0._R8P], [0.5_R8P, 1._R8P, 0._R8P]), &
            0, SELFINT_TOUCHING, ok, [0.5_R8P, 0.5_R8P, 0._R8P])
! 9. a point facet inside F, beside F
call expect(9, f, pnt([0.25_R8P, 0.25_R8P, 0._R8P]), 0, SELFINT_TOUCHING, ok, [0.25_R8P, 0.25_R8P, 0._R8P])
call expect(9, f, pnt([0.75_R8P, 0.75_R8P, 0._R8P]), 0, 0, ok)
! 10. sharing F's vertex O
call expect(10, f, tri(O, [0.5_R8P, 0.25_R8P, 0._R8P], [0.25_R8P, 0.125_R8P, 0._R8P]), SELFINT_CROSSING, SELFINT_CROSSING, ok)
call expect(10, f, tri([-0.25_R8P, -0.125_R8P, 0._R8P], [0.5_R8P, 0.25_R8P, 0._R8P], O), SELFINT_CROSSING, &
            SELFINT_CROSSING, ok)
call expect(10, f, tri(O, [0.5_R8P, 0._R8P, 0._R8P], [0.25_R8P, 0._R8P, 0._R8P]), 0, SELFINT_TOUCHING, ok)
call expect(10, f, tri(O, [-0.5_R8P, -0.25_R8P, 0._R8P], [-0.25_R8P, -0.125_R8P, 0._R8P]), 0, 0, ok)
call expect(10, f, tri(O, [0.25_R8P, 0.25_R8P, 1._R8P], [0.125_R8P, 0.125_R8P, 0.5_R8P]), 0, 0, ok)
call expect(10, sa, tri(O, [0.5_R8P, 0.5_R8P, 0._R8P], [1.5_R8P, 1.5_R8P, 0._R8P]), SELFINT_CROSSING, SELFINT_CROSSING, ok)
call expect(10, sa, tri(O, [1._R8P, -1._R8P, 0._R8P], [0.5_R8P, -0.5_R8P, 0._R8P]), 0, 0, ok)
! 11. the cube with a T-junction sliver
do k = 1, 2
   first = k == 1
   call selfint(cube_t_junction(first), .false., found, wt)
   ok(11) = ok(11) .and. wt .and. size(found) == 0
   call selfint(cube_t_junction(first), .true., found)
   ok(11) = ok(11) .and. size(found) == 2
   if (size(found) > 0) then
      ok(11) = ok(11) .and. all(found(:)%kind == SELFINT_TOUCHING)
      if (first) then
         ok(11) = ok(11) .and. all(found(:)%a /= 1)
      else
         ok(11) = ok(11) .and. all(found(:)%b /= 14)
      endif
   endif
   print '(A,L1,A,L1,A,I0)', '11. T-junction sliver (sliver first ', first, '): ', ok(11), '  pairs with touching ', &
                             size(found)
enddo
! 12. a balanced edge fan with a degenerate facet
call selfint(reshape([tri(O, X1, [0.5_R8P, 1._R8P, 0._R8P]), tri(X1, O, [0.5_R8P, 0._R8P, 1._R8P]), &
                      tri(X1, O, [0.5_R8P, 0._R8P, 0._R8P]), tri(O, X1, [0.5_R8P, 0._R8P, -1._R8P])], [9, 4]), &
             .true., found)
ok(12) = size(found) == 0
print '(A,L1,A,I0)', '12. edge fan with a degenerate facet: ', ok(12), '  pairs ', size(found)

do k = 1, 12
   if (.not. ok(k)) print '(A,I0,A)', 'case ', k, ' FAILED'
enddo
print '(A,L1)', 'Are all tests passed? ', all(ok)
endprogram fossil_test_selfint_degenerate
