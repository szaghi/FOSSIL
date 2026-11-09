!< FOSSIL test: body-body distance with exactly degenerate facets (issue #96).
!<
!< A degenerate facet (three collinear vertices: a segment; three equal vertices: a point) has no plane, so `pair_cut`
!< cannot classify it: `tri_tri_distance` decides contact by exact segment / point predicates instead. Each case is
!< checked three ways: the expected contact and distance by construction, an exact contact reference written here
!< (orient3d / orient2d, decomposed differently from the library: every non-degenerate projection, dominant-axis
!< interval order), and the distance reference: the minimum of the 15 sub-distances (6 vertex-triangle, 9 edge-edge)
!< by the existing kernels, in the library's order (equal bit for bit when apart).
!< 1. `tri_tri_distance`, both orders: the #96 probe (a segment 0.5 above a cube face's triangle), a segment crossing
!<    the face, segments touching it at an end (interior, a vertex of the face), a segment in the face's plane across its
!<    edge and beside it, a point inside, beside (in the plane) and off the face, two segments crossing, skew, collinear
!<    overlapping and apart, a point on and off a segment, two equal and two distinct points.
!< 2. `surface%distance_to` (tree descent) on bodies with one degenerate facet against the unit cube, and two bodies
!<    with one degenerate facet each: distance and contact against the reference over every facet pair.
!< 3. The 2-D twin `curve%distance_to` with zero-length segments (an open loop of two equal vertices; a square with a
!<    repeated vertex) against a unit square: a point off it, on an edge, on a corner, on an edge's line beyond it,
!<    inside; two equal and two distinct points: distance and contact by construction, both orders.
program fossil_test_body_distance_degenerate

use fossil,               only : surface_stl_object, facet_object, curve_object
use fossil_body_distance, only : tri_tri_distance, closest_on_triangle, closest_segment_segment
use fossil_predicates,    only : orient2d_xy, orient3d_xyz
use penf,                 only : I4P, R8P
use vecfor,               only : vector_R8P

implicit none

real(R8P)                :: face(3, 3), far_a(3, 3), far_b(3, 3)
type(surface_stl_object) :: cube, a, b
type(curve_object)       :: ca, cb
integer(I4P)             :: ncase, nbad
logical                  :: ok1, ok2, ok3

face = reshape([0._R8P, 1._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ! cube facet 2 (x = 0)
far_a = reshape([100._R8P, 100._R8P, 100._R8P, 101._R8P, 100._R8P, 100._R8P, 100._R8P, 101._R8P, 100._R8P], [3, 3])
far_b = -far_a
ncase = 0 ; nbad = 0
! segments against the face triangle (y <= z in the plane x = 0)
call pair(face, tri([0._R8P, 0._R8P, 1.5_R8P], [0.5_R8P, 0.5_R8P, 1.75_R8P], [1._R8P, 1._R8P, 2._R8P]), .false., &
          0.5_R8P)                                                                                  ! the #96 probe
call pair(face, tri([0._R8P, 0.25_R8P, 0.625_R8P], [-0.5_R8P, 0.125_R8P, 0.5_R8P], [0.5_R8P, 0.375_R8P, 0.75_R8P]), &
          .true., 0._R8P)                                                                           ! crossing
call pair(face, tri([0.5_R8P, 0.25_R8P, 0.5_R8P], [0._R8P, 0.25_R8P, 0.5_R8P], [1._R8P, 0.25_R8P, 0.5_R8P]), .true., &
          0._R8P)                                                                                   ! an end inside
call pair(face, tri([0._R8P, 0._R8P, 1._R8P], [0.5_R8P, 0.5_R8P, 1.5_R8P], [1._R8P, 1._R8P, 2._R8P]), .true., &
          0._R8P)                                                                                   ! an end on a vertex
call pair(face, tri([0._R8P, 0.75_R8P, 0.25_R8P], [0._R8P, 0.25_R8P, 0.75_R8P], [0._R8P, 0.5_R8P, 0.5_R8P]), .true., &
          0._R8P)                                                                                   ! in plane, across
call pair(face, tri([0._R8P, 0.75_R8P, 0.25_R8P], [0._R8P, 1._R8P, 0.5_R8P], [0._R8P, 0.875_R8P, 0.375_R8P]), &
          .false., -1._R8P)                                                                         ! in plane, beside
! points against the face
call pair(face, pnt([0._R8P, 0.25_R8P, 0.5_R8P]), .true., 0._R8P)                                ! inside
call pair(face, pnt([0._R8P, 0._R8P, 0.5_R8P]), .true., 0._R8P)                                  ! on an edge
call pair(face, pnt([0._R8P, 0.75_R8P, 0.25_R8P]), .false., -1._R8P)                             ! in plane, beside
call pair(face, pnt([0.5_R8P, 0.25_R8P, 0.5_R8P]), .false., 0.5_R8P)                             ! off the plane
! two degenerate facets
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0._R8P]), &
          tri([1._R8P, 0._R8P, 0._R8P], [0._R8P, 1._R8P, 0._R8P], [0.25_R8P, 0.75_R8P, 0._R8P]), .true., 0._R8P) ! crossing
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0._R8P]), &
          tri([1._R8P, 0._R8P, 1._R8P], [0._R8P, 1._R8P, 1._R8P], [0.25_R8P, 0.75_R8P, 1._R8P]), .false., 1._R8P) ! skew
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], [0.5_R8P, 0.5_R8P, 0.5_R8P]), &
          tri([2._R8P, 2._R8P, 2._R8P], [0.75_R8P, 0.75_R8P, 0.75_R8P], [1.5_R8P, 1.5_R8P, 1.5_R8P]), .true., 0._R8P) ! overlap
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], [0.5_R8P, 0.5_R8P, 0.5_R8P]), &
          tri([2._R8P, 2._R8P, 2._R8P], [3._R8P, 3._R8P, 3._R8P], [2.5_R8P, 2.5_R8P, 2.5_R8P]), .false., -1._R8P) ! apart
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0._R8P]), &
          pnt([0.25_R8P, 0.25_R8P, 0._R8P]), .true., 0._R8P)                                       ! point on a segment
call pair(tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0._R8P]), &
          pnt([2._R8P, 2._R8P, 0._R8P]), .false., sqrt(2._R8P))                                    ! on the line, apart
call pair(pnt([0.25_R8P, 0.5_R8P, 0._R8P]), pnt([0.25_R8P, 0.5_R8P, 0._R8P]), .true., 0._R8P)    ! equal points
call pair(pnt([0.25_R8P, 0.5_R8P, 0._R8P]), pnt([0.25_R8P, 0.5_R8P, 1._R8P]), .false., 1._R8P)   ! distinct points
ok1 = nbad == 0 .and. ncase == 36
print '(A,L1,A,I0,A,I0)', '1. tri_tri_distance with degenerate facets: ', ok1, '  cases ', ncase, '  failures ', nbad

ncase = 0 ; nbad = 0
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call body(b, tri([0._R8P, 0._R8P, 1.5_R8P], [0.5_R8P, 0.5_R8P, 1.75_R8P], [1._R8P, 1._R8P, 2._R8P]), far_b)
call compare(cube, b, .false., 0.5_R8P)                                                          ! the #96 probe
call body(b, tri([0._R8P, 0.25_R8P, 0.625_R8P], [-0.5_R8P, 0.125_R8P, 0.5_R8P], [0.5_R8P, 0.375_R8P, 0.75_R8P]), far_b)
call compare(cube, b, .true., 0._R8P)                                                            ! crossing a face
call body(b, tri([0._R8P, 0._R8P, 1._R8P], [-0.5_R8P, -0.5_R8P, 1.5_R8P], [-1._R8P, -1._R8P, 2._R8P]), far_b)
call compare(cube, b, .true., 0._R8P)                                                            ! touching a corner
call body(b, pnt([-0.5_R8P, 0.25_R8P, 0.5_R8P]), far_b)
call compare(cube, b, .false., 0.5_R8P)                                                          ! a point outside
call body(b, pnt([0._R8P, 0.25_R8P, 0.5_R8P]), far_b)
call compare(cube, b, .true., 0._R8P)                                                            ! a point on a face
call body(a, tri([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0._R8P]), far_a)
call body(b, tri([1._R8P, 0._R8P, 1._R8P], [0._R8P, 1._R8P, 1._R8P], [0.25_R8P, 0.75_R8P, 1._R8P]), far_b)
call compare(a, b, .false., 1._R8P)                                                              ! two skew segments
call body(b, tri([1._R8P, 0._R8P, 0._R8P], [0._R8P, 1._R8P, 0._R8P], [0.25_R8P, 0.75_R8P, 0._R8P]), far_b)
call compare(a, b, .true., 0._R8P)                                                               ! two crossing segments
ok2 = nbad == 0 .and. ncase == 14
print '(A,L1,A,I0,A,I0)', '2. distance_to with degenerate facets: ', ok2, '  cases ', ncase, '  failures ', nbad

ncase = 0 ; nbad = 0
call cb%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call point_curve(ca, 0.5_R8P, 2._R8P) ; call compare_2d(ca, cb, .false., 1._R8P)              ! a point off the square
call point_curve(ca, 0.5_R8P, 1._R8P) ; call compare_2d(ca, cb, .true., 0._R8P)               ! on an edge
call point_curve(ca, 1._R8P, 1._R8P) ; call compare_2d(ca, cb, .true., 0._R8P)                ! on a corner
call point_curve(ca, 2._R8P, 1._R8P) ; call compare_2d(ca, cb, .false., 1._R8P)               ! on an edge's line, beyond
call point_curve(ca, 0.5_R8P, 0.25_R8P) ; call compare_2d(ca, cb, .false., 0.25_R8P)          ! inside (distance to the curve)
call point_curve(cb, 0.5_R8P, 0.25_R8P) ; call compare_2d(ca, cb, .true., 0._R8P)             ! two equal points
call point_curve(cb, 0.5_R8P, 1.25_R8P) ; call compare_2d(ca, cb, .false., 1._R8P)            ! two distinct points
call ca%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P, 1._R8P], &
                    loop_sizes=[5_I4P])                                                          ! a repeated vertex
call cb%adopt_loops(u=[1._R8P, 2._R8P, 2._R8P], v=[2._R8P, 2._R8P, 3._R8P], loop_sizes=[3_I4P])
call compare_2d(ca, cb, .false., 1._R8P)                                                         ! apart from a corner
call point_curve(cb, 1._R8P, 1._R8P) ; call compare_2d(ca, cb, .true., 0._R8P)                ! on the repeated vertex
ok3 = nbad == 0 .and. ncase == 18
print '(A,L1,A,I0,A,I0)', '3. curve distance_to with zero-length segments: ', ok3, '  cases ', ncase, '  failures ', nbad
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3

contains
   subroutine pair(s, t, mexp, dexp)
   !< tri_tri_distance(s, t) and (t, s) against the construction (contact mexp; distance dexp if >= 0) and the references.
   real(R8P), intent(in) :: s(3, 3), t(3, 3)
   logical,   intent(in) :: mexp
   real(R8P), intent(in) :: dexp
   real(R8P)             :: x(3, 3), y(3, 3), d2, px(3), py(3), r2
   logical               :: m, mr
   integer(I4P)          :: o

   do o = 1, 2
      if (o == 1) then
         x = s ; y = t
      else
         x = t ; y = s
      endif
      ncase = ncase + 1
      call tri_tri_distance(x, y, d2, px, py, m)
      mr = ref_meet(x, y) ; r2 = ref_d2(x, y)
      if ((m .neqv. mexp) .or. (mr .neqv. mexp) .or. (m .and. d2 /= 0._R8P) .or. (.not. m .and. d2 /= r2) .or. &
          (dexp >= 0._R8P .and. abs(sqrt(d2) - dexp) > 1.e-15_R8P)) then
         nbad = nbad + 1
         print '(A,I0,A,2L2,A,L2,A,3ES24.16)', '   failure in case ', ncase, ': meet, ref ', m, mr, ', expected ', mexp, &
                                               ', d, ref, expected ', sqrt(d2), sqrt(r2), dexp
      endif
   enddo
   endsubroutine pair

   subroutine compare(s, t, mexp, dexp)
   !< s%distance_to(t) and t%distance_to(s) against the construction and the reference over every facet pair.
   type(surface_stl_object), intent(in) :: s, t
   logical,                  intent(in) :: mexp
   real(R8P),                intent(in) :: dexp

   call compare_one(s, t, mexp, dexp) ; call compare_one(t, s, mexp, dexp)
   endsubroutine compare

   subroutine compare_one(x, y, mexp, dexp)
   !< x%distance_to(y) against the construction and the reference over every facet pair.
   type(surface_stl_object), intent(in) :: x, y
   logical,                  intent(in) :: mexp
   real(R8P),                intent(in) :: dexp
   type(facet_object), pointer          :: f, g
   real(R8P)                            :: d, r2, u(3, 3), w(3, 3)
   logical                              :: h, mr
   integer(I4P)                         :: i, j, k

   ncase = ncase + 1
   call x%distance_to(y, d, intersecting=h)
   mr = .false. ; r2 = huge(1._R8P)
   do i = 1, x%get_facets_number()
      f => x%facet_at(i)
      do j = 1, y%get_facets_number()
         g => y%facet_at(j)
         do k = 1, 3
            u(:, k) = [f%vertex(k)%x, f%vertex(k)%y, f%vertex(k)%z] ; w(:, k) = [g%vertex(k)%x, g%vertex(k)%y, g%vertex(k)%z]
         enddo
         if (ref_meet(u, w)) mr = .true.
         r2 = min(r2, ref_d2(u, w))
      enddo
   enddo
   if (mr) r2 = 0._R8P
   if ((h .neqv. mexp) .or. (mr .neqv. mexp) .or. d /= sqrt(r2) .or. abs(d - dexp) > 1.e-15_R8P) then
      nbad = nbad + 1
      print '(A,I0,A,2L2,A,L2,A,3ES24.16)', '   failure in body case ', ncase, ': contact, ref ', h, mr, &
                                            ', expected ', mexp, ', d, ref, expected ', d, sqrt(r2), dexp
   endif
   endsubroutine compare_one

   subroutine compare_2d(s, t, mexp, dexp)
   !< s%distance_to(t) and t%distance_to(s) against the construction (contact mexp, distance dexp).
   type(curve_object), intent(in) :: s, t
   logical,            intent(in) :: mexp
   real(R8P),          intent(in) :: dexp
   real(R8P)                      :: d
   logical                        :: h
   integer(I4P)                   :: o

   do o = 1, 2
      ncase = ncase + 1
      if (o == 1) then
         call s%distance_to(t, d, intersecting=h)
      else
         call t%distance_to(s, d, intersecting=h)
      endif
      if ((h .neqv. mexp) .or. d /= dexp) then
         nbad = nbad + 1
         print '(A,I0,A,L2,A,L2,A,2ES24.16)', '   failure in curve case ', ncase, ': contact ', h, ', expected ', mexp, &
                                              ', d, expected ', d, dexp
      endif
   enddo
   endsubroutine compare_2d

   subroutine point_curve(c, u, v)
   !< An open loop of two equal vertices (one zero-length segment) at (u, v).
   type(curve_object), intent(inout) :: c
   real(R8P),          intent(in)    :: u, v

   call c%adopt_loops(u=[u, u], v=[v, v], loop_sizes=[2_I4P], closed=[.false.])
   endsubroutine point_curve

   ! references -------------------------------------------------------------------------------------------------------
   function ref_d2(x, y) result(d2)
   !< Minimum of the 6 vertex-triangle and 9 edge-edge squared distances, in the library's order (a NaN never wins).
   real(R8P), intent(in) :: x(3, 3), y(3, 3)
   real(R8P)             :: d2, c(3), c2(3), dd
   integer(I4P)          :: i, j

   d2 = huge(1._R8P)
   do i = 1, 3
      c = closest_on_triangle(x(:, i), y) ; dd = sum((x(:, i) - c)**2) ; if (dd < d2) d2 = dd
   enddo
   do i = 1, 3
      c = closest_on_triangle(y(:, i), x) ; dd = sum((y(:, i) - c)**2) ; if (dd < d2) d2 = dd
   enddo
   do i = 1, 3
      do j = 1, 3
         call closest_segment_segment(x(:, i), x(:, mod(i, 3) + 1), y(:, j), y(:, mod(j, 3) + 1), c, c2)
         dd = sum((c - c2)**2) ; if (dd < d2) d2 = dd
      enddo
   enddo
   endfunction ref_d2

   function ref_meet(x, y) result(meet)
   !< Exact contact of two closed triangles, each regular, a segment or a point (by the exact cross product).
   real(R8P), intent(in) :: x(3, 3), y(3, 3)
   logical               :: meet
   integer(I4P)          :: kx, ky
   real(R8P)             :: x0(3), x1(3), y0(3), y1(3)

   kx = shape_of(x, x0, x1) ; ky = shape_of(y, y0, y1)
   if (kx == 2 .and. ky == 2) then ! both regular: only far apart here (boxes disjoint), else no reference
      if (all(maxval(x, dim=2) >= minval(y, dim=2)) .and. all(maxval(y, dim=2) >= minval(x, dim=2))) &
         error stop 'ref_meet: two regular triangles with overlapping boxes'
      meet = .false.
   elseif (kx == 2) then
      meet = ref_seg_tri(y0, y1, ky, x)
   elseif (ky == 2) then
      meet = ref_seg_tri(x0, x1, kx, y)
   elseif (kx == 0 .and. ky == 0) then
      meet = all(x0 == y0)
   elseif (kx == 0) then
      meet = ref_on_seg(x0, y0, y1)
   elseif (ky == 0) then
      meet = ref_on_seg(y0, x0, x1)
   else
      meet = ref_seg_seg(x0, x1, y0, y1)
   endif
   endfunction ref_meet

   function shape_of(t, e0, e1) result(k)
   !< 2 regular (a nonzero exact cross product component), 1 segment (e0, e1: the two farthest vertices on the dominant
   !< axis), 0 point (e0).
   real(R8P), intent(in)  :: t(3, 3)
   real(R8P), intent(out) :: e0(3), e1(3)
   integer(I4P)           :: k, ax, i, lo, hi

   e0 = t(:, 1) ; e1 = t(:, 1)
   if (o2(t(:, 1), t(:, 2), t(:, 3), 1) /= 0 .or. o2(t(:, 1), t(:, 2), t(:, 3), 2) /= 0 .or. &
       o2(t(:, 1), t(:, 2), t(:, 3), 3) /= 0) then
      k = 2 ; return
   endif
   if (all(t(:, 1) == t(:, 2)) .and. all(t(:, 1) == t(:, 3))) then
      k = 0 ; return
   endif
   k = 1
   ax = maxloc(maxval(t, dim=2) - minval(t, dim=2), dim=1)
   lo = 1 ; hi = 1
   do i = 2, 3
      if (t(ax, i) < t(ax, lo)) lo = i
      if (t(ax, i) > t(ax, hi)) hi = i
   enddo
   e0 = t(:, lo) ; e1 = t(:, hi)
   endfunction shape_of

   function ref_seg_tri(p, q, k, t) result(meet)
   !< Segment (p, q) (a point if k = 0) against the regular triangle t.
   real(R8P),    intent(in) :: p(3), q(3), t(3, 3)
   integer(I4P), intent(in) :: k
   logical                  :: meet
   integer(I4P)             :: sp, sq, s(3), i

   sp = sgn(orient3d_xyz(t(:, 1), t(:, 2), t(:, 3), p)) ; sq = sgn(orient3d_xyz(t(:, 1), t(:, 2), t(:, 3), q))
   if (sp * sq > 0) then
      meet = .false.
   elseif (sp == 0 .and. sq == 0) then
      meet = ref_in_tri(p, t) .or. ref_in_tri(q, t)
      if (k == 1) then
         do i = 1, 3
            meet = meet .or. ref_seg_seg(p, q, t(:, i), t(:, mod(i, 3) + 1))
         enddo
      endif
   elseif (sp == 0) then
      meet = ref_in_tri(p, t)
   elseif (sq == 0) then
      meet = ref_in_tri(q, t)
   else
      do i = 1, 3
         s(i) = sgn(orient3d_xyz(p, q, t(:, i), t(:, mod(i, 3) + 1)))
      enddo
      meet = all(s >= 0) .or. all(s <= 0)
   endif
   endfunction ref_seg_tri

   function ref_in_tri(x, t) result(yes)
   !< x (on the plane of t) in the closed t: inside in every coordinate projection where t is not degenerate.
   real(R8P), intent(in) :: x(3), t(3, 3)
   logical               :: yes
   integer(I4P)          :: d, s1, s2, s3

   yes = .true.
   do d = 1, 3
      if (o2(t(:, 1), t(:, 2), t(:, 3), d) == 0) cycle
      s1 = o2(t(:, 1), t(:, 2), x, d) ; s2 = o2(t(:, 2), t(:, 3), x, d) ; s3 = o2(t(:, 3), t(:, 1), x, d)
      yes = yes .and. ((s1 >= 0 .and. s2 >= 0 .and. s3 >= 0) .or. (s1 <= 0 .and. s2 <= 0 .and. s3 <= 0))
   enddo
   endfunction ref_in_tri

   function ref_on_seg(x, p, q) result(yes)
   !< x on the closed segment (p, q), p /= q: collinear in every projection, between on the dominant axis.
   real(R8P), intent(in) :: x(3), p(3), q(3)
   logical               :: yes
   integer(I4P)          :: ax

   ax = maxloc(abs(q - p), dim=1)
   yes = o2(p, q, x, 1) == 0 .and. o2(p, q, x, 2) == 0 .and. o2(p, q, x, 3) == 0 .and. &
         x(ax) >= min(p(ax), q(ax)) .and. x(ax) <= max(p(ax), q(ax))
   endfunction ref_on_seg

   function ref_seg_seg(p, q, r, s) result(yes)
   !< Closed segments (p, q) and (r, s), both of positive length.
   real(R8P), intent(in) :: p(3), q(3), r(3), s(3)
   logical               :: yes
   integer(I4P)          :: d, d1, d2, d3, d4

   yes = ref_on_seg(p, r, s) .or. ref_on_seg(q, r, s) .or. ref_on_seg(r, p, q) .or. ref_on_seg(s, p, q)
   if (yes .or. orient3d_xyz(p, q, r, s) /= 0._R8P) return
   do d = 1, 3 ! a projection in which the common plane stays a plane: proper crossing there
      if (o2(p, q, r, d) == 0 .and. o2(p, q, s, d) == 0) cycle
      d1 = o2(r, s, p, d) ; d2 = o2(r, s, q, d) ; d3 = o2(p, q, r, d) ; d4 = o2(p, q, s, d)
      yes = d1 * d2 < 0 .and. d3 * d4 < 0
      return
   enddo
   endfunction ref_seg_seg

   function o2(p, q, r, d) result(s)
   !< Exact orientation sign of p, q, r projected along axis d.
   real(R8P),    intent(in) :: p(3), q(3), r(3)
   integer(I4P), intent(in) :: d
   integer(I4P)             :: s
   integer(I4P)             :: u, v

   u = mod(d, 3) + 1 ; v = mod(d + 1, 3) + 1
   s = sgn(orient2d_xy(p(u), p(v), q(u), q(v), r(u), r(v)))
   endfunction o2

   pure function sgn(x) result(s)
   real(R8P), intent(in) :: x
   integer(I4P)          :: s

   s = 0 ; if (x > 0._R8P) s = 1 ; if (x < 0._R8P) s = -1
   endfunction sgn

   ! builders ---------------------------------------------------------------------------------------------------------
   pure function tri(p, q, r) result(t)
   real(R8P), intent(in) :: p(3), q(3), r(3)
   real(R8P)             :: t(3, 3)

   t(:, 1) = p ; t(:, 2) = q ; t(:, 3) = r
   endfunction tri

   pure function pnt(p) result(t)
   real(R8P), intent(in) :: p(3)
   real(R8P)             :: t(3, 3)

   t(:, 1) = p ; t(:, 2) = p ; t(:, 3) = p
   endfunction pnt

   subroutine body(s, t, far)
   !< A body of two facets: t (degenerate) and a regular triangle far away.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: t(3, 3), far(3, 3)
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   integer(I4P)                            :: k

   s = fresh
   allocate(f(2))
   do k = 1, 3
      f(1)%vertex(k) = vector_R8P(t(1, k), t(2, k), t(3, k)) ; f(2)%vertex(k) = vector_R8P(far(1, k), far(2, k), far(3, k))
   enddo
   call f(1)%compute_normal ; call f(2)%compute_normal
   call s%adopt_facets(f)
   endsubroutine body
endprogram fossil_test_body_distance_degenerate
