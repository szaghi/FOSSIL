!< FOSSIL test: exact cuts on implicit points and the facet triangulation (issue #65, steps 2-4).
!<
!< 1. `pair_cut`: a proper crossing gives two line-plane points; the same edge crossing the same plane, from the two
!<    facets that share the edge, is the same identity (and the same rounded bits); an edge of B through an edge of A
!<    gives the edge-edge point from both sides of the pair (tie: the more specific construction); a vertex of B on
!<    A's plane is that explicit vertex; a vertex one ulp off the plane is not coplanar (no snap).
!< 2. `coplanar_cut`: a chord through a triangle (edge-edge ends, exact coordinates), the shared part of a collinear
!<    edge, nothing for a segment outside.
!< 3. `cdt_facet`: a triangle with interior points, points on its edges and constraints through collinear points:
!<    every constraint piece is an edge, the triangles tile the facet (rounded area sum), all keep the facet's
!<    orientation (exact orient2d).
!< 4. The arrangement of two unit cubes, the second moved by (0.5, 0, 0) and both rotated (coplanar faces off by
!<    round-off): each surface's sub-triangles close up by registry index (every edge used once in each direction),
!<    and distinct registry points have distinct rounded coordinates.
program fossil_test_exact_cuts

use fossil,                 only : surface_stl_object
use fossil_arrangement,     only : arrangement_t, arrangement_initialize, arrangement_collect_intersections, &
                                   arrangement_retriangulate, ARR_STATUS_OK
use fossil_cdt_implicit,    only : cdt_facet, CDTI_OK
use fossil_facet_object,    only : facet_object
use fossil_implicit_points, only : implicit_point_t, ip_explicit, ip_eei, ip_same, ip_approx, ip_orient2d, &
                                   IPK_LPI, IPK_EEI, IPK_EXPLICIT
use fossil_tritri_exact,    only : pair_cut, coplanar_cut, PAIR_CROSSING, PAIR_COPLANAR
use penf,                   only : I4P, R8P
use vecfor,                 only : vector_R8P

implicit none

type(vector_R8P)                :: ta(3), tb(3), tc(3)
type(implicit_point_t)          :: p, q, p2, q2
integer(I4P)                    :: kind, kind2
logical                         :: ok, ok1, ok2, ok3, ok4, found

! ---- 1. pair_cut
ta = [v(0._R8P, 0._R8P, 0._R8P), v(2._R8P, 0._R8P, 0._R8P), v(0._R8P, 2._R8P, 0._R8P)]
tb = [v(0.5_R8P, 0.5_R8P, -1._R8P), v(0.5_R8P, 0.5_R8P, 1._R8P), v(1.5_R8P, 0.25_R8P, 0.3_R8P)]
call pair_cut(ta, tb, kind, p, q)
ok1 = kind == PAIR_CROSSING .and. p%kind == IPK_LPI .and. q%kind == IPK_LPI
! the facet sharing tb's edge (1, 2) gives the same point for that edge
tc = [tb(2), tb(1), v(0._R8P, 1.5_R8P, 0.1_R8P)]
call pair_cut(ta, tc, kind2, p2, q2)
ok1 = ok1 .and. kind2 == PAIR_CROSSING .and. (ip_same(p, p2) .or. ip_same(p, q2) .or. ip_same(q, p2) .or. ip_same(q, q2))
! B's edge (1, 2) through A's edge (1, 2) (the x axis): the edge-edge point, identical from both sides of the tie
tb = [v(0.5_R8P, 0._R8P, -1._R8P), v(0.5_R8P, 0._R8P, 1._R8P), v(0.5_R8P, 1._R8P, 0.3_R8P)]
call pair_cut(ta, tb, kind, p, q)
ok1 = ok1 .and. kind == PAIR_CROSSING .and. &
      (ip_same(p, ip_eei(xyz(ta(1)), xyz(ta(2)), xyz(tb(1)), xyz(tb(2)))) .or. &
       ip_same(q, ip_eei(xyz(ta(1)), xyz(ta(2)), xyz(tb(1)), xyz(tb(2)))))
! a vertex of B on A's plane: the explicit vertex
tb = [v(0.5_R8P, 0.5_R8P, 0._R8P), v(0.5_R8P, 0.5_R8P, 1._R8P), v(1._R8P, 0.25_R8P, -1._R8P)]
call pair_cut(ta, tb, kind, p, q)
ok1 = ok1 .and. kind == PAIR_CROSSING .and. (ip_same(p, ip_explicit(xyz(tb(1)))) .or. ip_same(q, ip_explicit(xyz(tb(1)))))
! coplanar, and one ulp off: no longer coplanar
tb = [v(0.5_R8P, 0.5_R8P, 0._R8P), v(1.5_R8P, 0.5_R8P, 0._R8P), v(0.5_R8P, 1.5_R8P, 0._R8P)]
call pair_cut(ta, tb, kind, p, q)
ok1 = ok1 .and. kind == PAIR_COPLANAR
tb(2)%z = nearest(0._R8P, 1._R8P)
call pair_cut(ta, tb, kind, p, q)
ok1 = ok1 .and. kind /= PAIR_COPLANAR
print '(A,L1)', '1. pair_cut: ', ok1

! ---- 2. coplanar_cut
call coplanar_cut(v(-1._R8P, 0.5_R8P, 0._R8P), v(3._R8P, 0.5_R8P, 0._R8P), ta, found, p, q)
ok2 = found .and. p%kind == IPK_EEI .and. q%kind == IPK_EEI .and. &
      ends(ip_approx(p), ip_approx(q), [0._R8P, 0.5_R8P, 0._R8P], [1.5_R8P, 0.5_R8P, 0._R8P])
call coplanar_cut(v(0.5_R8P, 0._R8P, 0._R8P), v(2.5_R8P, 0._R8P, 0._R8P), ta, found, p, q)
ok2 = ok2 .and. found .and. ends(ip_approx(p), ip_approx(q), [0.5_R8P, 0._R8P, 0._R8P], [2._R8P, 0._R8P, 0._R8P])
call coplanar_cut(v(3._R8P, 3._R8P, 0._R8P), v(4._R8P, 3._R8P, 0._R8P), ta, found, p, q)
ok2 = ok2 .and. .not. found
print '(A,L1)', '2. coplanar_cut: ', ok2

! ---- 3. cdt_facet
block
   type(implicit_point_t)    :: pts(9)
   integer(I4P)              :: segs(2, 3), ntri, st, k, e, i
   integer(I4P), allocatable :: tris(:, :)
   real(R8P)                 :: area, x(3, 3)
   logical                   :: has

   pts(1) = ip_explicit([0._R8P, 0._R8P, 0._R8P]) ; pts(2) = ip_explicit([4._R8P, 0._R8P, 0._R8P])
   pts(3) = ip_explicit([0._R8P, 4._R8P, 0._R8P])
   pts(4) = ip_explicit([1._R8P, 1._R8P, 0._R8P]) ; pts(5) = ip_explicit([2._R8P, 1._R8P, 0._R8P])
   pts(6) = ip_explicit([3._R8P, 1._R8P, 0._R8P]) ; pts(7) = ip_explicit([2._R8P, 0._R8P, 0._R8P])
   pts(8) = ip_explicit([0._R8P, 2._R8P, 0._R8P]) ; pts(9) = ip_explicit([1._R8P, 2._R8P, 0._R8P])
   segs(:, 1) = [4, 6]   ! through 5; 6 = (3, 1) lies on the hypotenuse
   segs(:, 2) = [7, 5]   ! from the bottom edge
   segs(:, 3) = [8, 9]   ! from the left edge
   call cdt_facet(pts, 3_I4P, segs, tris, ntri, st)
   ok3 = st == CDTI_OK
   area = 0._R8P
   do k = 1, ntri
      do i = 1, 3
         x(:, i) = ip_approx(pts(tris(i, k)))
      enddo
      area = area + 0.5_R8P * ((x(1, 2) - x(1, 1)) * (x(2, 3) - x(2, 1)) - (x(2, 2) - x(2, 1)) * (x(1, 3) - x(1, 1)))
      ok3 = ok3 .and. ip_orient2d(pts(tris(1, k)), pts(tris(2, k)), pts(tris(3, k)), 3_I4P) > 0
   enddo
   ok3 = ok3 .and. area == 8._R8P
   do e = 1, 4 ! the pieces (4, 5), (5, 6), (7, 5), (8, 9) are edges
      has = .false.
      do k = 1, ntri
         do i = 1, 3
            if (piece(e, tris(i, k), tris(mod(i, 3) + 1, k))) has = .true.
         enddo
      enddo
      ok3 = ok3 .and. has
   enddo
   print '(A,I0,A,F6.2,A,L1)', '3. cdt_facet: triangles ', ntri, ', area ', area, '  ok ', ok3
endblock

! ---- 4. arrangement of rotated near-coplanar cubes
block
   type(surface_stl_object)        :: a, b
   type(facet_object), allocatable :: fa(:), fb(:), sub(:)
   type(arrangement_t)             :: arr
   integer(I4P), allocatable       :: owner(:), src(:), pid(:, :)
   integer(I4P)                    :: st, i, j, e, f, h, nbad, nsame
   type(facet_object), pointer     :: fp

   call a%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call b%translate(x=0.5_R8P, y=0._R8P, z=0._R8P)
   call a%rotate(axis=vector_R8P(1._R8P, 2._R8P, 3._R8P), angle=0.7_R8P)
   call b%rotate(axis=vector_R8P(1._R8P, 2._R8P, 3._R8P), angle=0.7_R8P)
   call a%analyze ; call b%analyze
   allocate(fa(a%get_facets_number()), fb(b%get_facets_number()))
   do i = 1, size(fa) ; fp => a%facet_at(i) ; fa(i) = fp ; enddo
   do i = 1, size(fb) ; fp => b%facet_at(i) ; fb(i) = fp ; enddo
   call arrangement_initialize(arr, fa, a%aabb, fb, b%aabb)
   call arrangement_collect_intersections(arr)
   call arrangement_retriangulate(arr, sub, owner, src, st, pid)
   nbad = 0
   do i = 1, size(sub)
      do e = 1, 3
         h = 0
         do j = 1, size(sub)
            if (owner(j) /= owner(i) .or. j == i) cycle
            do f = 1, 3
               if (pid(f, j) == pid(mod(e, 3) + 1, i) .and. pid(mod(f, 3) + 1, j) == pid(e, i)) h = h + 1
            enddo
         enddo
         if (h /= 1) nbad = nbad + 1
      enddo
   enddo
   nsame = 0
   do i = 1, arr%n_pts
      do j = i + 1, arr%n_pts
         if (all(arr%xyz(:, i) == arr%xyz(:, j))) nsame = nsame + 1
      enddo
   enddo
   ok4 = st == ARR_STATUS_OK .and. nbad == 0 .and. nsame == 0 .and. size(sub) > size(fa) + size(fb)
   print '(A,I0,A,I0,A,I0,A,I0,A,L1)', '4. rotated cubes: status ', st, ', sub-triangles ', size(sub), &
         ', unpaired half-edges ', nbad, ', registry points sharing rounded coordinates ', nsame, '  ok ', ok4
endblock

ok = ok1 .and. ok2 .and. ok3 .and. ok4
print '(A,L1)', 'Are all tests passed? ', ok

contains
   pure function v(x, y, z) result(r)
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: r

   r%x = x ; r%y = y ; r%z = z
   endfunction v

   pure function xyz(a) result(r)
   type(vector_R8P), intent(in) :: a
   real(R8P)                    :: r(3)

   r = [a%x, a%y, a%z]
   endfunction xyz

   pure function ends(a, b, e1, e2) result(s)
   real(R8P), intent(in) :: a(3), b(3), e1(3), e2(3)
   logical               :: s

   s = (all(a == e1) .and. all(b == e2)) .or. (all(a == e2) .and. all(b == e1))
   endfunction ends

   pure function piece(e, i, j) result(s)
   integer(I4P), intent(in) :: e, i, j
   logical                  :: s
   integer(I4P), parameter  :: PC(2, 4) = reshape([4, 5, 5, 6, 7, 5, 8, 9], [2, 4])

   s = (i == PC(1, e) .and. j == PC(2, e)) .or. (i == PC(2, e) .and. j == PC(1, e))
   endfunction piece
endprogram fossil_test_exact_cuts
