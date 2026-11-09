!< FOSSIL test: booleans and resolve_self_intersections with an exactly degenerate (zero-area) facet (issues #98, #99).
!<
!< The body: a tetrahedron whose edge (a, b) is split at its midpoint m in one of the two faces through it, the
!< T-junction closed by the zero-area sliver (a, b, m) (watertight, manifold: a common exporter artefact). The control
!< is the same tetrahedron unsplit. The other body is the unit cube (src/tests/cube.stl). The policy (#99 step 4):
!< a degenerate facet that does not meet the other body passes through the arrangement uncut (the pair is skipped,
!< exactly, by `degenerate_pair_meet`); one that touches or crosses it gives BOOL_STATUS_INVALID_INPUT, the surface
!< unchanged (`sanitize` the inputs first).
!< 1. #98's table: the tetrahedron 0.354 from the cube (x + y <= -0.5), the sliver's box overlapping the boxes of the
!<    cube facets on x = 0 and y = 0. UNION, INTERSECT, DIFFERENCE, both operand orders: status OK and the result is,
!<    bit for bit, the facets of the bodies it must contain (UNION both, INTERSECT none, DIFFERENCE the first); volume
!<    equal to the control's (round-off: the split face sums differently, the difference is printed). Base 1147c33:
!<    the orders with the sliver body first fail (BOOL_STATUS_CDT_FAILED).
!< 2. The tetrahedron crossing the cube with its vertex d (the sliver on the outside, its box overlapping cube facets,
!<    the split facets cut): the three ops, both orders: status OK, the result is the result of the same body without
!<    the sliver (the split left as a T-junction), plus the sliver where the op keeps the outside of the body (U, body
!<    minus cube), bit for bit, watertight; volume equal to the control's (round-off, printed). Base 1147c33: the
!<    orders with the sliver body first fail (BOOL_STATUS_CDT_FAILED).
!< 3. The sliver touching the cube at one point (a on the face x = 0) and crossing it (a inside): the three ops, both
!<    orders: BOOL_STATUS_INVALID_INPUT and the surface unchanged (the controls give OK). Base 1147c33: touching gives
!<    OK cube first (the sliver passed through uncut) and CDT_FAILED body first;
!<    crossing gives CDT_FAILED.
!< 4. `resolve_self_intersections`: two crossing cubes and the tetrahedron of 1 in one surface. The facets on either
!<    side of the T-junction share a segment of the split edge, so `find_self_intersections` reports them (a crossing
!<    along an edge) and the self-arrangement splits the facet across the T-junction at the sliver's middle vertex:
!<    that closes the T-junction and the sliver is dropped (kept, it would be a fin on a non-manifold edge: base +
!<    the pass-through alone gave BOOL_STATUS_NOT_CLOSED). OK, watertight, manifold, and equal bit for bit to the
!<    result of the same surface without the sliver (the T-junction open); volume equal to the control's (round-off,
!<    printed). The cube and the tetrahedron of 3 (crossing) in one surface: BOOL_STATUS_INVALID_INPUT, unchanged
!<    (the control resolves OK). Base 1147c33: BOOL_STATUS_CDT_FAILED for both.
program fossil_test_boolean_zero_area

use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                   BOOL_STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

integer(I4P), parameter  :: OPS(3) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE]
character(1), parameter  :: OPN(3) = ['U', 'I', 'D']
type(surface_stl_object) :: cube, cube2, sl, ctl, tj, r, rc, rt, soup, soup_c, cubes
real(R8P)                :: a(3), b(3), c(3), d(3)
integer(I4P)             :: k, o, st, stc, stt
logical                  :: ok(4), good

call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)

! 1. #98's table: apart
a = [-1._R8P, 0.5_R8P, 0.5_R8P] ; b = [0._R8P, -0.5_R8P, 0.5_R8P]
c = [-1._R8P, -0.5_R8P, 0._R8P] ; d = [-1._R8P, -0.5_R8P, 1._R8P]
call tetra(sl, a, b, c, d, 3_I4P, 1_I4P) ; call tetra(ctl, a, b, c, d, 3_I4P, 0_I4P)
ok(1) = sl%is_watertight() .and. sl%get_facets_number() == 6
do k = 1, 3
   do o = 1, 2
      call run(o, k, sl, ctl, r, rc, st, stc)
      good = st == BOOL_STATUS_OK .and. stc == BOOL_STATUS_OK
      select case (k)
      case (1)
         good = good .and. same_facets(r, [facets(cube), facets(sl)])
      case (2)
         good = good .and. r%get_facets_number() == 0
      case (3)
         if (o == 1) then
            good = good .and. same_facets(r, facets(cube))
         else
            good = good .and. same_facets(r, facets(sl))
         endif
      endselect
      call report('1. apart          ', o, k, st, stc, r, rc, good)
      ok(1) = ok(1) .and. good
   enddo
enddo
print '(A,L1)', '1. #98 table (apart, the sliver box overlapping): ', ok(1)

! 2. crossing the cube away from the sliver
d = [0.5_R8P, 0.25_R8P, 0.625_R8P]
call tetra(sl, a, b, c, d, 4_I4P, 1_I4P) ; call tetra(ctl, a, b, c, d, 4_I4P, 0_I4P)
call tetra(tj, a, b, c, d, 4_I4P, 2_I4P)
ok(2) = sl%is_watertight()
do k = 1, 3
   do o = 1, 2
      call run(o, k, sl, ctl, r, rc, st, stc)
      call run(o, k, tj, ctl, rt, rc, stt, stc)
      good = st == BOOL_STATUS_OK .and. stc == BOOL_STATUS_OK .and. stt == BOOL_STATUS_OK
      if (good) then ! the sliver is outside the cube: kept by U and by body \ cube
         if (k == 1 .or. (k == 3 .and. o == 2)) then
            good = same_facets(r, [facets(rt), sliver_of(sl)])
         else
            good = same_facets(r, facets(rt))
         endif
         good = good .and. r%is_watertight() .and. rc%is_watertight() .and. r%get_facets_number() > 0
      endif
      call report('2. crossing       ', o, k, st, stc, r, rc, good)
      ok(2) = ok(2) .and. good
   enddo
enddo
print '(A,L1)', '2. crossing the cube away from the sliver: ', ok(2)

! 3. the sliver touching the cube at a point, crossing it
ok(3) = .true.
a = [0._R8P, 0.5_R8P, 0.5_R8P] ; b = [-1._R8P, -0.5_R8P, 0.5_R8P]
c = [-1._R8P, 0.5_R8P, 0._R8P] ; d = [-1._R8P, 0.5_R8P, 1._R8P]
call tetra(sl, a, b, c, d, 3_I4P, 1_I4P) ; call tetra(ctl, a, b, c, d, 3_I4P, 0_I4P)
call invalid('3. touching       ', ok(3))
a = [0.25_R8P, 0.5_R8P, 0.5_R8P]
call tetra(sl, a, b, c, d, 3_I4P, 1_I4P) ; call tetra(ctl, a, b, c, d, 3_I4P, 0_I4P)
call invalid('3. crossing       ', ok(3))
print '(A,L1)', '3. the sliver meeting the cube: BOOL_STATUS_INVALID_INPUT, unchanged: ', ok(3)

! 4. resolve_self_intersections
cube2 = cube ; call cube2%translate(delta=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P)) ; call cube2%analyze
cubes = cube ; call cubes%merge_solids(other=cube2) ; call cubes%analyze
a = [-1._R8P, 0.5_R8P, 0.5_R8P] ; b = [0._R8P, -0.5_R8P, 0.5_R8P]
c = [-1._R8P, -0.5_R8P, 0._R8P] ; d = [-1._R8P, -0.5_R8P, 1._R8P]
call tetra(sl, a, b, c, d, 3_I4P, 1_I4P)
soup = cubes ; call soup%merge_solids(other=sl) ; call soup%analyze
r = soup ; call r%resolve_self_intersections(status=st)
call tetra(tj, a, b, c, d, 3_I4P, 2_I4P) ; call tetra(ctl, a, b, c, d, 3_I4P, 0_I4P)
soup_c = cubes ; call soup_c%merge_solids(other=tj) ; call soup_c%analyze
rt = soup_c ; call rt%resolve_self_intersections(status=stt)
soup_c = cubes ; call soup_c%merge_solids(other=ctl) ; call soup_c%analyze
rc = soup_c ; call rc%resolve_self_intersections(status=stc)
ok(4) = st == BOOL_STATUS_OK .and. stc == BOOL_STATUS_OK .and. stt == BOOL_STATUS_OK
if (ok(4)) ok(4) = same_facets(r, facets(rt)) .and. r%is_watertight() .and. r%is_manifold()
print '(A,I0,A,I0,A,F18.15,A,ES10.2,A,L1)', '4. resolve, far sliver: status ', st, ', facets ', &
   r%get_facets_number(), ', volume ', r%get_volume(), ', - control ', r%get_volume() - rc%get_volume(), &
   ', = the T-junction body resolved, bit for bit: ', ok(4)
a = [0.25_R8P, 0.5_R8P, 0.5_R8P] ; b = [-1._R8P, -0.5_R8P, 0.5_R8P]
c = [-1._R8P, 0.5_R8P, 0._R8P] ; d = [-1._R8P, 0.5_R8P, 1._R8P]
call tetra(sl, a, b, c, d, 3_I4P, 1_I4P) ; call tetra(ctl, a, b, c, d, 3_I4P, 0_I4P)
soup = cube ; call soup%merge_solids(other=sl) ; call soup%analyze
soup_c = cube ; call soup_c%merge_solids(other=ctl) ; call soup_c%analyze
r = soup ; call r%resolve_self_intersections(status=st)
rc = soup_c ; call rc%resolve_self_intersections(status=stc)
good = st == BOOL_STATUS_INVALID_INPUT .and. stc == BOOL_STATUS_OK .and. same_facets(r, facets(soup))
print '(A,I0,A,I0,A,L1)', '4. resolve, crossing sliver: status ', st, ', control ', stc, ', unchanged: ', good
ok(4) = ok(4) .and. good
print '(A,L1)', '4. resolve_self_intersections: ', ok(4)
print '(A,L1)', 'Are all tests passed? ', all(ok)

contains
   subroutine tetra(s, a, b, c, d, w, mode)
   !< The tetrahedron (a, b, c, d), outward. mode 0: unsplit; 1: the face (a, b, w-th vertex) split at the midpoint m of
   !< (a, b), the sliver (a, b, m) closing the T-junction (watertight); 2: split, no sliver (a T-junction).
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: a(3), b(3), c(3), d(3)
   integer(I4P),             intent(in)    :: w, mode
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:), h(:)
   real(R8P)                               :: v(3, 4), g(3), m(3), t(3, 3)
   integer(I4P), parameter                 :: FACE(3, 4) = reshape([1, 2, 3, 1, 2, 4, 1, 3, 4, 2, 3, 4], [3, 4])
   integer(I4P)                            :: i, n

   s = fresh
   v(:, 1) = a ; v(:, 2) = b ; v(:, 3) = c ; v(:, 4) = d
   g = 0.25_R8P * (a + b + c + d) ; m = 0.5_R8P * (a + b)
   allocate(f(6)) ; n = 0
   do i = 1, 4
      t = v(:, FACE(:, i))
      if (dot_product(cross3(t(:, 2) - t(:, 1), t(:, 3) - t(:, 1)), g - t(:, 1)) > 0._R8P) &
         t(:, [1, 2]) = t(:, [2, 1])
      if (mode > 0 .and. FACE(3, i) == w .and. FACE(1, i) == 1 .and. FACE(2, i) == 2) then
         call add(f, n, t(:, 1), m, t(:, 3))
         call add(f, n, m, t(:, 2), t(:, 3))
         if (mode == 1) call add(f, n, t(:, 1), t(:, 2), m)
      else
         call add(f, n, t(:, 1), t(:, 2), t(:, 3))
      endif
   enddo
   allocate(h(n)) ; h = f(1:n)
   call s%adopt_facets(h)
   endsubroutine tetra

   subroutine add(f, n, p, q, t)
   !< Append the facet (p, q, t).
   type(facet_object), intent(inout) :: f(:)
   integer(I4P),       intent(inout) :: n
   real(R8P),          intent(in)    :: p(3), q(3), t(3)

   n = n + 1
   f(n)%vertex(1) = vector_R8P(p(1), p(2), p(3))
   f(n)%vertex(2) = vector_R8P(q(1), q(2), q(3))
   f(n)%vertex(3) = vector_R8P(t(1), t(2), t(3))
   call f(n)%compute_metrix
   endsubroutine add

   pure function cross3(x, y) result(z)
   real(R8P), intent(in) :: x(3), y(3)
   real(R8P)             :: z(3)

   z = [x(2) * y(3) - x(3) * y(2), x(3) * y(1) - x(1) * y(3), x(1) * y(2) - x(2) * y(1)]
   endfunction cross3

   subroutine run(o, k, s, sc, r, rc, st, stc)
   !< Op k between the cube and s (o = 1: cube first, 2: s first), and the same with the control sc.
   integer(I4P),             intent(in)    :: o, k
   type(surface_stl_object), intent(in)    :: s, sc
   type(surface_stl_object), intent(inout) :: r, rc
   integer(I4P),             intent(out)   :: st, stc

   if (o == 1) then
      r = cube ; call r%boolean(other=s, op=OPS(k), status=st)
      rc = cube ; call rc%boolean(other=sc, op=OPS(k), status=stc)
   else
      r = s ; call r%boolean(other=cube, op=OPS(k), status=st)
      rc = sc ; call rc%boolean(other=cube, op=OPS(k), status=stc)
   endif
   endsubroutine run

   subroutine report(lab, o, k, st, stc, r, rc, good)
   !< One line per case.
   character(*),             intent(in) :: lab
   integer(I4P),             intent(in) :: o, k, st, stc
   type(surface_stl_object), intent(in) :: r, rc
   logical,                  intent(in) :: good
   real(R8P)                            :: vr, vc

   vr = 0._R8P ; vc = 0._R8P
   if (st == BOOL_STATUS_OK .and. r%get_facets_number() > 0) vr = r%get_volume()
   if (stc == BOOL_STATUS_OK .and. rc%get_facets_number() > 0) vc = rc%get_volume()
   print '(A,A,A,A,I0,A,I0,A,I0,A,F18.15,A,ES10.2,A,L1)', lab, OPN(k), merge(' cube first ', ' body first ', o == 1), &
      'status ', st, ' (control ', stc, '), facets ', r%get_facets_number(), ', volume ', vr, ', - control ', &
      vr - vc, ': ', good
   endsubroutine report

   subroutine invalid(lab, ok)
   !< The three ops, both orders, between the cube and sl: INVALID_INPUT and both surfaces unchanged; the control OK.
   character(*), intent(in)    :: lab
   logical,      intent(inout) :: ok
   integer(I4P)                :: k, o, st, stc
   logical                     :: good

   do k = 1, 3
      do o = 1, 2
         call run(o, k, sl, ctl, r, rc, st, stc)
         if (o == 1) then
            good = same_facets(r, facets(cube))
         else
            good = same_facets(r, facets(sl))
         endif
         good = good .and. st == BOOL_STATUS_INVALID_INPUT .and. stc == BOOL_STATUS_OK
         call report(lab, o, k, st, stc, r, rc, good)
         ok = ok .and. good
      enddo
   enddo
   endsubroutine invalid

   function facets(s) result(f)
   !< The facets of s.
   type(surface_stl_object), intent(in) :: s
   type(facet_object), allocatable      :: f(:)
   type(facet_object), pointer          :: p
   integer(I4P)                         :: i

   allocate(f(s%get_facets_number()))
   do i = 1, size(f)
      p => s%facet_at(i) ; f(i) = p
   enddo
   endfunction facets

   function sliver_of(s) result(f)
   !< The degenerate facets of s.
   type(surface_stl_object), intent(in) :: s
   type(facet_object), allocatable      :: f(:)
   type(facet_object), allocatable      :: fs(:), h(:)
   integer(I4P)                         :: i, n

   fs = facets(s)
   allocate(h(size(fs))) ; n = 0
   do i = 1, size(fs)
      if (any(abs(cross3(xyz(fs(i), 2) - xyz(fs(i), 1), xyz(fs(i), 3) - xyz(fs(i), 1))) > 0._R8P)) cycle
      n = n + 1 ; h(n) = fs(i)
   enddo
   allocate(f(n)) ; f = h(1:n)
   endfunction sliver_of

   pure function xyz(f, k) result(p)
   type(facet_object), intent(in) :: f
   integer(I4P),       intent(in) :: k
   real(R8P)                      :: p(3)

   p = [f%vertex(k)%x, f%vertex(k)%y, f%vertex(k)%z]
   endfunction xyz

   function same_facets(s, f) result(same)
   !< The facets of s are those of f (oriented vertex triples, up to rotation and order), bit for bit.
   type(surface_stl_object), intent(in) :: s
   type(facet_object),       intent(in) :: f(:)
   logical                              :: same
   type(facet_object), allocatable      :: g(:)
   real(R8P)                            :: kg(9), kf(9)
   logical, allocatable                 :: used(:)
   integer(I4P)                         :: i, j

   g = facets(s)
   same = size(g) == size(f)
   if (.not. same) return
   allocate(used(size(f))) ; used = .false.
   do i = 1, size(g)
      kg = key(g(i))
      same = .false.
      do j = 1, size(f)
         if (used(j)) cycle
         kf = key(f(j))
         if (all(kg == kf)) then
            used(j) = .true. ; same = .true. ; exit
         endif
      enddo
      if (.not. same) return
   enddo
   endfunction same_facets

   pure function key(f) result(k)
   !< The vertices of f rotated to start at the lexicographically smallest one (orientation kept).
   type(facet_object), intent(in) :: f
   real(R8P)                      :: k(9)
   real(R8P)                      :: v(3, 3)
   integer(I4P)                   :: i, s

   do i = 1, 3
      v(:, i) = xyz(f, i)
   enddo
   s = 1
   do i = 2, 3
      if (lex_less(v(:, i), v(:, s))) s = i
   enddo
   k = [v(:, s), v(:, mod(s, 3) + 1), v(:, mod(s + 1, 3) + 1)]
   endfunction key

   pure function lex_less(p, q) result(less)
   real(R8P), intent(in) :: p(3), q(3)
   logical               :: less
   integer(I4P)          :: i

   less = .false.
   do i = 1, 3
      if (p(i) /= q(i)) then
         less = p(i) < q(i) ; return
      endif
   enddo
   endfunction lex_less
endprogram fossil_test_boolean_zero_area
