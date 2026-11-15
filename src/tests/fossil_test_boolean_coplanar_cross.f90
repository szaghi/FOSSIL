!< FOSSIL test: a cut along a coplanar partner's edge crossing another cut inside a facet (roadmap #92 WP5 (a)).
!<
!< Boxes with dyadic corners (every coordinate exact), each face split along the diagonal from its lowest corner, merged
!< into one self-intersecting surface and resolved (`resolve_self_intersections`). A = [0, 1]^3; B = A + (0.25, 0.125,
!< 0) shares A's planes z = 0 and z = 1, so B's bottom facets are coplanar partners of A's bottom facet F (the triangle
!< y < x of z = 0) and cut it along their edges, among them B's bottom diagonal y = x - 0.125, whose facets on both
!< sides are coplanar partners too (no transversal facet runs along it).
!< 1. C = [0.625, 1.625] x [-0.25, 0.75] x [-0.5, 0.5] crosses z = 0: its faces x = 0.625 and y = 0.75 cut F
!<    transversally and cross B's diagonal inside F at (0.625, 0.5, 0) and (0.875, 0.75, 0), and B's edge y = 0.125 at
!<    (0.625, 0.125, 0) (a partner edge with a transversal facet along it, B's face y = 0.125).
!< 2. D = [0.5, 1.5] x [-0.3125, 0.6875] x [0, 1] shares z = 0 too: its bottom edges x = 0.5, y = 0.6875 and its
!<    diagonal y = x - 0.8125 cross B's diagonal and B's edge y = 0.125 inside F (cuts of two coplanar partners).
!< 3. A, B, C and D together.
!< 4. E = [0.375, 1.375] x [-0.5, 0.5] x [0, 1] shares z = 0, its bottom split along the other diagonal x + y = 0.875
!<    (both facets on it coplanar partners of F): it crosses B's diagonal inside F at (0.5, 0.375, 0), where no
!<    transversal facet runs along either cut (the edge-edge point of the two partners' edges, IPK_EEI).
!< Each self-union must equal the boolean chain of the boxes (A u B) u ... (no self-intersecting operand there) and the
!< analytic union volume (inclusion-exclusion over the boxes): status OK, watertight, manifold, one shell, no
!< self-intersection left, volume to 1e-12 of both, area to 1e-12 of the chain's; membership: at the points of a grid of
!< step 1/8 offset by 1/32 (never on a box plane, all at multiples of 1/16) inside(self-union) = inside(some box), by
!< the solid-angle sign. Then the boolean harness (`check_boolean`: volume and moment identities, deviation, membership
!< of the four ops) runs on the self-union against the last box (in contact: it lies on the self-union's boundary).
!< Negative control (base afe0425, before the crossing points of partners' cuts): every case returns
!< BOOL_STATUS_CDT_FAILED (cdt_facet: CDTI_CROSSING in F).
program fossil_test_boolean_coplanar_cross

use fossil,                 only : surface_stl_object, facet_object, intersection_pair_t, BOOL_UNION, BOOL_STATUS_OK, &
                                   SIGN_SOLID_ANGLE
use fossil_boolean_harness, only : boolean_report_t, check_boolean, boolean_pass, print_report
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none

real(R8P), parameter :: BOX_A(6) = [0._R8P, 1._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P]
real(R8P), parameter :: BOX_B(6) = [0.25_R8P, 1.25_R8P, 0.125_R8P, 1.125_R8P, 0._R8P, 1._R8P]
real(R8P), parameter :: BOX_C(6) = [0.625_R8P, 1.625_R8P, -0.25_R8P, 0.75_R8P, -0.5_R8P, 0.5_R8P]
real(R8P), parameter :: BOX_D(6) = [0.5_R8P, 1.5_R8P, -0.3125_R8P, 0.6875_R8P, 0._R8P, 1._R8P]
real(R8P), parameter :: BOX_E(6) = [0.375_R8P, 1.375_R8P, -0.5_R8P, 0.5_R8P, 0._R8P, 1._R8P]
logical              :: ok(4)

ok(1) = check('1. partner edges x transversal cuts', reshape([BOX_A, BOX_B, BOX_C], [6, 3]))
ok(2) = check('2. edges of two partners cross     ', reshape([BOX_A, BOX_B, BOX_D], [6, 3]))
ok(3) = check('3. both                            ', reshape([BOX_A, BOX_B, BOX_C, BOX_D], [6, 4]))
ok(4) = check('4. crossing diagonals of partners ', reshape([BOX_A, BOX_B, BOX_E], [6, 3]), anti=3)
print '(A,L1)', 'Are all tests passed? ', all(ok)

contains
   function check(label, boxes, anti) result(pass)
   !< Resolve the soup of the boxes and compare it with the boolean chain and the analytic union volume.
   character(*), intent(in)               :: label
   real(R8P),    intent(in)               :: boxes(:, :)
   integer(I4P), intent(in), optional     :: anti  !< Box whose bottom is split along the other diagonal.
   logical                                :: pass
   type(surface_stl_object)               :: soup, chain, box
   type(intersection_pair_t), allocatable :: pairs(:)
   type(boolean_report_t)                 :: rep
   integer(I4P)                           :: status, k, npairs, nbad, npts, i, j, l
   real(R8P)                              :: v_exact, lo(3), hi(3), p(3)
   logical                                :: flip

   call make_box(boxes(:, 1), soup) ; call make_box(boxes(:, 1), chain)
   do k = 2, size(boxes, 2)
      flip = .false.
      if (present(anti)) flip = k == anti
      call make_box(boxes(:, k), box, anti=flip)
      call soup%merge_solids(other=box)
      call chain%boolean(other=box, op=BOOL_UNION)
   enddo
   call soup%analyze
   call soup%find_self_intersections(pairs=pairs) ; npairs = size(pairs)
   call soup%resolve_self_intersections(status=status)
   call soup%find_self_intersections(pairs=pairs)
   v_exact = union_volume(boxes)
   print '(A,A,I0,A,I0,A,I0,A,F19.16,A,L1,L1,A,I0,A,I0)', label, ': input pairs ', npairs, ', status ', status, &
      ', facets ', soup%get_facets_number(), ', volume ', soup%get_volume(), ', watertight, manifold ', &
      soup%is_watertight(), soup%is_manifold(), ', shells ', soup%get_shells_number(), ', pairs ', size(pairs)
   print '(A,F19.16,A,F19.16,A,ES10.3,A,ES10.3)', '   exact volume ', v_exact, ', chain volume ', chain%get_volume(), &
      ', |dV| ', abs(soup%get_volume() - chain%get_volume()), ', |dA| ', abs(area(soup) - area(chain))
   lo = minval(boxes([1, 3, 5], :), dim=2) ; hi = maxval(boxes([2, 4, 6], :), dim=2)
   nbad = 0 ; npts = 0
   do i = 0, nint((hi(1) - lo(1)) * 8._R8P) - 1
      do j = 0, nint((hi(2) - lo(2)) * 8._R8P) - 1
         do l = 0, nint((hi(3) - lo(3)) * 8._R8P) - 1
            p = lo + [i, j, l] / 8._R8P + 1._R8P / 32._R8P
            npts = npts + 1
            if (soup%is_point_inside(point=vector_R8P(p(1), p(2), p(3)), sign_algorithm=SIGN_SOLID_ANGLE) .neqv. &
                any(all(spread(p, 2, size(boxes, 2)) >= boxes([1, 3, 5], :) .and. &
                        spread(p, 2, size(boxes, 2)) <= boxes([2, 4, 6], :), dim=1))) nbad = nbad + 1
         enddo
      enddo
   enddo
   print '(A,I0,A,I0)', '   membership mismatches ', nbad, ' of ', npts
   pass = status == BOOL_STATUS_OK .and. soup%is_watertight() .and. soup%is_manifold() .and. &
          soup%get_shells_number() == 1 .and. size(pairs) == 0 .and. &
          abs(soup%get_volume() - v_exact) < 1.e-12_R8P .and. abs(chain%get_volume() - v_exact) < 1.e-12_R8P .and. &
          abs(area(soup) - area(chain)) < 1.e-12_R8P .and. nbad == 0
   if (status == BOOL_STATUS_OK) then
      call check_boolean(soup, box, rep, 40_I4P, 20261115_I8P, area_identity=.false., contact=.true.)
      call print_report('   harness, self-union vs last box:', rep)
      pass = pass .and. boolean_pass(rep) .and. all(rep%status == BOOL_STATUS_OK)
   endif
   print '(A,A,L1)', label, ': self-union = boolean chain = exact union: ', pass
   endfunction check

   subroutine make_box(b, s, anti)
   !< The box [b1, b2] x [b3, b4] x [b5, b6], outward facets, each face split along the diagonal from its lowest corner
   !< (the bottom along the other one with `anti`).
   real(R8P),                intent(in)           :: b(6)
   type(surface_stl_object), intent(out)          :: s
   logical,                  intent(in), optional :: anti
   type(facet_object), allocatable       :: f(:)
   real(R8P)                             :: q(3, 4, 6)
   integer(I4P)                          :: k

   ! faces -z, +z, -x, +x, -y, +y: corners counter-clockwise seen from outside, the first one the lowest
   q(:, :, 1) = reshape([b(1), b(3), b(5), b(1), b(4), b(5), b(2), b(4), b(5), b(2), b(3), b(5)], [3, 4])
   q(:, :, 2) = reshape([b(1), b(3), b(6), b(2), b(3), b(6), b(2), b(4), b(6), b(1), b(4), b(6)], [3, 4])
   q(:, :, 3) = reshape([b(1), b(3), b(5), b(1), b(3), b(6), b(1), b(4), b(6), b(1), b(4), b(5)], [3, 4])
   q(:, :, 4) = reshape([b(2), b(3), b(5), b(2), b(4), b(5), b(2), b(4), b(6), b(2), b(3), b(6)], [3, 4])
   q(:, :, 5) = reshape([b(1), b(3), b(5), b(2), b(3), b(5), b(2), b(3), b(6), b(1), b(3), b(6)], [3, 4])
   q(:, :, 6) = reshape([b(1), b(4), b(5), b(1), b(4), b(6), b(2), b(4), b(6), b(2), b(4), b(5)], [3, 4])
   allocate(f(12))
   do k = 1, 6
      call set(f(2 * k - 1), q(:, 1, k), q(:, 2, k), q(:, 3, k))
      call set(f(2 * k), q(:, 1, k), q(:, 3, k), q(:, 4, k))
   enddo
   if (present(anti)) then
      if (anti) then
         call set(f(1), q(:, 1, 1), q(:, 2, 1), q(:, 4, 1))
         call set(f(2), q(:, 2, 1), q(:, 3, 1), q(:, 4, 1))
      endif
   endif
   call s%adopt_facets(facets=f)
   endsubroutine make_box

   subroutine set(f, p1, p2, p3)
   !< Facet (p1, p2, p3) with its normal.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_normal
   endsubroutine set

   pure function union_volume(boxes) result(v)
   !< Volume of the union of the boxes by inclusion-exclusion over their subsets.
   real(R8P), intent(in) :: boxes(:, :)
   real(R8P)             :: v
   real(R8P)             :: lo(3), hi(3)
   integer(I4P)          :: mask, k, n, bits

   n = size(boxes, 2) ; v = 0._R8P
   do mask = 1, 2 ** n - 1
      lo = -huge(1._R8P) ; hi = huge(1._R8P) ; bits = 0
      do k = 1, n
         if (.not. btest(mask, k - 1)) cycle
         bits = bits + 1
         lo = max(lo, boxes([1, 3, 5], k)) ; hi = min(hi, boxes([2, 4, 6], k))
      enddo
      if (any(hi <= lo)) cycle
      v = v + (-1._R8P) ** (bits + 1) * product(hi - lo)
   enddo
   endfunction union_volume

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
endprogram fossil_test_boolean_coplanar_cross
