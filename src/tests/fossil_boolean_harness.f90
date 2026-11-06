!< FOSSIL test-side module: identities every boolean result must satisfy (issue #63, roadmap #58 WP5).
module fossil_boolean_harness
!< FOSSIL test-side module: identities every boolean result must satisfy (issue #63, roadmap #58 WP5).
!<
!< `check_boolean(a, b, rep, ...)` runs the four ops of `surface%boolean` on copies of a and checks, with no analytic
!< answer:
!< - volume: V(U) + V(I) = V(A) + V(B), V(D) = V(A) - V(I), V(S) = V(U) - V(I);
!< - moments: the same identities on the raw second moments C = int x x^T about a fixed point o (the centre of the
!<   joint box), so on the inertia tensor I = tr(C) Id - C: where the volume is, not only how much; signed tetrahedra
!<   from o (`accumulate_tetrahedron`), valid on non-manifold results too;
!< - area (only when `area_identity`, i.e. no coplanar overlaps): area(U) + area(I) = area(A) + area(B);
!< - surface: every result lies on the inputs: the deviation lower bound (deterministic samples: vertices, edge
!<   midpoints, sub-triangle centroids) of the result from A and B merged is ~0; the certified upper bound is reported
!<   (capped work: a sub-triangle whose corners sit on different input facets does not settle at once);
!< - membership: at `npoints` fixed-seed points of the joint box farther than 1e-6 of its diagonal from A and B,
!<   inside(result) = op(inside(A), inside(B)), all by the solid-angle sign;
!< - validity (U, I, D of a configuration without contact, non-empty): watertight, manifold, V > 0, no
!<   self-intersecting facet pair.
!< Pass criterion (issue #33's): a result with status OK satisfies every identity that applies; a non-OK status is
!< allowed and counted. Residuals are relative to V(A) + V(B), |C(A)| + |C(B)|, area(A) + area(B), the joint diagonal.
use fossil,                 only : surface_stl_object, facet_object, intersection_pair_t, BOOL_UNION, BOOL_INTERSECT, &
                                   BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK, SIGN_SOLID_ANGLE
use fossil_mass_properties, only : moments_accumulator_t, accumulate_tetrahedron
use penf,                   only : I4P, I8P, R8P
use vecfor,                 only : vector_R8P

implicit none
private
public :: boolean_report_t, check_boolean, boolean_pass, print_report, moments, OPS, OPN, TOL_IDENTITY

integer(I4P), parameter :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF] !< The four ops.
character(1), parameter :: OPN(4) = ['U', 'I', 'D', 'S']                                         !< Their letters.
real(R8P),    parameter :: TOL_IDENTITY = 1e-9_R8P                                               !< Relative tolerance.

type :: boolean_report_t
   !< Outcome of the four ops on one pair.
   integer(I4P) :: status(4) = -1           !< Boolean status per op.
   integer(I4P) :: facets(4) = 0            !< Result facets per op.
   real(R8P)    :: volume(4) = 0._R8P       !< Result volume per op (signed tetrahedra).
   real(R8P)    :: va = 0._R8P, vb = 0._R8P !< Input volumes.
   real(R8P)    :: vol_res(3) = -1._R8P     !< Volume identity residuals: U+I, D, S (-1: not applicable).
   real(R8P)    :: mom_res(3) = -1._R8P     !< Moment identity residuals: U+I, D, S.
   real(R8P)    :: area_res = -1._R8P       !< Area identity residual.
   real(R8P)    :: dev_lower(4) = -1._R8P   !< Sampled deviation of each result from A u B (relative).
   real(R8P)    :: dev_upper(4) = -1._R8P   !< Certified deviation (relative), capped work.
   integer(I4P) :: member_bad(4) = 0        !< Membership mismatches per op.
   integer(I4P) :: member_n = 0             !< Membership points used.
   logical      :: validity_checked(4) = .false. !< Validity required for this op.
   logical      :: watertight(4) = .false., manifold(4) = .false. !< Topology per op.
   integer(I4P) :: self_pairs(4) = 0        !< Self-intersecting pairs per op (validity ops only).
endtype boolean_report_t

contains
   subroutine check_boolean(a, b, rep, npoints, seed, area_identity, contact, corrupt)
   !< Run the four ops of a against b and fill the report. `corrupt` (negative control only): drop the first facet of
   !< the result of that op before checking it.
   type(surface_stl_object), intent(in)  :: a, b          !< Inputs (closed, outward).
   type(boolean_report_t),   intent(out) :: rep           !< Report.
   integer(I4P),             intent(in)  :: npoints       !< Membership points to draw.
   integer(I8P),             intent(in)  :: seed          !< Park-Miller seed.
   logical,                  intent(in)  :: area_identity !< Check the area identity (no coplanar overlaps).
   logical,                  intent(in)  :: contact       !< The inputs touch (validity not required).
   integer(I4P), optional,   intent(in)  :: corrupt       !< Op whose result is damaged (negative control).
   type(surface_stl_object)              :: r, ab         !< Result, A and B merged.
   type(vector_R8P)                      :: bmin, bmax, p !< Joint box, a point.
   real(R8P)                             :: o(3), diag, ca(6), cb(6), cr(6, 4), ar(4), area_ab, dl, du !< Work.
   real(R8P),    allocatable             :: pts(:,:)      !< Membership points.
   logical,      allocatable             :: ina(:), inb(:) !< Input membership.
   type(intersection_pair_t), allocatable :: pairs(:)     !< Self-intersections.
   integer(I8P)                          :: s             !< Generator state.
   integer(I4P)                          :: k, i, n, st   !< Counters, status.
   logical                               :: want, inr     !< Op membership, result membership.

   call joint_box(a, b, bmin, bmax)
   o = 0.5_R8P * [bmin%x + bmax%x, bmin%y + bmax%y, bmin%z + bmax%z]
   diag = norm2([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
   call moments(a, o, rep%va, ca) ; call moments(b, o, rep%vb, cb)
   call merge_surfaces(a, b, ab)
   area_ab = a%get_area() + b%get_area()
   ! membership points, away from both surfaces
   allocate(pts(3, npoints), ina(npoints), inb(npoints))
   s = seed ; n = 0
   do i = 1, npoints
      p = vector_R8P(bmin%x + (bmax%x - bmin%x) * rnd(), bmin%y + (bmax%y - bmin%y) * rnd(), &
                     bmin%z + (bmax%z - bmin%z) * rnd())
      if (a%distance(point=p, is_signed=.false., is_square_root=.true.) < 1e-6_R8P * diag) cycle
      if (b%distance(point=p, is_signed=.false., is_square_root=.true.) < 1e-6_R8P * diag) cycle
      n = n + 1 ; pts(:, n) = [p%x, p%y, p%z]
      ina(n) = a%is_point_inside(point=p, sign_algorithm=SIGN_SOLID_ANGLE)
      inb(n) = b%is_point_inside(point=p, sign_algorithm=SIGN_SOLID_ANGLE)
   enddo
   rep%member_n = n
   do k = 1, 4
      r = a
      call r%boolean(other=b, op=OPS(k), status=st)
      rep%status(k) = st
      if (st /= BOOL_STATUS_OK) cycle
      if (present(corrupt)) then
         if (corrupt == k) call drop_first_facet(r)
      endif
      rep%facets(k) = r%get_facets_number()
      call moments(r, o, rep%volume(k), cr(:, k))
      ar(k) = r%get_area()
      do i = 1, n
         p = vector_R8P(pts(1, i), pts(2, i), pts(3, i))
         select case (k)
         case (1) ; want = ina(i) .or. inb(i)
         case (2) ; want = ina(i) .and. inb(i)
         case (3) ; want = ina(i) .and. .not. inb(i)
         case default ; want = ina(i) .neqv. inb(i)
         endselect
         inr = .false.
         if (rep%facets(k) > 0) inr = r%is_point_inside(point=p, sign_algorithm=SIGN_SOLID_ANGLE)
         if (inr .neqv. want) rep%member_bad(k) = rep%member_bad(k) + 1
      enddo
      if (rep%facets(k) > 0) then
         call r%deviation(ab, lower=dl, upper=du, tolerance=TOL_IDENTITY * diag, max_triangles=4 * rep%facets(k), status=st)
         rep%dev_lower(k) = dl / diag ; rep%dev_upper(k) = du / diag
      else
         rep%dev_lower(k) = 0._R8P ; rep%dev_upper(k) = 0._R8P
      endif
      rep%watertight(k) = r%is_watertight() ; rep%manifold(k) = r%is_manifold()
      rep%validity_checked(k) = k <= 3 .and. .not. contact .and. rep%facets(k) > 0
      if (rep%validity_checked(k)) then
         call r%find_self_intersections(pairs=pairs)
         rep%self_pairs(k) = size(pairs)
      endif
   enddo
   if (all(rep%status([1, 2]) == BOOL_STATUS_OK)) then
      rep%vol_res(1) = abs(rep%volume(1) + rep%volume(2) - rep%va - rep%vb) / (rep%va + rep%vb)
      rep%mom_res(1) = maxval(abs(cr(:, 1) + cr(:, 2) - ca - cb)) / maxval(abs(ca) + abs(cb))
      if (area_identity) rep%area_res = abs(ar(1) + ar(2) - area_ab) / area_ab
   endif
   if (all(rep%status([2, 3]) == BOOL_STATUS_OK)) then
      rep%vol_res(2) = abs(rep%volume(3) - rep%va + rep%volume(2)) / (rep%va + rep%vb)
      rep%mom_res(2) = maxval(abs(cr(:, 3) - ca + cr(:, 2))) / maxval(abs(ca) + abs(cb))
   endif
   if (all(rep%status([1, 2, 4]) == BOOL_STATUS_OK)) then
      rep%vol_res(3) = abs(rep%volume(4) - rep%volume(1) + rep%volume(2)) / (rep%va + rep%vb)
      rep%mom_res(3) = maxval(abs(cr(:, 4) - cr(:, 1) + cr(:, 2))) / maxval(abs(ca) + abs(cb))
   endif
   contains
      function rnd() result(x)
      !< Park-Miller in (0, 1).
      real(R8P) :: x

      s = mod(s * 48271_I8P, 2147483647_I8P)
      x = real(s, R8P) / 2147483647._R8P
      endfunction rnd
   endsubroutine check_boolean

   pure function boolean_pass(rep) result(ok)
   !< The pass criterion: every identity that applies holds (a non-OK status is allowed).
   type(boolean_report_t), intent(in) :: rep !< Report.
   logical                            :: ok  !< Pass.
   integer(I4P)                       :: k   !< Counter.

   ok = all(rep%vol_res <= TOL_IDENTITY) .and. all(rep%mom_res <= TOL_IDENTITY) .and. rep%area_res <= TOL_IDENTITY
   do k = 1, 4
      if (rep%status(k) /= BOOL_STATUS_OK) cycle
      ok = ok .and. rep%member_bad(k) == 0 .and. rep%dev_lower(k) <= TOL_IDENTITY
      if (rep%validity_checked(k)) ok = ok .and. rep%watertight(k) .and. rep%manifold(k) .and. &
                                         rep%volume(k) > 0._R8P .and. rep%self_pairs(k) == 0
   enddo
   endfunction boolean_pass

   subroutine print_report(label, rep)
   !< One line per pair: statuses, worst residuals, membership and validity.
   character(*),           intent(in) :: label !< Case name.
   type(boolean_report_t), intent(in) :: rep   !< Report.
   character(4)                       :: sts   !< Status letters.
   integer(I4P)                       :: k     !< Counter.

   do k = 1, 4
      sts(k:k) = merge(OPN(k), 'x', rep%status(k) == BOOL_STATUS_OK)
   enddo
   print '(A,1X,A4,A,ES8.1,A,ES8.1,A,ES8.1,A,ES8.1,A,I0,A,I0,A,L1)', label, sts, ' vol ', maxval(rep%vol_res), &
         ' mom ', maxval(rep%mom_res), ' area ', rep%area_res, ' dev ', maxval(rep%dev_lower), ' memb ', &
         sum(rep%member_bad), '/', rep%member_n, ' pass ', boolean_pass(rep)
   endsubroutine print_report

   subroutine joint_box(a, b, bmin, bmax)
   !< Bounding box of a and b, inflated by 5 %.
   type(surface_stl_object), intent(in)  :: a, b       !< Inputs.
   type(vector_R8P),         intent(out) :: bmin, bmax !< Box.
   type(vector_R8P)                      :: l1, h1, l2, h2 !< Input boxes.
   real(R8P)                             :: lo(3), hi(3), m(3) !< Work.

   l1 = a%get_bmin() ; h1 = a%get_bmax() ; l2 = b%get_bmin() ; h2 = b%get_bmax()
   lo = min([l1%x, l1%y, l1%z], [l2%x, l2%y, l2%z]) ; hi = max([h1%x, h1%y, h1%z], [h2%x, h2%y, h2%z])
   m = 0.05_R8P * (hi - lo) ; lo = lo - m ; hi = hi + m
   bmin = vector_R8P(lo(1), lo(2), lo(3)) ; bmax = vector_R8P(hi(1), hi(2), hi(3))
   endsubroutine joint_box

   subroutine moments(s, o, v, c)
   !< Volume and raw second moments [Cxx, Cyy, Czz, Cxy, Cxz, Cyz] about o, by signed tetrahedra.
   type(surface_stl_object), intent(in)  :: s    !< Surface.
   real(R8P),                intent(in)  :: o(3) !< Reference point.
   real(R8P),                intent(out) :: v    !< Volume.
   real(R8P),                intent(out) :: c(6) !< Second moments.
   type(moments_accumulator_t)           :: acc  !< Accumulator.
   type(facet_object), pointer           :: f    !< Facet.
   integer(I4P)                          :: i    !< Counter.

   do i = 1, s%get_facets_number()
      f => s%facet_at(i)
      call accumulate_tetrahedron(acc, [f%vertex(1)%x, f%vertex(1)%y, f%vertex(1)%z] - o, &
                                  [f%vertex(2)%x, f%vertex(2)%y, f%vertex(2)%z] - o, &
                                  [f%vertex(3)%x, f%vertex(3)%y, f%vertex(3)%z] - o)
   enddo
   v = acc%s(1) + acc%c(1) ; c = acc%s(5:10) + acc%c(5:10)
   endsubroutine moments

   subroutine drop_first_facet(r)
   !< Remove the first facet of r (a deliberately wrong result).
   type(surface_stl_object), intent(inout) :: r    !< Result.
   type(facet_object), allocatable         :: f(:) !< Facets.
   type(facet_object), pointer             :: fp   !< Facet.
   integer(I4P)                            :: i    !< Counter.

   allocate(f(r%get_facets_number() - 1))
   do i = 2, r%get_facets_number()
      fp => r%facet_at(i) ; f(i - 1) = fp
   enddo
   call r%adopt_facets(facets=f)
   call r%analyze
   endsubroutine drop_first_facet

   subroutine merge_surfaces(a, b, ab)
   !< A and B as one surface (their facets side by side).
   type(surface_stl_object), intent(in)  :: a, b !< Inputs.
   type(surface_stl_object), intent(out) :: ab   !< Merged.
   type(facet_object), allocatable       :: f(:) !< Facets.
   type(facet_object), pointer           :: fp   !< Facet.
   integer(I4P)                          :: i, na !< Counters.

   na = a%get_facets_number()
   allocate(f(na + b%get_facets_number()))
   do i = 1, na
      fp => a%facet_at(i) ; f(i) = fp
   enddo
   do i = 1, b%get_facets_number()
      fp => b%facet_at(i) ; f(na + i) = fp
   enddo
   call ab%adopt_facets(facets=f)
   call ab%analyze
   endsubroutine merge_surfaces
endmodule fossil_boolean_harness
