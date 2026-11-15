!< FOSSIL test: `ip_compare` on coordinates that a construction fixes to an input double (issue #92 WP5 d).
!<
!< `ip_compare` reads such a coordinate as that double (a line whose two points share it, an LPI plane normal to the
!< axis) instead of running the exact stage, which coordinate ties on axis-aligned faces and edges need otherwise. The
!< oracle: points built on small integers so that every construction is an integer point, which `ip_approx` returns
!< exactly (checked). 500 points (explicit, LPI, SSI, EEI; half of the lines axis-parallel, half of the LPI planes
!< axis-normal), every pair on every axis: `ip_compare` against the sign of the integer difference, and `ip_lex_compare`
!< against the lexicographic order. The counts of compares with the shortcut on both, one and no point, and of ties,
!< must all be positive.
program fossil_test_implicit_points_axis

use fossil_implicit_points, only : implicit_point_t, IPK_EXPLICIT, ip_explicit, ip_lpi, ip_ssi, ip_eei, ip_approx, &
                                   ip_compare, ip_lex_compare
use penf,                   only : I4P, I8P, R8P

implicit none

integer(I4P), parameter :: NP = 500
type(implicit_point_t)  :: pts(NP)
real(R8P)               :: xs(3, NP)
logical                 :: fixed(3, NP)
integer(I4P)            :: i, j, k, nbad, noracle, nboth, none_, nnone, nties, nlex, sg, ex
integer(I8P)            :: seed
logical                 :: ok

seed = 20261117_I8P ; noracle = 0
do i = 1, NP
   call make_point(i)
   if (any(ip_approx(pts(i)) /= xs(:, i))) noracle = noracle + 1
enddo
nbad = 0 ; nboth = 0 ; none_ = 0 ; nnone = 0 ; nties = 0 ; nlex = 0
do i = 1, NP
   do j = 1, NP
      do k = 1, 3
         sg = ip_compare(pts(i), pts(j), k)
         ex = isign3(xs(k, i) - xs(k, j))
         if (sg /= ex) nbad = nbad + 1
         if (ex == 0) nties = nties + 1
         if (fixed(k, i) .and. fixed(k, j)) then
            nboth = nboth + 1
         elseif (fixed(k, i) .or. fixed(k, j)) then
            none_ = none_ + 1
         else
            nnone = nnone + 1
         endif
      enddo
      if (mod(i + j, 7) == 0) then ! a sample of the lexicographic compares
         ex = 0
         do k = 1, 3
            ex = isign3(xs(k, i) - xs(k, j))
            if (ex /= 0) exit
         enddo
         if (ip_lex_compare(pts(i), pts(j)) /= ex) nlex = nlex + 1
      endif
   enddo
enddo
ok = nbad == 0 .and. nlex == 0 .and. noracle == 0 .and. nboth > 0 .and. none_ > 0 .and. nnone > 0 .and. nties > 0
print '(A,L1,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A,I0)', 'ip_compare vs integer oracle: ', ok, '  mismatches ', nbad, &
   '  lex mismatches ', nlex, '  oracle failures ', noracle, '  compares fixed on both ', nboth, '  on one ', none_, &
   '  on none ', nnone, '  ties ', nties
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine make_point(i)
   !< Point i: explicit (1/8), LPI (1/2), SSI (3/16), EEI (3/16), at an integer target in [-2, 2]^3; its target and the
   !< axes the construction fixes to an input double.
   integer(I4P), intent(in) :: i
   real(R8P)                :: x(3), d(3), e1(3), e2(3), n(3), r(1), u
   integer(I4P)             :: kind_, ax, drop, kk, kj

   do
      call ints(x, 2) ; call prng(seed, r) ; u = r(1)
      fixed(:, i) = .false.
      if (u < 0.125_R8P) then
         kind_ = 1
      elseif (u < 0.625_R8P) then
         kind_ = 2
      elseif (u < 0.8125_R8P) then
         kind_ = 3
      else
         kind_ = 4
      endif
      call ints(d, 2) ; call ints(e1, 2) ; call ints(e2, 2)
      call axis(ax)
      call prng(seed, r)
      if (r(1) < 0.5_R8P) d(ax) = 0._R8P          ! an axis-parallel line
      select case (kind_)
      case (1)
         pts(i) = ip_explicit(x) ; fixed(:, i) = .true. ; exit
      case (2) ! LPI: line x - d, x + 2 d; plane x + e1, x + e2, x - e1 - e2
         call axis(kk) ; call prng(seed, r)
         if (r(1) < 0.5_R8P) then                 ! a plane normal to axis kk
            e1(kk) = 0._R8P ; e2(kk) = 0._R8P
         endif
         n = cross(e1, e2)
         if (all(n == 0._R8P) .or. dot_product(n, d) == 0._R8P) cycle
         pts(i) = ip_lpi(x - d, x + 2._R8P * d, x + e1, x + e2, x - e1 - e2)
         fixed(:, i) = d == 0._R8P .or. (e1 == 0._R8P .and. e2 == 0._R8P)
         exit
      case (3, 4) ! two lines through x in one plane: x - d, x + 2 d and x - e1, x + e1
         n = cross(d, e1)
         if (all(n == 0._R8P)) cycle
         if (kind_ == 3) then ! the projection dropping the largest normal component keeps them crossing
            drop = maxloc(abs(n), dim=1)
            pts(i) = ip_ssi(x - d, x + 2._R8P * d, x - e1, x + e1, drop)
         else
            pts(i) = ip_eei(x - d, x + 2._R8P * d, x - e1, x + e1)
         endif
         do kj = 1, 3 ! the first segment of the canonical order is the line whose two points share a coordinate
            fixed(kj, i) = pts(i)%v(kj, 1) == pts(i)%v(kj, 2)
         enddo
         exit
      endselect
   enddo
   if (pts(i)%kind == IPK_EXPLICIT) fixed(:, i) = .true.
   xs(:, i) = x
   endsubroutine make_point

   subroutine ints(v, m)
   !< Integers in [-m, m].
   real(R8P),    intent(out) :: v(3)
   integer(I4P), intent(in)  :: m
   real(R8P)                 :: r(3)

   call prng(seed, r)
   v = real(floor(r * (2 * m + 1)) - m, R8P)
   endsubroutine ints

   subroutine axis(k)
   !< An axis 1..3.
   integer(I4P), intent(out) :: k
   real(R8P)                 :: r(1)

   call prng(seed, r)
   k = min(3_I4P, 1_I4P + int(3._R8P * r(1), I4P))
   endsubroutine axis

   pure function cross(a, b) result(c)
   real(R8P), intent(in) :: a(3), b(3)
   real(R8P)             :: c(3)

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   pure function isign3(x) result(s)
   real(R8P), intent(in) :: x
   integer(I4P)          :: s

   s = 0_I4P
   if (x > 0._R8P) s = 1_I4P
   if (x < 0._R8P) s = -1_I4P
   endfunction isign3

   pure subroutine prng(seed, r)
   !< Park-Miller draws in (0, 1), fixed seed.
   integer(I8P), intent(inout) :: seed
   real(R8P),    intent(out)   :: r(:)
   integer(I4P)                :: i

   do i = 1, size(r)
      seed = mod(seed * 48271_I8P, 2147483647_I8P)
      r(i) = real(seed, R8P) / 2147483647._R8P
   enddo
   endsubroutine prng
endprogram fossil_test_implicit_points_axis
