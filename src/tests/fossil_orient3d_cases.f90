!< FOSSIL test-side module: orient3d and line_crossing_before cases, and their batched loops (issue #100 step 3c).
module fossil_orient3d_cases
!< FOSSIL test-side module: orient3d and line_crossing_before cases, and their batched loops (issue #100 step 3c).
!<
!< The cases (`orient3d_quadruples`, `lcb_cases`) exercise every stage of `orient3d_xyz` and the perturbation of
!< `line_crossing_before`; `orient3d_many` / `lcb_many` evaluate them in one `!$acc parallel loop` (host code under
!< gfortran). The distance bench compares the device loops with the scalar host calls bit for bit (`orient3d device vs
!< serial mismatches`); `fossil_test_orient3d` runs the same loops on the host against the known answers.
!<
!< Quadruples [12, n] (a, b, c, d), with `expect` +1 / -1 (the exact sign is known), 0 (exactly zero) or 2 (unknown):
!<   1. the near-coplanar grid of `fossil_test_orient3d` at T = 2^0 and 2^40 (a moved by i, j, k ulps about the plane
!<      of b, c, d), in four orders (one even, two odd permutations): 24576;
!<   2. exactly coplanar dyadic points on z = 2 x - 3 y + 0.5 (exact differences: stage 2a): 500;
!<   3. four points on a shared coordinate plane (the exact zero), and the same with one point an ulp off it: 1800;
!<   4. coincident points (each pair, and all four): 700;
!<   5. near-coplanar random points far from the origin (inexact differences: stage 2b, 24 products), a moved by
!<      -2..2 ulps: 2000;
!<   6. well-separated random quadruples (the filter): 2000.
!< line_crossing_before cases: triangles on the planes w = w0 + p u + q v (dyadic, so the crossing is exact) of the
!< three axes, both orientations, at w0 = 0 and 2^40; points on the plane (the perturbation: the first nonzero normal
!< component), one ulp and one unit off it, entering and leaving; and collinear (degenerate) triangles. One ulp off a
!< crossing at 0 is 2^-52, not the subnormal `nearest(0, 1)`: orient3d is exact barring underflow of products below
!< ~1e-100 (its header), and on subnormals the nvfortran host (flush to zero, item 9b of #101) and the device differ.
use fossil_line_query, only : line_crossing_before
use fossil_predicates, only : orient3d_xyz
use penf,              only : I4P, R8P

implicit none
private
public :: orient3d_quadruples, lcb_cases, orient3d_many, lcb_many

contains
   subroutine orient3d_quadruples(q, expect)
   !< The orient3d cases (see the module header).
   real(R8P),    allocatable, intent(out) :: q(:,:)    !< Quadruples [12, n]: a, b, c, d.
   integer(I4P), allocatable, intent(out) :: expect(:) !< Known sign (+1, -1), exact zero (0) or unknown (2).
   real(R8P)    :: t, u, a(3), b(3), c(3), d(3), s0, r(12), w(4) !< Work.
   integer(I4P) :: m, i, j, k, n, ax, kk, l, ns !< Work.
   integer, allocatable :: seed(:) !< Seed.

   allocate(q(12, 24576 + 500 + 1800 + 700 + 2000 + 2000), expect(24576 + 500 + 1800 + 700 + 2000 + 2000))
   n = 0_I4P
   ! 1. near-coplanar grid (as fossil_test_orient3d)
   do m = 1, 2
      t = 0._R8P ; if (m == 2) t = 2._R8P**40
      u = spacing(t + 0.5_R8P)
      b = [t + 12._R8P, t + 3._R8P, 9._R8P] ; c = [t + 24._R8P, t + 7._R8P, 17._R8P] ; d = [t + 5._R8P, t + 30._R8P, -25._R8P]
      a = [t + 0.5_R8P, t + 0.5_R8P, 1._R8P] ! well separated: the fixed sign of a displacement in +z
      s0 = sign(1._R8P, naive(a, b, c, d))
      do k = -1, 1
         do j = 0, 31
            do i = 0, 31
               a = [t + 0.5_R8P + i * u, t + 0.5_R8P + j * u, (i - j) * u + k * u]
               call put(q, expect, n, a, b, c, d, int(s0, I4P) * k)
               call put(q, expect, n, b, c, a, d, int(s0, I4P) * k)
               call put(q, expect, n, b, a, c, d, -int(s0, I4P) * k)
               call put(q, expect, n, b, c, d, a, -int(s0, I4P) * k)
            enddo
         enddo
      enddo
   enddo
   call random_seed(size=ns)
   allocate(seed(ns)) ; seed = [(1003 + i, i = 1, ns)]
   call random_seed(put=seed)
   ! 2. exactly coplanar, dyadic coordinates on z = 2 x - 3 y + 0.5
   do kk = 1, 500
      call random_number(r)
      r = real(nint(r * 2._R8P**20), R8P) / 2._R8P**10 - 512._R8P
      call put(q, expect, n, plane_point(r(1), r(2)), plane_point(r(3), r(4)), plane_point(r(5), r(6)), &
               plane_point(r(7), r(8)), 0_I4P)
   enddo
   ! 3. a shared coordinate plane (exact zero), and one point an ulp off it (unknown sign)
   do ax = 1, 3
      do kk = 1, 300
         call random_number(r)
         r = 2._R8P * r - 1._R8P
         a = r(1:3) ; b = r(4:6) ; c = r(7:9) ; d = r(10:12)
         a(ax) = r(1) * 1.e3_R8P ; b(ax) = a(ax) ; c(ax) = a(ax) ; d(ax) = a(ax)
         call put(q, expect, n, a, b, c, d, 0_I4P)
         l = mod(kk, 4) + 1
         select case (l)
         case (1) ; a(ax) = nearest(a(ax), 1._R8P)
         case (2) ; b(ax) = nearest(b(ax), -1._R8P)
         case (3) ; c(ax) = nearest(c(ax), 1._R8P)
         case default ; d(ax) = nearest(d(ax), -1._R8P)
         endselect
         call put(q, expect, n, a, b, c, d, 2_I4P)
      enddo
   enddo
   ! 4. coincident points
   do kk = 1, 100
      call random_number(r)
      r = 4._R8P * r - 2._R8P
      a = r(1:3) ; b = r(4:6) ; c = r(7:9) ; d = r(10:12)
      call put(q, expect, n, a, a, c, d, 0_I4P) ; call put(q, expect, n, a, b, a, d, 0_I4P)
      call put(q, expect, n, a, b, c, a, 0_I4P) ; call put(q, expect, n, a, b, b, d, 0_I4P)
      call put(q, expect, n, a, b, c, b, 0_I4P) ; call put(q, expect, n, a, b, c, c, 0_I4P)
      call put(q, expect, n, d, d, d, d, 0_I4P)
   enddo
   ! 5. near-coplanar, far from the origin, inexact differences
   do kk = 1, 400
      call random_number(r)
      b = 1.e3_R8P + r(1:3) ; c = 1.e3_R8P + r(4:6) ; d = 1.e3_R8P + r(7:9)
      call random_number(w)
      do l = -2, 2
         a = b + w(1) * (c - b) + w(2) * (d - b)
         a(3) = a(3) + l * spacing(a(3))
         call put(q, expect, n, a, b, c, d, 2_I4P)
      enddo
   enddo
   ! 6. well separated
   do kk = 1, 2000
      call random_number(r)
      r = 2._R8P * r - 1._R8P
      call put(q, expect, n, r(1:3), r(4:6), r(7:9), r(10:12), 2_I4P)
   enddo
   if (n /= size(expect, kind=I4P)) error stop 'orient3d_quadruples: count'
   endsubroutine orient3d_quadruples

   subroutine put(q, expect, n, pa, pb, pc, pd, ex)
   !< Append one quadruple.
   real(R8P),    intent(inout) :: q(:,:)    !< Quadruples.
   integer(I4P), intent(inout) :: expect(:) !< Expected.
   integer(I4P), intent(inout) :: n         !< Count.
   real(R8P),    intent(in)    :: pa(3), pb(3), pc(3), pd(3) !< Points.
   integer(I4P), intent(in)    :: ex        !< Expected.

   n = n + 1_I4P
   q(1:3, n) = pa ; q(4:6, n) = pb ; q(7:9, n) = pc ; q(10:12, n) = pd ; expect(n) = ex
   endsubroutine put

   subroutine lcb_cases(tri, ax, ab, c, en, expect)
   !< The line_crossing_before cases (see the module header) with the exact answer.
   real(R8P),    allocatable, intent(out) :: tri(:,:)  !< Triangles [9, n] (A, B, C).
   integer(I4P), allocatable, intent(out) :: ax(:)     !< Line axes.
   real(R8P),    allocatable, intent(out) :: ab(:,:)   !< Line positions [2, n].
   real(R8P),    allocatable, intent(out) :: c(:)      !< Point coordinates along the line.
   logical,      allocatable, intent(out) :: en(:)     !< Entering flags.
   logical,      allocatable, intent(out) :: expect(:) !< Exact answers.
   real(R8P), parameter :: PQ(5) = [-1._R8P, -0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P] !< Plane slopes.
   real(R8P), parameter :: UV(2, 3) = reshape([1._R8P, 1._R8P, 0.5_R8P, 2.25_R8P, 2.75_R8P, 0.25_R8P], [2, 3]) !< Points.
   real(R8P)    :: w0, loc(3, 3), nl(3), nn(3), tt, cc, o !< Work.
   integer(I4P) :: k, i1, i2, ip, iq, io, iw, ipt, ic, ie, v, n, m, nmax !< Work.

   nmax = 3 * 2 * 5 * 5 * 3 * 5 * 2 * 2 + 3 * 2 * 3 * 2
   allocate(tri(9, nmax), ax(nmax), ab(2, nmax), c(nmax), en(nmax), expect(nmax))
   n = 0_I4P
   do k = 1, 3
      i1 = mod(k, 3_I4P) + 1_I4P ; i2 = mod(k + 1_I4P, 3_I4P) + 1_I4P
      do iw = 1, 2
         w0 = 0._R8P ; if (iw == 2) w0 = 2._R8P**40
         do ip = 1, 5
            do iq = 1, 5
               do io = 1, 2 ! orientation: (A, B, C) or (A, C, B)
                  ! local (u, v, w) vertices: A = (0, 0), B = (4, 0), C = (0, 4) on w = w0 + p u + q v
                  loc(:, 1) = [0._R8P, 0._R8P, w0]
                  loc(:, 2) = [4._R8P, 0._R8P, w0 + 4._R8P * PQ(ip)]
                  loc(:, 3) = [0._R8P, 4._R8P, w0 + 4._R8P * PQ(iq)]
                  if (io == 2) loc(:, 2:3) = loc(:, 3:2:-1)
                  ! local normal (B - A) x (C - A), exact (small dyadic)
                  nl = [(loc(2, 2) - loc(2, 1)) * (loc(3, 3) - loc(3, 1)) - (loc(3, 2) - loc(3, 1)) * (loc(2, 3) - loc(2, 1)), &
                        (loc(3, 2) - loc(3, 1)) * (loc(1, 3) - loc(1, 1)) - (loc(1, 2) - loc(1, 1)) * (loc(3, 3) - loc(3, 1)), &
                        (loc(1, 2) - loc(1, 1)) * (loc(2, 3) - loc(2, 1)) - (loc(2, 2) - loc(2, 1)) * (loc(1, 3) - loc(1, 1))]
                  ! global normal: local (u, v, w) = global (i1, i2, k); the cyclic relabelling keeps the orientation
                  nn(i1) = nl(1) ; nn(i2) = nl(2) ; nn(k) = nl(3)
                  do ipt = 1, 3
                     tt = w0 + PQ(ip) * UV(1, ipt) + PQ(iq) * UV(2, ipt) ! exact
                     do ic = 1, 5
                        select case (ic)
                        case (1) ; cc = tt
                        case (2) ; cc = nearest(tt, 1._R8P)
                           if (tt == 0._R8P) cc = epsilon(1._R8P) ! an ulp of 1, not a subnormal (see the header)
                        case (3) ; cc = nearest(tt, -1._R8P)
                           if (tt == 0._R8P) cc = -epsilon(1._R8P)
                        case (4) ; cc = tt + 1._R8P
                        case default ; cc = tt - 1._R8P
                        endselect
                        ! orient3d(A, B, C, P) = -n_axis (c - t), perturbed on the plane: -first nonzero of n (x, y, z)
                        if (cc /= tt) then
                           o = -sign(1._R8P, nn(k)) * sign(1._R8P, cc - tt)
                        else
                           o = 0._R8P
                           do m = 1, 3
                              if (nn(m) /= 0._R8P) then
                                 o = -sign(1._R8P, nn(m)) ; exit
                              endif
                           enddo
                        endif
                        do ie = 1, 2
                           n = n + 1_I4P
                           do v = 1, 3
                              tri(3 * (v - 1) + i1, n) = loc(1, v) ; tri(3 * (v - 1) + i2, n) = loc(2, v)
                              tri(3 * (v - 1) + k, n) = loc(3, v)
                           enddo
                           ax(n) = k ; ab(:, n) = UV(:, ipt) ; c(n) = cc ; en(n) = ie == 1
                           if (en(n)) then
                              expect(n) = o > 0._R8P
                           else
                              expect(n) = o < 0._R8P
                           endif
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
      ! degenerate (collinear) triangles: always before
      do iw = 1, 2
         w0 = 0._R8P ; if (iw == 2) w0 = 2._R8P**40
         do ipt = 1, 3
            do ie = 1, 2
               n = n + 1_I4P
               loc(:, 1) = [0._R8P, 0._R8P, w0] ; loc(:, 2) = [1._R8P, 1._R8P, w0 + 1._R8P]
               loc(:, 3) = [2._R8P, 2._R8P, w0 + 2._R8P]
               do v = 1, 3
                  tri(3 * (v - 1) + i1, n) = loc(1, v) ; tri(3 * (v - 1) + i2, n) = loc(2, v)
                  tri(3 * (v - 1) + k, n) = loc(3, v)
               enddo
               ax(n) = k ; ab(:, n) = UV(:, ipt) ; c(n) = w0 ; en(n) = ie == 1 ; expect(n) = .true.
            enddo
         enddo
      enddo
   enddo
   if (n /= nmax) error stop 'lcb_cases: count'
   endsubroutine lcb_cases

   subroutine orient3d_many(n, q, o)
   !< o(i) = orient3d_xyz of quadruple i, one `!$acc parallel loop` (the device path; host code without OpenACC).
   integer(I4P), intent(in)  :: n       !< Quadruples.
   real(R8P),    intent(in)  :: q(12, n) !< Quadruples.
   real(R8P),    intent(out) :: o(n)     !< Results.
   real(R8P)                 :: pa(3), pb(3), pc(3), pd(3) !< One quadruple.
   integer(I4P)              :: i, m !< Counters.

   !$acc parallel loop gang vector copyin(q) copyout(o) private(pa, pb, pc, pd)
   do i = 1, n
      do m = 1, 3
         pa(m) = q(m, i) ; pb(m) = q(m + 3, i) ; pc(m) = q(m + 6, i) ; pd(m) = q(m + 9, i)
      enddo
      o(i) = orient3d_xyz(pa, pb, pc, pd)
   enddo
   endsubroutine orient3d_many

   subroutine lcb_many(n, tri, ax, ab, c, en, before)
   !< before(i) = line_crossing_before of case i, one `!$acc parallel loop` (as `orient3d_many`).
   integer(I4P), intent(in)  :: n         !< Cases.
   real(R8P),    intent(in)  :: tri(9, n) !< Triangles.
   integer(I4P), intent(in)  :: ax(n)     !< Axes.
   real(R8P),    intent(in)  :: ab(2, n)  !< Line positions.
   real(R8P),    intent(in)  :: c(n)      !< Point coordinates.
   logical,      intent(in)  :: en(n)     !< Entering flags.
   logical,      intent(out) :: before(n) !< Results.
   real(R8P)                 :: t3(3, 3) !< One triangle.
   integer(I4P)              :: i, m, v  !< Counters.

   !$acc parallel loop gang vector copyin(tri, ax, ab, c, en) copyout(before) private(t3)
   do i = 1, n
      do v = 1, 3
         do m = 1, 3
            t3(m, v) = tri(3 * (v - 1) + m, i)
         enddo
      enddo
      before(i) = line_crossing_before(t3, ax(i), ab(1, i), ab(2, i), c(i), en(i))
   enddo
   endsubroutine lcb_many

   pure function plane_point(x, y) result(p)
   !< The point (x, y) on z = 2 x - 3 y + 0.5 (exact for the dyadic inputs).
   real(R8P), intent(in) :: x, y !< Coordinates.
   real(R8P)             :: p(3) !< Point.

   p = [x, y, 2._R8P * x - 3._R8P * y + 0.5_R8P]
   endfunction plane_point

   pure function naive(a, b, c, d) result(det)
   !< Naive determinant det[a - d; b - d; c - d].
   real(R8P), intent(in) :: a(3), b(3), c(3), d(3) !< Points.
   real(R8P)             :: det                    !< Determinant.
   real(R8P)             :: p(3), q(3), w(3)       !< Rows.

   p = a - d ; q = b - d ; w = c - d
   det = p(3) * (q(1) * w(2) - w(1) * q(2)) + q(3) * (w(1) * p(2) - p(1) * w(2)) + w(3) * (p(1) * q(2) - q(1) * p(2))
   endfunction naive
endmodule fossil_orient3d_cases
