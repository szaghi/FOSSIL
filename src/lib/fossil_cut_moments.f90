!< FOSSIL, geometric moments of the cut-cell pieces: monomials and exact integration kernels (issue #70).
module fossil_cut_moments
!< FOSSIL, geometric moments of the cut-cell pieces: monomials and exact integration kernels (issue #70).
!<
!< The moments of a piece are the integrals of the monomials X^a Y^b (Z^c) of degree a + b (+ c) <= p over it, in the
!< coordinates X = (x - x_c) / h_x of its cell (centre x_c, size h): M_ab(c) = int X^a Y^b (Z^c) dV, of the dimension
!< of a volume (an area in 2-D), O(|piece|) wherever the cell lies. Order: by total degree, then by decreasing a, then
!< by decreasing b (`mono_exponents`): 1, X, Y, (Z,) X^2, XY, ... The first ones are the volume and the volume times the
!< scaled centroid offset.
!<
!< Volumes come from the divergence theorem along x: int X^a Y^b Z^c dV = oint F n_x dA with
!< F = h_x X^(a+1) / (a+1) Y^b Z^c (dF/dx = X^a Y^b Z^c), so only the boundary parts with n_x /= 0 count: the walls
!< and the two faces x = x_c -+ h_x / 2 (X = -+1/2). 2-D areas by Green: int X^a Y^b dA = oint h_x X^(a+1) / (a+1) Y^b
!< dy along the boundary, counter-clockwise. The kernels are exact for these polynomials: Gauss-Legendre with 3 nodes on a
!< segment (degree 5) and a symmetric 7-point rule on a triangle (degree 5, Dunavant 1985 / Strang-Fix); the order is
!< capped at MOMENTS_MAX = 4, so the integrands have degree at most 5.
use penf, only : I4P, R8P

implicit none
private
public :: MOMENTS_MAX, mono_count, mono_exponents, segment_green_moments, triangle_moments, cell_moments
public :: tet_rule, triangle_rule, tensor_rule

integer(I4P), parameter :: MOMENTS_MAX = 4_I4P !< Highest moment order.
real(R8P),    parameter :: GX(3) = [0.5_R8P - sqrt(0.15_R8P), 0.5_R8P, 0.5_R8P + sqrt(0.15_R8P)] !< Gauss nodes on [0, 1].
real(R8P),    parameter :: GW(3) = [5._R8P / 18._R8P, 8._R8P / 18._R8P, 5._R8P / 18._R8P] !< Their weights.

contains
   pure integer(I4P) function mono_count(d, p) result(n)
   !< Monomials of degree <= p in d variables (2 or 3).
   integer(I4P), intent(in) :: d, p !< Dimension, order.

   if (d == 2_I4P) then
      n = (p + 1_I4P) * (p + 2_I4P) / 2_I4P
   else
      n = (p + 1_I4P) * (p + 2_I4P) * (p + 3_I4P) / 6_I4P
   endif
   endfunction mono_count

   pure function mono_exponents(d, p) result(e)
   !< Exponents [d, mono_count(d, p)] of the monomials, in the order of the module notes.
   integer(I4P), intent(in)  :: d, p !< Dimension, order.
   integer(I4P), allocatable :: e(:,:) !< Exponents.
   integer(I4P)              :: deg, a, b, k !< Work.

   allocate(e(d, mono_count(d, p))) ; k = 0_I4P
   do deg = 0, p
      do a = deg, 0, -1
         if (d == 2_I4P) then
            k = k + 1_I4P ; e(:, k) = [a, deg - a]
         else
            do b = deg - a, 0, -1
               k = k + 1_I4P ; e(:, k) = [a, b, deg - a - b]
            enddo
         endif
      enddo
   enddo
   endfunction mono_exponents

   pure subroutine segment_green_moments(p0, q0, e, sc, mom)
   !$acc routine seq
   !< Add the Green integrals of the straight path p0 -> q0 (scaled cell coordinates (X, Y)) to the 2-D moments:
   !< mom(k) += sc int X^(a+1) / (a+1) Y^b dY, sc = h_x h_y (the cell area), exact (3-point Gauss, degree <= 5). Summed
   !< counter-clockwise around a region it gives int X^a Y^b dx dy.
   real(R8P),    intent(in)    :: p0(2), q0(2) !< End points (scaled cell coordinates).
   integer(I4P), intent(in)    :: e(:,:)       !< Exponents [2, n].
   real(R8P),    intent(in)    :: sc           !< Scale: the cell area.
   real(R8P),    intent(inout) :: mom(:)       !< Moments [n].
   integer(I4P)                :: k, g          !< Work.
   real(R8P)                   :: pt(2), dy     !< Gauss point, dY.

   dy = q0(2) - p0(2)
   if (dy == 0._R8P) return
   do g = 1, 3
      pt = p0 + GX(g) * (q0 - p0)
      do k = 1, size(mom)
         mom(k) = mom(k) + sc * GW(g) * dy * pt(1) ** (e(1, k) + 1) / (e(1, k) + 1) * pt(2) ** e(2, k)
      enddo
   enddo
   endsubroutine segment_green_moments

   pure subroutine triangle_moments(t1, t2, t3, e, vc, ap, vol, wall)
   !< Add the wall terms of the triangle t1 t2 t3 (scaled cell coordinates, normal by the right-hand rule):
   !< vol(k) += vc int X^(a+1) / (a+1) Y^b Z^c N_X dA (scaled vector area: vc times it is the physical divergence term)
   !< and wall(k) += int X^a Y^b Z^c dA_phys (ap the physical area), exact for degree <= 5 (the symmetric 7-point rule,
   !< coordinates (6 -+ sqrt 15) / 21, (9 +- 2 sqrt 15) / 21 and weights (155 +- sqrt 15) / 1200, in closed form).
   real(R8P),    intent(in)    :: t1(3), t2(3), t3(3) !< Vertices.
   integer(I4P), intent(in)    :: e(:,:)              !< Exponents [3, n].
   real(R8P),    intent(in)    :: vc                  !< Scale: the cell volume.
   real(R8P),    intent(in)    :: ap                  !< Physical area of the triangle.
   real(R8P),    intent(inout) :: vol(:)              !< Volume terms (divergence along X) [n].
   real(R8P),    intent(inout) :: wall(:)             !< Wall moments [n].
   real(R8P), parameter :: S15 = sqrt(15._R8P)       !< sqrt 15.
   real(R8P), parameter :: A1 = (9._R8P - 2._R8P * S15) / 21._R8P, B1 = (6._R8P + S15) / 21._R8P !< Orbit 1.
   real(R8P), parameter :: A2 = (9._R8P + 2._R8P * S15) / 21._R8P, B2 = (6._R8P - S15) / 21._R8P !< Orbit 2.
   real(R8P), parameter :: W0 = 9._R8P / 40._R8P, W1 = (155._R8P + S15) / 1200._R8P, W2 = (155._R8P - S15) / 1200._R8P
   real(R8P) :: bar(3, 7), w(7), nv(3), pt(3)        !< Barycentric points, weights, scaled vector area, point.
   integer(I4P) :: g, k                               !< Work.

   bar(:, 1) = [1._R8P, 1._R8P, 1._R8P] / 3._R8P
   bar(:, 2) = [A1, B1, B1] ; bar(:, 3) = [B1, A1, B1] ; bar(:, 4) = [B1, B1, A1]
   bar(:, 5) = [A2, B2, B2] ; bar(:, 6) = [B2, A2, B2] ; bar(:, 7) = [B2, B2, A2]
   w = [W0, W1, W1, W1, W2, W2, W2]
   nv = 0.5_R8P * [(t2(2) - t1(2)) * (t3(3) - t1(3)) - (t2(3) - t1(3)) * (t3(2) - t1(2)), &
                   (t2(3) - t1(3)) * (t3(1) - t1(1)) - (t2(1) - t1(1)) * (t3(3) - t1(3)), &
                   (t2(1) - t1(1)) * (t3(2) - t1(2)) - (t2(2) - t1(2)) * (t3(1) - t1(1))] ! vector area
   do g = 1, 7
      pt = bar(1, g) * t1 + bar(2, g) * t2 + bar(3, g) * t3
      do k = 1, size(vol)
         vol(k) = vol(k) + vc * w(g) * nv(1) * pt(1) ** (e(1, k) + 1) / (e(1, k) + 1) * pt(2) ** e(2, k) * pt(3) ** e(3, k)
         wall(k) = wall(k) + w(g) * ap * pt(1) ** e(1, k) * pt(2) ** e(2, k) * pt(3) ** e(3, k)
      enddo
   enddo
   endsubroutine triangle_moments

   pure subroutine tet_rule(a, t1, t2, t3, pt, w)
   !< Quadrature on the signed tetrahedron (a, t1, t2, t3), exact for polynomials of degree <= 4 (MOMENTS_MAX): the
   !< Duffy map x = a + r (t1 - a + s (t2 - t1) + s t (t3 - t2)), Jacobian r^2 s det, Gauss-Legendre with 4 nodes in r
   !< (r^(2+p), degree 7) and 3 in s (s^(1+p), degree 5) and t (36 points). Weights carry the sign of det: the cones
   !< of a closed oriented boundary sum to the region.
   real(R8P), intent(in)  :: a(3), t1(3), t2(3), t3(3) !< Apex and base.
   real(R8P), intent(out) :: pt(3, 36), w(36)        !< Points, weights.
   real(R8P), parameter   :: R4(4) = [0.0694318442029737_R8P, 0.3300094782075719_R8P, 0.6699905217924281_R8P, &
                                      0.9305681557970263_R8P] !< Gauss nodes on [0, 1] (4).
   real(R8P), parameter   :: W4(4) = [0.1739274225687269_R8P, 0.3260725774312731_R8P, 0.3260725774312731_R8P, &
                                      0.1739274225687269_R8P] !< Their weights.
   real(R8P)              :: det, u1(3), u2(3), u3(3) !< Work.
   integer(I4P)           :: i, j, k, n              !< Work.

   u1 = t1 - a ; u2 = t2 - t1 ; u3 = t3 - t2
   det = u1(1) * (u2(2) * u3(3) - u2(3) * u3(2)) - u1(2) * (u2(1) * u3(3) - u2(3) * u3(1)) + &
         u1(3) * (u2(1) * u3(2) - u2(2) * u3(1))
   n = 0_I4P
   do i = 1, 4
      do j = 1, 3
         do k = 1, 3
            n = n + 1_I4P
            pt(:, n) = a + R4(i) * (u1 + GX(j) * (u2 + GX(k) * u3))
            w(n) = W4(i) * GW(j) * GW(k) * R4(i) * R4(i) * GX(j) * det
         enddo
      enddo
   enddo
   endsubroutine tet_rule

   pure subroutine triangle_rule(t1, t2, t3, pt, w, signed)
   !< The symmetric 7-point rule (degree 5) on the triangle t1 t2 t3 in 2-D (signed: weight sign = orientation) or 3-D
   !< (weights = area fractions, positive).
   real(R8P), intent(in)  :: t1(:), t2(:), t3(:)  !< Vertices (2 or 3 coordinates).
   real(R8P), intent(out) :: pt(:,:), w(7)         !< Points [size(t1), 7], weights.
   logical,   intent(in)  :: signed                !< 2-D signed area (else the 3-D area).
   real(R8P), parameter   :: S15 = sqrt(15._R8P)  !< sqrt 15.
   real(R8P), parameter   :: A1 = (9._R8P - 2._R8P * S15) / 21._R8P, B1 = (6._R8P + S15) / 21._R8P !< Orbit 1.
   real(R8P), parameter   :: A2 = (9._R8P + 2._R8P * S15) / 21._R8P, B2 = (6._R8P - S15) / 21._R8P !< Orbit 2.
   real(R8P), parameter   :: W0 = 9._R8P / 40._R8P, W1 = (155._R8P + S15) / 1200._R8P, W2 = (155._R8P - S15) / 1200._R8P
   real(R8P)              :: bar(3, 7), ar, c3(3), e1(3), e2(3) !< Work.
   integer(I4P)           :: g                     !< Work.

   bar(:, 1) = [1._R8P, 1._R8P, 1._R8P] / 3._R8P
   bar(:, 2) = [A1, B1, B1] ; bar(:, 3) = [B1, A1, B1] ; bar(:, 4) = [B1, B1, A1]
   bar(:, 5) = [A2, B2, B2] ; bar(:, 6) = [B2, A2, B2] ; bar(:, 7) = [B2, B2, A2]
   if (signed) then
      ar = 0.5_R8P * ((t2(1) - t1(1)) * (t3(2) - t1(2)) - (t2(2) - t1(2)) * (t3(1) - t1(1)))
   else
      e1 = t2 - t1 ; e2 = t3 - t1
      c3 = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
      ar = 0.5_R8P * norm2(c3)
   endif
   do g = 1, 7
      pt(:, g) = bar(1, g) * t1 + bar(2, g) * t2 + bar(3, g) * t3
   enddo
   w = ar * [W0, W1, W1, W1, W2, W2, W2]
   endsubroutine triangle_rule

   pure subroutine tensor_rule(lo, hi, pt, w)
   !< Tensor Gauss-Legendre rule (3 nodes per axis, degree 5) on the box [lo, hi] (2-D or 3-D).
   real(R8P), intent(in)  :: lo(:), hi(:) !< Corners.
   real(R8P), intent(out) :: pt(:,:), w(:) !< Points [d, 3^d], weights.
   integer(I4P)           :: i, j, k, n   !< Work.

   n = 0_I4P
   do k = 1, merge(3, 1, size(lo) == 3)
      do j = 1, 3
         do i = 1, 3
            n = n + 1_I4P
            pt(1, n) = lo(1) + GX(i) * (hi(1) - lo(1)) ; pt(2, n) = lo(2) + GX(j) * (hi(2) - lo(2))
            w(n) = GW(i) * GW(j) * product(hi - lo)
            if (size(lo) == 3) then
               pt(3, n) = lo(3) + GX(k) * (hi(3) - lo(3)) ; w(n) = w(n) * GW(k)
            endif
         enddo
      enddo
   enddo
   endsubroutine tensor_rule

   pure function cell_moments(e, vc) result(m)
   !$acc routine seq
   !< Moments of a whole cell: vc prod_i int_{-1/2}^{1/2} X^a_i dX (zero for an odd exponent, 1 / ((a + 1) 2^a) else).
   integer(I4P), intent(in) :: e(:,:) !< Exponents [d, n].
   real(R8P),    intent(in) :: vc     !< Cell volume (area).
   real(R8P)                :: m(size(e, 2)) !< Moments.
   integer(I4P)             :: k, i   !< Work.

   do k = 1, size(e, 2)
      m(k) = vc
      do i = 1, size(e, 1)
         if (mod(e(i, k), 2_I4P) == 1_I4P) then
            m(k) = 0._R8P
         else
            m(k) = m(k) / ((e(i, k) + 1) * 2._R8P ** e(i, k))
         endif
      enddo
   enddo
   endfunction cell_moments
endmodule fossil_cut_moments
