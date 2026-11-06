!< FOSSIL test-side module: closed bodies with known thickness and gaps (issue #62).
module fossil_thickness_bodies
!< FOSSIL test-side module: closed bodies with known thickness and gaps (issue #62).
!<
!< - `make_boxes`: n axis-aligned boxes, the k-th shifted along x by (k - 1)(width + gap); 12 facets each.
!< - `make_sphere`: a unit icosphere (the icosahedron subdivided `levels` times, vertices on the sphere).
!< - `make_wing`: a NACA 00tt section of unit chord (np points per side, cosine spacing, closed trailing edge)
!<   extruded over z in [0, span] with nz layers, the end caps fanned from the chord midpoint (the section is convex).
!< - `naca_half_thickness`: y_t(x) and its slope.
!< Every facet is oriented away from an inner point of its (convex) body.
use fossil, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: make_boxes, make_sphere, make_wing, naca_half_thickness

real(R8P), parameter :: PI = 3.141592653589793_R8P

contains
   subroutine make_boxes(lo, hi, n, gap, ff)
   !< n boxes [lo, hi], the k-th shifted by (k - 1)(hi_x - lo_x + gap) along x.
   real(R8P),                       intent(in)  :: lo(3), hi(3), gap
   integer(I4P),                    intent(in)  :: n
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: p(3, 8), q(3), sh
   integer(I4P), parameter :: Q4(4, 6) = reshape([1, 3, 7, 5,  2, 4, 8, 6,  1, 2, 6, 5,  3, 4, 8, 7,  1, 2, 4, 3, &
                                                  5, 6, 8, 7], [4, 6])
   integer(I4P)                                 :: b, j, m, kk

   allocate(ff(12 * n))
   m = 0
   do b = 1, n
      sh = (b - 1) * (hi(1) - lo(1) + gap)
      do kk = 1, 8
         q = [merge(hi(1), lo(1), btest(kk - 1, 0)), merge(hi(2), lo(2), btest(kk - 1, 1)), &
              merge(hi(3), lo(3), btest(kk - 1, 2))]
         p(:, kk) = q + [sh, 0._R8P, 0._R8P]
      enddo
      q = 0.5_R8P * (lo + hi) + [sh, 0._R8P, 0._R8P]
      do j = 1, 6
         m = m + 1 ; call tri(ff(m), p(:, Q4(1, j)), p(:, Q4(2, j)), p(:, Q4(3, j)), q)
         m = m + 1 ; call tri(ff(m), p(:, Q4(1, j)), p(:, Q4(3, j)), p(:, Q4(4, j)), q)
      enddo
   enddo
   endsubroutine make_boxes

   subroutine make_sphere(levels, ff)
   !< Unit icosphere.
   integer(I4P),                    intent(in)  :: levels
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P), allocatable                       :: tv(:,:,:), nv(:,:,:)
   real(R8P)                                    :: v(3, 12), a, b, m12(3), m23(3), m31(3)
   integer(I4P), parameter :: I20(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12,  2, 6, 10, &
      6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10,  5, 10, 6, &
      3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   integer(I4P)                                 :: j, l, nt0

   a = 1._R8P ; b = (1._R8P + sqrt(5._R8P)) / 2._R8P
   v = reshape([-a, b, 0._R8P,  a, b, 0._R8P,  -a, -b, 0._R8P,  a, -b, 0._R8P,  0._R8P, -a, b,  0._R8P, a, b, &
                0._R8P, -a, -b,  0._R8P, a, -b,  b, 0._R8P, -a,  b, 0._R8P, a,  -b, 0._R8P, -a,  -b, 0._R8P, a], [3, 12])
   allocate(tv(3, 3, 20))
   do j = 1, 20
      do l = 1, 3
         tv(:, l, j) = v(:, I20(l, j)) / norm2(v(:, I20(l, j)))
      enddo
   enddo
   do l = 1, levels
      nt0 = size(tv, 3)
      allocate(nv(3, 3, 4 * nt0))
      do j = 1, nt0
         m12 = tv(:, 1, j) + tv(:, 2, j) ; m12 = m12 / norm2(m12)
         m23 = tv(:, 2, j) + tv(:, 3, j) ; m23 = m23 / norm2(m23)
         m31 = tv(:, 3, j) + tv(:, 1, j) ; m31 = m31 / norm2(m31)
         nv(:, 1, 4*j-3) = tv(:, 1, j) ; nv(:, 2, 4*j-3) = m12 ; nv(:, 3, 4*j-3) = m31
         nv(:, 1, 4*j-2) = m12 ; nv(:, 2, 4*j-2) = tv(:, 2, j) ; nv(:, 3, 4*j-2) = m23
         nv(:, 1, 4*j-1) = m31 ; nv(:, 2, 4*j-1) = m23 ; nv(:, 3, 4*j-1) = tv(:, 3, j)
         nv(:, 1, 4*j)   = m12 ; nv(:, 2, 4*j)   = m23 ; nv(:, 3, 4*j)   = m31
      enddo
      call move_alloc(nv, tv)
   enddo
   allocate(ff(size(tv, 3)))
   do j = 1, size(tv, 3)
      call tri(ff(j), tv(:, 1, j), tv(:, 2, j), tv(:, 3, j), [0._R8P, 0._R8P, 0._R8P])
   enddo
   endsubroutine make_sphere

   pure subroutine naca_half_thickness(tt, x, yt, dyt)
   !< NACA 00tt half thickness y_t(x) (closed trailing edge) and its slope.
   real(R8P), intent(in)  :: tt, x
   real(R8P), intent(out) :: yt, dyt

   yt = 5._R8P * tt * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x ** 2 + 0.2843_R8P * x ** 3 - &
                       0.1036_R8P * x ** 4)
   dyt = 5._R8P * tt * (0.5_R8P * 0.2969_R8P / sqrt(x) - 0.1260_R8P - 2._R8P * 0.3516_R8P * x + 3._R8P * 0.2843_R8P * x ** 2 - &
                        4._R8P * 0.1036_R8P * x ** 3)
   endsubroutine naca_half_thickness

   subroutine make_wing(tt, np, span, nz, ff)
   !< NACA 00tt wing of unit chord, np points per side, extruded over [0, span] in nz layers.
   real(R8P),                       intent(in)  :: tt, span
   integer(I4P),                    intent(in)  :: np, nz
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: xs(0:np), ys(0:np), d, u(2 * np), v(2 * np), z(0:nz)
   integer(I4P)                                 :: i, k, n, m, i2

   do i = 0, np
      xs(i) = 0.5_R8P * (1._R8P - cos(PI * i / np))
      call naca_half_thickness(tt, max(xs(i), 0._R8P), ys(i), d)
   enddo
   ys(np) = 0._R8P
   n = 0
   do i = np, 1, -1                              ! upper side, trailing edge to leading edge
      n = n + 1 ; u(n) = xs(i) ; v(n) = ys(i)
   enddo
   do i = 0, np - 1                              ! lower side, leading edge to trailing edge (TE point once)
      n = n + 1 ; u(n) = xs(i) ; v(n) = -ys(i)
   enddo
   z = [(span * k / nz, k = 0, nz)]
   allocate(ff(2 * n * nz + 2 * n))
   m = 0
   do k = 0, nz - 1
      do i = 1, n
         i2 = mod(i, n) + 1
         m = m + 1 ; call tri(ff(m), [u(i), v(i), z(k)], [u(i2), v(i2), z(k)], [u(i2), v(i2), z(k + 1)], &
                              [0.5_R8P, 0._R8P, 0.5_R8P * (z(k) + z(k + 1))])
         m = m + 1 ; call tri(ff(m), [u(i), v(i), z(k)], [u(i2), v(i2), z(k + 1)], [u(i), v(i), z(k + 1)], &
                              [0.5_R8P, 0._R8P, 0.5_R8P * (z(k) + z(k + 1))])
      enddo
   enddo
   do i = 1, n
      i2 = mod(i, n) + 1
      m = m + 1 ; call tri(ff(m), [0.5_R8P, 0._R8P, 0._R8P], [u(i), v(i), 0._R8P], [u(i2), v(i2), 0._R8P], &
                           [0.5_R8P, 0._R8P, 0.5_R8P * span])
      m = m + 1 ; call tri(ff(m), [0.5_R8P, 0._R8P, span], [u(i), v(i), span], [u(i2), v(i2), span], &
                           [0.5_R8P, 0._R8P, 0.5_R8P * span])
   enddo
   endsubroutine make_wing

   subroutine tri(ft, p1, p2, p3, inner)
   !< A facet from three vertices, oriented with its normal away from `inner`.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3), inner(3)
   real(R8P)                         :: e1(3), e2(3), n(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(n, (p1 + p2 + p3) / 3._R8P - inner) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri
endmodule fossil_thickness_bodies
