!< FOSSIL, geometry kernels of the cells of curvilinear blocks (issue #71).
module fossil_block_geometry
!< FOSSIL, geometry kernels of the cells of curvilinear blocks (issue #71).
!<
!< Plain-real kernels on the corners of one cell: c(3, 8) in 3-D, corner (a, b, c) in {0, 1}^3 at 1 + a + 2 b + 4 c
!< (2-D: c(2, 4), corner 1 + a + 2 b). The cell is the trilinear (bilinear) map of the unit cube (square).
!<
!< - Volume: det J of the trilinear map has degree 2 in each variable, so the 2-point Gauss rule per direction
!<   integrates it exactly (8 determinants); the centroid integrand x det J has degree 3, exact with the same points.
!< - Face area vector of a bilinear face p00, p10, p11, p01: int x_s x x_t ds dt = 1/2 (p11 - p00) x (p01 - p10).
!< - 2-D: area 1/2 (p11 - p00) x (p01 - p10) (exact for the bilinear quadrilateral), centroid by 2 x 2 Gauss.
use penf, only : I4P, R8P

implicit none
private
public :: trilinear_point, trilinear_jacobian, det3, hex_volume, hex_centroid, bilinear_face_vector
public :: bilinear_point, bilinear_jacobian, quad_area, quad_centroid, GAUSS2

real(R8P), parameter :: GAUSS2(2) = [0.5_R8P - 0.5_R8P / sqrt(3._R8P), 0.5_R8P + 0.5_R8P / sqrt(3._R8P)] !< 2-point Gauss.

contains
   pure function trilinear_point(c, xi) result(x)
   !< Point of the trilinear map at xi.
   real(R8P), intent(in) :: c(3, 8) !< Corners.
   real(R8P), intent(in) :: xi(3)   !< Local coordinates.
   real(R8P)             :: x(3)    !< Point.
   real(R8P)             :: w(8)    !< Shape functions.

   call shape8(xi, w)
   x = matmul(c, w)
   endfunction trilinear_point

   pure function trilinear_jacobian(c, xi) result(jac)
   !< Jacobian [x_xi, x_eta, x_zeta] of the trilinear map at xi.
   real(R8P), intent(in) :: c(3, 8)    !< Corners.
   real(R8P), intent(in) :: xi(3)      !< Local coordinates.
   real(R8P)             :: jac(3, 3)  !< Columns: derivatives.
   real(R8P)             :: s(2, 3)    !< 1 - xi, xi per direction.

   s(1,:) = 1._R8P - xi ; s(2,:) = xi
   jac(:, 1) = (c(:, 2) - c(:, 1)) * s(1, 2) * s(1, 3) + (c(:, 4) - c(:, 3)) * s(2, 2) * s(1, 3) + &
               (c(:, 6) - c(:, 5)) * s(1, 2) * s(2, 3) + (c(:, 8) - c(:, 7)) * s(2, 2) * s(2, 3)
   jac(:, 2) = (c(:, 3) - c(:, 1)) * s(1, 1) * s(1, 3) + (c(:, 4) - c(:, 2)) * s(2, 1) * s(1, 3) + &
               (c(:, 7) - c(:, 5)) * s(1, 1) * s(2, 3) + (c(:, 8) - c(:, 6)) * s(2, 1) * s(2, 3)
   jac(:, 3) = (c(:, 5) - c(:, 1)) * s(1, 1) * s(1, 2) + (c(:, 6) - c(:, 2)) * s(2, 1) * s(1, 2) + &
               (c(:, 7) - c(:, 3)) * s(1, 1) * s(2, 2) + (c(:, 8) - c(:, 4)) * s(2, 1) * s(2, 2)
   endfunction trilinear_jacobian

   pure real(R8P) function det3(a) result(d)
   !< Determinant of a 3 x 3 matrix (columns a(:, 1..3)): a1 . (a2 x a3).
   real(R8P), intent(in) :: a(3, 3) !< Matrix.

   d = a(1, 1) * (a(2, 2) * a(3, 3) - a(3, 2) * a(2, 3)) + a(2, 1) * (a(3, 2) * a(1, 3) - a(1, 2) * a(3, 3)) + &
       a(3, 1) * (a(1, 2) * a(2, 3) - a(2, 2) * a(1, 3))
   endfunction det3

   pure real(R8P) function hex_volume(c) result(vol)
   !< Exact volume of the trilinear cell (signed: negative for a left-handed cell).
   real(R8P), intent(in) :: c(3, 8) !< Corners.
   integer(I4P)          :: a, b, e !< Gauss indices.

   vol = 0._R8P
   do e = 1, 2
      do b = 1, 2
         do a = 1, 2
            vol = vol + det3(trilinear_jacobian(c, [GAUSS2(a), GAUSS2(b), GAUSS2(e)]))
         enddo
      enddo
   enddo
   vol = vol / 8._R8P
   endfunction hex_volume

   pure function hex_centroid(c) result(x)
   !< Exact centroid of the trilinear cell (the corner mean if the volume is zero).
   real(R8P), intent(in) :: c(3, 8) !< Corners.
   real(R8P)             :: x(3)    !< Centroid.
   real(R8P)             :: d, vol, xi(3) !< Work.
   integer(I4P)          :: a, b, e !< Gauss indices.

   x = 0._R8P ; vol = 0._R8P
   do e = 1, 2
      do b = 1, 2
         do a = 1, 2
            xi = [GAUSS2(a), GAUSS2(b), GAUSS2(e)]
            d = det3(trilinear_jacobian(c, xi))
            x = x + d * trilinear_point(c, xi) ; vol = vol + d
         enddo
      enddo
   enddo
   if (vol /= 0._R8P) then
      x = x / vol
   else
      x = sum(c, dim=2) / 8._R8P
   endif
   endfunction hex_centroid

   pure function bilinear_face_vector(p00, p10, p11, p01) result(a)
   !< Exact area vector of the bilinear face p00, p10, p11, p01 (right-hand rule over that order).
   real(R8P), intent(in) :: p00(3), p10(3), p11(3), p01(3) !< Corners.
   real(R8P)             :: a(3)                           !< Area vector.
   real(R8P)             :: d1(3), d2(3)                   !< Diagonals.

   d1 = p11 - p00 ; d2 = p01 - p10
   a = 0.5_R8P * [d1(2) * d2(3) - d1(3) * d2(2), d1(3) * d2(1) - d1(1) * d2(3), d1(1) * d2(2) - d1(2) * d2(1)]
   endfunction bilinear_face_vector

   pure function bilinear_point(c, xi) result(x)
   !< Point of the bilinear map at xi.
   real(R8P), intent(in) :: c(2, 4) !< Corners.
   real(R8P), intent(in) :: xi(2)   !< Local coordinates.
   real(R8P)             :: x(2)    !< Point.

   x = c(:, 1) * (1._R8P - xi(1)) * (1._R8P - xi(2)) + c(:, 2) * xi(1) * (1._R8P - xi(2)) + &
       c(:, 3) * (1._R8P - xi(1)) * xi(2) + c(:, 4) * xi(1) * xi(2)
   endfunction bilinear_point

   pure function bilinear_jacobian(c, xi) result(jac)
   !< Jacobian [x_xi, x_eta] of the bilinear map at xi.
   real(R8P), intent(in) :: c(2, 4)   !< Corners.
   real(R8P), intent(in) :: xi(2)     !< Local coordinates.
   real(R8P)             :: jac(2, 2) !< Columns: derivatives.

   jac(:, 1) = (c(:, 2) - c(:, 1)) * (1._R8P - xi(2)) + (c(:, 4) - c(:, 3)) * xi(2)
   jac(:, 2) = (c(:, 3) - c(:, 1)) * (1._R8P - xi(1)) + (c(:, 4) - c(:, 2)) * xi(1)
   endfunction bilinear_jacobian

   pure real(R8P) function quad_area(c) result(area)
   !< Exact signed area of the bilinear quadrilateral.
   real(R8P), intent(in) :: c(2, 4) !< Corners.

   area = 0.5_R8P * ((c(1, 4) - c(1, 1)) * (c(2, 3) - c(2, 2)) - (c(2, 4) - c(2, 1)) * (c(1, 3) - c(1, 2)))
   endfunction quad_area

   pure function quad_centroid(c) result(x)
   !< Exact centroid of the bilinear quadrilateral (the corner mean if the area is zero).
   real(R8P), intent(in) :: c(2, 4) !< Corners.
   real(R8P)             :: x(2)    !< Centroid.
   real(R8P)             :: jac(2, 2), d, area, xi(2) !< Work.
   integer(I4P)          :: a, b    !< Gauss indices.

   x = 0._R8P ; area = 0._R8P
   do b = 1, 2
      do a = 1, 2
         xi = [GAUSS2(a), GAUSS2(b)]
         jac = bilinear_jacobian(c, xi) ; d = jac(1, 1) * jac(2, 2) - jac(2, 1) * jac(1, 2)
         x = x + d * bilinear_point(c, xi) ; area = area + d
      enddo
   enddo
   if (area /= 0._R8P) then
      x = x / area
   else
      x = sum(c, dim=2) / 4._R8P
   endif
   endfunction quad_centroid

   pure subroutine shape8(xi, w)
   !< Trilinear shape functions.
   real(R8P), intent(in)  :: xi(3) !< Local coordinates.
   real(R8P), intent(out) :: w(8)  !< Weights.
   real(R8P)              :: s(2, 3) !< 1 - xi, xi.
   integer(I4P)           :: a, b, e !< Corner bits.

   s(1,:) = 1._R8P - xi ; s(2,:) = xi
   do e = 1, 2
      do b = 1, 2
         do a = 1, 2
            w(a + 2 * (b - 1) + 4 * (e - 1)) = s(a, 1) * s(b, 2) * s(e, 3)
         enddo
      enddo
   enddo
   endsubroutine shape8
endmodule fossil_block_geometry
