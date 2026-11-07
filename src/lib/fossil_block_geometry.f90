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
!<
!< Validity (the sign of det J relative to the block's orientation `orient` = +-1):
!< - 3-D: det J has degree (2, 2, 2): 27 Bernstein coefficients per box from its values at {0, 1/2, 1}^3 of the box
!<   (B0 = f0, B2 = f2, B1 = 2 f1 - (f0 + f2) / 2 per direction), each with a floating error bound (the determinant's
!<   from the permanent of |J|, propagated through the transform). A box is positive when every coefficient exceeds its
!<   bound, non-negative when the others are exact zeros (coincident vertices: a collapsed edge), invalid as soon as a
!<   value is below minus its bound; otherwise it is cut in eight, down to VALIDITY_DEPTH levels.
!< - 2-D: det J has degree (1, 1): its coefficients are the four corner values, signed exactly by orient2d.
use penf, only : I4P, R8P
use fossil_predicates, only : orient2d_xy

implicit none
private
public :: trilinear_point, trilinear_jacobian, det3, hex_volume, hex_centroid, bilinear_face_vector
public :: bilinear_point, bilinear_jacobian, quad_area, quad_centroid, GAUSS2
public :: hex_validity, quad_validity, CELL_VALID, CELL_SINGULAR, CELL_INVALID, CELL_UNDETERMINED, VALIDITY_DEPTH

integer(I4P), parameter :: CELL_VALID        = 1_I4P !< det J has the block's sign everywhere.
integer(I4P), parameter :: CELL_SINGULAR     = 2_I4P !< det J >= 0, zero somewhere (collapsed edge or corner).
integer(I4P), parameter :: CELL_INVALID      = 3_I4P !< det J has the opposite sign somewhere, or vanishes everywhere.
integer(I4P), parameter :: CELL_UNDETERMINED = 4_I4P !< Not decided (3-D only).
integer(I4P), parameter :: VALIDITY_DEPTH    = 5_I4P !< Subdivision levels.
integer(I4P), parameter :: VALIDITY_BOXES    = 2048_I4P !< Boxes examined per cell at most.

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

   pure subroutine hex_validity(c, orient, status, jmin)
   !< Certified sign of det J over the trilinear cell (module notes).
   real(R8P),    intent(in)  :: c(3, 8) !< Corners.
   real(R8P),    intent(in)  :: orient  !< Block orientation, +1 or -1.
   integer(I4P), intent(out) :: status  !< CELL_*.
   real(R8P),    intent(out) :: jmin    !< Certified lower bound of orient det J (a negative value if invalid).
   real(R8P)                 :: blo(3, 64), bhi(3, 64), lo(3), hi(3), xi(3) !< Box stack, box, point.
   integer(I4P)              :: bdep(64), top, depth, a, b, e, used, ch !< Stack, counters.
   real(R8P)                 :: f(0:2, 0:2, 0:2), er(0:2, 0:2, 0:2), bc(0:2, 0:2, 0:2), eb(0:2, 0:2, 0:2) !< Values.
   logical                   :: singular, undecided !< Flags.

   status = CELL_VALID ; jmin = huge(1._R8P) ; singular = .false. ; undecided = .false.
   top = 1_I4P ; blo(:, 1) = 0._R8P ; bhi(:, 1) = 1._R8P ; bdep(1) = 0_I4P ; used = 0_I4P
   do while (top > 0)
      lo = blo(:, top) ; hi = bhi(:, top) ; depth = bdep(top) ; top = top - 1_I4P ; used = used + 1_I4P
      do e = 0, 2
         do b = 0, 2
            do a = 0, 2
               xi = lo + (hi - lo) * [a, b, e] * 0.5_R8P
               f(a, b, e) = orient * det3(trilinear_jacobian(c, xi))
               er(a, b, e) = 32._R8P * epsilon(1._R8P) * perm3(trilinear_jacobian_abs(c, xi))
               if (f(a, b, e) < -er(a, b, e)) then
                  status = CELL_INVALID ; jmin = f(a, b, e) ; return
               endif
            enddo
         enddo
      enddo
      call bernstein27(f, er, bc, eb)
      if (all(bc == 0._R8P .and. eb == 0._R8P)) then ! det J vanishes on the box, hence everywhere
         status = CELL_INVALID ; jmin = 0._R8P ; return
      elseif (all(bc > eb)) then
         jmin = min(jmin, minval(bc - eb))
      elseif (all(bc > eb .or. (bc == 0._R8P .and. eb == 0._R8P)) .and. any(bc > eb)) then
         singular = .true. ; jmin = min(jmin, 0._R8P)
      elseif (depth < VALIDITY_DEPTH .and. used + top + 8 <= VALIDITY_BOXES .and. top + 8 <= 64) then
         do ch = 0, 7
            top = top + 1_I4P ; bdep(top) = depth + 1_I4P
            do a = 1, 3
               if (btest(ch, a - 1)) then
                  blo(a, top) = 0.5_R8P * (lo(a) + hi(a)) ; bhi(a, top) = hi(a)
               else
                  blo(a, top) = lo(a) ; bhi(a, top) = 0.5_R8P * (lo(a) + hi(a))
               endif
            enddo
         enddo
      else
         undecided = .true. ; jmin = min(jmin, minval(bc - eb))
      endif
   enddo
   if (undecided) then
      status = CELL_UNDETERMINED
   elseif (singular) then
      status = CELL_SINGULAR
   endif
   endsubroutine hex_validity

   pure subroutine quad_validity(c, orient, status, jmin)
   !< Exact sign of det J over the bilinear cell: its four corner values by orient2d.
   real(R8P),    intent(in)  :: c(2, 4) !< Corners.
   real(R8P),    intent(in)  :: orient  !< Block orientation, +1 or -1.
   integer(I4P), intent(out) :: status  !< CELL_*.
   real(R8P),    intent(out) :: jmin    !< min of orient det J (at a corner).
   real(R8P)                 :: o(4)    !< Corner values.

   o(1) = orient2d_xy(c(1, 1), c(2, 1), c(1, 2), c(2, 2), c(1, 3), c(2, 3))
   o(2) = orient2d_xy(c(1, 1), c(2, 1), c(1, 2), c(2, 2), c(1, 4), c(2, 4))
   o(3) = orient2d_xy(c(1, 3), c(2, 3), c(1, 1), c(2, 1), c(1, 4), c(2, 4))
   o(4) = orient2d_xy(c(1, 4), c(2, 4), c(1, 3), c(2, 3), c(1, 2), c(2, 2))
   o = orient * o
   jmin = minval(o)
   if (all(o > 0._R8P)) then
      status = CELL_VALID
   elseif (any(o < 0._R8P) .or. all(o == 0._R8P)) then
      status = CELL_INVALID
   else
      status = CELL_SINGULAR
   endif
   endsubroutine quad_validity

   pure function trilinear_jacobian_abs(c, xi) result(jac)
   !< The trilinear Jacobian with every term in absolute value (bound for rounding errors).
   real(R8P), intent(in) :: c(3, 8)    !< Corners.
   real(R8P), intent(in) :: xi(3)      !< Local coordinates.
   real(R8P)             :: jac(3, 3)  !< Columns.
   real(R8P)             :: s(2, 3)    !< 1 - xi, xi per direction.

   s(1,:) = 1._R8P - xi ; s(2,:) = xi
   jac(:, 1) = abs(c(:, 2) - c(:, 1)) * s(1, 2) * s(1, 3) + abs(c(:, 4) - c(:, 3)) * s(2, 2) * s(1, 3) + &
               abs(c(:, 6) - c(:, 5)) * s(1, 2) * s(2, 3) + abs(c(:, 8) - c(:, 7)) * s(2, 2) * s(2, 3)
   jac(:, 2) = abs(c(:, 3) - c(:, 1)) * s(1, 1) * s(1, 3) + abs(c(:, 4) - c(:, 2)) * s(2, 1) * s(1, 3) + &
               abs(c(:, 7) - c(:, 5)) * s(1, 1) * s(2, 3) + abs(c(:, 8) - c(:, 6)) * s(2, 1) * s(2, 3)
   jac(:, 3) = abs(c(:, 5) - c(:, 1)) * s(1, 1) * s(1, 2) + abs(c(:, 6) - c(:, 2)) * s(2, 1) * s(1, 2) + &
               abs(c(:, 7) - c(:, 3)) * s(1, 1) * s(2, 2) + abs(c(:, 8) - c(:, 4)) * s(2, 1) * s(2, 2)
   endfunction trilinear_jacobian_abs

   pure real(R8P) function perm3(a) result(p)
   !< Permanent of a non-negative 3 x 3 matrix (the determinant with every sign +).
   real(R8P), intent(in) :: a(3, 3) !< Matrix.

   p = a(1, 1) * (a(2, 2) * a(3, 3) + a(3, 2) * a(2, 3)) + a(2, 1) * (a(3, 2) * a(1, 3) + a(1, 2) * a(3, 3)) + &
       a(3, 1) * (a(1, 2) * a(2, 3) + a(2, 2) * a(1, 3))
   endfunction perm3

   pure subroutine bernstein27(f, er, bc, eb)
   !< Bernstein coefficients of the tri-quadratic with values f at {0, 1/2, 1}^3, and their error bounds.
   real(R8P), intent(in)  :: f(0:2, 0:2, 0:2), er(0:2, 0:2, 0:2) !< Values, their error bounds.
   real(R8P), intent(out) :: bc(0:2, 0:2, 0:2), eb(0:2, 0:2, 0:2) !< Coefficients, their error bounds.
   real(R8P)              :: ab(0:2, 0:2, 0:2)                    !< Absolute transform (rounding of the transform).
   integer(I4P)           :: a, b                                 !< Counters.

   bc = f ; eb = er ; ab = abs(f)
   do b = 0, 2
      do a = 0, 2
         bc(1, a, b) = 2._R8P * bc(1, a, b) - 0.5_R8P * (bc(0, a, b) + bc(2, a, b))
         eb(1, a, b) = 2._R8P * eb(1, a, b) + 0.5_R8P * (eb(0, a, b) + eb(2, a, b))
         ab(1, a, b) = 2._R8P * ab(1, a, b) + 0.5_R8P * (ab(0, a, b) + ab(2, a, b))
      enddo
   enddo
   do b = 0, 2
      do a = 0, 2
         bc(a, 1, b) = 2._R8P * bc(a, 1, b) - 0.5_R8P * (bc(a, 0, b) + bc(a, 2, b))
         eb(a, 1, b) = 2._R8P * eb(a, 1, b) + 0.5_R8P * (eb(a, 0, b) + eb(a, 2, b))
         ab(a, 1, b) = 2._R8P * ab(a, 1, b) + 0.5_R8P * (ab(a, 0, b) + ab(a, 2, b))
      enddo
   enddo
   do b = 0, 2
      do a = 0, 2
         bc(a, b, 1) = 2._R8P * bc(a, b, 1) - 0.5_R8P * (bc(a, b, 0) + bc(a, b, 2))
         eb(a, b, 1) = 2._R8P * eb(a, b, 1) + 0.5_R8P * (eb(a, b, 0) + eb(a, b, 2))
         ab(a, b, 1) = 2._R8P * ab(a, b, 1) + 0.5_R8P * (ab(a, b, 0) + ab(a, b, 2))
      enddo
   enddo
   eb = 2._R8P * eb + 8._R8P * epsilon(1._R8P) * ab
   endsubroutine bernstein27

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
