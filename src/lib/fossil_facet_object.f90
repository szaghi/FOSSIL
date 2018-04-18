!< FOSSIL, facet class definition.

module fossil_facet_object
!< FOSSIL, facet class definition.

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : FR4P, I2P, I4P, R4P, R8P, str, ZeroR8P
use vecfor, only : face_normal3_R8P, sq_norm_R8P, normL2_R8P, vector_R4P, vector_R8P

implicit none
private
public :: facet_object
public :: FRLEN

integer(I4P), parameter :: FRLEN=80 !< Maximum length of facet record string.

type :: facet_object
   !< FOSSIL, facet class.
   type(vector_R8P) :: normal    !< Facet (outward) normal (versor), `(v2-v1).cross.(v3-v1)`.
   type(vector_R8P) :: vertex_1  !< Facet vertex 1.
   type(vector_R8P) :: vertex_2  !< Facet vertex 2.
   type(vector_R8P) :: vertex_3  !< Facet vertex 3.
   ! metrix
   ! triangle plane (parametric) equation: T(s,t) = B + s*E0 + t*E1
   type(vector_R8P) :: E0         !< Edge 2-1, `V2-V1`.
   type(vector_R8P) :: E1         !< Edge 3-1, `V3-V1`.
   real(R8P)        :: a=0._R8P   !< `E0.dot.E0`.
   real(R8P)        :: b=0._R8P   !< `E0.dot.E1`.
   real(R8P)        :: c=0._R8P   !< `E1.dot.E1`.
   real(R8P)        :: det=0._R8P !< `a*c - b*b`.
   ! triangle plane equation: ax + by + cz - d = 0, normal == [a, b, c]
   real(R8P) :: d=0._R8P !< `normal.dot.vertex_1`
   ! auxiliary
   type(vector_R8P) :: bb(2) !< Axis-aligned bounding box (AABB), bb(1)=min, bb(2)=max.
   contains
      ! public methods
      procedure, pass(self) :: check_normal          !< Check normal consistency.
      procedure, pass(self) :: compute_metrix        !< Compute local (plane) metrix.
      procedure, pass(self) :: distance              !< Compute the (unsigned, squared) distance from a point to the facet surface.
      procedure, pass(self) :: initialize            !< Initialize facet.
      procedure, pass(self) :: ray_intersect         !< Return true if facet is intersected by a ray.
      procedure, pass(self) :: load_from_file_ascii  !< Load facet from ASCII file.
      procedure, pass(self) :: load_from_file_binary !< Load facet from binary file.
      procedure, pass(self) :: sanitize_normal       !< Sanitize normal, make normal consistent with vertices.
      procedure, pass(self) :: save_into_file_ascii  !< Save facet into ASCII file.
      procedure, pass(self) :: save_into_file_binary !< Save facet into binary file.
      procedure, pass(self) :: solid_angle           !< Return the (projected) solid angle of the facet with respect the point.
      procedure, pass(self) :: winding_number        !< Return the winding number contribution of the facet with respect point.
      ! operators
      generic :: assignment(=) => facet_assign_facet !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: facet_assign_facet !< Operator `=`.
endtype facet_object

contains
   ! public methods
   elemental function check_normal(self) result(is_consistent)
   !< Check normal consistency.
   class(facet_object), intent(in) :: self          !< Facet.
   logical                         :: is_consistent !< Consistency check result.
   type(vector_R8P)                :: normal        !< Normal computed by means of vertices data.

   normal = face_normal3_R8P(pt1=self%vertex_1, pt2=self%vertex_2, pt3=self%vertex_3, norm='y')
   is_consistent = ((abs(normal%x - self%normal%x)<=2*ZeroR8P).and.&
                    (abs(normal%y - self%normal%y)<=2*ZeroR8P).and.&
                    (abs(normal%z - self%normal%z)<=2*ZeroR8P))
   endfunction check_normal

   elemental subroutine compute_metrix(self)
   !< Compute local (plane) metrix.
   class(facet_object), intent(inout) :: self !< Facet.

   call self%sanitize_normal

   self%E0  = self%vertex_2 - self%vertex_1
   self%E1  = self%vertex_3 - self%vertex_1
   self%a   = self%E0.dot.self%E0
   self%b   = self%E0.dot.self%E1
   self%c   = self%E1.dot.self%E1
   self%det = self%a * self%c - self%b * self%b

   self%d = self%normal.dot.self%vertex_1

   self%bb(1)%x = min(self%vertex_1%x, self%vertex_2%x, self%vertex_3%x)
   self%bb(1)%y = min(self%vertex_1%y, self%vertex_2%y, self%vertex_3%y)
   self%bb(1)%z = min(self%vertex_1%z, self%vertex_2%z, self%vertex_3%z)
   self%bb(2)%x = max(self%vertex_1%x, self%vertex_2%x, self%vertex_3%x)
   self%bb(2)%y = max(self%vertex_1%y, self%vertex_2%y, self%vertex_3%y)
   self%bb(2)%z = max(self%vertex_1%z, self%vertex_2%z, self%vertex_3%z)
   endsubroutine compute_metrix

   pure function distance(self, point)
   !< Compute the (unsigned, squared) distance from a point to the facet surface.
   !<
   !< @note Facet's metrix must be already computed.
   class(facet_object), intent(in) :: self                             !< Facet.
   type(vector_R8P),    intent(in) :: point                            !< Point.
   real(R8P)                       :: distance                         !< Closest distance from point to the facet.
   type(vector_R8P)                :: V1P                              !< `vertex_1-point`.
   real(R8P)                       :: d, e, f, s, t, sq, tq            !< Plane equation coefficients.
   real(R8P)                       :: tmp0, tmp1, numer, denom, invdet !< Temporary.

   associate(a=>self%a, b=>self%b, c=>self%c, det=>self%det)
   V1P = self%vertex_1 - point
   d = self%E0.dot.V1P
   e = self%E1.dot.V1P
   f = V1P.dot.V1P
   s = self%b * e - self%c * d
   t = self%b * d - self%a * e
   if (s+t <= det) then
      if (s < 0._R8P) then
         if (t < 0._R8P) then ! region 4
            if (e < 0._R8P) then
               sq = 0._R8P
               if (c >= -e) then
                  tq = -e / c
               else
                  tq = 1._R8P
               endif
            else
               if (d > 0._R8P) then
                  sq = 0._R8P
               else
                  if (a >= -d) then
                     sq = -d / a
                  else
                     sq = 1._R8P
                  endif
               endif
               tq = 0._R8P
            endif
         else ! region 3
            sq = 0._R8P
            if (e >= 0._R8P) then
               tq = 0._R8P
            else
               if (-e >= c) then
                  tq = 1._R8P
               else
                  tq = -e / c
               endif
            endif
         endif
      elseif (t < 0._R8P) then ! region 5
         if (d >= 0._R8P) then
            sq = 0._R8P
         else
            if (-d >= a) then
               sq = 1._R8P
            else
               sq = -d / a
            endif
         endif
         tq = 0._R8P
      else ! region 0
        invdet = 1._R8P / det
        sq = s * invdet
        tq = t * invdet
      endif
   else
      if (s < 0._R8P) then ! region 2
         tmp0 = b + d
         tmp1 = c + e
         if (tmp1 > tmp0) then
            numer = tmp1 - tmp0
            denom = a - 2._R8P * b + c
            if (numer >= denom) then
               sq = 1._R8P
            else
               sq = numer / denom
            endif
            tq = 1._R8P - sq
         else
            sq = 0._R8P
            if (tmp1 <= 0._R8P) then
               tq = 1._R8P
            else
               if (e >= 0._R8P) then
                  tq = 0._R8P
               else
                  tq = -e / c
               endif
            endif
         endif
      elseif (t < 0._R8P) then ! region 6
         tmp0 = a + d
         tmp1 = b + e
         if (tmp0 > tmp1) then
            numer = b + d - c - e
            denom = a - 2._R8P * b + c
            if (numer >= 0._R8P) then
               sq = 0._R8P
            else
               if (denom > -numer) then
                  sq = -numer / denom
               else
                  sq = 1._R8P
               endif
            endif
            tq = 1._R8P - sq
         else
            if (tmp0 <= 0._R8P) then
               sq = 1._R8P
            else
               if (d >= 0._R8P) then
                  sq = 0._R8P
               else
                  sq = -d / a
               endif
            endif
            tq = 0._R8P
         endif
      else ! region 1
         numer = c + e - b - d
         if (numer <= 0._R8P) then
            sq = 0._R8P
         else
            denom = a - 2._R8P * b + c
            if (numer >= denom) then
               sq = 1._R8P
            else
               sq = numer / denom
            endif
         endif
         tq = 1._R8P - sq
      endif
   endif
   distance = abs(a * sq * sq + 2._R8P * b * sq * tq + c * tq * tq + 2._R8P * d * sq + 2._R8P * e * tq + f)
   endassociate
   endfunction

   elemental subroutine initialize(self)
   !< Initialize facet.
   class(facet_object), intent(inout) :: self  !< Facet.
   type(facet_object)                 :: fresh !< Fresh instance of facet.

   self = fresh
   endsubroutine initialize

   pure function ray_intersect(self, ray_origin, ray_direction) result(intersect)
   !< Return true if facet is intersected by ray from point and oriented as ray direction vector.
   !<
   !< This based on Moller–Trumbore intersection algorithm.
   !<
   !< @note Facet's metrix must be already computed.
   class(facet_object), intent(in) :: self          !< Facet.
   type(vector_R8P),    intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),    intent(in) :: ray_direction !< Ray directio.
   logical                         :: intersect     !< Intersection test result.
   type(vector_R8P)                :: h, s, q       !< Projection vectors.
   real(R8P)                       :: a, f, u, v, t !< Baricentric abscissa.
   real(R8P), parameter            :: EPS=1e-7_R8P  !< Small espilon for round off errors control.

   intersect = .false.
   h = ray_direction.cross.self%E1
   a = self%E0.dot.h
   if ((a > -EPS).and.(a < EPS)) return
   f = 1._R8P / a
   s = ray_origin - self%vertex_1
   u = f * (s.dot.h)
   if ((u < 0._R8P).or.(u > 1._R8P)) return
   q = s.cross.self%E0
   v = f * ray_direction.dot.q
   if ((v < 0._R8P).or.(u + v > 1._R8P)) return
   t = f * self%E1.dot.q
   if (t > EPS) intersect = .true.
   endfunction ray_intersect

   subroutine load_from_file_ascii(self, file_unit)
   !< Load facet from ASCII file.
   class(facet_object), intent(inout) :: self      !< Facet.
   integer(I4P),        intent(in)    :: file_unit !< File unit.

   call load_facet_record(prefix='facet normal', record=self%normal)
   read(file_unit, *) ! outer loop
   call load_facet_record(prefix='vertex', record=self%vertex_1)
   call load_facet_record(prefix='vertex', record=self%vertex_2)
   call load_facet_record(prefix='vertex', record=self%vertex_3)
   read(file_unit, *) ! endloop
   read(file_unit, *) ! endfacet
   contains
      subroutine load_facet_record(prefix, record)
      !< Load a facet *record*, namely normal or vertex data.
      character(*),     intent(in)  :: prefix       !< Record prefix string.
      type(vector_R8P), intent(out) :: record       !< Record data.
      character(FRLEN)              :: facet_record !< Facet record string buffer.
      integer(I4P)                  :: i            !< Counter.

      read(file_unit, '(A)') facet_record
      i = index(string=facet_record, substring=prefix)
      if (i>0) then
         read(facet_record(i+len(prefix):), *) record%x, record%y, record%z
      else
         write(stderr, '(A)') 'error: impossible to read "'//prefix//'" from file unit "'//trim(str(file_unit))//'"!'
      endif
      endsubroutine load_facet_record
   endsubroutine load_from_file_ascii

   subroutine load_from_file_binary(self, file_unit)
   !< Load facet from binary file.
   class(facet_object), intent(inout) :: self       !< Facet.
   integer(I4P),        intent(in)    :: file_unit  !< File unit.
   integer(I2P)                       :: padding    !< Facet padding.
   real(R4P)                          :: triplet(3) !< Triplet record of R4P kind real.

   read(file_unit) triplet
   self%normal%x=real(triplet(1), R8P) ; self%normal%y=real(triplet(2), R8P) ; self%normal%z=real(triplet(3), R8P)
   read(file_unit) triplet
   self%vertex_1%x=real(triplet(1), R8P) ; self%vertex_1%y=real(triplet(2), R8P) ; self%vertex_1%z=real(triplet(3), R8P)
   read(file_unit) triplet
   self%vertex_2%x=real(triplet(1), R8P) ; self%vertex_2%y=real(triplet(2), R8P) ; self%vertex_2%z=real(triplet(3), R8P)
   read(file_unit) triplet
   self%vertex_3%x=real(triplet(1), R8P) ; self%vertex_3%y=real(triplet(2), R8P) ; self%vertex_3%z=real(triplet(3), R8P)
   read(file_unit) padding
   endsubroutine load_from_file_binary

   elemental subroutine sanitize_normal(self)
   !< Sanitize normal, make normal consistent with vertices.
   !<
   !<```fortran
   !< type(facet_object) :: facet
   !< facet%vertex_1 = -0.231369_R4P * ex_R4P + 0.0226865_R4P * ey_R4P + 1._R4P * ez_R4P
   !< facet%vertex_2 = -0.227740_R4P * ex_R4P + 0.0245457_R4P * ey_R4P + 0._R4P * ez_R4P
   !< facet%vertex_2 = -0.235254_R4P * ex_R4P + 0.0201881_R4P * ey_R4P + 0._R4P * ez_R4P
   !< call facet%sanitize_normal
   !< print "(3(F3.1,1X))", facet%normal%x, facet%normal%y, facet%normal%z
   !<```
   !=> -0.501673222 0.865057290 -2.12257713<<<
   class(facet_object), intent(inout) :: self !< Facet.

   self%normal = face_normal3_R8P(pt1=self%vertex_1, pt2=self%vertex_2, pt3=self%vertex_3, norm='y')
   endsubroutine sanitize_normal

   subroutine save_into_file_ascii(self, file_unit)
   !< Save facet into ASCII file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '  facet normal ', self%normal%x, ' ', self%normal%y, ' ', self%normal%z
   write(file_unit, '(A)')                            '    outer loop'
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex_1%x, ' ', self%vertex_1%y, ' ', self%vertex_1%z
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex_2%x, ' ', self%vertex_2%y, ' ', self%vertex_2%z
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex_3%x, ' ', self%vertex_3%y, ' ', self%vertex_3%z
   write(file_unit, '(A)')                            '    endloop'
   write(file_unit, '(A)')                            '  endfacet'
   endsubroutine save_into_file_ascii

   subroutine save_into_file_binary(self, file_unit)
   !< Save facet into binary file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.
   real(R4P)                       :: triplet(3) !< Triplet record of R4P kind real.

   triplet(1) = real(self%normal%x, R4P) ; triplet(2) = real(self%normal%y, R4P) ; triplet(3) = real(self%normal%z, R4P)
   write(file_unit) triplet
   triplet(1) = real(self%vertex_1%x, R4P) ; triplet(2) = real(self%vertex_1%y, R4P) ; triplet(3) = real(self%vertex_1%z, R4P)
   write(file_unit) triplet
   triplet(1) = real(self%vertex_2%x, R4P) ; triplet(2) = real(self%vertex_2%y, R4P) ; triplet(3) = real(self%vertex_2%z, R4P)
   write(file_unit) triplet
   triplet(1) = real(self%vertex_3%x, R4P) ; triplet(2) = real(self%vertex_3%y, R4P) ; triplet(3) = real(self%vertex_3%z, R4P)
   write(file_unit) triplet
   write(file_unit) 0_I2P
   endsubroutine save_into_file_binary

   pure function solid_angle(self, point)
   !< Return the (projected) solid angle of the facet with respect the point.
   class(facet_object), intent(in) :: self                      !< Facet.
   type(vector_R8P),    intent(in) :: point                     !< Point.
   real(R8P)                       :: solid_angle               !< Solid angle.
   type(vector_R8P)                :: R1, R2, R3                !< Edges from point to facet vertices.
   real(R8P)                       :: R1_norm, R2_norm, R3_norm !< Norms (L2) of edges from point to facet vertices.
   real(R8P)                       :: numerator                 !< Archtangent numerator.
   real(R8P)                       :: denominator               !< Archtangent denominator.

   R1 = self%vertex_1 - point ; R1_norm = R1%normL2()
   R2 = self%vertex_2 - point ; R2_norm = R2%normL2()
   R3 = self%vertex_3 - point ; R3_norm = R3%normL2()

   numerator = R1.dot.(R2.cross.R3)
   denominator = R1_norm * R2_norm * R3_norm + (R1.dot.R2) * R3_norm + &
                                               (R1.dot.R3) * R2_norm + &
                                               (R2.dot.R3) * R1_norm

   solid_angle = 2._R8P * atan2(numerator, denominator)
   endfunction solid_angle

   pure function winding_number(self, point)
   !< Return the winding number contribution of the facet with respect point.
   !<
   !< @note Raise error if the facet contains point.
   class(facet_object), intent(in) :: self           !< Facet.
   type(vector_R8P),    intent(in) :: point          !< Point.
   integer(I4P)                    :: winding_number !< Winding number contribution.
   integer(I4P)                    :: v1_sign        !< Vertex 1 integer sign.
   integer(I4P)                    :: v2_sign        !< Vertex 2 integer sign.
   integer(I4P)                    :: v3_sign        !< Vertex 3 integer sign.
   integer(I4P)                    :: face_boundary  !< Face boundaries count.

   v1_sign = vertex_sign(V=self%vertex_1, P=point)
   v2_sign = vertex_sign(V=self%vertex_2, P=point)
   v3_sign = vertex_sign(V=self%vertex_3, P=point)

                           face_boundary = 0
   if (v1_sign /= v2_sign) face_boundary = face_boundary + edge_sign(V1=self%vertex_1, V2=self%vertex_2, P=point)
   if (v2_sign /= v3_sign) face_boundary = face_boundary + edge_sign(V1=self%vertex_2, V2=self%vertex_3, P=point)
   if (v3_sign /= v1_sign) face_boundary = face_boundary + edge_sign(V1=self%vertex_3, V2=self%vertex_1, P=point)
   if (face_boundary == 0) then
      winding_number = 0
      return
   endif
   winding_number = triangle_sign(V1=self%vertex_1, V2=self%vertex_2, V3=self%vertex_3, P=point)
   contains
      pure function integer_sign(a)
      !< Return 1 if a is positive, -1 if it's negative and 0 if it's zero.
      real(R8P), intent(in) :: a            !< Generic coordinate.
      integer(I4P)          :: integer_sign !< Integer sign of "a".

      if     (a > 0._R8P) then
         integer_sign = 1
      elseif (a < 0._R8P) then
         integer_sign = -1
      else
         integer_sign = 0
      endif
      endfunction integer_sign

      pure function vertex_sign(V, P)
      !< Return the integer sign of the vertex V with respect to P.
      type(vector_R8P), intent(in) :: V           !< Vertex.
      type(vector_R8P), intent(in) :: P           !< Point of reference.
      integer(I4P)                 :: vertex_sign !< Integer sign of vertex "V" with respect "P".
      integer(I4P)                 :: signs(3)    !< Integer signs.

      signs(1) = integer_sign(V%x - P%x)
      signs(2) = integer_sign(V%y - P%y)
      signs(3) = integer_sign(V%z - P%z)
      if     (any(signs == -1)) then
         vertex_sign = -1
      elseif (any(signs == 1)) then
         vertex_sign = 1
      else
         vertex_sign = 0
         ! raise error: "V coincides with P"
      endif
      endfunction vertex_sign

      pure function edge_sign(V1, V2, P)
      !< Return the integer sign of the edge V1->V2 with respect to P.
      type(vector_R8P), intent(in) :: V1        !< Vertex 1.
      type(vector_R8P), intent(in) :: V2        !< Vertex 2.
      type(vector_R8P), intent(in) :: P         !< Point of reference.
      integer(I4P)                 :: edge_sign !< Integer sign of edge "V1->V2" with respect "P".
      integer(I4P)                 :: signs(3)  !< Integer signs.

      signs(1) = integer_sign((V1%y - P%y) * (V2%x - P%x) - (V1%x - P%x) * (V2%y - P%y))
      signs(2) = integer_sign((V1%z - P%z) * (V2%x - P%x) - (V1%x - P%x) * (V2%z - P%z))
      signs(3) = integer_sign((V1%z - P%z) * (V2%y - P%y) - (V1%y - P%y) * (V2%z - P%z))
      if     (any(signs == -1)) then
         edge_sign = -1
      elseif (any(signs == 1)) then
         edge_sign = 1
      else
         edge_sign = 0
         ! raise error: "V1->V2 is collinear with P"
      endif
      endfunction edge_sign

      pure function triangle_sign(V1, V2, V3, P)
      !< Return the integer sign of the triangle V1->V2->V3 with respect to P.
      type(vector_R8P), intent(in) :: V1            !< Vertex 1.
      type(vector_R8P), intent(in) :: V2            !< Vertex 2.
      type(vector_R8P), intent(in) :: V3            !< Vertex 3.
      type(vector_R8P), intent(in) :: P             !< Point of reference.
      integer(I4P)                 :: triangle_sign !< Integer sign of triangle "V1->V2->V3" with respect "P".
      real(R8P)                    :: m1_0, m1_1    !< Edge 1 2D coefficients.
      real(R8P)                    :: m2_0, m2_1    !< Edge 2 2D coefficients.
      real(R8P)                    :: m3_0, m3_1    !< Edge 3 2D coefficients.

      m1_0 = V1%x - P%x
      m1_1 = V1%y - P%y
      m2_0 = V2%x - P%x
      m2_1 = V2%y - P%y
      m3_0 = V3%x - P%x
      m3_1 = V3%y - P%y
      triangle_sign = integer_sign((m1_0 * m2_1 - m1_1 * m2_0) * (V3%z - P%z) + &
                                   (m2_0 * m3_1 - m2_1 * m3_0) * (V1%z - P%z) + &
                                   (m3_0 * m1_1 - m3_1 * m1_0) * (V2%z - P%z))
      ! if triangle sign == 0 raise errort: "V1->V2->V3 complanar with P"
      endfunction triangle_sign
   endfunction winding_number

   ! private methods
   ! `=` operator
   pure subroutine facet_assign_facet(lhs, rhs)
   !< Operator `=`.
   class(facet_object), intent(inout) :: lhs !< Left hand side.
   type(facet_object),  intent(in)    :: rhs !< Right hand side.

   lhs%normal = rhs%normal
   lhs%vertex_1 = rhs%vertex_1
   lhs%vertex_2 = rhs%vertex_2
   lhs%vertex_3 = rhs%vertex_3
   lhs%E0 = rhs%E0
   lhs%E1 = rhs%E1
   lhs%a = rhs%a
   lhs%b = rhs%b
   lhs%c = rhs%c
   lhs%d = rhs%d
   lhs%det = rhs%det
   lhs%bb = rhs%bb
   endsubroutine facet_assign_facet
endmodule fossil_facet_object
