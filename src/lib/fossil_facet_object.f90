!< FOSSIL, facet class definition.

module fossil_facet_object
!< FOSSIL, facet class definition.

use fossil_utils, only : EPS, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : FR4P, I2P, I4P, R4P, R8P, str, ZeroR8P
use vecfor, only : angle_R8P, face_normal3_R8P, normL2_R8P, vector_R8P

implicit none
private
public :: facet_object

type :: facet_object
   !< FOSSIL, facet class.
   type(vector_R8P) :: normal    !< Facet (outward) normal (versor), `(v2-v1).cross.(v3-v1)`.
   type(vector_R8P) :: vertex_1  !< Facet vertex 1.
   type(vector_R8P) :: vertex_2  !< Facet vertex 2.
   type(vector_R8P) :: vertex_3  !< Facet vertex 3.
   ! metrix
   ! triangle plane parametric equation: T(s,t) = B + s*E12 + t*E13
   type(vector_R8P) :: E12        !< Edge 1-2, `V2-V1`.
   type(vector_R8P) :: E13        !< Edge 1-3, `V3-V1`.
   real(R8P)        :: a=0._R8P   !< `E12.dot.E12`.
   real(R8P)        :: b=0._R8P   !< `E12.dot.E13`.
   real(R8P)        :: c=0._R8P   !< `E13.dot.E13`.
   real(R8P)        :: det=0._R8P !< `a*c - b*b`.
   ! triangle plane equation: nx*x + ny*y + nz*z - d = 0, normal == [nx, ny, nz]
   real(R8P) :: d=0._R8P !< `normal.dot.vertex_1`
   ! auxiliary
   type(vector_R8P) :: bb(2) !< Axis-aligned bounding box (AABB), bb(1)=min, bb(2)=max.
   ! connectivity
   integer(I4P)              :: id                     !< Facet global ID.
   integer(I4P)              :: fcon_edge_12=0_I4P     !< Connected face ID along edge 1-2.
   integer(I4P)              :: fcon_edge_23=0_I4P     !< Connected face ID along edge 2-3.
   integer(I4P)              :: fcon_edge_31=0_I4P     !< Connected face ID along edge 3-1.
   integer(I4P), allocatable :: vertex_1_occurrence(:) !< List of vertex 1 "occurrencies", list of facets global ID containing it.
   integer(I4P), allocatable :: vertex_2_occurrence(:) !< List of vertex 2 "occurrencies", list of facets global ID containing it.
   integer(I4P), allocatable :: vertex_3_occurrence(:) !< List of vertex 3 "occurrencies", list of facets global ID containing it.
   contains
      ! public methods
      procedure, pass(self) :: add_vertex_occurrence       !< Add vertex occurence.
      procedure, pass(self) :: check_normal                !< Check normal consistency.
      procedure, pass(self) :: check_vertices_occurrencies !< Check if vertices of facet are *identical* to the ones of other facet.
      procedure, pass(self) :: compute_metrix              !< Compute local (plane) metrix.
      procedure, pass(self) :: compute_normal              !< Compute normal by means of vertices data.
      procedure, pass(self) :: destroy                     !< Destroy facet.
      procedure, pass(self) :: distance                    !< Compute the (unsigned, squared) distance from a point to facet.
      procedure, pass(self) :: do_ray_intersect            !< Return true if facet is intersected by a ray.
      procedure, pass(self) :: initialize                  !< Initialize facet.
      procedure, pass(self) :: load_from_file_ascii        !< Load facet from ASCII file.
      procedure, pass(self) :: load_from_file_binary       !< Load facet from binary file.
      procedure, pass(self) :: make_normal_consistent      !< Make normal of other facet consistent with self.
      procedure, pass(self) :: reverse_normal              !< Reverse facet normal.
      procedure, pass(self) :: save_into_file_ascii        !< Save facet into ASCII file.
      procedure, pass(self) :: save_into_file_binary       !< Save facet into binary file.
      procedure, pass(self) :: solid_angle                 !< Return the (projected) solid angle of the facet with respect point.
      procedure, pass(self) :: tetrahedron_volume          !< Return the volume of tetrahedron built by facet and a given apex.
      procedure, pass(self) :: update_connectivity         !< Update facet connectivity.
      procedure, pass(self) :: vertex_global_id            !< Return the vertex global id given the local one.
      ! operators
      generic :: assignment(=) => facet_assign_facet !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: facet_assign_facet           !< Operator `=`.
      procedure, pass(self) :: edge_connection_in_other_ref !< Return the edge of connection in the other reference.
      procedure, pass(self) :: flip_edge                    !< Flip facet edge.
endtype facet_object

contains
   ! public methods
   elemental subroutine add_vertex_occurrence(self, vertex_id, facet_id)
   !< Add vertex occurrence.
   class(facet_object), intent(inout) :: self      !< Facet.
   integer(I4P),        intent(in)    :: vertex_id !< Vertex ID in local numeration, 1, 2 or 3.
   integer(I4P),        intent(in)    :: facet_id  !< Other facet ID containing vertex.

   select case(vertex_id)
   case(1)
      call add_occurrence(occurrence=self%vertex_1_occurrence)
   case(2)
      call add_occurrence(occurrence=self%vertex_2_occurrence)
   case(3)
      call add_occurrence(occurrence=self%vertex_3_occurrence)
   endselect
   contains
      pure subroutine add_occurrence(occurrence)
      !< Add new occurrence into a generic occurrencies array.
      integer(I4P), allocatable, intent(inout) :: occurrence(:)     !< Occurrences array.
      integer(I4P), allocatable                :: occurrence_tmp(:) !< Temporary occurences array.
      integer(I4P)                             :: no                !< Occurrences number.

      if (allocated(occurrence)) then
         no = size(occurrence, dim=1)
         allocate(occurrence_tmp(1:no+1))
         occurrence_tmp(1:no) = occurrence
         occurrence_tmp(no+1) = facet_id
         call move_alloc(from=occurrence_tmp, to=occurrence)
      else
         allocate(occurrence(1))
         occurrence(1) = facet_id
      endif
      endsubroutine add_occurrence
   endsubroutine add_vertex_occurrence

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

   pure subroutine check_vertices_occurrencies(self, other)
   !< Check if vertices of facet are *identical* (with tollerance) to the ones of other facet.
   !<
   !< If multiple occurrencies are found the counters are updated.
   class(facet_object), intent(inout) :: self  !< Facet.
   type(facet_object),  intent(inout) :: other !< Other facet.

   if     (check_pair(self%vertex_1, other%vertex_1)) then
      call self%add_vertex_occurrence( vertex_id=1, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=1, facet_id=self%id)
   elseif (check_pair(self%vertex_1, other%vertex_2)) then
      call self%add_vertex_occurrence( vertex_id=1, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=2, facet_id=self%id)
   elseif (check_pair(self%vertex_1, other%vertex_3)) then
      call self%add_vertex_occurrence( vertex_id=1, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=3, facet_id=self%id)
   endif
   if     (check_pair(self%vertex_2, other%vertex_1)) then
      call self%add_vertex_occurrence( vertex_id=2, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=1, facet_id=self%id)
   elseif (check_pair(self%vertex_2, other%vertex_2)) then
      call self%add_vertex_occurrence( vertex_id=2, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=2, facet_id=self%id)
   elseif (check_pair(self%vertex_2, other%vertex_3)) then
      call self%add_vertex_occurrence( vertex_id=2, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=3, facet_id=self%id)
   endif
   if     (check_pair(self%vertex_3, other%vertex_1)) then
      call self%add_vertex_occurrence( vertex_id=3, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=1, facet_id=self%id)
   elseif (check_pair(self%vertex_3, other%vertex_2)) then
      call self%add_vertex_occurrence( vertex_id=3, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=2, facet_id=self%id)
   elseif (check_pair(self%vertex_3, other%vertex_3)) then
      call self%add_vertex_occurrence( vertex_id=3, facet_id=other%id)
      call other%add_vertex_occurrence(vertex_id=3, facet_id=self%id)
   endif
   contains
      pure function check_pair(a, b)
      !< Check equality of vertices pair.
      type(vector_R8P), intent(in) :: a, b       !< Vertices pair.
      logical                      :: check_pair !< Check result.

      check_pair = ((abs(a%x - b%x) <= EPS).and.(abs(a%y - b%y) <= EPS).and.(abs(a%z - b%z) <= EPS))
      endfunction check_pair
   endsubroutine check_vertices_occurrencies

   elemental subroutine compute_metrix(self)
   !< Compute local (plane) metrix.
   class(facet_object), intent(inout) :: self !< Facet.

   call self%compute_normal

   self%E12 = self%vertex_2 - self%vertex_1
   self%E13 = self%vertex_3 - self%vertex_1
   self%a   = self%E12.dot.self%E12
   self%b   = self%E12.dot.self%E13
   self%c   = self%E13.dot.self%E13
   self%det = self%a * self%c - self%b * self%b

   self%d = self%normal.dot.self%vertex_1

   self%bb(1)%x = min(self%vertex_1%x, self%vertex_2%x, self%vertex_3%x)
   self%bb(1)%y = min(self%vertex_1%y, self%vertex_2%y, self%vertex_3%y)
   self%bb(1)%z = min(self%vertex_1%z, self%vertex_2%z, self%vertex_3%z)
   self%bb(2)%x = max(self%vertex_1%x, self%vertex_2%x, self%vertex_3%x)
   self%bb(2)%y = max(self%vertex_1%y, self%vertex_2%y, self%vertex_3%y)
   self%bb(2)%z = max(self%vertex_1%z, self%vertex_2%z, self%vertex_3%z)
   endsubroutine compute_metrix

   elemental subroutine compute_normal(self)
   !< Compute normal by means of vertices data.
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
   endsubroutine compute_normal

   elemental subroutine destroy(self)
   !< Destroy AABB.
   class(facet_object), intent(inout) :: self  !< Facet.
   type(facet_object)                 :: fresh !< Fresh instance of facet.

   self = fresh
   endsubroutine destroy

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
   d = self%E12.dot.V1P
   e = self%E13.dot.V1P
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

   pure function do_ray_intersect(self, ray_origin, ray_direction) result(intersect)
   !< Return true if facet is intersected by ray from origin and oriented as ray direction vector.
   !<
   !< This based on Moller–Trumbore intersection algorithm.
   !<
   !< @note Facet's metrix must be already computed.
   class(facet_object), intent(in) :: self          !< Facet.
   type(vector_R8P),    intent(in) :: ray_origin    !< Ray origin.
   type(vector_R8P),    intent(in) :: ray_direction !< Ray direction.
   logical                         :: intersect     !< Intersection test result.
   type(vector_R8P)                :: h, s, q       !< Projection vectors.
   real(R8P)                       :: a, f, u, v, t !< Baricentric abscissa.

   intersect = .false.
   h = ray_direction.cross.self%E13
   a = self%E12.dot.h
   if ((a > -EPS).and.(a < EPS)) return
   f = 1._R8P / a
   s = ray_origin - self%vertex_1
   u = f * (s.dot.h)
   if ((u < 0._R8P).or.(u > 1._R8P)) return
   q = s.cross.self%E12
   v = f * ray_direction.dot.q
   if ((v < 0._R8P).or.(u + v > 1._R8P)) return
   t = f * self%E13.dot.q
   if (t > EPS) intersect = .true.
   endfunction do_ray_intersect

   elemental subroutine initialize(self)
   !< Initialize facet.
   class(facet_object), intent(inout) :: self  !< Facet.
   type(facet_object)                 :: fresh !< Fresh instance of facet.

   self = fresh
   endsubroutine initialize

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

   pure subroutine make_normal_consistent(self, edge_dir, other)
   !< Make normal of other facet consistent with self.
   class(facet_object), intent(in)    :: self           !< Facet.
   character(*),        intent(in)    :: edge_dir       !< Edge (in self numeration) along which other is connected.
   type(facet_object),  intent(inout) :: other          !< Other facet to make consistent with self.
   character(len(edge_dir))           :: edge_dir_other !< Edge (in self numeration) along which other is connected.
   type(vector_R8P)                   :: edge           !< Edge of connection in the self reference.
   type(vector_R8P)                   :: edge_other     !< Edge of connection in the other reference.

   call self%edge_connection_in_other_ref(other=other, edge_dir=edge_dir_other, edge=edge_other)
   ! get self edge
   select case(edge_dir)
   case('edge_12')
      edge = self%vertex_2 - self%vertex_1
   case('edge_23')
      edge = self%vertex_3 - self%vertex_2
   case('edge_31')
      edge = self%vertex_1 - self%vertex_3
   endselect
   if (edge%dotproduct(edge_other)>0) then
      ! other numeration is consistent, normal has wrong orientation
      call other%flip_edge(edge_dir=edge_dir_other)
   endif
   endsubroutine make_normal_consistent

   elemental subroutine reverse_normal(self)
   !< Reverse facet normal.
   class(facet_object), intent(inout) :: self   !< Facet.
   type(vector_R8P)                   :: vertex !< Temporary vertex variable.

   call self%flip_edge(edge_dir='edge_23')
   endsubroutine reverse_normal

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

   pure function tetrahedron_volume(self, apex) result(volume)
   !< Return the volume of tetrahedron built by facet and a given apex.
   class(facet_object), intent(in) :: self   !< Facet.
   type(vector_R8P),    intent(in) :: apex   !< Tetrahedron apex.
   real(R8P)                       :: volume !< Tetrahedron volume.
   type(vector_R8P)                :: e12    !< Edge 1-2.
   type(vector_R8P)                :: e13    !< Edge 1-3.

   e12 = self%vertex_2 - self%vertex_1
   e13 = self%vertex_3 - self%vertex_1
   volume = 0.5_R8P * normL2_R8P(e12) * normL2_R8P(e13) * sin(angle_R8P(e12, e13)) * &
            apex%distance_to_plane(pt1=self%vertex_1, pt2=self%vertex_2, pt3=self%vertex_3) / 3._R8P
   endfunction

   pure subroutine update_connectivity(self)
   !< Update facet connectivity.
   !<
   !< @note Vertices occurrencies list must be already computed.
   class(facet_object), intent(inout) :: self !< Facet.

   self%fcon_edge_12 = facet_connected(occurrence_1=self%vertex_1_occurrence, occurrence_2=self%vertex_2_occurrence)
   self%fcon_edge_23 = facet_connected(occurrence_1=self%vertex_2_occurrence, occurrence_2=self%vertex_3_occurrence)
   self%fcon_edge_31 = facet_connected(occurrence_1=self%vertex_3_occurrence, occurrence_2=self%vertex_1_occurrence)
   contains
      pure function facet_connected(occurrence_1, occurrence_2)
      !< Return the facet ID connected by the edge. If no facet is found 0 is returned.
      !<
      !< @note Within two vertices occurrencies, namely one edge, there could be only two connected facets.
      integer(I4P), allocatable, intent(in) :: occurrence_1(:) !< Occurrences list of vertex 1.
      integer(I4P), allocatable, intent(in) :: occurrence_2(:) !< Occurrences list of vertex 2.
      integer(I4P)                          :: facet_connected !< ID of connected connected.
      integer(I4P)                          :: i1, i2          !< Counter.

      facet_connected = 0
      if (allocated(occurrence_1).and.allocated(occurrence_2)) then
         loop_1: do i1=1, size(occurrence_1, dim=1)
            do i2=1, size(occurrence_2, dim=1)
               if (occurrence_1(i1) == occurrence_2(i2)) then
                  facet_connected = occurrence_1(i1)
                  exit loop_1
               endif
            enddo
         enddo loop_1
      endif
      endfunction facet_connected
   endsubroutine update_connectivity

   pure function vertex_global_id(self, vertex_id)
   !< Return the vertex global id given the local one.
   class(facet_object), intent(in) :: self             !< Facet.
   integer(I4P),        intent(in) :: vertex_id        !< Local vertex id.
   integer(I4P)                    :: vertex_global_id !< Gloval vertex id.

   vertex_global_id = (self%id - 1) * 3 + vertex_id
   endfunction vertex_global_id

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
   lhs%E12 = rhs%E12
   lhs%E13 = rhs%E13
   lhs%a = rhs%a
   lhs%b = rhs%b
   lhs%c = rhs%c
   lhs%d = rhs%d
   lhs%det = rhs%det
   lhs%bb = rhs%bb
   lhs%id = rhs%id
   lhs%fcon_edge_12 = rhs%fcon_edge_12
   lhs%fcon_edge_23 = rhs%fcon_edge_23
   lhs%fcon_edge_31 = rhs%fcon_edge_31
   if (allocated(lhs%vertex_1_occurrence)) deallocate(lhs%vertex_1_occurrence)
   if (allocated(rhs%vertex_1_occurrence)) lhs%vertex_1_occurrence = rhs%vertex_1_occurrence
   if (allocated(lhs%vertex_2_occurrence)) deallocate(lhs%vertex_2_occurrence)
   if (allocated(rhs%vertex_2_occurrence)) lhs%vertex_2_occurrence = rhs%vertex_2_occurrence
   if (allocated(lhs%vertex_3_occurrence)) deallocate(lhs%vertex_3_occurrence)
   if (allocated(rhs%vertex_3_occurrence)) lhs%vertex_3_occurrence = rhs%vertex_3_occurrence
   endsubroutine facet_assign_facet

   pure subroutine edge_connection_in_other_ref(self, other, edge_dir, edge)
   !< Return the edge of connection in the other reference.
   class(facet_object), intent(in)  :: self     !< Facet.
   type(facet_object),  intent(in)  :: other    !< Other facet.
   character(*),        intent(out) :: edge_dir !< Edge (in other numeration) along which self is connected.
   type(vector_R8P),    intent(out) :: edge     !< Edge (in other numeration) along which self is connected.

   if     (other%fcon_edge_12 == self%id) then
      edge_dir = 'edge_12'
      edge = other%vertex_2 - other%vertex_1
   elseif (other%fcon_edge_23 == self%id) then
      edge_dir = 'edge_23'
      edge = other%vertex_3 - other%vertex_2
   elseif (other%fcon_edge_31 == self%id) then
      edge_dir = 'edge_31'
      edge = other%vertex_1 - other%vertex_3
   endif
   endsubroutine edge_connection_in_other_ref

   pure subroutine flip_edge(self, edge_dir)
   !< Flip facet edge.
   class(facet_object), intent(inout) :: self     !< Facet.
   character(*),        intent(in)    :: edge_dir !< Edge to be flipped.
   integer(I4P)                       :: fcon     !< Temporary facet connectiviy variable.

   select case(edge_dir)
   case('edge_12')
      call flip_vertices(a=self%vertex_1, b=self%vertex_2,                     &
                         fcon_bc=self%fcon_edge_23, fcon_ca=self%fcon_edge_31, &
                         vertex_a_occurrence=self%vertex_1_occurrence, vertex_b_occurrence=self%vertex_2_occurrence)
   case('edge_23')
      call flip_vertices(a=self%vertex_2, b=self%vertex_3,                     &
                         fcon_bc=self%fcon_edge_12, fcon_ca=self%fcon_edge_31, &
                         vertex_a_occurrence=self%vertex_2_occurrence, vertex_b_occurrence=self%vertex_3_occurrence)
   case('edge_31')
      call flip_vertices(a=self%vertex_3, b=self%vertex_1,                     &
                         fcon_bc=self%fcon_edge_12, fcon_ca=self%fcon_edge_23, &
                         vertex_a_occurrence=self%vertex_3_occurrence, vertex_b_occurrence=self%vertex_1_occurrence)
   endselect
   call self%compute_metrix
   contains
      pure subroutine flip_vertices(a, b, fcon_bc, fcon_ca, vertex_a_occurrence, vertex_b_occurrence)
      !< Flip two vertices of facet.
      type(vector_R8P),          intent(inout) :: a, b                   !< Vertices to be flipped.
      integer(I4P),              intent(inout) :: fcon_bc                !< Connected face ID along edge b-c.
      integer(I4P),              intent(inout) :: fcon_ca                !< Connected face ID along edge c-a.
      integer(I4P), allocatable, intent(inout) :: vertex_a_occurrence(:) !< List of vertex a "occurrencies".
      integer(I4P), allocatable, intent(inout) :: vertex_b_occurrence(:) !< List of vertex b "occurrencies".
      type(vector_R8P)                         :: vertex                 !< Temporary vertex variable.
      integer(I4P)                             :: fcon                   !< Temporary connected face ID.
      integer(I4P), allocatable                :: vertex_occurrence(:)   !< Temporary list of vertex "occurrencies".

      ! flip vertex
      vertex = a
      a = b
      b = vertex
      ! flip facet connectivity
      fcon = fcon_bc
      fcon_bc = fcon_ca
      fcon_ca = fcon
      ! flip vertex occurrences
      if (allocated(vertex_a_occurrence).and.allocated(vertex_a_occurrence)) then
         vertex_occurrence = vertex_a_occurrence
         vertex_a_occurrence = vertex_b_occurrence
         vertex_b_occurrence = vertex_occurrence
      elseif (allocated(vertex_a_occurrence)) then
         vertex_b_occurrence = vertex_a_occurrence
         deallocate(vertex_a_occurrence)
      elseif (allocated(vertex_b_occurrence)) then
         vertex_a_occurrence = vertex_b_occurrence
         deallocate(vertex_b_occurrence)
      endif
      endsubroutine flip_vertices
   endsubroutine flip_edge
endmodule fossil_facet_object
