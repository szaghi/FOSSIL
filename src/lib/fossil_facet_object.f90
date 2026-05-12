!< FOSSIL, facet class definition.

module fossil_facet_object
!< FOSSIL, facet class definition.

use fossil_list_id_object, only : list_id_object
use fossil_utils, only : EPS, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : FR4P, I2P, I4P, R4P, R8P, str, ZeroR8P
use vecfor, only : angle_R8P, face_normal3_R8P, mirror_matrix_R8P, normL2_R8P, rotation_matrix_R8P, vector_R8P

implicit none
private
public :: facet_object
public :: EDGE_12, EDGE_23, EDGE_31

! Edge indices for the connectivity API (fcon_edge, edge_pnormal, make_normal_consistent,
! flip_edge, edge_connection_in_other_ref). Convention:
!   EDGE_12 → vertex(1) → vertex(2)
!   EDGE_23 → vertex(2) → vertex(3)
!   EDGE_31 → vertex(3) → vertex(1)
integer(I4P), parameter :: EDGE_12 = 1_I4P
integer(I4P), parameter :: EDGE_23 = 2_I4P
integer(I4P), parameter :: EDGE_31 = 3_I4P

type :: facet_object
   !< FOSSIL, facet class.
   type(vector_R8P) :: normal    !< Facet (outward) normal (versor), `(v2-v1).cross.(v3-v1)`.
   type(vector_R8P) :: vertex(3) !< Facet vertices (cache of pool%coord(vertex_id) when pool is in use).
   integer(I4P)     :: vertex_id(3) = 0_I4P !< Pool ids for the three vertices (issue #5 stage 3a). 0 = unassigned.
   ! metrix
   type(vector_R8P) :: centroid !< Facet's centroid.
   ! triangle plane parametric equation: T(s,t) = B + s*E12 + t*E13
   type(vector_R8P) :: E12        !< Edge 1-2, `V2-V1`.
   type(vector_R8P) :: E13        !< Edge 1-3, `V3-V1`.
   real(R8P)        :: a=0._R8P   !< `E12.dot.E12`.
   real(R8P)        :: b=0._R8P   !< `E12.dot.E13`.
   real(R8P)        :: c=0._R8P   !< `E13.dot.E13`.
   real(R8P)        :: det=0._R8P !< Gram det of (E12, E13): `a*c - b*b` = |E12 x E13|^2 = (2*area)^2. Use `area()` for the facet area.
   ! triangle plane equation: nx*x + ny*y + nz*z - d = 0, normal == [nx, ny, nz]
   real(R8P) :: d=0._R8P !< `normal.dot.vertex(1)`
   ! auxiliary
   type(vector_R8P) :: bb(2) !< Axis-aligned bounding box (AABB), bb(1)=min, bb(2)=max.
   ! connectivity — edges are indexed 1..3 with the convention:
   !   edge 1 = vertex(1)→vertex(2)
   !   edge 2 = vertex(2)→vertex(3)
   !   edge 3 = vertex(3)→vertex(1)
   ! Named parameters EDGE_12, EDGE_23, EDGE_31 are exported for readability.
   integer(I4P)         :: id                   !< Facet global ID.
   integer(I4P)         :: fcon_edge(3)=0_I4P   !< Connected face ID along each edge (0 = disconnected).
   type(list_id_object) :: vertex_nearby(3)     !< List of vertices "nearby", list of vertices global ID nearby them.
   ! pseudo normals
   type(vector_R8P) :: edge_pnormal(3)   !< Edge pseudo-normals, indexed as fcon_edge.
   type(vector_R8P) :: vertex_pnormal(3) !< Vertices pseudo-normals.
   contains
      ! public methods
      procedure, pass(self) :: centroid_part                   !< Return facet's part to build up STL centroid.
      procedure, pass(self) :: check_normal                    !< Check normal consistency.
      procedure, pass(self) :: compute_distance                !< Compute the (unsigned, squared) distance from a point to facet.
      procedure, pass(self) :: compute_distance_with_region    !< Same as compute_distance, plus closest point and Voronoi region.
      procedure, pass(self) :: pseudo_normal_for_region        !< Return the pseudo-normal for a given closest-point region.
      procedure, pass(self) :: compute_metrix                  !< Compute local (plane) metrix.
      procedure, pass(self) :: compute_normal                  !< Compute normal by means of vertices data.
      procedure, pass(self) :: compute_vertices_nearby         !< Compute vertices nearby comparing to ones of other facet.
      procedure, pass(self) :: destroy                         !< Destroy facet.
      procedure, pass(self) :: destroy_connectivity            !< Destroy facet connectivity.
      procedure, pass(self) :: do_ray_intersect                !< Return true if facet is intersected by a ray.
      procedure, pass(self) :: initialize                      !< Initialize facet.
      procedure, pass(self) :: largest_edge_len                !< Return the largest edge length.
      procedure, pass(self) :: load_from_file_ascii            !< Load facet from ASCII file.
      procedure, pass(self) :: load_from_file_binary           !< Load facet from binary file.
      procedure, pass(self) :: make_normal_consistent          !< Make normal of other facet consistent with self.
      generic               :: mirror => mirror_by_normal, &
                                         mirror_by_matrix      !< Mirror facet.
      procedure, pass(self) :: reverse_normal                  !< Reverse facet normal.
      procedure, pass(self) :: resize                          !< Resize (scale) facet by x or y or z or vectorial factors.
      generic               :: rotate => rotate_by_axis_angle, &
                                         rotate_by_matrix      !< Rotate facet.
      procedure, pass(self) :: save_into_file_ascii            !< Save facet into ASCII file.
      procedure, pass(self) :: save_into_file_binary           !< Save facet into binary file.
      procedure, pass(self) :: smallest_edge_len               !< Return the smallest edge length.
      procedure, pass(self) :: solid_angle                     !< Return the (projected) solid angle of the facet with respect point.
      procedure, pass(self) :: tetrahedron_volume              !< Return the volume of tetrahedron built by facet and a given apex.
      procedure, pass(self) :: translate                       !< Translate facet given vectorial delta.
      procedure, pass(self) :: area                            !< Return the facet area (issue #7).
      procedure, pass(self) :: vertex_angle                    !< Return the subtended angle of given vertex.
      procedure, pass(self) :: vertex_global_id                !< Return the vertex global id given the local one.
      procedure, pass(self) :: set_vertex_ids                  !< Assign the three pool ids for this facet (issue #5 stage 3a).
      ! private methods
      procedure, pass(self), private :: edge_connection_in_other_ref !< Return the edge of connection in the other reference.
      procedure, pass(self), private :: flip_edge                    !< Flip facet edge.
      procedure, pass(self), private :: mirror_by_normal             !< Mirror facet given normal of mirroring plane.
      procedure, pass(self), private :: mirror_by_matrix             !< Mirror facet given matrix.
      procedure, pass(self), private :: rotate_by_axis_angle         !< Rotate facet given axis and angle.
      procedure, pass(self), private :: rotate_by_matrix             !< Rotate facet given matrix.
endtype facet_object

contains
   ! public methods
   pure function centroid_part(self)
   !< Return facet's part to build up STL centroid.
   !<
   !< @note Facet's normal should already computed/sanitized.
   class(facet_object), intent(in)  :: self          !< Facet.
   type(vector_R8P)                 :: centroid_part !< Facet's part of the STL centroid.

   associate(normal=>self%normal, vertex=>self%vertex)
      centroid_part%x = normal%x * ((vertex(1)%x + vertex(2)%x) * (vertex(1)%x + vertex(2)%x) + &
                                    (vertex(2)%x + vertex(3)%x) * (vertex(2)%x + vertex(3)%x) + &
                                    (vertex(3)%x + vertex(1)%x) * (vertex(3)%x + vertex(1)%x))
      centroid_part%y = normal%y * ((vertex(1)%y + vertex(2)%y) * (vertex(1)%y + vertex(2)%y) + &
                                    (vertex(2)%y + vertex(3)%y) * (vertex(2)%y + vertex(3)%y) + &
                                    (vertex(3)%y + vertex(1)%y) * (vertex(3)%y + vertex(1)%y))
      centroid_part%z = normal%z * ((vertex(1)%z + vertex(2)%z) * (vertex(1)%z + vertex(2)%z) + &
                                    (vertex(2)%z + vertex(3)%z) * (vertex(2)%z + vertex(3)%z) + &
                                    (vertex(3)%z + vertex(1)%z) * (vertex(3)%z + vertex(1)%z))
   endassociate
   endfunction centroid_part

   elemental function check_normal(self) result(is_consistent)
   !< Check normal consistency.
   class(facet_object), intent(in) :: self          !< Facet.
   logical                         :: is_consistent !< Consistency check result.
   type(vector_R8P)                :: normal        !< Normal computed by means of vertices data.

   normal = face_normal3_R8P(pt1=self%vertex(1), pt2=self%vertex(2), pt3=self%vertex(3), norm='y')
   is_consistent = ((abs(normal%x - self%normal%x)<=2*ZeroR8P).and.&
                    (abs(normal%y - self%normal%y)<=2*ZeroR8P).and.&
                    (abs(normal%z - self%normal%z)<=2*ZeroR8P))
   endfunction check_normal

   pure subroutine compute_distance(self, point, distance)
   !< Compute the (unsigned, squared) distance from a point to the facet surface.
   !<
   !< Thin wrapper around `compute_distance_with_region` that discards the
   !< closest-point and region outputs — kept for callers that only want d^2.
   class(facet_object), intent(in)  :: self      !< Facet.
   type(vector_R8P),    intent(in)  :: point     !< Point.
   real(R8P),           intent(out) :: distance  !< Closest squared distance from point to the facet.
   type(vector_R8P)                 :: closest   !< Discarded.
   integer(I4P)                     :: region    !< Discarded.

   call self%compute_distance_with_region(point=point, distance=distance, closest=closest, region=region)
   endsubroutine compute_distance

   pure subroutine compute_distance_with_region(self, point, distance, closest, region)
   !< Compute squared distance from a point to the facet, the closest point on the
   !< facet, and a tag identifying which Voronoi region of the triangle contains
   !< the closest point (face / one of three edges / one of three vertices).
   !<
   !< The closest point and region are needed by the signed-distance pseudo-normal
   !< sign test: the correct pseudo-normal at the closest point depends on whether
   !< it lies on the interior of a face, an edge, or a vertex.
   !<
   !< Region encoding:
   !<    0          : interior (face region) — use `self%normal`
   !<    1, 2, 3    : edge EDGE_12, EDGE_23, EDGE_31 respectively — use `self%edge_pnormal(k)`
   !<   -1, -2, -3  : vertex 1, 2, 3 respectively — use `self%vertex_pnormal(|k|)`
   !<
   !< @note Facet's metrix must be already computed.
   !<
   !< @note Algorithm by David Eberly, Geometric Tools LLC, http://www.geometrictools.com.
   class(facet_object), intent(in)  :: self                             !< Facet.
   type(vector_R8P),    intent(in)  :: point                            !< Point.
   real(R8P),           intent(out) :: distance                         !< Closest squared distance from point to the facet.
   type(vector_R8P),    intent(out) :: closest                          !< Closest point on the facet.
   integer(I4P),        intent(out) :: region                           !< Voronoi region tag (see encoding above).
   type(vector_R8P)                 :: V1P                              !< `vertex(1)-point`.
   real(R8P)                        :: d, e, f, s, t, sq, tq            !< Plane equation coefficients.
   real(R8P)                        :: tmp0, tmp1, numer, denom, invdet !< Temporary.
   real(R8P), parameter             :: BARY_TOL = 1.0e-12_R8P           !< Tolerance for classifying barycentric coords as 0/1.

   associate(a=>self%a, b=>self%b, c=>self%c, det=>self%det)
   V1P = self%vertex(1) - point
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
   closest = self%vertex(1) + sq * self%E12 + tq * self%E13
   ! Region classification from (sq, tq) barycentric-style coordinates.
   ! Vertex correspondence: V1 ↔ (0,0), V2 ↔ (1,0), V3 ↔ (0,1).
   ! Edge correspondence:   EDGE_12 (V1→V2) ↔ tq=0; EDGE_23 (V2→V3) ↔ sq+tq=1; EDGE_31 (V3→V1) ↔ sq=0.
   if (sq < BARY_TOL .and. tq < BARY_TOL) then
      region = -1_I4P                            ! vertex 1
   elseif (sq > 1._R8P - BARY_TOL .and. tq < BARY_TOL) then
      region = -2_I4P                            ! vertex 2
   elseif (tq > 1._R8P - BARY_TOL .and. sq < BARY_TOL) then
      region = -3_I4P                            ! vertex 3
   elseif (tq < BARY_TOL) then
      region = EDGE_12                           ! edge V1-V2
   elseif (sq < BARY_TOL) then
      region = EDGE_31                           ! edge V3-V1
   elseif (sq + tq > 1._R8P - BARY_TOL) then
      region = EDGE_23                           ! edge V2-V3
   else
      region = 0_I4P                             ! face interior
   endif
   endassociate
   endsubroutine compute_distance_with_region

   pure function pseudo_normal_for_region(self, region) result(n)
   !< Return the pseudo-normal of `self` corresponding to a closest-point region tag
   !< produced by `compute_distance_with_region`. See that routine for the encoding.
   class(facet_object), intent(in) :: self    !< Facet.
   integer(I4P),        intent(in) :: region  !< Region tag.
   type(vector_R8P)                :: n       !< Pseudo-normal at the closest point.

   select case (region)
   case (0_I4P)         ; n = self%normal
   case (1_I4P)         ; n = self%edge_pnormal(1)   ! EDGE_12
   case (2_I4P)         ; n = self%edge_pnormal(2)   ! EDGE_23
   case (3_I4P)         ; n = self%edge_pnormal(3)   ! EDGE_31
   case (-1_I4P)        ; n = self%vertex_pnormal(1)
   case (-2_I4P)        ; n = self%vertex_pnormal(2)
   case (-3_I4P)        ; n = self%vertex_pnormal(3)
   case default         ; n = self%normal             ! safety fallback
   end select
   endfunction pseudo_normal_for_region

   elemental subroutine compute_metrix(self)
   !< Compute local (plane) metrix.
   class(facet_object), intent(inout) :: self !< Facet.

   self%centroid = (self%vertex(1) + self%vertex(2) + self%vertex(3)) / 3._R8P
   call self%compute_normal

   self%E12 = self%vertex(2) - self%vertex(1)
   self%E13 = self%vertex(3) - self%vertex(1)
   self%a   = self%E12%dotproduct(rhs=self%E12)
   self%b   = self%E12%dotproduct(rhs=self%E13)
   self%c   = self%E13%dotproduct(rhs=self%E13)
   ! self%a   = self%E12.dot.self%E12
   ! self%b   = self%E12.dot.self%E13
   ! self%c   = self%E13.dot.self%E13
   self%det = self%a * self%c - self%b * self%b

   ! self%d = self%normal.dot.self%vertex(1)
   self%d = self%normal%dotproduct(rhs=self%vertex(1))

   self%bb(1)%x = min(self%vertex(1)%x, self%vertex(2)%x, self%vertex(3)%x)
   self%bb(1)%y = min(self%vertex(1)%y, self%vertex(2)%y, self%vertex(3)%y)
   self%bb(1)%z = min(self%vertex(1)%z, self%vertex(2)%z, self%vertex(3)%z)
   self%bb(2)%x = max(self%vertex(1)%x, self%vertex(2)%x, self%vertex(3)%x)
   self%bb(2)%y = max(self%vertex(1)%y, self%vertex(2)%y, self%vertex(3)%y)
   self%bb(2)%z = max(self%vertex(1)%z, self%vertex(2)%z, self%vertex(3)%z)
   endsubroutine compute_metrix

   elemental subroutine compute_normal(self)
   !< Compute normal by means of vertices data.
   !<
   !<```fortran
   !< type(facet_object) :: facet
   !< facet%vertex(1) = -0.231369_R4P * ex_R4P + 0.0226865_R4P * ey_R4P + 1._R4P * ez_R4P
   !< facet%vertex(2) = -0.227740_R4P * ex_R4P + 0.0245457_R4P * ey_R4P + 0._R4P * ez_R4P
   !< facet%vertex(2) = -0.235254_R4P * ex_R4P + 0.0201881_R4P * ey_R4P + 0._R4P * ez_R4P
   !< call facet%sanitize_normal
   !< print "(3(F3.1,1X))", facet%normal%x, facet%normal%y, facet%normal%z
   !<```
   !=> -0.501673222 0.865057290 -2.12257713<<<
   class(facet_object), intent(inout) :: self !< Facet.

   self%normal = face_normal3_R8P(pt1=self%vertex(1), pt2=self%vertex(2), pt3=self%vertex(3), norm='y')
   endsubroutine compute_normal

   pure subroutine compute_vertices_nearby(self, other, tolerance_to_be_nearby)
   !< Populate `vertex_nearby` (loose-tolerance) for each vertex of `self` and `other`.
   !<
   !< Strict-EPS coincidence is owned by `vertex_pool_object` (issue #5 stage 3c).
   !< This routine only records vertices within the looser sanitize tolerance, used
   !< downstream by `connect_nearby_vertices` to snap nearby vertices together.
   class(facet_object), intent(inout) :: self                   !< Facet.
   type(facet_object),  intent(inout) :: other                  !< Other facet.
   real(R8P),           intent(in)    :: tolerance_to_be_nearby !< Tolerance to identify nearby vertices.
   integer(I4P)                       :: vs, vo                 !< Counter.

   do vs=1, 3
      do vo=1, 3
         if (are_nearby(self%vertex(vs), other%vertex(vo), tolerance_to_be_nearby)) then
            call  self%vertex_nearby(vs)%put(id=other%vertex_global_id(vo))
            call other%vertex_nearby(vo)%put(id= self%vertex_global_id(vs))
         endif
      enddo
   enddo
   contains
      pure function are_nearby(a, b, tolerance)
      !< Check equality of vertices pair.
      type(vector_R8P), intent(in) :: a, b       !< Vertices pair.
      real(R8P),        intent(in) :: tolerance  !< Check tolerance.
      logical                      :: are_nearby !< Check result.

      are_nearby = ((abs(a%x - b%x) <= tolerance).and.&
                    (abs(a%y - b%y) <= tolerance).and.&
                    (abs(a%z - b%z) <= tolerance))
      endfunction are_nearby
   endsubroutine compute_vertices_nearby

   elemental subroutine destroy(self)
   !< Destroy facet.
   class(facet_object), intent(inout) :: self  !< Facet.

   self%normal         = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%vertex         = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%centroid       = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%E12            = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%E13            = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%a              = 0._R8P
   self%b              = 0._R8P
   self%c              = 0._R8P
   self%d              = 0._R8P
   self%det            = 0._R8P
   self%bb             = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%id             = 0_I4P
   self%fcon_edge      = 0_I4P
   self%edge_pnormal   = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%vertex_pnormal = vector_R8P(0._R8P, 0._R8P, 0._R8P)
   self%vertex_id      = 0_I4P
   call self%vertex_nearby%destroy
   endsubroutine destroy

   elemental subroutine destroy_connectivity(self)
   !< Destroy facet connectivity.
   class(facet_object), intent(inout) :: self  !< Facet.

   self%fcon_edge = 0_I4P
   call self%vertex_nearby%destroy
   endsubroutine destroy_connectivity

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
   ! h = ray_direction.cross.self%E13
   h = ray_direction%crossproduct(rhs=self%E13)
   ! a = self%E12.dot.h
   a = self%E12%dotproduct(rhs=h)
   if ((a > -EPS).and.(a < EPS)) return
   f = 1._R8P / a
   s = ray_origin - self%vertex(1)
   ! u = f * (s.dot.h)
   u = f * (s%dotproduct(rhs=h))
   if ((u < 0._R8P).or.(u > 1._R8P)) return
   ! q = s.cross.self%E12
   q = s%crossproduct(rhs=self%E12)
   v = f * ray_direction.dot.q
   if ((v < 0._R8P).or.(u + v > 1._R8P)) return
   t = f * self%E13.dot.q
   if (t > EPS) intersect = .true.
   endfunction do_ray_intersect

   elemental subroutine initialize(self)
   !< Initialize facet.
   class(facet_object), intent(inout) :: self  !< Facet.

   call self%destroy
   endsubroutine initialize

   pure function largest_edge_len(self) result(largest)
   !< Return the largest edge length.
   class(facet_object), intent(in) :: self    !< Facet.
   real(R8P)                       :: largest !< largest edge length.

   largest = max(normL2_R8P(self%vertex(2)-self%vertex(1)), &
                 normL2_R8P(self%vertex(3)-self%vertex(2)), &
                 normL2_R8P(self%vertex(1)-self%vertex(3)))
   endfunction largest_edge_len

   subroutine load_from_file_ascii(self, file_unit)
   !< Load facet from ASCII file.
   !<
   !< `iostat=` traps premature EOF / I/O errors. A corrupted STL (truncated mid-facet,
   !< header claiming more facets than the body contains, mixed line endings) now
   !< `error stop`s with the offending unit and operation instead of producing a generic
   !< runtime error with no context.
   class(facet_object), intent(inout) :: self      !< Facet.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I4P)                       :: ios       !< I/O status.

   call load_facet_record(prefix='facet normal', record=self%normal)
   read(file_unit, *, iostat=ios) ! outer loop
   call check_ios(ios, 'reading "outer loop" delimiter')
   call load_facet_record(prefix='vertex', record=self%vertex(1))
   call load_facet_record(prefix='vertex', record=self%vertex(2))
   call load_facet_record(prefix='vertex', record=self%vertex(3))
   read(file_unit, *, iostat=ios) ! endloop
   call check_ios(ios, 'reading "endloop" delimiter')
   read(file_unit, *, iostat=ios) ! endfacet
   call check_ios(ios, 'reading "endfacet" delimiter')
   contains
      subroutine check_ios(stat, op)
      !< Trap non-zero I/O status. Negative = EOF, positive = error.
      integer(I4P), intent(in) :: stat !< Status from a read.
      character(*), intent(in) :: op   !< Description of what was being read.
      if (stat == 0) return
      if (stat < 0) then
         error stop 'fossil_facet_object%load_from_file_ascii: unexpected end of file while '//op
      else
         error stop 'fossil_facet_object%load_from_file_ascii: I/O error while '//op
      endif
      endsubroutine check_ios

      subroutine load_facet_record(prefix, record)
      !< Load a facet *record*, namely normal or vertex data.
      !<
      !< Strict prefix check: the line must, after leading whitespace, *start with* the
      !< prefix. A previous substring-`index` check would accept any line containing the
      !< keyword anywhere (e.g. a stray "vertex count: 12" comment would be parsed as a
      !< coordinate triple), producing silent garbage.
      character(*),     intent(in)  :: prefix         !< Record prefix string.
      type(vector_R8P), intent(out) :: record         !< Record data.
      character(FRLEN)              :: facet_record   !< Facet record string buffer.
      character(len=:), allocatable :: trimmed        !< Line with leading whitespace removed.
      integer(I4P)                  :: prefix_len     !< Cached prefix length.
      integer(I4P)                  :: read_ios       !< I/O status for the line read.
      integer(I4P)                  :: parse_ios      !< I/O status for the coordinate parse.

      read(file_unit, '(A)', iostat=read_ios) facet_record
      if (read_ios /= 0) then
         if (read_ios < 0) then
            error stop 'fossil_facet_object%load_from_file_ascii: unexpected end of file before "'//prefix//'"'
         else
            error stop 'fossil_facet_object%load_from_file_ascii: I/O error reading "'//prefix//'" line'
         endif
      endif
      trimmed    = trim(adjustl(facet_record))
      prefix_len = len(prefix)
      if (len(trimmed) >= prefix_len) then
         if (trimmed(1:prefix_len) == prefix) then
            read(trimmed(prefix_len + 1:), *, iostat=parse_ios) record%x, record%y, record%z
            if (parse_ios /= 0) &
               error stop 'fossil_facet_object%load_from_file_ascii: bad coordinates after "'//prefix//'"'
            return
         endif
      endif
      write(stderr, '(A)') 'error: expected line to start with "'//prefix// &
                           '" on unit '//trim(str(file_unit))//', got: "'//trim(facet_record)//'"'
      error stop 'fossil_facet_object%load_from_file_ascii: malformed STL'
      endsubroutine load_facet_record
   endsubroutine load_from_file_ascii

   subroutine load_from_file_binary(self, file_unit)
   !< Load facet from binary file.
   !<
   !< `iostat=` traps premature EOF. A corrupted binary STL whose header claims more
   !< facets than the file actually contains (or that is truncated mid-facet) now
   !< `error stop`s with a specific message instead of producing a generic Fortran
   !< stream-read error.
   class(facet_object), intent(inout) :: self       !< Facet.
   integer(I4P),        intent(in)    :: file_unit  !< File unit.
   integer(I2P)                       :: padding    !< Facet padding.
   real(R4P)                          :: triplet(3) !< Triplet record of R4P kind real.
   integer(I4P)                       :: ios        !< I/O status.

   read(file_unit, iostat=ios) triplet
   call check_ios(ios, 'reading facet normal')
   self%normal%x=real(triplet(1), R8P) ; self%normal%y=real(triplet(2), R8P) ; self%normal%z=real(triplet(3), R8P)
   read(file_unit, iostat=ios) triplet
   call check_ios(ios, 'reading vertex 1')
   self%vertex(1)%x=real(triplet(1), R8P) ; self%vertex(1)%y=real(triplet(2), R8P) ; self%vertex(1)%z=real(triplet(3), R8P)
   read(file_unit, iostat=ios) triplet
   call check_ios(ios, 'reading vertex 2')
   self%vertex(2)%x=real(triplet(1), R8P) ; self%vertex(2)%y=real(triplet(2), R8P) ; self%vertex(2)%z=real(triplet(3), R8P)
   read(file_unit, iostat=ios) triplet
   call check_ios(ios, 'reading vertex 3')
   self%vertex(3)%x=real(triplet(1), R8P) ; self%vertex(3)%y=real(triplet(2), R8P) ; self%vertex(3)%z=real(triplet(3), R8P)
   read(file_unit, iostat=ios) padding
   call check_ios(ios, 'reading facet attribute padding')
   contains
      subroutine check_ios(stat, op)
      !< Trap non-zero I/O status. Negative = EOF, positive = error.
      integer(I4P), intent(in) :: stat !< Status from a read.
      character(*), intent(in) :: op   !< Description of what was being read.
      if (stat == 0) return
      if (stat < 0) then
         error stop 'fossil_facet_object%load_from_file_binary: unexpected end of file while '//op
      else
         error stop 'fossil_facet_object%load_from_file_binary: I/O error while '//op
      endif
      endsubroutine check_ios
   endsubroutine load_from_file_binary

   pure subroutine make_normal_consistent(self, edge, other)
   !< Make normal of other facet consistent with self.
   !<
   !< `edge` is the edge index in self numeration (use EDGE_12/EDGE_23/EDGE_31, values 1..3).
   class(facet_object), intent(in)    :: self        !< Facet.
   integer(I4P),        intent(in)    :: edge        !< Edge in self numeration (1..3).
   type(facet_object),  intent(inout) :: other       !< Other facet to make consistent with self.
   integer(I4P)                       :: edge_other  !< Edge in other numeration (1..3).
   type(vector_R8P)                   :: e_self      !< Edge vector in self reference.
   type(vector_R8P)                   :: e_other     !< Edge vector in other reference.

   if (edge < 1 .or. edge > 3) &
      error stop 'fossil_facet_object%make_normal_consistent: edge must be 1..3 (EDGE_12, EDGE_23, EDGE_31)'

   call self%edge_connection_in_other_ref(other=other, edge=edge_other, edge_vector=e_other)
   ! self edge vector: vertex(next) - vertex(edge), with next = mod(edge,3)+1
   e_self = self%vertex(mod(edge, 3) + 1) - self%vertex(edge)
   if (e_self%dotproduct(e_other) > 0) then
      ! other numeration is consistent, normal has wrong orientation
      call other%flip_edge(edge=edge_other)
   endif
   endsubroutine make_normal_consistent

   elemental subroutine resize(self, factor, center)
   !< Resize (scale) facet by x or y or z or vectorial factors.
   !<
   !< @note The name `scale` has not been used, it been a Fortran built-in.
   class(facet_object), intent(inout) :: self   !< Facet
   type(vector_R8P),    intent(in)    :: factor !< Vectorial factor.
   type(vector_R8P),    intent(in)    :: center !< Center of resize.

   self%vertex(1) = (self%vertex(1) - center) * factor + center
   self%vertex(2) = (self%vertex(2) - center) * factor + center
   self%vertex(3) = (self%vertex(3) - center) * factor + center
   endsubroutine resize

   elemental subroutine reverse_normal(self)
   !< Reverse facet normal.
   class(facet_object), intent(inout) :: self   !< Facet.
   type(vector_R8P)                   :: vertex !< Temporary vertex variable.

   call self%flip_edge(edge=EDGE_23)
   endsubroutine reverse_normal

   subroutine save_into_file_ascii(self, file_unit)
   !< Save facet into ASCII file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '  facet normal ', self%normal%x, ' ', self%normal%y, ' ', self%normal%z
   write(file_unit, '(A)')                            '    outer loop'
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex(1)%x, ' ', self%vertex(1)%y, ' ',self%vertex(1)%z
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex(2)%x, ' ', self%vertex(2)%y, ' ',self%vertex(2)%z
   write(file_unit, '(A,2('//FR4P//',A),'//FR4P//')') '      vertex ', self%vertex(3)%x, ' ', self%vertex(3)%y, ' ',self%vertex(3)%z
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
   triplet(1) = real(self%vertex(1)%x, R4P) ; triplet(2) = real(self%vertex(1)%y, R4P) ; triplet(3) = real(self%vertex(1)%z, R4P)
   write(file_unit) triplet
   triplet(1) = real(self%vertex(2)%x, R4P) ; triplet(2) = real(self%vertex(2)%y, R4P) ; triplet(3) = real(self%vertex(2)%z, R4P)
   write(file_unit) triplet
   triplet(1) = real(self%vertex(3)%x, R4P) ; triplet(2) = real(self%vertex(3)%y, R4P) ; triplet(3) = real(self%vertex(3)%z, R4P)
   write(file_unit) triplet
   write(file_unit) 0_I2P
   endsubroutine save_into_file_binary

   pure function smallest_edge_len(self) result(smallest)
   !< Return the smallest edge length.
   class(facet_object), intent(in) :: self     !< Facet.
   real(R8P)                       :: smallest !< Smallest edge length.

   smallest = min(normL2_R8P(self%vertex(2)-self%vertex(1)), &
                  normL2_R8P(self%vertex(3)-self%vertex(2)), &
                  normL2_R8P(self%vertex(1)-self%vertex(3)))
   endfunction smallest_edge_len

   pure function solid_angle(self, point)
   !< Return the (projected) solid angle of the facet with respect the point.
   class(facet_object), intent(in) :: self                      !< Facet.
   type(vector_R8P),    intent(in) :: point                     !< Point.
   real(R8P)                       :: solid_angle               !< Solid angle.
   type(vector_R8P)                :: R1, R2, R3                !< Edges from point to facet vertices.
   real(R8P)                       :: R1_norm, R2_norm, R3_norm !< Norms (L2) of edges from point to facet vertices.
   real(R8P)                       :: numerator                 !< Archtangent numerator.
   real(R8P)                       :: denominator               !< Archtangent denominator.

   R1 = self%vertex(1) - point ; R1_norm = R1%normL2()
   R2 = self%vertex(2) - point ; R2_norm = R2%normL2()
   R3 = self%vertex(3) - point ; R3_norm = R3%normL2()

   ! numerator = R1.dot.(R2.cross.R3)
   numerator = R1%dotproduct(rhs=R2%crossproduct(rhs=R3))
   ! denominator = R1_norm * R2_norm * R3_norm + (R1.dot.R2) * R3_norm + &
   !                                             (R1.dot.R3) * R2_norm + &
   !                                             (R2.dot.R3) * R1_norm
   denominator = R1_norm * R2_norm * R3_norm + (R1%dotproduct(rhs=R2)) * R3_norm + &
                                               (R1%dotproduct(rhs=R3)) * R2_norm + &
                                               (R2%dotproduct(rhs=R3)) * R1_norm

   solid_angle = 2._R8P * atan2(numerator, denominator)
   endfunction solid_angle

   pure function tetrahedron_volume(self, apex) result(volume)
   !< Return the signed volume of the tetrahedron (apex, v1, v2, v3) using the
   !< standard divergence-theorem convention: positive when the facet winding is
   !< outward-pointing relative to the apex.
   !<
   !< Derivation: the textbook signed tetrahedron volume is
   !<   V = (1/6) * (v1-apex) . ((v2-apex) x (v3-apex))
   !<     = (1/6) * ((v2-v1) x (v3-v1)) . (v1-apex)
   !<     = -(1/3) * area * (n . (apex - v1))      [n = unit winding normal]
   !< so summing over a closed surface with apex on (or inside) the surface gives
   !< a positive total iff windings are outward — the convention every other piece
   !< of mesh-processing code uses.
   class(facet_object), intent(in) :: self   !< Facet.
   type(vector_R8P),    intent(in) :: apex   !< Tetrahedron apex.
   real(R8P)                       :: volume !< Tetrahedron signed volume.
   type(vector_R8P)                :: e12    !< Edge 1-2.
   type(vector_R8P)                :: e13    !< Edge 1-3.

   e12 = self%vertex(2) - self%vertex(1)
   e13 = self%vertex(3) - self%vertex(1)
   volume = -0.5_R8P * normL2_R8P(e12) * normL2_R8P(e13) * sin(angle_R8P(e12, e13)) * &
            apex%distance_to_plane(pt1=self%vertex(1), pt2=self%vertex(2), pt3=self%vertex(3)) / 3._R8P
   endfunction tetrahedron_volume

   elemental subroutine translate(self, delta, recompute_metrix)
   !< Translate facet given vectorial delta.
   class(facet_object), intent(inout)        :: self             !< Facet.
   type(vector_R8P),    intent(in)           :: delta            !< Translation delta.
   logical,             intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.

   self%vertex(1) = self%vertex(1) + delta
   self%vertex(2) = self%vertex(2) + delta
   self%vertex(3) = self%vertex(3) + delta
   if (present(recompute_metrix)) then
      if (recompute_metrix) call self%compute_metrix
   endif
   endsubroutine translate

   pure function area(self) result(a)
   !< Return the facet area (issue #7).
   !<
   !< Area = (1/2) * |E12 x E13| = sqrt(self%det) / 2. Cheap closed form using
   !< the Gram-determinant cache already populated by `compute_metrix`. Returns
   !< zero if metrix has not been computed (`det` is zero-initialized).
   class(facet_object), intent(in) :: self !< Facet.
   real(R8P)                       :: a    !< Facet area.

   a = 0.5_R8P * sqrt(self%det)
   endfunction area

   pure function vertex_angle(self, vertex_id)
   !< Return the subtened angle of given vertex.
   class(facet_object), intent(in) :: self         !< Facet.
   integer(I4P),        intent(in) :: vertex_id    !< Local vertex id.
   real(R8P)                       :: vertex_angle !< Subtended angle.
   type(vector_R8P)                :: edge(2)      !< Two edge subtending the vertex.

   select case(vertex_id)
   case(1_I4P)
     edge(1) = self%vertex(2) - self%vertex(1)
     edge(2) = self%vertex(3) - self%vertex(1)
   case(2_I4P)
     edge(1) = self%vertex(3) - self%vertex(2)
     edge(2) = self%vertex(1) - self%vertex(2)
   case(3_I4P)
     edge(1) = self%vertex(1) - self%vertex(3)
     edge(2) = self%vertex(2) - self%vertex(3)
   endselect
   vertex_angle = edge(1)%angle(edge(2))
   endfunction vertex_angle

   pure function vertex_global_id(self, vertex_id)
   !< Return the vertex global id given the local one.
   class(facet_object), intent(in) :: self             !< Facet.
   integer(I4P),        intent(in) :: vertex_id        !< Local vertex id.
   integer(I4P)                    :: vertex_global_id !< Global vertex id.

   vertex_global_id = (self%id - 1) * 3 + vertex_id
   endfunction vertex_global_id

   pure subroutine set_vertex_ids(self, vid1, vid2, vid3)
   !< Assign the three pool ids for this facet (issue #5 stage 3a).
   !<
   !< Called by surface_stl_object after vertex_pool%initialize_from_facets has
   !< assigned ids. The id slots are the structural source of truth from stage 3b
   !< onward; `self%vertex(:)` remains as a coordinate cache for hot kernels.
   class(facet_object), intent(inout) :: self
   integer(I4P),        intent(in)    :: vid1, vid2, vid3

   self%vertex_id(1) = vid1
   self%vertex_id(2) = vid2
   self%vertex_id(3) = vid3
   endsubroutine set_vertex_ids

   ! private methods
   pure subroutine flip_edge(self, edge)
   !< Flip facet edge.
   !<
   !< Flipping edge `e` (1..3) swaps its two endpoint vertices and swaps the connectivity
   !< of the two other edges. The table BC_OF/CA_OF gives the indices of those two other
   !< edges (the two values in {1,2,3}\{e}, in ascending order).
   class(facet_object), intent(inout) :: self  !< Facet.
   integer(I4P),        intent(in)    :: edge  !< Edge to be flipped (1..3).
   integer(I4P)                       :: v1, v2, bc, ca !< Vertex indices and the two non-flipped edge indices.
   integer(I4P)                       :: tmp_vid        !< Temporary for pool id swap.
   integer(I4P), parameter            :: BC_OF(3) = [2_I4P, 1_I4P, 1_I4P]
   integer(I4P), parameter            :: CA_OF(3) = [3_I4P, 3_I4P, 2_I4P]

   if (edge < 1 .or. edge > 3) &
      error stop 'fossil_facet_object%flip_edge: edge must be 1..3 (EDGE_12, EDGE_23, EDGE_31)'

   v1 = edge
   v2 = mod(edge, 3) + 1
   bc = BC_OF(edge)
   ca = CA_OF(edge)

   ! Swap the two endpoint vertices: coordinates AND pool ids in lockstep so the
   ! cache stays consistent with the pool's facet_vid mapping (issue #5 stage 3a).
   block
      type(vector_R8P) :: tmp_v
      integer(I4P)     :: tmp_fcon
      tmp_v             = self%vertex(v1)
      self%vertex(v1)   = self%vertex(v2)
      self%vertex(v2)   = tmp_v
      tmp_vid           = self%vertex_id(v1)
      self%vertex_id(v1) = self%vertex_id(v2)
      self%vertex_id(v2) = tmp_vid
      ! Swap connectivity of the two non-flipped edges so neighbour links remain valid.
      tmp_fcon            = self%fcon_edge(bc)
      self%fcon_edge(bc)  = self%fcon_edge(ca)
      self%fcon_edge(ca)  = tmp_fcon
   end block
   call self%compute_metrix
   endsubroutine flip_edge

   pure subroutine mirror_by_normal(self, normal, recompute_metrix)
   !< Mirror facet given normal of mirroring plane.
   class(facet_object), intent(inout)        :: self             !< Facet.
   type(vector_R8P),    intent(in)           :: normal           !< Normal of mirroring plane.
   logical,             intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.

   call self%mirror_by_matrix(matrix=mirror_matrix_R8P(normal=normal))
   if (present(recompute_metrix)) then
      if (recompute_metrix) call self%compute_metrix
   endif
   endsubroutine mirror_by_normal

   pure subroutine mirror_by_matrix(self, matrix, recompute_metrix)
   !< Mirror facet given matrix (of mirroring).
   class(facet_object), intent(inout)        :: self             !< Facet.
   real(R8P),           intent(in)           :: matrix(3,3)      !< Mirroring matrix.
   logical,             intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.

   call self%vertex(1)%mirror(matrix=matrix)
   call self%vertex(2)%mirror(matrix=matrix)
   call self%vertex(3)%mirror(matrix=matrix)
   if (present(recompute_metrix)) then
      if (recompute_metrix) call self%compute_metrix
   endif
   endsubroutine mirror_by_matrix

   pure subroutine rotate_by_axis_angle(self, axis, angle, center, recompute_metrix)
   !< Rotate facet given axis and angle.
   !<
   !< Angle must be in radians. When `center` is supplied, the rotation pivots
   !< about that point (`v -> R*(v - center) + center`); otherwise it pivots
   !< about the world origin.
   class(facet_object), intent(inout)        :: self             !< Facet.
   type(vector_R8P),    intent(in)           :: axis             !< Axis of rotation.
   real(R8P),           intent(in)           :: angle            !< Angle of rotation.
   type(vector_R8P),    intent(in), optional :: center           !< Rotation centre (default: world origin).
   logical,             intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.

   call self%rotate_by_matrix(matrix=rotation_matrix_R8P(axis=axis, angle=angle), center=center)
   if (present(recompute_metrix)) then
      if (recompute_metrix) call self%compute_metrix
   endif
   endsubroutine rotate_by_axis_angle

   pure subroutine rotate_by_matrix(self, matrix, center, recompute_metrix)
   !< Rotate facet given matrix (of rotation).
   !<
   !< When `center` is supplied, the rotation pivots about that point
   !< (`v -> matrix*(v - center) + center`); otherwise it pivots about the
   !< world origin (`v -> matrix*v`).
   class(facet_object), intent(inout)        :: self             !< Facet.
   real(R8P),           intent(in)           :: matrix(3,3)      !< Rotation matrix.
   type(vector_R8P),    intent(in), optional :: center           !< Rotation centre (default: world origin).
   logical,             intent(in), optional :: recompute_metrix !< Sentinel to activate metrix recomputation.
   integer(I4P)                              :: v                !< Local vertex counter.

   if (present(center)) then
      do v = 1, 3
         self%vertex(v) = self%vertex(v) - center
         call self%vertex(v)%rotate(matrix=matrix)
         self%vertex(v) = self%vertex(v) + center
      enddo
   else
      call self%vertex(1)%rotate(matrix=matrix)
      call self%vertex(2)%rotate(matrix=matrix)
      call self%vertex(3)%rotate(matrix=matrix)
   endif
   if (present(recompute_metrix)) then
      if (recompute_metrix) call self%compute_metrix
   endif
   endsubroutine rotate_by_matrix

   pure subroutine edge_connection_in_other_ref(self, other, edge, edge_vector)
   !< Return the edge of connection in the other reference.
   !<
   !< Searches `other%fcon_edge(:)` for the edge that points back to `self%id` and
   !< returns its index (1..3) and the corresponding edge vector in `other`'s frame.
   class(facet_object), intent(in)  :: self        !< Facet.
   type(facet_object),  intent(in)  :: other       !< Other facet.
   integer(I4P),        intent(out) :: edge        !< Edge index in other numeration (1..3).
   type(vector_R8P),    intent(out) :: edge_vector !< Edge vector in other numeration.
   integer(I4P)                     :: e           !< Counter.

   do e=1, 3
      if (other%fcon_edge(e) == self%id) then
         edge = e
         edge_vector = other%vertex(mod(e, 3) + 1) - other%vertex(e)
         return
      endif
   enddo
   ! self is not connected to other along any edge — caller invariant broken
   error stop 'fossil_facet_object%edge_connection_in_other_ref: facets are not connected'
   endsubroutine edge_connection_in_other_ref

endmodule fossil_facet_object
