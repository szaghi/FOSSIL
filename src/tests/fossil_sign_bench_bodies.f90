!< FOSSIL test support: bodies with sharp features and query points of known side (issue #39).
module fossil_sign_bench_bodies
!< FOSSIL test support: bodies with sharp features and query points of known side (issue #39).
!<
!< Procedural, watertight, outward-oriented triangle meshes (a "soup" of vertex triples) with a
!< parameter for the sharpness, their analytic volumes, and the sharp features of each body:
!<
!< - prisms of a CCW polygon (cube, wedges, thin plate, L-shape): side faces, flat diagonals of the
!<   side quads, vertical edges (convex or reentrant), cap edges, prism vertices (convex or saddle);
!< - pyramid and faceted cone: the apex (a convex vertex where 4 or N faces meet), a lateral edge.
!<
!< Query points (decision 1 of #39) are x + delta d, with x a feature point and d a direction in its
!< inside or its outside sector (cone), so the side is known by construction:
!<
!< - face, flat edge: d = -n or +n;
!< - edge between faces 1 and 2: d(phi) = cos(phi) t1 - sin(phi) n1 (t1 in face 1, normal to the
!<   edge, pointing into it); the solid is phi in (0, theta), theta the dihedral angle measured
!<   through the solid; 5 directions at theta k / 6 inside and theta + (360 - theta) k / 6 outside;
!< - prism vertex (an edge ending on a cap normal to it): the edge directions tilted by 30 degrees
!<   away from the cap (inside and outside), plus the cap normal (outside);
!< - apex: inside -sum(n_i) and its averages with the edges leaving the apex; outside the face
!<   normals and +sum(n_i); kept only with a margin min|d . n_i| >= 1e-2 (the apex is the
!<   intersection of its face half-spaces).
!<
!< delta never exceeds a quarter of the distance from x to the nearest facet that does not contain
!< it (`lfs`), so the segment from x to the point crosses no facet.
!<
!< Dirty variants (decision 4): `soup_remove` (a hole), `soup_flip`, `soup_duplicate`,
!< `soup_split_edge` (needle triangles along an edge, kept watertight), `soup_tjunction`.
use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private
public :: soup_t, feature_t, body_t
public :: SB_FACE, SB_FLAT_EDGE, SB_EDGE_CONVEX, SB_EDGE_REENTRANT, SB_VERTEX_CONVEX, SB_VERTEX_SADDLE, SB_KIND_NAME
public :: bench_bodies, make_prism, make_apex_body, to_surface, make_points
public :: soup_remove, soup_flip, soup_duplicate, soup_split_edge, soup_tjunction

integer(I4P), parameter :: SB_FACE           = 1_I4P !< Face interior.
integer(I4P), parameter :: SB_FLAT_EDGE      = 2_I4P !< Diagonal between two coplanar triangles.
integer(I4P), parameter :: SB_EDGE_CONVEX    = 3_I4P !< Convex edge.
integer(I4P), parameter :: SB_EDGE_REENTRANT = 4_I4P !< Reentrant (concave) edge.
integer(I4P), parameter :: SB_VERTEX_CONVEX  = 5_I4P !< Convex vertex.
integer(I4P), parameter :: SB_VERTEX_SADDLE  = 6_I4P !< Saddle vertex (on a reentrant edge).
character(16), parameter :: SB_KIND_NAME(6) = ['face            ', 'flat edge       ', 'convex edge     ', &
                                               'reentrant edge  ', 'convex vertex   ', 'saddle vertex   '] !< Names.

real(R8P), parameter :: PI = acos(-1._R8P) !< Pi.

type :: soup_t
   !< Triangle soup: v(:, k, f) is vertex k of facet f (outward by the right-hand rule).
   real(R8P), allocatable :: v(:,:,:) !< Vertices [3, 3, facets].
endtype soup_t

type :: feature_t
   !< A feature point with the directions of its inside and outside sectors.
   real(R8P)              :: x(3) = 0._R8P   !< Feature point.
   integer(I4P)           :: kind = 0_I4P    !< SB_*.
   real(R8P)              :: angle = 0._R8P  !< Opening angle of the solid at the feature [deg] (180 on faces).
   real(R8P)              :: lfs = 0._R8P    !< Largest admissible delta (a quarter of the clearance).
   real(R8P), allocatable :: din(:,:)        !< Inside directions [3, m].
   real(R8P), allocatable :: dout(:,:)       !< Outside directions [3, m].
endtype feature_t

type :: body_t
   !< A benchmark body.
   character(32)                :: name = ''       !< Name.
   type(soup_t)                 :: soup            !< Clean mesh.
   type(feature_t), allocatable :: feat(:)         !< Features.
   real(R8P)                    :: volume = 0._R8P !< Analytic volume.
   real(R8P)                    :: diag = 0._R8P   !< Bounding-box diagonal L.
endtype body_t

contains
   subroutine bench_bodies(b)
   !< The benchmark set: cube, wedges 5/30/150 deg, thin plate (1e-3), pyramid (10 deg), faceted cone
   !< (15 deg, 32 faces), L-prism.
   type(body_t), allocatable, intent(out) :: b(:) !< Bodies.
   real(R8P) :: t !< Work.

   allocate(b(9))
   call make_prism(reshape([0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P], [2, 4]), 1._R8P, 1, 'cube', b(1))
   call make_prism(wedge(5._R8P), 1._R8P, 1, 'wedge 5', b(2))
   call make_prism(wedge(30._R8P), 1._R8P, 1, 'wedge 30', b(3))
   call make_prism(wedge(150._R8P), 1._R8P, 1, 'wedge 150', b(4))
   t = 1.e-3_R8P
   call make_prism(reshape([0._R8P, 0._R8P, 1._R8P, 0._R8P, 1._R8P, 1._R8P, 0._R8P, 1._R8P], [2, 4]), t, 1, 'plate 1e-3', b(5))
   call make_apex_body(4, 10._R8P, 'pyramid 10', b(6))
   call make_apex_body(32, 15._R8P, 'cone 15 (32)', b(7))
   call make_prism(reshape([1._R8P, 1._R8P, 1._R8P, 2._R8P, 0._R8P, 2._R8P, 0._R8P, 0._R8P, 2._R8P, 0._R8P, 2._R8P, 1._R8P], &
                           [2, 6]), 1._R8P, 1, 'L-prism', b(8))
   call make_prism(wedge(90._R8P), 1._R8P, 1, 'wedge 90', b(9))
   contains
      pure function wedge(theta) result(p)
      !< Triangle with apex angle theta at the origin, sides of length 1, symmetric about +x.
      real(R8P), intent(in) :: theta !< Apex angle [deg].
      real(R8P)             :: p(2, 3) !< Vertices (CCW).
      real(R8P)             :: h       !< Half angle [rad].

      h = 0.5_R8P * theta * PI / 180._R8P
      p(:, 1) = [0._R8P, 0._R8P] ; p(:, 2) = [cos(h), -sin(h)] ; p(:, 3) = [cos(h), sin(h)]
      endfunction wedge
   endsubroutine bench_bodies

   subroutine make_prism(p, h, fan, name, b)
   !< Prism of the CCW polygon p(2, n) over z in [0, h]; caps fanned from vertex `fan` (the polygon
   !< must be star-shaped from it). Features: the centroid of the first triangle of side 1 (face),
   !< the centre of side 1 (flat diagonal), the mid-height of every vertical edge, the mid-point of
   !< the top edge of side 1, the top end of every vertical edge, the cap centroid of the first fan
   !< triangle.
   real(R8P),    intent(in)  :: p(:,:) !< Polygon [2, n], CCW.
   real(R8P),    intent(in)  :: h      !< Height.
   integer,      intent(in)  :: fan    !< Fan vertex of the caps.
   character(*), intent(in)  :: name   !< Name.
   type(body_t), intent(out) :: b      !< Body.
   integer         :: n, k, k1, f, nf, m !< Work.
   real(R8P)       :: a(3), bb(3), c(3), d(3), area, nk(3), nkm(3), e(3), q1(3), q2(3), x(3) !< Work.
   type(feature_t) :: ft(64) !< Features buffer.

   n = size(p, 2)
   nf = 2 * n + 2 * (n - 2)
   allocate(b%soup%v(3, 3, nf))
   f = 0
   do k = 1, n ! side quads (a, b, c, d) = (Pk, Pk+1, Pk+1 + h, Pk + h), diagonal a-c
      k1 = mod(k, n) + 1
      a = [p(:, k), 0._R8P] ; bb = [p(:, k1), 0._R8P] ; c = [p(:, k1), h] ; d = [p(:, k), h]
      f = f + 1 ; b%soup%v(:, :, f) = reshape([a, bb, c], [3, 3])
      f = f + 1 ; b%soup%v(:, :, f) = reshape([a, c, d], [3, 3])
   enddo
   do k = 1, n - 2 ! caps, fan from vertex `fan`
      a  = [p(:, fan), 0._R8P]
      bb = [p(:, mod(fan + k - 1, n) + 1), 0._R8P]
      c  = [p(:, mod(fan + k, n) + 1), 0._R8P]
      f = f + 1 ; b%soup%v(:, :, f) = reshape([a, c, bb], [3, 3]) ! bottom, outward -z
      a(3) = h ; bb(3) = h ; c(3) = h
      f = f + 1 ; b%soup%v(:, :, f) = reshape([a, bb, c], [3, 3]) ! top, outward +z
   enddo
   area = 0._R8P
   do k = 1, n
      k1 = mod(k, n) + 1
      area = area + 0.5_R8P * (p(1, k) * p(2, k1) - p(1, k1) * p(2, k))
   enddo
   b%name = name ; b%volume = area * h
   b%diag = norm2([maxval(p(1, :)) - minval(p(1, :)), maxval(p(2, :)) - minval(p(2, :)), h])
   m = 0
   ! side 1: face and flat diagonal
   nk = side_normal(1)
   m = m + 1 ; ft(m) = flat(sum(b%soup%v(:, :, 1), dim=2) / 3._R8P, nk, SB_FACE)
   m = m + 1 ; ft(m) = flat(0.5_R8P * (b%soup%v(:, 1, 1) + b%soup%v(:, 3, 1)), nk, SB_FLAT_EDGE)
   ! top cap face
   m = m + 1 ; ft(m) = flat(sum(b%soup%v(:, :, 2 * n + 2), dim=2) / 3._R8P, [0._R8P, 0._R8P, 1._R8P], SB_FACE)
   ! vertical edges and prism vertices at the top
   do k = 1, n
      nkm = side_normal(mod(k + n - 2, n) + 1) ! side k-1 (ends at Pk)
      nk  = side_normal(k)                     ! side k (starts at Pk)
      x = [p(:, k), 0.5_R8P * h] ; e = [0._R8P, 0._R8P, 1._R8P]
      q1 = [0.5_R8P * (p(:, k) + p(:, mod(k + n - 2, n) + 1)), 0.5_R8P * h] ! inside side k-1
      q2 = [0.5_R8P * (p(:, k) + p(:, mod(k, n) + 1)), 0.5_R8P * h]         ! inside side k
      m = m + 1 ; ft(m) = edge(x, e, nkm, q1, q2)
      m = m + 1 ; ft(m) = prism_vertex(ft(m - 1), [p(:, k), h], [0._R8P, 0._R8P, 1._R8P])
   enddo
   ! top edge of side 1 (side face and top cap)
   x = [0.5_R8P * (p(:, 1) + p(:, 2)), h] ; e = [p(:, 2) - p(:, 1), 0._R8P] ; e = e / norm2(e)
   q1 = [0.5_R8P * (p(:, 1) + p(:, 2)), 0.5_R8P * h]
   q2 = [sum(b%soup%v(1:2, :, 2 * n + 2), dim=2) / 3._R8P, h]
   m = m + 1 ; ft(m) = edge(x, e, side_normal(1), q1, q2)
   b%feat = ft(1:m)
   call set_lfs(b)
   contains
      pure function side_normal(s) result(nn)
      !< Outward normal of side s (CCW polygon: (dy, -dx)).
      integer, intent(in) :: s     !< Side.
      real(R8P)           :: nn(3) !< Normal.
      real(R8P)           :: dd(2) !< Direction.

      dd = p(:, mod(s, n) + 1) - p(:, s)
      nn = [dd(2), -dd(1), 0._R8P] / norm2(dd)
      endfunction side_normal
   endsubroutine make_prism

   subroutine make_apex_body(n, alpha, name, b)
   !< Pyramid (n = 4) or faceted cone: regular n-gon base of circumradius 0.5 at z = 0, apex on the
   !< axis at the height where the faces make the angle alpha with the axis. Features: the apex, the
   !< mid-point of the lateral edge 1, the centroid of lateral face 1.
   integer,      intent(in)  :: n     !< Base vertices.
   real(R8P),    intent(in)  :: alpha !< Half-angle between the axis and a face [deg].
   character(*), intent(in)  :: name  !< Name.
   type(body_t), intent(out) :: b     !< Body.
   real(R8P)              :: r, ri, h, apex(3), q(3, n), nn(3, n), s(3), d(3), g(3), dl(3, 2 * n + 1), w(3) !< Work.
   integer                :: k, k1, f, mi, mo !< Work.
   type(feature_t)        :: fa       !< Apex feature.

   r = 0.5_R8P ; ri = r * cos(PI / n) ; h = ri / tan(alpha * PI / 180._R8P)
   apex = [0._R8P, 0._R8P, h]
   do k = 1, n
      q(:, k) = [r * cos(2 * PI * (k - 1) / n), r * sin(2 * PI * (k - 1) / n), 0._R8P]
   enddo
   allocate(b%soup%v(3, 3, 2 * n - 2))
   f = 0
   do k = 1, n
      k1 = mod(k, n) + 1
      f = f + 1 ; b%soup%v(:, :, f) = reshape([q(:, k), q(:, k1), apex], [3, 3])
      nn(:, k) = cross(q(:, k1) - q(:, k), apex - q(:, k)) ; nn(:, k) = nn(:, k) / norm2(nn(:, k))
   enddo
   do k = 2, n - 1
      f = f + 1 ; b%soup%v(:, :, f) = reshape([q(:, 1), q(:, k + 1), q(:, k)], [3, 3])
   enddo
   b%name = name ; b%volume = 0.5_R8P * n * r * r * sin(2 * PI / n) * h / 3._R8P
   b%diag = norm2([2 * r, 2 * r, h])
   ! apex: half-space intersection of the lateral faces
   fa%x = apex ; fa%kind = SB_VERTEX_CONVEX ; fa%angle = 2._R8P * alpha
   s = sum(nn, dim=2) ; s = s / norm2(s)
   mi = 0 ; mo = 0
   call try_in(-s)
   do k = 1, n
      g = q(:, k) - apex ; g = g / norm2(g)
      d = -s + g ; call try_in(d / norm2(d))
   enddo
   allocate(fa%din(3, mi)) ; fa%din = dl(:, 1:mi)
   mo = 0
   call try_out(s)
   do k = 1, n
      call try_out(nn(:, k))
   enddo
   allocate(fa%dout(3, mo)) ; fa%dout = dl(:, 1:mo)
   allocate(b%feat(3))
   b%feat(1) = fa
   ! lateral edge 1 (q1 -> apex) between faces n and 1
   w = apex - q(:, 1) ; w = w / norm2(w)
   b%feat(2) = edge(0.5_R8P * (q(:, 1) + apex), w, nn(:, n), sum(b%soup%v(:, :, n), dim=2) / 3._R8P, &
                    sum(b%soup%v(:, :, 1), dim=2) / 3._R8P)
   b%feat(3) = flat(sum(b%soup%v(:, :, 1), dim=2) / 3._R8P, nn(:, 1), SB_FACE)
   call set_lfs(b)
   contains
      subroutine try_in(dd)
      !< Keep dd as an inside direction if every face is at least 1e-2 away.
      real(R8P), intent(in) :: dd(3) !< Direction.

      if (maxval(matmul(dd, nn)) <= -1.e-2_R8P) then
         mi = mi + 1 ; dl(:, mi) = dd
      endif
      endsubroutine try_in

      subroutine try_out(dd)
      !< Keep dd as an outside direction if some face half-space is violated by at least 1e-2.
      real(R8P), intent(in) :: dd(3) !< Direction.

      if (maxval(matmul(dd, nn)) >= 1.e-2_R8P) then
         mo = mo + 1 ; dl(:, mo) = dd
      endif
      endsubroutine try_out
   endsubroutine make_apex_body

   function flat(x, n, kind) result(ft)
   !< Face or flat-edge feature: inside -n, outside +n.
   real(R8P),    intent(in) :: x(3), n(3) !< Point, outward normal.
   integer(I4P), intent(in) :: kind       !< SB_FACE or SB_FLAT_EDGE.
   type(feature_t)          :: ft         !< Feature.

   ft%x = x ; ft%kind = kind ; ft%angle = 180._R8P
   allocate(ft%din(3, 1), ft%dout(3, 1))
   ft%din(:, 1) = -n ; ft%dout(:, 1) = n
   endfunction flat

   function edge(x, e, n1, q1, q2) result(ft)
   !< Edge feature at x (unit direction e) between face 1 (normal n1, interior point q1) and face 2.
   real(R8P), intent(in) :: x(3), e(3)   !< Point, edge direction.
   real(R8P), intent(in) :: n1(3), q1(3) !< Face 1: outward normal, interior point.
   real(R8P), intent(in) :: q2(3)        !< Face 2: interior point.
   type(feature_t)       :: ft           !< Feature.
   real(R8P)             :: t1(3), t2(3), theta, phi !< Work.
   integer               :: k            !< Counter.

   t1 = q1 - x ; t1 = t1 - dot_product(t1, e) * e ; t1 = t1 / norm2(t1)
   t2 = q2 - x ; t2 = t2 - dot_product(t2, e) * e ; t2 = t2 / norm2(t2)
   theta = atan2(-dot_product(t2, n1), dot_product(t2, t1)) ! d(phi) = cos(phi) t1 - sin(phi) n1
   if (theta <= 0._R8P) theta = theta + 2._R8P * PI
   ft%x = x ; ft%angle = theta * 180._R8P / PI
   ft%kind = SB_EDGE_CONVEX ; if (theta > PI) ft%kind = SB_EDGE_REENTRANT
   allocate(ft%din(3, 5), ft%dout(3, 5))
   do k = 1, 5
      phi = theta * k / 6._R8P
      ft%din(:, k) = cos(phi) * t1 - sin(phi) * n1
      phi = theta + (2._R8P * PI - theta) * k / 6._R8P
      ft%dout(:, k) = cos(phi) * t1 - sin(phi) * n1
   enddo
   endfunction edge

   function prism_vertex(fe, x, nc) result(ft)
   !< Vertex at the end x of an edge feature fe, on a cap with outward normal nc (normal to the edge):
   !< the edge directions tilted by 30 degrees away from the cap, plus nc (outside).
   type(feature_t), intent(in) :: fe     !< Edge feature.
   real(R8P),       intent(in) :: x(3)   !< Vertex.
   real(R8P),       intent(in) :: nc(3)  !< Cap outward normal.
   type(feature_t)             :: ft     !< Feature.
   real(R8P), parameter        :: TILT = 30._R8P * PI / 180._R8P !< Tilt.
   integer                     :: k      !< Counter.

   ft%x = x ; ft%angle = fe%angle
   ft%kind = SB_VERTEX_CONVEX ; if (fe%kind == SB_EDGE_REENTRANT) ft%kind = SB_VERTEX_SADDLE
   allocate(ft%din(3, size(fe%din, 2)), ft%dout(3, size(fe%dout, 2) + 1))
   do k = 1, size(fe%din, 2)
      ft%din(:, k) = cos(TILT) * fe%din(:, k) - sin(TILT) * nc
   enddo
   do k = 1, size(fe%dout, 2)
      ft%dout(:, k) = cos(TILT) * fe%dout(:, k) - sin(TILT) * nc
   enddo
   ft%dout(:, size(fe%dout, 2) + 1) = nc
   endfunction prism_vertex

   subroutine set_lfs(b)
   !< lfs of each feature: a quarter of the distance from its point to the nearest facet that does
   !< not contain it (distance > 1e-13 L).
   type(body_t), intent(inout) :: b !< Body.
   integer   :: i, f !< Counters.
   real(R8P) :: d2, dmin !< Work.

   do i = 1, size(b%feat)
      dmin = huge(1._R8P)
      do f = 1, size(b%soup%v, 3)
         d2 = point_triangle_d2(b%feat(i)%x, b%soup%v(:, 1, f), b%soup%v(:, 2, f), b%soup%v(:, 3, f))
         if (sqrt(d2) > 1.e-13_R8P * b%diag) dmin = min(dmin, sqrt(d2))
      enddo
      b%feat(i)%lfs = 0.25_R8P * dmin
   enddo
   endsubroutine set_lfs

   subroutine to_surface(soup, s)
   !< Build a surface from a soup (adopt_facets: metrics, connectivity, tree, pseudo-normals).
   type(soup_t),             intent(in)    :: soup !< Soup.
   type(surface_stl_object), intent(inout) :: s    !< Surface.
   type(facet_object), allocatable :: f(:) !< Facets.
   integer                         :: i, k !< Counters.

   allocate(f(size(soup%v, 3)))
   do i = 1, size(f)
      do k = 1, 3
         f(i)%vertex(k) = vector_R8P(soup%v(1, k, i), soup%v(2, k, i), soup%v(3, k, i))
      enddo
   enddo
   call s%adopt_facets(facets=f)
   endsubroutine to_surface

   subroutine make_points(b, decades, p, inside, kind, decade, angle, feature)
   !< Query points of every feature: x + delta d over the decades delta = 10^-k L (k in `decades`)
   !< that do not exceed the feature's lfs, every inside and outside direction.
   type(body_t),              intent(in)  :: b          !< Body.
   integer,                   intent(in)  :: decades(:) !< Exponents k.
   real(R8P),    allocatable, intent(out) :: p(:,:)     !< Points [3, m].
   logical,      allocatable, intent(out) :: inside(:)  !< Side by construction [m].
   integer(I4P), allocatable, intent(out) :: kind(:)    !< Feature kind [m].
   integer,      allocatable, intent(out) :: decade(:)  !< Decade k [m].
   real(R8P),    allocatable, intent(out) :: angle(:)   !< Feature opening angle [m].
   integer,      allocatable, intent(out) :: feature(:) !< Feature index [m].
   integer   :: i, j, k, m, n !< Counters.
   real(R8P) :: delta !< Distance.

   n = 0
   do i = 1, size(b%feat)
      do k = 1, size(decades)
         if (10._R8P**(-decades(k)) * b%diag <= b%feat(i)%lfs) n = n + size(b%feat(i)%din, 2) + size(b%feat(i)%dout, 2)
      enddo
   enddo
   allocate(p(3, n), inside(n), kind(n), decade(n), angle(n), feature(n))
   m = 0
   do i = 1, size(b%feat)
      do k = 1, size(decades)
         delta = 10._R8P**(-decades(k)) * b%diag
         if (delta > b%feat(i)%lfs) cycle
         do j = 1, size(b%feat(i)%din, 2)
            m = m + 1 ; p(:, m) = b%feat(i)%x + delta * b%feat(i)%din(:, j) ; inside(m) = .true.
            kind(m) = b%feat(i)%kind ; decade(m) = decades(k) ; angle(m) = b%feat(i)%angle ; feature(m) = i
         enddo
         do j = 1, size(b%feat(i)%dout, 2)
            m = m + 1 ; p(:, m) = b%feat(i)%x + delta * b%feat(i)%dout(:, j) ; inside(m) = .false.
            kind(m) = b%feat(i)%kind ; decade(m) = decades(k) ; angle(m) = b%feat(i)%angle ; feature(m) = i
         enddo
      enddo
   enddo
   endsubroutine make_points

   subroutine soup_remove(s, f)
   !< Remove facet f (a hole).
   type(soup_t), intent(inout) :: s !< Soup.
   integer,      intent(in)    :: f !< Facet.
   integer                     :: j !< Counter.

   s%v = s%v(:, :, pack([(j, j = 1, size(s%v, 3))], [(j /= f, j = 1, size(s%v, 3))]))
   endsubroutine soup_remove

   subroutine soup_flip(s, f)
   !< Reverse the orientation of facet f.
   type(soup_t), intent(inout) :: s !< Soup.
   integer,      intent(in)    :: f !< Facet.
   real(R8P)                   :: t(3) !< Swap.

   t = s%v(:, 2, f) ; s%v(:, 2, f) = s%v(:, 3, f) ; s%v(:, 3, f) = t
   endsubroutine soup_flip

   subroutine soup_duplicate(s, f)
   !< Append a copy of facet f.
   type(soup_t), intent(inout) :: s !< Soup.
   integer,      intent(in)    :: f !< Facet.
   real(R8P), allocatable      :: w(:,:,:) !< Work.

   allocate(w(3, 3, size(s%v, 3) + 1))
   w(:, :, 1:size(s%v, 3)) = s%v ; w(:, :, size(w, 3)) = s%v(:, :, f)
   call move_alloc(w, s%v)
   endsubroutine soup_duplicate

   subroutine soup_split_edge(s, a, b, t)
   !< Split the edge (a, b) at a + t(j) (b - a) in every facet that has it (needles, watertight).
   type(soup_t), intent(inout) :: s      !< Soup.
   real(R8P),    intent(in)    :: a(3), b(3) !< Edge end points.
   real(R8P),    intent(in)    :: t(:)   !< Increasing parameters in (0, 1).
   real(R8P), allocatable      :: w(:,:,:) !< New soup.
   real(R8P)                   :: q(3, size(t) + 2), c(3) !< Points along the edge, opposite vertex.
   integer                     :: f, k, m, i1, i2 !< Work.
   logical                     :: fwd  !< Edge found as a -> b in the facet orientation.

   allocate(w(3, 3, size(s%v, 3) * (size(t) + 1)))
   m = 0
   do f = 1, size(s%v, 3)
      call find_edge(s%v(:, :, f), i1, i2, fwd)
      if (i1 == 0) then
         m = m + 1 ; w(:, :, m) = s%v(:, :, f)
         cycle
      endif
      c = s%v(:, 6 - i1 - i2, f)
      q(:, 1) = a ; q(:, size(t) + 2) = b
      do k = 1, size(t)
         q(:, k + 1) = a + t(k) * (b - a)
      enddo
      do k = 1, size(t) + 1
         m = m + 1
         if (fwd) then
            w(:, :, m) = reshape([q(:, k), q(:, k + 1), c], [3, 3])
         else
            w(:, :, m) = reshape([q(:, k + 1), q(:, k), c], [3, 3])
         endif
      enddo
   enddo
   s%v = w(:, :, 1:m)
   contains
      subroutine find_edge(v, j1, j2, forward)
      !< Local indices of a and b in facet v (0 if absent) and whether a -> b follows the orientation.
      real(R8P), intent(in)  :: v(3, 3)  !< Facet.
      integer,   intent(out) :: j1, j2   !< Local indices of a and b.
      logical,   intent(out) :: forward  !< a -> b in the facet orientation.
      integer                :: l        !< Counter.

      j1 = 0 ; j2 = 0
      do l = 1, 3
         if (all(v(:, l) == a)) j1 = l
         if (all(v(:, l) == b)) j2 = l
      enddo
      if (j1 == 0 .or. j2 == 0) then
         j1 = 0 ; j2 = 0 ; forward = .false.
         return
      endif
      forward = mod(j1, 3) + 1 == j2
      endsubroutine find_edge
   endsubroutine soup_split_edge

   subroutine soup_tjunction(s, f, t)
   !< Split edge 1-2 of facet f alone at parameter t (a T-junction: the neighbour is not split).
   type(soup_t), intent(inout) :: s !< Soup.
   integer,      intent(in)    :: f !< Facet.
   real(R8P),    intent(in)    :: t !< Parameter in (0, 1).
   real(R8P)                   :: q(3), v(3, 3) !< Split point, facet.
   real(R8P), allocatable      :: w(:,:,:) !< Work.

   v = s%v(:, :, f)
   q = v(:, 1) + t * (v(:, 2) - v(:, 1))
   allocate(w(3, 3, size(s%v, 3) + 1))
   w(:, :, 1:size(s%v, 3)) = s%v
   w(:, :, f) = reshape([v(:, 1), q, v(:, 3)], [3, 3])
   w(:, :, size(w, 3)) = reshape([q, v(:, 2), v(:, 3)], [3, 3])
   call move_alloc(w, s%v)
   endsubroutine soup_tjunction

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   pure function point_triangle_d2(p, a, b, c) result(d2)
   !< Squared distance from p to the triangle abc (Ericson, Real-Time Collision Detection 5.1.5).
   real(R8P), intent(in) :: p(3), a(3), b(3), c(3) !< Point, triangle.
   real(R8P)             :: d2                     !< Squared distance.
   real(R8P)             :: ab(3), ac(3), ap(3), bp(3), cp(3), d1, d2_, d3, d4, d5, d6, va, vb, vc, v, w, q(3) !< Work.

   ab = b - a ; ac = c - a ; ap = p - a
   d1 = dot_product(ab, ap) ; d2_ = dot_product(ac, ap)
   if (d1 <= 0._R8P .and. d2_ <= 0._R8P) then
      q = a
   else
      bp = p - b ; d3 = dot_product(ab, bp) ; d4 = dot_product(ac, bp)
      cp = p - c ; d5 = dot_product(ab, cp) ; d6 = dot_product(ac, cp)
      vc = d1 * d4 - d3 * d2_ ; vb = d5 * d2_ - d1 * d6 ; va = d3 * d6 - d5 * d4
      if (d3 >= 0._R8P .and. d4 <= d3) then
         q = b
      elseif (vc <= 0._R8P .and. d1 >= 0._R8P .and. d3 <= 0._R8P) then
         q = a + d1 / (d1 - d3) * ab
      elseif (d6 >= 0._R8P .and. d5 <= d6) then
         q = c
      elseif (vb <= 0._R8P .and. d2_ >= 0._R8P .and. d6 <= 0._R8P) then
         q = a + d2_ / (d2_ - d6) * ac
      elseif (va <= 0._R8P .and. (d4 - d3) >= 0._R8P .and. (d5 - d6) >= 0._R8P) then
         q = b + (d4 - d3) / ((d4 - d3) + (d5 - d6)) * (c - b)
      else
         v = vb / (va + vb + vc) ; w = vc / (va + vb + vc)
         q = a + ab * v + ac * w
      endif
   endif
   d2 = sum((p - q)**2)
   endfunction point_triangle_d2
endmodule fossil_sign_bench_bodies
