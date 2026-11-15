!< FOSSIL test: grid-adapted surface markers (issue #42 step 4).
program fossil_test_surface_markers
!< FOSSIL test: grid-adapted surface markers (issue #42 step 4).
!<
!< - A unit icosphere (2 subdivisions), dx = 0.06: the marker areas sum to the input area within 2 %,
!<   the normals are unit and outward, and the nearest-marker spacings lie in [0.8, 1.2] dx for more
!<   than 85 % of the markers (0.912 on the first run).
!< - A unit cube (faces split in 2 x 2) with sharp features kept, dx = 0.1: area 6 to round-off, normals
!<   on the six axis directions, spacing fraction at least 0.65 (0.723 on the first run; 0.528 after issues #79 and
!<   #91, when the bound was 0.48; 0.726 since issue #92 WP5 c). The spacing is the distance in space between the
!<   centroids of facets sharing an edge: across a 90-degree crease two equilateral facets of edge sqrt(3) dx on a
!<   crease edge are 0.707 dx apart, below the window, and they are a share 3 dx of the facets of a unit cube, so an
!<   ideal marker set (equilateral, rows along the creases) gets 1 - 3 dx = 0.70. The spacings are printed in bins of
!<   0.1 dx, for the markers within 0.75 dx of a crease and for the others.
!< - The crease vertices in the edge-length window (issue #92 WP5 c): `isotropic_remesh(preserve_features=.true.)` of
!<   that cube at L = sqrt(3) 0.1 splits every edge above 4L/3 each iteration, collapses the straight crease runs along
!<   their line and relaxes their vertices along it: every cube edge keeps 7 or 8 vertices (nint(1 / L) = 6 segments,
!<   or 7) and every crease spacing lies in [4L/5, 4L/3] (0.807 L to 1.019 L measured; 11 edges at 6 segments, one at
!<   7). Before, 9 vertices on every edge at 0.125 = 0.722 L: the crease vertices stayed where the splits put them.
!< - The input surface is unchanged (intent(in)); dx <= 0 is refused (negative control).
!< - Curves: a unit square and a 64-gon resampled exactly: every marker on the polyline (distance
!<   below 1e-15), lengths all equal to L / n with n = nint(L / dx), summing to the perimeter; an open
!<   line gets n = nint(L / dx) markers; the chord spacing fraction of the 64-gon is 1.

use fossil, only : surface_stl_object, facet_object, curve_object, STATUS_INVALID_INPUT, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P)
type(surface_stl_object)        :: s
type(curve_object)              :: c
type(facet_object), allocatable :: f(:)
real(R8P), allocatable          :: cen(:,:), area(:), nrm(:,:)
real(R8P)                       :: frac, a0, u(64), v(64), d, l0
integer(I4P)                    :: i, st, nf0, m
real(R8P), allocatable          :: sp(:)
real(R8P)                       :: de, ce(3), smin, smax
integer(I4P)                    :: hc(13), hf(13), k, j, ncv(12)
logical                         :: ok(6)

! sphere
call icosphere(2_I4P, f) ; a0 = 0._R8P
do i = 1, size(f)
   a0 = a0 + f(i)%area()
enddo
call s%adopt_facets(facets=f) ; nf0 = s%get_facets_number()
call s%surface_markers(dx=0.06_R8P, centroids=cen, areas=area, normals=nrm, spacing_fraction=frac, status=st)
m = size(area)
ok(1) = st == STATUS_OK .and. abs(sum(area) - a0) < 0.02_R8P * a0 .and. s%get_facets_number() == nf0 .and. &
        all(abs(norm2(nrm, dim=1) - 1._R8P) < 1.e-12_R8P) .and. all(sum(nrm * cen, dim=1) > 0.9_R8P * norm2(cen, dim=1)) .and. &
        frac > 0.85_R8P
print '(A,I0,A,F8.5,A,F8.5,A,F6.3,A,L1)', 'sphere dx=0.06: markers ', m, '  area ', sum(area), ' (input ', a0, &
      ')  spacing in [0.8, 1.2] dx: ', frac, '  ok ', ok(1)
! cube
call box(2_I4P, f) ; call s%adopt_facets(facets=f)
call s%surface_markers(dx=0.1_R8P, centroids=cen, areas=area, normals=nrm, spacing_fraction=frac, spacings=sp, &
                       preserve_features=.true., status=st)
ok(2) = st == STATUS_OK .and. abs(sum(area) - 6._R8P) < 1.e-12_R8P .and. all(maxval(abs(nrm), dim=1) > 1._R8P - 1.e-9_R8P) .and. &
        frac >= 0.65_R8P
print '(A,I0,A,F8.5,A,F6.3,A,L1)', 'cube dx=0.1: markers ', size(area), '  area ', sum(area), '  spacing in [0.8, 1.2] dx: ', &
      frac, '  ok ', ok(2)
hc = 0 ; hf = 0
do i = 1, size(area)
   ce = min(cen(:, i), 1._R8P - cen(:, i))
   de = minval(ce, mask=[(j /= minloc(ce, 1), j = 1, 3)]) ! distance to the nearer crease (the smallest is the face's)
   k = min(13, int(sp(i) / 0.01_R8P) + 1)
   if (de < 0.075_R8P) then
      hc(k) = hc(k) + 1
   else
      hf(k) = hf(k) + 1
   endif
enddo
print '(A)', '  spacing / dx in bins of 0.1 (0-0.1, ..., 1.1-1.2, >= 1.2):'
print '(A,13I4)', '  near a crease ', hc
print '(A,13I4)', '  on the faces  ', hf
! the crease vertices in the edge-length window
call box(2_I4P, f) ; call s%adopt_facets(facets=f)
call s%isotropic_remesh(target_length=sqrt(3._R8P) * 0.1_R8P, preserve_features=.true., status=st)
call crease_vertices(s, ncv, smin, smax)
smin = smin / (sqrt(3._R8P) * 0.1_R8P) ; smax = smax / (sqrt(3._R8P) * 0.1_R8P)
ok(6) = st == STATUS_OK .and. all(ncv == 7_I4P .or. ncv == 8_I4P) .and. smin >= 0.8_R8P .and. smax <= 4._R8P / 3._R8P
print '(A,12I3,A,2F7.3,A,L1)', 'cube crease vertices: per cube edge ', ncv, '  spacings / L ', smin, smax, &
      '  ok ', ok(6)
call s%surface_markers(dx=0._R8P, centroids=cen, areas=area, normals=nrm, status=st)
ok(3) = st == STATUS_INVALID_INPUT
call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
call c%curve_markers(dx=-1._R8P, centroids=cen, lengths=area, normals=nrm, status=st)
ok(3) = ok(3) .and. st == STATUS_INVALID_INPUT
print '(A,L1)', 'dx <= 0 refused: ', ok(3)
! curves
call c%curve_markers(dx=0.07_R8P, centroids=cen, lengths=area, normals=nrm, status=st)
d = 0._R8P
do i = 1, size(area)
   d = max(d, abs(c%distance(cen(:, i))))
enddo
ok(4) = st == STATUS_OK .and. size(area) == nint(4._R8P / 0.07_R8P) .and. d < 1.e-15_R8P .and. &
        maxval(abs(area - 4._R8P / size(area))) < 1.e-15_R8P .and. abs(sum(area) - 4._R8P) < 1.e-13_R8P .and. &
        all(abs(norm2(nrm, dim=1) - 1._R8P) < 1.e-15_R8P)
print '(A,I0,A,ES9.1,A,L1)', 'square dx=0.07: markers ', size(area), '  max distance to the polyline ', d, '  ok ', ok(4)
do i = 1, 64
   u(i) = cos(2 * PI * (i - 1) / 64) ; v(i) = sin(2 * PI * (i - 1) / 64)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
l0 = 64 * 2 * sin(PI / 64)
call c%curve_markers(dx=0.05_R8P, centroids=cen, lengths=area, normals=nrm, spacing_fraction=frac, status=st)
d = 0._R8P
do i = 1, size(area)
   d = max(d, abs(c%distance(cen(:, i))))
enddo
ok(5) = st == STATUS_OK .and. size(area) == nint(l0 / 0.05_R8P) .and. d < 1.e-15_R8P .and. abs(sum(area) - l0) < 1.e-13_R8P .and. &
        frac == 1._R8P
call c%adopt_loops(u=[0._R8P, 1._R8P, 3._R8P], v=[0._R8P, 0._R8P, 0._R8P], loop_sizes=[3_I4P], closed=[.false.])
call c%curve_markers(dx=0.4_R8P, centroids=cen, lengths=area, normals=nrm, status=st)
ok(5) = ok(5) .and. size(area) == 8 .and. maxval(abs(area - 0.375_R8P)) < 1.e-15_R8P .and. abs(cen(1, 1) - 0.1875_R8P) < 1.e-15_R8P
print '(A,I0,A,ES9.1,A,F6.3,A,I0,A,L1)', '64-gon dx=0.05: markers ', nint(l0 / 0.05_R8P), '  max distance ', d, &
      '  spacing fraction ', frac, '  open line markers ', size(area), '  ok ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine crease_vertices(su, nc, smin, smax)
   !< Number of distinct vertices on each of the 12 edges of the unit cube (corners included), and the shortest and
   !< longest spacing between consecutive ones.
   type(surface_stl_object), intent(in)  :: su         !< Surface.
   integer(I4P),             intent(out) :: nc(12)     !< Vertices per cube edge.
   real(R8P),                intent(out) :: smin, smax !< Spacings.
   type(facet_object), pointer           :: fq
   real(R8P)                             :: p(3), t(400, 12), x
   integer(I4P)                          :: kk, jj, a, e, i2, j2
   logical                               :: on(3)

   nc = 0
   do kk = 1, su%get_facets_number()
      fq => su%facet_at(kk)
      do jj = 1, 3
         p = [fq%vertex(jj)%x, fq%vertex(jj)%y, fq%vertex(jj)%z]
         on = p == 0._R8P .or. p == 1._R8P
         if (count(on) < 2) cycle
         do a = 1, 3 ! the cube edges along axis a through p: the other two coordinates fixed
            if (.not. (on(mod(a, 3) + 1) .and. on(mod(a + 1, 3) + 1))) cycle
            e = 4 * (a - 1) + 2 * nint(p(mod(a, 3) + 1)) + nint(p(mod(a + 1, 3) + 1)) + 1
            if (any(t(1:nc(e), e) == p(a))) cycle
            nc(e) = nc(e) + 1 ; t(nc(e), e) = p(a)
         enddo
      enddo
   enddo
   smin = huge(1._R8P) ; smax = 0._R8P
   do e = 1, 12
      do i2 = 2, nc(e) ! insertion sort along the edge
         x = t(i2, e) ; j2 = i2 - 1
         do while (j2 >= 1)
            if (t(j2, e) <= x) exit
            t(j2 + 1, e) = t(j2, e) ; j2 = j2 - 1
         enddo
         t(j2 + 1, e) = x
      enddo
      smin = min(smin, minval(t(2:nc(e), e) - t(1:nc(e) - 1, e)))
      smax = max(smax, maxval(t(2:nc(e), e) - t(1:nc(e) - 1, e)))
   enddo
   endsubroutine crease_vertices

   subroutine put(fo, k, a, b, cc)
   !< Facet k from three points.
   type(facet_object), intent(inout) :: fo(:)     !< Facets.
   integer(I4P),       intent(in)    :: k         !< Index.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo(k)%vertex(1) = vector_R8P(a(1), a(2), a(3))
   fo(k)%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo(k)%vertex(3) = vector_R8P(cc(1), cc(2), cc(3))
   call fo(k)%compute_metrix
   endsubroutine put

   subroutine icosphere(levels, fo)
   !< Unit icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere.
   integer(I4P),                    intent(in)  :: levels   !< Subdivisions.
   type(facet_object), allocatable, intent(out) :: fo(:)    !< Facets.
   real(R8P), allocatable  :: tri(:,:,:), nxt(:,:,:)
   real(R8P)               :: vv(3, 12), g
   integer(I4P), parameter :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                    2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                    4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                    5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   integer(I4P)            :: k, l, mm

   g = (1._R8P + sqrt(5._R8P)) / 2._R8P
   vv = reshape([-1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, 0._R8P, &
                 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P, -g, &
                 g, 0._R8P, -1._R8P, g, 0._R8P, 1._R8P, -g, 0._R8P, -1._R8P, -g, 0._R8P, 1._R8P], [3, 12])
   do k = 1, 12
      vv(:, k) = vv(:, k) / norm2(vv(:, k))
   enddo
   allocate(tri(3, 3, 20))
   do k = 1, 20
      tri(:, :, k) = vv(:, IDX(:, k))
   enddo
   do l = 1, levels
      mm = size(tri, dim=3)
      allocate(nxt(3, 3, 4 * mm))
      do k = 1, mm
         nxt(:, :, 4 * k - 3) = reshape([tri(:, 1, k), mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
         nxt(:, :, 4 * k - 2) = reshape([mid(tri(:, 1, k), tri(:, 2, k)), tri(:, 2, k), mid(tri(:, 2, k), tri(:, 3, k))], [3, 3])
         nxt(:, :, 4 * k - 1) = reshape([mid(tri(:, 3, k), tri(:, 1, k)), mid(tri(:, 2, k), tri(:, 3, k)), tri(:, 3, k)], [3, 3])
         nxt(:, :, 4 * k)     = reshape([mid(tri(:, 1, k), tri(:, 2, k)), mid(tri(:, 2, k), tri(:, 3, k)), &
                                         mid(tri(:, 3, k), tri(:, 1, k))], [3, 3])
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   mm = size(tri, dim=3)
   allocate(fo(mm))
   do k = 1, mm
      call put(fo, k, tri(:, 1, k), tri(:, 2, k), tri(:, 3, k))
   enddo
   endsubroutine icosphere

   pure function mid(a, b) result(mm)
   !< Unit midpoint.
   real(R8P), intent(in) :: a(3), b(3) !< Points.
   real(R8P)             :: mm(3)      !< Midpoint on the unit sphere.

   mm = (a + b) / norm2(a + b)
   endfunction mid

   subroutine box(n, fo)
   !< Unit cube, every face an n x n grid of squares split in two triangles (outward normals).
   integer(I4P),                    intent(in)  :: n  !< Squares per side.
   type(facet_object), allocatable, intent(out) :: fo(:) !< Facets.
   real(R8P)    :: a(3), b(3), cc(3), dd(3), s0, s1, t0, t1
   integer(I4P) :: ax, side, i, j, mm, ia, ib

   allocate(fo(12 * n * n))
   mm = 0
   do ax = 1, 3
      ia = mod(ax, 3) + 1 ; ib = mod(ax + 1, 3) + 1
      do side = 0, 1
         do j = 0, n - 1
            do i = 0, n - 1
               s0 = real(i, R8P) / n ; s1 = real(i + 1, R8P) / n ; t0 = real(j, R8P) / n ; t1 = real(j + 1, R8P) / n
               a(ax) = side ; b(ax) = side ; cc(ax) = side ; dd(ax) = side
               a(ia) = s0 ; a(ib) = t0 ; b(ia) = s1 ; b(ib) = t0 ; cc(ia) = s1 ; cc(ib) = t1 ; dd(ia) = s0 ; dd(ib) = t1
               if (side == 1) then
                  mm = mm + 1 ; call put(fo, mm, a, b, cc)
                  mm = mm + 1 ; call put(fo, mm, a, cc, dd)
               else
                  mm = mm + 1 ; call put(fo, mm, a, cc, b)
                  mm = mm + 1 ; call put(fo, mm, a, dd, cc)
               endif
            enddo
         enddo
      enddo
   enddo
   endsubroutine box
endprogram fossil_test_surface_markers
