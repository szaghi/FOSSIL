!< FOSSIL test: the load samples of a surface (#109 WP2 step 01).
program fossil_test_load_samples_3d
!< FOSSIL test: the load samples of a surface (#109 WP2 step 01).
!<
!< 1. Weights: on an icosphere (320 facets) the weights of every facet sum to its area, all of them to the surface area,
!<    to round-off, for the three orders; q = 1, 3, 6 points per facet; normals outward (n . p > 0), points on the facet
!<    plane.
!< 2. Exactness: on one triangle the rule of order 1, 2, 3 integrates every monomial x^a y^b z^c of degree 1, 2, 4 to
!<    round-off (the exact integral from the barycentric moments), and misses one of the next degree.
!< 3. Visibility: two unit boxes 0.25 apart; probes 0.1 out all visible, 0.5 out blocked exactly on the two facing faces,
!<    1.5 in all blocked (the opposite face); a zero length is the point itself, visible.
!< 4. Cells: every probe's cell on stretched lines against a brute-force search (low-side rule, 0 outside).
!< 5. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host in every field, bit for bit
!<    (`load_samples_mismatches`), with lines, posed too.
!< 6. Posed: points and probes = `to_world` of the unposed ones, normals rotated (to round-off), weights and visibility
!<    bit for bit; against a copy moved by `move` to round-off, the same visibility.
!< 7. Resident lifecycle, refusals (order, partial or bad lines, unmapped surface, probes on an octree), the empty surface,
!<    a degenerate facet (zero weight and normal, probes not visible).
use fossil,                  only : surface_stl_object, load_samples_t, load_samples_3d, load_samples_reset, &
                                    load_samples_mismatches, rigid_pose_t, pose_from_axis_angle, to_world, rotate, &
                                    STATUS_OK, STATUS_INVALID_INPUT, AABB_TREE_OCTREE
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_boxes, make_sphere
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: EPS = epsilon(1._R8P)
type(surface_stl_object)        :: sphere, boxes, tri, moved, octree, empty, degen
type(facet_object), allocatable :: fc(:), fd(:)
type(load_samples_t)            :: sh, sd, sr, sc, sp
type(rigid_pose_t)              :: pose
real(R8P), allocatable          :: x(:), y(:), z(:), bad(:)
real(R8P)                       :: e, area, a1, v(3, 3)
integer(I4P)                    :: i, f, s, k, o, st, nbad, nblk
integer(I4P), parameter         :: QQ(3) = [1_I4P, 3_I4P, 6_I4P], DEG(3) = [1_I4P, 2_I4P, 4_I4P]
logical                         :: ok, miss

ok = .true.
call make_sphere(2_I4P, fc) ; call sphere%adopt_facets(facets=fc)
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, 0.25_R8P, fc) ; call boxes%adopt_facets(facets=fc)
! 1. weights, normals, points
do o = 1, 3
   call load_samples_3d(sphere, samples=sh, status=st, order=o)
   call check('sphere order '//char(48 + o)//': status, q, n', st == STATUS_OK .and. sh%q == QQ(o) .and. &
              sh%n == 320_I4P * QQ(o) .and. sh%ne == 320_I4P .and. sh%m == 0_I4P .and. .not. sh%with_cells)
   e = 0._R8P ; area = 0._R8P
   do f = 1, 320
      a1 = facet_area(sphere, f) ; area = area + a1
      e = max(e, abs(sum(sh%weight((f - 1) * sh%q + 1:f * sh%q)) - a1) / a1)
   enddo
   call check('sphere order '//char(48 + o)//': weights of a facet sum to its area (8 eps)', e <= 8._R8P * EPS)
   call check('sphere order '//char(48 + o)//': all weights sum to the area', abs(sum(sh%weight) - area) <= 1.e-13_R8P * area)
   nbad = 0
   do s = 1, sh%n
      f = sh%element(s)
      v = facet_vertices(sphere, f)
      if (dot_product(sh%normal(:, s), sh%points(:, s)) <= 0._R8P) nbad = nbad + 1
      if (abs(norm2(sh%normal(:, s)) - 1._R8P) > 4._R8P * EPS) nbad = nbad + 1
      if (abs(dot_product(sh%points(:, s) - v(:, 1), sh%normal(:, s))) > 8._R8P * EPS) nbad = nbad + 1
      if (sh%element(s) /= (s - 1) / sh%q + 1 .or. sh%body(s) /= 1 .or. sh%patch(s) /= 1) nbad = nbad + 1
   enddo
   call check('sphere order '//char(48 + o)//': unit outward normals, points on the facet, ids', nbad == 0)
enddo
! 2. exactness on one triangle
allocate(fc(1))
fc(1)%vertex(1) = vector_R8P(0.3_R8P, -0.2_R8P, 0.5_R8P) ; fc(1)%vertex(2) = vector_R8P(1.1_R8P, 0.4_R8P, 0.2_R8P)
fc(1)%vertex(3) = vector_R8P(0.1_R8P, 0.9_R8P, 1.3_R8P) ; call fc(1)%compute_metrix
call tri%adopt_facets(facets=fc)
v = facet_vertices(tri, 1_I4P)
do o = 1, 3
   call load_samples_3d(tri, samples=sh, status=st, order=o)
   e = 0._R8P ; miss = .false.
   do i = 0, DEG(o) + 1
      do k = 0, DEG(o) + 1 - i
         do f = 0, DEG(o) + 1 - i - k
            if (i + k + f > DEG(o) + 1) cycle
            a1 = 0._R8P
            do s = 1, sh%n
               a1 = a1 + sh%weight(s) * sh%points(1, s)**i * sh%points(2, s)**k * sh%points(3, s)**f
            enddo
            if (i + k + f <= DEG(o)) then
               e = max(e, abs(a1 - tri_integral(v, i, k, f)))
            elseif (abs(a1 - tri_integral(v, i, k, f)) > 1.e-8_R8P) then
               miss = .true.
            endif
         enddo
      enddo
   enddo
   call check('triangle order '//char(48 + o)//': exact to degree '//char(48 + DEG(o))//' (1e-14), not to the next', &
              e < 1.e-14_R8P .and. miss)
enddo
! 3. visibility blocked by a second body
call load_samples_3d(boxes, samples=sh, status=st, lengths=[0.1_R8P, 0.5_R8P, -1.5_R8P, 0._R8P])
nbad = 0 ; nblk = 0
do s = 1, sh%n
   if (.not. sh%visible(1, s)) nbad = nbad + 1
   if (sh%visible(3, s)) nbad = nbad + 1
   if (.not. sh%visible(4, s) .or. any(sh%probe(:, 4, s) /= sh%points(:, s))) nbad = nbad + 1
   if (.not. sh%visible(2, s)) then
      nblk = nblk + 1
      if (.not. ((sh%body(s) == 1 .and. sh%normal(1, s) == 1._R8P) .or. (sh%body(s) == 2 .and. sh%normal(1, s) == -1._R8P))) &
         nbad = nbad + 1
   endif
   if (any(sh%probe(:, 2, s) /= sh%points(:, s) + 0.5_R8P * sh%normal(:, s))) nbad = nbad + 1
enddo
call check('boxes: 0.1 out visible, 1.5 in blocked, length 0 the point', nbad == 0)
call check('boxes: 0.5 out blocked on the facing faces only (4 facets x 3 points)', nblk == 12)
! 4. cells against brute force
x = [(-0.3_R8P + 0.13_R8P * i + 0.004_R8P * i * i, i = 0, 22)] ; y = [(-0.6_R8P + 0.2_R8P * i, i = 0, 9)]
z = [(-0.25_R8P + 0.25_R8P * i, i = 0, 6)]                   ! node lines through the faces z = 0 and 1
call load_samples_3d(boxes, x, y, z, sh, status=st, order=3_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P])
call check('cells: lines, flags', st == STATUS_OK .and. sh%with_cells .and. sh%ni == 23 .and. sh%nk == 7)
call check('cells: brute force on every probe', count_bad_cells(sh) == 0)
! 5. twins
call sphere%enter_device ; call boxes%enter_device
call load_samples_3d(boxes, x, y, z, sd, status=st, order=3_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P], device=.true.)
call check('twin: device = host, every field', st == STATUS_OK .and. load_samples_mismatches(sh, sd) == 0)
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.05_R8P, -0.02_R8P, 0.03_R8P])
call load_samples_3d(boxes, x, y, z, sh, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P], pose=pose)
call load_samples_3d(boxes, x, y, z, sd, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P], pose=pose, &
                     device=.true.)
call check('twin posed: device = host, every field', st == STATUS_OK .and. load_samples_mismatches(sh, sd) == 0)
call load_samples_3d(sphere, x, y, z, sh, status=st, order=1_I4P, lengths=[0.2_R8P])
call load_samples_3d(sphere, x, y, z, sd, status=st, order=1_I4P, lengths=[0.2_R8P], device=.true.)
call check('twin sphere: device = host', load_samples_mismatches(sh, sd) == 0)
! 6. posed against unposed and against a moved copy
call load_samples_3d(boxes, x, y, z, sh, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P], pose=pose)
call load_samples_3d(boxes, samples=sp, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P])
e = 0._R8P
do s = 1, sh%n
   e = max(e, norm2(sh%points(:, s) - to_world(pose, sp%points(:, s))), norm2(sh%normal(:, s) - rotate(pose, sp%normal(:, s))))
   do k = 1, 3
      e = max(e, norm2(sh%probe(:, k, s) - to_world(pose, sp%probe(:, k, s))))
   enddo
enddo
call check('posed: points, probes = to_world, normals rotated (1e-14)', e < 1.e-14_R8P)
call check('posed: weights and visibility = unposed, bit for bit', all(sh%weight == sp%weight) .and. &
           all(sh%visible .eqv. sp%visible) .and. count(.not. sh%visible(2, :)) == 2 * 2 * 3)
call check('posed: cells brute force', count_bad_cells(sh) == 0)
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 2_I4P, 0.25_R8P, fc) ; call moved%adopt_facets(facets=fc)
call moved%move(pose)
call load_samples_3d(moved, samples=sc, status=st, order=2_I4P, lengths=[0.1_R8P, 0.5_R8P, -0.4_R8P])
call check('posed: = a moved copy (1e-14), the same visibility', maxval(abs(sc%points - sh%points)) < 1.e-14_R8P .and. &
           maxval(abs(sc%probe - sh%probe)) < 1.e-14_R8P .and. maxval(abs(sc%weight - sh%weight)) < 1.e-14_R8P .and. &
           all(sc%visible .eqv. sh%visible))
! 7. resident lifecycle
call load_samples_3d(boxes, x, y, z, sh, status=st, order=3_I4P, lengths=[0.1_R8P, 0.5_R8P])
call load_samples_3d(boxes, x, y, z, sr, status=st, order=3_I4P, lengths=[0.1_R8P, 0.5_R8P], resident=.true.)
call check('resident: flagged', st == STATUS_OK .and. sr%on_device)
call sr%update_host
call check('resident: update_host gives the host bits', load_samples_mismatches(sh, sr) == 0 .and. sr%on_device)
sc = sr
call sc%update_host                                     ! a copy: warns on stderr, drops the flag
call check('resident: a copy is not mapped', .not. sc%on_device .and. load_samples_mismatches(sh, sc) == 0)
call sr%exit_device
call check('resident: exit_device', .not. sr%on_device)
call sr%enter_device
call sr%enter_device
call check('resident: enter_device, twice a no-op', sr%on_device)
call load_samples_3d(boxes, samples=sr, status=st, device=.true.)
call check('resident: a new call releases it first', st == STATUS_OK .and. .not. sr%on_device .and. sr%n == 72)
call load_samples_3d(boxes, samples=sr, status=st, resident=.true.)
call load_samples_reset(sr)
call check('resident: load_samples_reset', .not. sr%on_device .and. sr%n == 0 .and. .not. allocated(sr%points))
! refusals, the empty surface, a degenerate facet
call load_samples_3d(boxes, samples=sd, status=st, order=4_I4P)
call check('refused: order 4', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_3d(boxes, samples=sd, status=st, order=0_I4P)
call check('refused: order 0', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_3d(boxes, x, y, samples=sd, status=st, lengths=[0.1_R8P])
call check('refused: lines in part', st == STATUS_INVALID_INPUT .and. sd%n == 0)
bad = [0._R8P, 1._R8P, 1._R8P]
call load_samples_3d(boxes, x, bad, z, sd, status=st, lengths=[0.1_R8P])
call check('refused: a line not increasing', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_3d(boxes, x, [0._R8P], z, sd, status=st, lengths=[0.1_R8P])
call check('refused: a line of one point', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call boxes%exit_device
call load_samples_3d(boxes, samples=sd, status=st, device=.true.)
call check('refused: device path of an unmapped surface (stderr message)', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 1_I4P, 0._R8P, fc)
call octree%adopt_facets(facets=fc, aabb_tree_kind=AABB_TREE_OCTREE)
call load_samples_3d(octree, samples=sd, status=st, lengths=[0.1_R8P])
call check('refused: probes on an octree (stderr message)', st == STATUS_INVALID_INPUT .and. sd%n == 0)
call load_samples_3d(octree, samples=sd, status=st, order=3_I4P)
call check('octree: samples without probes', st == STATUS_OK .and. sd%n == 72)
call load_samples_3d(empty, x, y, z, sd, status=st, lengths=[0.1_R8P], device=.true.)
call check('empty surface: no sample, both paths', st == STATUS_OK .and. sd%n == 0 .and. allocated(sd%points))
call make_boxes([0._R8P, 0._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], 1_I4P, 0._R8P, fc)
allocate(fd(13))
do i = 1, 12
   fd(i) = fc(i)
enddo
fd(13)%vertex(1) = vector_R8P(0.2_R8P, 0.2_R8P, 0.5_R8P) ; fd(13)%vertex(2) = vector_R8P(0.4_R8P, 0.4_R8P, 0.5_R8P)
fd(13)%vertex(3) = vector_R8P(0.6_R8P, 0.6_R8P, 0.5_R8P) ; call fd(13)%compute_metrix
call degen%adopt_facets(facets=fd)
call load_samples_3d(degen, samples=sd, status=st, order=2_I4P, lengths=[0.1_R8P])
call check('degenerate facet: zero weights and normals, probes at the points, not visible', st == STATUS_OK .and. &
           all(sd%weight(37:39) == 0._R8P) .and. all(sd%normal(:, 37:39) == 0._R8P) .and. .not. any(sd%visible(1, 37:39)) .and. &
           all(sd%probe(:, 1, 37:39) == sd%points(:, 37:39)) .and. all(sd%visible(1, 1:36)))
call sphere%exit_device

print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(label, cond)
   !< Report one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   if (.not. cond) then
      print '(A)', 'FAILED: '//label
      ok = .false.
   else
      print '(A)', 'ok: '//label
   endif
   endsubroutine check

   function facet_vertices(surf, f) result(vv)
   !< The three vertices of facet f [3, 3].
   type(surface_stl_object), intent(in) :: surf    !< Surface.
   integer(I4P),             intent(in) :: f       !< Facet.
   real(R8P)                            :: vv(3, 3) !< Vertices.

   vv = reshape(surf%tri_vertices(:, f), [3, 3])
   endfunction facet_vertices

   real(R8P) function facet_area(surf, f) result(a)
   !< Area of facet f.
   type(surface_stl_object), intent(in) :: surf !< Surface.
   integer(I4P),             intent(in) :: f    !< Facet.
   real(R8P) :: vv(3, 3), c(3)                  !< Work.

   vv = facet_vertices(surf, f)
   c = cross(vv(:, 2) - vv(:, 1), vv(:, 3) - vv(:, 1))
   a = 0.5_R8P * norm2(c)
   endfunction facet_area

   pure function cross(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

   real(R8P) function tri_integral(vv, ia, ib, ic) result(r)
   !< Exact integral of x^ia y^ib z^ic over the triangle vv: the monomial as a polynomial in the barycentric coordinates,
   !< integrated term by term (int l1^i l2^j l3^k = 2 A i! j! k! / (i + j + k + 2)!).
   real(R8P),    intent(in) :: vv(3, 3)   !< Vertices.
   integer(I4P), intent(in) :: ia, ib, ic !< Exponents.
   real(R8P)    :: p(0:6, 0:6, 0:6), pn(0:6, 0:6, 0:6) !< Polynomials.
   integer(I4P) :: d, j, a, b, c                       !< Counters.

   p = 0._R8P ; p(0, 0, 0) = 1._R8P
   do d = 1, 3
      do j = 1, merge(ia, merge(ib, ic, d == 2), d == 1)
         pn = 0._R8P
         do c = 0, 5 ; do b = 0, 5 ; do a = 0, 5
            pn(a + 1, b, c) = pn(a + 1, b, c) + vv(d, 1) * p(a, b, c)
            pn(a, b + 1, c) = pn(a, b + 1, c) + vv(d, 2) * p(a, b, c)
            pn(a, b, c + 1) = pn(a, b, c + 1) + vv(d, 3) * p(a, b, c)
         enddo ; enddo ; enddo
         p = pn
      enddo
   enddo
   r = 0._R8P
   do c = 0, 6 ; do b = 0, 6 ; do a = 0, 6
      if (p(a, b, c) /= 0._R8P) r = r + p(a, b, c) * 2._R8P * fact(a) * fact(b) * fact(c) / fact(a + b + c + 2)
   enddo ; enddo ; enddo
   r = r * 0.5_R8P * norm2(cross(vv(:, 2) - vv(:, 1), vv(:, 3) - vv(:, 1)))
   endfunction tri_integral

   pure real(R8P) function fact(n)
   !< n!.
   integer(I4P), intent(in) :: n !< Argument.
   integer(I4P) :: i             !< Counter.

   fact = 1._R8P
   do i = 2, n
      fact = fact * i
   enddo
   endfunction fact

   integer(I4P) function count_bad_cells(smp) result(nb)
   !< Probes whose cells differ from a brute-force search on x, y, z (low-side rule, 0 outside).
   type(load_samples_t), intent(in) :: smp !< Samples.
   integer(I4P) :: s, k, a                 !< Counters.

   nb = 0
   do s = 1, smp%n
      do k = 1, smp%m
         do a = 1, 3
            select case (a)
            case (1) ; if (smp%probe_cell(a, k, s) /= brute(x, smp%probe(a, k, s))) nb = nb + 1
            case (2) ; if (smp%probe_cell(a, k, s) /= brute(y, smp%probe(a, k, s))) nb = nb + 1
            case default ; if (smp%probe_cell(a, k, s) /= brute(z, smp%probe(a, k, s))) nb = nb + 1
            endselect
         enddo
      enddo
   enddo
   endfunction count_bad_cells

   pure integer(I4P) function brute(e, c) result(idx)
   !< The cell i with e(i) < c <= e(i+1) by a linear scan, 0 if none.
   real(R8P), intent(in) :: e(:) !< Lines.
   real(R8P), intent(in) :: c    !< Coordinate.
   integer(I4P) :: i             !< Counter.

   idx = 0
   do i = 1, size(e) - 1
      if (e(i) < c .and. c <= e(i + 1)) idx = i
   enddo
   endfunction brute
endprogram fossil_test_load_samples_3d
