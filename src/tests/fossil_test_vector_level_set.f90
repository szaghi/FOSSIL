!< FOSSIL test: vector level set, closest point / id / normal stored in the field (issue #37 step 5).
program fossil_test_vector_level_set
!< FOSSIL test: vector level set, closest point / id / normal stored in the field (issue #37 step 5).
!<
!< 3-D (faceted cone `make_apex_body(256, 15)`: clean, closed, curved; 20^3 grid) and 2-D (256-gon of
!< radius 0.5, 64^2 grid), `compute(store_closest=.true.)`:
!< - on every band cell: phi, closest point, id and normal equal query_many at that cell, bit for
!<   bit; |phi| = |x - closest| to 1e-12 relative; the sign of phi is the side of x along the normal;
!< - every non-band cell holds the sentinel (closest = FAR_VALUE, id 0, normal 0);
!< - cube-02 at 16^3: for band cells whose closest point is inside a facet, the image point
!<   closest + |phi| n equals the mirror 2 closest - x to 1e-12 (no grad phi involved);
!< - without store_closest the arrays are not allocated.

use fossil, only : surface_stl_object, curve_object, distance_field_object, FAR_VALUE, FEATURE_FACE
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
use fossil_sign_bench_bodies, only : body_t, make_apex_body, to_surface

implicit none

real(R8P), parameter        :: PI = acos(-1._R8P)
type(surface_stl_object)    :: cone, cube
type(body_t)                :: cb
type(curve_object)          :: circle
type(distance_field_object) :: f
real(R8P)                   :: x(64), u(256), v(256)
integer(I4P)                :: k, st
logical                     :: ok(4)

call make_apex_body(256, 15._R8P, 'cone 15 (256)', cb)
call to_surface(cb%soup, cone)
call grid3(cone, 20_I4P)
ok(1) = check3(cone, 'cone 15 (256)')

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call grid3(cube, 16_I4P)
ok(2) = check3(cube, 'cube-02') .and. image_check(cube)

do k = 1, 256
   u(k) = 0.5_R8P * cos(2._R8P * PI * (k - 1) / 256) ; v(k) = 0.5_R8P * sin(2._R8P * PI * (k - 1) / 256)
enddo
call circle%adopt_loops(u=u, v=v, loop_sizes=[256_I4P])
do k = 1, 64
   x(k) = -1._R8P + (2._R8P * k - 1._R8P) / 64._R8P
enddo
call f%initialize(x=x, y=x)
call f%compute(curve=circle, store_closest=.true., status=st)
ok(3) = st == 0 .and. check2(circle, x)

call f%initialize(x=x, y=x)
call f%compute(curve=circle, status=st)
ok(4) = .not. allocated(f%closest) .and. .not. allocated(f%closest_id) .and. .not. allocated(f%closest_normal)
print '(A,L1)', 'no store_closest, no arrays: ', ok(4)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine grid3(s, n)
   !< n^3 cell-centre grid over 1.4 x the surface box; field computed with store_closest.
   type(surface_stl_object), intent(inout) :: s !< Surface.
   integer(I4P),             intent(in)    :: n !< Cells per axis.
   type(vector_R8P)                        :: lo, hi !< Box.
   real(R8P)                               :: gx(n), gy(n), gz(n) !< Coordinates.
   integer(I4P)                            :: q !< Counter.

   lo = s%get_bmin() ; hi = s%get_bmax()
   do q = 1, n
      gx(q) = lo%x + (hi%x - lo%x) * (-0.2_R8P + 1.4_R8P * (q - 0.5_R8P) / n)
      gy(q) = lo%y + (hi%y - lo%y) * (-0.2_R8P + 1.4_R8P * (q - 0.5_R8P) / n)
      gz(q) = lo%z + (hi%z - lo%z) * (-0.2_R8P + 1.4_R8P * (q - 0.5_R8P) / n)
   enddo
   call f%initialize(x=gx, y=gy, z=gz)
   call f%compute(surface=s, store_closest=.true.)
   endsubroutine grid3

   function check3(s, label) result(good)
   !< Band cells vs query_many, |phi| = |x - cp|, sign along the normal, non-band sentinel.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Label.
   logical                              :: good  !< Result.
   real(R8P),        allocatable :: gx(:), gy(:), gz(:) !< Coordinates (from the field's dims).
   type(vector_R8P)              :: p(1) !< Point.
   real(R8P)                     :: d(1), cp(3, 1), nn(3, 1), r !< query_many outputs.
   integer(I4P)                  :: fid(1), ni, nj, nk, i, j, k, bad, nb !< Work.

   call f%get_dims(ni, nj, nk)
   bad = 0 ; nb = 0
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            if (.not. f%in_band(i, j, k)) then
               if (any(f%closest(:, i, j, k) /= FAR_VALUE) .or. f%closest_id(i, j, k) /= 0 .or. &
                   any(f%closest_normal(:, i, j, k) /= 0._R8P)) bad = bad + 1
               cycle
            endif
            nb = nb + 1
            p(1) = cell_point(i, j, k)
            call s%query_many(points=p, distances=d, closest_points=cp, facets=fid, normals=nn)
            if (f%phi(i, j, k) /= d(1) .or. any(f%closest(:, i, j, k) /= cp(:, 1)) .or. f%closest_id(i, j, k) /= fid(1) .or. &
                any(f%closest_normal(:, i, j, k) /= nn(:, 1))) bad = bad + 1
            r = norm2([p(1)%x, p(1)%y, p(1)%z] - cp(:, 1))
            if (abs(abs(f%phi(i, j, k)) - r) > 1.e-12_R8P * max(r, 1._R8P)) bad = bad + 1
            if (r > 1.e-9_R8P .and. ((f%phi(i, j, k) < 0._R8P) .neqv. &
                (dot_product([p(1)%x, p(1)%y, p(1)%z] - cp(:, 1), nn(:, 1)) < 0._R8P))) bad = bad + 1
         enddo
      enddo
   enddo
   good = bad == 0 .and. nb > 0
   print '(A,A,I0,A,I0)', label, ': band cells ', nb, '  mismatches ', bad
   endfunction check3

   function image_check(s) result(good)
   !< Face-region band cells: closest + |phi| n == 2 closest - x (to 1e-12).
   type(surface_stl_object), intent(in) :: s    !< Surface.
   logical                              :: good !< Result.
   type(vector_R8P)                     :: p(1) !< Point.
   real(R8P)                            :: d(1), img(3), mir(3) !< Work.
   integer(I4P)                         :: ft(1), ni, nj, nk, i, j, k, bad, nf !< Work.

   call f%get_dims(ni, nj, nk)
   bad = 0 ; nf = 0
   do k = 1, nk
      do j = 1, nj
         do i = 1, ni
            if (.not. f%in_band(i, j, k)) cycle
            p(1) = cell_point(i, j, k)
            call s%query_many(points=p, distances=d, features=ft)
            if (ft(1) /= FEATURE_FACE .or. f%phi(i, j, k) >= 0._R8P) cycle ! ghost points: inside, face region
            nf = nf + 1
            img = f%closest(:, i, j, k) + abs(f%phi(i, j, k)) * f%closest_normal(:, i, j, k)
            mir = 2._R8P * f%closest(:, i, j, k) - [p(1)%x, p(1)%y, p(1)%z]
            if (any(abs(img - mir) > 1.e-12_R8P)) bad = bad + 1
         enddo
      enddo
   enddo
   good = bad == 0 .and. nf > 0
   print '(A,I0,A,I0)', 'cube-02 image points (inside, face region): ', nf, '  mismatches ', bad
   endfunction image_check

   function check2(cur, xc) result(good)
   !< 2-D band cells vs query_many, |phi| = |x - cp|, closest(3) = 0, sentinel elsewhere.
   type(curve_object), intent(in) :: cur   !< Curve.
   real(R8P),          intent(in) :: xc(:) !< Coordinates (both axes).
   logical                        :: good  !< Result.
   real(R8P)                      :: d(1), cp(2, 1), nn(2, 1), r !< Work.
   integer(I4P)                   :: sg(1), i, j, bad, nb !< Work.

   bad = 0 ; nb = 0
   do j = 1, size(xc)
      do i = 1, size(xc)
         if (.not. f%in_band(i, j, 1)) then
            if (any(f%closest(:, i, j, 1) /= FAR_VALUE) .or. f%closest_id(i, j, 1) /= 0) bad = bad + 1
            cycle
         endif
         nb = nb + 1
         call cur%query_many(points=reshape([xc(i), xc(j)], [2, 1]), distances=d, closest_points=cp, segments=sg, normals=nn)
         if (f%phi(i, j, 1) /= d(1) .or. any(f%closest(1:2, i, j, 1) /= cp(:, 1)) .or. f%closest(3, i, j, 1) /= 0._R8P .or. &
             f%closest_id(i, j, 1) /= sg(1) .or. any(f%closest_normal(1:2, i, j, 1) /= nn(:, 1))) bad = bad + 1
         r = hypot(xc(i) - cp(1, 1), xc(j) - cp(2, 1))
         if (abs(abs(f%phi(i, j, 1)) - r) > 1.e-12_R8P) bad = bad + 1
      enddo
   enddo
   good = bad == 0 .and. nb > 0
   print '(A,I0,A,I0)', '256-gon (2-D): band cells ', nb, '  mismatches ', bad
   endfunction check2

   function cell_point(i, j, k) result(pt)
   !< Coordinates of cell (i, j, k) of the current field, rebuilt from its grid (3-D).
   integer(I4P), intent(in) :: i, j, k !< Cell.
   type(vector_R8P)         :: pt      !< Point.
   real(R8P)                :: c(3)    !< Coordinates.

   call f%get_coordinates(i, j, k, c)
   pt%x = c(1) ; pt%y = c(2) ; pt%z = c(3)
   endfunction cell_point
endprogram fossil_test_vector_level_set
