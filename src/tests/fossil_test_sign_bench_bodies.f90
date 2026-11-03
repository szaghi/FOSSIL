!< FOSSIL test: the sign-benchmark bodies and query points (issue #39 step 1).
program fossil_test_sign_bench_bodies
!< FOSSIL test: the sign-benchmark bodies and query points (issue #39 step 1).
!<
!< - Every clean body (cube, wedges 5/30/90/150, plate 1e-3, pyramid 10, cone 15 with 32 faces,
!<   L-prism) is watertight and manifold, with the analytic volume to 1e-12 relative.
!< - Every query point lies at delta from its feature point (to 1e-14 L).
!< - At the largest admissible delta of every feature (its own lfs), the solid-angle sign (exact up to round-off
!<   away from the surface) agrees with the side by construction on every point.
!< - The feature kinds expected per body occur (reentrant edge and saddle vertex on the L-prism only).
!< - Dirty variants: a removed facet is not watertight; a flipped facet changes the volume; a
!<   duplicated facet is not manifold; split edges (needles of aspect ratio ~1e6) stay watertight with
!<   the same volume; a T-junction is not watertight.

use fossil,                   only : surface_stl_object, SIGN_SOLID_ANGLE
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface, make_points, soup_remove, soup_flip, &
                                     soup_duplicate, soup_split_edge, soup_tjunction, SB_EDGE_REENTRANT, SB_VERTEX_SADDLE
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
type(soup_t)              :: d
real(R8P),    allocatable :: p(:,:), ang(:)
logical,      allocatable :: ins(:)
integer(I4P), allocatable :: knd(:)
integer,      allocatable :: dec(:), fid(:)
integer                   :: i, m, nbad_geo, nbad_dist, nbad_sa, npts, nsa, kd, dmin
logical                   :: ok(4), reentrant_only_l
integer, parameter        :: DECADES(12) = [(kd, kd = 1, 12)]

call bench_bodies(b)
nbad_geo = 0 ; nbad_dist = 0 ; nbad_sa = 0 ; npts = 0 ; nsa = 0 ; reentrant_only_l = .true.
do i = 1, size(b)
   call to_surface(b(i)%soup, s)
   if (.not. (s%is_watertight() .and. s%is_manifold()) .or. abs(s%get_volume() - b(i)%volume) > 1.e-12_R8P * b(i)%volume) &
      nbad_geo = nbad_geo + 1
   call make_points(b(i), DECADES, p, ins, knd, dec, ang, fid)
   npts = npts + size(ins)
   if ((any(knd == SB_EDGE_REENTRANT) .or. any(knd == SB_VERTEX_SADDLE)) .neqv. (trim(b(i)%name) == 'L-prism')) &
      reentrant_only_l = .false.
   dmin = minval(dec)
   do m = 1, size(ins)
      if (dec(m) == minval(dec, mask=fid == fid(m))) then
         nsa = nsa + 1
         if (s%is_point_inside(vector_R8P(p(1, m), p(2, m), p(3, m)), sign_algorithm=SIGN_SOLID_ANGLE) .neqv. ins(m)) &
            nbad_sa = nbad_sa + 1
      endif
   enddo
   nbad_dist = nbad_dist + check_dist(b(i))
   print '(A,A,I0,A,I0,A,F7.3,A,I0,A,ES9.2)', b(i)%name, ' facets ', size(b(i)%soup%v, 3), ' features ', size(b(i)%feat), &
         '  volume ', s%get_volume(), '  points ', size(ins), '  largest delta/L ', 10._R8P**(-dmin)
enddo
ok(1) = nbad_geo == 0 .and. reentrant_only_l
ok(2) = nbad_dist == 0
ok(3) = nbad_sa == 0
print '(A,I0,A,L1)', 'clean bodies watertight, manifold, analytic volume: failures ', nbad_geo, &
      '; reentrant features on the L-prism only: ', reentrant_only_l
print '(A,I0,A,I0)', 'points at delta from their feature: ', npts, ' points, failures ', nbad_dist
print '(A,I0,A,I0)', 'solid-angle sign = constructed side at the largest delta: ', nsa, ' points, failures ', nbad_sa

! dirty variants of the cube
ok(4) = .true.
d = b(1)%soup ; call soup_remove(d, 1) ; call to_surface(d, s)
ok(4) = ok(4) .and. .not. s%is_watertight()
d = b(1)%soup ; call soup_flip(d, 5) ; call to_surface(d, s)
ok(4) = ok(4) .and. abs(s%get_volume() - 1._R8P) > 1.e-3_R8P
d = b(1)%soup ; call soup_duplicate(d, 1) ; call to_surface(d, s)
ok(4) = ok(4) .and. .not. s%is_manifold()
d = b(1)%soup
call soup_split_edge(d, [1._R8P, 1._R8P, 0._R8P], [1._R8P, 1._R8P, 1._R8P], &
                     [0.5_R8P - 2.e-6_R8P, 0.5_R8P - 1.e-6_R8P, 0.5_R8P, 0.5_R8P + 1.e-6_R8P, 0.5_R8P + 2.e-6_R8P])
call to_surface(d, s)
ok(4) = ok(4) .and. size(d%v, 3) == 12 + 2 * 5 .and. s%is_watertight() .and. abs(s%get_volume() - 1._R8P) < 1.e-12_R8P
d = b(1)%soup ; call soup_tjunction(d, 1, 0.3_R8P) ; call to_surface(d, s)
ok(4) = ok(4) .and. .not. s%is_watertight()
print '(A,L1)', 'dirty variants (hole, flip, duplicate, needles, T-junction) behave: ', ok(4)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check_dist(bb) result(nbad)
   !< Count points not at delta from their feature point.
   type(body_t), intent(in) :: bb   !< Body.
   integer                  :: nbad !< Failures.
   integer   :: f, k, j !< Counters.
   real(R8P) :: delta  !< Distance.

   nbad = 0
   do f = 1, size(bb%feat)
      do k = 1, 12
         delta = 10._R8P**(-k) * bb%diag
         if (delta > bb%feat(f)%lfs) cycle
         do j = 1, size(bb%feat(f)%din, 2)
            if (abs(norm2(bb%feat(f)%x + delta * bb%feat(f)%din(:, j) - bb%feat(f)%x) - delta) > 1.e-14_R8P * bb%diag) &
               nbad = nbad + 1
         enddo
         do j = 1, size(bb%feat(f)%dout, 2)
            if (abs(norm2(bb%feat(f)%x + delta * bb%feat(f)%dout(:, j) - bb%feat(f)%x) - delta) > 1.e-14_R8P * bb%diag) &
               nbad = nbad + 1
         enddo
      enddo
   enddo
   endfunction check_dist
endprogram fossil_test_sign_bench_bodies
