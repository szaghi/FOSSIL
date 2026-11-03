!< FOSSIL test: point-triangle distance accurate down to 1e-12 of the size (issue #39).
program fossil_test_distance_tiny
!< FOSSIL test: point-triangle distance accurate down to 1e-12 of the size (issue #39).
!<
!< The squared distance used to be the expanded quadratic form of the Eberly solve, whose absolute
!< error is ~eps L^2: below d ~ 1e-8 L it was noise, could be negative, and then flipped
!< `is_point_inside(SIGN_PSEUDO_NORMAL)`. On the procedural cube and the 5-degree wedge of
!< `fossil_sign_bench_bodies`, points at delta = 1e-6 ... 1e-12 L (side known by construction):
!< - the pseudo-normal inside test equals the constructed side on every point;
!< - the squared distance is never negative, and it is never zero (no point lies on the surface);
!< - on face points the distance equals delta to 1e-14 L.

use fossil,                   only : surface_stl_object, SIGN_PSEUDO_NORMAL
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface, make_points, SB_FACE
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable :: b(:)
type(surface_stl_object)  :: s
real(R8P),    allocatable :: p(:,:), ang(:)
logical,      allocatable :: ins(:)
integer(I4P), allocatable :: knd(:)
integer,      allocatable :: dec(:), fid(:)
integer                   :: i, m, k, nsign, nneg, nface, npts
real(R8P)                 :: d2, d
logical                   :: ok
integer, parameter        :: BODY(2) = [1, 2] ! cube, wedge 5

nsign = 0 ; nneg = 0 ; nface = 0 ; npts = 0
call bench_bodies(b)
do k = 1, size(BODY)
   i = BODY(k)
   call to_surface(b(i)%soup, s)
   call make_points(b(i), [6, 7, 8, 9, 10, 11, 12], p, ins, knd, dec, ang, fid)
   npts = npts + size(ins)
   do m = 1, size(ins)
      if (s%is_point_inside(vector_R8P(p(1, m), p(2, m), p(3, m)), sign_algorithm=SIGN_PSEUDO_NORMAL) .neqv. ins(m)) &
         nsign = nsign + 1
      call s%compute_distance(point=vector_R8P(p(1, m), p(2, m), p(3, m)), distance=d2)
      if (d2 <= 0._R8P) nneg = nneg + 1
      if (knd(m) == SB_FACE) then
         call s%compute_distance(point=vector_R8P(p(1, m), p(2, m), p(3, m)), distance=d, is_square_root=.true.)
         if (abs(d - 10._R8P**(-dec(m)) * b(i)%diag) > 1.e-14_R8P * b(i)%diag) nface = nface + 1
      endif
   enddo
enddo
ok = nsign == 0 .and. nneg == 0 .and. nface == 0
print '(A,I0,A,I0,A,I0,A,I0)', 'points ', npts, ': pseudo-normal sign errors ', nsign, '  non-positive d^2 ', nneg, &
      '  face distance errors ', nface
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_distance_tiny
