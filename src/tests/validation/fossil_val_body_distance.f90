!< FOSSIL validation: body-body distance against brute force (issue #67, step 1).
!<
!< The bunny (decimated to 20 000 facets) against the dragon (decimated to 5 000) under 20 fixed-seed rigid motions that
!< keep them apart, touching or crossing: distance_to (best-first descent of the two BVHs) must give the distance, the
!< facets and the witness points of the brute-force loop over all facet pairs, bit for bit. The speedup is reported.
program fossil_val_body_distance

use fossil,                only : surface_stl_object, facet_object
use fossil_body_distance,  only : body_distance_brute
use fossil_boolean_bodies, only : move_body
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

type(surface_stl_object)        :: a, b, c
type(facet_object), allocatable :: xa(:), xb(:)
type(vector_R8P)                :: ca, cb
real(R8P)                       :: d, d2, pa(3), pb(3), qa(3), qb(3), t0, t1, tt, tb, axis(3), ang, sh(3)
integer(I4P)                    :: fa, fb, ia, ib, k, nbad, nhit
logical                         :: hit, meet
integer(I8P)                    :: seed

call a%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call a%decimate(target_facets=20000_I4P)
call c%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call c%decimate(target_facets=5000_I4P)
call scale_to_unit(a) ; call scale_to_unit(c)
call facets_of(a, xa)
seed = 20261106_I8P ; nbad = 0 ; nhit = 0 ; tt = 0._R8P ; tb = 0._R8P
do k = 1, 20
   b = c
   axis = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] ; ang = 6.28_R8P * rnd()
   sh = (0.3_R8P + 0.9_R8P * rnd()) * [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P]
   call move_body(b, axis, ang, sh)
   call facets_of(b, xb)
   call cpu_time(t0)
   call a%distance_to(b, d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb, intersecting=hit)
   call cpu_time(t1) ; tt = tt + (t1 - t0)
   call body_distance_brute(xa, xb, d2, ia, ib, qa, qb, meet)
   call cpu_time(t0) ; tb = tb + (t0 - t1)
   if (d /= sqrt(d2) .or. fa /= ia .or. fb /= ib .or. any(pa /= qa) .or. any(pb /= qb) .or. (hit .neqv. meet)) &
      nbad = nbad + 1
   if (hit) nhit = nhit + 1
   print '(A,I3,A,ES12.4,A,L1)', 'motion ', k, ': distance ', d, ', intersecting ', hit
enddo
print '(A,I0,A,I0,A,F8.3,A,F8.3,A,F8.1)', 'mismatches ', nbad, ', intersecting ', nhit, ', time tree ', tt, &
      ' s, brute ', tb, ' s, speedup ', tb / max(tt, 1e-9_R8P)
print '(A,L1)', 'Are all tests passed? ', nbad == 0

contains
   subroutine scale_to_unit(s)
   type(surface_stl_object), intent(inout) :: s
   type(vector_R8P)                        :: lo, hi
   real(R8P)                               :: f

   lo = s%get_bmin() ; hi = s%get_bmax()
   call s%translate(delta=-0.5_R8P * (lo + hi))
   f = 1._R8P / norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call s%resize(factor=vector_R8P(f, f, f)) ; call s%analyze
   endsubroutine scale_to_unit

   subroutine facets_of(s, f)
   type(surface_stl_object),        intent(in)  :: s
   type(facet_object), allocatable, intent(out) :: f(:)
   type(facet_object), pointer                  :: p
   integer(I4P)                                 :: i

   allocate(f(s%get_facets_number()))
   do i = 1, size(f)
      p => s%facet_at(i) ; f(i) = p
   enddo
   endsubroutine facets_of

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_val_body_distance
