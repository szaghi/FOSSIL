!< FOSSIL test: distance between two bodies (issue #67, step 1).
!<
!< 1. Two unit cubes apart by 1 along x, and along the diagonal by d (each axis): distance 1 and sqrt(3) d.
!< 2. Touching at a vertex, along an edge, on a face: 0, intersecting. Overlapping: 0, intersecting.
!< 3. A 10^3 sphere and a 10^3 torus at 6 fixed-seed rigid motions: distance, facets and witness points bit-identical to
!<    the brute-force loop (body_distance_brute); the witness points are at that distance.
!< 4. Two parallel faces 1e-12 apart (a cube and its copy moved by 1 + 1e-12): distance 1e-12 to round-off.
program fossil_test_body_distance

use fossil,                only : surface_stl_object, facet_object
use fossil_body_distance,  only : body_distance_brute
use fossil_boolean_bodies, only : BODY_SPHERE, BODY_TORUS, mc_body, move_body
use penf,                  only : I4P, I8P, R8P

implicit none

real(R8P), parameter            :: OFFS(3, 4) = reshape([1._R8P, 1._R8P, 1._R8P, 1._R8P, 1._R8P, 0._R8P, &
                                                         1._R8P, 0.2_R8P, 0.3_R8P, 0.5_R8P, 0.5_R8P, 0.5_R8P], [3, 4])
type(surface_stl_object)        :: a, b
real(R8P)                       :: d, pa(3), pb(3), d2, qa(3), qb(3), axis(3), ang, t(3)
integer(I4P)                    :: fa, fb, ia, ib, k, nbad
logical                         :: hit, meet, ok1, ok2, ok3, ok4
integer(I8P)                    :: seed
type(facet_object), allocatable :: xa(:), xb(:)

call cube(a, [0._R8P, 0._R8P, 0._R8P]) ; call cube(b, [2._R8P, 0._R8P, 0._R8P])
call a%distance_to(b, d, intersecting=hit)
ok1 = abs(d - 1._R8P) < 1e-15_R8P .and. .not. hit
call cube(b, [1.5_R8P, 1.5_R8P, 1.5_R8P])
call a%distance_to(b, d)
ok1 = ok1 .and. abs(d - sqrt(3._R8P) * 0.5_R8P) < 1e-15_R8P
print '(A,L1)', '1. apart: ', ok1

ok2 = .true.
do k = 1, 4
   call cube(b, OFFS(:, k))
   call a%distance_to(b, d, intersecting=hit)
   ok2 = ok2 .and. d == 0._R8P .and. hit
enddo
print '(A,L1)', '2. touching and overlapping: ', ok2

call mc_body(BODY_SPHERE, 10_I4P, a)
seed = 20261106_I8P ; nbad = 0
do k = 1, 6
   call mc_body(BODY_TORUS, 10_I4P, b)
   axis = [rnd() - 0.5_R8P, rnd() - 0.5_R8P, rnd() - 0.5_R8P] ; ang = 6.28_R8P * rnd()
   t = [rnd() * 3._R8P - 1.5_R8P, rnd() * 3._R8P - 1.5_R8P, rnd() * 3._R8P - 1.5_R8P]
   call move_body(b, axis, ang, t)
   call a%distance_to(b, d, point_self=pa, point_other=pb, facet_self=fa, facet_other=fb)
   call facets_of(a, xa) ; call facets_of(b, xb)
   call body_distance_brute(xa, xb, d2, ia, ib, qa, qb, meet)
   if (d /= sqrt(d2) .or. fa /= ia .or. fb /= ib .or. any(pa /= qa) .or. any(pb /= qb)) nbad = nbad + 1
   if (abs(norm2(pa - pb) - d) > 1e-14_R8P) nbad = nbad + 1
enddo
ok3 = nbad == 0
print '(A,L1,A,I0)', '3. sphere x torus vs brute force: ', ok3, '  mismatches ', nbad

call cube(a, [0._R8P, 0._R8P, 0._R8P]) ; call cube(b, [1._R8P + 1.e-12_R8P, 0._R8P, 0._R8P])
call a%distance_to(b, d, intersecting=hit)
ok4 = abs(d - 1.e-12_R8P) < 1e-15_R8P .and. .not. hit
print '(A,L1,A,ES12.4)', '4. faces 1e-12 apart: ', ok4, '  distance ', d

print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4

contains
   subroutine cube(s, o)
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: o(3)
   type(surface_stl_object)                :: fresh

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call s%translate(x=o(1), y=o(2), z=o(3)) ; call s%analyze
   endsubroutine cube

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
endprogram fossil_test_body_distance
