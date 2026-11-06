!< FOSSIL test: named patches in the IB queries (issue #66, step 6).
!<
!< A marching-cubes sphere (12^3) split into `top` (facet centroid z > 0) and `bottom`; 400 fixed-seed points in its box.
!< 1. query_many(patches=) is the patch of the closest facet (facets=), and points with |z| > 0.35 get their side's patch.
!< 2. query_many_device(patches=) equals query_many's (here the host build of the device path; the GPU bench has the
!<    device line `patch device vs serial mismatches`).
!< 3. probe_many(patches=) and image_points(patches=) equal query_many's.
program fossil_test_patches_queries

use fossil,                only : surface_stl_object, facet_object
use fossil_boolean_bodies, only : BODY_SPHERE, mc_body
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

integer(I4P), parameter      :: N = 400
type(surface_stl_object)     :: s
type(vector_R8P)             :: pts(N)
real(R8P)                    :: d(N), dd(N), probes(3, N), images(3, N), lengths(N)
integer(I4P)                 :: fac(N), pat(N), pat2(N), pat3(N), pat4(N), i, itop
logical                      :: vis(N), ok1, ok2, ok3
integer(I8P)                 :: seed
type(facet_object), pointer  :: f

call mc_body(BODY_SPHERE, 12_I4P, s)
call split(s)
itop = 0
do i = 1, s%get_patches_number()
   if (s%patch_name(i) == 'top') itop = i
enddo
seed = 20261106_I8P
do i = 1, N
   pts(i) = vector_R8P(rnd() * 1.4_R8P - 0.7_R8P, rnd() * 1.4_R8P - 0.7_R8P, rnd() * 1.4_R8P - 0.7_R8P)
enddo

call s%query_many(points=pts, distances=d, facets=fac, patches=pat)
ok1 = .true.
do i = 1, N
   f => s%facet_at(fac(i))
   ok1 = ok1 .and. pat(i) == f%patch
   if (pts(i)%z > 0.35_R8P) ok1 = ok1 .and. pat(i) == itop
   if (pts(i)%z < -0.35_R8P) ok1 = ok1 .and. pat(i) /= itop
enddo
print '(A,L1)', '1. query_many patches: ', ok1

call s%enter_device
call s%query_many_device(points=pts, distances=dd, patches=pat2)
call s%exit_device
ok2 = all(pat2 == pat)
print '(A,L1)', '2. query_many_device patches equal: ', ok2

lengths = 0.05_R8P
call s%probe_many(points=pts, lengths=lengths, probe_points=probes, visible=vis, patches=pat3)
call s%image_points(points=pts, images=images, patches=pat4)
ok3 = all(pat3 == pat) .and. all(pat4 == pat)
print '(A,L1)', '3. probe_many and image_points patches equal: ', ok3

print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3

contains
   subroutine split(x)
   type(surface_stl_object), intent(inout) :: x
   type(facet_object), pointer             :: g
   integer(I4P), allocatable               :: t(:), b(:)
   integer(I4P)                            :: k

   allocate(t(0), b(0))
   do k = 1, x%get_facets_number()
      g => x%facet_at(k)
      if (g%centroid%z > 0._R8P) then
         t = [t, k]
      else
         b = [b, k]
      endif
   enddo
   call x%set_patch(t, 'top')
   call x%set_patch(b, 'bottom')
   endsubroutine split

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_patches_queries
