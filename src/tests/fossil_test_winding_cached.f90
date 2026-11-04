!< FOSSIL, test `surface%winding_number` on the persistent cache (issue #29, step 4).
!<
!< Regression variant (issue #45): dragon.stl (6588 facets) only, no timing. The full version (bunny and dragon-fine,
!< plus the time per fast query on bunny) is src/tests/validation/fossil_val_winding_cached.f90.
!<
!< For the SAH BVH and the octree, on dragon, over 10^3 grid points:
!<   1. fast (default beta = 2) vs exact (beta <= 0): max |dw| < 1e-2 (bound set on bunny; measured on
!<      dragon: BVH 5.7e-4, octree 6.3e-4, BVH vs octree 7.2e-4). The former implementation failed this
!<      on both tree kinds: far-field sign error (|dw| ~ 2 on bunny), and on octrees every facet counted
!<      once per tree level;
!<   2. the default equals an explicit beta = 2 bit for bit;
!<   3. BVH and octree fast results agree within 2e-2 (both within 1e-2 of exact).

program fossil_test_winding_cached

use fossil, only : surface_stl_object, AABB_TREE_SAH_BVH, AABB_TREE_OCTREE
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

logical :: ok(1)

ok(1) = check('src/tests/dragon.stl', 'dragon')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(file_name, label) result(passed)
   character(*), intent(in) :: file_name, label
   logical                  :: passed
   type(surface_stl_object) :: sb, so
   type(vector_R8P)         :: p, bmin, h
   real(R8P)                :: wb, wo, we, eb, eo, ebo
   integer(I4P)             :: i, j, k
   logical                  :: default_ok

   call sb%load_from_file(file_name=file_name, guess_format=.true., aabb_tree_kind=AABB_TREE_SAH_BVH)
   call so%load_from_file(file_name=file_name, guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
   bmin = sb%get_bmin() ; h = sb%get_bmax() - bmin
   bmin = bmin - 0.1_R8P * h ; h = 1.2_R8P * h / 9._R8P
   eb = 0._R8P ; eo = 0._R8P ; ebo = 0._R8P ; default_ok = .true.
   do k = 0, 9 ; do j = 0, 9 ; do i = 0, 9
      p  = vector_R8P(bmin%x + i * h%x, bmin%y + j * h%y, bmin%z + k * h%z)
      we = sb%winding_number(point=p, beta=-1._R8P)
      wb = sb%winding_number(point=p)
      wo = so%winding_number(point=p)
      if (wb /= sb%winding_number(point=p, beta=2._R8P)) default_ok = .false.
      eb  = max(eb,  abs(wb - we))
      eo  = max(eo,  abs(wo - we))
      ebo = max(ebo, abs(wb - wo))
   enddo ; enddo ; enddo
   passed = (eb < 1.e-2_R8P) .and. (eo < 1.e-2_R8P) .and. (ebo < 2.e-2_R8P) .and. default_ok
   print '(A,A,ES10.3,A,ES10.3,A,ES10.3,A,L1,A,L1)', label, ': max |w_fast - w_exact| BVH ', eb, ' octree ', eo, &
         '  BVH vs octree ', ebo, '  default==beta2 ', default_ok, ' ok=', passed
   endfunction check
endprogram fossil_test_winding_cached
