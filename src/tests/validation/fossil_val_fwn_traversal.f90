!< VALIDATION version (issue #45): the full, slow test, run by scripts/run_validation.sh at -O2.
!< The fast regression variant run per push is src/tests/fossil_test_fwn_traversal.f90.
!<
!< FOSSIL, test the iterative fast-winding-number traversal `bvh_winding_number_oac` (issue #29, step 3).
!<
!< On 12^3 grids over the bounding box grown by 10%:
!<   1. exact path (beta2 <= 0) equals the per-facet sum of `facet%solid_angle` to 1e-10 in w
!<      (the traversal rebuilds vertices as v1 + e12, v1 + e13 from the payload);
!<   2. fast path vs exact: |w_fast - w_exact| < 1e-2 at beta = 2 and < 2e-3 at beta = 4, and
!<      the error decreases with beta, on bunny and dragon-fine (measured with the second-order
!<      expansion: bunny 8.7e-3 / 1.7e-3, dragon-fine 9.8e-4 / 7.2e-5). The bound also catches
!<      the far-field sign error of the former implementation (error ~2.0 on bunny);
!<   3. watertight cube-05 (12288 facets): w within 1e-2 of 1 at interior points and of 0 at
!<      exterior points, with the far field active.

program fossil_val_fwn_traversal

use fossil,                  only : surface_stl_object, facet_object
use fossil_aabb_tree_object, only : bvh_winding_number_oac
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter :: FOURPI = 4._R8P * acos(-1._R8P)
logical              :: ok(5)

ok(1) = check_exact('src/tests/bunny.stl',       'bunny      ')
ok(2) = check_exact('src/tests/dragon-fine.stl', 'dragon-fine')
ok(3) = check_fast('src/tests/bunny.stl',        'bunny      ')
ok(4) = check_fast('src/tests/dragon-fine.stl',  'dragon-fine')
ok(5) = check_cube()

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function w_of(s, p, beta) result(w)
   type(surface_stl_object), intent(in) :: s
   type(vector_R8P),         intent(in) :: p
   real(R8P),                intent(in) :: beta
   real(R8P)                            :: w
   w = bvh_winding_number_oac(s%aabb%node, s%aabb%payload, s%wn_payload, p%x, p%y, p%z, sign(beta*beta, beta)) / FOURPI
   endfunction w_of

   function grid_point(s, i, j, k) result(p)
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: i, j, k
   type(vector_R8P)                     :: p, bmin, h
   bmin = s%get_bmin() ; h = s%get_bmax() - bmin
   bmin = bmin - 0.1_R8P * h ; h = 1.2_R8P * h / 11._R8P
   p = vector_R8P(bmin%x + i * h%x, bmin%y + j * h%y, bmin%z + k * h%z)
   endfunction grid_point

   function check_exact(file_name, label) result(passed)
   character(*), intent(in)    :: file_name, label
   logical                     :: passed
   type(surface_stl_object)    :: s
   type(facet_object), pointer :: f
   type(vector_R8P)            :: p
   real(R8P)                   :: we, wf, err
   integer(I4P)                :: i, j, k, m

   call s%load_from_file(file_name=file_name, guess_format=.true.)
   err = 0._R8P
   do k = 0, 11, 3 ; do j = 0, 11, 3 ; do i = 0, 11, 3
      p = grid_point(s, i, j, k)
      we = w_of(s, p, -1._R8P)
      wf = 0._R8P
      do m = 1, s%get_facets_number()
         f => s%facet_at(m)
         wf = wf + f%solid_angle(point=p)
      enddo
      err = max(err, abs(we - wf / FOURPI))
   enddo ; enddo ; enddo
   passed = err < 1.e-10_R8P
   print '(A,A,ES10.3,A,L1)', label, ': exact path vs per-facet sum, max |dw| = ', err, ' ok=', passed
   endfunction check_exact

   function check_fast(file_name, label) result(passed)
   character(*), intent(in) :: file_name, label
   logical                  :: passed
   type(surface_stl_object) :: s
   type(vector_R8P)         :: p
   real(R8P)                :: we, e2, e4
   integer(I4P)             :: i, j, k

   call s%load_from_file(file_name=file_name, guess_format=.true.)
   e2 = 0._R8P ; e4 = 0._R8P
   do k = 0, 11 ; do j = 0, 11 ; do i = 0, 11
      p  = grid_point(s, i, j, k)
      we = w_of(s, p, -1._R8P)
      e2 = max(e2, abs(w_of(s, p, 2._R8P) - we))
      e4 = max(e4, abs(w_of(s, p, 4._R8P) - we))
   enddo ; enddo ; enddo
   passed = (e2 < 1.e-2_R8P) .and. (e4 < 2.e-3_R8P) .and. (e4 < e2)
   print '(A,A,ES10.3,A,ES10.3,A,L1)', label, ': max |w_fast - w_exact|, beta=2: ', e2, '  beta=4: ', e4, ' ok=', passed
   endfunction check_fast

   function check_cube() result(passed)
   logical                  :: passed
   type(surface_stl_object) :: s
   real(R8P)                :: w_in, w_out, e_in, e_out
   integer(I4P)             :: i, j, k

   call s%load_from_file(file_name='src/tests/cube-05.stl', guess_format=.true.)
   e_in = 0._R8P ; e_out = 0._R8P
   do k = 1, 3 ; do j = 1, 3 ; do i = 1, 3
      w_in  = w_of(s, vector_R8P(0.25_R8P * i, 0.25_R8P * j, 0.25_R8P * k), 2._R8P)
      w_out = w_of(s, vector_R8P(1.5_R8P + 0.5_R8P * i, 0.25_R8P * j, 0.25_R8P * k), 2._R8P)
      e_in  = max(e_in,  abs(w_in - 1._R8P))
      e_out = max(e_out, abs(w_out))
   enddo ; enddo ; enddo
   passed = (e_in < 1.e-2_R8P) .and. (e_out < 1.e-2_R8P)
   print '(A,ES10.3,A,ES10.3,A,L1)', 'cube-05: max |w - 1| inside = ', e_in, '  max |w| outside = ', e_out, ' ok=', passed
   endfunction check_cube
endprogram fossil_val_fwn_traversal
