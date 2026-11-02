!< FOSSIL, assert the ray-intersection sign agrees with the pseudo-normal sign on watertight meshes.
!<
!< SIGN_RAY_INTERSECTIONS decides inside/outside by the parity of ray-facet crossings
!< (majority of three rays). With axis-aligned rays, a ray through a shared edge was
!< counted once per incident triangle, so the parity flipped: the unit cube's centre
!< came out "outside". On a watertight, outward-oriented mesh the pseudo-normal sign
!< is the reference. Grids over cube, cube-01, cube-02 and a cube-02 rotated by an
!< arbitrary angle, including lattice-aligned points (the worst case for axis rays);
!< points within 1e-9 of the surface are skipped (both signs are ambiguous there).

program fossil_test_ray_sign

use fossil, only : surface_stl_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

logical :: ok(5) !< Per-case results.

ok(1) = check_centre()
ok(2) = check_mesh('src/tests/cube.stl',    .false., 'cube           ')
ok(3) = check_mesh('src/tests/cube-01.stl', .false., 'cube-01        ')
ok(4) = check_mesh('src/tests/cube-02.stl', .false., 'cube-02        ')
ok(5) = check_mesh('src/tests/cube-02.stl', .true.,  'cube-02 rotated')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check_centre() result(passed)
   !< The unit cube's centre must be inside by ray intersections.
   logical                  :: passed !< Check result.
   type(surface_stl_object) :: cube   !< Unit cube.
   real(R8P)                :: d      !< Signed distance by ray intersections.

   call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call cube%sanitize
   d = cube%distance(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), is_signed=.true., is_square_root=.true., &
                     sign_algorithm=SIGN_RAY_INTERSECTIONS)
   passed = (d < 0._R8P) .and. cube%is_point_inside(point=vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), &
                                                    sign_algorithm=SIGN_RAY_INTERSECTIONS)
   print '(A,F8.4,A,L1)', 'cube centre by ray intersections: d=', d, ' ok=', passed
   endfunction check_centre

   function check_mesh(file_name, rotated, label) result(passed)
   !< Compare ray-intersection and pseudo-normal signs on a 13^3 grid over the grown bbox.
   character(*), intent(in) :: file_name  !< Mesh.
   logical,      intent(in) :: rotated    !< Rotate the mesh by an arbitrary angle first.
   character(*), intent(in) :: label      !< Case label.
   logical                  :: passed     !< Check result.
   type(surface_stl_object) :: surface    !< Surface.
   type(vector_R8P)         :: bmin, bmax, h, p !< Grid geometry and query point.
   real(R8P)                :: d_pn, d_ri !< Signed distances (pseudo-normal, ray intersections).
   integer(I4P), parameter  :: N = 13_I4P !< Grid points per direction.
   integer(I4P)             :: i, j, k    !< Counters.
   integer(I4P)             :: n_checked, n_mismatch !< Statistics.

   call surface%load_from_file(file_name=file_name, guess_format=.true.)
   call surface%sanitize
   if (rotated) then
      call surface%rotate(axis=ex_R8P + 0.7_R8P * ey_R8P + 0.3_R8P * ez_R8P, angle=0.61_R8P)
      call surface%analyze
   endif
   bmin = surface%get_bmin()
   bmax = surface%get_bmax()
   h    = bmax - bmin
   bmin = bmin - 0.25_R8P * h
   h    = 1.5_R8P * h / real(N - 1, R8P)   ! N = 13 puts grid points on the cube's half/quarter lattice
   n_checked = 0 ; n_mismatch = 0
   do k = 0, N - 1
      do j = 0, N - 1
         do i = 0, N - 1
            p%x = bmin%x + real(i, R8P) * h%x
            p%y = bmin%y + real(j, R8P) * h%y
            p%z = bmin%z + real(k, R8P) * h%z
            d_pn = surface%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_PSEUDO_NORMAL)
            if (abs(d_pn) < 1.e-9_R8P) cycle
            d_ri = surface%distance(point=p, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_RAY_INTERSECTIONS)
            n_checked = n_checked + 1
            if (sign(1._R8P, d_pn) /= sign(1._R8P, d_ri)) n_mismatch = n_mismatch + 1
         enddo
      enddo
   enddo
   passed = (n_mismatch == 0) .and. (n_checked > 0)
   print '(A,A,I0,A,I0,A,L1)', label, ': points checked=', n_checked, ' sign mismatches=', n_mismatch, ' ok=', passed
   endfunction check_mesh
endprogram fossil_test_ray_sign
