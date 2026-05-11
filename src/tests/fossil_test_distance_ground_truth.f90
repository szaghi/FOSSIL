!< FOSSIL, signed-distance ground-truth benchmark.
!<
!< Compares the AABB-accelerated distance against a brute-force O(N) facet scan
!< on identical points and asserts maximum absolute error is within tolerance.
!< Brute force is treated as ground truth because it touches every facet — any
!< discrepancy indicates the tree traversal pruned the true closest facet.
!<
!< The test sweeps a regular grid of points covering the mesh bounding box plus
!< a small outside margin and a few deliberately near-surface points, exercising
!< the three regimes the tree must handle: deep-inside, deep-outside, near-surface.

program fossil_test_distance_ground_truth

use fossil, only : surface_stl_object
use fossil_aabb_tree_object, only : AABB_USE_INDEX, AABB_USE_BRUTE_FORCE
use penf, only : I4P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P

implicit none

! tolerance for "equal" unsigned squared distance — brute force and tree must
! reach the same closest facet, so the only allowed slack is floating-point
! roundoff in the O(N) reduction order.
real(R8P), parameter :: DIST_TOL = 1.0e-10_R8P
! tolerance for the analytical cube SDF check (linear distance, not squared).
real(R8P), parameter :: SDF_TOL  = 1.0e-12_R8P

! Reference meshes — each row probes a different stress on the AABB traversal.
! cube           : closed, axis-aligned, symmetric  — sanity check.
! naca0012-ascii : thin elongated surface  — long, narrow children stress the
!                  "closest-box != contains-nearest" pathology of greedy descent.
! dragon         : detailed organic surface — many siblings with overlapping
!                  bounding boxes near the surface.
character(*), parameter :: STL_FILES(*) = [character(64) :: &
   'src/tests/cube.stl',                                   &
   'src/tests/naca0012-ascii.stl',                         &
   'src/tests/dragon.stl']
integer(I4P), parameter :: REF_LEVELS = 2
integer(I4P), parameter :: NI = 9, NJ = 9, NK = 9   ! 9^3 = 729 query points per mesh
real(R8P),    parameter :: MARGIN = 0.5_R8P         ! outside margin as fraction of bbox

type(surface_stl_object)      :: surface
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_brute(:), d_tree(:)
real(R8P)                     :: max_abs_err, abs_err
integer(I4P)                  :: sign_mismatches
integer(I4P)                  :: n_points, p, worst_p, m
logical                       :: tests_passed
logical                       :: overall_passed

overall_passed = .true.

do m = 1, size(STL_FILES)
   call surface%destroy
   call surface%load_from_file(file_name=trim(STL_FILES(m)), guess_format=.true., aabb_refinement_levels=REF_LEVELS)
   call surface%sanitize
   call surface%analyze(aabb_refinement_levels=REF_LEVELS)

   call build_query_points(points, n_points)
   if (allocated(d_brute)) deallocate(d_brute)
   if (allocated(d_tree))  deallocate(d_tree)
   allocate(d_brute(n_points), d_tree(n_points))

   ! Pass 1 — brute force (ground truth)
   call surface%aabb%set_use_index(AABB_USE_BRUTE_FORCE)
   do p = 1, n_points
      d_brute(p) = surface%distance(point=points(p), is_signed=.true.)
   enddo

   ! Pass 2 — AABB tree
   call surface%aabb%set_use_index(AABB_USE_INDEX)
   do p = 1, n_points
      d_tree(p) = surface%distance(point=points(p), is_signed=.true.)
   enddo

   ! Compare
   max_abs_err = 0._R8P
   sign_mismatches = 0
   worst_p = 0
   do p = 1, n_points
      abs_err = abs(abs(d_tree(p)) - abs(d_brute(p)))
      if (abs_err > max_abs_err) then
         max_abs_err = abs_err
         worst_p = p
      endif
      if (sign(1._R8P, d_tree(p)) /= sign(1._R8P, d_brute(p))) sign_mismatches = sign_mismatches + 1
   enddo

   tests_passed = (max_abs_err < DIST_TOL) .and. (sign_mismatches == 0)
   overall_passed = overall_passed .and. tests_passed

   print '(A,A)',      'mesh:                   ', trim(STL_FILES(m))
   print '(A,I0)',     '  query points:           ', n_points
   print '(A,ES12.5)', '  max |d_tree - d_brute|: ', max_abs_err
   print '(A,I0)',     '  sign mismatches:        ', sign_mismatches
   if (worst_p > 0 .and. max_abs_err >= DIST_TOL) then
      print '(A,3(ES12.5,1X))', '  worst point (x,y,z):    ', &
         points(worst_p)%x, points(worst_p)%y, points(worst_p)%z
      print '(A,ES12.5,A,ES12.5)', '    d_brute=', d_brute(worst_p), '  d_tree=', d_tree(worst_p)
   endif
   print '(A,L1)', '  tree==brute (d^2):      ', tests_passed

   ! Analytical SDF check — applicable only to the unit cube [0,1]^3, which has a
   ! closed-form signed distance. This validates pseudo-normal sign decisions
   ! against an oracle independent of the mesh data structure. Requires
   ! sanitize_normals to produce consistently outward normals — that is what
   ! gives the pseudo-normal sign convention "inside = negative" meaning.
   if (trim(STL_FILES(m)) == 'src/tests/cube.stl') then
      call check_cube_analytical(surface, points, n_points, tests_passed)
      overall_passed = overall_passed .and. tests_passed
      print '(A,L1)', '  cube SDF (analytical):  ', tests_passed
   endif
enddo

print '(A,L1)', 'Are all tests passed? ', overall_passed
if (.not. overall_passed) error stop 1

contains

   subroutine build_query_points(pts, n)
   !< Build a query set covering inside/outside/near-surface regimes.
   type(vector_R8P), allocatable, intent(out) :: pts(:)
   integer(I4P),                  intent(out) :: n
   real(R8P)                                  :: dx, dy, dz, xmin, ymin, zmin, xspan, yspan, zspan
   integer(I4P)                               :: i, j, k, idx

   associate(bmin => surface%get_bmin(), bmax => surface%get_bmax())
      xspan = bmax%x - bmin%x; yspan = bmax%y - bmin%y; zspan = bmax%z - bmin%z
      xmin = bmin%x - MARGIN * xspan
      ymin = bmin%y - MARGIN * yspan
      zmin = bmin%z - MARGIN * zspan
      dx = (1._R8P + 2 * MARGIN) * xspan / (NI - 1)
      dy = (1._R8P + 2 * MARGIN) * yspan / (NJ - 1)
      dz = (1._R8P + 2 * MARGIN) * zspan / (NK - 1)
   end associate

   n = NI * NJ * NK
   allocate(pts(n))
   idx = 0
   do k = 1, NK
      do j = 1, NJ
         do i = 1, NI
            idx = idx + 1
            pts(idx) = (xmin + (i - 1) * dx) * ex_R8P + &
                       (ymin + (j - 1) * dy) * ey_R8P + &
                       (zmin + (k - 1) * dz) * ez_R8P
         enddo
      enddo
   enddo
   endsubroutine build_query_points

   subroutine check_cube_analytical(surf, qpts, nq, passed)
   !< Compare FOSSIL's pseudo-normal signed distance against the closed-form SDF
   !< of the unit cube [0,1]^3. Fails the assertion on any mismatch beyond SDF_TOL.
   type(surface_stl_object), intent(inout) :: surf
   type(vector_R8P),         intent(in)    :: qpts(:)
   integer(I4P),             intent(in)    :: nq
   logical,                  intent(out)   :: passed
   real(R8P)                               :: d_fossil, d_true, err, max_err
   integer(I4P)                            :: q, sign_diff, worst

   max_err = 0._R8P
   sign_diff = 0
   worst = 0
   do q = 1, nq
      d_fossil = surf%distance(point=qpts(q), is_signed=.true., is_square_root=.true.)
      d_true   = cube_sdf(qpts(q))
      err = abs(d_fossil - d_true)
      if (err > max_err) then
         max_err = err
         worst = q
      endif
      if (sign(1._R8P, d_fossil) /= sign(1._R8P, d_true)) sign_diff = sign_diff + 1
   enddo
   passed = (max_err < SDF_TOL) .and. (sign_diff == 0)
   print '(A,ES12.5)', '    max |sdf_fossil - sdf_true|: ', max_err
   print '(A,I0)',     '    cube sign mismatches:        ', sign_diff
   if (worst > 0 .and. max_err >= SDF_TOL) then
      d_fossil = surf%distance(point=qpts(worst), is_signed=.true., is_square_root=.true.)
      d_true   = cube_sdf(qpts(worst))
      print '(A,3(F8.4,1X),A,F10.5,A,F10.5)', '    worst pt=', qpts(worst)%x, qpts(worst)%y, qpts(worst)%z, &
          ' fossil=', d_fossil, ' true=', d_true
   endif
   endsubroutine check_cube_analytical

   pure function cube_sdf(p) result(sd)
   !< Closed-form signed distance to the unit cube [0,1]^3, outside positive.
   !< Standard SDF construction: separate the "outside leg" (clamped distance to
   !< the box) from the "inside leg" (max signed projection on each axis, capped
   !< at 0). Their sum gives a smooth signed distance everywhere except on edges
   !< and corners, where it is the L2 distance.
   type(vector_R8P), intent(in) :: p
   real(R8P)                    :: sd
   real(R8P)                    :: qx, qy, qz
   real(R8P)                    :: ax, ay, az
   real(R8P)                    :: outside_len, inside

   ax = abs(p%x - 0.5_R8P)
   ay = abs(p%y - 0.5_R8P)
   az = abs(p%z - 0.5_R8P)
   qx = max(ax - 0.5_R8P, 0._R8P)
   qy = max(ay - 0.5_R8P, 0._R8P)
   qz = max(az - 0.5_R8P, 0._R8P)
   outside_len = sqrt(qx * qx + qy * qy + qz * qz)
   inside = min(max(ax - 0.5_R8P, ay - 0.5_R8P, az - 0.5_R8P), 0._R8P)
   sd = outside_len + inside
   endfunction cube_sdf

endprogram fossil_test_distance_ground_truth
