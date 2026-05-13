!< FOSSIL test: per-facet Shape Diameter Function (issue #18 §1.9, step 1).
!<
!< Public contract for `compute_sdf` that step-N implementations must honour:
!<   - returns `sdf(1:nf)` aligned to `facet(1:nf)`,
!<   - facets with too few cone-ray hits get `SDF_SENTINEL = -1`,
!<   - status = SDF_STATUS_OK on valid input, SDF_STATUS_BAD_INPUT otherwise.
!<
!< Step 1 covers invariants 1..3 (per-facet SDF). Step 2 adds invariant 4
!< for the dual-graph Laplacian smoothing pass. Step 3 adds invariants 5..6
!< for the GMM clustering and `surface%segment_sdf` capstone TBP:
!<   1. Cube SDF is roughly constant. cube.stl spans [0,1]^3, so every
!<      interior cone ray hits the opposite face at distance ~1.0 (with some
!<      spread from cone aperture and triangulation discretization). Assert:
!<      no sentinels, mean SDF ≈ 1.0, stddev/mean ≤ 0.30 (loose: cone rays
!<      to corners traverse longer paths and inflate spread).
!<   2. Sphere SDF is roughly constant ≈ 2R. From the center each radial ray
!<      crosses to the antipode at distance 2R. Cone rays hit shorter chords,
!<      so the median is ≤ 2R but should still be tightly clustered.
!<      (Sphere built from §1.5 marching cubes on an analytic SDF.)
!<   3. Two-scale composite (small cube inside a larger one, disjoint) →
!<      bimodal SDF distribution: small-cube facets cluster near 0.3, big-cube
!<      facets cluster near 1.0. Assert: at least 80% of facets are within
!<      ±20% of one of the two expected values.
!<   4. Smoothing reduces variance. On the bunny.stl SDF, two passes of
!<      Laplacian smoothing (lambda = 0.5, the default) must strictly reduce
!<      the per-facet variance — the whole point of the pass is to attenuate
!<      cone-ray Monte-Carlo noise without erasing real geometric structure.
!<      Assert: stddev(smoothed) < stddev(raw), with mean approximately
!<      preserved (mass-conservation property of symmetric Laplacian
!<      averaging on a closed mesh).
!<   5. Sphere + k=1 → all facets get label 1. Trivial cluster-path test:
!<      with one cluster GMM degenerates to "every valid facet to label 1"
!<      regardless of SDF distribution. (We avoid the obvious "uniform-input
!<      collapses to one cluster with k=2" assertion: GMM is a mode-finder,
!<      not a mode-merger, so it WILL split a tightly-clustered distribution
!<      into two near-identical components if asked. That's intended GMM
!<      behaviour, not a bug, and a graph-cut post-pass would be needed to
!<      enforce true cluster collapse — see the deferred §1.9b follow-up.)
!<   6. Bunny + k=4 capstone via `surface%segment_sdf`. Asserts: status OK,
!<      all labels in [0, 4], no surprise sentinels (bunny is closed so the
!<      cone-ray hit rate should be ~100%; a high sentinel fraction would
!<      indicate something broken in the ray cast). Verifies the public TBP
!<      end-to-end (compute_sdf → smooth_sdf → GMM → labels).

program fossil_test_sdf

use fossil, only : surface_stl_object, extract_isosurface, &
                   SDF_STATUS_OK, SDF_LABEL_UNASSIGNED
use fossil_facet_object, only : facet_object
use fossil_sdf, only : compute_sdf, smooth_sdf, SDF_SENTINEL
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, sphere, composite, small_cube, bunny
type(facet_object), allocatable :: sphere_facets(:)
real(R8P), allocatable :: sdf(:), sdf_smoothed(:), values(:,:,:)
integer(I4P), allocatable :: labels(:)
type(vector_R8P) :: bmin, bmax
real(R8P) :: x, y, z, dx
real(R8P) :: mean_sdf, stddev_sdf, frac_sentinel
real(R8P) :: mean_raw, stddev_raw, mean_smooth, stddev_smooth
real(R8P) :: frac_in_cluster_a, frac_in_cluster_b, target_a, target_b
integer(I4P) :: status, i, j, k, f, nf, n_sentinel, n_a, n_b, n_dummy
integer(I4P) :: cube_first_label, n_label_diff, lab_min, lab_max, n_unassigned
logical :: are_tests_passed(6)

integer(I4P), parameter :: N_GRID = 32_I4P
real(R8P),    parameter :: TOL_BIMODAL = 0.20_R8P  ! ±20% of expected SDF value

are_tests_passed = .false.

! ---- Invariant 1: cube SDF ≈ 1, low spread.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call compute_sdf(facet=ref_facets(cube), tree=cube%aabb, &
                 bmin=cube%get_bmin(), bmax=cube%get_bmax(), &
                 sdf=sdf, status=status)
nf = size(sdf, kind=I4P)
call summary(sdf=sdf, mean=mean_sdf, stddev=stddev_sdf, n_sentinel=n_sentinel)
frac_sentinel = real(n_sentinel, R8P) / real(nf, R8P)
print '(A,I0,A,F8.4,A,F8.4,A,F6.3)', 'inv 1 (cube): nf=', nf, ' mean=', mean_sdf, &
      ' stddev=', stddev_sdf, ' sentinel_frac=', frac_sentinel
are_tests_passed(1) = (status == SDF_STATUS_OK)         .and. &
                      (frac_sentinel == 0._R8P)         .and. &
                      (abs(mean_sdf - 1._R8P) < 0.20_R8P) .and. &
                      (stddev_sdf / mean_sdf < 0.30_R8P)

! ---- Invariant 2: sphere SDF ≈ 2R = 2.0, low relative spread.
allocate(values(N_GRID, N_GRID, N_GRID))
dx = 4._R8P / real(N_GRID - 1, R8P)
do k = 1_I4P, N_GRID
   z = -2._R8P + (k - 1) * dx
   do j = 1_I4P, N_GRID
      y = -2._R8P + (j - 1) * dx
      do i = 1_I4P, N_GRID
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1._R8P  ! unit sphere
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=sphere_facets, status=status)
call sphere%adopt_facets(facets=sphere_facets)
deallocate(values)
call compute_sdf(facet=ref_facets(sphere), tree=sphere%aabb, &
                 bmin=sphere%get_bmin(), bmax=sphere%get_bmax(), &
                 sdf=sdf, status=status)
nf = size(sdf, kind=I4P)
call summary(sdf=sdf, mean=mean_sdf, stddev=stddev_sdf, n_sentinel=n_sentinel)
frac_sentinel = real(n_sentinel, R8P) / real(nf, R8P)
print '(A,I0,A,F8.4,A,F8.4,A,F6.3)', 'inv 2 (sphere): nf=', nf, ' mean=', mean_sdf, &
      ' stddev=', stddev_sdf, ' sentinel_frac=', frac_sentinel
! Median chord through a sphere (uniformly distributed on the cap) is
! significantly less than the diameter; expect mean SDF in [1.0, 2.0].
are_tests_passed(2) = (status == SDF_STATUS_OK)              .and. &
                      (frac_sentinel < 0.05_R8P)             .and. &
                      (mean_sdf > 1.0_R8P)                   .and. &
                      (mean_sdf < 2.1_R8P)                   .and. &
                      (stddev_sdf / mean_sdf < 0.40_R8P)

! ---- Invariant 3: two-scale composite (big cube + far-translated small cube).
!      Two disjoint solids in one surface. SDF queries on big-cube facets
!      should report ~1.0; on small-cube facets ~0.3.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call small_cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call small_cube%resize(factor=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P))
call small_cube%translate(delta=vector_R8P(3._R8P, 0._R8P, 0._R8P))  ! disjoint, outside big cube
composite = cube
call composite%merge_solids(other=small_cube)
call compute_sdf(facet=ref_facets(composite), tree=composite%aabb, &
                 bmin=composite%get_bmin(), bmax=composite%get_bmax(), &
                 sdf=sdf, status=status)
nf = size(sdf, kind=I4P)
target_a = 1.0_R8P
target_b = 0.3_R8P
n_a = 0_I4P; n_b = 0_I4P
do f = 1_I4P, nf
   if (sdf(f) == SDF_SENTINEL) cycle
   if (abs(sdf(f) - target_a) <= TOL_BIMODAL * target_a) n_a = n_a + 1_I4P
   if (abs(sdf(f) - target_b) <= TOL_BIMODAL * target_b) n_b = n_b + 1_I4P
enddo
frac_in_cluster_a = real(n_a, R8P) / real(nf, R8P)
frac_in_cluster_b = real(n_b, R8P) / real(nf, R8P)
print '(A,I0,A,F6.3,A,F6.3)', 'inv 3 (bimodal): nf=', nf, &
      ' frac_near_1.0=', frac_in_cluster_a, ' frac_near_0.3=', frac_in_cluster_b
are_tests_passed(3) = (status == SDF_STATUS_OK) .and. &
                      ((frac_in_cluster_a + frac_in_cluster_b) >= 0.80_R8P) .and. &
                      (frac_in_cluster_a > 0._R8P) .and. &
                      (frac_in_cluster_b > 0._R8P)

! ---- Invariant 4: smoothing strictly reduces variance on bunny.
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call compute_sdf(facet=ref_facets(bunny), tree=bunny%aabb, &
                 bmin=bunny%get_bmin(), bmax=bunny%get_bmax(), &
                 sdf=sdf, status=status)
allocate(sdf_smoothed(size(sdf, kind=I4P)))
sdf_smoothed = sdf
call smooth_sdf(facet=ref_facets(bunny), sdf=sdf_smoothed, status=status)
call summary(sdf=sdf,          mean=mean_raw,    stddev=stddev_raw,    n_sentinel=n_dummy)
call summary(sdf=sdf_smoothed, mean=mean_smooth, stddev=stddev_smooth, n_sentinel=n_dummy)
print '(A,F8.4,A,F8.4,A,F8.4,A,F8.4)', &
      'inv 4 (smoothing): raw mean=', mean_raw, ' stddev=', stddev_raw, &
      '  smooth mean=', mean_smooth, ' stddev=', stddev_smooth
are_tests_passed(4) = (status == SDF_STATUS_OK)                         .and. &
                      (stddev_smooth < stddev_raw)                      .and. &
                      (abs(mean_smooth - mean_raw) < 0.05_R8P * mean_raw)
deallocate(sdf_smoothed)

! ---- Invariant 5: sphere + k=1 → every valid facet gets label 1.
!      Trivial cluster-path test (with k=1, GMM degenerates to one label).
call sphere%segment_sdf(facet_labels=labels, num_clusters=1_I4P, status=status)
nf = size(labels, kind=I4P)
n_a = count(labels == 1_I4P)
n_unassigned = count(labels == SDF_LABEL_UNASSIGNED)
print '(A,I0,A,I0,A,I0)', 'inv 5 (sphere k=1): nf=', nf, ' n_label_1=', n_a, &
      ' n_unassigned=', n_unassigned
are_tests_passed(5) = (status == SDF_STATUS_OK) .and. &
                      (n_a + n_unassigned == nf) .and. &
                      (n_a > 0_I4P)

! ---- Invariant 6: bunny + k=4 capstone via TBP, all labels in [0, 4].
call bunny%segment_sdf(facet_labels=labels, num_clusters=4_I4P, status=status)
nf = size(labels, kind=I4P)
lab_min = minval(labels)
lab_max = maxval(labels)
n_unassigned = count(labels == SDF_LABEL_UNASSIGNED)
print '(A,I0,A,I0,A,I0,A,I0)', 'inv 6 (bunny k=4): nf=', nf, ' lab_min=', lab_min, &
      ' lab_max=', lab_max, ' n_unassigned=', n_unassigned
are_tests_passed(6) = (status == SDF_STATUS_OK)             .and. &
                      (lab_min >= 0_I4P)                    .and. &
                      (lab_max <= 4_I4P)                    .and. &
                      (real(n_unassigned, R8P) / real(nf, R8P) < 0.05_R8P)

print '(A,6L2)', 'per-case results: ', are_tests_passed
print '(A,L1)', 'Are all tests passed? ', all(are_tests_passed)

contains

   function ref_facets(s) result(arr)
   !< Read-only handle to a surface's facet array.
   type(surface_stl_object), target, intent(in) :: s
   type(facet_object),       pointer            :: arr(:)

   arr => s%facets_ref()
   endfunction ref_facets

   pure subroutine summary(sdf, mean, stddev, n_sentinel)
   !< Summary statistics of a SDF array, excluding sentinels.
   real(R8P),    intent(in)  :: sdf(:)
   real(R8P),    intent(out) :: mean, stddev
   integer(I4P), intent(out) :: n_sentinel
   integer(I4P)              :: ii, n_valid
   real(R8P)                 :: s, ss

   n_sentinel = 0_I4P
   n_valid = 0_I4P
   s = 0._R8P
   ss = 0._R8P
   do ii = 1_I4P, size(sdf, kind=I4P)
      if (sdf(ii) == SDF_SENTINEL) then
         n_sentinel = n_sentinel + 1_I4P
         cycle
      endif
      n_valid = n_valid + 1_I4P
      s  = s + sdf(ii)
      ss = ss + sdf(ii) * sdf(ii)
   enddo
   if (n_valid == 0_I4P) then
      mean = 0._R8P
      stddev = 0._R8P
   else
      mean = s / real(n_valid, R8P)
      stddev = sqrt(max(0._R8P, ss / real(n_valid, R8P) - mean * mean))
   endif
   endsubroutine summary

endprogram fossil_test_sdf
