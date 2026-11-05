!< FOSSIL test: Marching Cubes isosurface extraction (issue #18 §1.5).
!<
!< Five invariants:
!<
!<   1. Sphere SDF — sample `f(x,y,z) = sqrt(x^2+y^2+z^2) - 1.0` on a 32^3
!<      grid in [-2, 2]^3, extract iso=0. Result volume should equal
!<      4/3 π = 4.189 within ~5% (grid-resolution slack).
!<   2. Cube SDF — sample `f(x,y,z) = max(|x|, |y|, |z|) - 0.5` on a 32^3
!<      grid in [-1, 1]^3, extract iso=0. Result volume = 1.0.
!<   3. Off-iso sphere — same sphere SDF, extract at iso=0.3. The level
!<      surface is at radius 0.7, so volume should equal 4/3 π * 0.7^3.
!<   4. Roundtrip via distance field — load cube.stl, call
!<      `surface%resample_via_distance_field(32)`, assert the resulting
!<      surface has volume ≈ 1.0. Tests integration with the existing
!<      `distance` API and the §1.5 step 5 helper.
!<   5. Watertight straight from the extraction (issue #48): spheres of radius 1, 0.96, 0.8137, 1.0731 and 0.5 on
!<      31^3 and 32^3 grids over [-1.5, 1.5]^3, adopted with no repair, have no boundary edge, no non-manifold edge
!<      and no zero-area facet. Before #48 radius 1 on 31^3 had 336 boundary edges, 12 non-manifold edges and 26
!<      zero-area facets: each cell interpolated a shared grid edge in its own direction, so the copies of a crossing
!<      differed in the last bits and the exact weld left them apart.

program fossil_test_marching_cubes

use fossil, only : surface_stl_object, &
                   extract_isosurface, MC_STATUS_OK
use fossil_facet_object, only : facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: sphere_surface, cube_surface, offset_surface, resampled
type(facet_object), allocatable :: facets(:)
real(R8P), allocatable   :: values(:, :, :)
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: vol, vol_target
integer(I4P)             :: status, i, j, k, n
real(R8P)                :: x, y, z, dx
type(surface_stl_object) :: mc
type(facet_object), pointer :: fp
real(R8P)                :: radius(5)
integer(I4P)             :: ir, ig, nbad, f
logical                  :: are_tests_passed(5)

integer(I4P), parameter :: N_GRID         = 32_I4P
real(R8P),    parameter :: PI             = 3.141592653589793_R8P
real(R8P),    parameter :: TOL_VOLUME_5P  = 0.05_R8P  ! 5% relative tolerance
real(R8P),    parameter :: TOL_VOLUME_10P = 0.10_R8P  ! 10% for sharp-edged cube SDF (MC stair-steps the corners)

are_tests_passed = .false.

! ----- 1. Sphere SDF: f = r - 1.0; iso=0 → unit sphere.
n  = N_GRID
dx = 4._R8P / real(n - 1, R8P)
allocate(values(n, n, n))
do k = 1, n
   z = -2._R8P + (k - 1) * dx
   do j = 1, n
      y = -2._R8P + (j - 1) * dx
      do i = 1, n
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1.0_R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=facets, status=status)
call sphere_surface%adopt_facets(facets=facets)
vol = sphere_surface%get_volume()
vol_target = 4._R8P / 3._R8P * PI
are_tests_passed(1) = (status == MC_STATUS_OK .and. &
                       abs(vol - vol_target) <= TOL_VOLUME_5P * vol_target)
print '(A,I0,A,ES14.6,A,ES14.6)', 'sphere: status=', status, '  vol=', vol, '  target=', vol_target
deallocate(values)

! ----- 2. Cube SDF: f = max(|x|, |y|, |z|) - 0.5; iso=0 → unit cube.
n  = N_GRID
dx = 2._R8P / real(n - 1, R8P)
allocate(values(n, n, n))
do k = 1, n
   z = -1._R8P + (k - 1) * dx
   do j = 1, n
      y = -1._R8P + (j - 1) * dx
      do i = 1, n
         x = -1._R8P + (i - 1) * dx
         values(i, j, k) = max(abs(x), abs(y), abs(z)) - 0.5_R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-1._R8P, -1._R8P, -1._R8P)
bmax = vector_R8P( 1._R8P,  1._R8P,  1._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                        surface=facets, status=status)
call cube_surface%adopt_facets(facets=facets)
vol = cube_surface%get_volume()
vol_target = 1._R8P
are_tests_passed(2) = (status == MC_STATUS_OK .and. &
                       abs(vol - vol_target) <= TOL_VOLUME_10P * vol_target)
print '(A,I0,A,ES14.6,A,ES14.6)', 'cube SDF: status=', status, '  vol=', vol, '  target=', vol_target
deallocate(values)

! ----- 3. Off-iso sphere: same sphere SDF, iso=0.3 → sphere of radius 0.7.
n  = N_GRID
dx = 4._R8P / real(n - 1, R8P)
allocate(values(n, n, n))
do k = 1, n
   z = -2._R8P + (k - 1) * dx
   do j = 1, n
      y = -2._R8P + (j - 1) * dx
      do i = 1, n
         x = -2._R8P + (i - 1) * dx
         values(i, j, k) = sqrt(x**2 + y**2 + z**2) - 1.0_R8P
      enddo
   enddo
enddo
bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=-0.3_R8P, &
                        surface=facets, status=status)
call offset_surface%adopt_facets(facets=facets)
vol = offset_surface%get_volume()
vol_target = 4._R8P / 3._R8P * PI * 0.7_R8P**3
are_tests_passed(3) = (status == MC_STATUS_OK .and. &
                       abs(vol - vol_target) <= TOL_VOLUME_5P * vol_target)
print '(A,I0,A,ES14.6,A,ES14.6)', 'off-iso: status=', status, '  vol=', vol, '  target=', vol_target
deallocate(values)

! ----- 4. Roundtrip: cube.stl → distance field → MC at iso=0 → cube_resampled.
call cube_surface%destroy
call cube_surface%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call cube_surface%resample_via_distance_field(resolution=32_I4P, surface_out=resampled, status=status)
vol = resampled%get_volume()
vol_target = 1._R8P
are_tests_passed(4) = (status == MC_STATUS_OK .and. &
                       abs(vol - vol_target) <= TOL_VOLUME_10P * vol_target)
print '(A,I0,A,ES14.6,A,ES14.6)', 'roundtrip: status=', status, '  vol=', vol, '  target=', vol_target

! ----- 5. Watertight, manifold, no zero-area facet, with no repair.
radius = [1._R8P, 0.96_R8P, 0.8137_R8P, 1.0731_R8P, 0.5_R8P]
are_tests_passed(5) = .true.
do ig = 31, 32
   if (allocated(values)) deallocate(values)
   allocate(values(ig, ig, ig))
   dx = 3._R8P / (ig - 1)
   do ir = 1, size(radius)
      do k = 1, ig ; do j = 1, ig ; do i = 1, ig
         values(i, j, k) = norm2(-1.5_R8P + dx * ([i, j, k] - 1)) - radius(ir)
      enddo ; enddo ; enddo
      call extract_isosurface(values=values, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                              bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
      call mc%adopt_facets(facets=facets)
      nbad = 0
      do f = 1, mc%get_facets_number()
         fp => mc%facet_at(f)
         nbad = nbad + count(fp%fcon_edge == 0_I4P)
         if (.not. (fp%area() > 0._R8P)) nbad = nbad + 1
      enddo
      nbad = nbad + mc%get_non_manifold_edges_number()
      if (nbad /= 0 .or. .not. mc%is_watertight()) then
         are_tests_passed(5) = .false.
         print '(A,I0,A,F7.4,A,I0)', 'NOT watertight: grid ', ig, ' radius ', radius(ir), ' defects ', nbad
      endif
   enddo
enddo
print '(A,L1)', 'spheres watertight with no repair (2 grids x 5 radii): ', are_tests_passed(5)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

endprogram fossil_test_marching_cubes
