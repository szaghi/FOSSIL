!< Shape Diameter Function — per-facet thickness scalar (issue #18 §1.9, step 1).
!<
!< Reference: Shapira, Shamir & Cohen-Or, *Consistent Mesh Partitioning and
!< Skeletonisation Using the Shape Diameter Function* (Visual Computer 2008).
!< CGAL package: Surface_mesh_segmentation.
!<
!< Pipeline (this module covers step 1 only — smoothing in step 2, GMM in step 3):
!<   For each facet, cast `num_rays` rays from the centroid into a cone of half-
!<   angle `cone_angle_deg/2` around the *inward* normal, collect the first-hit
!<   distance from the opposite surface, and take a robust statistic (median).
!<   This scalar — the "shape diameter" — is approximately the local thickness
!<   of the solid at that facet. Thin features (fillets, fins) get small SDF;
!<   thick body regions get large SDF.

module fossil_sdf
!< Shape Diameter Function — per-facet thickness scalar (issue #18 §1.9).

use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_ray_query, only : ray_hit_t
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private

public :: compute_sdf
public :: SDF_STATUS_OK, SDF_STATUS_BAD_INPUT
public :: SDF_SENTINEL
public :: SDF_DEFAULT_NUM_RAYS, SDF_DEFAULT_CONE_DEG
public :: SDF_MIN_HIT_FRACTION

integer(I4P), parameter :: SDF_STATUS_OK         = 0_I4P
integer(I4P), parameter :: SDF_STATUS_BAD_INPUT  = 1_I4P  !< Empty surface, or num_rays < 1, or invalid cone angle.

real(R8P),    parameter :: SDF_SENTINEL          = -1._R8P  !< SDF for facets with too few hits (degenerate shells).

integer(I4P), parameter :: SDF_DEFAULT_NUM_RAYS  = 30_I4P   !< Shapira's default.
real(R8P),    parameter :: SDF_DEFAULT_CONE_DEG  = 120._R8P !< Shapira's default cone half-aperture is 60°, full angle 120°.

real(R8P),    parameter :: SDF_MIN_HIT_FRACTION  = 0.5_R8P  !< Below this hit-rate, the SDF for that facet is the sentinel.

real(R8P),    parameter :: PI_R8P                = 3.14159265358979323846_R8P
real(R8P),    parameter :: GOLDEN_ANGLE          = PI_R8P * (3._R8P - sqrt(5._R8P))  !< Fibonacci-spiral angular step.
real(R8P),    parameter :: ORIGIN_OFFSET_REL     = 1.0e-6_R8P  !< Ray origin offset = this * bbox_diagonal.

contains

   subroutine compute_sdf(facet, tree, bmin, bmax, sdf, num_rays, cone_angle_deg, status)
   !< Compute the per-facet Shape Diameter Function (issue #18 §1.9 step 1).
   !<
   !< Operates on a `facet(:)` array + `aabb_tree` (no `surface_stl_object`
   !< import — feature modules in this codebase don't depend on the surface
   !< type, to avoid a use cycle; the surface module wires the high-level TBP).
   !<
   !< Output `sdf(1:nf)` corresponds to `facet(1:nf)`. A facet whose ray cone
   !< returns hits on fewer than `SDF_MIN_HIT_FRACTION` of the rays gets
   !< `SDF_SENTINEL` — downstream clustering must exclude those. Documented
   !< contract for "degenerate shell" facets where interior thickness is
   !< ill-defined (single-face shells, near-holes, isolated patches).
   type(facet_object),     intent(in)              :: facet(:)         !< Surface facets.
   type(aabb_tree_object), intent(in)              :: tree             !< Built AABB tree over `facet`.
   type(vector_R8P),       intent(in)              :: bmin             !< Surface bbox minimum.
   type(vector_R8P),       intent(in)              :: bmax             !< Surface bbox maximum.
   real(R8P), allocatable, intent(out)             :: sdf(:)           !< Per-facet SDF scalar (length = nf).
   integer(I4P),           intent(in),    optional :: num_rays         !< Rays per facet (default SDF_DEFAULT_NUM_RAYS).
   real(R8P),              intent(in),    optional :: cone_angle_deg   !< Full cone aperture, degrees (default SDF_DEFAULT_CONE_DEG).
   integer(I4P),           intent(out),   optional :: status           !< Status code.
   integer(I4P)                                    :: nf               !< Facet count.
   integer(I4P)                                    :: nrays            !< Resolved num_rays.
   real(R8P)                                       :: cone_deg         !< Resolved cone angle.
   real(R8P)                                       :: bdiag            !< Bbox diagonal (for offset and max ray length).
   real(R8P)                                       :: offset           !< Ray-origin offset along -normal.
   real(R8P)                                       :: max_t            !< Cap on ray distance.
   type(vector_R8P), allocatable                   :: dirs(:)          !< Cone-ray directions for current facet.
   type(ray_hit_t)                                 :: hit              !< Per-ray first-hit record.
   type(vector_R8P)                                :: ray_origin       !< Per-ray origin (centroid - offset * normal).
   type(vector_R8P)                                :: inward           !< -normal versor.
   logical                                         :: has_hit          !< Per-ray flag.
   real(R8P), allocatable                          :: hits_t(:)        !< Per-facet hits' t buffer (size up to nrays).
   integer(I4P)                                    :: f, r, n_hits     !< Counters.

   if (present(status)) status = SDF_STATUS_OK

   nf = size(facet, kind=I4P)
   if (nf == 0_I4P) then
      allocate(sdf(0))
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif
   nrays    = SDF_DEFAULT_NUM_RAYS;  if (present(num_rays))       nrays    = num_rays
   cone_deg = SDF_DEFAULT_CONE_DEG;  if (present(cone_angle_deg)) cone_deg = cone_angle_deg
   if (nrays < 1_I4P .or. cone_deg <= 0._R8P .or. cone_deg > 180._R8P) then
      allocate(sdf(nf))
      sdf = SDF_SENTINEL
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif

   bdiag = sqrt((bmax%x - bmin%x)**2 + (bmax%y - bmin%y)**2 + (bmax%z - bmin%z)**2)
   offset = ORIGIN_OFFSET_REL * bdiag
   max_t  = bdiag  ! ray can't traverse further than the bbox diagonal in a convex sense

   allocate(sdf(nf))
   allocate(hits_t(nrays))

   do f = 1_I4P, nf
      inward = -1._R8P * facet(f)%normal
      ray_origin = facet(f)%centroid + offset * inward    ! step just inside the surface
      call cone_directions(axis=inward, cone_deg=cone_deg, n=nrays, dirs=dirs)
      n_hits = 0_I4P
      do r = 1_I4P, nrays
         call tree%intersect_ray_first_tree(facet=facet, ray_origin=ray_origin, &
                                            ray_direction=dirs(r), hit=hit, has_hit=has_hit)
         if (.not. has_hit) cycle
         if (hit%t < 0._R8P) cycle
         if (hit%t > max_t)  cycle
         n_hits = n_hits + 1_I4P
         hits_t(n_hits) = hit%t
      enddo
      if (real(n_hits, R8P) / real(nrays, R8P) < SDF_MIN_HIT_FRACTION) then
         sdf(f) = SDF_SENTINEL
      else
         sdf(f) = median(hits_t(1:n_hits))
      endif
   enddo

   deallocate(hits_t)
   if (allocated(dirs)) deallocate(dirs)
   endsubroutine compute_sdf

   pure subroutine cone_directions(axis, cone_deg, n, dirs)
   !< Generate `n` directions uniformly distributed in a cone of full aperture
   !< `cone_deg` around `axis` via Fibonacci-spiral sampling.
   !<
   !< Deterministic (no PRNG, same inputs → same outputs), low-discrepancy
   !< (no clumping at the pole or the rim), unit-length. Caller must pass a
   !< unit-length `axis`; we don't normalize defensively — Shapira's pipeline
   !< always passes the facet normal which is already unit.
   type(vector_R8P),              intent(in)  :: axis      !< Cone central axis (unit).
   real(R8P),                     intent(in)  :: cone_deg  !< Full cone aperture in degrees.
   integer(I4P),                  intent(in)  :: n         !< Number of directions.
   type(vector_R8P), allocatable, intent(out) :: dirs(:)   !< Output cone directions (unit, length n).
   real(R8P)                                  :: half_rad  !< cone_deg/2 in radians.
   real(R8P)                                  :: cos_half  !< Cosine of half-aperture (caps polar angle).
   real(R8P)                                  :: t, phi, ct, st  !< Sampling parameters.
   real(R8P)                                  :: cos_th, sin_th  !< Polar cosine, sine.
   type(vector_R8P)                           :: e1, e2          !< Tangent basis perpendicular to axis.
   integer(I4P)                               :: i

   if (allocated(dirs)) deallocate(dirs)
   allocate(dirs(n))
   half_rad = cone_deg * 0.5_R8P * PI_R8P / 180._R8P
   cos_half = cos(half_rad)
   call ortho_basis(axis=axis, e1=e1, e2=e2)
   do i = 1_I4P, n
      ! Map index i ∈ [1, n] to t ∈ [0, 1]; use stratified midpoints to avoid the pole.
      t = (real(i, R8P) - 0.5_R8P) / real(n, R8P)
      cos_th = 1._R8P - t * (1._R8P - cos_half)   ! cos uniformly in [cos_half, 1] → uniform on cone cap area
      sin_th = sqrt(max(0._R8P, 1._R8P - cos_th * cos_th))
      phi = real(i, R8P) * GOLDEN_ANGLE
      ct = cos(phi)
      st = sin(phi)
      dirs(i) = cos_th * axis + sin_th * (ct * e1 + st * e2)
   enddo
   endsubroutine cone_directions

   pure subroutine ortho_basis(axis, e1, e2)
   !< Build a right-handed orthonormal basis (axis, e1, e2). Uses Hughes-Möller
   !< stable construction (no branch on near-pole degeneracy).
   type(vector_R8P), intent(in)  :: axis  !< Unit input axis.
   type(vector_R8P), intent(out) :: e1    !< First tangent (unit, perpendicular to axis).
   type(vector_R8P), intent(out) :: e2    !< Second tangent (unit, perpendicular to both).
   real(R8P)                     :: ax, ay, az
   real(R8P)                     :: nrm

   ax = abs(axis%x); ay = abs(axis%y); az = abs(axis%z)
   ! Pick the axis-aligned vector least parallel to `axis` to seed e1.
   if (ax <= ay .and. ax <= az) then
      e1 = vector_R8P(0._R8P, -axis%z, axis%y)
   elseif (ay <= ax .and. ay <= az) then
      e1 = vector_R8P(-axis%z, 0._R8P, axis%x)
   else
      e1 = vector_R8P(-axis%y, axis%x, 0._R8P)
   endif
   nrm = sqrt(e1%dotproduct(rhs=e1))
   if (nrm > 0._R8P) e1 = e1 / nrm
   e2 = axis%crossproduct(rhs=e1)
   endsubroutine ortho_basis

   pure function median(x) result(m)
   !< Exact median of `x` via a copy-and-sort. n is at most ~30 (one cone per
   !< facet), so insertion sort beats heapsort's overhead. For larger n switch
   !< to nth_element-style partition (out of scope here).
   real(R8P), intent(in) :: x(:)
   real(R8P)             :: m
   real(R8P), allocatable :: tmp(:)
   real(R8P)              :: key
   integer(I4P)           :: n, i, j

   n = size(x, kind=I4P)
   if (n == 0_I4P) then
      m = 0._R8P
      return
   endif
   allocate(tmp(n))
   tmp = x
   do i = 2_I4P, n
      key = tmp(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (tmp(j) <= key) exit
         tmp(j + 1_I4P) = tmp(j)
         j = j - 1_I4P
      enddo
      tmp(j + 1_I4P) = key
   enddo
   if (mod(n, 2_I4P) == 1_I4P) then
      m = tmp((n + 1_I4P) / 2_I4P)
   else
      m = 0.5_R8P * (tmp(n / 2_I4P) + tmp(n / 2_I4P + 1_I4P))
   endif
   deallocate(tmp)
   endfunction median

endmodule fossil_sdf
