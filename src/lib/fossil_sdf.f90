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
public :: smooth_sdf
public :: segment_sdf
public :: SDF_STATUS_OK, SDF_STATUS_BAD_INPUT
public :: SDF_SENTINEL
public :: SDF_LABEL_UNASSIGNED
public :: SDF_DEFAULT_NUM_RAYS, SDF_DEFAULT_CONE_DEG
public :: SDF_DEFAULT_SMOOTHING_LAMBDA, SDF_DEFAULT_SMOOTHING_ITERATIONS
public :: SDF_DEFAULT_NUM_CLUSTERS
public :: SDF_MIN_HIT_FRACTION

integer(I4P), parameter :: SDF_STATUS_OK         = 0_I4P
integer(I4P), parameter :: SDF_STATUS_BAD_INPUT  = 1_I4P  !< Empty surface, or num_rays < 1, or invalid cone angle.

real(R8P),    parameter :: SDF_SENTINEL          = -1._R8P  !< SDF for facets with too few hits (degenerate shells).
integer(I4P), parameter :: SDF_LABEL_UNASSIGNED  = 0_I4P    !< Label for sentinel-SDF facets in `segment_sdf` output.

integer(I4P), parameter :: SDF_DEFAULT_NUM_CLUSTERS = 4_I4P  !< Mechanical-CAD-friendly default; CGAL uses 5.

integer(I4P), parameter :: SDF_DEFAULT_NUM_RAYS  = 30_I4P   !< Shapira's default.
real(R8P),    parameter :: SDF_DEFAULT_CONE_DEG  = 120._R8P !< Shapira's default cone half-aperture is 60°, full angle 120°.

real(R8P),    parameter :: SDF_DEFAULT_SMOOTHING_LAMBDA     = 0.5_R8P  !< Per-iteration blend weight toward neighbour mean.
integer(I4P), parameter :: SDF_DEFAULT_SMOOTHING_ITERATIONS = 2_I4P    !< Shapira's default smoothing pass count.

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

   subroutine smooth_sdf(facet, sdf, lambda, iterations, status)
   !< Laplacian smoothing of a per-facet SDF over the dual graph (issue #18 §1.9 step 2).
   !<
   !< Update rule per pass: `sdf_new[f] = (1 - lambda) * sdf[f] + lambda * mean(sdf[neighbours])`.
   !< Adjacency is read from `facet%fcon_edge` (0 = boundary). Sentinel facets
   !< are passed through unchanged AND excluded from neighbour averages — a
   !< sentinel neighbour does not contaminate the mean. Facets whose only
   !< neighbours are all sentinels (or boundary) keep their pre-pass SDF.
   !<
   !< Two passes (the default) are sufficient for the canonical Shapira pipeline.
   !< With `lambda = 0.5` (default) and 2 passes, a feature spans ~3 facets in
   !< its support — the typical fillet width. Larger `lambda` or more passes
   !< washes out genuine SDF discontinuities; smaller values leave per-facet
   !< noise from the cone-ray Monte Carlo. The defaults are not knobs the user
   !< usually wants to tune.
   type(facet_object),     intent(in)              :: facet(:)    !< Facet array (only `fcon_edge` is read).
   real(R8P),              intent(inout)           :: sdf(:)      !< Per-facet SDF; mutated in place.
   real(R8P),              intent(in),    optional :: lambda      !< Blend weight in [0, 1] (default SDF_DEFAULT_SMOOTHING_LAMBDA).
   integer(I4P),           intent(in),    optional :: iterations  !< Pass count (default SDF_DEFAULT_SMOOTHING_ITERATIONS).
   integer(I4P),           intent(out),   optional :: status      !< Status code.
   real(R8P), allocatable                          :: sdf_new(:)  !< Double-buffer output of one pass.
   real(R8P)                                       :: lam, sum_nb !< Resolved lambda + per-facet neighbour sum.
   integer(I4P)                                    :: it, f, e, n_nb, n_iter, nf
   integer(I4P)                                    :: nb_id

   if (present(status)) status = SDF_STATUS_OK
   nf = size(sdf, kind=I4P)
   if (nf == 0_I4P) return
   if (size(facet, kind=I4P) /= nf) then
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif
   lam = SDF_DEFAULT_SMOOTHING_LAMBDA;     if (present(lambda))     lam = lambda
   n_iter = SDF_DEFAULT_SMOOTHING_ITERATIONS; if (present(iterations)) n_iter = iterations
   if (lam < 0._R8P .or. lam > 1._R8P .or. n_iter < 0_I4P) then
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif
   if (n_iter == 0_I4P) return  ! caller asked for no smoothing — no-op

   allocate(sdf_new(nf))
   do it = 1_I4P, n_iter
      do f = 1_I4P, nf
         if (sdf(f) == SDF_SENTINEL) then
            sdf_new(f) = SDF_SENTINEL
            cycle
         endif
         sum_nb = 0._R8P
         n_nb = 0_I4P
         do e = 1_I4P, 3_I4P
            nb_id = facet(f)%fcon_edge(e)
            if (nb_id <= 0_I4P)      cycle  ! boundary edge
            if (nb_id > nf)          cycle  ! defensive: shouldn't happen with valid connectivity
            if (sdf(nb_id) == SDF_SENTINEL) cycle
            sum_nb = sum_nb + sdf(nb_id)
            n_nb = n_nb + 1_I4P
         enddo
         if (n_nb == 0_I4P) then
            sdf_new(f) = sdf(f)  ! isolated valid facet: no smoothing source
         else
            sdf_new(f) = (1._R8P - lam) * sdf(f) + lam * (sum_nb / real(n_nb, R8P))
         endif
      enddo
      sdf = sdf_new
   enddo
   deallocate(sdf_new)
   endsubroutine smooth_sdf

   subroutine segment_sdf(facet, tree, bmin, bmax, facet_labels, sdf, &
                          num_clusters, smoothing_lambda, smoothing_iterations, &
                          num_rays, cone_angle_deg, status)
   !< Capstone driver: compute SDF, smooth it, fit a 1D GMM, return per-facet
   !< hard labels via posterior argmax (issue #18 §1.9 step 3).
   !<
   !< Output `facet_labels(1:nf)`:
   !<   - `SDF_LABEL_UNASSIGNED` (== 0) for facets with sentinel SDF (degenerate
   !<     shells where the cone-of-rays returned too few hits).
   !<   - `1..num_clusters` for valid facets, with the cluster index ordered by
   !<     ascending GMM mean (cluster 1 = thinnest features, cluster k = thickest).
   !<
   !< Optional output `sdf(1:nf)` exposes the smoothed per-facet SDF for users
   !< who want to inspect or post-process the underlying scalar field.
   !<
   !< If the user requests `num_clusters` larger than the data supports (e.g.
   !< `k=4` on a cube where only 1 SDF value exists), GMM may converge with
   !< several near-identical means; argmax labelling still yields a valid
   !< partition (just into fewer effective groups). This is intended behaviour
   !< — we don't try to merge clusters automatically because a small bimodal
   !< feature might genuinely warrant its own label that is dimensionally close
   !< to the bulk.
   type(facet_object),       intent(in)              :: facet(:)            !< Surface facets.
   type(aabb_tree_object),   intent(in)              :: tree                !< Built AABB tree.
   type(vector_R8P),         intent(in)              :: bmin                !< Bbox min.
   type(vector_R8P),         intent(in)              :: bmax                !< Bbox max.
   integer(I4P), allocatable, intent(out)            :: facet_labels(:)     !< Per-facet hard label ∈ [0, num_clusters].
   real(R8P), allocatable,   intent(out),   optional :: sdf(:)              !< Per-facet smoothed SDF (optional output).
   integer(I4P),             intent(in),    optional :: num_clusters        !< Default SDF_DEFAULT_NUM_CLUSTERS (4).
   real(R8P),                intent(in),    optional :: smoothing_lambda    !< Default SDF_DEFAULT_SMOOTHING_LAMBDA.
   integer(I4P),             intent(in),    optional :: smoothing_iterations!< Default SDF_DEFAULT_SMOOTHING_ITERATIONS.
   integer(I4P),             intent(in),    optional :: num_rays            !< Default SDF_DEFAULT_NUM_RAYS.
   real(R8P),                intent(in),    optional :: cone_angle_deg      !< Default SDF_DEFAULT_CONE_DEG.
   integer(I4P),             intent(out),   optional :: status              !< Status code.
   real(R8P), allocatable                            :: sdf_local(:)        !< Working SDF (raw → smoothed).
   real(R8P), allocatable                            :: sdf_valid(:)        !< Non-sentinel SDF values, packed.
   integer(I4P), allocatable                         :: valid_idx(:)        !< Index map: sdf_valid(i) → original facet id.
   real(R8P), allocatable                            :: gmm_mean(:), gmm_var(:), gmm_w(:)  !< GMM parameters.
   integer(I4P), allocatable                         :: order(:)            !< Sort permutation of GMM means.
   integer(I4P)                                      :: nf, n_valid, k, i, f, lab
   integer(I4P)                                      :: ierr

   if (present(status)) status = SDF_STATUS_OK

   nf = size(facet, kind=I4P)
   if (nf == 0_I4P) then
      allocate(facet_labels(0))
      if (present(sdf)) allocate(sdf(0))
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif
   k = SDF_DEFAULT_NUM_CLUSTERS;  if (present(num_clusters)) k = num_clusters
   if (k < 1_I4P) then
      allocate(facet_labels(nf))
      facet_labels = SDF_LABEL_UNASSIGNED
      if (present(sdf)) then
         allocate(sdf(nf))
         sdf = SDF_SENTINEL
      endif
      if (present(status)) status = SDF_STATUS_BAD_INPUT
      return
   endif

   call compute_sdf(facet=facet, tree=tree, bmin=bmin, bmax=bmax, &
                    sdf=sdf_local, num_rays=num_rays, cone_angle_deg=cone_angle_deg, &
                    status=ierr)
   if (ierr /= SDF_STATUS_OK) then
      allocate(facet_labels(nf))
      facet_labels = SDF_LABEL_UNASSIGNED
      if (present(sdf)) call move_alloc(from=sdf_local, to=sdf)
      if (present(status)) status = ierr
      return
   endif
   call smooth_sdf(facet=facet, sdf=sdf_local, lambda=smoothing_lambda, &
                   iterations=smoothing_iterations, status=ierr)
   if (ierr /= SDF_STATUS_OK) then
      allocate(facet_labels(nf))
      facet_labels = SDF_LABEL_UNASSIGNED
      if (present(sdf)) call move_alloc(from=sdf_local, to=sdf)
      if (present(status)) status = ierr
      return
   endif

   ! Pack the valid (non-sentinel) SDF values for GMM.
   n_valid = 0_I4P
   do f = 1_I4P, nf
      if (sdf_local(f) /= SDF_SENTINEL) n_valid = n_valid + 1_I4P
   enddo
   allocate(facet_labels(nf))
   facet_labels = SDF_LABEL_UNASSIGNED
   if (n_valid == 0_I4P) then
      if (present(sdf)) call move_alloc(from=sdf_local, to=sdf)
      return  ! all facets were sentinels — every label stays UNASSIGNED, status OK
   endif
   allocate(sdf_valid(n_valid))
   allocate(valid_idx(n_valid))
   i = 0_I4P
   do f = 1_I4P, nf
      if (sdf_local(f) /= SDF_SENTINEL) then
         i = i + 1_I4P
         sdf_valid(i) = sdf_local(f)
         valid_idx(i) = f
      endif
   enddo

   ! Cap k at n_valid (can't have more clusters than data points).
   k = min(k, n_valid)

   call gmm_fit_1d(x=sdf_valid, k=k, mean=gmm_mean, var=gmm_var, weight=gmm_w)

   ! Sort cluster indices by ascending mean → stable label ordering for tests.
   allocate(order(k))
   call argsort_ascending(gmm_mean, order)

   ! Assign each valid facet to argmax-posterior cluster.
   do i = 1_I4P, n_valid
      lab = gmm_argmax_posterior(x=sdf_valid(i), mean=gmm_mean, var=gmm_var, weight=gmm_w)
      ! Map raw cluster index → ascending-mean rank (so label 1 = thinnest).
      facet_labels(valid_idx(i)) = inv_perm(order, lab)
   enddo

   if (present(sdf)) call move_alloc(from=sdf_local, to=sdf)
   endsubroutine segment_sdf

   ! ---------- 1D Gaussian Mixture Model via EM with k-means++ init ----------

   subroutine gmm_fit_1d(x, k, mean, var, weight)
   !< Fit a 1D `k`-component Gaussian mixture to `x(:)` via EM, with k-means++
   !< initialization. Returns (mean, var, weight) for each component.
   !<
   !< Hard-coded safety: variance floored at `1e-12 * total_variance` (kills the
   !< empty-cluster collapse where σ² → 0 → log-likelihood → −∞), max 50 EM
   !< iterations (well over the typical 5-10 for 1D), convergence on log-
   !< likelihood Δ < 1e-9 * |log-likelihood|.
   real(R8P),              intent(in)  :: x(:)        !< Data points.
   integer(I4P),           intent(in)  :: k           !< Number of components.
   real(R8P), allocatable, intent(out) :: mean(:)     !< (k) component means.
   real(R8P), allocatable, intent(out) :: var(:)      !< (k) component variances.
   real(R8P), allocatable, intent(out) :: weight(:)   !< (k) component mixture weights, sum to 1.
   real(R8P), allocatable              :: resp(:,:)   !< Responsibilities r(n,k).
   real(R8P)                           :: total_var, var_floor
   real(R8P)                           :: ll, ll_prev
   integer(I4P)                        :: n, it, j

   integer(I4P), parameter :: MAX_EM_ITER = 50_I4P
   real(R8P),    parameter :: LL_REL_TOL  = 1.0e-9_R8P

   n = size(x, kind=I4P)
   allocate(mean(k), var(k), weight(k), resp(n, k))

   total_var = data_variance(x)
   var_floor = max(1.0e-12_R8P * total_var, 1.0e-30_R8P)  ! absolute floor for total_var = 0 case

   ! ---- k-means++ initialization for means; uniform weights; total-variance for vars.
   call kmeanspp_init(x=x, k=k, centers=mean)
   weight = 1._R8P / real(k, R8P)
   var = max(total_var, var_floor)

   ll_prev = -huge(1._R8P)
   do it = 1_I4P, MAX_EM_ITER
      ! E step: responsibilities r(n, k) ∝ weight(k) * N(x_n | mean(k), var(k)).
      call gmm_e_step(x=x, mean=mean, var=var, weight=weight, resp=resp, ll=ll)

      ! Convergence check (after at least one E step).
      if (it > 1_I4P) then
         if (abs(ll - ll_prev) <= LL_REL_TOL * max(1._R8P, abs(ll))) exit
      endif
      ll_prev = ll

      ! M step: update mean / var / weight from responsibilities.
      call gmm_m_step(x=x, resp=resp, mean=mean, var=var, weight=weight, var_floor=var_floor)
   enddo

   deallocate(resp)
   ! Cluster ordering is left as-fit; the caller (segment_sdf) sorts by mean.
   ! Suppress unused-loop-var warning in some compilers:
   j = it
   endsubroutine gmm_fit_1d

   subroutine kmeanspp_init(x, k, centers)
   !< Stateless k-means++ seeding for 1D data. Deterministic: first center is the
   !< median (stable, not noisy), subsequent centers picked greedily by maximum
   !< squared distance from the current center set. Skips duplicates so a
   !< constant-data input still produces k distinct centers (with vanishing
   !< variance — handled by the variance floor in `gmm_fit_1d`).
   real(R8P),    intent(in)  :: x(:)
   integer(I4P), intent(in)  :: k
   real(R8P),    intent(out) :: centers(:)
   real(R8P), allocatable    :: d2(:)         !< Squared distance to nearest existing center.
   integer(I4P)              :: n, i, c, picked
   real(R8P)                 :: best_d2

   n = size(x, kind=I4P)
   allocate(d2(n))
   centers(1) = median(x)
   do i = 1_I4P, n
      d2(i) = (x(i) - centers(1))**2
   enddo
   do c = 2_I4P, k
      ! Pick the point farthest from any existing center.
      picked = 1_I4P
      best_d2 = -1._R8P
      do i = 1_I4P, n
         if (d2(i) > best_d2) then
            best_d2 = d2(i)
            picked = i
         endif
      enddo
      centers(c) = x(picked)
      do i = 1_I4P, n
         d2(i) = min(d2(i), (x(i) - centers(c))**2)
      enddo
   enddo
   deallocate(d2)
   endsubroutine kmeanspp_init

   pure subroutine gmm_e_step(x, mean, var, weight, resp, ll)
   !< Compute responsibilities r(n,k) and total log-likelihood of `x` under the
   !< current GMM parameters. Numerically stable via log-sum-exp on log-
   !< responsibilities.
   real(R8P), intent(in)  :: x(:)
   real(R8P), intent(in)  :: mean(:), var(:), weight(:)
   real(R8P), intent(out) :: resp(:,:)
   real(R8P), intent(out) :: ll
   integer(I4P)           :: n, k, i, j
   real(R8P)              :: log_max, denom
   real(R8P), allocatable :: log_p(:)         !< log[ w_j * N(x_i | mean_j, var_j) ]
   real(R8P), parameter   :: LOG_2PI = 1.8378770664093454835606594728112_R8P

   n = size(x, kind=I4P)
   k = size(mean, kind=I4P)
   ll = 0._R8P
   allocate(log_p(k))
   do i = 1_I4P, n
      do j = 1_I4P, k
         log_p(j) = log(weight(j)) - 0.5_R8P * (LOG_2PI + log(var(j)) + (x(i) - mean(j))**2 / var(j))
      enddo
      log_max = maxval(log_p)
      denom = 0._R8P
      do j = 1_I4P, k
         denom = denom + exp(log_p(j) - log_max)
      enddo
      ll = ll + log_max + log(denom)
      do j = 1_I4P, k
         resp(i, j) = exp(log_p(j) - log_max) / denom
      enddo
   enddo
   deallocate(log_p)
   endsubroutine gmm_e_step

   pure subroutine gmm_m_step(x, resp, mean, var, weight, var_floor)
   !< Update GMM parameters from responsibilities. Variance floored at
   !< `var_floor` to prevent collapse on duplicate points.
   real(R8P), intent(in)    :: x(:)
   real(R8P), intent(in)    :: resp(:,:)
   real(R8P), intent(inout) :: mean(:), var(:), weight(:)
   real(R8P), intent(in)    :: var_floor
   real(R8P)                :: nk, mu_j, vs
   integer(I4P)             :: n, k, i, j

   n = size(x, kind=I4P)
   k = size(mean, kind=I4P)
   do j = 1_I4P, k
      nk = 0._R8P
      do i = 1_I4P, n
         nk = nk + resp(i, j)
      enddo
      if (nk <= 0._R8P) then
         ! Empty cluster: keep its mean, set tiny variance and weight.
         var(j) = var_floor
         weight(j) = 1.0e-12_R8P
         cycle
      endif
      mu_j = 0._R8P
      do i = 1_I4P, n
         mu_j = mu_j + resp(i, j) * x(i)
      enddo
      mu_j = mu_j / nk
      vs = 0._R8P
      do i = 1_I4P, n
         vs = vs + resp(i, j) * (x(i) - mu_j)**2
      enddo
      mean(j)   = mu_j
      var(j)    = max(vs / nk, var_floor)
      weight(j) = nk / real(n, R8P)
   enddo
   endsubroutine gmm_m_step

   pure function gmm_argmax_posterior(x, mean, var, weight) result(j_star)
   !< Return the cluster index j ∈ [1, k] with maximum posterior at scalar `x`.
   !< No need to compute the normalizer — comparing log-numerators is enough.
   real(R8P),    intent(in) :: x
   real(R8P),    intent(in) :: mean(:), var(:), weight(:)
   integer(I4P)             :: j_star
   integer(I4P)             :: j, k
   real(R8P)                :: lp, best_lp

   k = size(mean, kind=I4P)
   j_star = 1_I4P
   best_lp = log(weight(1)) - 0.5_R8P * (log(var(1)) + (x - mean(1))**2 / var(1))
   do j = 2_I4P, k
      lp = log(weight(j)) - 0.5_R8P * (log(var(j)) + (x - mean(j))**2 / var(j))
      if (lp > best_lp) then
         best_lp = lp
         j_star = j
      endif
   enddo
   endfunction gmm_argmax_posterior

   pure function data_variance(x) result(v)
   !< Sample variance (unbiased denom n; we just want a scale, not an estimator).
   real(R8P), intent(in) :: x(:)
   real(R8P)             :: v, mu
   integer(I4P)          :: n, i

   n = size(x, kind=I4P)
   if (n <= 1_I4P) then
      v = 0._R8P
      return
   endif
   mu = 0._R8P
   do i = 1_I4P, n
      mu = mu + x(i)
   enddo
   mu = mu / real(n, R8P)
   v = 0._R8P
   do i = 1_I4P, n
      v = v + (x(i) - mu)**2
   enddo
   v = v / real(n, R8P)
   endfunction data_variance

   pure subroutine argsort_ascending(x, order)
   !< Insertion sort returning the permutation `order` such that
   !< `x(order(1)) <= x(order(2)) <= ...`. n is small (= num_clusters).
   real(R8P),    intent(in)  :: x(:)
   integer(I4P), intent(out) :: order(:)
   integer(I4P)              :: i, j, k, key

   k = size(x, kind=I4P)
   do i = 1_I4P, k
      order(i) = i
   enddo
   do i = 2_I4P, k
      key = order(i)
      j = i - 1_I4P
      do while (j >= 1_I4P)
         if (x(order(j)) <= x(key)) exit
         order(j + 1_I4P) = order(j)
         j = j - 1_I4P
      enddo
      order(j + 1_I4P) = key
   enddo
   endsubroutine argsort_ascending

   pure function inv_perm(order, j) result(rank)
   !< Return rank of cluster index `j` in the ascending-mean order (1-based).
   !< Equivalent to: find i such that order(i) == j, return i.
   integer(I4P), intent(in) :: order(:)
   integer(I4P), intent(in) :: j
   integer(I4P)             :: rank, i

   rank = 0_I4P
   do i = 1_I4P, size(order, kind=I4P)
      if (order(i) == j) then
         rank = i
         return
      endif
   enddo
   endfunction inv_perm

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
