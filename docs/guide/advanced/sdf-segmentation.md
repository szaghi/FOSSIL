---
title: SDF segmentation
---

# SDF segmentation

## What it does

Synthesizes per-facet **semantic labels** from geometry alone,
partitioning the surface into regions of similar local thickness. CAD-
exported STL almost never carries patch / surface-group information —
every facet is in one anonymous bag. This TBP gives you a per-facet
integer in `[0, k]` that downstream code (BC tagging, CFD setup) can
use as a substitute for the `.cgns` patch hierarchy that's missing.

Three-stage pipeline (Shapira, Shamir & Cohen-Or 2008):

1. **Per-facet Shape Diameter Function** — for each facet, cast a
   cone of rays inward (along `-normal`), measure first-hit distance
   on the opposite surface, take the median. The **shape diameter** is
   approximately the local thickness of the solid at that facet.
2. **Dual-graph Laplacian smoothing** — average each facet's SDF with
   its 3 neighbours via `facet%fcon_edge`. Two passes by default
   (Shapira's value).
3. **1D Gaussian-mixture clustering** — fit a `k`-component GMM to
   the smoothed SDF values via EM with k-means++ initialization.
   Per-facet hard label = argmax posterior, with clusters re-ordered
   ascending by mean so label 1 = thinnest features, label `k` =
   thickest.

## Pipeline detail

For step 1: each facet emits `num_rays` rays (default 30) in a
cone of half-aperture `cone_angle_deg/2` (default 60° = 120° full
aperture) around the inward normal, distributed via Fibonacci-spiral
sampling (deterministic, low-discrepancy, no clumping at pole or
rim). Rays start at `centroid - ε·normal` (offset trick) to skip the
self-hit. If fewer than 50 % of rays hit the opposite surface, the
facet is flagged with sentinel SDF (-1) and excluded from clustering.

For step 2: `sdf_new[f] = (1-λ)·sdf[f] + λ·mean(sdf[neighbours])`.
Sentinels propagate (a sentinel facet stays sentinel) AND are
excluded from neighbour averages (sentinel neighbours don't
contaminate the mean).

For step 3: variance-floored EM (`var ≥ 1e-12 · total_var` prevents
empty-cluster collapse), log-sum-exp stable E step, 50-iteration cap
with log-likelihood Δ < 1e-9 convergence. K-means++ init is
**deterministic** (first center = median, subsequent = max-d²
greedy) for reproducibility.

## API

```fortran
call surface%segment_sdf(facet_labels, sdf, num_clusters, smoothing_lambda, &
                         smoothing_iterations, num_rays, cone_angle_deg, status)
```

| Argument                | Intent       | Type                       | Default       | Notes                                                                |
| ----------------------- | ------------ | -------------------------- | ------------- | -------------------------------------------------------------------- |
| `facet_labels`          | `out`        | `integer(I4P), allocatable`| —             | Per-facet labels in `[0, num_clusters]`.                             |
| `sdf`                   | `out`, opt.  | `real(R8P), allocatable`   | —             | Smoothed per-facet SDF for inspection.                               |
| `num_clusters`          | `in`, opt.   | `integer(I4P)`             | `4`           | GMM component count. Mechanical-CAD-friendly; CGAL uses 5.           |
| `smoothing_lambda`      | `in`, opt.   | `real(R8P)`                | `0.5`         | Laplacian blend weight ∈ `[0, 1]`.                                   |
| `smoothing_iterations`  | `in`, opt.   | `integer(I4P)`             | `2`           | Smoothing pass count. Shapira's default.                             |
| `num_rays`              | `in`, opt.   | `integer(I4P)`             | `30`          | Cone rays per facet. Shapira's default.                              |
| `cone_angle_deg`        | `in`, opt.   | `real(R8P)`                | `120.0`       | Full cone aperture, degrees.                                         |
| `status`                | `out`, opt.  | `integer(I4P)`             | —             | `SDF_STATUS_*`.                                                      |

`facet_labels(i) == 0` (`SDF_LABEL_UNASSIGNED`) means the facet had
sentinel SDF (degenerate shell, isolated patch). Downstream code
should treat unassigned facets as their own implicit "unknown" group.

## Example

```fortran
program ex_segment_sdf
use fossil
use penf, only : I4P, R8P
implicit none

type(surface_stl_object)  :: bunny
integer(I4P), allocatable :: labels(:)
real(R8P),    allocatable :: sdf(:)
integer(I4P)              :: status, k, lo, hi, n_unassigned

call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
call bunny%segment_sdf(facet_labels=labels, sdf=sdf, num_clusters=4_I4P, status=status)

lo = minval(labels)
hi = maxval(labels)
n_unassigned = count(labels == 0_I4P)
print '(A,I0,A,I0,A,I0,A,I0,A,I0)', 'bunny labels: nf=', size(labels), &
      '  range=[', lo, ',', hi, ']  unassigned=', n_unassigned
do k = 1_I4P, 4_I4P
   print '(A,I0,A,I0)', '  cluster ', k, ': nf = ', count(labels == k)
enddo
endprogram ex_segment_sdf
```

On the bunny fixture: 69 451 facets, labels in `[1, 4]`, zero
unassigned (closed manifold, ray hit-rate ~100 %). The four cluster
counts roughly correspond to "thinnest features" (ears, tail tips —
label 1) through "thickest body region" (label 4).

## Known limitations

- **GMM is a mode-finder, not a mode-merger.** With `k` larger than
  the data supports, GMM will split a tightly-clustered distribution
  into two near-identical components — and the argmax labelling
  treats those as distinct labels. On a uniform input (sphere with
  `k = 2`), the label distribution is **not** trivially "all label
  1"; it'll split ~62 / 38 between two labels of essentially the
  same SDF cluster. **This is intended GMM behaviour, not a bug.**
- **Cube + `k = 2` is a misleading test fixture.** The cube's 12-facet
  diagonal-split tessellation produces a *genuinely bimodal* SDF
  (some facets see longer cone-ray paths to the opposite face), so
  GMM correctly splits — even though intuitively "the cube has only
  one feature scale." Use sphere with `k = 1` to test the trivial
  cluster path.
- **No alpha-expansion graph cut (Shapira §5)**. The full Shapira
  pipeline runs a graph-cut post-pass over the dual graph to enforce
  label coherence across low-dihedral-angle edges (segmentation
  respects sharp edges). FOSSIL's MVP omits this — labels can flip
  across edges that smoothing didn't already align. **Deferred to
  §1.9b follow-up.** For typical AMR-CFD STL the cone-ray + smoothing
  signal is dominant and the labels are usable as-is.
- **Sentinel SDF facets are excluded from clustering.** They get
  label 0 (`SDF_LABEL_UNASSIGNED`). On clean closed meshes this is
  rare (< 0.1 % of facets); on dirty meshes with isolated patches it
  can be substantial.
- **Reproducibility**: k-means++ init is deterministic but **not
  parallel-safe**. Don't multi-thread `segment_sdf`.

## See also

- [ray queries](/guide/advanced/ray-queries) — the `intersect_ray_first`
  TBP that the cone-of-rays SDF query uses internally.
- [`get_facets_number`](/guide/api/surface-stl-object) — to size your
  `labels` allocation if pre-allocating.
- [Constants](/guide/api/constants) — `SDF_STATUS_*`,
  `SDF_LABEL_UNASSIGNED`, `SDF_SENTINEL`.

## References

- Shapira, Shamir & Cohen-Or, *Consistent Mesh Partitioning and
  Skeletonisation Using the Shape Diameter Function*, The Visual
  Computer 24(4), 2008. The reference algorithm.
- CGAL `Surface_mesh_segmentation` package documentation. The
  reference implementation; FOSSIL matches its API shape but omits
  the alpha-expansion post-pass.
- Boykov & Kolmogorov, *An Experimental Comparison of Min-Cut/Max-
  Flow Algorithms for Energy Minimization in Vision*, PAMI 2004.
  The graph-cut algorithm that §1.9b would adopt.
