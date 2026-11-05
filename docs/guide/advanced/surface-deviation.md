---
title: Surface deviation
---

# Surface deviation (Hausdorff distance)

## What it does

Measures how far one surface lies from another (issue #60, roadmap #58):

- the **one-sided deviation** $h(A, B) = \sup_{x \in A} d(x, B)$: how far the farthest point of A is from B;
- the **Hausdorff distance** $H(A, B) = \max\big(h(A, B), h(B, A)\big)$;
- the **mean** and **rms** of $d(x, B)$ over A, and $d$ at every vertex of A.

The two one-sided deviations differ in general: a small body inside a large one is close to it ($h$ = the gap), while
the large body's far corners are far from the small one.

Remeshing, decimation, smoothing, hole filling, alpha wrapping and booleans all move a surface. The deviation measures
by how much, so a test or a user can bound it instead of trusting it.

## Method

$d(x, B)$ comes from B's bounding-volume hierarchy (`distance_many`): exact closest-point distances, bit-identical
across tree kinds and thread counts.

**Lower bound.** A is sampled deterministically:

- every vertex;
- every edge midpoint;
- the centroids of the $n_f^2$ sub-triangles of a regular subdivision of each facet f, with
  $n_f = \lceil \sqrt{\text{area}_f \cdot \text{samples\_per\_area}} \rceil$. The default density is 4 samples per
  mean facet area, about 4 per facet.

The largest sampled distance is a **lower bound** of $h(A, B)$: the true maximum can only be larger.

**Mean and rms.** The sub-triangle centroids, each weighted by $\text{area}_f / n_f^2$, are a quadrature over A: the
centroid rule on each sub-triangle, exact for a function linear over it. The mean and rms of $d$ are therefore
**estimates**, exact when $d$ is constant (two parallel planes) and converging as the density grows.

## API

```fortran
call surface%deviation(other, lower, mean, rms, per_vertex, samples_per_area, samples, status)
   class(surface_stl_object), intent(in)            :: self             ! A
   type(surface_stl_object),  intent(in)            :: other            ! B
   real(R8P),                 intent(out), optional :: lower            ! lower bound of h(A, B)
   real(R8P),                 intent(out), optional :: mean, rms        ! estimates of d over A
   real(R8P), allocatable,    intent(out), optional :: per_vertex(:)    ! d at the vertex-pool vertices of A
   real(R8P),                 intent(in),  optional :: samples_per_area ! default 4 / mean facet area of A
   integer(I4P),              intent(out), optional :: samples          ! points evaluated
   integer(I4P),              intent(out), optional :: status           ! DEV_STATUS_OK, DEV_STATUS_BAD_INPUT

call surface%hausdorff(other, lower, samples_per_area, status)          ! H(A, B): both one-sided runs
```

- `status = DEV_STATUS_BAD_INPUT` when either surface is empty.
- Neither surface needs to be closed: the deviation is defined for any triangle soup.
- The work is one batched distance query per call. The result does not depend on the thread count: the queries are
  independent and the maximum and sums run in a fixed order.

## Example: how far did the repair move the dragon?

<<< @/examples/snippets/ck_deviation.f90

<<< @/examples/output/ck_deviation.ansi{ansi}

- **Repaired against original:** 0.076. The 54 hole patches span the holes, and their middles are that far from the
  original surface.
- **Original against repaired:** 3e-3, at the facets `sanitize` removed (degenerate and duplicate ones).
- The means are tiny: almost every point of either surface lies on the other.

## Validation

`src/tests/fossil_test_deviation.f90`:

| Case | Checked | Measured |
|---|---|---|
| `cube-02` against itself | 0 at the vertices, round-off inside | 2.8e-17 |
| `cube-02` against itself translated by 0.3 along x | h = 0.3 both ways (attained at vertices) | exact |
| Marching-cubes unit sphere against itself scaled by 1.05 | 0.98 · 0.05 ≤ h ≤ 0.05 | 0.04995 |
| The cube [0.4, 0.6]³ inside the unit cube | h(small, unit) = 0.4, h(unit, small) = 0.4√3, H = 0.4√3 | exact |
| Two parallel unit squares 0.3 apart | lower, mean, rms = 0.3 | exact |
| 1 and 4 OpenMP threads (`tests-gnu-openmp`) | lower, mean, per-vertex bit for bit | identical |
| An empty surface | `DEV_STATUS_BAD_INPUT` | |

## See also

- [Mass properties](/guide/advanced/mass-properties), [decimation](/guide/advanced/decimate),
  [isotropic remeshing](/guide/advanced/isotropic-remesh): operations whose result the deviation bounds.
- [`surface_stl_object`](/guide/api/surface-stl-object): `distance_many`, the query behind it.

## References

- M. Tang, M. Lee, Y. J. Kim, *Interactive Hausdorff distance computation for general polygonal models*, ACM SIGGRAPH
  2009.
- P. Cignoni, C. Rocchini, R. Scopigno, *Metro: measuring error on simplified surfaces*, Computer Graphics Forum 17(2),
  1998.
