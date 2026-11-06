---
title: Thickness and gaps
---

# Thickness and gaps

## What it does

Measures, for every facet of a closed surface (issue #62, roadmap #58):

- the **thickness**: how thick the body is behind the facet;
- the **gap**: how far the next wall is in front of the facet, which is the clearance between two parts of a body or
  between two bodies merged into one surface. It is `huge` when nothing is in front.

On real geometry, immersed-boundary methods fail where a wall or a gap is thinner than the cell. Such a wall gives
the wrong sign, cells with no fluid, and image points that land across the wall. The
[ambiguous foot points](/guide/advanced/ib-queries) of a query find that afterwards, point by point. Thickness and gap
find it before the grid exists, for the whole surface.

![a plate tapering from 0.05 to 0.45 along its length, coloured by its ray thickness from teal (thin) to ember (thick)](../../examples/images/ref-thickness-ray.svg){.fossil-figure}

![a flat plate under a tilted one, the upper face of the flat plate coloured by the gap to the tilted plate, from teal (0.15) to ember (0.35); faces with nothing in front are plain](../../examples/images/ref-thickness-gap.svg){.fossil-figure}

## Ray method

The ray method is the *shape diameter function* of Shapira, Shamir and Cohen-Or (2008), measured rather than used
for segmentation.

- **Rays.** From the facet centroid, `rays` rays (default 31) inside a cone of full aperture `cone_angle` (default
  60°), around $-n$ for the thickness and $+n$ for the gap. The directions are fixed: ray 1 is the axis, and ray $k$
  has the polar cosine
  $\cos\theta_k = 1 - (1 - \cos\tfrac{\alpha}{2})\,\tfrac{k-1}{m-1}$ (equal-area steps up to the rim) and the
  azimuth $(k-1)$ golden angles. No random numbers, so the same facet gives the same value on every run.
- **The source facet is excluded**, not stepped over: `intersect_ray_first(exclude=f)` skips it, so no origin offset
  enters the lengths.
- **Thickness rays** count when their first hit **leaves** the solid ($n_\text{hit}\cdot d > 0$). A miss, or a hit
  through a facet facing the wrong way, does not count.
- **Gap rays** count when their first hit **enters** the solid ($n_\text{hit}\cdot d < 0$), and when they miss, with
  length $\infty$.
- **Value.** The lower median of the counted lengths. A facet with fewer than half of its rays counted gets `huge`
  and is counted in `invalid`: an open or misoriented surface is reported, never silently averaged. A gap whose
  median is a miss is `huge` but valid: nothing is in front.
- **Vertex fields** (`vertex_thickness`, `vertex_gap`) are area-weighted averages of the finite values of the
  incident facets.
- **Parallel.** OpenMP over the facets. Each facet is computed on its own, so the result does not depend on the
  thread count.

### Bias of the cone

A cone measures oblique lengths, and the median keeps that bias. On a plate and a sphere it is known in closed form:

| Body | Ray thickness | Default (60°, 31 rays) | Shapira's 120° |
|---|---|---|---|
| plate of thickness $t$ | $t \cdot \mathrm{median}(1/\cos\theta_k)$ | $1.072\,t$ | $1.333\,t$ |
| sphere of radius $R$ | $2R \cdot \mathrm{median}(\cos\theta_k)$ | $0.933 \cdot 2R$ | $0.750 \cdot 2R$ |
| any, `rays=1` | the axial chord | $t$, $2R$ | |

The default cone is narrower than Shapira's 120°, which serves segmentation. For a measurement it would overestimate
a plate by a third. A cone still matters: one axial ray can graze an edge, while the median of 31 rays does not
depend on any single one.

[`segment_sdf`](/guide/advanced/sdf-segmentation) keeps its own shape diameter (120°, an origin offset, no
orientation test), unchanged: it is a segmentation feature, not a measurement.

## API

```fortran
call surface%thickness(thickness, gap, method, cone_angle, rays, vertex_thickness, vertex_gap, invalid, status)
!  thickness(nf), gap(nf)         per facet, huge where undefined or nothing is in front
!  method                         THICKNESS_RAY (default)
!  cone_angle, rays               the cone: full aperture [degrees] (default 60), rays per facet (default 31)
!  vertex_thickness, vertex_gap   per vertex of the vertex pool
!  invalid                        facets whose rays gave no value
```

Every argument is optional, and at least one output is required. The constants `THICKNESS_RAY`, `THICKNESS_MEDIAL`,
`THICKNESS_DEFAULT_CONE` and `THICKNESS_DEFAULT_RAYS` come from `use fossil`.

## Validation

`src/tests/fossil_test_thickness_ray.f90`:

| Case | Result |
|---|---|
| a plate 10 × 10 × 0.1, facets of the large faces | $t/\cos\theta_{16}$ to 1.4e-16; with `rays=1`, $t$ to 2.8e-16; gap `huge` everywhere |
| a faceted unit sphere (1280 facets) | $2\cos\theta_{16}$ to 6.1e-3, axial $2$ to 4.5e-3: the faceting |
| two unit boxes 0.2 apart | gap $0.2/\cos\theta_{16}$ on the facing faces to 2.8e-16, `huge` on the others |
| a box without one face | the facets whose rays escape are `invalid`, thickness `huge` there |
| `exclude` | equal, facet and $t$ bit for bit, to the first other hit of `intersect_ray_all` |
| 1 and 4 OpenMP threads | bit-identical |

## Limits

- **An estimate.** Ray thickness carries the cone bias above.
- **Closed, outward-oriented input.** Rays escape from open surfaces and hit misoriented facets from the wrong side.
  Both are counted in `invalid`.

## References

- L. Shapira, A. Shamir, D. Cohen-Or, *Consistent mesh partitioning and skeletonisation using the shape diameter
  function*, The Visual Computer 24 (2008).

## See also

- [Ray queries](/guide/advanced/ray-queries): the first-hit traversal the rays use.
- [Immersed-boundary queries](/guide/advanced/ib-queries): ambiguous foot points, the per-query symptom of a thin wall.
