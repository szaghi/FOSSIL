---
title: Isotropic remeshing (§1.7)
---

# §1.7 — Isotropic remeshing

## What it does

Equalizes **edge length** across a triangulated surface while preserving
its geometry. CAD-exported STL routinely mixes 0.01 mm and 10 mm
triangles on the same body — perfectly fine for visualization, **fatal**
for AABB-tree query balance and for cut-cell quality near refinement
bands. Isotropic remeshing rewrites the tessellation so every edge is
within ~30% of a target length `L`, producing a mesh whose every facet
has roughly the same shape and size.

Complementary to [decimation (§1.3)](https://github.com/szaghi/FOSSIL/issues/18):
decimation reduces facet *count* while preserving the silhouette;
isotropic remeshing equalizes edge *length* without preferring either
direction. Use decimation to make a million-facet model tractable; use
isotropic remeshing to make any model's tessellation regular.

## Pipeline

Each outer iteration runs four passes (Botsch & Kobbelt 2004):

1. **Split** edges longer than `4L/3`. New vertex at the midpoint;
   incident triangles re-triangulated.
2. **Collapse** edges shorter than `4L/5`. Reuses §1.3's safety
   predicates: reject the collapse if it would flip a triangle's
   normal, create a duplicate vertex, or open a non-manifold edge.
3. **Flip** interior edges to drive vertex valence toward 6 (the
   regular-triangulation ideal). Rejects flips that would create a
   duplicate edge or a degenerate triangle.
4. **Tangential relaxation**: move each vertex toward the area-weighted
   centroid of its one-ring, then **project** back onto the original
   surface via the AABB tree's closest-point query. Without projection
   the mesh would shrink toward its volumetric centroid.

Three to five outer iterations is typical. Connectivity is rebuilt from
scratch between passes — chosen after the first attempt's packed-array
incremental design produced uncatchable heap corruption; the
"rebuild every pass" approach trades speed for debuggability.

**Feature preservation** is optional and off by default. When enabled,
edges with dihedral angle > 30° are flagged as **sharp**, and their
endpoint vertices are **locked**: collapse skips them, flip skips any
edge with a locked vertex, relaxation skips them entirely. This is what
keeps the 12 edges of a cube sharp through remeshing — without it, the
cube rounds out into a ball over a few iterations.

## API

```fortran
call surface%isotropic_remesh(target_length, iterations, preserve_features, status)
   real(R8P),    intent(in), optional :: target_length     ! ≤ 0 → median input edge
   integer(I4P), intent(in), optional :: iterations        ! default 5
   logical,      intent(in), optional :: preserve_features ! default .false.
   integer(I4P), intent(out), optional :: status
```

- **`target_length`**: desired uniform edge length. Pass `≤ 0` (or
  omit) to use the input's *median* edge length — produces a near-
  identity remesh that just smooths the distribution. Pass a smaller
  value to densify, larger to coarsen.
- **`iterations`**: outer-loop count. Default `5`. Convergence is
  typically visible after 3 passes; further iterations smooth the edge-
  length distribution further but also accumulate the inward Laplacian
  bias (see Limitations).
- **`preserve_features`**: lock vertices on edges with dihedral > 30°.
  **Essential for mechanical CAD inputs**; should usually be `.true.`
  on anything that isn't an organic shape.

Mutates `self` in place.

## Example

Remesh a cube with sharp-feature preservation, then a sphere without
(the sphere has no sharp features so the flag is moot):

```fortran
program ex_isotropic_remesh
use fossil
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object) :: cube, sphere
integer(I4P)             :: status

! Cube — preserve the 12 sharp edges through remeshing.
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
print '(A,I0,A,F8.4)', 'cube before: n_facets = ', cube%get_facets_number(), &
                       '   smallest edge = ', cube%smallest_edge_len()
call cube%isotropic_remesh(target_length=0.2_R8P, iterations=3_I4P, &
                            preserve_features=.true., status=status)
print '(A,I0,A,F8.4)', 'cube after : n_facets = ', cube%get_facets_number(), &
                       '   smallest edge = ', cube%smallest_edge_len()

! Sphere — no sharp features, smooth remesh OK.
call sphere%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
call sphere%isotropic_remesh(iterations=3_I4P, status=status)
print '(A,I0)', 'sphere remeshed: n_facets = ', sphere%get_facets_number()
endprogram ex_isotropic_remesh
```

## Bunny remesh — visual reference

Input bunny (Stanford fixture, 69451 facets, mixed-scale tessellation
typical of scanned data):

![Stanford bunny — input](/pictures/advanced/remesh-bunny-input.png)

After 3 iterations of remeshing at `target_length = 0.05`:

![Stanford bunny — after 3 iterations](/pictures/advanced/remesh-bunny-output.png)

The silhouette is preserved; the per-facet edge-length distribution
tightens substantially (median unchanged, stddev drops by >5×). Facet
count is roughly preserved (69451 → 69259) — remeshing isn't
decimation, it just redistributes triangles.

## Known limitations

- **Inward Laplacian bias on convex regions.** The area-weighted
  Laplacian centroid of a vertex's one-ring sits *inside* the surface
  by ~`r · (1 - cos θ)` per relax step on a region of radius `r` with
  half-aperture `θ`. The projection step after relaxation pulls the
  vertex back toward the input surface but doesn't fully compensate.
  **Net effect**: a sphere loses ~5% volume per outer iteration even
  with projection enabled. For convex inputs where preservation
  matters, run fewer iterations (2-3) or use the `target_length`
  parameter conservatively. Tighter preservation requires
  normal-direction projection or scaled-step heuristics; deferred to a
  future refinement.
- **Feature preservation uses a hard 30° dihedral cutoff.** Edges below
  the threshold get smoothed even if they were intended as features
  (subtle creases on a fillet). The cutoff isn't user-tunable in the
  current API; if your input has sub-30° features that must survive,
  pre-sharpen them by inserting a tiny fillet before remeshing.
- **`target_length` is a guideline, not a contract.** The split (`>
  4L/3`) and collapse (`< 4L/5`) thresholds let edges in `[0.8L,
  1.33L]` survive untouched. Ratio fixed in the literature; not a
  parameter.
- **Connectivity is rebuilt every pass**, not incrementally maintained.
  Per-pass cost is `O(N log N)` instead of `O(δ_changed)`. Acceptable
  for typical inputs; a 70k-facet remesh runs in ~1 s. For million-
  facet meshes consider running fewer iterations.
- **Volume drift compounds across iterations.** A sphere remeshed for
  10 iterations loses ~40% volume — verifiable via
  `surface%get_volume()` before and after. The current MVP-grade
  implementation should not be used for volume-critical pipelines
  without an explicit post-remesh volume-restoration pass.

## See also

- [`decimate`](/guide/api/surface-stl-object) (§1.3 — facet count
  reduction; complementary primitive).
- [`statistics`](/guide/api/surface-stl-object#statistics) — print
  edge-length distribution before and after to verify the remesh
  worked.
- The [§1.7 implementation issue](https://github.com/szaghi/FOSSIL/issues/18)
  has the full multi-commit history and design notes (rectangular vs.
  packed connectivity, the heap-corruption pivot, the relax-bias
  diagnosis).

## References

- Botsch & Kobbelt, *A Remeshing Approach to Multiresolution Modeling*,
  Symposium on Geometry Processing 2004. The four-pass pipeline.
- Botsch et al., *Polygon Mesh Processing*, AK Peters 2010, Chapter 6.
  Textbook treatment with implementation guidance.
- Alliez, Cohen-Steiner, Devillers, Lévy, & Desbrun, *Anisotropic
  Polygonal Remeshing*, SIGGRAPH 2003. The anisotropic
  generalization (out of FOSSIL scope; for surfaces with directional
  feature preservation requirements).
