---
title: Self-intersection
---

# Self-intersection detection and resolution

## What it does

Two TBPs that bracket the self-intersection problem:

- **`surface%find_self_intersections(pairs, status, include_touching)`** —
  *detection*. Returns the list of every facet pair `(a, b)` whose
  triangles cross or overlap in a common plane, decided exactly, with the
  kind of contact and the 3D segment endpoints of the intersection.
- **`surface%resolve_self_intersections(status)`** — *resolution*.
  Cleans the surface up by collapsing every self-crossing through a
  self-boolean union (`boolean(self, self, BOOL_UNION)`).

Real-world STL exports **routinely** have self-intersections. They
slip through `is_watertight()` (the topological-edge test is satisfied
even when triangles geometrically cross), so any downstream operation
that assumes a clean solid — booleans, signed distance, voxelization
for cut-cell — silently produces wrong answers. Detection is the first
step; resolution is what you reach for when sanitize alone can't fix
the input.

## Pipeline (detection)

1. **Tree-vs-tree broad phase.** Walk the AABB tree against itself,
   collecting every leaf-pair whose AABBs overlap. Reuses the existing
   `enumerate_children` dispatcher so SAH BVH and octree both share
   the traversal code.
2. **Pool-aware adjacency filter.** The shared vertices of the pair are
   counted by vertex-pool ID (exact coordinate equality when the pool is
   not in use). Two facets sharing an **edge** meet along it, as expected;
   they are a defect only when exactly coplanar and folded over it (both
   free vertices on the same side of the edge, exact `orient2d`):
   `SELFINT_EDGE_FOLD` (issue #69).
3. **Narrow phase.** No shared vertex (issue #69): the exact `pair_cut`
   of the boolean (`fossil_tritri_exact`, exact `orient3d` signs and
   implicit section points). A common segment of positive length is a
   crossing (`SELFINT_CROSSING`), reported with its endpoints rounded to
   doubles. A single common point is a touch (`SELFINT_TOUCHING`),
   reported only with `include_touching=.true.`. Two exactly coplanar
   triangles are compared in the projection that keeps them non-degenerate
   (exact `orient2d`, separating-axis test on the six edge lines):
   overlapping interiors are a defect (`SELFINT_COPLANAR`, `p` and `q`
   the two centroids), boundaries meeting without area a touch. Before
   #69 this was the floating Möller test, which never reported a coplanar
   overlap. One shared **vertex** (issue #46): the exact `pair_classify`
   (`orient3d`); a common segment of positive length is a crossing, the
   vertex alone is a touch and is not reported. Pairs whose other
   vertices lie strictly on one side of the other plane can only touch
   and are skipped. Exactly coplanar pairs sharing a vertex (issue #69;
   skipped to 1e-12 rad before, so crossings tilted by less went unseen)
   are decided in the projection: overlapping interiors are reported
   (`SELFINT_COPLANAR`), a common edge segment from the vertex is a touch
   (`include_touching`).

4. **Non-manifold edges** (issue #69). Sheets crossing along a shared
   edge (e.g. two cubes merged so that one passes through the other along a
   common edge) meet only on the edge, so no pair test sees them. A pass
   over the edges with 4 or more facets orders the facets around each
   edge by the angle of their half-planes, exactly (`orient3d` on input
   vertices), and notes whether each traverses the edge forwards or
   backwards. Embedded oriented sheets alternate around the edge (the
   winding number of the wedges stays in {w, w + 1}); two adjacent facets
   with the same direction bound a wedge wound twice, where the sheets
   cross. Each such pair is reported as `SELFINT_SHEETS`, with the edge as
   `p`–`q`. Fans with unequal numbers of forward and backward facets (an
   orientation or boundary defect, not a crossing) and fans with two
   facets in the same half-plane (folds, reported pairwise) are skipped.
   The criterion assumes consistently oriented sheets: on a misoriented
   input a reported pair can be an orientation defect rather than a
   crossing (`sanitize` orients the normals first).

Not reported: duplicate facets (three shared vertices), which
`remove_duplicate_facets` handles.

## Pipeline (resolution)

`resolve_self_intersections` is a **one-line wrapper** over the
boolean engine:

```fortran
call self%boolean(other=self, op=BOOL_UNION, status=status)
```

The §1.1 arrangement code treats `self vs self` symmetrically: every
self-crossing becomes an arrangement segment, every input triangle
gets retriangulated, the per-sub-triangle winding-number tagging
is meant to keep only the outer-manifold sub-triangles. It does not yet:
see [Known limitations](#known-limitations).

## API

### Detection

```fortran
call surface%find_self_intersections(pairs=pairs, status=status, include_touching=touching)
   type(intersection_pair_t), allocatable, intent(out) :: pairs(:)
   integer(I4P),                           intent(out), optional :: status
   logical,                                intent(in),  optional :: include_touching ! default .false.
```

`pairs` is always allocated (zero-length if no intersections found).
Each record carries:

- `pairs(i)%a`, `pairs(i)%b` — facet IDs with `a < b`.
- `pairs(i)%kind` — `SELFINT_CROSSING` (a common segment of positive
  length), `SELFINT_COPLANAR` (coplanar, overlapping interiors),
  `SELFINT_EDGE_FOLD` (sharing an edge, coplanar, folded over it),
  `SELFINT_SHEETS` (adjacent around a non-manifold edge with the same
  direction: the sheets cross along the edge) or, with
  `include_touching`, `SELFINT_TOUCHING` (one common point, or coplanar
  boundaries meeting without area).
- `pairs(i)%p`, `pairs(i)%q` — segment endpoints in world coordinates
  (the common point twice for a touch, the two centroids for a coplanar
  overlap).

### Resolution

```fortran
call surface%resolve_self_intersections(status=status)
   integer(I4P), intent(out), optional :: status  ! BOOL_STATUS_*
```

Mutates `self` in place. It does not currently remove the crossings:
see [Known limitations](#known-limitations).

## Example

Two unit cubes, the second moved by (0.3, 0.4, 0.45), merged into one
surface by `merge_solids`: each shell is closed, so the surface is
watertight, yet the shells pass through each other. Detection finds the
crossings, and removing them goes through the union of the two bodies
as separate surfaces (`report` prints the facets, closure, shells,
volume and remaining pairs; `split_shells` copies the facets of each
shell, by `facet_shell`, into a surface of its own):

<<< @/examples/snippets/ref_self_intersection-merge.f90

<<< @/examples/snippets/ref_self_intersection-detect.f90

<<< @/examples/snippets/ref_self_intersection-split.f90

<<< @/examples/output/ref_self_intersection.ansi{ansi}

- **Watertight is not enough.** The merged surface passes
  `is_watertight()` with two shells and volume 2, the overlap counted
  twice.
- **Detection.** 32 pairs, all `SELFINT_CROSSING`: no facets of the two
  cubes are coplanar at this offset. Their segments join into one closed
  curve, the seam where the two boundaries cross (figure below).
- **Removal.** The [boolean](/guide/advanced/booleans) union of the two
  shells is watertight, one shell, of volume 1.7690 = 2 − 0.7 · 0.6 ·
  0.55, with no self-intersection left. It needs the bodies as separate,
  individually clean surfaces; for a self-intersecting single shell see
  [Known limitations](#known-limitations).

![two unit cubes merged into one surface, drawn translucent, the second moved by (0.3, 0.4, 0.45); the 32 crossing segments find_self_intersections reports form one closed ember curve where the two boundaries cross](../../examples/images/ref-self-intersection-cubes.svg){.fossil-figure}

On the dragon fixture (6588 facets) the same call reports 1191 pairs:
1167 crossings and, since issue #69, 24 pairs of sheets crossing at its
non-manifold edges. For a clean mesh (`cube-01.stl`, `bunny.stl`) it
returns no pair: shared edges and the fans around shared vertices never
register as intersections.

## Visual reference

The dragon fixture (6588 facets), drawn translucent, with the segment
`p`–`q` of each of the 1191 pairs `find_self_intersections`
reports:

![the dragon drawn translucent, the 1191 self-intersection segments find_self_intersections reports in ember, most of them at the head and the feet](../../examples/images/ref-self-intersection-dragon.svg){.fossil-figure}

Drawn opaque, the dragon looks fine: the crossings are small folds
the eye does not catch. That's the *whole point*
of `find_self_intersections` — these are the defects that escape eye
inspection and that `is_watertight()` doesn't catch.

## Known limitations

- **`resolve_self_intersections` does not resolve.** It currently
  returns the arrangement unselected: on the two crossing cubes of the
  [Example](#example) it returns `BOOL_STATUS_OK` with 224 facets,
  volume 2.0 (the overlap counted twice) and a surface that is not
  watertight, and on `dragon.stl` it returns status 1 and leaves the
  surface unchanged; a library issue is pending. Until it is fixed,
  unite separate bodies with `boolean`, as the Example does, or wrap
  the input with [`alpha_wrap`](/guide/advanced/alpha-wrap).
- **Intended behaviour of the resolution.** Once fixed, `resolve_self_intersections` returns the topology of
  `BOOL_UNION` applied to the surface against itself (see [Pipeline](#pipeline-resolution)): the outer manifold only,
  so interior cavities formed by the self-intersections are dropped, not re-meshed.
- **Detection is `O(n_pairs)` in the worst case**, where `n_pairs` is
  the number of overlapping leaf-pairs found in the broad phase. For
  the dragon fixture (~6.6k facets) this is ~10× faster than brute-
  force `O(n²)` pair enumeration; for cleaner meshes the speedup is
  much larger because few leaves overlap at all.

## See also

- [booleans](/guide/advanced/booleans) — the resolution path is a
  one-liner over the boolean engine; same coplanar limitation
  applies.
- [`surface%sanitize`](/guide/api/surface-stl-object#sanitize) — local
  defect repair (degenerate / duplicate facets, normal orientation).
  Run this *before* `find_self_intersections` to filter out
  intersection records that exist only because of degenerate-facet
  noise.
- [alpha wrapping](/guide/advanced/alpha-wrap) — the heavyweight
  alternative for inputs too dirty for `resolve` to handle (large
  holes, missing facets that aren't just self-intersection
  artifacts).

## References

- Möller, *A Fast Triangle-Triangle Intersection Test*, Journal of
  Graphics Tools 2(2), 1997. The narrow-phase test until issue #69,
  replaced by the exact `pair_cut` (see [booleans](/guide/advanced/booleans)).
- Zhou, Grinspun, Zorin & Jacobson, *Mesh Arrangements for Solid
  Geometry*, SIGGRAPH 2016. The reference algorithm for the §1.1
  boolean that `resolve_self_intersections` reuses.
- libigl `igl::SelfIntersectMesh` — the closest open-source
  reference; FOSSIL's detector matches its broad-phase + narrow-phase
  decomposition.
