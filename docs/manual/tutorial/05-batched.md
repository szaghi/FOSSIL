# 5. Batched queries

A grid asks for the distance of every one of its points, not of five. `distance_many` takes the whole batch in one
call. This chapter runs it on 16³ points over the domain, on both kinds of AABB tree and with both device-ready sign
algorithms.

<<< @/examples/snippets/tut05_batched-grid.f90

<<< @/examples/snippets/tut05_batched-trees.f90

- `distance_many(points=, distances=)` is a loop of `compute_distance` over `points`, with the same `is_signed`,
  `sign_algorithm` and `is_square_root` for every point (unsigned and squared by default, as `distance`).
  `distances` is allocated by the caller, with the size of `points`; a size mismatch gives
  `status = STATUS_INVALID_INPUT` and leaves `distances` untouched.
- The loop is an OpenMP `parallel do`: the tree is read-only and each iteration writes its own slot. Built with
  `-fopenmp` it runs on all threads, without it serially, with the same results.
- `aabb_tree_kind=` on `load_from_file` (also on `analyze` and `adopt_facets`) selects the tree:
  `AABB_TREE_SAH_BVH` (the default), a binary tree split by the surface-area heuristic, or `AABB_TREE_OCTREE`, the
  older 8-way space partition. The octree is kept for comparison; the BVH is the faster one and the only one with a
  device path (chapter 12).

<<< @/examples/snippets/tut05_batched-signs.f90

## Running it

<<< @/examples/output/tut05.ansi{ansi}

- The two trees give the same distance magnitudes bit for bit: the traversal never prunes the true closest facet.
  One sign differs. At that point several facets are equally close, the trees report different ones, and on a
  non-manifold surface their pseudo-normals disagree.
- The two sign algorithms give the same magnitudes and different signs at 50 of the 4096 points. The pseudo-normal
  sign puts 152 points inside, one of them at d = -0.647, deeper than the dragon is wide: a sign read from the
  closest feature is wrong where the surface is open. The fast winding number puts 108 points inside, none deeper
  than 0.070. On a surface that is not closed, use `SIGN_FAST_WINDING_NUMBER`; the rest of the tutorial does.

The timing is not printed: it depends on the machine. The figures are in
[GPU distance queries](/guide/advanced/device-distance) (the BVH is about 110 times faster than the octree on a mesh
of 24k facets).

::: tip What you learned
`distance_many` answers a whole grid in one OpenMP call; the SAH BVH and the octree agree on the distances; the
pseudo-normal sign needs a closed surface, the fast winding number does not. Reference:
[surface_stl_object](/guide/api/surface-stl-object).
:::

Next: [6. A distance field](./06-distance-field).
