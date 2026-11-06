---
title: Named patches
---

# Named patches

## What it does

A CFD body carries several boundary conditions: walls, inlets, outlets, porous or heated regions. Exporters write
them as named `solid <name> … endsolid` blocks of one ASCII STL, or as one STL file per patch. FOSSIL keeps them
(issue #66): every facet carries a patch id (`facet%patch`), an index into the surface's table of patch names.

## Reading and writing

- `load_from_file` reads one patch per `solid` block of an ASCII file, named after the block, in file order; blocks
  with the same name merge. A file with one block, or a binary file, is one patch named after the header.
- `load_from_files(file_names)` reads several files as one surface, one patch per file (or per block), named after
  the `solid` name or, when that is empty, the file name.
- `save_into_file(is_ascii=.true.)` writes one block per patch; a single-patch surface writes exactly the file it
  wrote before. Binary STL has no names: all patches go into one.
- `get_patches_number()`, `patch_name(i)`, `facet_patch(f)` read the table; `set_patch(facets, name)` moves facets
  to a named patch, created if new.

Before issue #66 the ASCII reader stopped with "malformed STL" at the first `endsolid` of a multi-solid file.

## Through the pipeline

The patch id travels with the facet, so every operation that keeps or copies facets keeps it:

- `sanitize`: degenerate and duplicate facets are removed with their ids; of two duplicates in different patches the
  first one in the file wins.
- `rotate`, `translate`, `mirror`, `resize` move vertices only.
- `clip` keeps the patch of every facet, and the remainder keeps the name table.
- `merge_solids(other)` maps `other`'s patches into the surface's table by name: equal names merge, new names are
  appended.
- `fill_holes` gives the new facets the patch of the hole's boundary when every boundary facet is in one patch, and a
  patch named `filled` otherwise (a hole across a wall/inlet edge has no single right answer; the caller can
  `set_patch` it).

## Through booleans

`boolean` and `resolve_self_intersections` keep the patch of every facet: a sub-triangle of the arrangement takes
its parent facet's patch, so the area of a patch in the result is the area of its parent patch outside (or inside) the
other body, as the operation selects. The other surface's patches join the table by name before the arrangement:
`A ∪ B` of a body with patches `low`, `high` and one with `b_low`, `b_high` has all four, and equal names merge.

## Through remeshing

Patch boundaries are kept where they are, like the sharp edges of [feature edges](/guide/advanced/feature-edges), and
every facet keeps a patch:

- `isotropic_remesh`: the vertices of the edges between two patches are locked (whatever `preserve_features` says),
  so those edges are never flipped or collapsed; a split facet's halves keep its patch, and a midpoint on a boundary
  edge becomes a boundary vertex.
- `decimate`: two boundary vertices never merge, and an interior vertex merges onto a boundary one at its place; the
  boundary vertices are the input ones exactly.
- `smooth`: the boundary vertices do not move.
- `alpha_wrap`: each facet of the wrap takes the patch of the input facet closest to its centroid.

## API

See [Patches](/guide/api/surface-stl-object#patches) in the API reference.
