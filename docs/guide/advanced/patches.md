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

## API

See [Patches](/guide/api/surface-stl-object#patches) in the API reference.
