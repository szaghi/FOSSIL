---
title: Features
---

# Features

## STL File I/O

- Load ASCII or binary STL files from disk
- Automatic format detection (`guess_format=.true.`) — no need to know the file format upfront
- Load with on-the-fly clipping: `clip_min` / `clip_max` bounding box arguments discard facets outside the box during load
- Save surfaces to ASCII or binary STL

## Surface Analysis

After calling `surface%analize`, the following data are computed and stored:

- Axis-aligned bounding box (`bmin`, `bmax`)
- Enclosed volume (`volume`) and centroid (`centroid`)
- Facet connectivity — which facets share each edge
- Disconnected-edge detection: facets with 1, 2, or 3 disconnected edges are catalogued separately
- AABB octree built automatically for acceleration

## Surface Manipulation

All manipulation methods operate directly on `surface_stl_object`:

- **Translate** — by a 3D vector (`delta=`) or by scalar components (`x=`, `y=`, `z=`)
- **Rotate** — around an arbitrary axis by angle in radians, or by a given rotation matrix
- **Mirror** — with respect to a plane defined by its normal, or by a given mirror matrix
- **Resize (scale)** — by a 3D vector factor or by scalar per-axis factors; optionally scale about the surface centroid
- **Clip** — discard facets outside an axis-aligned bounding box; the cut-off remainder is returned as a separate surface
- **Merge** — combine two STL surfaces into one (`merge_solids`)
- **Sanitize normals** — make all facet normals consistent (outward or inward)
- **Reverse normals** — flip all facet normals
- **Connect nearby vertices** — repair disconnected edges by snapping nearby vertices together

## Distance and Point-in-Polyhedra Queries

- Minimum unsigned (squared or Euclidean) distance from a point to the surface
- Signed distance (negative inside, positive outside) via:
  - Solid angle computation
  - Ray intersection counting
- Point-in-polyhedron test via solid angle or ray intersection
- `compute_mesh_distance` — compute distance from each node of a structured mesh block to the STL surface
- All distance queries benefit from AABB octree acceleration (up to 7× speedup over brute force)

## AABB Octree Acceleration

The Axis-Aligned Bounding Box (AABB) octree is built automatically during `analize`. It uses 8-child (octree) subdivision and is embedded directly in `surface_stl_object%aabb`. The refinement depth is user-configurable and reported in `statistics()`.

## `fossilizer` CLI

A command-line utility (`src/app/fossilizer.f90`) wraps the library for interactive STL processing:

- Load one or more STL files
- Apply clip, merge, mirror, rotate, translate, resize, sanitize, and connectivity repair operations
- Print per-file statistics to stdout
- Save the result to a named output file

## Compiler Support

| Compiler | Status |
|----------|--------|
| GNU gfortran ≥ 5.3 | Supported |
| Intel Fortran ≥ 16.x | Supported |
| NVIDIA HPC Fortran (nvfortran) | Supported |
| IBM XL Fortran | Not tested |
| g95 | Not tested |
| NAG Fortran | Not tested |

## Design Principles

- **Pure Fortran** — no C extensions, no system calls beyond standard I/O
- **OOP** — all functionality exposed as type-bound procedures on three well-defined types
- **TDD** — every public operation is exercised by automated tests in `src/tests/`
- **KISS** — simple, focused API without unnecessary abstractions
- **Free & Open Source** — multi-licensed for both FOSS and commercial use
