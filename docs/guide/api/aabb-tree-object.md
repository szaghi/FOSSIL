---
title: aabb_tree_object
---

# `aabb_tree_object`

The AABB (axis-aligned bounding-box) tree is the acceleration structure that
makes distance queries and inside tests run in O(log N) instead of O(N). Every
[`surface_stl_object`](/guide/api/surface-stl-object) owns one as
`surface%aabb` — you usually don't interact with it directly. You configure
it by passing the optional `aabb_tree_kind` and `aabb_refinement_levels`
arguments to [`load_from_file`](/guide/api/surface-stl-object#load_from_file)
or [`analyze`](/guide/api/surface-stl-object#analyze).

This page documents the **inspection and benchmarking** API on the tree
itself: useful for tooling, regression tests, and the occasional debugging
session, but not for normal use.

## Two kinds, one interface

| Kind                  | Build cost | Query speed (typical)              | When picked                                              |
| --------------------- | ---------- | ---------------------------------- | -------------------------------------------------------- |
| `AABB_TREE_SAH_BVH`   | O(N log N) | ~110× the octree on dragon-fine    | **Default.** Always, unless you have a benchmarking reason. |
| `AABB_TREE_OCTREE`    | O(N)       | Slower; depth-limited              | Benchmarking the BVH, debugging, teaching.                |

Both kinds use the same traversal algorithm — best-first descent with
squared-distance pruning. As a consequence, **both kinds produce bit-exact
distance results** on a given mesh and a given query point. The choice is
purely a performance trade-off.

[[toc]]

---

## Accessors

### `get_tree_kind`

```fortran
kind = surface%aabb%get_tree_kind()
```

**Purpose.** Return the kind currently configured, as one of
[`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh) or
[`AABB_TREE_OCTREE`](/guide/api/constants#aabb_tree_octree).

**Arguments.** None.

**Returns.** `integer(I4P)`.

**Example.**

```fortran
program ex_aabb_get_tree_kind
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
select case (surface%aabb%get_tree_kind())
case (AABB_TREE_SAH_BVH) ; print '(A)', 'using SAH BVH'
case (AABB_TREE_OCTREE)  ; print '(A)', 'using octree'
end select
end program ex_aabb_get_tree_kind
```

**See also.** [`set_tree_kind`](#set_tree_kind),
[`get_refinement_levels`](#get_refinement_levels).

---

### `set_tree_kind`

```fortran
call surface%aabb%set_tree_kind(tree_kind=AABB_TREE_OCTREE)
call surface%analyze
```

**Purpose.** Change the kind in-place. **Does not rebuild the tree** —
`set_tree_kind` invalidates the build state; the actual rebuild happens at
the next call to [`surface%analyze`](/guide/api/surface-stl-object#analyze).

**Arguments.**

| Argument | Intent | Type           | Notes                                                              |
| -------- | ------ | -------------- | ------------------------------------------------------------------ |
| `kind`   | `in`   | `integer(I4P)` | (required) `AABB_TREE_SAH_BVH` or `AABB_TREE_OCTREE`.              |

**Example.** Swap to octree to benchmark.

```fortran
program ex_aabb_set_tree_kind
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%aabb%set_tree_kind(tree_kind=AABB_TREE_OCTREE)
call surface%analyze(aabb_refinement_levels=AABB_AUTO_REFINEMENT)
print '(A, I0)', 'octree depth = ', surface%aabb%get_refinement_levels()
end program ex_aabb_set_tree_kind
```

**Notes.**

- The simpler idiom is to pass `aabb_tree_kind=` directly to
  `load_from_file` or `analyze` and skip this method altogether.

**See also.** [`get_tree_kind`](#get_tree_kind),
[`surface%analyze`](/guide/api/surface-stl-object#analyze).

---

### `get_refinement_levels`

```fortran
n = surface%aabb%get_refinement_levels()
```

**Purpose.** Return the octree depth. Meaningless when the kind is
`AABB_TREE_SAH_BVH` (the BVH ignores the field).

**Returns.** `integer(I4P)`.

**Example.** See [`set_tree_kind`](#set_tree_kind).

**See also.** [`set_refinement_levels`](#set_refinement_levels),
[`AABB_AUTO_REFINEMENT`](/guide/api/constants#aabb_auto_refinement).

---

### `set_refinement_levels`

```fortran
call surface%aabb%set_refinement_levels(refinement_levels=4)
call surface%analyze
```

**Purpose.** Set the octree depth manually, or pass `AABB_AUTO_REFINEMENT` to
ask the heuristic to pick. Invalidates the build state; the next
`surface%analyze` does the actual rebuild.

**Arguments.**

| Argument            | Intent | Type           | Notes                                                                              |
| ------------------- | ------ | -------------- | ---------------------------------------------------------------------------------- |
| `refinement_levels` | `in`   | `integer(I4P)` | (required) A non-negative integer depth, or the sentinel `AABB_AUTO_REFINEMENT`.   |

**Notes.**

- Has no effect when the kind is `AABB_TREE_SAH_BVH`.
- Depth 0 corresponds to a single-node tree (effectively brute force).
- The auto-pick heuristic is `ceil(log8(N/64))`, clamped to `[1, 6]`.

**See also.** [`get_refinement_levels`](#get_refinement_levels),
[`AABB_AUTO_REFINEMENT`](/guide/api/constants#aabb_auto_refinement).

---

### `get_nodes_number`

```fortran
n = surface%aabb%get_nodes_number()
```

**Purpose.** Total number of tree nodes (internal + leaves) after the last
build.

**Returns.** `integer(I4P)`.

**Example.** Compare BVH vs octree size on the same mesh.

```fortran
program ex_aabb_get_nodes_number
use fossil
use penf, only : I4P, R8P
implicit none
type(surface_stl_object) :: surface
integer(I4P)             :: n_bvh, n_oct

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
n_bvh = surface%aabb%get_nodes_number()
call surface%aabb%set_tree_kind(tree_kind=AABB_TREE_OCTREE)
call surface%analyze(aabb_refinement_levels=AABB_AUTO_REFINEMENT)
n_oct = surface%aabb%get_nodes_number()
print '(A, I0, A, I0)', 'nodes: BVH=', n_bvh, ' octree=', n_oct
end program ex_aabb_get_nodes_number
```

**See also.** [`node_at`](#node_at).

---

### `get_is_initialized`

```fortran
ok = surface%aabb%get_is_initialized()
```

**Purpose.** Return `.true.` once the tree has been built (after
`surface%analyze` or `surface%load_from_file`). Useful as a guard in code
that may run before a surface has been loaded.

**Returns.** `logical`.

---

### `node_at`

```fortran
node => surface%aabb%node_at(i)
```

**Purpose.** Return a pointer to node `i` for tree-walking diagnostics
(visualisation, regression tests). Returns `null()` for out-of-range `i`.

**Arguments.**

| Argument | Intent | Type           | Notes                                |
| -------- | ------ | -------------- | ------------------------------------ |
| `i`      | `in`   | `integer(I4P)` | (required) 1-based node index.       |

**Returns.** `type(aabb_node_object), pointer` — or `null()` if out of range.

**Notes.**

- The `aabb_node_object` type holds bounding-box extents and an integer list
  of contained facet ids. It is **not** re-exported by `use fossil`; access
  it through the pointer you get here.
- The pointer is read-only by convention; mutating the node invalidates the
  tree.

**See also.** [`get_nodes_number`](#get_nodes_number).

---

## Dispatch knob *(benchmarking only)*

Two named `logical` constants choose between tree-accelerated and brute-force
query dispatch. They live in `fossil_aabb_tree_object` and are **not**
re-exported by `use fossil`; to use them, import the submodule explicitly.

### `get_use_index` / `set_use_index`

```fortran
use fossil_aabb_tree_object, only : AABB_USE_INDEX, AABB_USE_BRUTE_FORCE
...
call surface%aabb%set_use_index(use_index=AABB_USE_BRUTE_FORCE)
```

**Purpose.** Switch between the tree-accelerated query path
(`AABB_USE_INDEX`, the default — `.true.`) and a brute-force O(N) scan over
every facet (`AABB_USE_BRUTE_FORCE` — `.false.`). The brute-force mode is
strictly for benchmarking the speed-up the tree provides; both modes produce
bit-identical results.

**Arguments (`set_use_index`).**

| Argument     | Intent | Type      | Notes                                       |
| ------------ | ------ | --------- | ------------------------------------------- |
| `use_index`  | `in`   | `logical` | (required) `AABB_USE_INDEX` (`.true.`) or `AABB_USE_BRUTE_FORCE` (`.false.`). |

**Returns (`get_use_index`).** `logical`.

**Example.** Time-compare tree vs brute force.

```fortran
program ex_aabb_brute_force
use fossil
use penf, only : I4P, R8P
use fossil_aabb_tree_object, only : AABB_USE_INDEX, AABB_USE_BRUTE_FORCE
use vecfor_R8P, only : vector_R8P, ex_R8P
implicit none
type(surface_stl_object) :: surface
real(R8P)                :: d_tree, d_brute

call surface%load_from_file(file_name='cube.stl', guess_format=.true.)
call surface%aabb%set_use_index(use_index=AABB_USE_INDEX)
d_tree  = surface%distance(point=2._R8P * ex_R8P, is_square_root=.true.)
call surface%aabb%set_use_index(use_index=AABB_USE_BRUTE_FORCE)
d_brute = surface%distance(point=2._R8P * ex_R8P, is_square_root=.true.)
print '(A, ES20.12)', 'tree   = ', d_tree
print '(A, ES20.12)', 'brute  = ', d_brute
print '(A, ES20.12)', 'diff   = ', abs(d_tree - d_brute)
end program ex_aabb_brute_force
```

**Notes.**

- The two constants are `.true.` and `.false.` respectively; the named form
  exists purely for readability at call sites.
- Brute force is faster than tree only on degenerate inputs (e.g., a surface
  with one facet) — the regression suite includes these to verify the
  pruning predicate stays exact.

**See also.** [`AABB_TREE_SAH_BVH`](/guide/api/constants#aabb_tree_sah_bvh),
[`surface%distance`](/guide/api/surface-stl-object#distance).
