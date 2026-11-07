# Cut cells: solver support

FOSSIL's [cut cells](./cut-cells) give a block's fluid fractions, apertures,
centroids, wall measures and fluid **pieces**. A conservative cut-cell solver
needs more geometry before its first step: the neighbourhoods its small-cell
cure works on, merge targets for small, fresh and dead cells, moments above the
first, coarse levels, and the fluid fractions of the cell edges. Issue #70
(roadmap #68, WP1) adds them as pure geometry and graph work; the solver keeps
the arithmetic on its states and fluxes.

::: info Design note (issue #70, step 2)
This page fixes the definitions, the tie rules and the API shapes before the
code. Each section says what the literature prescribes, with the source, and
what FOSSIL chooses where the literature leaves room. Statements from general
knowledge rather than a cited source are marked *(gk)*. Later steps of #70
replace the plans here with the measured results.
:::

## Principles

1. **The piece is the unit.** Every neighbourhood, link and target is indexed by
   the pieces `1..cut%np`, not by cells. A cell with one piece is that piece; a
   regular cell is one piece with κ = 1. AMReX's embedded boundaries are
   single-valued per cell *(gk)*; a thin plate across a cell gives FOSSIL two
   pieces with two states, and the graph keeps them apart.
2. **Trusted cells only.** A cell with `pieces = -1` has no pieces, so it takes
   no part in the graph; its neighbours see it as a wall. Every output reports
   such gaps through a status, never through a silent fallback.
3. **Zero-volume pieces carry no state.** A wall on the high plane of a cell
   leaves a fluid face bounding no fluid volume ([limitations](./cut-cells#limitations)).
   These pieces get links (their face is open) but are never members, centres
   or targets of a neighbourhood; a solver treats their faces as walls.
4. **Deterministic.** Every choice follows a fixed order on (cell index in
   x-fastest order, piece index): the same input gives the same output, with
   or without OpenMP.
5. **CPU only.** The cut cells are computed on the CPU, and this is
   preprocessing per block, not a batched point query: no device path.
6. **New free modules**, the cut-cell types unchanged in shape:
   `fossil_cut_graph` (links, neighbourhoods, targets), `fossil_cut_moments`
   (moments, quadrature), `fossil_cut_coarsen` (coarse levels). The links and
   moments need the transient wall polygons of `cut_cells_3d`, so it takes
   optional arguments that compute them only when asked.

## The piece graph

A **link** (p, q, f, A) says that piece p of cell c and piece q of the cell
across face f of c share a fluid region of that face, of area A. It is stored
from p's side in CSR form: `link_start(np + 1)`, `link_piece(:)`,
`link_face(:)` (1 to 6: −x, +x, −y, +y, −z, +z; 1 to 4 in 2-D) and
`link_area(:)`. Both directions are stored, with the same area.

- **Which pieces meet.** A connected fluid region of a face lies next to one
  piece on each side: the fluid next to a connected region is connected.
  A face region is matched to its pieces through the fluid intervals of the
  face's four edges, which both cells compute exactly and identically (the
  edge crossings under the global perturbation of issue #69).
- **Areas.** When the pieces of one side that touch the face are a single
  piece, every link of the face takes the face area of the piece on the other
  side (exact: the same numbers as `piece_face`). When both sides have several
  pieces on the face, the area of each region comes from the plane section's
  2-D pieces; a region with no fluid edge interval (a fluid island inside a
  solid trace) is matched by its area, and an ambiguous match is reported.
- **Invariants (the gate):** the links of a face sum to the face fraction ×
  face area, symmetry (p → q ⇔ q → p, same area), and per piece the outward
  face areas minus the wall normal integral vanish per axis.

## Neighbourhoods for state redistribution

State redistribution ([Berger & Giuliani 2021](https://doi.org/10.1016/j.jcp.2020.109820),
§4; llm-wiki `high-order-capturing-schemes/wiki/methods/state-redistribution.md`)
merges each small cell with neighbours "preferably in the direction of the
wall normal" until the merged volume reaches ½ hᵈ, falling back to a centred
3ᵈ, then 5ᵈ, tile. The rule for choosing neighbours inside the tile is left
open; FOSSIL's is:

- A piece p is **small** if |I_p| < τ hᵈ (τ = ½ by default, an argument; the
  threshold is a 1-D stability result, unproved in 2-D and 3-D, by the
  authors' own account). A piece that is not small has N_p = {p}.
- **Candidates** are the pieces of the 3ᵈ tile around p's cell reachable from p
  through links between pieces of the tile, zero-volume pieces excluded.
- **Order:** by the alignment of the cell offset o ∈ {−1, 0, 1}ᵈ with the unit
  wall normal n̂ of p (into the fluid), s = n̂ · o / |o|, largest first; ties by
  smaller |o|₁, then lexicographic o, then piece index.
- **Growth:** candidates join in that order until Σ |I_j| ≥ τ hᵈ. If the 3ᵈ
  tile is exhausted, the 5ᵈ tile is used with the same rule; if that is
  exhausted too, the neighbourhood is kept with status *incomplete*.
- **Outputs:** the neighbourhoods in CSR form (p first), the overlap counts
  N^ov_j (the neighbourhoods containing j, its own included), the weighted
  volumes V̂_p = Σ_{j∈N_p} |I_j| / N^ov_j and centroids x̂_p, and the weighted
  factors of [Giuliani et al. 2022](https://doi.org/10.1016/j.jcp.2022.111305)
  (§3.1): η^nb_p = (τ hᵈ − |I_p|) / Σ_{j∈N_p∖p} |I_j| for small p, else 0, and
  η^self_p = 1 − (1/N^ov_p) Σ_{k≠p, p∈N_k} η^nb_k.
- **Gate:** unit column sums of the redistribution matrix (the conservation
  criterion of Giuliani et al.: Σ_{p ∋ j} 1/N^ov_j = 1 for every j), the
  threshold reached by every complete neighbourhood, x̂ the weighted centroid
  of its members, and the same result with the cells traversed in reverse.

## Neighbourhoods for flux redistribution

Flux redistribution (llm-wiki `methods/flux-redistribution.md`; Natarajan et
al. 2022 §2; Barrio Sanchez et al. 2024 App. A) averages over the fluid cells
of the 3ᵈ block around a cut cell "that can be reached from it", made precise
by Barrio Sanchez et al. as a monotone index path through valid cells.

- **Path:** from p, a sequence of links each of which moves one index by one
  step towards the target cell's offset, never away, through pieces that are
  trusted and have volume.
- **Weights:** κ_j / Σ κ over the reached set, p included (κ = |I| / hᵈ).
- **Gate:** equality with brute-force enumeration of the monotone paths,
  weights summing to 1, and a wall between two cells of the block cutting the
  path.

## Merge targets

Two families, both reported per piece:

- **One target** (cell merging; Falcovitz et al. 1997 in llm-wiki
  `methods/cut-cell-finite-volume.md`: covered, newly exposed and small cells
  merge "with an adjacent wet cell"). The target is the linked piece, not
  small, that maximises n̂ · o / |o|, ties by the larger link area, then the
  lower index. When every linked piece is small, the chain is followed to the
  first piece that is not; a cycle is reported.
- **Normal-directed targets** (conservative mixing; Bezgin et al. 2023 §3.6 in
  llm-wiki `weno-teno/wiki/methods/level-set-method.md`: "up to seven target
  cells in 3D chosen from the interface normal"). The targets are the pieces at
  the offsets o with o_a ∈ {0, sign n_a}, o ≠ 0 (3 in 2-D, 7 in 3-D), reached
  by a monotone path. The weights of the cited scheme are not in the wiki;
  FOSSIL reports β_o ∝ (Π_{a: o_a ≠ 0} |n_a|)^{2/|supp o|} *(gk, after Hu et
  al. 2006)*, normalised over the targets found, and leaves the volume-fraction
  factor to the solver.
- **Fresh and dead pieces** of a moving body come from `swept_cells_2d/3d`
  (`cleared`, `covered`): the caller passes the pieces to serve as a mask and
  the cells that may not serve as targets (covered at the end of the step for a
  dead piece), and the same rules apply on the cut cells at the right pose.

## Moments and quadrature

- **Moments** of a piece: ∫ ((x − x_c)/h_x)^a ((y − y_c)/h_y)^b ((z − z_c)/h_z)^c dV
  for a + b + c ≤ p (p ≤ 4), about the cell centre x_c, scaled by the cell size
  so they stay O(|I|) far from the origin; the same for the face fluid regions
  and the wall polygons. Volume moments by the divergence theorem along x,
  ∫_V x^α dV = (1/(α_x + 1)) ∮ x^{α + e_x} n_x dA: the wall polygons (fan
  triangles, a rule exact to degree p + 1) and the +x face only (local
  coordinates put the −x face at 0); face moments by Green's theorem over the
  same boundary arcs as the face fractions; 2-D by Green. First moments must
  equal the existing centroids.
- **Wall quadrature:** fan triangles × a symmetric triangle rule, positive
  weights, exact to degree p.
- **Volume quadrature:** cones from the piece centroid over its boundary
  polygons, split into tetrahedra (triangles in 2-D), × a simplex rule *(gk)*.
  Exact to degree p; weights can be negative where the piece is not
  star-shaped about its centroid, and the minimum weight is reported.

## Coarsening

Coarse–fine cut cells are where re-redistribution fails at κ < 10⁻⁴ (llm-wiki
`methods/high-order-amr.md`). The coarse geometry at ratio r (default 2,
nested edges checked) is built from the fine one, never recomputed:

- coarse pieces are the components of the fine pieces linked inside the coarse
  cell, and the fine → coarse piece map is an output;
- volumes, face fractions, wall measures, moments and links are sums of the
  fine ones in a fixed order, so fine fractions add up to coarse ones by
  construction.
- **Gate:** against `cut_cells_3d` run directly on the coarse grid, the same
  piece count in every cell and the same links, volumes and areas to 10⁻¹⁴ of
  the cell.

## Edge fractions

The fluid length of each cell edge (3-D: `edge_x(nx, ny + 1, nz + 1)`,
`edge_y`, `edge_z`; 2-D: the node inside flags), from the exact edge
crossings the pieces already use. Embedded-boundary constrained transport needs
them (llm-wiki `div-free-amr/wiki/methods/embedded-boundary-constrained-transport.md`).

## Left out

- **Redistribution itself**: applying the weights to states or fluxes, the
  slope limiter, re-redistribution registers. The graph and W in sparse form
  are what such a register needs.
- **Piece identity over time** for moving bodies beyond fresh and dead pieces:
  decided with the moving overset grids (WP3 of roadmap #68).
