# Cut cells: solver support

FOSSIL's [cut cells](./cut-cells) give a block's fluid fractions, apertures,
centroids, wall measures and fluid **pieces**. A conservative cut-cell solver
needs more geometry before its first step: the neighbourhoods its small-cell
cure works on, merge targets for small, fresh and dead cells, moments above the
first, coarse levels, and the fluid fractions of the cell edges. Issue #70
(roadmap #68, WP1) adds them as pure geometry and graph work; the solver keeps
the arithmetic on its states and fluxes.

::: info Design and reference (issue #70)
Each section says what the literature prescribes, with the source, what FOSSIL
chooses where the literature leaves room, the interface, and what the tests
measured. Statements from general knowledge rather than a cited source are
marked *(gk)*. A worked example is in the cookbook:
[Neighbourhoods for a small-cell cure](/manual/cookbook#neighbourhoods-for-a-small-cell-cure).
:::

## At a glance

![an ellipse translating and turning through 24 by 12 cells: per frame, the cut cells shaded by their fluid fraction, each small piece in teal joined to its merge target in yellow, and the cells of the smallest piece's state-redistribution neighbourhood outlined in ember](../../examples/images/anim-cut-solver.svg){.fossil-figure}

Each frame recomputes the cut cells of a moving ellipse, their links, the
state-redistribution neighbourhoods and the merge targets: 8 to 13 small
pieces per frame (area below half the cell), the smallest down to 2.61e-05 of
its cell (`anim_cut_solver.f90`).

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
from p's side in CSR form (`cut_links_t`): `start(np + 1)`, `piece(:)`,
`face(:)` (1 to 6: −x, +x, −y, +y, −z, +z; 1 to 4 in 2-D) and `area(:)`,
sorted by face, then by the piece across. Both directions are stored, with the
same area.

<<< @/examples/snippets/ref_cut_solver-links.f90

![a thin plate 0.03 thick, tilted, across 8 by 6 cells of 0.125: a dot at the centroid of every piece, ember for the 7 cells the plate splits into two pieces, teal for the others; gold lines join each piece to the middle of the fluid edge interval it shares with a neighbour, and on to that neighbour: the pieces above the plate link only to pieces above it, those below only to those below](../../examples/images/ref-cut-solver-links.svg){.fossil-figure}

The plate splits 7 of the 48 cells: 55 pieces, 176 links counted in both directions, none across the plate
(`fig_ref_cut_solver_links.f90`).

- **2-D.** The fluid intervals of a cell edge are bounded by the crossings of
  the walk and the cell corners, the same points in the two cells that share
  the edge, so the intervals of two pieces across an edge are equal or
  disjoint: one link per pair of equal intervals, of their length. The cut
  cells now also give each piece's apertures (`piece_aperture(4, np)`) and
  its fluid intervals (`piece_iv_start`, `piece_iv_edge`, `piece_iv_lo`,
  `piece_iv_hi`).
- **3-D.** A connected fluid region of a face lies next to one piece on each
  side. The regions are the pieces of the plane section's 2-D cut cell (their
  areas are the same numbers as the face fractions), matched to the 3-D pieces
  on both sides through the fluid edge intervals they share, the edges of the
  face being cell edges. A one-piece side takes every region; two one-piece
  cells get one link of the face's fluid area. A region with no fluid
  interval (a fluid island inside a solid trace on the face) goes to the only
  piece of its side whose face area it leaves unexplained, else it is counted
  in `links%unmatched` and gets no link.
- **Untrusted cells** (`pieces = -1`) have no pieces and no links: their
  neighbours see a wall there.

Measured (`fossil_test_cut_links`, `fossil_test_cut_links_3d`):

| Case | Result |
|---|---|
| 2-D: circle on a generic and an aligned grid, thin plate | link sums = piece aperture × edge length to 1·10⁻¹⁶, exact symmetry, per-piece apertures add up to the cell's exactly; plate pieces never link across the plate |
| 3-D: tilted wedge, cone, L-prism (generic, aligned), thin plate, square tube (hole as a face island), 8 lattice tetrahedra; 113 298 links, 193 multi-piece cells | link sums = piece face fraction × face area to 1.2·10⁻¹⁴, exact symmetry, 0 unmatched regions; plate pieces never link across the plate |

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
- Each neighbourhood depends on its own piece only, so the result does not
  depend on the order in which pieces are visited.

<<< @/examples/snippets/ref_cut_solver-srd.f90

Measured (`fossil_test_cut_srd`; 2-D circle, thin plate and half-plane, 3-D tilted
wedge, cone and thin plate, τ = ½, 327 small pieces): every piece with volume is
in exactly N^ov_j neighbourhoods and the column sums are 1 exactly; Σ V̂ is the
total fluid volume to 10⁻¹⁴; x̂ is the weighted centroid of the members to
round-off; no neighbourhood is incomplete; 0 ≤ η^nb ≤ 1 and 0 < η^self ≤ 1;
on the half-plane x < 0.58 (cells of 0.1, fraction 0.2) every small cell merges
with the cell across its +x face only, the wall normal; no neighbourhood crosses
a plate.

## Neighbourhoods for flux redistribution

Flux redistribution (llm-wiki `methods/flux-redistribution.md`; Natarajan et
al. 2022 §2; Barrio Sanchez et al. 2024 App. A) averages over the fluid cells
of the 3ᵈ block around a cut cell "that can be reached from it", made precise
by Barrio Sanchez et al. as a monotone index path through valid cells.

- **Path:** from p, a sequence of links each of which moves one index by one
  step towards the target cell's offset, never away, through pieces that are
  trusted and have volume.
- **Weights:** κ_j / Σ κ over the reached set, p included (κ = |I| / hᵈ).
- A piece filling its cell has itself only, weight 1; a zero-volume piece has
  no set.
- `frd_neighbourhoods(links, piece_cell, piece_volume, x, y, frd, z=)` gives
  `frd_t` (`start`, `member` with the piece first, `weight`).

<<< @/examples/snippets/ref_cut_solver-frd.f90

- **Measured** (`fossil_test_cut_frd`; 2-D circle and thin plate, 3-D tilted
  wedge and thin plate, 654 cut pieces): every set equals a brute-force
  fixpoint that pulls pieces in from the cell grid (a formulation independent
  of the library's breadth-first search); the weights sum to 1 to 8.9·10⁻¹⁶;
  no set crosses a plate where the block lies within the plate (near its ends
  a monotone path goes around).

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
  the pieces that may not serve as targets (for a dead piece, the pieces of the
  covered cells), and the same rules apply on the cut cells at the right pose.
- A piece is **suitable** as a target if it has volume, is not small and is
  not excluded. An unsuitable single target (a small one, or a dying neighbour
  of a dying piece) starts a chain followed to the first suitable piece
  (`final_target`, status `MERGE_CHAIN`); a chain that loops or ends is
  reported (`MERGE_CYCLE`), a piece with no linked piece with volume too
  (`MERGE_ALONE`).
- `merge_targets(links, piece_cell, piece_volume, wet_normal, x, y, mt, z=, tau=, mask=, exclude=)`
  gives `merge_t` (`status`, `target`, `final_target`, `mix_start`,
  `mix_piece`, `mix_weight`).

<<< @/examples/snippets/ref_cut_solver-merge.f90

- **Measured** (`fossil_test_cut_merge`; 2-D half-plane, circle and a
  translating circle, 3-D tilted wedge, 250 served pieces): on the
  half-plane x < 0.58 every small cell targets the cell across its +x face,
  also as its only normal-directed target, weight 1; every target is linked
  and suitable or ends a chain at a suitable piece; normal-directed targets lie
  at the offsets with o_a ∈ {0, sign n_a}, at most 2ᵈ − 1, in distinct cells,
  weights summing to 1 to 2.2·10⁻¹⁶; the 15 fresh pieces of the moving circle
  (at the end pose) and its 18 dead ones (at the start pose, the covered cells
  excluded) all get targets, no cycle and no lone piece.

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
- **Interface:** `cut_cells_3d(..., moments_order=p)` gives `piece_moments(:, np)`,
  `piece_wall_moments(:, np)`, `piece_face_moments(:, 6, np)` and the face
  totals `face_x/y/z_moments`, in the order of `mono_exponents` (`fossil_cut_moments`:
  by degree, then decreasing exponent of X, then of Y); 2-D
  `curve%cut_cells(..., moments_order=p)` gives `piece_moments` and
  `piece_wall_moments`. The faces of a one-piece cell come from the plane
  sections (the same numbers as the face fractions), those of a multi-piece
  cell from the same boundary arcs as its face areas; a fluid cell is
  analytic. Kernels: 3-point Gauss on segments and the symmetric 7-point
  triangle rule in closed form, both exact to degree 5.
- **Measured** (`fossil_test_cut_moments`, order 4): summed over the pieces
  and re-centred to global monomials, the moments equal the block's minus the
  body's (integrated independently over tetrahedra from the origin) to
  9.9·10⁻¹⁴ relative on the tilted wedge, the aligned L-prism and the thin
  plate, and to 3.6·10⁻¹⁵ in 2-D (a 64-gon); the zeroth and first moments
  equal the piece volumes and centroids, the wall ones the wetted areas and
  centroids (2.0·10⁻¹⁵ of the cell); the face moments of the pieces add up to
  the plane sections' (5.7·10⁻¹⁵ of the face).
- **Quadrature** (`cut_cells_3d(..., quadrature=.true.)`, 2-D
  `curve%cut_cells(..., quadrature=.true.)`): per piece `quad_start`,
  `quad_point`, `quad_weight`, and `wall_quad_start`, `wall_quad_point`,
  `wall_quad_weight`; 3-D `quad_status` (0: the weights add up to the piece
  volume).
  - Wall rules: the symmetric 7-point rule on the wall fan triangles (3-D),
    3-point Gauss on the wall segments (2-D); positive weights, exact to
    degree 5.
  - Volume rules: signed cones from the cell centre over the oriented fluid
    boundary of the piece *(gk: the divergence-theorem cone construction)*:
    the wall fans reversed and, on every face, the triangles from the face
    centre to the face region's oriented boundary, from the plane section for
    a one-piece cell (exact, as the face fractions) and from the face arcs of
    the pieces otherwise; 36 points per cone (Duffy map, Gauss 4 × 3 × 3),
    exact to degree 4; a fluid cell takes the 3ᵈ Gauss rule. Weights are
    negative where a piece is not star-shaped about the cell centre: they are
    counted, not hidden. 2-D: signed triangles from the cell centre over the
    exported oriented boundary (`boundaries=.true.`: `piece_seg_start`,
    `piece_seg`, `piece_seg_wall`).
  - **Measured** (`fossil_test_cut_quadrature`; L-prism generic and
    aligned, thin plate, 4 lattice tetrahedra, a 64-gon; 5028 pieces): the
    rules integrate every monomial of degree ≤ 4 to the piece's moments to
    4.8·10⁻¹⁵ of the cell, the wall rules to 1.4·10⁻¹⁵; 360 pieces have
    negative volume weights; 1 piece has `quad_status` 1, a two-piece cell of
    a lattice tetrahedron whose face arcs (the clipped traces, as its face
    areas) do not close.

![a five-pointed star, the solid, on 4 by 4 cells of 0.25: the volume quadrature points of the cut pieces in teal where their weight is positive and in ember where it is negative, many of them inside the star's arms, and the wall points in gold along the star's edges](../../examples/images/ref-cut-solver-quadrature.svg){.fossil-figure}

On the star, 14 of the 20 pieces are cut; their rules have 544 volume points, 98 of them with a negative weight, and
72 wall points. The weights of every piece add up to its area within 3.82e-17 (`fig_ref_cut_solver_quadrature.f90`).
A signed triangle from the cell centre covers what lies between the centre and the boundary, solid included: its
points can lie outside the piece, and the triangles over the solid side carry the negative weights that cancel them.

## Coarsening

Coarse–fine cut cells are where re-redistribution fails at κ < 10⁻⁴ (llm-wiki
`methods/high-order-amr.md`). The coarse geometry at ratio r (nested edges, the
cell counts multiples of r) is built from the fine one, never recomputed:
`cut_coarsen_3d(fine, links, x, y, z, ratio, coarse, coarse_links, fine_to_coarse, status)`
(2-D `cut_coarsen_2d`, without z).

- Coarse pieces are the connected components of the fine pieces of a coarse
  cell, joined by the fine links inside it; `fine_to_coarse` maps every fine
  piece to its coarse one.
- Volumes, centroids, wetted measures, the fluid areas of the coarse faces,
  the moments (re-centred and rescaled to the coarse cell) and the links are
  sums of the fine ones; node flags and edge fractions come from the fine
  ones. A coarse cell with an untrusted fine cell is untrusted. Face piece
  counts are not derived (−1); quadrature rules are not coarsened
  (concatenating the fine rules is exact).

<<< @/examples/snippets/ref_cut_solver-coarsen.f90

- **Measured** (`fossil_test_cut_coarsen`; tilted wedge, L-prism, thin plate,
  a 64-gon and a thin plate in 2-D; 2612 coarse cells at ratio 2): against
  the cut cells computed directly on the coarse grid, the same pieces in every
  cell, piece volumes to 2.5·10⁻¹⁴ of the cell, fractions and wetted areas to
  2.5·10⁻¹⁴, node flags exactly, edge fractions to 1.1·10⁻¹⁶, moments to
  1.2·10⁻¹⁵, link areas per coarse face to 7.8·10⁻¹⁵.

## Edge fractions

`cut_cells_3d(..., edges=.true.)` gives the fluid fraction of every cell edge, `edge_x(nx, ny + 1, nz + 1)`,
`edge_y(nx + 1, ny, nz + 1)` and `edge_z(nx + 1, ny + 1, nz)` (one query per grid line, off by default:
the swept geometry calls the cut cells many times and needs no edges). Every call gives the inside flag of
every grid node, `node_in(nx + 1, ny + 1, nz + 1)` (the z lines); the 2-D cut cells give
`node_in(nx + 1, ny + 1)` (their edges are the apertures). Embedded-boundary constrained transport needs the edge fractions (llm-wiki
`div-free-amr/wiki/methods/embedded-boundary-constrained-transport.md`).

- **Exact topology.** Every grid line is queried once (`line_crossings`). A
  crossing belongs to the segment between two nodes when it lies after the
  perturbed lower node and before the perturbed upper one (`line_crossing_before`,
  one exact `orient3d`, ties by the global perturbation of issue #69). The
  first node's flag is the parity of the crossings before it; the next ones
  flip at each crossing. Lengths use the rounded crossing positions.
- The same node flags are the corner states of the 3-D pieces.
- **Measured** (`fossil_test_cut_edges`): against the apertures of the plane
  sections' 2-D cut cells (a separate code path: `slice_surface` and the 2-D
  walk), on the tilted wedge, the L-prism (generic and aligned), cube-01
  through its faces and 6 lattice tetrahedra: 209 132 edges agree to
  4.3·10⁻¹⁵. Of 60 231 node flags, 7 differ, all at nodes exactly on the
  surface of a lattice tetrahedron, where the sections decide on rounded
  crossing points ([limitations](./cut-cells#limitations)); the 3-D flags are
  exact predicates on the input triangles.

## Example

The regions shown above, as one program: the cube-01 body on a 3-D block of 8³ cells of 0.25, then the half-plane
x < 0.58 on 10 × 10 cells of 0.1 (the solid on the left, the fluid fraction of the cut column 0.2). The cookbook recipe
[Neighbourhoods for a small-cell cure](/manual/cookbook#neighbourhoods-for-a-small-cell-cure) does the same on a disc.

<<< @/examples/snippets/ref_cut_solver.f90

<<< @/examples/output/ref_cut_solver.ansi{ansi}

- **3-D links.** cube-01 lies on cell faces: 485 pieces, 36 of them of zero volume, where a face of the cube is the
  high plane of a cell. The first piece, in the corner cell, links to its three neighbours over whole faces (0.0625).
  `edges=.true.` gives the edge fractions: 584 of the 648 z edges hold fluid. No fluid region is unmatched.
- **SRD.** The 10 pieces of the cut column are small (fraction 0.20 < τ = ½). The one in cell (6, 5) merges with the
  piece of cell (7, 5) only, across its +x face, the wall normal; η^nb = (½ − 0.2) / 1 = 0.3.
- **FRD.** Its flux set holds the 6 pieces of the 3 × 3 block that are fluid; its own weight is κ_p / Σ κ =
  0.2 / 3.6 = 0.056.
- **Merge.** Its single target is the piece across +x; it is also its only normal-directed target, weight 1.
- **Moments.** In the cell's scaled coordinates the piece spans X ∈ [0.3, 0.5] and Y ∈ [−½, ½]: ∫1 = 0.2,
  ∫X = 0.08, ∫X² = 0.03267, ∫Y² = 0.01667 (per cell area), in the order of `mono_exponents`: 1, X, Y, X², XY, Y².
- **Coarsening.** At ratio 2 the 50 fine pieces become 15 coarse ones, the cut column's fraction 0.10, and the fine
  piece maps to the coarse piece 7, in coarse column 3.

## Left out

- **Redistribution itself**: applying the weights to states or fluxes, the
  slope limiter, re-redistribution registers. The graph and W in sparse form
  are what such a register needs.
- **Piece identity over time** for moving bodies beyond fresh and dead pieces:
  decided with the moving overset grids (WP3 of roadmap #68).
