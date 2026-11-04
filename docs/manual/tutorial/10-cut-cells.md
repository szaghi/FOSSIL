# 10. Cut cells

A conservative cut-cell solver needs, for every cell the wall crosses, the fluid fraction of the cell, the fluid
fraction of each face, and the separate fluid pieces when the wall splits a cell. FOSSIL computes them exactly from
the geometry. This chapter does it in 2-D on the section of chapter 7 and in 3-D on the dragon.

<<< @/examples/snippets/tut10_cut_cells-cut2d.f90

- `cut_cells(x=, y=, cut=)` takes the cell edges of a block, strictly increasing in each direction, and fills a
  `cut_cells_2d_t`. The fluid is outside the solid. Per cell (nx × ny arrays):
  - `alpha`: the fluid area fraction, 1 in the fluid, 0 in the solid;
  - `aperture(4, :, :)`: the fluid fraction of the four edges, in the order -u, +u, -v, +v;
  - `centroid`, `wet_length`, `wet_centroid`, `wet_normal`: the fluid centroid and the wetted wall;
  - `pieces`: the number of separate fluid pieces of the cell, with their measures in `piece_area`,
    `piece_centroid`, ... from index `first_piece`.
- The topology is exact (the crossings are decided with the exact `orient2d`, a point on a grid line belongs to the
  low cell), the measures exact to round-off. A cell whose boundary does not close (an open loop) is flagged in
  `closed`, with `alpha = -1`.

<<< @/examples/snippets/tut10_cut_cells-cut3d.f90

- `cut_cells_3d(surface=, x=, y=, z=, cut=)` is a free routine (from `use fossil`), with a `cut_cells_3d_t` result:
  `alpha` (volume fraction), `face_x`, `face_y`, `face_z` (face fractions), the wetted area, centroid and normal, and
  the fluid pieces: `pieces(i, j, k)`, numbered from `first_piece(i, j, k)`, each with `piece_volume`,
  `piece_centroid`, its wall terms and `piece_face(6, :)`, the fluid fraction of its six faces.
- `pieces = -1` flags a cell that cannot be trusted: a boundary edge of the mesh crosses it, or its faces do not
  close. Such a cell has no measures.

## Running it

<<< @/examples/output/tut10.ansi{ansi}

- 2-D: 91 cells are cut, 2 of them hold two fluid pieces. Summed over the block, the solid area equals the area of the curve to the printed digits.
  The cell with the longest wall, (23, 41), is 33 % fluid; its -v edge is entirely solid.
- 3-D: 895 cells are cut and 75 have two or more fluid pieces. 17 are flagged: the dragon is open (chapter 2), and
  the cells at its open edges get no measures.
- In cell (13, 7, 6) the fluid is 6.856 % of the cell, in two pieces of 6.533 % and 0.323 %: their volumes add up
  to the cell's.

::: tip What you learned
`curve%cut_cells` and `cut_cells_3d` give exact fluid fractions, apertures and fluid pieces per cell, and flag the
cells an open surface makes untrustworthy. Reference: [Cut-cell geometry](/guide/advanced/cut-cells).
:::

Next: [11. Moving bodies](./11-moving-bodies).
