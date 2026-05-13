!< FOSSIL, Marching Cubes isosurface extraction (issue #18 §1.5).

module fossil_marching_cubes
!< FOSSIL, Marching Cubes isosurface extraction.
!<
!< Standard Lorensen-Cline 1987 Marching Cubes (the original — not MC33). For
!< each cube of the input scalar grid, compute an 8-bit case index from the
!< sign of `value(Vi) - iso` at each of the 8 corners; look up which of the
!< 12 cube edges have a sign change; for each active edge linearly interpolate
!< the iso-crossing point; assemble the resulting points into the triangles
!< prescribed by the per-case triangle table.
!<
!< Canonical conventions used (Paul Bourke's tables):
!<
!<   Vertex numbering of a unit cube (origin at V0):
!<     V0 = (0, 0, 0), V1 = (1, 0, 0), V2 = (1, 1, 0), V3 = (0, 1, 0)
!<     V4 = (0, 0, 1), V5 = (1, 0, 1), V6 = (1, 1, 1), V7 = (0, 1, 1)
!<
!<   Edge numbering (each edge connects two corners):
!<     E0=V0-V1, E1=V1-V2, E2=V2-V3, E3=V3-V0  (bottom face)
!<     E4=V4-V5, E5=V5-V6, E6=V6-V7, E7=V7-V4  (top face)
!<     E8=V0-V4, E9=V1-V5, E10=V2-V6, E11=V3-V7 (vertical edges)
!<
!<   Sign convention: bit `i` of the case index is set iff `value(Vi) < iso`
!<   ("vertex is inside the iso-surface"). Note: Bourke's tables are
!<   defined so that triangles are wound with normals pointing toward the
!<   "outside" (vertices with bit clear) when this convention is used, which
!<   for an SDF (negative inside, positive outside) means normals point
!<   **outward** of the bounded region — the convention FOSSIL uses for
!<   `compute_volume` to return a positive volume.
!<
!< Output is a flat list of independent triangles — adjacent cubes that share
!< an iso-crossing on a common edge produce duplicate vertices. Calling
!< `surface%adopt_facets(...)` on the output runs `analyze` (which includes
!< `connect_nearby_vertices`) and deduplicates near-coincident vertices,
!< producing a watertight mesh. Hash-based per-edge dedup inside the MC loop
!< would cut the output size by ~3x; deferred as future work.
!<
!< Known ambiguity (Lorensen-Cline limitation): on cases 105 and 150 (and
!< their rotations) the case-table choice can leave small holes in the iso-
!< surface where two diagonally-opposite corners are inside. Marching Cubes
!< 33 (Chernyaev 1995) resolves this with face- and body-tests; not
!< implemented here. For typical AMR-CFD use cases (smooth SDFs sampled on
!< well-resolved grids) the ambiguity rarely fires.

use fossil_facet_object,   only : facet_object
use penf,                  only : I4P, R8P
use vecfor,                only : vector_R8P

implicit none
private
public :: extract_isosurface
public :: MC_STATUS_OK, MC_STATUS_BAD_DIMENSIONS

integer(I4P), parameter :: MC_STATUS_OK             = 0_I4P
integer(I4P), parameter :: MC_STATUS_BAD_DIMENSIONS = 1_I4P  !< Grid has fewer than 2 cells along some axis.

! Edge endpoints — for each of the 12 edges, the indices (0..7) of its two corners.
! Used in the linear interpolation step to recover the 3D position of the iso-crossing.
integer(I4P), parameter :: EDGE_VERTEX(2, 0:11) = reshape([ &
   0, 1, &  !< E0
   1, 2, &  !< E1
   2, 3, &  !< E2
   3, 0, &  !< E3
   4, 5, &  !< E4
   5, 6, &  !< E5
   6, 7, &  !< E6
   7, 4, &  !< E7
   0, 4, &  !< E8
   1, 5, &  !< E9
   2, 6, &  !< E10
   3, 7  &  !< E11
   ], shape=[2, 12])

! Cube vertex offsets in (i, j, k) cell-relative coordinates.
! The cube at grid cell (i, j, k) has corners at (i+VOX(0,n), j+VOX(1,n), k+VOX(2,n)).
integer(I4P), parameter :: VOX(0:2, 0:7) = reshape([ &
   0, 0, 0, &  !< V0 = (i,   j,   k)
   1, 0, 0, &  !< V1 = (i+1, j,   k)
   1, 1, 0, &  !< V2 = (i+1, j+1, k)
   0, 1, 0, &  !< V3 = (i,   j+1, k)
   0, 0, 1, &  !< V4 = (i,   j,   k+1)
   1, 0, 1, &  !< V5 = (i+1, j,   k+1)
   1, 1, 1, &  !< V6 = (i+1, j+1, k+1)
   0, 1, 1  &  !< V7 = (i,   j+1, k+1)
   ], shape=[3, 8])

! 256-entry edge table: for each case index (0..255), a 12-bit mask of which
! cube edges have a sign change. Bit `e` set ⇔ edge E_e has a sign-flip between
! its two endpoints. Constants below are the canonical Lorensen-Cline /
! Bourke values (`http://paulbourke.net/geometry/polygonise/`, public domain).
integer(I4P), parameter :: MC_EDGE_TABLE(0:255) = [ &
   int(z'000', I4P), int(z'109', I4P), int(z'203', I4P), int(z'30a', I4P), &
   int(z'406', I4P), int(z'50f', I4P), int(z'605', I4P), int(z'70c', I4P), &
   int(z'80c', I4P), int(z'905', I4P), int(z'a0f', I4P), int(z'b06', I4P), &
   int(z'c0a', I4P), int(z'd03', I4P), int(z'e09', I4P), int(z'f00', I4P), &
   int(z'190', I4P), int(z'099', I4P), int(z'393', I4P), int(z'29a', I4P), &
   int(z'596', I4P), int(z'49f', I4P), int(z'795', I4P), int(z'69c', I4P), &
   int(z'99c', I4P), int(z'895', I4P), int(z'b9f', I4P), int(z'a96', I4P), &
   int(z'd9a', I4P), int(z'c93', I4P), int(z'f99', I4P), int(z'e90', I4P), &
   int(z'230', I4P), int(z'339', I4P), int(z'033', I4P), int(z'13a', I4P), &
   int(z'636', I4P), int(z'73f', I4P), int(z'435', I4P), int(z'53c', I4P), &
   int(z'a3c', I4P), int(z'b35', I4P), int(z'83f', I4P), int(z'936', I4P), &
   int(z'e3a', I4P), int(z'f33', I4P), int(z'c39', I4P), int(z'd30', I4P), &
   int(z'3a0', I4P), int(z'2a9', I4P), int(z'1a3', I4P), int(z'0aa', I4P), &
   int(z'7a6', I4P), int(z'6af', I4P), int(z'5a5', I4P), int(z'4ac', I4P), &
   int(z'bac', I4P), int(z'aa5', I4P), int(z'9af', I4P), int(z'8a6', I4P), &
   int(z'faa', I4P), int(z'ea3', I4P), int(z'da9', I4P), int(z'ca0', I4P), &
   int(z'460', I4P), int(z'569', I4P), int(z'663', I4P), int(z'76a', I4P), &
   int(z'066', I4P), int(z'16f', I4P), int(z'265', I4P), int(z'36c', I4P), &
   int(z'c6c', I4P), int(z'd65', I4P), int(z'e6f', I4P), int(z'f66', I4P), &
   int(z'86a', I4P), int(z'963', I4P), int(z'a69', I4P), int(z'b60', I4P), &
   int(z'5f0', I4P), int(z'4f9', I4P), int(z'7f3', I4P), int(z'6fa', I4P), &
   int(z'1f6', I4P), int(z'0ff', I4P), int(z'3f5', I4P), int(z'2fc', I4P), &
   int(z'dfc', I4P), int(z'cf5', I4P), int(z'fff', I4P), int(z'ef6', I4P), &
   int(z'9fa', I4P), int(z'8f3', I4P), int(z'bf9', I4P), int(z'af0', I4P), &
   int(z'650', I4P), int(z'759', I4P), int(z'453', I4P), int(z'55a', I4P), &
   int(z'256', I4P), int(z'35f', I4P), int(z'055', I4P), int(z'15c', I4P), &
   int(z'e5c', I4P), int(z'f55', I4P), int(z'c5f', I4P), int(z'd56', I4P), &
   int(z'a5a', I4P), int(z'b53', I4P), int(z'859', I4P), int(z'950', I4P), &
   int(z'7c0', I4P), int(z'6c9', I4P), int(z'5c3', I4P), int(z'4ca', I4P), &
   int(z'3c6', I4P), int(z'2cf', I4P), int(z'1c5', I4P), int(z'0cc', I4P), &
   int(z'fcc', I4P), int(z'ec5', I4P), int(z'dcf', I4P), int(z'cc6', I4P), &
   int(z'bca', I4P), int(z'ac3', I4P), int(z'9c9', I4P), int(z'8c0', I4P), &
   int(z'8c0', I4P), int(z'9c9', I4P), int(z'ac3', I4P), int(z'bca', I4P), &
   int(z'cc6', I4P), int(z'dcf', I4P), int(z'ec5', I4P), int(z'fcc', I4P), &
   int(z'0cc', I4P), int(z'1c5', I4P), int(z'2cf', I4P), int(z'3c6', I4P), &
   int(z'4ca', I4P), int(z'5c3', I4P), int(z'6c9', I4P), int(z'7c0', I4P), &
   int(z'950', I4P), int(z'859', I4P), int(z'b53', I4P), int(z'a5a', I4P), &
   int(z'd56', I4P), int(z'c5f', I4P), int(z'f55', I4P), int(z'e5c', I4P), &
   int(z'15c', I4P), int(z'055', I4P), int(z'35f', I4P), int(z'256', I4P), &
   int(z'55a', I4P), int(z'453', I4P), int(z'759', I4P), int(z'650', I4P), &
   int(z'af0', I4P), int(z'bf9', I4P), int(z'8f3', I4P), int(z'9fa', I4P), &
   int(z'ef6', I4P), int(z'fff', I4P), int(z'cf5', I4P), int(z'dfc', I4P), &
   int(z'2fc', I4P), int(z'3f5', I4P), int(z'0ff', I4P), int(z'1f6', I4P), &
   int(z'6fa', I4P), int(z'7f3', I4P), int(z'4f9', I4P), int(z'5f0', I4P), &
   int(z'b60', I4P), int(z'a69', I4P), int(z'963', I4P), int(z'86a', I4P), &
   int(z'f66', I4P), int(z'e6f', I4P), int(z'd65', I4P), int(z'c6c', I4P), &
   int(z'36c', I4P), int(z'265', I4P), int(z'16f', I4P), int(z'066', I4P), &
   int(z'76a', I4P), int(z'663', I4P), int(z'569', I4P), int(z'460', I4P), &
   int(z'ca0', I4P), int(z'da9', I4P), int(z'ea3', I4P), int(z'faa', I4P), &
   int(z'8a6', I4P), int(z'9af', I4P), int(z'aa5', I4P), int(z'bac', I4P), &
   int(z'4ac', I4P), int(z'5a5', I4P), int(z'6af', I4P), int(z'7a6', I4P), &
   int(z'0aa', I4P), int(z'1a3', I4P), int(z'2a9', I4P), int(z'3a0', I4P), &
   int(z'd30', I4P), int(z'c39', I4P), int(z'f33', I4P), int(z'e3a', I4P), &
   int(z'936', I4P), int(z'83f', I4P), int(z'b35', I4P), int(z'a3c', I4P), &
   int(z'53c', I4P), int(z'435', I4P), int(z'73f', I4P), int(z'636', I4P), &
   int(z'13a', I4P), int(z'033', I4P), int(z'339', I4P), int(z'230', I4P), &
   int(z'e90', I4P), int(z'f99', I4P), int(z'c93', I4P), int(z'd9a', I4P), &
   int(z'a96', I4P), int(z'b9f', I4P), int(z'895', I4P), int(z'99c', I4P), &
   int(z'69c', I4P), int(z'795', I4P), int(z'49f', I4P), int(z'596', I4P), &
   int(z'29a', I4P), int(z'393', I4P), int(z'099', I4P), int(z'190', I4P), &
   int(z'f00', I4P), int(z'e09', I4P), int(z'd03', I4P), int(z'c0a', I4P), &
   int(z'b06', I4P), int(z'a0f', I4P), int(z'905', I4P), int(z'80c', I4P), &
   int(z'70c', I4P), int(z'605', I4P), int(z'50f', I4P), int(z'406', I4P), &
   int(z'30a', I4P), int(z'203', I4P), int(z'109', I4P), int(z'000', I4P)  &
   ]

! 256 × 16 triangle table: for each case, up to 5 triangles each described by
! 3 edge indices (0..11). Sentinel value -1 marks the end of the list within
! a row. Layout is row-major: `MC_TRI_TABLE(slot, case)` where slot ∈ 0..15.
! Same Bourke convention as the edge table; the two must come from the same
! source for the case index → triangles mapping to be consistent.
integer(I4P), parameter :: MC_TRI_TABLE(0:15, 0:255) = reshape([ &
   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  0
    0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  1
    0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  2
    1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  3
    1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  4
    0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  5
    9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  6
    2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1, &  !  7
    3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  8
    0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  !  9
    1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 10
    1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1, &  ! 11
    3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 12
    0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1, &  ! 13
    3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 14
    9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 15
    4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 16
    4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 17
    0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 18
    4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 19
    1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 20
    3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1, &  ! 21
    9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 22
    2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1, &  ! 23
    8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 24
   11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 25
    9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1, &  ! 26
    4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1, &  ! 27
    3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 28
    1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1, &  ! 29
    4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1, &  ! 30
    4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1, &  ! 31
    9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 32
    9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 33
    0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 34
    8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 35
    1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 36
    3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 37
    5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 38
    2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1, &  ! 39
    9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 40
    0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 41
    0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1, &  ! 42
    2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1, &  ! 43
   10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 44
    4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1, &  ! 45
    5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1, &  ! 46
    5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1, &  ! 47
    9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 48
    9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 49
    0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 50
    1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 51
    9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 52
   10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1, &  ! 53
    8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1, &  ! 54
    2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 55
    7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 56
    9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1, &  ! 57
    2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1, &  ! 58
   11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 59
    9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1, &  ! 60
    5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1, &  ! 61
   11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1, &  ! 62
   11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 63
   10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 64
    0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 65
    9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 66
    1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 67
    1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 68
    1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1, &  ! 69
    9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 70
    5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1, &  ! 71
    2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 72
   11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 73
    0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 74
    5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1, &  ! 75
    6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 76
    0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1, &  ! 77
    3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1, &  ! 78
    6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1, &  ! 79
    5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 80
    4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1, &  ! 81
    1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 82
   10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1, &  ! 83
    6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1, &  ! 84
    1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1, &  ! 85
    8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1, &  ! 86
    7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1, &  ! 87
    3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 88
    5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1, &  ! 89
    0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1, &  ! 90
    9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1, &  ! 91
    8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1, &  ! 92
    5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1, &  ! 93
    0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1, &  ! 94
    6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1, &  ! 95
   10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 96
    4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 97
   10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 98
    8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1, &  ! 99
    1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 100
    3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1, &  ! 101
    0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 102
    8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 103
   10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 104
    0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1, &  ! 105
    3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1, &  ! 106
    6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1, &  ! 107
    9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1, &  ! 108
    8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1, &  ! 109
    3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 110
    6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 111
    7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1, &  ! 112
    0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1, &  ! 113
   10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1, &  ! 114
   10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 115
    1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1, &  ! 116
    2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1, &  ! 117
    7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 118
    7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 119
    2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1, &  ! 120
    2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1, &  ! 121
    1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1, &  ! 122
   11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1, &  ! 123
    8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1, &  ! 124
    0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 125
    7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1, &  ! 126
    7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 127
    7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 128
    3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 129
    0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 130
    8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 131
   10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 132
    1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 133
    2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 134
    6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1, &  ! 135
    7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 136
    7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 137
    2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 138
    1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1, &  ! 139
   10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 140
   10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1, &  ! 141
    0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1, &  ! 142
    7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 143
    6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 144
    3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 145
    8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 146
    9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1, &  ! 147
    6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 148
    1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1, &  ! 149
    4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1, &  ! 150
   10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1, &  ! 151
    8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 152
    0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 153
    1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1, &  ! 154
    1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 155
    8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1, &  ! 156
   10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1, &  ! 157
    4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1, &  ! 158
   10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 159
    4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 160
    0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1, &  ! 161
    5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1, &  ! 162
   11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1, &  ! 163
    9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1, &  ! 164
    6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1, &  ! 165
    7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1, &  ! 166
    3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1, &  ! 167
    7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 168
    9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1, &  ! 169
    3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1, &  ! 170
    6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1, &  ! 171
    9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1, &  ! 172
    1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1, &  ! 173
    4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1, &  ! 174
    7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1, &  ! 175
    6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 176
    3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1, &  ! 177
    0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1, &  ! 178
    6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 179
    1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1, &  ! 180
    0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1, &  ! 181
   11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1, &  ! 182
    6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1, &  ! 183
    5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1, &  ! 184
    9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1, &  ! 185
    1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1, &  ! 186
    1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 187
    1, 3, 6, 1, 6,10, 3, 8, 6, 5, 6, 9, 8, 9, 6,-1, &  ! 188
   10, 1, 0,10, 0, 6, 9, 5, 0, 5, 6, 0,-1,-1,-1,-1, &  ! 189
    0, 3, 8, 5, 6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 190
   10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 191
   11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 192
   11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 193
    5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 194
   10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1, &  ! 195
   11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 196
    0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1, &  ! 197
    9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1, &  ! 198
    7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1, &  ! 199
    2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 200
    8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1, &  ! 201
    9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1, &  ! 202
    9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1, &  ! 203
    1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 204
    0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1, &  ! 205
    9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1, &  ! 206
    9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 207
    5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1, &  ! 208
    5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1, &  ! 209
    0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1, &  ! 210
   10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1, &  ! 211
    2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1, &  ! 212
    0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1, &  ! 213
    0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1, &  ! 214
    9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 215
    2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1, &  ! 216
    5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 217
    3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1, &  ! 218
    5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1, &  ! 219
    8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1, &  ! 220
    0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 221
    8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1, &  ! 222
    9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 223
    4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1, &  ! 224
    0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1, &  ! 225
    1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1, &  ! 226
    3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1, &  ! 227
    4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1, &  ! 228
    9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1, &  ! 229
   11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1, &  ! 230
   11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1, &  ! 231
    2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1, &  ! 232
    9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1, &  ! 233
    3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1, &  ! 234
    1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 235
    4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1, &  ! 236
    4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1, &  ! 237
    4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 238
    4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 239
    9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 240
    3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1, &  ! 241
    0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1, &  ! 242
    3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 243
    1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1, &  ! 244
    3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1, &  ! 245
    0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 246
    3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 247
    2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1, &  ! 248
    9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 249
    2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1, &  ! 250
    1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 251
    1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 252
    0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 253
    0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, &  ! 254
   -1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1  &  ! 255
   ], shape=[16, 256])

contains

   subroutine extract_isosurface(values, bmin, bmax, iso, surface, status)
   !< Extract the iso=`iso` surface from a 3D scalar field `values` defined on
   !< a regular grid spanning the bounding box `[bmin, bmax]`.
   !<
   !< `values(1:nx, 1:ny, 1:nz)` is the scalar at grid corners; the cell at
   !< (i, j, k) has corners at indices (i..i+1, j..j+1, k..k+1) and spatial
   !< coordinates linearly interpolated within `[bmin, bmax]`.
   !<
   !< Returns the iso-surface as a flat array of `facet_object` (one allocate
   !< per call). Adjacent cubes may produce duplicate vertices on shared
   !< edges; downstream `surface%adopt_facets` runs `analyze` which dedupes
   !< via `connect_nearby_vertices`, producing a watertight mesh.
   !<
   !< Sign convention: a vertex with `value < iso` is "inside" the iso-surface,
   !< matching Bourke's table convention.
   real(R8P),                       intent(in)            :: values(:, :, :) !< Scalar field at grid corners.
   type(vector_R8P),                intent(in)            :: bmin, bmax      !< Spatial extent of the grid.
   real(R8P),                       intent(in)            :: iso             !< Iso-value to extract.
   type(facet_object), allocatable, intent(out)           :: surface(:)      !< Output triangle list.
   integer(I4P),                    intent(out), optional :: status          !< MC_STATUS_*.
   integer(I4P)                                           :: nx, ny, nz, i, j, k, c, t
   integer(I4P)                                           :: cube_case, edge_mask, idx
   real(R8P)                                              :: corner_val(0:7)
   type(vector_R8P)                                       :: corner_pos(0:7)
   type(vector_R8P)                                       :: edge_pt(0:11)
   real(R8P)                                              :: dx, dy, dz
   ! Output accumulator with geometric growth.
   type(facet_object), allocatable                        :: buf(:), tmp(:)
   integer(I4P)                                           :: cap, used
   type(facet_object)                                     :: f
   integer(I4P)                                           :: ev0, ev1
   real(R8P)                                              :: t_interp
   integer(I4P), parameter :: INITIAL_CAP = 1024_I4P

   if (present(status)) status = MC_STATUS_OK
   nx = size(values, dim=1); ny = size(values, dim=2); nz = size(values, dim=3)
   if (nx < 2 .or. ny < 2 .or. nz < 2) then
      if (present(status)) status = MC_STATUS_BAD_DIMENSIONS
      allocate(surface(0))
      return
   endif

   dx = (bmax%x - bmin%x) / real(nx - 1, R8P)
   dy = (bmax%y - bmin%y) / real(ny - 1, R8P)
   dz = (bmax%z - bmin%z) / real(nz - 1, R8P)

   cap  = INITIAL_CAP
   used = 0_I4P
   allocate(buf(cap))

   do k = 1, nz - 1
      do j = 1, ny - 1
         do i = 1, nx - 1
            ! Gather the 8 corner values + positions for this cell.
            do c = 0, 7
               corner_val(c) = values(i + VOX(0, c), j + VOX(1, c), k + VOX(2, c))
               corner_pos(c) = vector_R8P(bmin%x + (i - 1 + VOX(0, c)) * dx, &
                                          bmin%y + (j - 1 + VOX(1, c)) * dy, &
                                          bmin%z + (k - 1 + VOX(2, c)) * dz)
            enddo
            ! Build the 8-bit case index.
            cube_case = 0_I4P
            do c = 0, 7
               if (corner_val(c) < iso) cube_case = ior(cube_case, ishft(1_I4P, c))
            enddo
            edge_mask = MC_EDGE_TABLE(cube_case)
            if (edge_mask == 0_I4P) cycle  ! cell entirely inside or entirely outside
            ! Linear interpolation along each active edge.
            do c = 0, 11
               if (iand(edge_mask, ishft(1_I4P, c)) == 0_I4P) cycle
               ev0 = EDGE_VERTEX(1, c) ; ev1 = EDGE_VERTEX(2, c)
               t_interp = interp_t(corner_val(ev0), corner_val(ev1), iso)
               edge_pt(c)%x = corner_pos(ev0)%x + t_interp * (corner_pos(ev1)%x - corner_pos(ev0)%x)
               edge_pt(c)%y = corner_pos(ev0)%y + t_interp * (corner_pos(ev1)%y - corner_pos(ev0)%y)
               edge_pt(c)%z = corner_pos(ev0)%z + t_interp * (corner_pos(ev1)%z - corner_pos(ev0)%z)
            enddo
            ! Emit triangles per the case-specific table.
            t = 0
            do
               if (t > 15_I4P - 2_I4P) exit
               idx = MC_TRI_TABLE(t, cube_case)
               if (idx < 0_I4P) exit
               ! Vertex order: Bourke's tables list (a, b, c) with normal
               ! pointing toward the "outside" (vertices with bit clear,
               ! i.e. value > iso) under the (b-a) × (c-a) right-hand rule.
               ! For SDFs (value < iso = inside), this gives an outward
               ! normal which is what compute_volume requires for positive
               ! volumes — but `vecfor`'s `face_normal3_R8P` uses the
               ! `(p2-p1) × (p3-p1)` convention which when combined with
               ! Bourke's listing produces the **opposite** sign in our
               ! integration. Reverse v(2) ↔ v(3) to flip the normal.
               f%vertex(1) = edge_pt(idx)
               f%vertex(3) = edge_pt(MC_TRI_TABLE(t + 1, cube_case))
               f%vertex(2) = edge_pt(MC_TRI_TABLE(t + 2, cube_case))
               call f%compute_metrix
               if (used == cap) then
                  cap = 2 * cap
                  allocate(tmp(cap))
                  tmp(1:used) = buf(1:used)
                  call move_alloc(from=tmp, to=buf)
               endif
               used = used + 1
               buf(used) = f
               t = t + 3
            enddo
         enddo
      enddo
   enddo

   allocate(surface(used))
   surface(1:used) = buf(1:used)
   endsubroutine extract_isosurface

   pure function interp_t(va, vb, iso) result(t)
   !< Parameter t ∈ [0, 1] along the edge (a → b) where the iso-crossing lies:
   !<   `value(a) + t * (value(b) - value(a)) = iso`.
   !< Defensive against a == b (no crossing): returns 0.5.
   real(R8P), intent(in) :: va, vb, iso
   real(R8P)             :: t

   if (abs(vb - va) < tiny(1._R8P)) then
      t = 0.5_R8P
   else
      t = (iso - va) / (vb - va)
   endif
   endfunction interp_t

endmodule fossil_marching_cubes
