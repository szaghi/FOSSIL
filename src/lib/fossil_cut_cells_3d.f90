!< FOSSIL, 3-D cut cells of a surface on a rectilinear block (issue #40 steps 4-5).
module fossil_cut_cells_3d
!< FOSSIL, 3-D cut cells of a surface on a rectilinear block (issue #40 steps 4-5).
!<
!< Faces (step 4). The solid part of every cell face of the block is a 2-D problem: the section of
!< the surface by the grid plane (`slice_surface` with `on_plane_low`: the plane at x_i + eps, a
!< vertex on it counts as below, the rule of the line crossings and of the 2-D cut cells), cut by
!< the 2-D cut-cell machinery on the block of the two other axes. Each grid plane is sliced once,
!< so the two cells that share a face see the same fluid fraction and centroid.
!<
!< In-plane axes follow `slice_surface` (cyclic): planes x = x_i are cut on (y, z), planes y = y_j on
!< (z, x), planes z = z_k on (x, y). The crossing points of the sections are canonical rounded
!< constructions; decisions on them are made on those rounded coordinates, the same in every
!< kernel. A face whose section does not close (an open surface) is flagged.
!<
!< Cells (step 5). Every triangle is clipped to every cell it overlaps (Sutherland-Hodgman with the
!< same rule: a point on the plane x_i belongs to the cell below it, the crossing points canonical).
!< The solid volume and first moments of S n C follow from the divergence theorem over the boundary
!< of S n C: the clipped triangles (outward normal of the solid) and the solid parts of the six cell
!< faces (outward normal of the cell, areas from the faces above): V = 1/3 oint x . n dA,
!< M_x = 1/2 oint x^2 n_x dA, ..., each cell in coordinates local to its lower corner (no loss of
!< accuracy far from the origin). The fluid volume, centroid and the wetted area, centroid and
!< normal integral follow. A cell with a face not closed is flagged and gets no measures.
!<
!< Pieces (issue #40). The connected fluid pieces of every cell, by a union-find over the parts of
!< the boundary of the fluid: the clipped wall polygons (with on-plane flags from the coordinates: a
!< vertex on a cell plane lies on that face, two flags of two axes put it on a cell edge) and the fluid
!< intervals of the 12 cell edges. Since issue #70 the corner states are the exact parities of the
!< perturbed nodes (`edge_fractions`) and the piercings of an edge its exact crossings under the global
!< perturbation (`tri_line_hit_oac`, kept when `line_crossing_before` puts them between the perturbed
!< corners), so the states alternate from one corner to the other; crossings at one rounded position are
!< reordered by that alternation. Before #70 the piercings were the polygon vertices on the edge, one per
!< incident facet at a mesh edge or vertex on it, and lattice bodies flagged up to 38 cells each.
!< Joins: polygons sharing a vertex (the clipped points of a shared mesh edge are bit-identical),
!< a piercing with the fluid interval next to it, the fluid intervals at a fluid corner, and a solid
!< island traced on a face (a loop of wall polygon edges on the face that touches no cell edge, the
!< fluid outside it) with the nearest trace, or the fluid interval of the face edge, a ray in +u
!< from its rightmost point reaches. Every set that touches the cell boundary is the outer boundary
!< of one fluid piece; a closed wall set inside the cell adds a piece only if it bounds a void
!< (negative solid volume). So a passage narrower than the cell (the fluid on the cell boundary in
!< two disks joined inside) is one piece, and a thin plate across the cell gives two. The count needs
!< a closed surface in the cell: a cell that a boundary edge of the mesh (an edge of one facet)
!< overlaps, a cell whose fraction is outside [0, 1] (an open surface corrupts the parity of the
!< faces beyond the cells it flags), a cell without triangles whose fraction is not 0 or 1, and a
!< cell whose edge states do not close get -1. A surface on the high plane of a cell lies inside that cell, so the open face
!< bounds a zero-volume fluid piece there (as in 2-D).
!<
!< Piece measures (issue #40). Per piece: fluid volume and centroid, wetted area, centroid and
!< normal integral, and the fluid fraction of each of the six cell faces (`piece_*`, the pieces of
!< a cell at first_piece .. first_piece + pieces - 1). Walls: the divergence-theorem terms of every
!< clipped polygon go to its set. Faces: Green's theorem per set over the fluid edge intervals and
!< the trace segments, oriented by the polygon winding (fluid on the left on a high face, on the right
!< on a low face). Volume and moments follow from the high faces and the walls. An inclusion (a closed
!< wall set of positive solid volume) goes to the set around it, found by a +x ray. A single-piece
!< cell takes the cell totals; zero-volume pieces are kept, with their volume (0), for a solver to drop.
!< Traces are read on the perturbed face plane (issue #70): a polygon lying in the face plane bounds the
!< fluid only if its solid lies beyond the face, and then not along the cell edges (the edge intervals
!< carry the boundary there); a facet ending on the face plane counts there only if the facet across its
!< edge on the plane crosses it. A multi-piece cell is trusted only if its split closes (`split_closes`:
!< face sums equal the plane sections, volumes add up, every piece closes per axis), else -1: on 1.57
!< million cells of lattice bodies 19 remain (947 before #70), where a facet touches the cell along a
!< segment of a cell edge only, so its trace on the face is missing.
!<
!< Kernels (issue #100). The cells are per-cell kernels on plain arrays (the `!$acc routine seq` form of #20, for the
!< device in WP2 step 3c): `cell_totals_3d` (clipped triangles and totals, a gather over the cell's candidates),
!< `cell_pieces_3d` (pieces and measures, count then fill, on a workspace slice sized by `cell_pieces_sizes` from the
!< cell's candidates and edge piercings), `single_piece_3d`. The candidates are binned by bounding box, then kept only if
!< the box-triangle test says they touch the cell (`cell_candidates_3d`); the exact edge piercings come from the host
!< (`cell_piercings`: exact orient3d). Moments, quadrature and the links' intervals take the same kernels.
use fossil_box_predicates,     only : box_triangle_overlap
use fossil_box_query,          only : box_tolerance
use fossil_curve_object,       only : curve_object, cut_cells_2d_t, cut_cells_2d_blocks_t, cut_cells_device_many
use fossil_curve_slice,        only : slice_surface
use fossil_cut_graph,          only : cut_links_t, build_csr
use fossil_cut_moments,        only : MOMENTS_MAX, mono_count, mono_exponents, triangle_moments, segment_green_moments, &
                                      cell_moments, tet_rule, triangle_rule, tensor_rule
use fossil_facet_object,       only : facet_object
use fossil_line_query,         only : line_crossing_before, tri_line_hit_oac
use fossil_predicates,         only : orient2d_xy
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P

implicit none
private
public :: cut_cells_3d_t, cut_cells_3d, cut_cells_3d_sections

integer(I4P), parameter :: MONO2_MAX = 15_I4P !< mono_count(2, MOMENTS_MAX): moments per face, at most.
integer(I4P), parameter :: MONO3_MAX = 35_I4P !< mono_count(3, MOMENTS_MAX): moments per piece, at most.

type :: cut_cells_3d_t
   !< Cut-cell geometry of a 3-D block (fluid = outside the solid).
   integer(I4P)              :: nx = 0_I4P, ny = 0_I4P, nz = 0_I4P !< Cells.
   real(R8P),    allocatable :: face_x(:,:,:)      !< Fluid fraction of the faces x = x_i [nx+1, ny, nz].
   real(R8P),    allocatable :: face_y(:,:,:)      !< Fluid fraction of the faces y = y_j [nx, ny+1, nz].
   real(R8P),    allocatable :: face_z(:,:,:)      !< Fluid fraction of the faces z = z_k [nx, ny, nz+1].
   real(R8P),    allocatable :: face_x_centroid(:,:,:,:) !< Fluid centroid of the x faces [3, nx+1, ny, nz].
   real(R8P),    allocatable :: face_y_centroid(:,:,:,:) !< Fluid centroid of the y faces [3, nx, ny+1, nz].
   real(R8P),    allocatable :: face_z_centroid(:,:,:,:) !< Fluid centroid of the z faces [3, nx, ny, nz+1].
   integer(I4P), allocatable :: face_x_pieces(:,:,:) !< Fluid pieces of the x faces.
   integer(I4P), allocatable :: face_y_pieces(:,:,:) !< Fluid pieces of the y faces.
   integer(I4P), allocatable :: face_z_pieces(:,:,:) !< Fluid pieces of the z faces.
   logical,      allocatable :: face_x_closed(:,:,:) !< The x face section closed.
   logical,      allocatable :: face_y_closed(:,:,:) !< The y face section closed.
   logical,      allocatable :: face_z_closed(:,:,:) !< The z face section closed.
   real(R8P),    allocatable :: alpha(:,:,:)        !< Fluid volume fraction [nx, ny, nz] (-1 if not closed).
   real(R8P),    allocatable :: centroid(:,:,:,:)   !< Fluid centroid [3, nx, ny, nz].
   real(R8P),    allocatable :: wet_area(:,:,:)     !< Wetted area [nx, ny, nz].
   real(R8P),    allocatable :: wet_centroid(:,:,:,:) !< Wetted centroid [3, nx, ny, nz].
   real(R8P),    allocatable :: wet_normal(:,:,:,:) !< Integral of the outward (solid to fluid) normal over the wall [3, ...].
   logical,      allocatable :: closed(:,:,:)       !< All six faces of the cell closed.
   integer(I4P), allocatable :: pieces(:,:,:)       !< Fluid pieces per cell [nx, ny, nz] (-1 if not closed).
   integer(I4P), allocatable :: first_piece(:,:,:)  !< Index of the first piece of the cell [nx, ny, nz].
   integer(I4P)              :: np = 0_I4P          !< Pieces in the block.
   integer(I4P), allocatable :: piece_cell(:,:)     !< Cell of each piece [3, np].
   real(R8P),    allocatable :: piece_volume(:)     !< Fluid volume of each piece [np].
   real(R8P),    allocatable :: piece_centroid(:,:) !< Fluid centroid of each piece [3, np] (cell centre if no volume).
   real(R8P),    allocatable :: piece_wet_area(:)   !< Wetted area of each piece [np].
   real(R8P),    allocatable :: piece_wet_centroid(:,:) !< Wetted centroid of each piece [3, np].
   real(R8P),    allocatable :: piece_wet_normal(:,:) !< Outward (solid to fluid) normal integral of each piece [3, np].
   real(R8P),    allocatable :: piece_face(:,:)     !< Fluid fraction of the cell faces per piece [6 (-x, +x, -y, +y, -z, +z), np].
   real(R8P),    allocatable :: edge_x(:,:,:)       !< Fluid fraction of the cell edges along x [nx, ny+1, nz+1] (issue #70,
                                                    !< with `edges=.true.`; edge_y, edge_z likewise).
   real(R8P),    allocatable :: edge_y(:,:,:)       !< Fluid fraction of the cell edges along y [nx+1, ny, nz+1].
   real(R8P),    allocatable :: edge_z(:,:,:)       !< Fluid fraction of the cell edges along z [nx+1, ny+1, nz].
   logical,      allocatable :: node_in(:,:,:)      !< Grid node inside the solid [nx+1, ny+1, nz+1] (perturbed, exact).
   integer(I4P)              :: moments_order = 0_I4P !< Order of the moments computed (issue #70; 0: none).
   real(R8P),    allocatable :: piece_moments(:,:)  !< int X^a Y^b Z^c dV of each piece [mono_count(3, p), np].
   real(R8P),    allocatable :: piece_wall_moments(:,:) !< int X^a Y^b Z^c dA over the wall of each piece [mono_count(3, p), np].
   real(R8P),    allocatable :: piece_face_moments(:,:,:) !< Fluid moments of the six faces per piece [mono_count(2, p), 6, np]
                                                    !< (in-face coordinates in cyclic order: x faces (Y, Z), y (Z, X), z (X, Y)).
   real(R8P),    allocatable :: face_x_moments(:,:,:,:) !< Fluid moments of the x faces [mono_count(2, p), nx+1, ny, nz].
   real(R8P),    allocatable :: face_y_moments(:,:,:,:) !< Fluid moments of the y faces [mono_count(2, p), nx, ny+1, nz].
   real(R8P),    allocatable :: face_z_moments(:,:,:,:) !< Fluid moments of the z faces [mono_count(2, p), nx, ny, nz+1].
   integer(I4P), allocatable :: quad_start(:)       !< Volume rule of piece p: quad_start(p) .. quad_start(p + 1) - 1 (#70).
   real(R8P),    allocatable :: quad_point(:,:)     !< Its points [3, nq].
   real(R8P),    allocatable :: quad_weight(:)      !< Its weights (signed: cones over the boundary).
   integer(I4P), allocatable :: quad_status(:)      !< 0: the weights add up to the piece volume (1e-12 of the cell), 1: not.
   integer(I4P), allocatable :: wall_quad_start(:)  !< Wall rule of piece p (positive weights).
   real(R8P),    allocatable :: wall_quad_point(:,:) !< Its points [3, nw].
   real(R8P),    allocatable :: wall_quad_weight(:) !< Its weights.
endtype cut_cells_3d_t

contains
   subroutine cut_cells_3d(surface, x, y, z, cut, status, links, edges, moments_order, quadrature, device_sections)
   !< Cut cells of `surface` on the block with cell edges x(nx+1), y(ny+1), z(nz+1) (issue #40).
   !<
   !< Step 4: the face fluid fractions, centroids and pieces of every grid face (see the module
   !< notes). `status = STATUS_INVALID_INPUT` on bad edges.
   !<
   !< `device_sections = .true.` (issue #100 step 3b; default .false.) cuts the plane sections on the device, batched
   !< per axis (`sections_device`), with the same bits; it applies only without `links`, `moments_order` > 0 and
   !< `quadrature` (else the sections stay on the host).
   type(surface_stl_object), intent(in)            :: surface    !< Surface.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(cut_cells_3d_t),     intent(out)           :: cut        !< Cut cells.
   integer(I4P),             intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   type(cut_links_t),        intent(out), optional :: links      !< Links between the pieces (issue #70).
   logical,                  intent(in),  optional :: edges      !< Also the fluid fractions of the cell edges (issue #70).
   integer(I4P),             intent(in),  optional :: moments_order !< Moments of the pieces up to this order (0-4, #70).
   logical,                  intent(in),  optional :: quadrature !< Quadrature rules on the pieces and their walls (#70).
   logical,                  intent(in),  optional :: device_sections !< The plane sections on the device (#100 step 3b).
   type(cut_cells_2d_t), allocatable :: px(:), py(:), pz(:) !< Plane sections kept for the links.
   logical              :: all_edges !< Edge fractions asked.
   integer(I4P)         :: mo, nm2   !< Moment order, 2-D moments per face.
   logical              :: qd        !< Quadrature asked.
   integer(I4P),         allocatable :: piv_piece(:), piv_edge(:) !< Fluid edge intervals of the pieces of multi-piece
   real(R8P),            allocatable :: piv_lo(:), piv_hi(:)      !< cells: piece, cell edge (1-12), range.
   integer(I4P)         :: nx, ny, nz, i, j, k, st, f, m, c, nf !< Work.
   real(R8P),    allocatable :: tri(:,:,:), nrm(:,:) !< Triangles [3, 3, nf], normals [3, nf].
   integer(I4P), allocatable :: cnt(:), kstart(:), klst(:) !< Candidates of the cells (see `cell_candidates_3d`).
   type(facet_object), pointer :: fp      !< Facet.
   logical              :: dsec   !< The sections on the device.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   if (nx < 1 .or. ny < 1 .or. nz < 1) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny)) .or. any(z(2:) <= z(:nz))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   cut%nx = nx ; cut%ny = ny ; cut%nz = nz
   allocate(cut%face_x(nx + 1, ny, nz), cut%face_y(nx, ny + 1, nz), cut%face_z(nx, ny, nz + 1))
   allocate(cut%face_x_centroid(3, nx + 1, ny, nz), cut%face_y_centroid(3, nx, ny + 1, nz), &
            cut%face_z_centroid(3, nx, ny, nz + 1))
   allocate(cut%face_x_pieces(nx + 1, ny, nz), cut%face_y_pieces(nx, ny + 1, nz), cut%face_z_pieces(nx, ny, nz + 1))
   allocate(cut%face_x_closed(nx + 1, ny, nz), cut%face_y_closed(nx, ny + 1, nz), cut%face_z_closed(nx, ny, nz + 1))
   mo = 0_I4P ; if (present(moments_order)) mo = max(0_I4P, min(MOMENTS_MAX, moments_order))
   qd = .false. ; if (present(quadrature)) qd = quadrature
   if (present(links) .or. qd) allocate(px(nx + 1), py(ny + 1), pz(nz + 1))
   cut%moments_order = mo ; nm2 = 0_I4P ; if (mo > 0_I4P) nm2 = mono_count(2_I4P, mo)
   allocate(cut%face_x_moments(nm2, nx + 1, ny, nz), cut%face_y_moments(nm2, nx, ny + 1, nz), &
            cut%face_z_moments(nm2, nx, ny, nz + 1))
   ! the plane sections: faces of every grid plane
   dsec = .false. ; if (present(device_sections)) dsec = device_sections
   if (dsec .and. .not. (present(links) .or. qd .or. mo > 0_I4P)) then
      call sections_device(surface, x, y, z, cut)
   else
      call sections_host(surface, x, y, z, cut, present(links), mo, qd, nm2, px, py, pz)
   endif
   ! grid nodes (issue #70), and the cell edges if asked: exact crossings along the grid lines
   all_edges = .false. ; if (present(edges)) all_edges = edges ! no short-circuit in Fortran: never edges if absent
   call edge_fractions(surface, x, y, z, cut, all_edges)
   ! cells (issue #100): the candidate facets of every cell (`cell_candidates_3d`: bounding-box binning, then the box
   ! test), then per cell the clipped triangles, the faces and the totals (`cell_totals_3d`). Before #100 a facet loop
   ! scattered the clipped triangles into the cells; the gather adds the same terms to each cell in the same order (its
   ! facets in increasing order, the fan of each in order), so the bits are the same; a candidate the box test drops has
   ! no clipped polygon in the cell (see `cell_candidates_3d`).
   nf = surface%get_facets_number()
   allocate(tri(3, 3, nf), nrm(3, nf))
   do f = 1, nf
      fp => surface%facet_at(f)
      do m = 1, 3
         tri(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
      nrm(:, f) = cross3(tri(:, 2, f) - tri(:, 1, f), tri(:, 3, f) - tri(:, 1, f))
   enddo
   call cell_candidates_3d(x, y, z, nf, tri, nrm, cnt, kstart, klst)
   allocate(cut%alpha(nx, ny, nz), cut%centroid(3, nx, ny, nz), cut%wet_area(nx, ny, nz), cut%wet_centroid(3, nx, ny, nz), &
            cut%wet_normal(3, nx, ny, nz), cut%closed(nx, ny, nz))
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1_I4P) + ny * (k - 1_I4P))
            call cell_totals_3d(nf, tri, kstart(c + 1) - kstart(c), klst(kstart(c)), x(i), x(i + 1), y(j), y(j + 1), &
                                z(k), z(k + 1), cut%face_x(i + 1, j, k), cut%face_y(i, j + 1, k), cut%face_z(i, j, k + 1), &
                                cut%face_x_closed(i, j, k), cut%face_x_closed(i + 1, j, k), cut%face_y_closed(i, j, k), &
                                cut%face_y_closed(i, j + 1, k), cut%face_z_closed(i, j, k), cut%face_z_closed(i, j, k + 1), &
                                cut%alpha(i, j, k), cut%centroid(1, i, j, k), cut%wet_area(i, j, k), &
                                cut%wet_centroid(1, i, j, k), cut%wet_normal(1, i, j, k), cut%closed(i, j, k))
         enddo
      enddo
   enddo
   if (present(links)) then
      call cut_pieces_3d(x, y, z, cut, nf, tri, nrm, cnt, kstart, klst, mo, qd, piv_piece, piv_edge, piv_lo, piv_hi, px, py, pz)
      call links_3d(cut, x, y, z, px, py, pz, piv_piece, piv_edge, piv_lo, piv_hi, links)
   else
      if (qd) then
         call cut_pieces_3d(x, y, z, cut, nf, tri, nrm, cnt, kstart, klst, mo, qd, px=px, py=py, pz=pz)
      else
         call cut_pieces_3d(x, y, z, cut, nf, tri, nrm, cnt, kstart, klst, mo, qd)
      endif
   endif
   endsubroutine cut_cells_3d

   subroutine cell_candidates_3d(x, y, z, nf, tri, nrm, cnt, kstart, klst)
   !< Candidate facets of every cell (issue #100): `cnt(c)` facets whose bounding box overlaps cell c (the binning of
   !< #40, kept for the cells' classification), and among them, in increasing order, the ones that touch the closed cell
   !< box (`klst(kstart(c) : kstart(c + 1) - 1)`). The test is `box_triangle_overlap`, conservative: inflated by `tol`
   !< it never misses a contact. A dropped facet has no clipped polygon in the cell (the clip keeps points of the closed
   !< box only, up to round-off far below tol), no piercing of a cell edge (a crossing between the perturbed corners lies
   !< on the closed edge) and no edge shared with a trace (`across_crosses`: the facet would contain the trace). Its atom
   !< in `cell_pieces_core` would be a singleton never counted, and the root of every set is its least atom, so the
   !< pieces, their order and every sum (over the kept polygons in list order) are those of the full list. The
   !< inflation covers the round-off of the projections over the grid and the surface; a facet whose computed normal is
   !< small next to its edges (|n|_1 < 1e-2 e_max^2: the spread of its vertices along that axis could exceed the
   !< inflation) is kept without the test.
   real(R8P),                 intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),              intent(in)  :: nf               !< Facets.
   real(R8P),                 intent(in)  :: tri(3, 3, nf)    !< Triangles.
   real(R8P),                 intent(in)  :: nrm(3, nf)       !< Their normals (not unit).
   integer(I4P), allocatable, intent(out) :: cnt(:)           !< Bounding-box candidates per cell [ncell + 1].
   integer(I4P), allocatable, intent(out) :: kstart(:)        !< Kept candidates of cell c: kstart(c) .. kstart(c + 1) - 1.
   integer(I4P), allocatable, intent(out) :: klst(:)          !< Kept candidates [their number + 1].
   integer(I4P), allocatable :: start(:), lst(:)              !< Bounding-box candidates (CSR).
   logical,      allocatable :: thin(:)                       !< Facets kept without the box test.
   real(R8P)                 :: tol, emax2                    !< Box test inflation, longest squared edge.
   integer(I4P)              :: nx, ny, nz, f, m, i, j, k, il, ih, jl, jh, kl, kh, c, w !< Work.

   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   allocate(cnt(nx * ny * nz + 1)) ; cnt = 0_I4P
   do m = 1, 2 ! count, then fill
      if (m == 2) then
         allocate(start(nx * ny * nz + 1), lst(sum(cnt)))
         start(1) = 1_I4P
         do c = 1, nx * ny * nz
            start(c + 1) = start(c) + cnt(c)
         enddo
         cnt = 0_I4P
      endif
      do f = 1, nf
         il = max(1_I4P, col_of(minval(tri(1, :, f)), x)) ; ih = min(nx, col_of(maxval(tri(1, :, f)), x))
         jl = max(1_I4P, col_of(minval(tri(2, :, f)), y)) ; jh = min(ny, col_of(maxval(tri(2, :, f)), y))
         kl = max(1_I4P, col_of(minval(tri(3, :, f)), z)) ; kh = min(nz, col_of(maxval(tri(3, :, f)), z))
         do k = kl, kh
            do j = jl, jh
               do i = il, ih
                  c = i + nx * ((j - 1) + ny * (k - 1))
                  if (m == 2) lst(start(c) + cnt(c)) = f
                  cnt(c) = cnt(c) + 1_I4P
               enddo
            enddo
         enddo
      enddo
   enddo
   allocate(thin(nf))
   tol = box_tolerance([min(x(1), minval(tri(1, :, :))), min(y(1), minval(tri(2, :, :))), min(z(1), minval(tri(3, :, :)))], &
                       [max(x(nx + 1), maxval(tri(1, :, :))), max(y(ny + 1), maxval(tri(2, :, :))), &
                        max(z(nz + 1), maxval(tri(3, :, :)))])
   do f = 1, nf
      emax2 = max(sum((tri(:, 2, f) - tri(:, 1, f)) ** 2), sum((tri(:, 3, f) - tri(:, 2, f)) ** 2), &
                  sum((tri(:, 1, f) - tri(:, 3, f)) ** 2))
      thin(f) = sum(abs(nrm(:, f))) < 1.e-2_R8P * emax2
   enddo
   allocate(kstart(nx * ny * nz + 1), klst(size(lst) + 1)) ; klst = 0_I4P
   w = 0_I4P
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1) + ny * (k - 1))
            kstart(c) = w + 1_I4P
            do m = start(c), start(c + 1) - 1_I4P
               f = lst(m)
               if (thin(f) .or. box_triangle_overlap([x(i), y(j), z(k)], [x(i + 1), y(j + 1), z(k + 1)], tri(:, 1, f), &
                                                     tri(:, 2, f), tri(:, 3, f), tol)) then
                  w = w + 1_I4P ; klst(w) = f
               endif
            enddo
         enddo
      enddo
   enddo
   kstart(nx * ny * nz + 1) = w + 1_I4P
   endsubroutine cell_candidates_3d

   subroutine cell_totals_3d(nf, tri, nt, lst, x0, x1, y0, y1, z0, z1, fx1, fy1, fz1, cx0, cx1, cy0, cy1, cz0, cz1, alpha, &
                             centroid, wet_area, wet_centroid, wet_normal, closed)
   !< Kernel (issue #100): the cell [x0, x1] x [y0, y1] x [z0, z1] from its facets lst (increasing): the clipped
   !< triangles (Sutherland-Hodgman, a point on a plane belongs to the cell below it), their divergence-theorem terms in
   !< coordinates local to the low corner, then with the solid parts of the +x, +y, +z faces (the low faces are at local
   !< 0) the solid volume and moments; the fluid fraction and centroid, the wetted area, centroid and normal integral.
   !< A cell with a face not closed gets alpha = -1 (see the module notes).
   !$acc routine seq
   integer(I4P), intent(in)  :: nf                     !< Facets.
   real(R8P),    intent(in)  :: tri(3, 3, nf)          !< Triangles.
   integer(I4P), intent(in)  :: nt                     !< Facets of the cell.
   integer(I4P), intent(in)  :: lst(nt)                !< Their indexes (increasing).
   real(R8P),    intent(in)  :: x0, x1, y0, y1, z0, z1 !< Cell.
   real(R8P),    intent(in)  :: fx1, fy1, fz1          !< Fluid fractions of the faces +x, +y, +z.
   logical,      intent(in)  :: cx0, cx1, cy0, cy1, cz0, cz1 !< The faces closed.
   real(R8P),    intent(out) :: alpha                  !< Fluid fraction (-1 if not closed).
   real(R8P),    intent(out) :: centroid(3)            !< Fluid centroid.
   real(R8P),    intent(out) :: wet_area               !< Wetted area.
   real(R8P),    intent(out) :: wet_centroid(3)        !< Wetted centroid.
   real(R8P),    intent(out) :: wet_normal(3)          !< Wall normal integral.
   logical,      intent(out) :: closed                 !< All six faces closed.
   real(R8P)                 :: acc(11), poly(3, 16), c0(3), h(3), sq(3), vc, vs, ms(3), d1(3), d2(3), d3(3), s !< Work.
   integer(I4P)              :: q, np, m, a            !< Work.

   acc = 0._R8P
   c0(1) = x0 ; c0(2) = y0 ; c0(3) = z0
   do q = 1, nt
      poly(:, 1:3) = tri(:, :, lst(q)) ; np = 3_I4P
      call clip_poly(poly, np, 1_I4P, x0, .true.) ; call clip_poly(poly, np, 1_I4P, x1, .false.)
      call clip_poly(poly, np, 2_I4P, y0, .true.) ; call clip_poly(poly, np, 2_I4P, y1, .false.)
      call clip_poly(poly, np, 3_I4P, z0, .true.) ; call clip_poly(poly, np, 3_I4P, z1, .false.)
      if (np < 3_I4P) cycle
      do m = 2, np - 1 ! fan, local coordinates
         d1 = poly(:, 1) - c0 ; d2 = poly(:, m) - c0 ; d3 = poly(:, m + 1) - c0
         call tri_terms(d1, d2, d3, acc)
      enddo
   enddo
   h(1) = x1 - x0 ; h(2) = y1 - y0 ; h(3) = z1 - z0
   vc = h(1) * h(2) * h(3)
   closed = cx0 .and. cx1 .and. cy0 .and. cy1 .and. cz0 .and. cz1
   ! +x, +y, +z faces (local coordinate h); the -x, -y, -z faces are at local 0 and add nothing
   sq(1) = (1._R8P - fx1) * h(2) * h(3) ; sq(2) = (1._R8P - fy1) * h(3) * h(1) ; sq(3) = (1._R8P - fz1) * h(1) * h(2)
   s = 0._R8P
   do a = 1, 3 ! sum(h * sq), in its order
      s = s + h(a) * sq(a)
   enddo
   vs = acc(1) + s / 3._R8P
   do a = 1, 3
      ms(a) = acc(1 + a) + 0.5_R8P * h(a) * h(a) * sq(a)
   enddo
   wet_area = acc(5)
   do a = 1, 3
      wet_normal(a) = acc(8 + a)
      wet_centroid(a) = c0(a) + 0.5_R8P * h(a)
      if (acc(5) > 0._R8P) wet_centroid(a) = c0(a) + acc(5 + a) / acc(5)
      centroid(a) = c0(a) + 0.5_R8P * h(a)
   enddo
   if (.not. closed) then
      alpha = -1._R8P
      return
   endif
   alpha = (vc - vs) / vc
   if (vc - vs > 0._R8P) then
      do a = 1, 3
         centroid(a) = c0(a) + (vc * 0.5_R8P * h(a) - ms(a)) / (vc - vs)
      enddo
   endif
   endsubroutine cell_totals_3d

   subroutine cut_pieces_3d(x, y, z, cut, nf, tri, nrm, cnt, kstart, klst, mo, qd, piv_piece, piv_edge, piv_lo, piv_hi, &
                            px, py, pz)
   !< Fluid pieces of every cell (see the module notes); with piv_*, the fluid edge intervals of the pieces of every
   !< multi-piece cell (issue #70, for the links).
   !<
   !< Count, then fill (issue #100). Every cell is classified: no count (-1), a cell without candidates (fluid or solid
   !< by its fraction), a cell whose candidates all miss it (the box test) with equal corner states (one piece or none,
   !< what `cell_pieces_core` gives on an empty list), or a cut cell. A cut cell gets its exact edge piercings
   !< (`cell_piercings`, host: exact orient3d), a workspace sized from its candidates and piercings
   !< (`cell_pieces_sizes`) and the count pass of `cell_pieces_3d`: its pieces and the numbers of its intervals and rule
   !< points. Prefix sums over the cells in order give every cell its first piece, interval and rule points; the fill
   !< writes the outputs in place: the cell totals for a single piece, the fill pass of `cell_pieces_3d` (the same
   !< kernel on the same inputs, the same bits) where the kernel's pieces are needed (several pieces, or one with moments
   !< or quadrature), and the plane-section cones of the quadrature. No buffer grows; the outputs are those of the
   !< cell-by-cell loop with growing buffers that this replaces, bit for bit and in the same order.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(inout) :: cut        !< Cut cells (alpha, closed set).
   integer(I4P),             intent(in)    :: nf         !< Facets.
   real(R8P),                intent(in)    :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P),             intent(in)    :: cnt(*), kstart(*), klst(*) !< Candidates (see `cell_candidates_3d`).
   integer(I4P),             intent(in)    :: mo         !< Moment order (0: none).
   logical,                  intent(in)    :: qd         !< Quadrature rules.
   integer(I4P), allocatable, intent(out), optional :: piv_piece(:), piv_edge(:) !< Intervals: piece, cell edge.
   real(R8P),    allocatable, intent(out), optional :: piv_lo(:), piv_hi(:)      !< Intervals: range.
   type(cut_cells_2d_t),     intent(in),  optional :: px(:), py(:), pz(:) !< Plane sections (with boundaries, for qd).
   integer(I4P), parameter     :: K_NONE = 0_I4P, K_EMPTY = 1_I4P, K_QUIET = 2_I4P, K_CUT = 3_I4P !< Cell kinds.
   integer(I4P), allocatable   :: kind_(:)                !< Kind of each cell.
   integer(I4P), allocatable   :: civ(:), cq(:), cw(:), cq1(:), cns(:) !< Per cell: intervals, rule points, wall rule
                                                          !< points, wall cones of a single piece, section cones.
   integer(I4P), allocatable   :: ivo(:), qo(:), wo(:)    !< Intervals, rule points, wall rule points before each cell.
   integer(I4P)                :: nm2, nm3                !< Moments per face, per piece.
   integer(I4P), allocatable   :: e2(:,:), e3(:,:)        !< Their exponents.
   real(R8P),    allocatable   :: secm(:,:)               !< Section moments of the six faces of a cell [nm2, 6].
   real(R8P),    allocatable   :: pmom(:,:), pwmom(:,:), pfmom(:,:,:) !< Piece moments: volume, wall, faces.
   integer(I4P), allocatable   :: ivp(:), ive(:)          !< Intervals of the multi-piece cells: piece, cell edge.
   real(R8P),    allocatable   :: ivl(:), ivh(:)          !< Their ranges.
   real(R8P),    allocatable   :: qpt(:,:), qwt(:), wpt(:,:), wwt(:) !< Rules: volume points, weights; wall.
   integer(I4P), allocatable   :: qcnt(:), wcnt(:)        !< Rule points per piece.
   integer(I4P), allocatable   :: pce(:), pcq(:), pcen(:) !< Piercings of a cell: edge, candidate, entering (1).
   real(R8P),    allocatable   :: pct(:)                  !< Their positions.
   real(R8P),    allocatable   :: rw(:)                   !< Real workspace of a cell.
   integer(I4P), allocatable   :: iw(:)                   !< Integer workspace of a cell.
   integer(I4P)                :: nx, ny, nz, m, i, j, k, il, ih, jl, jh, kl, kh, c, np, nt, npc, nrw, niw, q, kd !< Work.
   integer(I4P)                :: ncell, npt, nivt, nqt, nwt, nivc, nqc, nwc, qwn1, p0, m3, ns, f_, ntmax !< Work.
   real(R8P)                   :: c0(3), c1(3)            !< Cell corners.
   logical,      allocatable   :: open_cell(:,:,:)        !< Cells a mesh boundary edge overlaps.
   real(R8P),    allocatable   :: ek(:,:)                 !< Edge keys [6, 3 nf] (end points, lexicographic order).
   integer(I4P), allocatable   :: idx(:)                  !< Sorted edges.
   real(R8P)                   :: pa(3), pb(3)            !< Edge ends.
   real(R8P)                   :: h(3), tpt(3, 27), tw(27), faces(6), zero(MONO3_MAX) !< Work.
   integer(I4P)                :: cs(8)                   !< Corner states.
   logical                     :: want_iv                 !< Intervals asked.

   nx = cut%nx ; ny = cut%ny ; nz = cut%nz ; ncell = nx * ny * nz
   want_iv = present(piv_piece)
   allocate(cut%pieces(nx, ny, nz), cut%first_piece(nx, ny, nz))
   nm2 = 0_I4P ; nm3 = 0_I4P
   if (mo > 0_I4P) then
      nm2 = mono_count(2_I4P, mo) ; nm3 = mono_count(3_I4P, mo) ; e2 = mono_exponents(2_I4P, mo) ; e3 = mono_exponents(3_I4P, mo)
   else
      allocate(e2(2, 0), e3(3, 0))
   endif
   allocate(secm(nm2, 6)) ; zero = 0._R8P
   ! mesh boundary edges (one facet): the cells they overlap get no count
   allocate(open_cell(nx, ny, nz), ek(6, 3 * nf), idx(3 * nf)) ; open_cell = .false.
   do f_ = 1, nf
      do m = 1, 3
         pa = tri(:, m, f_) ; pb = tri(:, mod(m, 3) + 1, f_)
         if (lexless(pb, pa)) then
            ek(:, 3 * (f_ - 1) + m) = [pb, pa]
         else
            ek(:, 3 * (f_ - 1) + m) = [pa, pb]
         endif
      enddo
   enddo
   call sort_keys(6_I4P, 3_I4P * nf, ek, idx)
   m = 1
   do while (m <= 3 * nf)
      f_ = m
      do while (f_ < 3 * nf)
         if (any(ek(:, idx(f_ + 1)) /= ek(:, idx(m)))) exit
         f_ = f_ + 1
      enddo
      if (f_ == m) then ! an edge of one facet
         pa = ek(1:3, idx(m)) ; pb = ek(4:6, idx(m))
         il = max(1_I4P, col_of(min(pa(1), pb(1)), x)) ; ih = min(nx, col_of(max(pa(1), pb(1)), x))
         jl = max(1_I4P, col_of(min(pa(2), pb(2)), y)) ; jh = min(ny, col_of(max(pa(2), pb(2)), y))
         kl = max(1_I4P, col_of(min(pa(3), pb(3)), z)) ; kh = min(nz, col_of(max(pa(3), pb(3)), z))
         if (il <= ih .and. jl <= jh .and. kl <= kh) open_cell(il:ih, jl:jh, kl:kh) = .true.
      endif
      m = f_ + 1
   enddo
   ! 1. classify and count; corner states: the exact inside flags of the perturbed nodes (issue #70, `edge_fractions`)
   allocate(kind_(ncell), civ(ncell), cq(ncell), cw(ncell), cq1(ncell), cns(ncell))
   kind_ = K_NONE ; civ = 0_I4P ; cq = 0_I4P ; cw = 0_I4P ; cq1 = 0_I4P ; cns = 0_I4P
   ntmax = 1_I4P
   do c = 1, ncell
      ntmax = max(ntmax, kstart(c + 1) - kstart(c))
   enddo
   allocate(pce(12 * ntmax), pcq(12 * ntmax), pcen(12 * ntmax), pct(12 * ntmax)) ! at most 12 nt piercings
   ! zero-size outputs for the count pass (it writes none)
   allocate(rw(1), iw(1), pmom(nm3, 0), pwmom(nm3, 0), pfmom(nm2, 6, 0), ivp(0), ive(0), ivl(0), ivh(0), qpt(3, 0), qwt(0), &
            wpt(3, 0), wwt(0), qcnt(0), wcnt(0))
   allocate(cut%piece_volume(0), cut%piece_centroid(3, 0), cut%piece_wet_area(0), cut%piece_wet_centroid(3, 0), &
            cut%piece_wet_normal(3, 0), cut%piece_face(6, 0))
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1) + ny * (k - 1))
            if (.not. cut%closed(i, j, k) .or. open_cell(i, j, k) .or. cut%alpha(i, j, k) < -1.e-12_R8P .or. &
                cut%alpha(i, j, k) > 1._R8P + 1.e-12_R8P) then
               cut%pieces(i, j, k) = -1_I4P
            elseif (cnt(c) == 0_I4P .and. abs(cut%alpha(i, j, k) - anint(cut%alpha(i, j, k))) > 1.e-12_R8P) then
               cut%pieces(i, j, k) = -1_I4P
            elseif (cnt(c) == 0_I4P) then
               cut%pieces(i, j, k) = merge(1_I4P, 0_I4P, cut%alpha(i, j, k) > 0.5_R8P) ; kind_(c) = K_EMPTY
            else
               do q = 1, 8
                  cs(q) = merge(0_I4P, 1_I4P, cut%node_in(i + mod(q - 1, 2), j + mod((q - 1) / 2, 2), k + (q - 1) / 4))
               enddo
               nt = kstart(c + 1) - kstart(c)
               if (nt == 0_I4P .and. all(cs == cs(1))) then ! what the kernel gives on an empty list
                  cut%pieces(i, j, k) = cs(1) ; kind_(c) = K_QUIET
               else
                  kind_(c) = K_CUT ; c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
                  call cell_piercings(nf, tri, nt, klst(kstart(c):kstart(c + 1) - 1), c0, c1, npc, pce, pcq, pct, pcen)
                  call cell_pieces_sizes(nt, npc, nm2, nm3, qd, nrw, niw)
                  if (nrw > size(rw)) then
                     deallocate(rw) ; allocate(rw(nrw))
                  endif
                  if (niw > size(iw)) then
                     deallocate(iw) ; allocate(iw(niw))
                  endif
                  faces = [cut%face_x(i, j, k), cut%face_x(i + 1, j, k), cut%face_y(i, j, k), cut%face_y(i, j + 1, k), &
                           cut%face_z(i, j, k), cut%face_z(i, j, k + 1)]
                  call cell_pieces_3d(.false., nf, tri, nrm, nt, klst(kstart(c)), x(i), x(i + 1), y(j), y(j + 1), z(k), &
                                      z(k + 1), cs, cut%alpha(i, j, k), faces, npc, pce, pcq, pct, pcen, nm2, nm3, e2, e3, &
                                      qd, want_iv, secm, nrw, rw, niw, iw, np, nivc, nqc, nwc, qwn1, &
                                      0_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P, 0_I4P, cut%piece_volume, &
                                      cut%piece_centroid, cut%piece_wet_area, cut%piece_wet_centroid, cut%piece_wet_normal, &
                                      cut%piece_face, pmom, pwmom, pfmom, ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
                  cut%pieces(i, j, k) = np
                  if (np > 1_I4P) then
                     if (want_iv) civ(c) = nivc
                     cq(c) = nqc ; cw(c) = nwc
                  elseif (np == 1_I4P) then
                     cq1(c) = qwn1 ; cw(c) = nwc
                  endif
               endif
            endif
            np = cut%pieces(i, j, k)
            if (np < 1_I4P .or. .not. qd) cycle
            ! section cones of a single piece; the rule points of the cell
            if (np == 1_I4P .and. kind_(c) /= K_EMPTY) cns(c) = section_cones_count(px, py, pz, i, j, k)
            select case (kind_(c))
            case (K_EMPTY) ; cq(c) = 27_I4P
            case (K_QUIET) ; cq(c) = cns(c)
            case default
               if (np == 1_I4P) cq(c) = cq1(c) + cns(c)
            endselect
         enddo
      enddo
   enddo
   ! 2. prefix sums over the cells in order
   allocate(ivo(ncell), qo(ncell), wo(ncell))
   npt = 0_I4P ; nivt = 0_I4P ; nqt = 0_I4P ; nwt = 0_I4P
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1) + ny * (k - 1))
            cut%first_piece(i, j, k) = npt + 1_I4P ; ivo(c) = nivt ; qo(c) = nqt ; wo(c) = nwt
            npt = npt + max(0_I4P, cut%pieces(i, j, k))
            nivt = nivt + civ(c) ; nqt = nqt + cq(c) ; nwt = nwt + cw(c)
         enddo
      enddo
   enddo
   ! 3. fill
   cut%np = npt
   deallocate(cut%piece_volume, cut%piece_centroid, cut%piece_wet_area, cut%piece_wet_centroid, cut%piece_wet_normal, &
              cut%piece_face)
   allocate(cut%piece_cell(3, npt), cut%piece_volume(npt), cut%piece_centroid(3, npt), cut%piece_wet_area(npt), &
            cut%piece_wet_centroid(3, npt), cut%piece_wet_normal(3, npt), cut%piece_face(6, npt))
   deallocate(pmom, pwmom, pfmom, ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
   allocate(pmom(nm3, npt), pwmom(nm3, npt), pfmom(nm2, 6, npt), ivp(nivt), ive(nivt), ivl(nivt), ivh(nivt))
   if (qd) then
      allocate(qpt(3, nqt), qwt(nqt), wpt(3, nwt), wwt(nwt), qcnt(npt), wcnt(npt))
   else
      allocate(qpt(3, 0), qwt(0), wpt(3, 0), wwt(0), qcnt(0), wcnt(0))
   endif
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1) + ny * (k - 1))
            np = cut%pieces(i, j, k)
            if (np < 1_I4P) cycle
            p0 = cut%first_piece(i, j, k) - 1_I4P ; kd = kind_(c)
            do m = 1, np
               cut%piece_cell(1, p0 + m) = i ; cut%piece_cell(2, p0 + m) = j ; cut%piece_cell(3, p0 + m) = k
            enddo
            h = [x(i + 1) - x(i), y(j + 1) - y(j), z(k + 1) - z(k)]
            if (np == 1_I4P) & ! measures: the cell totals for a single piece (exact), else the per-piece ones (kernel)
               call single_piece_3d(h(1), h(2), h(3), cut%alpha(i, j, k), cut%centroid(1, i, j, k), cut%wet_area(i, j, k), &
                                    cut%wet_centroid(1, i, j, k), cut%wet_normal(1, i, j, k), cut%face_x(i, j, k), &
                                    cut%face_x(i + 1, j, k), cut%face_y(i, j, k), cut%face_y(i, j + 1, k), &
                                    cut%face_z(i, j, k), cut%face_z(i, j, k + 1), cut%piece_volume(p0 + 1), &
                                    cut%piece_centroid(1, p0 + 1), cut%piece_wet_area(p0 + 1), &
                                    cut%piece_wet_centroid(1, p0 + 1), cut%piece_wet_normal(1, p0 + 1), &
                                    cut%piece_face(1, p0 + 1))
            ! moments (issue #70): faces from the plane sections for a single piece, else per piece; volumes by the
            ! divergence theorem along x from the wall terms and the two x faces
            if (nm3 > 0_I4P) then
               if (kd == K_EMPTY) then ! a fluid cell
                  pmom(:, p0 + 1) = cell_moments(e3, product(h)) ; pwmom(:, p0 + 1) = 0._R8P
                  pfmom(:, 1, p0 + 1) = cell_moments(e2, h(2) * h(3)) ; pfmom(:, 2, p0 + 1) = pfmom(:, 1, p0 + 1)
                  pfmom(:, 3, p0 + 1) = cell_moments(e2, h(3) * h(1)) ; pfmom(:, 4, p0 + 1) = pfmom(:, 3, p0 + 1)
                  pfmom(:, 5, p0 + 1) = cell_moments(e2, h(1) * h(2)) ; pfmom(:, 6, p0 + 1) = pfmom(:, 5, p0 + 1)
               elseif (np == 1_I4P) then
                  call section_face_moments(cut, i, j, k, nm2, secm)
                  pfmom(:, :, p0 + 1) = secm
                  if (kd == K_QUIET) then ! no wall
                     pwmom(:, p0 + 1) = 0._R8P
                     call volume_moments_k(nm3, nm2, zero, secm(:, 1), secm(:, 2), h(1), e2, e3, pmom(:, p0 + 1))
                  endif
               endif
            endif
            ! quadrature (issue #70): a fluid cell takes the tensor rule, a cut cell its cones
            if (qd) then
               if (kd == K_EMPTY) then
                  call tensor_rule([x(i), y(j), z(k)], [x(i + 1), y(j + 1), z(k + 1)], tpt, tw)
                  qpt(:, qo(c) + 1:qo(c) + 27) = tpt ; qwt(qo(c) + 1:qo(c) + 27) = tw
                  qcnt(p0 + 1) = 27_I4P ; wcnt(p0 + 1) = 0_I4P
               elseif (np == 1_I4P) then ! the wall cones (kernel), then cones over the exact section faces
                  ns = qo(c) + cq1(c)
                  do m3 = 1, 6
                     select case ((m3 + 1) / 2)
                     case (1) ; call face_cones(px(i + mod(m3 + 1, 2)), x, y, z, [i, j, k], m3, qpt, qwt, ns)
                     case (2) ; call face_cones(py(j + mod(m3 + 1, 2)), x, y, z, [i, j, k], m3, qpt, qwt, ns)
                     case default ; call face_cones(pz(k + mod(m3 + 1, 2)), x, y, z, [i, j, k], m3, qpt, qwt, ns)
                     endselect
                  enddo
                  qcnt(p0 + 1) = cq1(c) + cns(c) ; wcnt(p0 + 1) = cw(c)
               endif
            endif
            ! the kernel's pieces: several pieces, or one with moments or quadrature
            if (kd == K_CUT .and. (np > 1_I4P .or. nm3 > 0_I4P .or. qd)) then
               do q = 1, 8
                  cs(q) = merge(0_I4P, 1_I4P, cut%node_in(i + mod(q - 1, 2), j + mod((q - 1) / 2, 2), k + (q - 1) / 4))
               enddo
               nt = kstart(c + 1) - kstart(c) ; c0 = [x(i), y(j), z(k)] ; c1 = [x(i + 1), y(j + 1), z(k + 1)]
               call cell_piercings(nf, tri, nt, klst(kstart(c):kstart(c + 1) - 1), c0, c1, npc, pce, pcq, pct, pcen)
               call cell_pieces_sizes(nt, npc, nm2, nm3, qd, nrw, niw)
               if (nrw > size(rw)) then
                  deallocate(rw) ; allocate(rw(nrw))
               endif
               if (niw > size(iw)) then
                  deallocate(iw) ; allocate(iw(niw))
               endif
               faces = [cut%face_x(i, j, k), cut%face_x(i + 1, j, k), cut%face_y(i, j, k), cut%face_y(i, j + 1, k), &
                        cut%face_z(i, j, k), cut%face_z(i, j, k + 1)]
               call cell_pieces_3d(.true., nf, tri, nrm, nt, klst(kstart(c)), x(i), x(i + 1), y(j), y(j + 1), z(k), &
                                   z(k + 1), cs, cut%alpha(i, j, k), faces, npc, pce, pcq, pct, pcen, nm2, nm3, e2, e3, &
                                   qd, want_iv, secm, nrw, rw, niw, iw, np, nivc, nqc, nwc, qwn1, &
                                   p0, ivo(c), qo(c), wo(c), npt, nivt, nqt, nwt, size(qcnt, kind=I4P), cut%piece_volume, &
                                   cut%piece_centroid, cut%piece_wet_area, cut%piece_wet_centroid, cut%piece_wet_normal, &
                                   cut%piece_face, pmom, pwmom, pfmom, ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
            endif
         enddo
      enddo
   enddo
   if (nm3 > 0_I4P) then
      call move_alloc(pmom, cut%piece_moments) ; call move_alloc(pwmom, cut%piece_wall_moments)
      call move_alloc(pfmom, cut%piece_face_moments)
   endif
   if (qd) then
      allocate(cut%quad_start(npt + 1), cut%wall_quad_start(npt + 1), cut%quad_status(npt))
      cut%quad_start(1) = 1_I4P ; cut%wall_quad_start(1) = 1_I4P
      do m3 = 1, npt
         cut%quad_start(m3 + 1) = cut%quad_start(m3) + qcnt(m3) ; cut%wall_quad_start(m3 + 1) = cut%wall_quad_start(m3) + wcnt(m3)
      enddo
      call move_alloc(qpt, cut%quad_point) ; call move_alloc(qwt, cut%quad_weight)
      call move_alloc(wpt, cut%wall_quad_point) ; call move_alloc(wwt, cut%wall_quad_weight)
      do m3 = 1, npt
         h = [x(cut%piece_cell(1, m3) + 1) - x(cut%piece_cell(1, m3)), y(cut%piece_cell(2, m3) + 1) - y(cut%piece_cell(2, m3)), &
              z(cut%piece_cell(3, m3) + 1) - z(cut%piece_cell(3, m3))]
         cut%quad_status(m3) = merge(0_I4P, 1_I4P, abs(sum(cut%quad_weight(cut%quad_start(m3):cut%quad_start(m3 + 1) - 1)) - &
                                                     cut%piece_volume(m3)) <= 1.e-12_R8P * product(h))
      enddo
   endif
   if (want_iv) then
      call move_alloc(ivp, piv_piece) ; call move_alloc(ive, piv_edge) ; call move_alloc(ivl, piv_lo)
      call move_alloc(ivh, piv_hi)
   endif
   endsubroutine cut_pieces_3d

   subroutine cell_piercings(nf, tri, nt, lst, c0, c1, npc, pce, pcq, pct, pcen)
   !< Exact crossings of the 12 edges of the cell [c0, c1] by its candidates lst (issue #70), host only (exact orient3d):
   !< a crossing belongs to the edge when it lies after the perturbed low corner and before the perturbed high one
   !< (`line_crossing_before`), so the corner states and the crossings of every edge agree. Raw, in the order edge, then
   !< candidate: edge, candidate (position in lst), the kernel's position along the edge, entering (1) or leaving (0);
   !< `cell_pieces_core` clamps the position to the edge or takes the clipped polygon's vertex on it (issue #100: the
   !< exact tests before the kernel, which does the rest).
   integer(I4P), intent(in)    :: nf            !< Facets.
   real(R8P),    intent(in)    :: tri(3, 3, nf) !< Triangles.
   integer(I4P), intent(in)    :: nt            !< Candidates of the cell.
   integer(I4P), intent(in)    :: lst(nt)       !< Their indexes.
   real(R8P),    intent(in)    :: c0(3), c1(3)  !< Cell corners.
   integer(I4P), intent(out)   :: npc           !< Piercings (at most 12 nt).
   integer(I4P), intent(inout) :: pce(:), pcq(:), pcen(:) !< Edge, candidate, entering.
   real(R8P),    intent(inout) :: pct(:)        !< Position.
   integer(I4P)                :: e, a, kk, q   !< Work.
   real(R8P)                   :: la, lb, tc    !< Edge line, kernel parameter.
   logical                     :: hit, en       !< Crossing, entering.

   npc = 0_I4P
   do e = 1, 12
      a = (e - 1_I4P) / 4_I4P + 1_I4P ; kk = mod(e - 1_I4P, 4_I4P) ! edge axis, sides
      la = merge(c1(mod(a, 3_I4P) + 1_I4P), c0(mod(a, 3_I4P) + 1_I4P), mod(kk, 2_I4P) == 1_I4P)
      lb = merge(c1(mod(a + 1_I4P, 3_I4P) + 1_I4P), c0(mod(a + 1_I4P, 3_I4P) + 1_I4P), kk / 2_I4P == 1_I4P)
      do q = 1, nt
         call tri_line_hit_oac(reshape(tri(:, :, lst(q)), [9]), a, la, lb, hit, tc, en)
         if (.not. hit) cycle
         if (line_crossing_before(tri(:, :, lst(q)), a, la, lb, c0(a), en)) cycle
         if (.not. line_crossing_before(tri(:, :, lst(q)), a, la, lb, c1(a), en)) cycle
         npc = npc + 1_I4P
         pce(npc) = e ; pcq(npc) = q ; pct(npc) = tc ; pcen(npc) = merge(1_I4P, 0_I4P, en)
      enddo
   enddo
   endsubroutine cell_piercings

   pure subroutine cell_pieces_sizes(nt, npc, nm2, nm3, qd, nrw, niw)
   !< Workspace of `cell_pieces_3d` for a cell with nt candidates, npc piercings, nm2 / nm3 moments per face / piece and
   !< the quadrature flag qd (issue #100). With the atoms NA = nt + 8 + npc + 12 (polygons, corners, at most npc + 12
   !< fluid edge intervals NI), the polygon vertices MM = 16 nt and, with qd, the face segments NS = 48 nt + 2 NI (a
   !< polygon edge lies on at most three faces, an interval on two):
   !<    nrw = 48 nt + npc + 2 NI + NA (35 + 2 nm3 + 6 nm2) + 4 MM + 6 NS + 1 reals,
   !<    niw = 17 nt + 3 npc + 2 NI + 4 NA + 3 MM + 2 NS + 2 + (qd: 5 NA) + 1 integers
   !< (layout in `cell_pieces_layout`; the + 1 keeps every slice start inside the slice when the last parts are empty).
   !< Every count of the cell is bounded by its part, so nothing can overflow.
   !$acc routine seq
   integer(I4P), intent(in)  :: nt, npc, nm2, nm3 !< Candidates, piercings, moments per face, per piece.
   logical,      intent(in)  :: qd                !< Quadrature.
   integer(I4P), intent(out) :: nrw, niw          !< Reals, integers.
   integer(I4P)              :: ro(14), io(20)    !< Layout.

   call cell_pieces_layout(nt, npc, nm2, nm3, qd, ro, io)
   nrw = ro(14) ; niw = io(20)
   endsubroutine cell_pieces_sizes

   pure subroutine cell_pieces_layout(nt, npc, nm2, nm3, qd, ro, io)
   !< Start of each part of a cell's workspace (sizes in `cell_pieces_sizes`; ro(14), io(20): the sizes).
   !< Reals: polygons pv(3, 16, nt), piercing positions(npc), interval ranges lo(NI), hi(NI), solid volumes(NA), sort
   !< keys and face segments (4, MM), wall terms(11, NA), face areas(6, NA), measures(17, NA), wall divergence moments
   !< (nm3, NA), wall moments(nm3, NA), face moments(nm2, 6, NA), oriented face segments(6, NS).
   !< Integers: vertex flags(16, nt), vertex counts(nt), piercing edge, polygon, entering(npc), union-find parents(NA),
   !< interval edge(NI), atom(NI), touch flags(NA), inclusion targets(NA), piece ids(NA), sort and segment work(3, MM),
   !< segment face(NS), set(NS), rule starts(NA + 1), wall rule starts(NA + 1), wall cones(NA), fill pointers(NA, NA).
   !$acc routine seq
   integer(I4P), intent(in)  :: nt, npc, nm2, nm3 !< Candidates, piercings, moments per face, per piece.
   logical,      intent(in)  :: qd                !< Quadrature.
   integer(I4P), intent(out) :: ro(14), io(20)    !< Starts (the last: total size).
   integer(I4P)              :: na, ni, mm, ns, nr !< Atoms, intervals, vertices, segments, rule parts.

   na = nt + 20_I4P + npc ; ni = npc + 12_I4P ; mm = 16_I4P * nt
   ns = 0_I4P ; nr = 0_I4P
   if (qd) then
      ns = 48_I4P * nt + 2_I4P * ni ; nr = na
   endif
   ro(1) = 1_I4P
   ro(2) = ro(1) + 48_I4P * nt ; ro(3) = ro(2) + npc ; ro(4) = ro(3) + ni ; ro(5) = ro(4) + ni ; ro(6) = ro(5) + na
   ro(7) = ro(6) + 4_I4P * mm ; ro(8) = ro(7) + 11_I4P * na ; ro(9) = ro(8) + 6_I4P * na ; ro(10) = ro(9) + 17_I4P * na
   ro(11) = ro(10) + nm3 * na ; ro(12) = ro(11) + nm3 * na ; ro(13) = ro(12) + 6_I4P * nm2 * na
   ro(14) = ro(13) + 6_I4P * ns + 1_I4P
   io(1) = 1_I4P
   io(2) = io(1) + 16_I4P * nt ; io(3) = io(2) + nt ; io(4) = io(3) + npc ; io(5) = io(4) + npc ; io(6) = io(5) + npc
   io(7) = io(6) + na ; io(8) = io(7) + ni ; io(9) = io(8) + ni ; io(10) = io(9) + na ; io(11) = io(10) + na
   io(12) = io(11) + na ; io(13) = io(12) + 3_I4P * mm ; io(14) = io(13) + ns ; io(15) = io(14) + ns
   io(16) = io(15) + nr + 1_I4P ; io(17) = io(16) + nr + 1_I4P
   io(18) = io(17) + nr ; io(19) = io(18) + nr ; io(20) = io(19) + nr + 1_I4P
   endsubroutine cell_pieces_layout

   subroutine cell_pieces_3d(fill, nf, tri, nrm, nt, lst, x0, x1, y0, y1, z0, z1, cs, alpha, faces, npc, pce, pcq, pct, &
                             pcen, nm2, nm3, e2, e3, qd, want_iv, secm, nrw, rw, niw, iw, npieces, nivc, nqc, nwc, qwn1, &
                             p0, iv0, q0, w0, npo, nivo, nqo, nwo, nqcd, pvol, pcen_, pwa, pwc, pwn, pfc, pmom, pwmom, pfmom, &
                             ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
   !< Kernel (issue #100): the fluid pieces of the cut cell [x0, x1] x [y0, y1] x [z0, z1] with candidates lst and
   !< piercings (pce, pcq, pct, pcen) from `cell_piercings`, on its workspace (rw, iw; `cell_pieces_sizes`). Count pass
   !< (fill false): npieces (-1 if the boundary states do not close or a multi-piece split does not close, see
   !< `split_closes`), the cell's fluid edge intervals for the links (nivc, several pieces), its rule points nqc and wall
   !< rule points nwc, and the wall cones of a single piece qwn1. Fill pass: the same, and the outputs written at the
   !< cell's offsets (pieces p0 + 1.., intervals iv0 + 1.., rule points q0 + 1.., wall rule points w0 + 1..): with
   !< several pieces their measures, moments, intervals and rules; with one (its measures are the cell totals, its faces
   !< the plane sections) its volume and wall moments, its wall cones and wall rule.
   !$acc routine seq
   logical,      intent(in)    :: fill                   !< Fill pass.
   integer(I4P), intent(in)    :: nf                     !< Facets.
   real(R8P),    intent(in)    :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P), intent(in)    :: nt                     !< Candidates of the cell.
   integer(I4P), intent(in)    :: lst(nt)                !< Their indexes.
   real(R8P),    intent(in)    :: x0, x1, y0, y1, z0, z1 !< Cell.
   integer(I4P), intent(in)    :: cs(8)                  !< Corner states (0 solid, 1 fluid), exact.
   real(R8P),    intent(in)    :: alpha, faces(6)        !< Cell fraction, face fractions (-x, +x, -y, +y, -z, +z).
   integer(I4P), intent(in)    :: npc                    !< Piercings.
   integer(I4P), intent(in)    :: pce(npc), pcq(npc), pcen(npc) !< Their edges, candidates, entering flags.
   real(R8P),    intent(in)    :: pct(npc)               !< Their positions.
   integer(I4P), intent(in)    :: nm2, nm3               !< Moments per face, per piece (0: none).
   integer(I4P), intent(in)    :: e2(2, nm2), e3(3, nm3) !< Their exponents.
   logical,      intent(in)    :: qd, want_iv            !< Quadrature, intervals.
   real(R8P),    intent(in)    :: secm(nm2, 6)           !< Section moments of the six faces (a single piece).
   integer(I4P), intent(in)    :: nrw, niw               !< Workspace sizes.
   real(R8P),    intent(inout) :: rw(nrw)                !< Real workspace.
   integer(I4P), intent(inout) :: iw(niw)                !< Integer workspace.
   integer(I4P), intent(out)   :: npieces, nivc, nqc, nwc, qwn1 !< Pieces, intervals, rule points, wall rule points, wall cones.
   integer(I4P), intent(in)    :: p0, iv0, q0, w0        !< Pieces, intervals, rule points, wall rule points before the cell's.
   integer(I4P), intent(in)    :: npo, nivo, nqo, nwo, nqcd !< Sizes of the outputs.
   real(R8P),    intent(inout) :: pvol(npo), pcen_(3, npo), pwa(npo), pwc(3, npo), pwn(3, npo), pfc(6, npo) !< Measures.
   real(R8P),    intent(inout) :: pmom(nm3, npo), pwmom(nm3, npo), pfmom(nm2, 6, npo) !< Moments.
   integer(I4P), intent(inout) :: ivp(nivo), ive(nivo)   !< Intervals: piece, cell edge.
   real(R8P),    intent(inout) :: ivl(nivo), ivh(nivo)   !< Intervals: range.
   real(R8P),    intent(inout) :: qpt(3, nqo), qwt(nqo), wpt(3, nwo), wwt(nwo) !< Rules.
   integer(I4P), intent(inout) :: qcnt(nqcd), wcnt(nqcd) !< Rule points per piece.
   integer(I4P)                :: ro(14), io(20), na, ni, mm, ns, nr !< Layout.

   call cell_pieces_layout(nt, npc, nm2, nm3, qd, ro, io)
   na = nt + 20_I4P + npc ; ni = npc + 12_I4P ; mm = 16_I4P * nt
   ns = 0_I4P ; nr = 0_I4P
   if (qd) then
      ns = 48_I4P * nt + 2_I4P * ni ; nr = na
   endif
   call cell_pieces_core(fill, nf, tri, nrm, nt, lst, x0, x1, y0, y1, z0, z1, cs, alpha, faces, npc, pce, pcq, pct, pcen, &
                         nm2, nm3, e2, e3, qd, want_iv, secm, na, ni, mm, ns, nr, &
                         rw(ro(1)), rw(ro(2)), rw(ro(3)), rw(ro(4)), rw(ro(5)), rw(ro(6)), rw(ro(7)), rw(ro(8)), rw(ro(9)), &
                         rw(ro(10)), rw(ro(11)), rw(ro(12)), rw(ro(13)), &
                         iw(io(1)), iw(io(2)), iw(io(3)), iw(io(4)), iw(io(5)), iw(io(6)), iw(io(7)), iw(io(8)), iw(io(9)), &
                         iw(io(10)), iw(io(11)), iw(io(12)), iw(io(13)), iw(io(14)), iw(io(15)), iw(io(16)), iw(io(17)), &
                         iw(io(18)), iw(io(19)), &
                         npieces, nivc, nqc, nwc, qwn1, p0, iv0, q0, w0, npo, nivo, nqo, nwo, nqcd, pvol, pcen_, pwa, pwc, &
                         pwn, pfc, pmom, pwmom, pfmom, ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
   endsubroutine cell_pieces_3d

   subroutine cell_pieces_core(fill, nf, tri, nrm, nt, lst, x0, x1, y0, y1, z0, z1, cs, alpha, faces, npc, pce, pcq, pct, &
                               pcen, nm2, nm3, e2, e3, qd, want_iv, secm, na, ni, mm, ns, nr, &
                               pv, pt, ivl_, ivh_, vol, keyw, wt, fa, pm, wvm, wwm, fmm, fsg, &
                               pfl, pn, pe, pp, pen, parent, ive_, iva, touch, tgt, pid, iwk, fsf, fsr, qst, wst, qwn, fq, fw, &
                               npieces, nivc, nqc, nwc, qwn1, p0, iv0, q0, w0, npo, nivo, nqo, nwo, nqcd, pvol, pcen_, pwa, &
                               pwc, pwn, pfc, pmom, pwmom, pfmom, ivp, ive, ivl, ivh, qpt, qwt, wpt, wwt, qcnt, wcnt)
   !< `cell_pieces_3d` on the workspace parts (see `cell_pieces_layout`). Fluid pieces of the cell cut by the candidates
   !< lst (see the module notes): atoms are the clipped polygons 1..nt, the corners nt+1..nt+8 and the fluid edge
   !< intervals after them, joined by a union-find whose roots are the least atoms of their sets.
   !<
   !< The corner states cs and the crossings of the 12 cell edges are exact under the global perturbation (issue #70):
   !< a crossing belongs to the edge when it lies after the perturbed low corner and before the perturbed high one
   !< (`cell_piercings`), so the corner states and the crossings of every edge agree. Before #70 the piercings were the
   !< clipped polygon vertices on two planes, and a mesh edge or vertex on a cell edge (lattice bodies) gave one
   !< piercing per incident facet: the states did not alternate.
   !$acc routine seq
   logical,      intent(in)    :: fill                   !< Fill pass.
   integer(I4P), intent(in)    :: nf                     !< Facets.
   real(R8P),    intent(in)    :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P), intent(in)    :: nt                     !< Candidates.
   integer(I4P), intent(in)    :: lst(nt)                !< Their indexes.
   real(R8P),    intent(in)    :: x0, x1, y0, y1, z0, z1 !< Cell.
   integer(I4P), intent(in)    :: cs(8)                  !< Corner states.
   real(R8P),    intent(in)    :: alpha, faces(6)        !< Cell and face fractions.
   integer(I4P), intent(in)    :: npc                    !< Piercings.
   integer(I4P), intent(in)    :: pce(npc), pcq(npc), pcen(npc) !< Edges, candidates, entering.
   real(R8P),    intent(in)    :: pct(npc)               !< Positions.
   integer(I4P), intent(in)    :: nm2, nm3               !< Moments per face, per piece.
   integer(I4P), intent(in)    :: e2(2, nm2), e3(3, nm3) !< Exponents.
   logical,      intent(in)    :: qd, want_iv            !< Quadrature, intervals.
   real(R8P),    intent(in)    :: secm(nm2, 6)           !< Section moments of the faces.
   integer(I4P), intent(in)    :: na, ni, mm, ns, nr     !< Atoms, intervals, polygon vertices, face segments, rule pieces.
   real(R8P),    intent(inout) :: pv(3, 16, nt)          !< Polygons.
   real(R8P),    intent(inout) :: pt(npc)                !< Piercing positions.
   real(R8P),    intent(inout) :: ivl_(ni), ivh_(ni)     !< Fluid intervals: range.
   real(R8P),    intent(inout) :: vol(na)                !< Solid volume of closed wall sets.
   real(R8P),    intent(inout) :: keyw(4 * mm)           !< Sort keys, face segment end points.
   real(R8P),    intent(inout) :: wt(11, na), fa(6, na)  !< Wall terms, face fluid areas (per root).
   real(R8P),    intent(inout) :: pm(17, na)             !< Measures (per piece).
   real(R8P),    intent(inout) :: wvm(nm3, na), wwm(nm3, na), fmm(nm2, 6, na) !< Moments: wall divergence, wall, faces.
   real(R8P),    intent(inout) :: fsg(6, ns)             !< Quadrature: oriented face segments.
   integer(I4P), intent(inout) :: pfl(16, nt), pn(nt)    !< Vertex flags, vertex counts.
   integer(I4P), intent(inout) :: pe(npc), pp(npc), pen(npc) !< Piercings: edge, polygon (0: none in the cell), entering.
   integer(I4P), intent(inout) :: parent(na)             !< Union-find.
   integer(I4P), intent(inout) :: ive_(ni), iva(ni)      !< Fluid intervals: edge, atom.
   integer(I4P), intent(inout) :: touch(na), tgt(na), pid(na) !< Touch flags, inclusion targets, piece ids.
   integer(I4P), intent(inout) :: iwk(3 * mm)            !< Sort and face segment work.
   integer(I4P), intent(inout) :: fsf(ns), fsr(ns)       !< Segment face and set.
   integer(I4P), intent(inout) :: qst(nr + 1), wst(nr + 1) !< Rule starts (per piece).
   integer(I4P), intent(inout) :: qwn(nr), fq(nr), fw(nr) !< Wall cones per piece, fill pointers.
   integer(I4P), intent(out)   :: npieces, nivc, nqc, nwc, qwn1 !< Counts.
   integer(I4P), intent(in)    :: p0, iv0, q0, w0, npo, nivo, nqo, nwo, nqcd !< Offsets, output sizes.
   real(R8P),    intent(inout) :: pvol(npo), pcen_(3, npo), pwa(npo), pwc(3, npo), pwn(3, npo), pfc(6, npo) !< Measures.
   real(R8P),    intent(inout) :: pmom(nm3, npo), pwmom(nm3, npo), pfmom(nm2, 6, npo) !< Moments.
   integer(I4P), intent(inout) :: ivp(nivo), ive(nivo)   !< Intervals: piece, edge.
   real(R8P),    intent(inout) :: ivl(nivo), ivh(nivo)   !< Intervals: range.
   real(R8P),    intent(inout) :: qpt(3, nqo), qwt(nqo), wpt(3, nwo), wwt(nwo) !< Rules.
   integer(I4P), intent(inout) :: qcnt(nqcd), wcnt(nqcd) !< Rule points per piece.
   integer(I4P), parameter     :: NCORN = 8_I4P          !< Corner atoms follow the polygons.
   real(R8P)                   :: c0(3), c1(3), d0(3), e1(3), e2v(3), cr(3), h3(3) !< Cell corners, work.
   integer(I4P)                :: cstate(8)              !< Corner states: -1 unknown, 0 solid, 1 fluid.
   integer(I4P)                :: q, v, w, e, a, natom, niv, nc, npv, cpp, nn, nfs !< Work.
   real(R8P)                   :: tt                     !< Crossing position.
   logical                     :: ok, any_fl, measures   !< Consistency, work, measures needed.

   c0(1) = x0 ; c0(2) = y0 ; c0(3) = z0 ; c1(1) = x1 ; c1(2) = y1 ; c1(3) = z1
   nivc = 0_I4P ; nqc = 0_I4P ; nwc = 0_I4P ; qwn1 = 0_I4P
   do q = 1, nt
      pv(:, 1:3, q) = tri(:, :, lst(q)) ; pfl(1:3, q) = 0_I4P ; npv = 3_I4P
      do a = 1, 3
         call clip_flags(pv(:, :, q), pfl(:, q), npv, a, c0(a), .true., 2_I4P * (a - 1_I4P))
         call clip_flags(pv(:, :, q), pfl(:, q), npv, a, c1(a), .false., 2_I4P * (a - 1_I4P) + 1_I4P)
      enddo
      pn(q) = npv ; if (npv < 3_I4P) pn(q) = 0_I4P
      ! on-plane flags from the coordinates (issue #70): a vertex made by one clip can lie exactly on a plane
      ! clipped before or after it (lattice bodies), and the traces and edge piercings need that flag too
      do v = 1, pn(q)
         pfl(v, q) = 0_I4P
         do a = 1, 3
            if (pv(a, v, q) == c0(a)) pfl(v, q) = ibset(pfl(v, q), 2_I4P * (a - 1_I4P))
            if (pv(a, v, q) == c1(a)) pfl(v, q) = ibset(pfl(v, q), 2_I4P * (a - 1_I4P) + 1_I4P)
         enddo
      enddo
   enddo
   ! piercings (issue #70): the exact crossings of `cell_piercings`, at the clipped polygon's vertex on the edge if it
   ! has one (same rounding), else at the kernel's position clamped to the edge
   do w = 1, npc
      e = pce(w) ; q = pcq(w) ; a = (e - 1_I4P) / 4_I4P + 1_I4P
      tt = min(max(pct(w), c0(a)), c1(a))
      do v = 1, pn(q)
         if (edge_of(pfl(v, q)) == e) then
            tt = pv(a, v, q) ; exit
         endif
      enddo
      pe(w) = e ; pt(w) = tt ; pen(w) = pcen(w) ; pp(w) = merge(q, 0_I4P, pn(q) > 0_I4P)
   enddo
   call sort_pierce(npc, pe, pt, pen, pp)
   ! corner states from the pierced edges, then propagation along the others
   ok = .true. ; cstate = cs ; natom = nt + NCORN ; niv = 0_I4P
   do e = 1, 12 ! check: the crossings alternate from the low corner's state to the high one's
      call walk_edge(e, .false., c0, c1, nt, npc, pe, pt, pen, pp, cstate, ok, na, parent, natom, ni, niv, ive_, iva, &
                     ivl_, ivh_)
      if (.not. ok) exit
   enddo
   do a = 1, 8
      if (cstate(a) == -1_I4P) ok = .false.
   enddo
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   ! atoms: polygons 1..nt, corners nt+1..nt+8, intervals after
   do a = 1, na
      parent(a) = a
   enddo
   do e = 1, 12
      call walk_edge(e, .true., c0, c1, nt, npc, pe, pt, pen, pp, cstate, ok, na, parent, natom, ni, niv, ive_, iva, &
                     ivl_, ivh_)
      if (.not. ok) exit
   enddo
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   if (nt > 0_I4P) call join_shared_vertices(nt, pv, pn, mm, keyw, iwk(1), iwk(mm + 1), na, parent)
   do a = 1, merge(3_I4P, 0_I4P, nt > 0_I4P) ! (no polygon: no trace)
      call faces_islands(a, 0_I4P, nf, nrm, nt, lst, pv, pfl, pn, mm, keyw(1), keyw(2 * mm + 1), iwk(1), iwk(mm + 1), &
                         iwk(2 * mm + 1), na, parent, ni, niv, ive_, iva, ivl_, ivh_, ok)
      if (.not. ok) exit
      call faces_islands(a, 1_I4P, nf, nrm, nt, lst, pv, pfl, pn, mm, keyw(1), keyw(2 * mm + 1), iwk(1), iwk(mm + 1), &
                         iwk(2 * mm + 1), na, parent, ni, niv, ive_, iva, ivl_, ivh_, ok)
      if (.not. ok) exit
   enddo
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   ! count: sets touching the boundary, plus closed wall sets that bound a void
   do a = 1, natom
      touch(a) = 0_I4P ; vol(a) = 0._R8P
   enddo
   do q = 1, nt
      if (pn(q) == 0_I4P) cycle
      w = uf_find(na, parent, q)
      any_fl = .false.
      do v = 1, pn(q)
         if (pfl(v, q) /= 0_I4P) any_fl = .true.
      enddo
      if (any_fl) touch(w) = 1_I4P
      do v = 2, pn(q) - 1
         d0 = pv(:, 1, q) - c0 ; e1 = pv(:, v, q) - pv(:, 1, q) ; e2v = pv(:, v + 1, q) - pv(:, 1, q)
         call cross3_into(e1, e2v, cr)
         vol(w) = vol(w) + dot_product(d0, cr) / 6._R8P
      enddo
   enddo
   do a = 1, NCORN
      if (cstate(a) == 1_I4P) touch(uf_find(na, parent, nt + a)) = 1_I4P
   enddo
   do a = 1, niv
      touch(uf_find(na, parent, iva(a))) = 1_I4P
   enddo
   nc = 0_I4P
   do a = 1, natom
      if (uf_find(na, parent, a) /= a) cycle
      if (touch(a) == 1_I4P) then
         nc = nc + 1_I4P
      elseif (a <= nt) then
         if (pn(a) > 0_I4P .and. vol(a) < 0._R8P) nc = nc + 1_I4P
      endif
   enddo
   npieces = nc
   ! inclusions, piece ids (always: a missing target is no count); the measures where the pass needs them
   call piece_ids(nf, tri, nrm, nt, lst, c1, pv, pfl, pn, na, natom, parent, touch, vol, ni, niv, ive_, iva, ivl_, ivh_, &
                  tgt, pid, nn)
   if (nn < 0_I4P) then
      npieces = -1_I4P ; return
   endif
   measures = fill .or. npieces > 1_I4P .or. qd
   if (.not. measures) return
   call piece_measures(fill, nf, tri, nrm, nt, lst, c0, c1, nm2, nm3, e2, e3, qd, pv, pfl, pn, na, natom, parent, tgt, &
                       pid, ni, niv, ive_, iva, ivl_, ivh_, wt, fa, pm, wvm, wwm, fmm, ns, nfs, fsg, fsf, fsr, nn, qst, wst, &
                       qwn, fw, w0, nwo, wpt, wwt)
   if (npieces > 1_I4P) then ! (issue #70) trust the split only if it closes
      h3 = c1 - c0
      if (.not. split_closes(npieces, pm, h3, alpha, faces)) then
         npieces = -1_I4P ; return
      endif
   endif
   ! counts
   if (want_iv .and. npieces > 1_I4P) then
      do a = 1, niv
         if (pid(uf_find(na, parent, iva(a))) >= 1_I4P) nivc = nivc + 1_I4P
      enddo
   endif
   if (qd) then
      nqc = qst(nn + 1) - 1_I4P ; nwc = wst(nn + 1) - 1_I4P
      if (nn >= 1_I4P) qwn1 = qwn(1)
   endif
   if (.not. fill) return
   ! fill: several pieces, their measures, moments and intervals; one, its volume and wall moments; the rules
   if (npieces > 1_I4P) then
      do a = 1, nn
         pvol(p0 + a) = pm(1, a)
         pcen_(:, p0 + a) = pm(2:4, a) ; pwa(p0 + a) = pm(5, a) ; pwc(:, p0 + a) = pm(6:8, a)
         pwn(:, p0 + a) = pm(9:11, a) ; pfc(:, p0 + a) = pm(12:17, a)
      enddo
      if (want_iv) then
         w = iv0
         do a = 1, niv
            cpp = pid(uf_find(na, parent, iva(a)))
            if (cpp < 1_I4P) cycle
            w = w + 1_I4P
            ivp(w) = p0 + cpp ; ive(w) = ive_(a) ; ivl(w) = ivl_(a) ; ivh(w) = ivh_(a)
         enddo
      endif
   endif
   if (nm3 > 0_I4P) then
      do q = 1, natom
         if (pid(q) == 0_I4P) cycle
         a = pid(q)
         pwmom(:, p0 + a) = wwm(:, q)
         if (npieces > 1_I4P) then
            pfmom(:, :, p0 + a) = fmm(:, :, q)
            call volume_moments_k(nm3, nm2, wvm(:, q), fmm(:, 1, q), fmm(:, 2, q), x1 - x0, e2, e3, pmom(:, p0 + a))
         else
            call volume_moments_k(nm3, nm2, wvm(:, q), secm(:, 1), secm(:, 2), x1 - x0, e2, e3, pmom(:, p0 + a))
         endif
      enddo
   endif
   if (qd) call piece_rules(npieces > 1_I4P, nt, c0, c1, pv, pn, na, parent, tgt, pid, nn, ns, nfs, fsg, fsf, fsr, qst, &
                            wst, fq, p0, q0, nqo, qpt, qwt, nqcd, qcnt, wcnt)
   endsubroutine cell_pieces_core

   pure integer(I4P) function uf_find(na, parent, i) result(r)
   !< Root of atom i (the least atom of its set; no path compression: the sets are small).
   !$acc routine seq
   integer(I4P), intent(in) :: na         !< Atoms.
   integer(I4P), intent(in) :: parent(na) !< Union-find.
   integer(I4P), intent(in) :: i          !< Atom.

   r = i
   do while (parent(r) /= r)
      r = parent(r)
   enddo
   endfunction uf_find

   pure subroutine uf_unite(na, parent, i, j)
   !< Join the sets of i and j: the larger root points to the smaller.
   !$acc routine seq
   integer(I4P), intent(in)    :: na         !< Atoms.
   integer(I4P), intent(inout) :: parent(na) !< Union-find.
   integer(I4P), intent(in)    :: i, j       !< Atoms.
   integer(I4P)                :: ri, rj     !< Roots.

   ri = uf_find(na, parent, i) ; rj = uf_find(na, parent, j)
   if (ri /= rj) parent(max(ri, rj)) = min(ri, rj)
   endsubroutine uf_unite

   pure subroutine walk_edge(ed, build, c0, c1, nt, npc, pe, pt, pen, pp, cstate, ok, na, parent, natom, ni, niv, ive, iva, &
                             ivl, ivh)
   !< Walk the piercings of edge ed from its low corner's exact state. Without `build`: check that they alternate and
   !< end at the high corner's state, reordering crossings at one rounded position (a mesh edge or vertex on the cell
   !< edge: their exact order is the alternation). With `build`: the fluid interval atoms and their joins.
   !$acc routine seq
   integer(I4P), intent(in)    :: ed                   !< Edge.
   logical,      intent(in)    :: build                !< Build the intervals.
   real(R8P),    intent(in)    :: c0(3), c1(3)         !< Cell corners.
   integer(I4P), intent(in)    :: nt, npc              !< Polygons, piercings.
   integer(I4P), intent(inout) :: pe(npc), pen(npc), pp(npc) !< Piercings: edge, entering, polygon.
   real(R8P),    intent(inout) :: pt(npc)              !< Piercings: position.
   integer(I4P), intent(inout) :: cstate(8)            !< Corner states.
   logical,      intent(inout) :: ok                   !< Consistency.
   integer(I4P), intent(in)    :: na                   !< Atoms (bound).
   integer(I4P), intent(inout) :: parent(na)           !< Union-find.
   integer(I4P), intent(inout) :: natom                !< Atoms so far.
   integer(I4P), intent(in)    :: ni                   !< Intervals (bound).
   integer(I4P), intent(inout) :: niv                  !< Intervals so far.
   integer(I4P), intent(inout) :: ive(ni), iva(ni)     !< Intervals: edge, atom.
   real(R8P),    intent(inout) :: ivl(ni), ivh(ni)     !< Intervals: range.
   integer(I4P)                :: st, ii, i1, i2, prev, ax, jt, kt, cl, ch, ti !< Work.
   real(R8P)                   :: tprev, lo, hi, tr   !< Previous position, edge range, swap.

   ax = (ed - 1_I4P) / 4_I4P + 1_I4P
   lo = c0(ax) ; hi = c1(ax)
   cl = corner_of(ed, 0_I4P) ; ch = corner_of(ed, 1_I4P)
   i1 = 0_I4P ; i2 = -1_I4P
   do ii = 1, npc
      if (pe(ii) /= ed) cycle
      if (i1 == 0_I4P) i1 = ii
      i2 = ii
   enddo
   if (i1 == 0_I4P) then ! no piercing
      if (build .and. cstate(cl) == 1_I4P) then
         natom = natom + 1_I4P ; niv = niv + 1_I4P
         ive(niv) = ed ; ivl(niv) = lo ; ivh(niv) = hi ; iva(niv) = natom
         call uf_unite(na, parent, natom, nt + cl) ; call uf_unite(na, parent, natom, nt + ch)
      endif
      return
   endif
   st = cstate(cl) ! exact (issue #70); crossings at one rounded position are reordered below
   prev = nt + cl ; tprev = lo
   do ii = i1, i2
      if ((pen(ii) == 1_I4P) .neqv. (st == 1_I4P)) then ! alternation broken: repair a tie, else not closed
         jt = 0_I4P
         do kt = ii + 1_I4P, i2
            if (abs(pt(kt) - pt(ii)) > 64._R8P * spacing(max(abs(lo), abs(hi)))) exit
            if ((pen(kt) == 1_I4P) .eqv. (st == 1_I4P)) then
               jt = kt ; exit
            endif
         enddo
         if (jt == 0_I4P) then
            ok = .false. ; return
         endif
         ti = pe(ii) ; pe(ii) = pe(jt) ; pe(jt) = ti
         tr = pt(ii) ; pt(ii) = pt(jt) ; pt(jt) = tr
         ti = pen(ii) ; pen(ii) = pen(jt) ; pen(jt) = ti
         ti = pp(ii) ; pp(ii) = pp(jt) ; pp(jt) = ti
      endif
      if (build .and. st == 1_I4P) then
         natom = natom + 1_I4P ; niv = niv + 1_I4P
         ive(niv) = ed ; ivl(niv) = tprev ; ivh(niv) = pt(ii) ; iva(niv) = natom
         if (prev > 0_I4P) call uf_unite(na, parent, natom, prev) ! 0: a crossing without a polygon in the cell
         if (pp(ii) > 0_I4P) call uf_unite(na, parent, natom, pp(ii))
      endif
      st = 1_I4P - st ; prev = pp(ii) ; tprev = pt(ii)
   enddo
   if (.not. build) then ! set (or check) the high corner's state
      if (cstate(ch) == -1_I4P) then
         cstate(ch) = st
      elseif (cstate(ch) /= st) then
         ok = .false.
      endif
   endif
   if (build .and. st == 1_I4P) then
      natom = natom + 1_I4P ; niv = niv + 1_I4P
      ive(niv) = ed ; ivl(niv) = tprev ; ivh(niv) = hi ; iva(niv) = natom
      if (prev > 0_I4P) call uf_unite(na, parent, natom, prev)
      call uf_unite(na, parent, natom, nt + ch)
   endif
   endsubroutine walk_edge

   pure subroutine join_shared_vertices(nt, pv, pn, mm, key, own, idx, na, parent)
   !< Join the polygons that share a vertex (bit-identical points).
   !$acc routine seq
   integer(I4P), intent(in)    :: nt            !< Polygons.
   real(R8P),    intent(in)    :: pv(3, 16, nt) !< Polygons.
   integer(I4P), intent(in)    :: pn(nt)        !< Vertex counts.
   integer(I4P), intent(in)    :: mm            !< Vertices (bound).
   real(R8P),    intent(inout) :: key(3, mm)    !< Vertices.
   integer(I4P), intent(inout) :: own(mm), idx(mm) !< Owner polygon, sorted order.
   integer(I4P), intent(in)    :: na            !< Atoms.
   integer(I4P), intent(inout) :: parent(na)    !< Union-find.
   integer(I4P)                :: m, ii, qq, vv !< Work.

   m = 0_I4P
   do qq = 1, nt
      do vv = 1, pn(qq)
         m = m + 1_I4P ; key(:, m) = pv(:, vv, qq) ; own(m) = qq
      enddo
   enddo
   call sort_keys(3_I4P, m, key, idx)
   do ii = 2, m
      if (key(1, idx(ii)) == key(1, idx(ii - 1)) .and. key(2, idx(ii)) == key(2, idx(ii - 1)) .and. &
          key(3, idx(ii)) == key(3, idx(ii - 1))) call uf_unite(na, parent, own(idx(ii)), own(idx(ii - 1)))
   enddo
   endsubroutine join_shared_vertices

   pure subroutine faces_islands(ax, sd, nf, nrm, nt, lst, pv, pfl, pn, mm, su, sv, sp, sc, stouch, na, parent, ni, niv, &
                                 ive, iva, ivl, ivh, ok)
   !< Join every solid island traced on the face (axis ax, side sd) to the fluid around it: a loop of wall polygon edges
   !< on the face that touches no cell edge, the fluid outside it, joins the nearest trace, or the fluid interval of the
   !< face edge, a ray in +u from its rightmost point reaches.
   !$acc routine seq
   integer(I4P), intent(in)    :: ax, sd           !< Face.
   integer(I4P), intent(in)    :: nf               !< Facets.
   real(R8P),    intent(in)    :: nrm(3, nf)       !< Normals.
   integer(I4P), intent(in)    :: nt               !< Polygons.
   integer(I4P), intent(in)    :: lst(nt)          !< Their facets.
   real(R8P),    intent(in)    :: pv(3, 16, nt)    !< Polygons.
   integer(I4P), intent(in)    :: pfl(16, nt), pn(nt) !< Vertex flags, counts.
   integer(I4P), intent(in)    :: mm               !< Segments (bound).
   real(R8P),    intent(inout) :: su(2, mm), sv(2, mm) !< Segment end points in face coordinates.
   integer(I4P), intent(inout) :: sp(mm), sc(mm), stouch(mm) !< Segment polygon, component, touch flag.
   integer(I4P), intent(in)    :: na               !< Atoms.
   integer(I4P), intent(inout) :: parent(na)       !< Union-find.
   integer(I4P), intent(in)    :: ni, niv          !< Intervals: bound, number.
   integer(I4P), intent(in)    :: ive(ni), iva(ni) !< Intervals: edge, atom.
   real(R8P),    intent(in)    :: ivl(ni), ivh(ni) !< Intervals: range.
   logical,      intent(inout) :: ok               !< Consistency.
   integer(I4P)                :: bt, a1, a2, ns, ii, jj, qq, vv, best, comp, sides(3), ed, iv, l !< Work.
   real(R8P)                   :: uv(2), vv0, uc, ubest, nbest !< Work.
   logical                     :: found            !< Work.

   bt = 2_I4P * (ax - 1_I4P) + sd
   a1 = mod(ax, 3_I4P) + 1_I4P ; a2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P
   ns = 0_I4P
   do qq = 1, nt
      do vv = 1, pn(qq)
         jj = mod(vv, pn(qq)) + 1
         if (.not. (btest(pfl(vv, qq), bt) .and. btest(pfl(jj, qq), bt))) cycle
         ns = ns + 1_I4P
         su(1, ns) = pv(a1, vv, qq) ; su(2, ns) = pv(a1, jj, qq) ; sv(1, ns) = pv(a2, vv, qq) ; sv(2, ns) = pv(a2, jj, qq)
         sp(ns) = qq
         stouch(ns) = 0_I4P
         if (edge_of(pfl(vv, qq)) /= 0_I4P .or. edge_of(pfl(jj, qq)) /= 0_I4P .or. &
             popcnt(pfl(vv, qq)) > 2 .or. popcnt(pfl(jj, qq)) > 2) stouch(ns) = 1_I4P
      enddo
   enddo
   if (ns == 0_I4P) return
   ! segment components by shared end points (small sets: quadratic)
   do ii = 1, ns
      sc(ii) = ii
   enddo
   do ii = 1, ns
      do jj = ii + 1, ns
         if (sc(ii) == sc(jj)) cycle
         if ((su(1, ii) == su(1, jj) .and. sv(1, ii) == sv(1, jj)) .or. (su(1, ii) == su(2, jj) .and. sv(1, ii) == sv(2, jj)) &
             .or. (su(2, ii) == su(1, jj) .and. sv(2, ii) == sv(1, jj)) .or. &
             (su(2, ii) == su(2, jj) .and. sv(2, ii) == sv(2, jj))) then
            iv = max(sc(ii), sc(jj)) ; best = min(sc(ii), sc(jj))
            do l = 1, ns
               if (sc(l) == iv) sc(l) = best
            enddo
         endif
      enddo
   enddo
   do ii = 1, ns
      found = .false.
      do l = 1, ns
         if (stouch(l) == 1_I4P .and. sc(l) == sc(ii)) found = .true.
      enddo
      if (found) stouch(ii) = 1_I4P
   enddo
   do comp = 1, ns
      ii = 0_I4P
      do l = 1, ns
         if (sc(l) == comp) then
            ii = l ; exit
         endif
      enddo
      if (ii == 0_I4P) cycle
      if (stouch(ii) == 1_I4P) cycle
      ! rightmost point of the loop and the incident segment with the largest |n_u|
      uv(1) = -huge(1._R8P) ; uv(2) = -huge(1._R8P)
      do jj = 1, ns
         if (sc(jj) /= comp) cycle
         do iv = 1, 2
            if (su(iv, jj) > uv(1) .or. (su(iv, jj) == uv(1) .and. sv(iv, jj) > uv(2))) then
               uv(1) = su(iv, jj) ; uv(2) = sv(iv, jj)
            endif
         enddo
      enddo
      best = 0_I4P ; nbest = 0._R8P
      do jj = 1, ns
         if (sc(jj) /= comp) cycle
         if (.not. ((su(1, jj) == uv(1) .and. sv(1, jj) == uv(2)) .or. (su(2, jj) == uv(1) .and. sv(2, jj) == uv(2)))) cycle
         if (abs(nrm(a1, lst(sp(jj)))) > abs(nbest)) then
            best = jj ; nbest = nrm(a1, lst(sp(jj)))
         endif
      enddo
      if (best == 0_I4P .or. nbest <= 0._R8P) cycle ! a hole (fluid inside the loop): nothing to join
      ! ray in +u from the rightmost point: the nearest trace of another set, else the face edge
      vv0 = uv(2) ; ubest = huge(1._R8P) ; qq = 0_I4P
      do jj = 1, ns
         if (sc(jj) == comp) cycle
         if ((sv(1, jj) > vv0) .eqv. (sv(2, jj) > vv0)) cycle
         uc = su(1, jj) + (vv0 - sv(1, jj)) / (sv(2, jj) - sv(1, jj)) * (su(2, jj) - su(1, jj))
         if (uc > uv(1) .and. uc < ubest) then
            ubest = uc ; qq = sp(jj)
         endif
      enddo
      if (qq /= 0_I4P) then
         call uf_unite(na, parent, sp(best), qq)
         cycle
      endif
      sides = 0_I4P ; sides(ax) = sd ; sides(a1) = 1_I4P
      ed = edge_index(a2, sides)
      qq = 0_I4P
      do jj = 1, niv
         if (ive(jj) == ed .and. ivl(jj) <= vv0 .and. vv0 <= ivh(jj)) qq = iva(jj)
      enddo
      if (qq == 0_I4P) then
         ok = .false. ; return
      endif
      call uf_unite(na, parent, sp(best), qq)
   enddo
   endsubroutine faces_islands

   pure subroutine piece_ids(nf, tri, nrm, nt, lst, c1, pv, pfl, pn, na, natom, parent, touch, vol, ni, niv, ive, iva, &
                             ivl, ivh, tgt, pid, nn)
   !< Inclusion targets and piece ids. An inclusion (a closed set of positive solid volume) goes to the set around it
   !< (`inclusion_target`; nn = -1 if one is not found); the counted sets (touching the boundary, or closed voids) get
   !< ids 1..nn in the order of their roots.
   !$acc routine seq
   integer(I4P), intent(in)    :: nf               !< Facets.
   real(R8P),    intent(in)    :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P), intent(in)    :: nt               !< Polygons.
   integer(I4P), intent(in)    :: lst(nt)          !< Their facets.
   real(R8P),    intent(in)    :: c1(3)            !< Cell high corner.
   real(R8P),    intent(in)    :: pv(3, 16, nt)    !< Polygons.
   integer(I4P), intent(in)    :: pfl(16, nt), pn(nt) !< Vertex flags, counts.
   integer(I4P), intent(in)    :: na, natom        !< Atoms: bound, number.
   integer(I4P), intent(in)    :: parent(na)       !< Union-find.
   integer(I4P), intent(in)    :: touch(na)        !< Touch flags.
   real(R8P),    intent(in)    :: vol(na)          !< Solid volumes.
   integer(I4P), intent(in)    :: ni, niv          !< Intervals: bound, number.
   integer(I4P), intent(in)    :: ive(ni), iva(ni) !< Intervals: edge, atom.
   real(R8P),    intent(in)    :: ivl(ni), ivh(ni) !< Intervals: range.
   integer(I4P), intent(inout) :: tgt(na), pid(na) !< Inclusion targets, piece ids.
   integer(I4P), intent(out)   :: nn               !< Pieces (-1: an inclusion without a target).
   integer(I4P)                :: r, kk, guard     !< Work.

   do r = 1, natom
      tgt(r) = 0_I4P ; pid(r) = 0_I4P
   enddo
   nn = -1_I4P
   do r = 1, nt
      if (uf_find(na, parent, r) /= r .or. touch(r) == 1_I4P .or. pn(r) == 0_I4P .or. vol(r) < 0._R8P) cycle
      tgt(r) = inclusion_target(r, nf, tri, nrm, nt, lst, c1, pv, pfl, pn, na, parent, ni, niv, ive, iva, ivl, ivh)
      if (tgt(r) == 0_I4P) return
   enddo
   do r = 1, nt
      if (tgt(r) == 0_I4P) cycle
      kk = tgt(r) ; guard = 0_I4P
      do while (tgt(kk) /= 0_I4P .and. guard < natom)
         kk = tgt(kk) ; guard = guard + 1_I4P
      enddo
      tgt(r) = kk
   enddo
   nn = 0_I4P
   do r = 1, natom
      if (uf_find(na, parent, r) /= r) cycle
      if (touch(r) == 1_I4P .or. (r <= nt .and. tgt(r) == 0_I4P .and. vol(r) < 0._R8P)) then
         nn = nn + 1_I4P ; pid(r) = nn
      endif
   enddo
   endsubroutine piece_ids

   pure integer(I4P) function inclusion_target(ri, nf, tri, nrm, nt, lst, c1, pv, pfl, pn, na, parent, ni, niv, ive, iva, &
                                               ivl, ivh) result(tr)
   !< Root of the set around the inclusion ri: the nearest wall a ray in +x from its rightmost vertex hits (exact
   !< inclusive test on the projected triangle), else the fluid region where the ray leaves the cell through the +x face
   !< (the ray of `faces_islands` on that face); 0 if not found.
   !$acc routine seq
   integer(I4P), intent(in) :: ri               !< Root of the inclusion.
   integer(I4P), intent(in) :: nf               !< Facets.
   real(R8P),    intent(in) :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P), intent(in) :: nt               !< Polygons.
   integer(I4P), intent(in) :: lst(nt)          !< Their facets.
   real(R8P),    intent(in) :: c1(3)            !< Cell high corner.
   real(R8P),    intent(in) :: pv(3, 16, nt)    !< Polygons.
   integer(I4P), intent(in) :: pfl(16, nt), pn(nt) !< Vertex flags, counts.
   integer(I4P), intent(in) :: na               !< Atoms.
   integer(I4P), intent(in) :: parent(na)       !< Union-find.
   integer(I4P), intent(in) :: ni, niv          !< Intervals: bound, number.
   integer(I4P), intent(in) :: ive(ni), iva(ni) !< Intervals: edge, atom.
   real(R8P),    intent(in) :: ivl(ni), ivh(ni) !< Intervals: range.
   real(R8P)                :: pr(3), t(3, 3), nn(3), xh, xb, o(3), ub, uc !< Work.
   integer(I4P)             :: q2, v2, j2, qb, ed, sides(3) !< Work.

   pr = -huge(1._R8P)
   do q2 = 1, nt
      if (pn(q2) == 0_I4P) cycle
      if (uf_find(na, parent, q2) /= ri) cycle
      do v2 = 1, pn(q2)
         if (pv(1, v2, q2) > pr(1)) pr = pv(:, v2, q2)
      enddo
   enddo
   xb = huge(1._R8P) ; qb = 0_I4P
   do q2 = 1, nt
      if (pn(q2) == 0_I4P) cycle
      if (uf_find(na, parent, q2) == ri) cycle
      t = tri(:, :, lst(q2)) ; nn = nrm(:, lst(q2))
      if (nn(1) == 0._R8P) cycle
      if (orient2d_xy(t(2, 1), t(3, 1), t(2, 2), t(3, 2), t(2, 3), t(3, 3)) == 0._R8P) cycle
      o(1) = orient2d_xy(t(2, 1), t(3, 1), t(2, 2), t(3, 2), pr(2), pr(3))
      o(2) = orient2d_xy(t(2, 2), t(3, 2), t(2, 3), t(3, 3), pr(2), pr(3))
      o(3) = orient2d_xy(t(2, 3), t(3, 3), t(2, 1), t(3, 1), pr(2), pr(3))
      if (.not. ((o(1) >= 0._R8P .and. o(2) >= 0._R8P .and. o(3) >= 0._R8P) .or. &
                 (o(1) <= 0._R8P .and. o(2) <= 0._R8P .and. o(3) <= 0._R8P))) cycle
      xh = t(1, 1) - (nn(2) * (pr(2) - t(2, 1)) + nn(3) * (pr(3) - t(3, 1))) / nn(1)
      if (xh > pr(1) .and. xh <= c1(1) .and. xh < xb) then
         xb = xh ; qb = q2
      endif
   enddo
   if (qb /= 0_I4P) then
      tr = uf_find(na, parent, qb) ; return
   endif
   ! the +x face, (u, v) = (y, z): the nearest trace in +u, else the fluid interval of the edge u = high
   ub = huge(1._R8P) ; qb = 0_I4P
   do q2 = 1, nt
      do v2 = 1, pn(q2)
         j2 = mod(v2, pn(q2)) + 1_I4P
         if (.not. (btest(pfl(v2, q2), 1) .and. btest(pfl(j2, q2), 1))) cycle
         if ((pv(3, v2, q2) > pr(3)) .eqv. (pv(3, j2, q2) > pr(3))) cycle
         uc = pv(2, v2, q2) + (pr(3) - pv(3, v2, q2)) / (pv(3, j2, q2) - pv(3, v2, q2)) * (pv(2, j2, q2) - pv(2, v2, q2))
         if (uc > pr(2) .and. uc < ub) then
            ub = uc ; qb = q2
         endif
      enddo
   enddo
   if (qb /= 0_I4P) then
      tr = uf_find(na, parent, qb) ; return
   endif
   sides(1) = 1_I4P ; sides(2) = 1_I4P ; sides(3) = 0_I4P ; ed = edge_index(3_I4P, sides)
   tr = 0_I4P
   do j2 = 1, niv
      if (ive(j2) == ed .and. ivl(j2) <= pr(3) .and. pr(3) <= ivh(j2)) tr = uf_find(na, parent, iva(j2))
   enddo
   endfunction inclusion_target

   pure logical function across_crosses(fi, ax, lev, nf, tri, nt, lst) result(yes)
   !< The facet of the cell across the edge of facet fi that lies on the plane x(ax) = lev has a vertex beyond the
   !< plane (x(ax) > lev); .false. if fi has no edge on the plane or no facet of the cell shares it.
   !$acc routine seq
   integer(I4P), intent(in) :: fi            !< Facet.
   integer(I4P), intent(in) :: ax            !< Plane axis.
   real(R8P),    intent(in) :: lev           !< Plane level.
   integer(I4P), intent(in) :: nf            !< Facets.
   real(R8P),    intent(in) :: tri(3, 3, nf) !< Triangles.
   integer(I4P), intent(in) :: nt            !< Candidates.
   integer(I4P), intent(in) :: lst(nt)       !< Their facets.
   integer(I4P)             :: m1, m2, g, mm, nshare, fg !< Work.

   yes = .false.
   do m1 = 1, 3
      m2 = mod(m1, 3_I4P) + 1_I4P
      if (tri(ax, m1, fi) /= lev .or. tri(ax, m2, fi) /= lev) cycle
      do g = 1, nt
         fg = lst(g)
         if (fg == fi) cycle
         nshare = 0_I4P
         do mm = 1, 3
            if ((tri(1, mm, fg) == tri(1, m1, fi) .and. tri(2, mm, fg) == tri(2, m1, fi) .and. &
                 tri(3, mm, fg) == tri(3, m1, fi)) .or. (tri(1, mm, fg) == tri(1, m2, fi) .and. &
                 tri(2, mm, fg) == tri(2, m2, fi) .and. tri(3, mm, fg) == tri(3, m2, fi))) nshare = nshare + 1_I4P
         enddo
         if (nshare == 2_I4P .and. max(tri(ax, 1, fg), tri(ax, 2, fg), tri(ax, 3, fg)) > lev) then
            yes = .true. ; return
         endif
      enddo
   enddo
   endfunction across_crosses

   subroutine piece_measures(fill, nf, tri, nrm, nt, lst, c0, c1, nm2, nm3, e2, e3, qd, pv, pfl, pn, na, natom, parent, &
                             tgt, pid, ni, niv, ive, iva, ivl, ivh, wt, fa, pm, wvm, wwm, fmm, ns, nfs, fsg, fsf, fsr, nn, &
                             qst, wst, qwn, fw, w0, nwo, wpt, wwt)
   !< Measures of the counted sets (ids 1..nn): pm(1) fluid volume, pm(2:4) centroid, pm(5) wetted area, pm(6:8) wetted
   !< centroid, pm(9:11) outward normal integral, pm(12:17) fluid fraction of the faces -x, +x, -y, +y, -z, +z.
   !<
   !< Walls: the divergence-theorem terms of every polygon go to its set (the outward normal of the fluid is minus the
   !< triangle's). Faces: the fluid area of each face by Green's theorem in face coordinates local to c0,
   !< A = oint u dv counter-clockwise seen from +e_a, over the boundary arcs of the fluid region of the face: the fluid
   !< edge intervals (an arc on a low edge, u = 0 or v = 0, adds nothing) and the trace segments. A clipped polygon lies
   !< inside the cell, so the fluid is on the left of its trace as wound on a high face and on the right on a low face:
   !< the orientation comes from the winding, no sign test. Volume and moments need the high faces only (local
   !< coordinates). With fill, the moments and the wall rule points (written at w0 + their place); with qd, the face
   !< segments of the rules (fsg, nfs) and the rule starts per piece (qst: wall cones qwn, then the segment cones; wst).
   !$acc routine seq
   logical,      intent(in)    :: fill                   !< Fill pass.
   integer(I4P), intent(in)    :: nf                     !< Facets.
   real(R8P),    intent(in)    :: tri(3, 3, nf), nrm(3, nf) !< Triangles, normals.
   integer(I4P), intent(in)    :: nt                     !< Polygons.
   integer(I4P), intent(in)    :: lst(nt)                !< Their facets.
   real(R8P),    intent(in)    :: c0(3), c1(3)           !< Cell corners.
   integer(I4P), intent(in)    :: nm2, nm3               !< Moments per face, per piece.
   integer(I4P), intent(in)    :: e2(2, nm2), e3(3, nm3) !< Exponents.
   logical,      intent(in)    :: qd                     !< Quadrature.
   real(R8P),    intent(in)    :: pv(3, 16, nt)          !< Polygons.
   integer(I4P), intent(in)    :: pfl(16, nt), pn(nt)    !< Vertex flags, counts.
   integer(I4P), intent(in)    :: na, natom              !< Atoms: bound, number.
   integer(I4P), intent(in)    :: parent(na), tgt(na), pid(na) !< Union-find, inclusion targets, piece ids.
   integer(I4P), intent(in)    :: ni, niv                !< Intervals: bound, number.
   integer(I4P), intent(in)    :: ive(ni), iva(ni)       !< Intervals: edge, atom.
   real(R8P),    intent(in)    :: ivl(ni), ivh(ni)       !< Intervals: range.
   real(R8P),    intent(inout) :: wt(11, na), fa(6, na), pm(17, na) !< Wall terms, face areas (per root), measures.
   real(R8P),    intent(inout) :: wvm(nm3, na), wwm(nm3, na), fmm(nm2, 6, na) !< Moments (per root).
   integer(I4P), intent(in)    :: ns                     !< Face segments (bound).
   integer(I4P), intent(out)   :: nfs                    !< Face segments.
   real(R8P),    intent(inout) :: fsg(6, ns)             !< Face segments.
   integer(I4P), intent(inout) :: fsf(ns), fsr(ns)       !< Their faces and sets.
   integer(I4P), intent(in)    :: nn                     !< Pieces.
   integer(I4P), intent(inout) :: qst(*), wst(*), qwn(*), fw(*) !< Rule starts, wall cones, wall fill pointers (qd).
   integer(I4P), intent(in)    :: w0, nwo                !< Wall rule points before the cell's, size of the output.
   real(R8P),    intent(inout) :: wpt(3, nwo), wwt(nwo)  !< Wall rule.
   real(R8P)                   :: cc(3), ts(3, 3), tp(3), mtmp(MONO2_MAX), h(3), pl(3), ql(3), sg, vp, mp(3), hf(3) !< Work.
   real(R8P)                   :: d1(3), d2(3), d3(3), e1(3), e2v(3), ga(2), gb(2), tp7(3, 7), tw7(7), sa(3), sb(3), s !< Work.
   integer(I4P)                :: qq, vv, jj, r, rr, ax, sd, bt, b1, b2, kk, ee, k1, k2, s1, s2, gv, kq, l, a !< Work.
   logical                     :: inplane, outside       !< Polygon in the face plane, its outward normal leaves the cell.

   h = c1 - c0
   do r = 1, natom
      wt(:, r) = 0._R8P ; fa(:, r) = 0._R8P
      if (nm3 > 0_I4P) then
         wvm(:, r) = 0._R8P ; wwm(:, r) = 0._R8P ; fmm(:, :, r) = 0._R8P
      endif
   enddo
   cc = c0 + 0.5_R8P * h
   nfs = 0_I4P
   if (qd) then ! wall cones (36 per fan triangle) and wall rule points (7) per piece, before the segments
      do kq = 1, nn + 1
         qst(kq) = 0_I4P ; wst(kq) = 0_I4P
      enddo
      do kq = 1, nn
         qwn(kq) = 0_I4P
      enddo
      do qq = 1, nt
         if (pn(qq) == 0_I4P) cycle
         rr = uf_find(na, parent, qq) ; if (tgt(rr) /= 0_I4P) rr = tgt(rr)
         if (pid(rr) /= 0_I4P) then
            qst(pid(rr) + 1) = qst(pid(rr) + 1) + 36_I4P * (pn(qq) - 2_I4P)
            qwn(pid(rr)) = qwn(pid(rr)) + 36_I4P * (pn(qq) - 2_I4P)
         endif
         if (tgt(rr) /= 0_I4P) rr = tgt(rr) ! the wall points of an inclusion follow its target once more (as before #100)
         if (pid(rr) /= 0_I4P) wst(pid(rr) + 1) = wst(pid(rr) + 1) + 7_I4P * (pn(qq) - 2_I4P)
      enddo
      wst(1) = 1_I4P
      do kq = 1, nn
         wst(kq + 1) = wst(kq + 1) + wst(kq) ; fw(kq) = wst(kq)
      enddo
   endif
   ! walls and traces
   do qq = 1, nt
      if (pn(qq) == 0_I4P) cycle
      r = uf_find(na, parent, qq) ; if (tgt(r) /= 0_I4P) r = tgt(r)
      do vv = 2, pn(qq) - 1
         d1 = pv(:, 1, qq) - c0 ; d2 = pv(:, vv, qq) - c0 ; d3 = pv(:, vv + 1, qq) - c0
         call tri_terms(d1, d2, d3, wt(:, r))
         if (qd .and. fill) then ! wall rule on the fan triangle
            rr = r ; if (tgt(rr) /= 0_I4P) rr = tgt(rr)
            if (pid(rr) /= 0_I4P) then
               call triangle_rule(pv(:, 1, qq), pv(:, vv, qq), pv(:, vv + 1, qq), tp7, tw7, .false.)
               kq = pid(rr)
               do gv = 1, 7
                  wpt(:, w0 + fw(kq)) = tp7(:, gv) ; wwt(w0 + fw(kq)) = tw7(gv) ; fw(kq) = fw(kq) + 1_I4P
               enddo
            endif
         endif
         if (fill .and. nm3 > 0_I4P) then
            ts(:, 1) = (pv(:, 1, qq) - cc) / h ; ts(:, 2) = (pv(:, vv, qq) - cc) / h ; ts(:, 3) = (pv(:, vv + 1, qq) - cc) / h
            e1 = pv(:, vv, qq) - pv(:, 1, qq) ; e2v = pv(:, vv + 1, qq) - pv(:, 1, qq)
            call cross3_into(e1, e2v, tp)
            call triangle_moments(ts(:, 1), ts(:, 2), ts(:, 3), e3, h(1) * h(2) * h(3), 0.5_R8P * norm2_3(tp), wvm(:, r), &
                                  wwm(:, r))
         endif
      enddo
      do ax = 1, 3
         b1 = mod(ax, 3_I4P) + 1_I4P ; b2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P
         do sd = 0, 1
            bt = 2_I4P * (ax - 1_I4P) + sd ; sg = merge(1._R8P, -1._R8P, sd == 1_I4P)
            inplane = .true.
            do vv = 1, pn(qq)
               if (.not. btest(pfl(vv, qq), bt)) inplane = .false.
            enddo
            outside = (sd == 1_I4P .and. nrm(ax, lst(qq)) > 0._R8P) .or. (sd == 0_I4P .and. nrm(ax, lst(qq)) < 0._R8P)
            do vv = 1, pn(qq)
               jj = mod(vv, pn(qq)) + 1_I4P
               if (.not. (btest(pfl(vv, qq), bt) .and. btest(pfl(jj, qq), bt))) cycle
               ! traces on the perturbed face (issue #70): the face plane is taken at c + eps (#69)
               if (inplane) then ! the polygon lies in the face plane, just below the perturbed one
                  if (outside) exit ! solid on the cell side: the perturbed face is fluid there, no boundary
                  ! solid beyond the face: its outline bounds the fluid (it stands for the facets around it,
                  ! which cross the perturbed plane), except along the cell edges, where the fluid edge
                  ! intervals carry the boundary
                  if (iand(iand(pfl(vv, qq), pfl(jj, qq)), not(ibset(0_I4P, bt))) /= 0_I4P) cycle
               elseif (max(tri(ax, 1, lst(qq)), tri(ax, 2, lst(qq)), tri(ax, 3, lst(qq))) <= &
                       merge(c1(ax), c0(ax), sd == 1_I4P)) then
                  ! the facet ends on the face plane, which the perturbed plane does not cut: its edge on the plane
                  ! bounds the fluid only if the facet across that edge crosses the plane (that facet's polygon in
                  ! the cell is the edge alone, so this one stands for it)
                  if (.not. across_crosses(lst(qq), ax, merge(c1(ax), c0(ax), sd == 1_I4P), nf, tri, nt, lst)) cycle
               endif
               pl = pv(:, vv, qq) - c0 ; ql = pv(:, jj, qq) - c0
               fa(bt + 1, r) = fa(bt + 1, r) + sg * 0.5_R8P * (pl(b1) + ql(b1)) * (ql(b2) - pl(b2))
               if (qd) then ! the trace, counter-clockwise seen from +e_ax
                  nfs = nfs + 1_I4P ; fsf(nfs) = bt + 1_I4P ; fsr(nfs) = r
                  if (sg > 0._R8P) then
                     fsg(1:3, nfs) = pv(:, vv, qq) ; fsg(4:6, nfs) = pv(:, jj, qq)
                  else
                     fsg(1:3, nfs) = pv(:, jj, qq) ; fsg(4:6, nfs) = pv(:, vv, qq)
                  endif
               endif
               if (fill .and. nm2 > 0_I4P) then
                  mtmp(1:nm2) = 0._R8P
                  ga(1) = pl(b1) / h(b1) - 0.5_R8P ; ga(2) = pl(b2) / h(b2) - 0.5_R8P
                  gb(1) = ql(b1) / h(b1) - 0.5_R8P ; gb(2) = ql(b2) / h(b2) - 0.5_R8P
                  call segment_green_moments(ga, gb, e2, h(b1) * h(b2), mtmp(1:nm2))
                  fmm(:, bt + 1, r) = fmm(:, bt + 1, r) + sg * mtmp(1:nm2)
               endif
            enddo
         enddo
      enddo
   enddo
   ! fluid edge intervals: edge along k at sides s1 (axis k1), s2 (axis k2); on the face (k1, s1) it
   ! runs along v at u = s2 (counter-clockwise: +v at u = h); on the face (k2, s2) along u (no area)
   do jj = 1, niv
      ee = ive(jj) ; kk = (ee - 1_I4P) / 4_I4P + 1_I4P
      s1 = mod(mod(ee - 1_I4P, 4_I4P), 2_I4P) ; s2 = mod(ee - 1_I4P, 4_I4P) / 2_I4P
      k1 = mod(kk, 3_I4P) + 1_I4P ; k2 = mod(kk + 1_I4P, 3_I4P) + 1_I4P
      r = uf_find(na, parent, iva(jj))
      if (s2 == 1_I4P) fa(2 * (k1 - 1) + s1 + 1, r) = fa(2 * (k1 - 1) + s1 + 1, r) + h(k2) * (ivh(jj) - ivl(jj))
      if (fill .and. nm2 > 0_I4P) then ! moments: centred coordinates, so the edge at U = -1/2 counts too (going -V there)
         mtmp(1:nm2) = 0._R8P
         if (s2 == 1_I4P) then
            ga(1) = 0.5_R8P ; ga(2) = (ivl(jj) - cc(kk)) / h(kk) ; gb(1) = 0.5_R8P ; gb(2) = (ivh(jj) - cc(kk)) / h(kk)
         else
            ga(1) = -0.5_R8P ; ga(2) = (ivh(jj) - cc(kk)) / h(kk) ; gb(1) = -0.5_R8P ; gb(2) = (ivl(jj) - cc(kk)) / h(kk)
         endif
         call segment_green_moments(ga, gb, e2, h(k2) * h(kk), mtmp(1:nm2))
         fmm(:, 2 * (k1 - 1) + s1 + 1, r) = fmm(:, 2 * (k1 - 1) + s1 + 1, r) + mtmp(1:nm2)
      endif
      if (qd) then ! the interval bounds two faces: (k1, s1) along V at U = s2, (k2, s2) along U at V = s1
         sa = 0._R8P ; sa(kk) = ivl(jj) ; sa(k1) = merge(c1(k1), c0(k1), s1 == 1_I4P)
         sa(k2) = merge(c1(k2), c0(k2), s2 == 1_I4P) ; sb = sa ; sb(kk) = ivh(jj)
         nfs = nfs + 1_I4P ; fsf(nfs) = 2_I4P * (k1 - 1_I4P) + s1 + 1_I4P ; fsr(nfs) = r
         if (s2 == 1_I4P) then
            fsg(1:3, nfs) = sa ; fsg(4:6, nfs) = sb
         else
            fsg(1:3, nfs) = sb ; fsg(4:6, nfs) = sa
         endif
         nfs = nfs + 1_I4P ; fsf(nfs) = 2_I4P * (k2 - 1_I4P) + s2 + 1_I4P ; fsr(nfs) = r
         if (s1 == 0_I4P) then
            fsg(1:3, nfs) = sa ; fsg(4:6, nfs) = sb
         else
            fsg(1:3, nfs) = sb ; fsg(4:6, nfs) = sa
         endif
      endif
   enddo
   if (qd) then ! the segment cones per piece, then the starts
      do l = 1, nfs
         rr = uf_find(na, parent, fsr(l)) ; if (pid(rr) /= 0_I4P) qst(pid(rr) + 1) = qst(pid(rr) + 1) + 36_I4P
      enddo
      qst(1) = 1_I4P
      do kq = 1, nn
         qst(kq + 1) = qst(kq + 1) + qst(kq)
      enddo
   endif
   hf(1) = h(2) * h(3) ; hf(2) = h(3) * h(1) ; hf(3) = h(1) * h(2)
   do r = 1, natom
      if (pid(r) == 0_I4P) cycle
      kk = pid(r)
      s = 0._R8P
      do a = 1, 3 ! sum(h * fa([2, 4, 6], r)), in its order
         s = s + h(a) * fa(2 * a, r)
      enddo
      vp = s / 3._R8P - wt(1, r)
      do a = 1, 3
         mp(a) = 0.5_R8P * h(a) * h(a) * fa(2 * a, r) - wt(1 + a, r)
      enddo
      pm(1, kk) = vp
      pm(2:4, kk) = c0 + 0.5_R8P * h
      if (vp > 0._R8P) pm(2:4, kk) = c0 + mp / vp
      pm(5, kk) = wt(5, r)
      pm(6:8, kk) = c0 + 0.5_R8P * h
      if (wt(5, r) > 0._R8P) pm(6:8, kk) = c0 + wt(6:8, r) / wt(5, r)
      pm(9:11, kk) = wt(9:11, r)
      do a = 1, 6
         pm(11 + a, kk) = fa(a, r) / hf((a + 1) / 2)
      enddo
   enddo
   endsubroutine piece_measures

   pure subroutine piece_rules(multi, nt, c0, c1, pv, pn, na, parent, tgt, pid, nn, ns, nfs, fsg, fsf, fsr, qst, wst, fq, &
                               p0, q0, nqo, qpt, qwt, nqcd, qcnt, wcnt)
   !< Volume rules of the pieces (issue #70): cones from the cell centre over the fluid boundary (the wall fans reversed:
   !< their normal is the solid's; on every face the triangles from the face centre to the oriented segments), 36 points
   !< per cone (`tet_rule`), written at q0 + the piece's place (qst). A single piece takes its wall cones only (its faces:
   !< the plane sections, by the caller), several their segment cones too and their point counts.
   !$acc routine seq
   logical,      intent(in)    :: multi                  !< Several pieces.
   integer(I4P), intent(in)    :: nt                     !< Polygons.
   real(R8P),    intent(in)    :: c0(3), c1(3)           !< Cell corners.
   real(R8P),    intent(in)    :: pv(3, 16, nt)          !< Polygons.
   integer(I4P), intent(in)    :: pn(nt)                 !< Vertex counts.
   integer(I4P), intent(in)    :: na                     !< Atoms.
   integer(I4P), intent(in)    :: parent(na), tgt(na), pid(na) !< Union-find, inclusion targets, piece ids.
   integer(I4P), intent(in)    :: nn, ns, nfs            !< Pieces, segments: bound, number.
   real(R8P),    intent(in)    :: fsg(6, ns)             !< Face segments.
   integer(I4P), intent(in)    :: fsf(ns), fsr(ns)       !< Their faces and sets.
   integer(I4P), intent(in)    :: qst(nn + 1), wst(nn + 1) !< Rule starts.
   integer(I4P), intent(inout) :: fq(nn)                 !< Fill pointers.
   integer(I4P), intent(in)    :: p0, q0, nqo, nqcd      !< Pieces, rule points before the cell's; output sizes.
   real(R8P),    intent(inout) :: qpt(3, nqo), qwt(nqo)  !< Volume rules.
   integer(I4P), intent(inout) :: qcnt(nqcd), wcnt(nqcd) !< Rule points per piece.
   real(R8P)                   :: cc(3), pf(3), rp(3, 36), rw(36) !< Work.
   integer(I4P)                :: qq, vv, rr, kq, m, ft, ax_, sd_ !< Work.

   cc = c0 + 0.5_R8P * (c1 - c0)
   do kq = 1, nn
      fq(kq) = qst(kq)
   enddo
   do qq = 1, nt
      if (pn(qq) == 0_I4P) cycle
      rr = uf_find(na, parent, qq) ; if (tgt(rr) /= 0_I4P) rr = tgt(rr)
      if (pid(rr) == 0_I4P) cycle
      kq = pid(rr)
      do vv = 2, pn(qq) - 1 ! reversed: the fluid's outward normal
         call tet_rule(cc, pv(:, 1, qq), pv(:, vv + 1, qq), pv(:, vv, qq), rp, rw)
         qpt(:, q0 + fq(kq):q0 + fq(kq) + 35) = rp ; qwt(q0 + fq(kq):q0 + fq(kq) + 35) = rw ; fq(kq) = fq(kq) + 36_I4P
      enddo
   enddo
   if (.not. multi) return
   do m = 1, nfs
      rr = uf_find(na, parent, fsr(m)) ; if (pid(rr) == 0_I4P) cycle
      kq = pid(rr) ; ft = fsf(m) ; ax_ = (ft + 1_I4P) / 2_I4P ; sd_ = mod(ft + 1_I4P, 2_I4P)
      pf = cc ; pf(ax_) = merge(c1(ax_), c0(ax_), sd_ == 1_I4P)
      if (sd_ == 1_I4P) then ! high face: outward +e_ax, as the segments turn
         call tet_rule(cc, pf, fsg(1:3, m), fsg(4:6, m), rp, rw)
      else
         call tet_rule(cc, pf, fsg(4:6, m), fsg(1:3, m), rp, rw)
      endif
      qpt(:, q0 + fq(kq):q0 + fq(kq) + 35) = rp ; qwt(q0 + fq(kq):q0 + fq(kq) + 35) = rw ; fq(kq) = fq(kq) + 36_I4P
   enddo
   do kq = 1, nn
      qcnt(p0 + kq) = qst(kq + 1) - qst(kq) ; wcnt(p0 + kq) = wst(kq + 1) - wst(kq)
   enddo
   endsubroutine piece_rules

   pure subroutine sort_pierce(n, pe, pt, pen, pp)
   !< Sort piercings by edge, then position (insertion, stable).
   !$acc routine seq
   integer(I4P), intent(in)    :: n                  !< Piercings.
   integer(I4P), intent(inout) :: pe(n), pen(n), pp(n) !< Edge, entering, polygon.
   real(R8P),    intent(inout) :: pt(n)              !< Position.
   integer(I4P)                :: i, j, ve, ven, vp  !< Work.
   real(R8P)                   :: vt                 !< Work.

   do i = 2, n
      ve = pe(i) ; vt = pt(i) ; ven = pen(i) ; vp = pp(i) ; j = i - 1_I4P
      do while (j >= 1_I4P)
         if (pe(j) < ve .or. (pe(j) == ve .and. pt(j) <= vt)) exit
         pe(j + 1) = pe(j) ; pt(j + 1) = pt(j) ; pen(j + 1) = pen(j) ; pp(j + 1) = pp(j) ; j = j - 1_I4P
      enddo
      pe(j + 1) = ve ; pt(j + 1) = vt ; pen(j + 1) = ven ; pp(j + 1) = vp
   enddo
   endsubroutine sort_pierce

   pure subroutine sort_keys(nr, m, key, idx)
   !< Index heapsort of the keys key(:, 1:m), lexicographic over the nr rows.
   !$acc routine seq
   integer(I4P), intent(in)  :: nr, m       !< Rows, keys.
   real(R8P),    intent(in)  :: key(nr, m)  !< Keys.
   integer(I4P), intent(out) :: idx(m)      !< Sorted order.
   integer(I4P)              :: i, t        !< Work.

   do i = 1, m
      idx(i) = i
   enddo
   do i = m / 2_I4P, 1_I4P, -1_I4P
      call sift_keys(nr, m, key, idx, i, m)
   enddo
   do i = m, 2_I4P, -1_I4P
      t = idx(1) ; idx(1) = idx(i) ; idx(i) = t
      call sift_keys(nr, m, key, idx, 1_I4P, i - 1_I4P)
   enddo
   endsubroutine sort_keys

   pure subroutine sift_keys(nr, m, key, idx, lo, hi)
   !< Sift down idx(lo) in the heap idx(1:hi).
   !$acc routine seq
   integer(I4P), intent(in)    :: nr, m      !< Rows, keys.
   real(R8P),    intent(in)    :: key(nr, m) !< Keys.
   integer(I4P), intent(inout) :: idx(m)     !< Heap.
   integer(I4P), intent(in)    :: lo, hi     !< Root, heap end.
   integer(I4P)                :: r, c, tt   !< Work.

   r = lo
   do while (2_I4P * r <= hi)
      c = 2_I4P * r
      if (c < hi) then
         if (key_less(nr, m, key, idx(c), idx(c + 1))) c = c + 1_I4P
      endif
      if (.not. key_less(nr, m, key, idx(r), idx(c))) exit
      tt = idx(r) ; idx(r) = idx(c) ; idx(c) = tt ; r = c
   enddo
   endsubroutine sift_keys

   pure logical function key_less(nr, m, key, a, b) result(less)
   !< Lexicographic order of the keys a, b.
   !$acc routine seq
   integer(I4P), intent(in) :: nr, m      !< Rows, keys.
   real(R8P),    intent(in) :: key(nr, m) !< Keys.
   integer(I4P), intent(in) :: a, b       !< Keys compared.
   integer(I4P)             :: r          !< Row.

   less = .false.
   do r = 1, nr
      if (key(r, a) /= key(r, b)) then
         less = key(r, a) < key(r, b)
         return
      endif
   enddo
   endfunction key_less

   pure subroutine cross3_into(a, b, c)
   !< Cross product c = a x b (the formula of `cross3`).
   !$acc routine seq
   real(R8P), intent(in)  :: a(3), b(3) !< Vectors.
   real(R8P), intent(out) :: c(3)       !< a x b.

   c(1) = a(2) * b(3) - a(3) * b(2) ; c(2) = a(3) * b(1) - a(1) * b(3) ; c(3) = a(1) * b(2) - a(2) * b(1)
   endsubroutine cross3_into

   pure function norm2_3(v) result(r)
   !< norm2(v) of a 3-vector by the scaled algorithm of gfortran's norm2, written out (issue #100; the 3-D twin of
   !< `norm2_xy` in fossil_cut_cells_2d: nvfortran's device and host norm2 differ by 1 ulp; this one gives the bits of
   !< gfortran's intrinsic, checked on 4 million vectors from 1e-180 to 1e170, zeros included).
   !$acc routine seq
   real(R8P), intent(in) :: v(3)               !< Vector.
   real(R8P)             :: r                  !< Norm.
   real(R8P)             :: scale, res, a, val !< Work.
   integer(I4P)          :: i                  !< Work.

   scale = 1._R8P ; res = 0._R8P
   do i = 1, 3
      if (v(i) /= 0._R8P) then
         a = abs(v(i))
         if (a > scale) then
            val = scale / a ; res = (val * val) * res + 1._R8P ; scale = a
         else
            val = a / scale ; res = val * val + res
         endif
      endif
   enddo
   r = sqrt(res) * scale
   endfunction norm2_3

   pure subroutine volume_moments_k(nm3, nm2, wv, fm, fp, hx, e2, e3, m)
   !< Volume moments of a piece (issue #70): int X^a Y^b Z^c dV = oint h_x X^(a+1) / (a+1) Y^b Z^c n_x dA over its
   !< boundary: minus the wall terms wv (taken with the solid's normal) plus the faces x = -+ h_x / 2 (fluid moments fm,
   !< fp in (Y, Z), the outward normal -+1).
   !$acc routine seq
   integer(I4P), intent(in)  :: nm3, nm2      !< Moments per piece, per face.
   real(R8P),    intent(in)  :: wv(nm3)       !< Wall divergence terms (solid normal).
   real(R8P),    intent(in)  :: fm(nm2), fp(nm2) !< Moments of the -x and +x faces.
   real(R8P),    intent(in)  :: hx            !< Cell width along x.
   integer(I4P), intent(in)  :: e2(2, nm2), e3(3, nm3) !< Exponents.
   real(R8P),    intent(out) :: m(nm3)        !< Moments.
   integer(I4P)              :: k, l          !< Work.

   do k = 1, nm3
      m(k) = -wv(k)
      do l = 1, nm2
         if (e2(1, l) /= e3(2, k) .or. e2(2, l) /= e3(3, k)) cycle
         m(k) = m(k) + hx / (e3(1, k) + 1) * (0.5_R8P ** (e3(1, k) + 1) * fp(l) - (-0.5_R8P) ** (e3(1, k) + 1) * fm(l))
      enddo
   enddo
   endsubroutine volume_moments_k

   pure logical function split_closes(np, pm, h, alpha, faces) result(ok)
   !< The pieces pm [17, np] of a cell of size h split it consistently (issue #70): their face fractions add up to the
   !< exact face totals (the plane sections), their volumes to alpha |cell|, and every piece closes, (f+ - f-) A_a = the
   !< a component of its wall normal integral, per axis. Tolerance 1e-12 of the face area or cell volume, as the other
   !< trust tests. A face split by the traces of the clipped polygons can miss a boundary on the perturbed face plane
   !< (a facet lying in it or ending on it next to one crossing it): such a cell is not trusted.
   !$acc routine seq
   integer(I4P), intent(in) :: np        !< Pieces.
   real(R8P),    intent(in) :: pm(17, np) !< Piece measures (see `piece_measures`).
   real(R8P),    intent(in) :: h(3)      !< Cell size.
   real(R8P),    intent(in) :: alpha     !< Cell volume fraction.
   real(R8P),    intent(in) :: faces(6)  !< Face fractions (-x, +x, -y, +y, -z, +z).
   real(R8P)                :: af(3), vc, s !< Face areas, cell volume, sum.
   integer(I4P)             :: a, p      !< Work.
   real(R8P), parameter     :: TOL = 1.e-12_R8P !< Relative tolerance.

   af(1) = h(2) * h(3) ; af(2) = h(3) * h(1) ; af(3) = h(1) * h(2) ; vc = h(1) * h(2) * h(3)
   s = 0._R8P
   do p = 1, np
      s = s + pm(1, p)
   enddo
   ok = abs(s - alpha * vc) <= TOL * vc
   do a = 1, 6
      s = 0._R8P
      do p = 1, np
         s = s + pm(11 + a, p)
      enddo
      ok = ok .and. abs(s - faces(a)) <= TOL
   enddo
   do p = 1, np
      do a = 1, 3
         ok = ok .and. abs((pm(13 + 2 * (a - 1), p) - pm(12 + 2 * (a - 1), p)) * af(a) - pm(8 + a, p)) <= TOL * af(a)
      enddo
   enddo
   endfunction split_closes

   pure subroutine single_piece_3d(hx, hy, hz, alpha, centroid, wet_area, wet_centroid, wet_normal, fx0, fx1, fy0, fy1, &
                                   fz0, fz1, pvol, pcen, pwa, pwc, pwn, pfc)
   !< Kernel (issue #100): the measures of the single piece of a cell: the cell totals (exact).
   !$acc routine seq
   real(R8P), intent(in)  :: hx, hy, hz                   !< Cell size.
   real(R8P), intent(in)  :: alpha, centroid(3), wet_area, wet_centroid(3), wet_normal(3) !< Cell totals.
   real(R8P), intent(in)  :: fx0, fx1, fy0, fy1, fz0, fz1 !< Face fractions.
   real(R8P), intent(out) :: pvol, pcen(3), pwa, pwc(3), pwn(3), pfc(6) !< Piece measures.

   pvol = alpha * (hx * hy * hz)
   pcen = centroid ; pwa = wet_area ; pwc = wet_centroid ; pwn = wet_normal
   pfc(1) = fx0 ; pfc(2) = fx1 ; pfc(3) = fy0 ; pfc(4) = fy1 ; pfc(5) = fz0 ; pfc(6) = fz1
   endsubroutine single_piece_3d

   subroutine section_face_moments(cut, i, j, k, nm2, secm)
   !< The fluid moments of the six faces of cell (i, j, k) from the plane sections (-x, +x, -y, +y, -z, +z).
   type(cut_cells_3d_t), intent(in)  :: cut          !< Cut cells.
   integer(I4P),         intent(in)  :: i, j, k, nm2 !< Cell, moments per face.
   real(R8P),            intent(out) :: secm(nm2, 6) !< Moments.

   secm(:, 1) = cut%face_x_moments(:, i, j, k) ; secm(:, 2) = cut%face_x_moments(:, i + 1, j, k)
   secm(:, 3) = cut%face_y_moments(:, i, j, k) ; secm(:, 4) = cut%face_y_moments(:, i, j + 1, k)
   secm(:, 5) = cut%face_z_moments(:, i, j, k) ; secm(:, 6) = cut%face_z_moments(:, i, j, k + 1)
   endsubroutine section_face_moments

   integer(I4P) function section_cones_count(px, py, pz, i, j, k) result(n)
   !< Points of the cones over the exact section faces of cell (i, j, k) (see `face_cones`): 36 per boundary segment.
   type(cut_cells_2d_t), intent(in) :: px(:), py(:), pz(:) !< Plane sections.
   integer(I4P),         intent(in) :: i, j, k             !< Cell.
   integer(I4P)                     :: f, ic(3), ax, b1, b2, iu, iv, ip, pl !< Work.

   n = 0_I4P ; ic = [i, j, k]
   do f = 1, 6
      ax = (f + 1_I4P) / 2_I4P ; b1 = mod(ax, 3_I4P) + 1_I4P ; b2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P
      iu = ic(b1) ; iv = ic(b2) ; pl = ic(ax) + mod(f + 1_I4P, 2_I4P)
      select case (ax)
      case (1)
         if (px(pl)%pieces(iu, iv) < 1_I4P) cycle
         do ip = px(pl)%first_piece(iu, iv), px(pl)%first_piece(iu, iv) + px(pl)%pieces(iu, iv) - 1_I4P
            n = n + 36_I4P * (px(pl)%piece_seg_start(ip + 1) - px(pl)%piece_seg_start(ip))
         enddo
      case (2)
         if (py(pl)%pieces(iu, iv) < 1_I4P) cycle
         do ip = py(pl)%first_piece(iu, iv), py(pl)%first_piece(iu, iv) + py(pl)%pieces(iu, iv) - 1_I4P
            n = n + 36_I4P * (py(pl)%piece_seg_start(ip + 1) - py(pl)%piece_seg_start(ip))
         enddo
      case default
         if (pz(pl)%pieces(iu, iv) < 1_I4P) cycle
         do ip = pz(pl)%first_piece(iu, iv), pz(pl)%first_piece(iu, iv) + pz(pl)%pieces(iu, iv) - 1_I4P
            n = n + 36_I4P * (pz(pl)%piece_seg_start(ip + 1) - pz(pl)%piece_seg_start(ip))
         enddo
      endselect
   enddo
   endfunction section_cones_count

   subroutine face_cones(sec, x, y, z, ic, f, pts, ws, n)
   !< Cones from the centre of cell ic over the fluid region of its face f, from the oriented boundary of the plane
   !< section sec (exact, as the face fractions): triangles from the face centre (issue #70), written at n + 1.. (36
   !< points per boundary segment; the caller sized the rule with `section_cones_count`).
   type(cut_cells_2d_t), intent(in)    :: sec              !< Plane section.
   real(R8P),            intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   integer(I4P),         intent(in)    :: ic(3), f         !< Cell, face (1-6).
   real(R8P),            intent(inout) :: pts(:,:), ws(:)  !< Points, weights.
   integer(I4P),         intent(inout) :: n                !< Points so far.
   real(R8P)    :: c0(3), c1(3), cc(3), pf(3), a3(3), b3(3), rp(3, 36), rw(36) !< Work.
   integer(I4P) :: ax, sd, ip, kk, iu, iv, b1, b2      !< Work.

   ax = (f + 1_I4P) / 2_I4P ; sd = mod(f + 1_I4P, 2_I4P)
   b1 = mod(ax, 3_I4P) + 1_I4P ; b2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P ; iu = ic(b1) ; iv = ic(b2)
   c0 = [x(ic(1)), y(ic(2)), z(ic(3))] ; c1 = [x(ic(1) + 1), y(ic(2) + 1), z(ic(3) + 1)] ; cc = 0.5_R8P * (c0 + c1)
   pf = cc ; pf(ax) = merge(c1(ax), c0(ax), sd == 1_I4P)
   if (sec%pieces(iu, iv) < 1_I4P) return
   do ip = sec%first_piece(iu, iv), sec%first_piece(iu, iv) + sec%pieces(iu, iv) - 1_I4P
      do kk = sec%piece_seg_start(ip), sec%piece_seg_start(ip + 1) - 1_I4P
         a3(ax) = pf(ax) ; a3(b1) = sec%piece_seg(1, kk) ; a3(b2) = sec%piece_seg(2, kk)
         b3(ax) = pf(ax) ; b3(b1) = sec%piece_seg(3, kk) ; b3(b2) = sec%piece_seg(4, kk)
         if (sd == 1_I4P) then ! high face: outward +e_ax, counter-clockwise in (u, v)
            call tet_rule(cc, pf, a3, b3, rp, rw)
         else
            call tet_rule(cc, pf, b3, a3, rp, rw)
         endif
         pts(:, n + 1:n + 36) = rp ; ws(n + 1:n + 36) = rw ; n = n + 36_I4P
      enddo
   enddo
   endsubroutine face_cones

   subroutine links_3d(cut, x, y, z, px, py, pz, piv_piece, piv_edge, piv_lo, piv_hi, links)
   !< Links between the pieces of the block (issue #70; `fossil_cut_graph`). A connected fluid region of a face lies
   !< next to one piece on each side; the regions are the pieces of the plane section's 2-D cut cell (areas from the
   !< section, the same numbers as the face fractions), matched to the 3-D pieces through the fluid edge intervals the
   !< two share (the edges of the face are cell edges). A one-piece side takes every region. A region with no fluid
   !< interval (a fluid island inside a solid trace) goes to the only piece of that side whose face area it leaves
   !< unexplained; otherwise it is counted in `unmatched` and gets no link.
   type(cut_cells_3d_t), intent(in)  :: cut              !< Cut cells (pieces).
   real(R8P),            intent(in)  :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_2d_t), intent(in)  :: px(:), py(:), pz(:) !< Plane sections.
   integer(I4P),         intent(in)  :: piv_piece(:), piv_edge(:) !< Fluid edge intervals of multi-piece cells.
   real(R8P),            intent(in)  :: piv_lo(:), piv_hi(:)      !< Their ranges.
   type(cut_links_t),    intent(out) :: links            !< Links.
   integer(I4P), allocatable :: lp(:), lq(:), lf(:), ivs(:), ivo(:) !< Raw links; intervals by piece (CSR).
   real(R8P),    allocatable :: la(:)                    !< Raw link areas.
   integer(I4P) :: nl, a, i, ic(3), jc(3), np2, ua, va, n(3) !< Work.
   real(R8P)    :: af, fr                                !< Face area, face fraction.
   real(R8P)    :: e(3, 2)                               !< Low and high edges of the cell.

   links%np = cut%np ; n = [cut%nx, cut%ny, cut%nz]
   ! intervals by piece
   allocate(ivs(cut%np + 1), ivo(size(piv_piece))) ; ivs = 0_I4P
   do i = 1, size(piv_piece)
      ivs(piv_piece(i) + 1) = ivs(piv_piece(i) + 1) + 1_I4P
   enddo
   ivs(1) = 1_I4P
   do i = 1, cut%np
      ivs(i + 1) = ivs(i + 1) + ivs(i)
   enddo
   allocate(lp(64), lq(64), lf(64), la(64)) ; nl = 0_I4P
   block
      integer(I4P) :: fill(cut%np)
      fill = 0_I4P
      do i = 1, size(piv_piece)
         ivo(ivs(piv_piece(i)) + fill(piv_piece(i))) = i ; fill(piv_piece(i)) = fill(piv_piece(i)) + 1_I4P
      enddo
   endblock
   do a = 1, 3
      ua = mod(a, 3_I4P) + 1_I4P ; va = mod(a + 1_I4P, 3_I4P) + 1_I4P
      do np2 = 2, n(a) ! interior planes: the face between cells np2 - 1 and np2 along a
         select case (a)
         case (1) ; call plane_links(px(np2))
         case (2) ; call plane_links(py(np2))
         case default ; call plane_links(pz(np2))
         endselect
      enddo
   enddo
   call build_csr(nl, lp, lq, lf, la, links)
   contains
      subroutine plane_links(sc)
      !< Links over the faces of one interior grid plane (axis a, index np2), with its section sc.
      type(cut_cells_2d_t), intent(in) :: sc !< Plane section.
      integer(I4P) :: i1, i2, rr, m, nreg   !< Work.
      integer(I4P), allocatable :: pr(:,:)  !< Pieces of each region on the two sides.
      real(R8P),    allocatable :: ar(:)    !< Region areas.

      do i2 = 1, n(va)
         do i1 = 1, n(ua)
            ic(a) = np2 - 1_I4P ; ic(ua) = i1 ; ic(va) = i2 ; jc = ic ; jc(a) = np2
            if (cut%pieces(ic(1), ic(2), ic(3)) < 1_I4P .or. cut%pieces(jc(1), jc(2), jc(3)) < 1_I4P) cycle
            select case (a)
            case (1) ; fr = cut%face_x(jc(1), jc(2), jc(3))
            case (2) ; fr = cut%face_y(jc(1), jc(2), jc(3))
            case default ; fr = cut%face_z(jc(1), jc(2), jc(3))
            endselect
            if (fr <= 0._R8P) cycle
            e(:, 1) = [x(ic(1)), y(ic(2)), z(ic(3))] ; e(:, 2) = [x(ic(1) + 1), y(ic(2) + 1), z(ic(3) + 1)]
            af = (e(ua, 2) - e(ua, 1)) * (e(va, 2) - e(va, 1))
            if (cut%pieces(ic(1), ic(2), ic(3)) == 1_I4P .and. cut%pieces(jc(1), jc(2), jc(3)) == 1_I4P) then
               call add_link(cut%first_piece(ic(1), ic(2), ic(3)), cut%first_piece(jc(1), jc(2), jc(3)), fr * af)
               cycle
            endif
            if (sc%pieces(i1, i2) < 1_I4P) cycle
            nreg = sc%pieces(i1, i2)
            allocate(pr(2, nreg), ar(nreg))
            do rr = 1, nreg
               m = sc%first_piece(i1, i2) + rr - 1_I4P
               ar(rr) = sc%piece_area(m)
               pr(1, rr) = side_piece(sc, m, ic, 1_I4P) ! the low cell: its +a face
               pr(2, rr) = side_piece(sc, m, jc, 0_I4P) ! the high cell: its -a face
            enddo
            call resolve(pr, ar)
            do rr = 1, nreg
               if (ar(rr) <= 0._R8P) cycle
               if (pr(1, rr) == 0_I4P .or. pr(2, rr) == 0_I4P) then
                  links%unmatched = links%unmatched + 1_I4P ; cycle
               endif
               call add_link(pr(1, rr), pr(2, rr), ar(rr))
            enddo
            deallocate(pr, ar)
         enddo
      enddo
      endsubroutine plane_links

      integer(I4P) function side_piece(sc, m, cc, sd) result(pq)
      !< The piece of cell cc whose fluid edge intervals overlap those of region m of the section (0: none). A
      !< one-piece cell gives its piece.
      type(cut_cells_2d_t), intent(in) :: sc           !< Section.
      integer(I4P),         intent(in) :: m, cc(3), sd !< Region, cell, side of the face (axis a) in cc.
      integer(I4P) :: t, ed, ll, kax, sides(3), le, nu, nv, w, pw !< Work.
      real(R8P)    :: tol                              !< Overlap tolerance.

      pq = 0_I4P
      if (cut%pieces(cc(1), cc(2), cc(3)) == 1_I4P) then
         pq = cut%first_piece(cc(1), cc(2), cc(3)) ; return
      endif
      nu = sc%nx ; nv = sc%ny
      do t = sc%piece_iv_start(m), sc%piece_iv_start(m + 1) - 1
         ed = sc%piece_iv_edge(t) ; sides = 0_I4P ; sides(a) = sd
         if (ed <= (nu + 1_I4P) * nv) then ! along v, at u line ll
            ll = mod(ed - 1_I4P, nu + 1_I4P) + 1_I4P ; kax = va
            sides(ua) = merge(0_I4P, 1_I4P, ll == cc(ua))
         else ! along u, at v line ll
            ll = (ed - (nu + 1_I4P) * nv - 1_I4P) / nu + 1_I4P ; kax = ua
            sides(va) = merge(0_I4P, 1_I4P, ll == cc(va))
         endif
         le = edge_index(kax, sides)
         tol = 1.e-12_R8P * (sc%piece_iv_hi(t) - sc%piece_iv_lo(t) + abs(sc%piece_iv_hi(t)))
         do pw = cut%first_piece(cc(1), cc(2), cc(3)), cut%first_piece(cc(1), cc(2), cc(3)) + cut%pieces(cc(1), cc(2), cc(3)) - 1
            do w = ivs(pw), ivs(pw + 1) - 1
               if (piv_edge(ivo(w)) /= le) cycle
               if (min(piv_hi(ivo(w)), sc%piece_iv_hi(t)) - max(piv_lo(ivo(w)), sc%piece_iv_lo(t)) > tol) then
                  pq = pw ; return
               endif
            enddo
         enddo
      enddo
      endfunction side_piece

      subroutine resolve(pr, ar)
      !< Regions left without a piece on a side: the only piece of that side whose face area is not explained by its
      !< matched regions takes them, one at a time.
      integer(I4P), intent(inout) :: pr(:,:)       !< Pieces of the regions.
      real(R8P),    intent(in)    :: ar(:)         !< Region areas.
      integer(I4P) :: side, rr, pw, cand, ncand, cc(3), fidx !< Work.
      real(R8P)    :: res, afc                     !< Residual, face area.

      do side = 1, 2
         if (side == 1) then
            cc = ic ; fidx = 2_I4P * a
         else
            cc = jc ; fidx = 2_I4P * a - 1_I4P
         endif
         afc = af
         do rr = 1, size(ar)
            if (pr(side, rr) /= 0_I4P .or. ar(rr) <= 0._R8P) cycle
            ncand = 0_I4P ; cand = 0_I4P
            do pw = cut%first_piece(cc(1), cc(2), cc(3)), cut%first_piece(cc(1), cc(2), cc(3)) + cut%pieces(cc(1), cc(2), cc(3)) - 1
               res = cut%piece_face(fidx, pw) * afc - sum(ar, mask=pr(side, :) == pw)
               if (res > 1.e-12_R8P * afc) then
                  ncand = ncand + 1_I4P ; cand = pw
               endif
            enddo
            if (ncand == 1_I4P) pr(side, rr) = cand
         enddo
      enddo
      endsubroutine resolve

      subroutine add_link(pa, pb, ar_)
      !< The two directed links of piece pa (low cell, face +a) and pb (high cell, face -a), area ar_.
      integer(I4P), intent(in) :: pa, pb !< Pieces.
      real(R8P),    intent(in) :: ar_    !< Area.

      if (nl + 2_I4P > size(lp)) then
         lp = [lp, lp] ; lq = [lq, lq] ; lf = [lf, lf] ; la = [la, la]
      endif
      nl = nl + 1_I4P ; lp(nl) = pa ; lq(nl) = pb ; lf(nl) = 2_I4P * a ; la(nl) = ar_
      nl = nl + 1_I4P ; lp(nl) = pb ; lq(nl) = pa ; lf(nl) = 2_I4P * a - 1_I4P ; la(nl) = ar_
      endsubroutine add_link
   endsubroutine links_3d

   function section_moments(c2, iu, iv, nm2) result(m)
   !< Fluid moments of the cell (iu, iv) of a plane section: the sum over its pieces (issue #70).
   type(cut_cells_2d_t), intent(in) :: c2      !< Section.
   integer(I4P),         intent(in) :: iu, iv  !< Cell.
   integer(I4P),         intent(in) :: nm2     !< Moments.
   real(R8P)                        :: m(nm2)  !< Moments.
   integer(I4P)                     :: p       !< Work.

   m = 0._R8P
   if (c2%pieces(iu, iv) < 1_I4P) return
   do p = c2%first_piece(iu, iv), c2%first_piece(iu, iv) + c2%pieces(iu, iv) - 1_I4P
      m = m + c2%piece_moments(:, p)
   enddo
   endfunction section_moments

   subroutine edge_fractions(surface, x, y, z, cut, all_edges)
   !< Inside flags of the grid nodes and fluid fractions of the cell edges along z (issue #70), from the exact crossings
   !< of the z lines (`line_fractions`); with all_edges, the x and y lines too, for edge_x and edge_y.
   type(surface_stl_object), intent(in)    :: surface          !< Surface.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(inout) :: cut              !< Cut cells.
   logical,                  intent(in)    :: all_edges        !< Edges along x and y too.
   integer(I4P)              :: nx, ny, nz, i, j, k            !< Work.
   real(R8P),    allocatable :: fr(:)                          !< Fractions of one line.
   logical,      allocatable :: nd(:)                          !< Node flags of one line.

   nx = cut%nx ; ny = cut%ny ; nz = cut%nz
   allocate(cut%node_in(nx + 1, ny + 1, nz + 1), cut%edge_z(nx + 1, ny + 1, nz))
   do j = 1, ny + 1 ! z lines at (x, y): every node, and the edges along z
      do i = 1, nx + 1
         call line_fractions(surface, 3_I4P, x(i), y(j), z, fr, nd)
         cut%edge_z(i, j, :) = fr ; cut%node_in(i, j, :) = nd
      enddo
   enddo
   if (.not. all_edges) then
      deallocate(cut%edge_z) ; return
   endif
   allocate(cut%edge_x(nx, ny + 1, nz + 1), cut%edge_y(nx + 1, ny, nz + 1))
   do k = 1, nz + 1 ! x lines at (y, z)
      do j = 1, ny + 1
         call line_fractions(surface, 1_I4P, y(j), z(k), x, fr, nd)
         cut%edge_x(:, j, k) = fr
      enddo
   enddo
   do k = 1, nz + 1 ! y lines at (z, x)
      do i = 1, nx + 1
         call line_fractions(surface, 2_I4P, z(k), x(i), y, fr, nd)
         cut%edge_y(i, :, k) = fr
      enddo
   enddo
   endsubroutine edge_fractions

   subroutine line_fractions(surface, axis, a, b, e, frac, inside)
   !< Fluid fraction of the segments [e(s), e(s+1)] of the grid line along `axis` at (a, b) and the inside flags of
   !< its nodes e(s) (issue #70). The nodes are the perturbed grid points of #69: a crossing lies in segment s when it
   !< is after the perturbed e(s) and before the perturbed e(s+1) (`line_crossing_before`, exact); the first node's
   !< flag is the parity of the crossings before it, the next ones flip at each crossing. Lengths use the rounded
   !< crossing positions, clamped to the segment.
   type(surface_stl_object), intent(in)  :: surface   !< Surface.
   integer(I4P),             intent(in)  :: axis      !< Line axis.
   real(R8P),                intent(in)  :: a, b      !< Line position (cyclic order, see `fossil_line_query`).
   real(R8P),                intent(in)  :: e(:)      !< Nodes along the line (strictly increasing).
   real(R8P), allocatable,   intent(out) :: frac(:)   !< Fluid fractions [size(e) - 1].
   logical,   allocatable,   intent(out) :: inside(:) !< Node inside flags [size(e)].
   real(R8P),    allocatable :: t(:)                  !< Crossings.
   integer(I4P), allocatable :: fc(:), seg(:)         !< Facets, segment of each crossing (0: before e(1), ns+1: after).
   logical,      allocatable :: en(:)                 !< Entering flags.
   integer(I4P)              :: n, ns, c, sg, m       !< Work.
   real(R8P)                 :: tri(3, 3), pos, tc, fl !< Work.
   logical                   :: st                    !< Inside state.

   ns = size(e, kind=I4P) - 1_I4P
   allocate(frac(ns), inside(ns + 1))
   call surface%line_crossings(axis, a, b, -huge(1._R8P), huge(1._R8P), t, fc, en, n)
   allocate(seg(n))
   m = 0_I4P
   do c = 1, n
      tri = reshape(surface%tri_vertices(:, fc(c)), [3, 3])
      sg = max(0_I4P, min(ns + 1_I4P, col_index(t(c), e))) ! from the rounded position, then exactly
      do while (sg >= 1_I4P)
         if (.not. line_crossing_before(tri, axis, a, b, e(sg), en(c))) exit
         sg = sg - 1_I4P
      enddo
      do while (sg <= ns)
         if (line_crossing_before(tri, axis, a, b, e(sg + 1), en(c))) exit
         sg = sg + 1_I4P
      enddo
      seg(c) = sg
      if (sg == 0_I4P) m = m + 1_I4P
   enddo
   st = mod(m, 2_I4P) == 1_I4P
   inside(1) = st
   c = 1_I4P
   do sg = 1, ns
      pos = e(sg) ; fl = 0._R8P
      do c = 1, n
         if (seg(c) /= sg) cycle
         tc = min(max(t(c), e(sg)), e(sg + 1))
         if (.not. st) fl = fl + (tc - pos)
         pos = tc ; st = .not. st
      enddo
      if (.not. st) fl = fl + (e(sg + 1) - pos)
      frac(sg) = fl / (e(sg + 1) - e(sg))
      inside(sg + 1) = st
   enddo
   endsubroutine line_fractions

   pure integer(I4P) function col_index(w, e) result(k)
   !< Segment of w: 0 if w <= e(1), s if e(s) < w <= e(s+1), size(e) beyond.
   real(R8P), intent(in) :: w    !< Coordinate.
   real(R8P), intent(in) :: e(:) !< Nodes.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   if (w <= e(1)) then
      k = 0_I4P ; return
   endif
   if (w > e(size(e))) then
      k = size(e, kind=I4P) ; return
   endif
   lo = 1_I4P ; hi = size(e, kind=I4P)
   do while (hi - lo > 1_I4P)
      mid = (lo + hi) / 2_I4P
      if (w > e(mid)) then
         lo = mid
      else
         hi = mid
      endif
   enddo
   k = lo
   endfunction col_index

   pure integer(I4P) function edge_index(k, sides)
   !< Cell edge (1-12) along axis k at the sides (0 low, 1 high) of the two other axes.
   !$acc routine seq
   integer(I4P), intent(in) :: k        !< Edge axis.
   integer(I4P), intent(in) :: sides(3) !< Sides per axis (sides(k) ignored).

   edge_index = (k - 1_I4P) * 4_I4P + sides(mod(k, 3_I4P) + 1_I4P) + 2_I4P * sides(mod(k + 1_I4P, 3_I4P) + 1_I4P) + 1_I4P
   endfunction edge_index

   pure integer(I4P) function corner_of(e, side)
   !< Corner (1-8, sides x + 2 y + 4 z + 1) at the low (side 0) or high (1) end of edge e.
   !$acc routine seq
   integer(I4P), intent(in) :: e, side !< Edge, end.
   integer(I4P)             :: k, r, sides(3) !< Work.

   k = (e - 1_I4P) / 4_I4P + 1_I4P ; r = mod(e - 1_I4P, 4_I4P)
   sides(k) = side ; sides(mod(k, 3_I4P) + 1_I4P) = mod(r, 2_I4P) ; sides(mod(k + 1_I4P, 3_I4P) + 1_I4P) = r / 2_I4P
   corner_of = sides(1) + 2_I4P * sides(2) + 4_I4P * sides(3) + 1_I4P
   endfunction corner_of

   pure integer(I4P) function edge_of(fl)
   !< Cell edge of a vertex with on-plane flags fl (bit 2 (a - 1) + side): exactly two axes flagged,
   !< else 0.
   !$acc routine seq
   integer(I4P), intent(in) :: fl !< Flags.
   integer(I4P)             :: a, na, k, sides(3) !< Work.

   na = 0_I4P ; k = 0_I4P ; sides = 0_I4P
   do a = 1, 3
      if (btest(fl, 2 * (a - 1)) .or. btest(fl, 2 * (a - 1) + 1)) then
         na = na + 1_I4P
         if (btest(fl, 2 * (a - 1) + 1)) sides(a) = 1_I4P
      else
         k = a
      endif
   enddo
   edge_of = 0_I4P
   if (na == 2_I4P) edge_of = edge_index(k, sides)
   endfunction edge_of

   pure logical function lexless(a, b)
   !< Lexicographic order of two points.
   real(R8P), intent(in) :: a(3), b(3) !< Points.

   lexless = a(1) < b(1) .or. (a(1) == b(1) .and. (a(2) < b(2) .or. (a(2) == b(2) .and. a(3) < b(3))))
   endfunction lexless

   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   pure subroutine clip_flags(p, fl, n, ax, level, keep_above, bit)
   !< `clip_poly` (same kept points and crossings) carrying on-plane flags: a crossing gets the bit of
   !< the plane plus the flags its two end points share; kept points keep theirs.
   !$acc routine seq
   real(R8P),    intent(inout) :: p(3, 16)   !< Polygon.
   integer(I4P), intent(inout) :: fl(16)     !< Vertex flags.
   integer(I4P), intent(inout) :: n          !< Vertices.
   integer(I4P), intent(in)    :: ax         !< Axis.
   real(R8P),    intent(in)    :: level      !< Plane.
   logical,      intent(in)    :: keep_above !< Side kept.
   integer(I4P), intent(in)    :: bit        !< Flag bit of the plane.
   real(R8P)                   :: o(3, 16)   !< Output.
   integer(I4P)                :: of(16)     !< Output flags.
   integer(I4P)                :: k, kk, no  !< Work.
   logical                     :: ik, ikk    !< Inside flags.
   real(R8P)                   :: lo(3), hi(3) !< Ordered end points.

   if (n == 0) return
   no = 0
   do k = 1, n
      kk = mod(k, n) + 1
      if (keep_above) then ! kept side
         ik = p(ax, k) > level ; ikk = p(ax, kk) > level
      else
         ik = p(ax, k) <= level ; ikk = p(ax, kk) <= level
      endif
      if (ik) then
         no = no + 1 ; o(:, no) = p(:, k) ; of(no) = fl(k)
      endif
      if (ik .neqv. ikk) then
         if (p(1, k) < p(1, kk) .or. (p(1, k) == p(1, kk) .and. (p(2, k) < p(2, kk) .or. &
             (p(2, k) == p(2, kk) .and. p(3, k) <= p(3, kk))))) then
            lo = p(:, k) ; hi = p(:, kk)
         else
            lo = p(:, kk) ; hi = p(:, k)
         endif
         no = no + 1
         o(:, no) = lo + (level - lo(ax)) / (hi(ax) - lo(ax)) * (hi - lo)
         o(ax, no) = level
         of(no) = ibset(iand(fl(k), fl(kk)), bit)
      endif
   enddo
   n = no ; p(:, 1:n) = o(:, 1:n) ; fl(1:n) = of(1:n)
   endsubroutine clip_flags

   pure subroutine tri_terms(p1, p2, p3, a)
   !< Divergence-theorem terms of the triangle p1 p2 p3 (local coordinates, outward orientation):
   !< a(1) += 1/3 p1 . A (volume), a(2:4) += 1/2 A_k int x_k^2 dA / area (moments), a(5) += area,
   !< a(6:8) += area * centroid, a(9:11) += A (vector area).
   !$acc routine seq
   real(R8P), intent(in)    :: p1(3), p2(3), p3(3) !< Vertices.
   real(R8P), intent(inout) :: a(11)               !< Accumulators.
   real(R8P)                :: va(3), ar, q(3)     !< Vector area, area, quadratic sums.

   va = 0.5_R8P * [(p2(2) - p1(2)) * (p3(3) - p1(3)) - (p2(3) - p1(3)) * (p3(2) - p1(2)), &
                   (p2(3) - p1(3)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(3) - p1(3)), &
                   (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))]
   ar = norm2_3(va) ! gfortran's norm2, bit for bit (issue #100)
   q = (p1 * p1 + p2 * p2 + p3 * p3 + p1 * p2 + p2 * p3 + p3 * p1) / 6._R8P ! int x_k^2 dA / area
   a(1) = a(1) + dot_product(p1, va) / 3._R8P
   a(2:4) = a(2:4) + 0.5_R8P * va * q
   a(5) = a(5) + ar
   a(6:8) = a(6:8) + ar * (p1 + p2 + p3) / 3._R8P
   a(9:11) = a(9:11) + va
   endsubroutine tri_terms

   pure subroutine clip_poly(p, n, ax, level, keep_above)
   !< Clip the polygon p(:, 1:n) by the plane x(ax) = level: keep x > level (keep_above) or x <= level.
   !< Crossing points are canonical (from the edge end points in lexicographic order, the axis
   !< coordinate set to level), so neighbouring cells get the same points.
   !$acc routine seq
   real(R8P),    intent(inout) :: p(3, 16)   !< Polygon.
   integer(I4P), intent(inout) :: n          !< Vertices.
   integer(I4P), intent(in)    :: ax         !< Axis.
   real(R8P),    intent(in)    :: level      !< Plane.
   logical,      intent(in)    :: keep_above !< Side kept.
   real(R8P)                   :: o(3, 16)   !< Output.
   integer(I4P)                :: k, kk, no  !< Work.
   logical                     :: ik, ikk    !< Inside flags.
   real(R8P)                   :: lo(3), hi(3) !< Ordered end points.

   if (n == 0) return
   no = 0
   do k = 1, n
      kk = mod(k, n) + 1
      if (keep_above) then ! kept side
         ik = p(ax, k) > level ; ikk = p(ax, kk) > level
      else
         ik = p(ax, k) <= level ; ikk = p(ax, kk) <= level
      endif
      if (ik) then
         no = no + 1 ; o(:, no) = p(:, k)
      endif
      if (ik .neqv. ikk) then ! canonical crossing of the segment with the plane
         if (p(1, k) < p(1, kk) .or. (p(1, k) == p(1, kk) .and. (p(2, k) < p(2, kk) .or. &
             (p(2, k) == p(2, kk) .and. p(3, k) <= p(3, kk))))) then
            lo = p(:, k) ; hi = p(:, kk)
         else
            lo = p(:, kk) ; hi = p(:, k)
         endif
         no = no + 1
         o(:, no) = lo + (level - lo(ax)) / (hi(ax) - lo(ax)) * (hi - lo)
         o(ax, no) = level
      endif
   enddo
   n = no ; p(:, 1:n) = o(:, 1:n)
   endsubroutine clip_poly

   pure integer(I4P) function col_of(w, e)
   !< Cell column of w: 0 if w <= e(1), k if e(k) < w <= e(k+1), n+1 beyond.
   real(R8P), intent(in) :: w    !< Coordinate.
   real(R8P), intent(in) :: e(:) !< Edges.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   if (w <= e(1)) then
      col_of = 0 ; return
   endif
   if (w > e(size(e))) then
      col_of = size(e, kind=I4P) ; return
   endif
   lo = 1 ; hi = size(e, kind=I4P)
   do while (hi - lo > 1)
      mid = (lo + hi) / 2
      if (w > e(mid)) then
         lo = mid
      else
         hi = mid
      endif
   enddo
   col_of = lo
   endfunction col_of

   subroutine sections_host(surface, x, y, z, cut, wlinks, mo, qd, nm2, px, py, pz)
   !< The plane sections of `cut_cells_3d` on the host: every grid plane sliced (`slice_surface` with `on_plane_low`)
   !< and cut by `curve%cut_cells`, its cells copied into the face arrays of `cut` (allocated by the caller); with
   !< links or quadrature the sections are kept in px, py, pz (allocated by the caller).
   type(surface_stl_object),          intent(in)    :: surface    !< Surface.
   real(R8P),                         intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),              intent(inout) :: cut        !< Cut cells (face arrays allocated).
   logical,                           intent(in)    :: wlinks     !< Links asked: fluid edge intervals, sections kept.
   integer(I4P),                      intent(in)    :: mo, nm2    !< Moment order, 2-D moments per face.
   logical,                           intent(in)    :: qd         !< Quadrature asked: boundaries, sections kept.
   type(cut_cells_2d_t), allocatable, intent(inout) :: px(:), py(:), pz(:) !< Plane sections kept.
   type(curve_object)   :: sec    !< Section.
   type(cut_cells_2d_t) :: c2     !< Its cut cells.
   integer(I4P)         :: nx, ny, nz, i, j, k, st !< Work.

   nx = cut%nx ; ny = cut%ny ; nz = cut%nz
   ! planes x = x_i: (u, v) = (y, z)
   do i = 1, nx + 1
      call slice_surface(surface=surface, normal_axis=1_I4P, level=x(i), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=y, y=z, cut=c2, status=st, intervals=wlinks, moments_order=mo, boundaries=qd)
      if (wlinks .or. qd) px(i) = c2
      do k = 1, nz
         do j = 1, ny
            cut%face_x(i, j, k) = c2%alpha(j, k) ; cut%face_x_pieces(i, j, k) = c2%pieces(j, k)
            cut%face_x_closed(i, j, k) = c2%closed(j, k)
            cut%face_x_centroid(:, i, j, k) = [x(i), c2%centroid(1, j, k), c2%centroid(2, j, k)]
            if (nm2 > 0_I4P) cut%face_x_moments(:, i, j, k) = section_moments(c2, j, k, nm2)
         enddo
      enddo
   enddo
   ! planes y = y_j: (u, v) = (z, x)
   do j = 1, ny + 1
      call slice_surface(surface=surface, normal_axis=2_I4P, level=y(j), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=z, y=x, cut=c2, status=st, intervals=wlinks, moments_order=mo, boundaries=qd)
      if (wlinks .or. qd) py(j) = c2
      do k = 1, nz
         do i = 1, nx
            cut%face_y(i, j, k) = c2%alpha(k, i) ; cut%face_y_pieces(i, j, k) = c2%pieces(k, i)
            cut%face_y_closed(i, j, k) = c2%closed(k, i)
            cut%face_y_centroid(:, i, j, k) = [c2%centroid(2, k, i), y(j), c2%centroid(1, k, i)]
            if (nm2 > 0_I4P) cut%face_y_moments(:, i, j, k) = section_moments(c2, k, i, nm2)
         enddo
      enddo
   enddo
   ! planes z = z_k: (u, v) = (x, y)
   do k = 1, nz + 1
      call slice_surface(surface=surface, normal_axis=3_I4P, level=z(k), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=x, y=y, cut=c2, status=st, intervals=wlinks, moments_order=mo, boundaries=qd)
      if (wlinks .or. qd) pz(k) = c2
      do j = 1, ny
         do i = 1, nx
            cut%face_z(i, j, k) = c2%alpha(i, j) ; cut%face_z_pieces(i, j, k) = c2%pieces(i, j)
            cut%face_z_closed(i, j, k) = c2%closed(i, j)
            cut%face_z_centroid(:, i, j, k) = [c2%centroid(1, i, j), c2%centroid(2, i, j), z(k)]
            if (nm2 > 0_I4P) cut%face_z_moments(:, i, j, k) = section_moments(c2, i, j, nm2)
         enddo
      enddo
   enddo
   endsubroutine sections_host

   subroutine cut_cells_3d_sections(surface, x, y, z, cut, device, status)
   !< The faces of `cut_cells_3d` only (issue #100 step 3b, for the benches): nx, ny, nz and the face fractions,
   !< centroids, pieces and closed flags of every grid face, moments order 0, by the host sections (`sections_host`,
   !< the path of `cut_cells_3d`) or with `device = .true.` by the batched device sections (`sections_device`): the
   !< same bits. `status = STATUS_INVALID_INPUT` on bad edges.
   type(surface_stl_object), intent(in)            :: surface    !< Surface.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(cut_cells_3d_t),     intent(out)           :: cut        !< Faces (the other outputs unallocated).
   logical,                  intent(in),  optional :: device     !< The sections on the device (default .false.).
   integer(I4P),             intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   type(cut_cells_2d_t), allocatable :: px(:), py(:), pz(:) !< Not kept.
   integer(I4P)                      :: nx, ny, nz            !< Cells.
   logical                           :: dev                   !< On the device.

   if (present(status)) status = STATUS_OK
   nx = size(x, kind=I4P) - 1_I4P ; ny = size(y, kind=I4P) - 1_I4P ; nz = size(z, kind=I4P) - 1_I4P
   if (nx < 1 .or. ny < 1 .or. nz < 1) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   if (any(x(2:) <= x(:nx)) .or. any(y(2:) <= y(:ny)) .or. any(z(2:) <= z(:nz))) then
      if (present(status)) status = STATUS_INVALID_INPUT
      return
   endif
   dev = .false. ; if (present(device)) dev = device
   cut%nx = nx ; cut%ny = ny ; cut%nz = nz
   if (dev) then
      call sections_device(surface, x, y, z, cut)
   else
      allocate(cut%face_x(nx + 1, ny, nz), cut%face_y(nx, ny + 1, nz), cut%face_z(nx, ny, nz + 1))
      allocate(cut%face_x_centroid(3, nx + 1, ny, nz), cut%face_y_centroid(3, nx, ny + 1, nz), &
               cut%face_z_centroid(3, nx, ny, nz + 1))
      allocate(cut%face_x_pieces(nx + 1, ny, nz), cut%face_y_pieces(nx, ny + 1, nz), cut%face_z_pieces(nx, ny, nz + 1))
      allocate(cut%face_x_closed(nx + 1, ny, nz), cut%face_y_closed(nx, ny + 1, nz), cut%face_z_closed(nx, ny, nz + 1))
      call sections_host(surface, x, y, z, cut, .false., 0_I4P, .false., 0_I4P, px, py, pz)
   endif
   endsubroutine cut_cells_3d_sections

   subroutine sections_device(surface, x, y, z, cut)
   !< The plane sections of `cut_cells_3d` on the device (issue #100 step 3b), moments order 0, no intervals, no
   !< boundaries: per axis, the N + 1 planes are sliced on the host (`slice_surface` with `on_plane_low`, as
   !< `sections_host`) and cut as the blocks of one batched launch set (`cut_cells_device_many`, resident, cell arrays
   !< only); a device kernel (`scatter_sections`) writes their cells into the face arrays, which alone come back. The
   !< face arrays of `cut` are replaced (built in locals, moved in). Every value is a copy of the 2-D cut cell's, which
   !< equals `curve%cut_cells` bit for bit: the faces of `sections_host`.
   type(surface_stl_object), intent(in)    :: surface    !< Surface.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(inout) :: cut        !< Cut cells (nx, ny, nz set).

   call sections_axis_device(surface, 1_I4P, x, y, z, cut%nx + 1_I4P, cut%ny, cut%nz, cut%face_x, cut%face_x_centroid, &
                             cut%face_x_pieces, cut%face_x_closed)
   call sections_axis_device(surface, 2_I4P, y, z, x, cut%nx, cut%ny + 1_I4P, cut%nz, cut%face_y, cut%face_y_centroid, &
                             cut%face_y_pieces, cut%face_y_closed)
   call sections_axis_device(surface, 3_I4P, z, x, y, cut%nx, cut%ny, cut%nz + 1_I4P, cut%face_z, cut%face_z_centroid, &
                             cut%face_z_pieces, cut%face_z_closed)
   endsubroutine sections_device

   subroutine sections_axis_device(surface, axis, lv, ga, gb, n1, n2, n3, face, fcen, fpcs, fcls)
   !< The faces of the planes x(axis) = lv(b) (in-plane grid ga, gb in the cyclic order of `slice_surface`) on the
   !< device: the sections as the blocks of `cut_cells_device_many` (resident), scattered by `scatter_sections` into
   !< the face arrays [n1, n2, n3] (the 3-D shape of the axis), copied back.
   type(surface_stl_object),  intent(in)  :: surface      !< Surface.
   integer(I4P),              intent(in)  :: axis         !< Plane normal axis.
   real(R8P),                 intent(in)  :: lv(:)        !< Plane levels.
   real(R8P),                 intent(in)  :: ga(:), gb(:) !< In-plane cell edges (u, v).
   integer(I4P),              intent(in)  :: n1, n2, n3   !< Face array shape.
   real(R8P),    allocatable, intent(out) :: face(:,:,:)  !< Fluid fractions.
   real(R8P),    allocatable, intent(out) :: fcen(:,:,:,:) !< Fluid centroids [3, n1, n2, n3].
   integer(I4P), allocatable, intent(out) :: fpcs(:,:,:)  !< Fluid pieces.
   logical,      allocatable, intent(out) :: fcls(:,:,:)  !< Closed flags.
   type(curve_object),   allocatable :: secs(:)           !< Sections.
   type(cut_cells_2d_blocks_t)       :: cb                !< Their cut cells.
   real(R8P),    allocatable :: fa(:,:,:), fc(:,:,:,:)    !< Face arrays (built here).
   integer(I4P), allocatable :: fp(:,:,:)                 !< Face pieces (built here).
   logical,      allocatable :: fl(:,:,:)                 !< Face closed flags (built here).
   integer(I4P)              :: nb, b, st                 !< Work.

   nb = size(lv, kind=I4P)
   allocate(secs(nb))
   do b = 1_I4P, nb
      call slice_surface(surface=surface, normal_axis=axis, level=lv(b), curve=secs(b), on_plane_low=.true.)
   enddo
   call cut_cells_device_many(secs, ga, gb, cb, st, intervals=.false., pieces=.false., resident=.true.)
   deallocate(secs)
   allocate(fa(n1, n2, n3), fc(3, n1, n2, n3), fp(n1, n2, n3), fl(n1, n2, n3))
   !$acc data copyout(fa, fc, fp, fl) copyin(lv)
   call scatter_sections(axis, cb%nx, cb%ny, nb, lv, cb%alpha, cb%centroid, cb%pieces, cb%closed, n1, n2, n3, fa, fc, &
                         fp, fl)
   !$acc end data
   call cb%exit_device
   call move_alloc(fa, face) ; call move_alloc(fc, fcen) ; call move_alloc(fp, fpcs) ; call move_alloc(fl, fcls)
   endsubroutine sections_axis_device

   subroutine scatter_sections(axis, nu, nv, nb, lv, alpha, cen, pcs, cls, n1, n2, n3, fa, fc, fp, fl)
   !< Device kernel (issue #100 step 3b): cell (iu, iv) of the section of plane b into the face arrays of the axis, the
   !< copy of `sections_host` (x planes: face (b, iu, iv), centroid (lv, c1, c2); y planes: (iv, b, iu), (c2, lv, c1);
   !< z planes: (iu, iv, b), (c1, c2, lv)). Every array must be present on the device.
   integer(I4P), intent(in)  :: axis, nu, nv, nb       !< Axis, section cells, planes.
   real(R8P),    intent(in)  :: lv(nb)                 !< Plane levels.
   real(R8P),    intent(in)  :: alpha(nu, nv, nb)      !< Section fluid fractions.
   real(R8P),    intent(in)  :: cen(2, nu, nv, nb)     !< Section fluid centroids.
   integer(I4P), intent(in)  :: pcs(nu, nv, nb)        !< Section pieces.
   logical,      intent(in)  :: cls(nu, nv, nb)        !< Section closed flags.
   integer(I4P), intent(in)  :: n1, n2, n3             !< Face array shape.
   real(R8P),    intent(out) :: fa(n1, n2, n3)         !< Face fractions.
   real(R8P),    intent(out) :: fc(3, n1, n2, n3)      !< Face centroids.
   integer(I4P), intent(out) :: fp(n1, n2, n3)         !< Face pieces.
   logical,      intent(out) :: fl(n1, n2, n3)         !< Face closed flags.
   integer(I4P)              :: b, iu, iv              !< Counters.

   select case (axis)
   case (1)
      !$acc parallel loop gang vector collapse(3) default(present) present(lv, alpha, cen, pcs, cls, fa, fc, fp, fl)
      do iv = 1_I4P, nv
         do iu = 1_I4P, nu
            do b = 1_I4P, nb
               fa(b, iu, iv) = alpha(iu, iv, b) ; fp(b, iu, iv) = pcs(iu, iv, b) ; fl(b, iu, iv) = cls(iu, iv, b)
               fc(1, b, iu, iv) = lv(b) ; fc(2, b, iu, iv) = cen(1, iu, iv, b) ; fc(3, b, iu, iv) = cen(2, iu, iv, b)
            enddo
         enddo
      enddo
   case (2)
      !$acc parallel loop gang vector collapse(3) default(present) present(lv, alpha, cen, pcs, cls, fa, fc, fp, fl)
      do iu = 1_I4P, nu
         do b = 1_I4P, nb
            do iv = 1_I4P, nv
               fa(iv, b, iu) = alpha(iu, iv, b) ; fp(iv, b, iu) = pcs(iu, iv, b) ; fl(iv, b, iu) = cls(iu, iv, b)
               fc(1, iv, b, iu) = cen(2, iu, iv, b) ; fc(2, iv, b, iu) = lv(b) ; fc(3, iv, b, iu) = cen(1, iu, iv, b)
            enddo
         enddo
      enddo
   case default
      !$acc parallel loop gang vector collapse(3) default(present) present(lv, alpha, cen, pcs, cls, fa, fc, fp, fl)
      do b = 1_I4P, nb
         do iv = 1_I4P, nv
            do iu = 1_I4P, nu
               fa(iu, iv, b) = alpha(iu, iv, b) ; fp(iu, iv, b) = pcs(iu, iv, b) ; fl(iu, iv, b) = cls(iu, iv, b)
               fc(1, iu, iv, b) = cen(1, iu, iv, b) ; fc(2, iu, iv, b) = cen(2, iu, iv, b) ; fc(3, iu, iv, b) = lv(b)
            enddo
         enddo
      enddo
   endselect
   endsubroutine scatter_sections
endmodule fossil_cut_cells_3d
