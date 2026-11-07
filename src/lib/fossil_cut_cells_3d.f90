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
!< perturbed nodes (`column_parity`) and the piercings of an edge its exact crossings under the global
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
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_curve_slice,        only : slice_surface
use fossil_facet_object,       only : facet_object
use fossil_line_query,         only : line_crossing_before, tri_line_hit_oac
use fossil_predicates,         only : orient2d_xy
use fossil_surface_stl_object, only : surface_stl_object, STATUS_OK, STATUS_INVALID_INPUT
use penf,                      only : I4P, R8P

implicit none
private
public :: cut_cells_3d_t, cut_cells_3d

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
endtype cut_cells_3d_t

type :: pierce_t
   !< A wall polygon vertex on a cell edge.
   integer(I4P) :: edge = 0_I4P     !< Cell edge (1-12).
   real(R8P)    :: t = 0._R8P       !< Position along the edge.
   logical      :: enter = .false.  !< The edge, going in +axis, enters the solid there.
   integer(I4P) :: poly = 0_I4P     !< Polygon.
endtype pierce_t

contains
   subroutine cut_cells_3d(surface, x, y, z, cut, status)
   !< Cut cells of `surface` on the block with cell edges x(nx+1), y(ny+1), z(nz+1) (issue #40).
   !<
   !< Step 4: the face fluid fractions, centroids and pieces of every grid face (see the module
   !< notes). `status = STATUS_INVALID_INPUT` on bad edges.
   type(surface_stl_object), intent(in)            :: surface    !< Surface.
   real(R8P),                intent(in)            :: x(:), y(:), z(:) !< Cell edges (strictly increasing).
   type(cut_cells_3d_t),     intent(out)           :: cut        !< Cut cells.
   integer(I4P),             intent(out), optional :: status     !< STATUS_OK or STATUS_INVALID_INPUT.
   integer(I4P)         :: nx, ny, nz, i, j, k, st, f !< Work.
   integer(I4P)         :: il, ih, jl, jh, kl, kh, m, np !< Work.
   real(R8P)            :: tri(3, 3), poly(3, 16), c0(3), h(3), sq(3), vc, vs, ms(3) !< Work.
   real(R8P), allocatable :: acc(:,:,:,:)  !< Per cell: solid volume, moments (3), wet area, moments (3), normal (3).
   type(facet_object), pointer :: fp      !< Facet.
   type(curve_object)   :: sec    !< Section.
   type(cut_cells_2d_t) :: c2     !< Its cut cells.

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
   ! planes x = x_i: (u, v) = (y, z)
   do i = 1, nx + 1
      call slice_surface(surface=surface, normal_axis=1_I4P, level=x(i), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=y, y=z, cut=c2, status=st)
      do k = 1, nz
         do j = 1, ny
            cut%face_x(i, j, k) = c2%alpha(j, k) ; cut%face_x_pieces(i, j, k) = c2%pieces(j, k)
            cut%face_x_closed(i, j, k) = c2%closed(j, k)
            cut%face_x_centroid(:, i, j, k) = [x(i), c2%centroid(1, j, k), c2%centroid(2, j, k)]
         enddo
      enddo
   enddo
   ! planes y = y_j: (u, v) = (z, x)
   do j = 1, ny + 1
      call slice_surface(surface=surface, normal_axis=2_I4P, level=y(j), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=z, y=x, cut=c2, status=st)
      do k = 1, nz
         do i = 1, nx
            cut%face_y(i, j, k) = c2%alpha(k, i) ; cut%face_y_pieces(i, j, k) = c2%pieces(k, i)
            cut%face_y_closed(i, j, k) = c2%closed(k, i)
            cut%face_y_centroid(:, i, j, k) = [c2%centroid(2, k, i), y(j), c2%centroid(1, k, i)]
         enddo
      enddo
   enddo
   ! planes z = z_k: (u, v) = (x, y)
   do k = 1, nz + 1
      call slice_surface(surface=surface, normal_axis=3_I4P, level=z(k), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=x, y=y, cut=c2, status=st)
      do j = 1, ny
         do i = 1, nx
            cut%face_z(i, j, k) = c2%alpha(i, j) ; cut%face_z_pieces(i, j, k) = c2%pieces(i, j)
            cut%face_z_closed(i, j, k) = c2%closed(i, j)
            cut%face_z_centroid(:, i, j, k) = [c2%centroid(1, i, j), c2%centroid(2, i, j), z(k)]
         enddo
      enddo
   enddo
   ! cells: clipped triangles
   allocate(acc(11, nx, ny, nz)) ; acc = 0._R8P
   do f = 1, surface%get_facets_number()
      fp => surface%facet_at(f)
      do m = 1, 3
         tri(:, m) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
      il = max(1_I4P, col_of(minval(tri(1, :)), x)) ; ih = min(nx, col_of(maxval(tri(1, :)), x))
      jl = max(1_I4P, col_of(minval(tri(2, :)), y)) ; jh = min(ny, col_of(maxval(tri(2, :)), y))
      kl = max(1_I4P, col_of(minval(tri(3, :)), z)) ; kh = min(nz, col_of(maxval(tri(3, :)), z))
      do k = kl, kh
         do j = jl, jh
            do i = il, ih
               poly(:, 1:3) = tri ; np = 3
               call clip_plane(poly, np, 1, x(i), .true.) ; call clip_plane(poly, np, 1, x(i + 1), .false.)
               call clip_plane(poly, np, 2, y(j), .true.) ; call clip_plane(poly, np, 2, y(j + 1), .false.)
               call clip_plane(poly, np, 3, z(k), .true.) ; call clip_plane(poly, np, 3, z(k + 1), .false.)
               if (np < 3) cycle
               c0 = [x(i), y(j), z(k)]
               do m = 2, np - 1 ! fan, local coordinates
                  call tri_terms(poly(:, 1) - c0, poly(:, m) - c0, poly(:, m + 1) - c0, acc(:, i, j, k))
               enddo
            enddo
         enddo
      enddo
   enddo
   ! cells: faces and totals
   allocate(cut%alpha(nx, ny, nz), cut%centroid(3, nx, ny, nz), cut%wet_area(nx, ny, nz), cut%wet_centroid(3, nx, ny, nz), &
            cut%wet_normal(3, nx, ny, nz), cut%closed(nx, ny, nz))
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c0 = [x(i), y(j), z(k)] ; h = [x(i + 1) - x(i), y(j + 1) - y(j), z(k + 1) - z(k)]
            vc = h(1) * h(2) * h(3)
            cut%closed(i, j, k) = cut%face_x_closed(i, j, k) .and. cut%face_x_closed(i + 1, j, k) .and. &
                                  cut%face_y_closed(i, j, k) .and. cut%face_y_closed(i, j + 1, k) .and. &
                                  cut%face_z_closed(i, j, k) .and. cut%face_z_closed(i, j, k + 1)
            ! +x, +y, +z faces (local coordinate h); the -x, -y, -z faces are at local 0 and add nothing
            sq = [(1._R8P - cut%face_x(i + 1, j, k)) * h(2) * h(3), (1._R8P - cut%face_y(i, j + 1, k)) * h(3) * h(1), &
                  (1._R8P - cut%face_z(i, j, k + 1)) * h(1) * h(2)]
            vs = acc(1, i, j, k) + sum(h * sq) / 3._R8P
            ms = acc(2:4, i, j, k) + 0.5_R8P * h * h * sq
            cut%wet_area(i, j, k) = acc(5, i, j, k)
            cut%wet_normal(:, i, j, k) = acc(9:11, i, j, k)
            cut%wet_centroid(:, i, j, k) = c0 + 0.5_R8P * h
            if (acc(5, i, j, k) > 0._R8P) cut%wet_centroid(:, i, j, k) = c0 + acc(6:8, i, j, k) / acc(5, i, j, k)
            if (.not. cut%closed(i, j, k)) then
               cut%alpha(i, j, k) = -1._R8P ; cut%centroid(:, i, j, k) = c0 + 0.5_R8P * h
               cycle
            endif
            cut%alpha(i, j, k) = (vc - vs) / vc
            cut%centroid(:, i, j, k) = c0 + 0.5_R8P * h
            if (vc - vs > 0._R8P) cut%centroid(:, i, j, k) = c0 + (vc * 0.5_R8P * h - ms) / (vc - vs)
         enddo
      enddo
   enddo
   call cut_pieces_3d(surface, x, y, z, cut)
   endsubroutine cut_cells_3d

   subroutine cut_pieces_3d(surface, x, y, z, cut)
   !< Fluid pieces of every cell (see the module notes).
   type(surface_stl_object), intent(in)    :: surface    !< Surface.
   real(R8P),                intent(in)    :: x(:), y(:), z(:) !< Cell edges.
   type(cut_cells_3d_t),     intent(inout) :: cut        !< Cut cells (alpha, closed set).
   type(facet_object), pointer :: fp                      !< Facet.
   real(R8P),    allocatable   :: tri(:,:,:), nrm(:,:)    !< Triangles [3, 3, nf], normals [3, nf].
   integer(I4P), allocatable   :: cnt(:), start(:), lst(:) !< Cell -> triangles (CSR).
   integer(I4P)                :: nx, ny, nz, nf, f, m, i, j, k, il, ih, jl, jh, kl, kh, c, np !< Work.
   logical,      allocatable   :: open_cell(:,:,:)        !< Cells a mesh boundary edge overlaps.
   real(R8P),    allocatable   :: ek(:,:)                 !< Edge keys [6, 3 nf] (end points, lexicographic order).
   integer(I4P), allocatable   :: idx(:)                  !< Sorted edges.
   real(R8P)                   :: pa(3), pb(3)            !< Edge ends.
   real(R8P),    allocatable   :: pm(:,:)                 !< Measures of the pieces of a cell [17, np].
   real(R8P),    allocatable   :: pbuf(:,:)               !< Measures of all pieces [17, capacity].
   integer(I4P), allocatable   :: pcel(:,:)               !< Cells of all pieces [3, capacity].
   integer(I4P)                :: npb                     !< Pieces so far.
   real(R8P)                   :: tot(17), h(3)           !< Cell totals, cell size.
   logical,      allocatable   :: node_in(:,:,:)          !< Exact inside flags of the (perturbed) grid nodes.
   logical,      allocatable   :: need(:,:)               !< Node columns that a cut cell needs.
   integer(I4P)                :: cs(8), q                !< Corner states, corner.

   nx = cut%nx ; ny = cut%ny ; nz = cut%nz
   nf = surface%get_facets_number()
   allocate(tri(3, 3, nf), nrm(3, nf), cnt(nx * ny * nz + 1), cut%pieces(nx, ny, nz), cut%first_piece(nx, ny, nz))
   allocate(pbuf(17, 64), pcel(3, 64)) ; npb = 0_I4P
   cnt = 0_I4P
   do f = 1, nf
      fp => surface%facet_at(f)
      do m = 1, 3
         tri(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
      enddo
      nrm(:, f) = cross3(tri(:, 2, f) - tri(:, 1, f), tri(:, 3, f) - tri(:, 1, f))
   enddo
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
   ! mesh boundary edges (one facet): the cells they overlap get no count
   allocate(open_cell(nx, ny, nz), ek(6, 3 * nf), idx(3 * nf)) ; open_cell = .false.
   do f = 1, nf
      do m = 1, 3
         pa = tri(:, m, f) ; pb = tri(:, mod(m, 3) + 1, f)
         if (lexless(pb, pa)) then
            ek(:, 3 * (f - 1) + m) = [pb, pa]
         else
            ek(:, 3 * (f - 1) + m) = [pa, pb]
         endif
      enddo
   enddo
   call sort_points(ek, idx)
   m = 1
   do while (m <= 3 * nf)
      f = m
      do while (f < 3 * nf)
         if (any(ek(:, idx(f + 1)) /= ek(:, idx(m)))) exit
         f = f + 1
      enddo
      if (f == m) then ! an edge of one facet
         pa = ek(1:3, idx(m)) ; pb = ek(4:6, idx(m))
         il = max(1_I4P, col_of(min(pa(1), pb(1)), x)) ; ih = min(nx, col_of(max(pa(1), pb(1)), x))
         jl = max(1_I4P, col_of(min(pa(2), pb(2)), y)) ; jh = min(ny, col_of(max(pa(2), pb(2)), y))
         kl = max(1_I4P, col_of(min(pa(3), pb(3)), z)) ; kh = min(nz, col_of(max(pa(3), pb(3)), z))
         if (il <= ih .and. jl <= jh .and. kl <= kh) open_cell(il:ih, jl:jh, kl:kh) = .true.
      endif
      m = f + 1
   enddo
   ! corner states: the exact inside flags of the perturbed nodes (issue #70), one z line per needed column
   allocate(need(nx + 1, ny + 1), node_in(nx + 1, ny + 1, nz + 1)) ; need = .false. ; node_in = .false.
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            if (cnt(i + nx * ((j - 1) + ny * (k - 1))) > 0_I4P) need(i:i + 1, j:j + 1) = .true.
         enddo
      enddo
   enddo
   do j = 1, ny + 1
      do i = 1, nx + 1
         if (need(i, j)) call column_parity(surface, x(i), y(j), z, node_in(i, j, :))
      enddo
   enddo
   do k = 1, nz
      do j = 1, ny
         do i = 1, nx
            c = i + nx * ((j - 1) + ny * (k - 1))
            cut%first_piece(i, j, k) = npb + 1_I4P
            if (.not. cut%closed(i, j, k) .or. open_cell(i, j, k) .or. cut%alpha(i, j, k) < -1.e-12_R8P .or. &
                cut%alpha(i, j, k) > 1._R8P + 1.e-12_R8P) then
               cut%pieces(i, j, k) = -1_I4P
            elseif (cnt(c) == 0_I4P .and. abs(cut%alpha(i, j, k) - anint(cut%alpha(i, j, k))) > 1.e-12_R8P) then
               cut%pieces(i, j, k) = -1_I4P
            elseif (cnt(c) == 0_I4P) then
               cut%pieces(i, j, k) = merge(1_I4P, 0_I4P, cut%alpha(i, j, k) > 0.5_R8P)
            else
               do q = 1, 8
                  cs(q) = merge(0_I4P, 1_I4P, node_in(i + mod(q - 1, 2), j + mod((q - 1) / 2, 2), k + (q - 1) / 4))
               enddo
               call cell_pieces(tri, nrm, lst(start(c):start(c + 1) - 1), &
                                [x(i), y(j), z(k)], [x(i + 1), y(j + 1), z(k + 1)], cs, np, pm)
               if (np > 1_I4P) then ! (issue #70) trust the split only if it closes
                  h = [x(i + 1) - x(i), y(j + 1) - y(j), z(k + 1) - z(k)]
                  if (.not. split_closes(pm(:, 1:np), h, cut%alpha(i, j, k), [cut%face_x(i, j, k), cut%face_x(i + 1, j, k), &
                      cut%face_y(i, j, k), cut%face_y(i, j + 1, k), cut%face_z(i, j, k), cut%face_z(i, j, k + 1)])) np = -1_I4P
               endif
               cut%pieces(i, j, k) = np
            endif
            if (cut%pieces(i, j, k) < 1_I4P) cycle
            ! measures: the cell totals for a single piece (exact), else the per-piece ones
            if (cut%pieces(i, j, k) == 1_I4P) then
               h = [x(i + 1) - x(i), y(j + 1) - y(j), z(k + 1) - z(k)]
               tot = [cut%alpha(i, j, k) * product(h), cut%centroid(:, i, j, k), cut%wet_area(i, j, k), &
                      cut%wet_centroid(:, i, j, k), cut%wet_normal(:, i, j, k), cut%face_x(i, j, k), cut%face_x(i + 1, j, k), &
                      cut%face_y(i, j, k), cut%face_y(i, j + 1, k), cut%face_z(i, j, k), cut%face_z(i, j, k + 1)]
               call add_pieces(reshape(tot, [17, 1]))
            else
               call add_pieces(pm)
            endif
         enddo
      enddo
   enddo
   cut%np = npb
   allocate(cut%piece_cell(3, npb), cut%piece_volume(npb), cut%piece_centroid(3, npb), cut%piece_wet_area(npb), &
            cut%piece_wet_centroid(3, npb), cut%piece_wet_normal(3, npb), cut%piece_face(6, npb))
   cut%piece_cell = pcel(:, 1:npb)
   cut%piece_volume = pbuf(1, 1:npb) ; cut%piece_centroid = pbuf(2:4, 1:npb) ; cut%piece_wet_area = pbuf(5, 1:npb)
   cut%piece_wet_centroid = pbuf(6:8, 1:npb) ; cut%piece_wet_normal = pbuf(9:11, 1:npb) ; cut%piece_face = pbuf(12:17, 1:npb)
   contains
      subroutine add_pieces(pp)
      !< Append the pieces pp [17, n] of the cell (i, j, k).
      real(R8P), intent(in)     :: pp(:,:) !< Measures.
      real(R8P),    allocatable :: tb(:,:) !< Growth.
      integer(I4P), allocatable :: tc(:,:) !< Growth.
      integer(I4P)              :: n, cap  !< Work.

      n = size(pp, 2, kind=I4P)
      if (npb + n > size(pbuf, 2)) then
         cap = max(2_I4P * size(pbuf, 2, kind=I4P), npb + n)
         allocate(tb(17, cap), tc(3, cap)) ; tb(:, 1:npb) = pbuf(:, 1:npb) ; tc(:, 1:npb) = pcel(:, 1:npb)
         call move_alloc(tb, pbuf) ; call move_alloc(tc, pcel)
      endif
      pbuf(:, npb + 1:npb + n) = pp
      pcel(1, npb + 1:npb + n) = i ; pcel(2, npb + 1:npb + n) = j ; pcel(3, npb + 1:npb + n) = k
      npb = npb + n
      endsubroutine add_pieces
   endsubroutine cut_pieces_3d

   subroutine cell_pieces(tri, nrm, lst, c0, c1, cs, npieces, pm)
   !< Fluid pieces of the cell [c0, c1] cut by the triangles lst (see the module notes); -1 if the
   !< boundary states do not close. pm: the measures of the pieces (see `piece_measures`).
   !<
   !< The corner states cs and the crossings of the 12 cell edges are exact under the global perturbation
   !< (issue #70): a crossing belongs to the edge when it lies after the perturbed low corner and before the
   !< perturbed high one (`line_crossing_before`), so the corner states and the crossings of every edge agree.
   !< Before #70 the piercings were the clipped polygon vertices on two planes, and a mesh edge or vertex on a
   !< cell edge (lattice bodies) gave one piercing per incident facet: the states did not alternate.
   real(R8P),                intent(in)  :: tri(:,:,:), nrm(:,:) !< Triangles, normals.
   integer(I4P),             intent(in)  :: lst(:)    !< Triangles of the cell.
   real(R8P),                intent(in)  :: c0(3), c1(3) !< Cell corners.
   integer(I4P),             intent(in)  :: cs(8)     !< Corner states (0 solid, 1 fluid), exact.
   integer(I4P),             intent(out) :: npieces   !< Pieces.
   real(R8P), allocatable,   intent(out) :: pm(:,:)   !< Measures [17, npieces] (allocated if npieces >= 0).
   integer(I4P), parameter     :: NCORN = 8_I4P       !< Corner atoms follow the polygons.
   real(R8P),    allocatable   :: pv(:,:,:)           !< Polygons [3, 16, nt].
   integer(I4P), allocatable   :: pfl(:,:), pn(:)     !< Vertex flags, vertex counts.
   type(pierce_t), allocatable :: pc(:)               !< Piercings.
   integer(I4P), allocatable   :: parent(:)           !< Union-find.
   integer(I4P), allocatable   :: iv_edge(:), iv_atom(:) !< Fluid intervals: edge, atom.
   real(R8P),    allocatable   :: iv_lo(:), iv_hi(:)  !< Fluid intervals: range.
   logical,      allocatable   :: touch(:)            !< Atom touches the cell boundary.
   real(R8P),    allocatable   :: vol(:)              !< Solid volume of closed wall sets.
   integer(I4P)                :: cstate(8)           !< Corner states: -1 unknown, 0 solid, 1 fluid.
   integer(I4P)                :: nt, np, q, v, w, e, a, natom, niv, npc, nc, kk !< Work.
   logical                     :: ok, changed         !< Consistency, propagation.
   real(R8P)                   :: tt, la, lb, tc      !< Crossing position, edge line, kernel parameter.
   logical                     :: hit, en             !< Crossing, entering.

   nt = size(lst, kind=I4P)
   allocate(pv(3, 16, nt), pfl(16, nt), pn(nt))
   do q = 1, nt
      pv(:, 1:3, q) = tri(:, :, lst(q)) ; pfl(1:3, q) = 0_I4P ; np = 3_I4P
      do a = 1, 3
         call clip_flags(pv(:, :, q), pfl(:, q), np, a, c0(a), .true., 2_I4P * (a - 1_I4P))
         call clip_flags(pv(:, :, q), pfl(:, q), np, a, c1(a), .false., 2_I4P * (a - 1_I4P) + 1_I4P)
      enddo
      pn(q) = np ; if (np < 3_I4P) pn(q) = 0_I4P
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
   ! piercings: the exact crossings of the cell edges (issue #70)
   allocate(pc(12 * nt)) ; npc = 0_I4P
   do e = 1, 12
      a = (e - 1_I4P) / 4_I4P + 1_I4P ; kk = mod(e - 1_I4P, 4_I4P) ! edge axis, sides
      la = merge(c1(mod(a, 3_I4P) + 1_I4P), c0(mod(a, 3_I4P) + 1_I4P), mod(kk, 2_I4P) == 1_I4P)
      lb = merge(c1(mod(a + 1_I4P, 3_I4P) + 1_I4P), c0(mod(a + 1_I4P, 3_I4P) + 1_I4P), kk / 2_I4P == 1_I4P)
      do q = 1, nt
         call tri_line_hit_oac(reshape(tri(:, :, lst(q)), [9]), a, la, lb, hit, tc, en)
         if (.not. hit) cycle
         if (line_crossing_before(tri(:, :, lst(q)), a, la, lb, c0(a), en)) cycle
         if (.not. line_crossing_before(tri(:, :, lst(q)), a, la, lb, c1(a), en)) cycle
         tt = min(max(tc, c0(a)), c1(a)) ! the clipped polygon's vertex on the edge if it has one (same rounding)
         do v = 1, pn(q)
            if (edge_of(pfl(v, q)) == e) then
               tt = pv(a, v, q) ; exit
            endif
         enddo
         npc = npc + 1_I4P
         pc(npc) = pierce_t(e, tt, en, merge(q, 0_I4P, pn(q) > 0_I4P))
      enddo
   enddo
   call sort_pierce(pc, npc)
   ! corner states from the pierced edges, then propagation along the others
   ok = .true. ; cstate = cs
   do e = 1, 12 ! check: the crossings alternate from the low corner's state to the high one's
      call walk_edge(e, .false.)
   enddo
   if (any(cstate == -1_I4P)) ok = .false.
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   ! atoms: polygons 1..nt, corners nt+1..nt+8, intervals after
   allocate(parent(nt + NCORN + npc + 12), iv_edge(npc + 12), iv_atom(npc + 12), iv_lo(npc + 12), iv_hi(npc + 12))
   parent = [(a, a = 1, size(parent))]
   niv = 0_I4P ; natom = nt + NCORN
   do e = 1, 12
      call walk_edge(e, .true.)
   enddo
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   call join_shared_vertices
   do a = 1, 3
      call faces_islands(a, 0_I4P)
      call faces_islands(a, 1_I4P)
   enddo
   if (.not. ok) then
      npieces = -1_I4P ; return
   endif
   ! count: sets touching the boundary, plus closed wall sets that bound a void
   allocate(touch(natom), vol(natom)) ; touch = .false. ; vol = 0._R8P
   do q = 1, nt
      if (pn(q) == 0_I4P) cycle
      w = find(q)
      if (any(pfl(1:pn(q), q) /= 0_I4P)) touch(w) = .true.
      do v = 2, pn(q) - 1
         vol(w) = vol(w) + dot_product(pv(:, 1, q) - c0, &
                  cross3(pv(:, v, q) - pv(:, 1, q), pv(:, v + 1, q) - pv(:, 1, q))) / 6._R8P
      enddo
   enddo
   do a = 1, NCORN
      if (cstate(a) == 1_I4P) touch(find(nt + a)) = .true.
   enddo
   do a = 1, niv
      touch(find(iv_atom(a))) = .true.
   enddo
   nc = 0_I4P
   do a = 1, natom
      if (find(a) /= a) cycle
      if (touch(a)) then
         nc = nc + 1_I4P
      elseif (a <= nt) then
         if (pn(a) > 0_I4P .and. vol(a) < 0._R8P) nc = nc + 1_I4P
      endif
   enddo
   npieces = nc
   call piece_measures
   contains
      integer(I4P) function find(i) result(r)
      !< Root of i (path halving).
      integer(I4P), intent(in) :: i !< Atom.
      integer(I4P)             :: j !< Work.

      j = i
      do while (parent(j) /= j)
         parent(j) = parent(parent(j)) ; j = parent(j)
      enddo
      r = j
      endfunction find

      subroutine unite(i, j)
      !< Join the sets of i and j.
      integer(I4P), intent(in) :: i, j !< Atoms.
      integer(I4P)             :: ri, rj !< Roots.

      ri = find(i) ; rj = find(j)
      if (ri /= rj) parent(max(ri, rj)) = min(ri, rj)
      endsubroutine unite

      subroutine propagate(ca, cb)
      !< Equal states at the two corners of an edge without piercings.
      integer(I4P), intent(in) :: ca, cb !< Corners.

      if (cstate(ca) == -1_I4P .and. cstate(cb) /= -1_I4P) then
         cstate(ca) = cstate(cb) ; changed = .true.
      elseif (cstate(cb) == -1_I4P .and. cstate(ca) /= -1_I4P) then
         cstate(cb) = cstate(ca) ; changed = .true.
      elseif (cstate(ca) /= cstate(cb)) then
         ok = .false.
      endif
      endsubroutine propagate

      subroutine set_corner(cc, st)
      !< Set (or check) a corner state.
      integer(I4P), intent(in) :: cc, st !< Corner, state.

      if (cstate(cc) == -1_I4P) then
         cstate(cc) = st
      elseif (cstate(cc) /= st) then
         ok = .false.
      endif
      endsubroutine set_corner

      subroutine walk_edge(ed, build)
      !< Walk the piercings of edge ed from its low corner's exact state. Without `build`: check that they
      !< alternate and end at the high corner's state, reordering crossings at one rounded position (a mesh
      !< edge or vertex on the cell edge: their exact order is the alternation). With `build`: the fluid
      !< interval atoms and their joins.
      integer(I4P), intent(in) :: ed    !< Edge.
      logical,      intent(in) :: build !< Build the intervals.
      integer(I4P)             :: st, ii, i1, i2, prev, ax, jt, kt !< State, piercings, previous endpoint atom.
      real(R8P)                :: tprev, lo, hi !< Previous position, edge range.
      type(pierce_t)           :: tmp   !< Swap.

      ax = (ed - 1_I4P) / 4_I4P + 1_I4P
      lo = c0(ax) ; hi = c1(ax)
      i1 = 0_I4P ; i2 = -1_I4P
      do ii = 1, npc
         if (pc(ii)%edge /= ed) cycle
         if (i1 == 0_I4P) i1 = ii
         i2 = ii
      enddo
      if (i1 == 0_I4P) then ! no piercing
         if (build .and. cstate(corner_of(ed, 0_I4P)) == 1_I4P) then
            call new_interval(ed, lo, hi)
            call unite(natom, nt + corner_of(ed, 0_I4P)) ; call unite(natom, nt + corner_of(ed, 1_I4P))
         endif
         return
      endif
      st = cstate(corner_of(ed, 0_I4P)) ! exact (issue #70); crossings at one rounded position are reordered below
      prev = nt + corner_of(ed, 0_I4P) ; tprev = lo
      do ii = i1, i2
         if (pc(ii)%enter .neqv. (st == 1_I4P)) then ! alternation broken: repair a tie, else not closed
            jt = 0_I4P
            do kt = ii + 1_I4P, i2
               if (abs(pc(kt)%t - pc(ii)%t) > 64._R8P * spacing(max(abs(lo), abs(hi)))) exit
               if (pc(kt)%enter .eqv. (st == 1_I4P)) then
                  jt = kt ; exit
               endif
            enddo
            if (jt == 0_I4P) then
               ok = .false. ; return
            endif
            tmp = pc(ii) ; pc(ii) = pc(jt) ; pc(jt) = tmp
         endif
         if (build .and. st == 1_I4P) then
            call new_interval(ed, tprev, pc(ii)%t)
            if (prev > 0_I4P) call unite(natom, prev) ! 0: a crossing without a polygon in the cell
            if (pc(ii)%poly > 0_I4P) call unite(natom, pc(ii)%poly)
         endif
         st = 1_I4P - st ; prev = pc(ii)%poly ; tprev = pc(ii)%t
      enddo
      if (.not. build) call set_corner(corner_of(ed, 1_I4P), st)
      if (build .and. st == 1_I4P) then
         call new_interval(ed, tprev, hi)
         if (prev > 0_I4P) call unite(natom, prev)
         call unite(natom, nt + corner_of(ed, 1_I4P))
      endif
      endsubroutine walk_edge

      subroutine new_interval(ed, tlo, thi)
      !< A new fluid interval atom (natom).
      integer(I4P), intent(in) :: ed       !< Edge.
      real(R8P),    intent(in) :: tlo, thi !< Range.

      natom = natom + 1_I4P ; niv = niv + 1_I4P
      iv_edge(niv) = ed ; iv_lo(niv) = tlo ; iv_hi(niv) = thi ; iv_atom(niv) = natom
      endsubroutine new_interval

      subroutine join_shared_vertices
      !< Join the polygons that share a vertex (bit-identical points).
      real(R8P),    allocatable :: key(:,:) !< Vertices [3, m].
      integer(I4P), allocatable :: own(:), idx(:) !< Owner polygon, sorted order.
      integer(I4P)              :: mm, ii, qq, vv !< Work.

      mm = sum(pn)
      allocate(key(3, mm), own(mm), idx(mm))
      mm = 0_I4P
      do qq = 1, nt
         do vv = 1, pn(qq)
            mm = mm + 1_I4P ; key(:, mm) = pv(:, vv, qq) ; own(mm) = qq
         enddo
      enddo
      call sort_points(key, idx)
      do ii = 2, mm
         if (all(key(:, idx(ii)) == key(:, idx(ii - 1)))) call unite(own(idx(ii)), own(idx(ii - 1)))
      enddo
      endsubroutine join_shared_vertices

      subroutine faces_islands(ax, sd)
      !< Join every solid island traced on the face (axis ax, side sd) to the fluid around it.
      integer(I4P), intent(in) :: ax, sd !< Face.
      integer(I4P)             :: bt, a1, a2, ns, ii, jj, qq, vv, best, comp, sides(3), ed !< Work.
      real(R8P),    allocatable :: su(:,:), sv(:,:)  !< Segment end points in face coordinates [2, ns].
      integer(I4P), allocatable :: sp(:), sc(:)      !< Segment polygon, component.
      logical,      allocatable :: stouch(:)         !< Segment component touches a cell edge.
      real(R8P)                 :: uv(2), vv0, uc, ubest, nbest !< Work.
      integer(I4P)              :: iv                 !< Vertex of the rightmost point.

      bt = 2_I4P * (ax - 1_I4P) + sd
      a1 = mod(ax, 3_I4P) + 1_I4P ; a2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P
      ns = 0_I4P
      do qq = 1, nt
         do vv = 1, pn(qq)
            if (btest(pfl(vv, qq), bt) .and. btest(pfl(mod(vv, pn(qq)) + 1, qq), bt)) ns = ns + 1_I4P
         enddo
      enddo
      if (ns == 0_I4P) return
      allocate(su(2, ns), sv(2, ns), sp(ns), sc(ns), stouch(ns))
      ns = 0_I4P
      do qq = 1, nt
         do vv = 1, pn(qq)
            jj = mod(vv, pn(qq)) + 1
            if (.not. (btest(pfl(vv, qq), bt) .and. btest(pfl(jj, qq), bt))) cycle
            ns = ns + 1_I4P
            su(:, ns) = [pv(a1, vv, qq), pv(a1, jj, qq)] ; sv(:, ns) = [pv(a2, vv, qq), pv(a2, jj, qq)]
            sp(ns) = qq
            stouch(ns) = edge_of(pfl(vv, qq)) /= 0_I4P .or. edge_of(pfl(jj, qq)) /= 0_I4P .or. &
                         popcnt(pfl(vv, qq)) > 2 .or. popcnt(pfl(jj, qq)) > 2
         enddo
      enddo
      ! segment components by shared end points (small sets: quadratic)
      sc = [(ii, ii = 1, ns)]
      do ii = 1, ns
         do jj = ii + 1, ns
            if (sc(ii) == sc(jj)) cycle
            if ((su(1, ii) == su(1, jj) .and. sv(1, ii) == sv(1, jj)) .or. (su(1, ii) == su(2, jj) .and. sv(1, ii) == sv(2, jj)) &
                .or. (su(2, ii) == su(1, jj) .and. sv(2, ii) == sv(1, jj)) .or. &
                (su(2, ii) == su(2, jj) .and. sv(2, ii) == sv(2, jj))) then
               iv = max(sc(ii), sc(jj)) ; best = min(sc(ii), sc(jj))
               where (sc == iv) sc = best
            endif
         enddo
      enddo
      do ii = 1, ns
         if (any(stouch .and. sc == sc(ii))) stouch(ii) = .true.
      enddo
      do comp = 1, ns
         if (.not. any(sc == comp)) cycle
         ii = findloc(sc, comp, 1)
         if (stouch(ii)) cycle
         ! rightmost point of the loop and the incident segment with the largest |n_u|
         uv = [-huge(1._R8P), -huge(1._R8P)]
         do jj = 1, ns
            if (sc(jj) /= comp) cycle
            do iv = 1, 2
               if (su(iv, jj) > uv(1) .or. (su(iv, jj) == uv(1) .and. sv(iv, jj) > uv(2))) uv = [su(iv, jj), sv(iv, jj)]
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
            call unite(sp(best), qq)
            cycle
         endif
         sides = 0_I4P ; sides(ax) = sd ; sides(a1) = 1_I4P
         ed = edge_index(a2, sides)
         qq = 0_I4P
         do jj = 1, niv
            if (iv_edge(jj) == ed .and. iv_lo(jj) <= vv0 .and. vv0 <= iv_hi(jj)) qq = iv_atom(jj)
         enddo
         if (qq == 0_I4P) then
            ok = .false. ; return
         endif
         call unite(sp(best), qq)
      enddo
      endsubroutine faces_islands

      subroutine piece_measures
      !< Measures of the counted sets, in the order of their roots: pm(1) fluid volume, pm(2:4) centroid,
      !< pm(5) wetted area, pm(6:8) wetted centroid, pm(9:11) outward normal integral, pm(12:17) fluid
      !< fraction of the faces -x, +x, -y, +y, -z, +z.
      !<
      !< Walls: the divergence-theorem terms of every polygon go to its set (the outward normal of the
      !< fluid is minus the triangle's). Faces: the fluid area of each face by Green's theorem in face
      !< coordinates local to c0, A = oint u dv counter-clockwise seen from +e_a, over the boundary arcs
      !< of the fluid region of the face: the fluid edge intervals (an arc on a low edge, u = 0 or v = 0,
      !< adds nothing) and the trace segments. A clipped polygon lies inside the cell, so the fluid is on
      !< the left of its trace as wound on a high face and on the right on a low face: the orientation
      !< comes from the winding, no sign test. Volume and moments need the high faces only (local
      !< coordinates). An inclusion (closed set, positive solid volume) goes to the set around it: the
      !< nearest wall a ray in +x from its rightmost vertex hits (exact inclusive test on the projected
      !< triangle), else the fluid region where the ray leaves the cell through the +x face (the ray of
      !< `faces_islands` on that face); -1 if neither is found.
      real(R8P),    allocatable :: wt(:,:), fa(:,:) !< Wall terms [11, natom], face fluid areas [6, natom].
      integer(I4P), allocatable :: tgt(:), pid(:)   !< Inclusion target, piece of a root.
      real(R8P)                 :: h(3), pl(3), ql(3), sg, vp, mp(3), hf(3) !< Work.
      integer(I4P)              :: qq, vv, jj, r, ax, sd, bt, b1, b2, kk, ee, k1, k2, s1, s2, guard, npc_ !< Work.
      logical                   :: inplane, outside !< Polygon in the face plane, its outward normal leaves the cell.

      h = c1 - c0
      allocate(wt(11, natom), fa(6, natom), tgt(natom), pid(natom))
      wt = 0._R8P ; fa = 0._R8P ; tgt = 0_I4P ; pid = 0_I4P
      do r = 1, nt
         if (find(r) /= r .or. touch(r) .or. pn(r) == 0_I4P .or. vol(r) < 0._R8P) cycle
         tgt(r) = inclusion_target(r)
         if (tgt(r) == 0_I4P) then
            npieces = -1_I4P ; return
         endif
      enddo
      do r = 1, nt
         if (tgt(r) == 0_I4P) cycle
         kk = tgt(r) ; guard = 0_I4P
         do while (tgt(kk) /= 0_I4P .and. guard < natom)
            kk = tgt(kk) ; guard = guard + 1_I4P
         enddo
         tgt(r) = kk
      enddo
      npc_ = 0_I4P
      do r = 1, natom
         if (find(r) /= r) cycle
         if (touch(r) .or. (r <= nt .and. tgt(r) == 0_I4P .and. vol(r) < 0._R8P)) then
            npc_ = npc_ + 1_I4P ; pid(r) = npc_
         endif
      enddo
      ! walls and traces
      do qq = 1, nt
         if (pn(qq) == 0_I4P) cycle
         r = find(qq) ; if (tgt(r) /= 0_I4P) r = tgt(r)
         do vv = 2, pn(qq) - 1
            call tri_terms(pv(:, 1, qq) - c0, pv(:, vv, qq) - c0, pv(:, vv + 1, qq) - c0, wt(:, r))
         enddo
         do ax = 1, 3
            b1 = mod(ax, 3_I4P) + 1_I4P ; b2 = mod(ax + 1_I4P, 3_I4P) + 1_I4P
            do sd = 0, 1
               bt = 2_I4P * (ax - 1_I4P) + sd ; sg = merge(1._R8P, -1._R8P, sd == 1_I4P)
               inplane = all(btest(pfl(1:pn(qq), qq), bt))
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
                  elseif (maxval(tri(ax, :, lst(qq))) <= merge(c1(ax), c0(ax), sd == 1_I4P)) then
                     ! the facet ends on the face plane, which the perturbed plane does not cut: its edge on the plane
                     ! bounds the fluid only if the facet across that edge crosses the plane (that facet's polygon in
                     ! the cell is the edge alone, so this one stands for it)
                     if (.not. across_crosses(lst(qq), ax, merge(c1(ax), c0(ax), sd == 1_I4P))) cycle
                  endif
                  pl = pv(:, vv, qq) - c0 ; ql = pv(:, jj, qq) - c0
                  fa(bt + 1, r) = fa(bt + 1, r) + sg * 0.5_R8P * (pl(b1) + ql(b1)) * (ql(b2) - pl(b2))
               enddo
            enddo
         enddo
      enddo
      ! fluid edge intervals: edge along k at sides s1 (axis k1), s2 (axis k2); on the face (k1, s1) it
      ! runs along v at u = s2 (counter-clockwise: +v at u = h); on the face (k2, s2) along u (no area)
      do jj = 1, niv
         ee = iv_edge(jj) ; kk = (ee - 1_I4P) / 4_I4P + 1_I4P
         s1 = mod(mod(ee - 1_I4P, 4_I4P), 2_I4P) ; s2 = mod(ee - 1_I4P, 4_I4P) / 2_I4P
         k1 = mod(kk, 3_I4P) + 1_I4P ; k2 = mod(kk + 1_I4P, 3_I4P) + 1_I4P
         if (s2 == 1_I4P) then
            r = find(iv_atom(jj))
            fa(2 * (k1 - 1) + s1 + 1, r) = fa(2 * (k1 - 1) + s1 + 1, r) + h(k2) * (iv_hi(jj) - iv_lo(jj))
         endif
      enddo
      allocate(pm(17, npc_))
      hf = [h(2) * h(3), h(3) * h(1), h(1) * h(2)]
      do r = 1, natom
         if (pid(r) == 0_I4P) cycle
         kk = pid(r)
         vp = sum(h * fa([2, 4, 6], r)) / 3._R8P - wt(1, r)
         mp = 0.5_R8P * h * h * fa([2, 4, 6], r) - wt(2:4, r)
         pm(1, kk) = vp
         pm(2:4, kk) = c0 + 0.5_R8P * h
         if (vp > 0._R8P) pm(2:4, kk) = c0 + mp / vp
         pm(5, kk) = wt(5, r)
         pm(6:8, kk) = c0 + 0.5_R8P * h
         if (wt(5, r) > 0._R8P) pm(6:8, kk) = c0 + wt(6:8, r) / wt(5, r)
         pm(9:11, kk) = wt(9:11, r)
         pm(12:17, kk) = fa(:, r) / [hf(1), hf(1), hf(2), hf(2), hf(3), hf(3)]
      enddo
      endsubroutine piece_measures

      logical function across_crosses(fi, ax, lev) result(yes)
      !< The facet of the cell across the edge of facet fi that lies on the plane x(ax) = lev has a vertex beyond the
      !< plane (x(ax) > lev); .false. if fi has no edge on the plane or no facet of the cell shares it.
      integer(I4P), intent(in) :: fi  !< Facet.
      integer(I4P), intent(in) :: ax  !< Plane axis.
      real(R8P),    intent(in) :: lev !< Plane level.
      integer(I4P)             :: m1, m2, g, mm, nshare !< Work.

      yes = .false.
      do m1 = 1, 3
         m2 = mod(m1, 3_I4P) + 1_I4P
         if (tri(ax, m1, fi) /= lev .or. tri(ax, m2, fi) /= lev) cycle
         do g = 1, nt
            if (lst(g) == fi) cycle
            nshare = 0_I4P
            do mm = 1, 3
               if (all(tri(:, mm, lst(g)) == tri(:, m1, fi)) .or. all(tri(:, mm, lst(g)) == tri(:, m2, fi))) &
                  nshare = nshare + 1_I4P
            enddo
            if (nshare == 2_I4P .and. maxval(tri(ax, :, lst(g))) > lev) then
               yes = .true. ; return
            endif
         enddo
      enddo
      endfunction across_crosses

      integer(I4P) function inclusion_target(ri) result(tr)
      !< Root of the set around the inclusion ri (see `piece_measures`); 0 if not found.
      integer(I4P), intent(in) :: ri !< Root of the inclusion.
      real(R8P)                :: pr(3), t(3, 3), nn(3), xh, xb, o(3), ub, uc !< Work.
      integer(I4P)             :: q2, v2, j2, qb, ed, sides(3) !< Work.

      pr = -huge(1._R8P)
      do q2 = 1, nt
         if (pn(q2) == 0_I4P .or. find(q2) /= ri) cycle
         do v2 = 1, pn(q2)
            if (pv(1, v2, q2) > pr(1)) pr = pv(:, v2, q2)
         enddo
      enddo
      xb = huge(1._R8P) ; qb = 0_I4P
      do q2 = 1, nt
         if (pn(q2) == 0_I4P .or. find(q2) == ri) cycle
         t = tri(:, :, lst(q2)) ; nn = nrm(:, lst(q2))
         if (nn(1) == 0._R8P) cycle
         if (orient2d_xy(t(2, 1), t(3, 1), t(2, 2), t(3, 2), t(2, 3), t(3, 3)) == 0._R8P) cycle
         o(1) = orient2d_xy(t(2, 1), t(3, 1), t(2, 2), t(3, 2), pr(2), pr(3))
         o(2) = orient2d_xy(t(2, 2), t(3, 2), t(2, 3), t(3, 3), pr(2), pr(3))
         o(3) = orient2d_xy(t(2, 3), t(3, 3), t(2, 1), t(3, 1), pr(2), pr(3))
         if (.not. (all(o >= 0._R8P) .or. all(o <= 0._R8P))) cycle
         xh = t(1, 1) - (nn(2) * (pr(2) - t(2, 1)) + nn(3) * (pr(3) - t(3, 1))) / nn(1)
         if (xh > pr(1) .and. xh <= c1(1) .and. xh < xb) then
            xb = xh ; qb = q2
         endif
      enddo
      if (qb /= 0_I4P) then
         tr = find(qb) ; return
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
         tr = find(qb) ; return
      endif
      sides = [1_I4P, 1_I4P, 0_I4P] ; ed = edge_index(3_I4P, sides)
      tr = 0_I4P
      do j2 = 1, niv
         if (iv_edge(j2) == ed .and. iv_lo(j2) <= pr(3) .and. pr(3) <= iv_hi(j2)) tr = find(iv_atom(j2))
      enddo
      endfunction inclusion_target
   endsubroutine cell_pieces

   pure logical function split_closes(pm, h, alpha, faces) result(ok)
   !< The pieces pm [17, n] of a cell of size h split it consistently (issue #70): their face fractions add up to the
   !< exact face totals (the plane sections), their volumes to alpha |cell|, and every piece closes, (f+ - f-) A_a = the
   !< a component of its wall normal integral, per axis. Tolerance 1e-12 of the face area or cell volume, as the other
   !< trust tests. A face split by the traces of the clipped polygons can miss a boundary on the perturbed face plane
   !< (a facet lying in it or ending on it next to one crossing it): such a cell is not trusted.
   real(R8P), intent(in) :: pm(:,:)   !< Piece measures (see `piece_measures`).
   real(R8P), intent(in) :: h(3)      !< Cell size.
   real(R8P), intent(in) :: alpha     !< Cell volume fraction.
   real(R8P), intent(in) :: faces(6)  !< Face fractions (-x, +x, -y, +y, -z, +z).
   real(R8P)             :: af(3), vc  !< Face areas, cell volume.
   integer(I4P)          :: a, p      !< Work.
   real(R8P), parameter  :: TOL = 1.e-12_R8P !< Relative tolerance.

   af = [h(2) * h(3), h(3) * h(1), h(1) * h(2)] ; vc = product(h)
   ok = abs(sum(pm(1, :)) - alpha * vc) <= TOL * vc
   do a = 1, 6
      ok = ok .and. abs(sum(pm(11 + a, :)) - faces(a)) <= TOL
   enddo
   do p = 1, size(pm, 2)
      do a = 1, 3
         ok = ok .and. abs((pm(13 + 2 * (a - 1), p) - pm(12 + 2 * (a - 1), p)) * af(a) - pm(8 + a, p)) <= TOL * af(a)
      enddo
   enddo
   endfunction split_closes

   subroutine column_parity(surface, xa, yb, z, inside)
   !< Exact inside flags of the nodes (xa, yb, z(k)) of a grid column (issue #70): the parity of the crossings of
   !< the z line before each node, decided by `line_crossing_before` under the global perturbation, a node on the
   !< wall taken at P + (eps1, eps2, eps3). The same test as the node parities of #69, so the corner states of
   !< every cell agree with the crossings of its edges.
   type(surface_stl_object), intent(in)  :: surface   !< Surface.
   real(R8P),                intent(in)  :: xa, yb    !< Column.
   real(R8P),                intent(in)  :: z(:)      !< Nodes.
   logical,                  intent(out) :: inside(:) !< Inside flags.
   real(R8P),    allocatable :: t(:)     !< Crossings.
   integer(I4P), allocatable :: fc(:)    !< Facets.
   logical,      allocatable :: en(:)    !< Entering flags.
   integer(I4P)              :: n, k, m, kk !< Work.

   call surface%line_crossings(3_I4P, xa, yb, -huge(1._R8P), huge(1._R8P), t, fc, en, n)
   do k = 1, size(z)
      m = 0_I4P
      do kk = 1_I4P, n
         if (line_crossing_before(reshape(surface%tri_vertices(:, fc(kk)), [3, 3]), 3_I4P, xa, yb, z(k), en(kk))) &
            m = m + 1_I4P
      enddo
      inside(k) = mod(m, 2_I4P) == 1_I4P
   enddo
   endsubroutine column_parity

   pure integer(I4P) function edge_index(k, sides)
   !< Cell edge (1-12) along axis k at the sides (0 low, 1 high) of the two other axes.
   integer(I4P), intent(in) :: k        !< Edge axis.
   integer(I4P), intent(in) :: sides(3) !< Sides per axis (sides(k) ignored).

   edge_index = (k - 1_I4P) * 4_I4P + sides(mod(k, 3_I4P) + 1_I4P) + 2_I4P * sides(mod(k + 1_I4P, 3_I4P) + 1_I4P) + 1_I4P
   endfunction edge_index

   pure integer(I4P) function corner_of(e, side)
   !< Corner (1-8, sides x + 2 y + 4 z + 1) at the low (side 0) or high (1) end of edge e.
   integer(I4P), intent(in) :: e, side !< Edge, end.
   integer(I4P)             :: k, r, sides(3) !< Work.

   k = (e - 1_I4P) / 4_I4P + 1_I4P ; r = mod(e - 1_I4P, 4_I4P)
   sides(k) = side ; sides(mod(k, 3_I4P) + 1_I4P) = mod(r, 2_I4P) ; sides(mod(k + 1_I4P, 3_I4P) + 1_I4P) = r / 2_I4P
   corner_of = sides(1) + 2_I4P * sides(2) + 4_I4P * sides(3) + 1_I4P
   endfunction corner_of

   pure integer(I4P) function edge_of(fl)
   !< Cell edge of a vertex with on-plane flags fl (bit 2 (a - 1) + side): exactly two axes flagged,
   !< else 0.
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

   pure subroutine sort_pierce(pc, n)
   !< Sort piercings by edge, then position (insertion).
   type(pierce_t), intent(inout) :: pc(:) !< Piercings.
   integer(I4P),   intent(in)    :: n     !< Their number.
   type(pierce_t)                :: v     !< Work.
   integer(I4P)                  :: i, j  !< Work.

   do i = 2, n
      v = pc(i) ; j = i - 1_I4P
      do while (j >= 1_I4P)
         if (pc(j)%edge < v%edge .or. (pc(j)%edge == v%edge .and. pc(j)%t <= v%t)) exit
         pc(j + 1) = pc(j) ; j = j - 1_I4P
      enddo
      pc(j + 1) = v
   enddo
   endsubroutine sort_pierce

   pure logical function lexless(a, b)
   !< Lexicographic order of two points.
   real(R8P), intent(in) :: a(3), b(3) !< Points.

   lexless = a(1) < b(1) .or. (a(1) == b(1) .and. (a(2) < b(2) .or. (a(2) == b(2) .and. a(3) < b(3))))
   endfunction lexless

   subroutine sort_points(key, idx)
   !< Index heapsort of the keys key(:, m), lexicographic over the rows.
   real(R8P),    intent(in)  :: key(:,:) !< Points [3, m].
   integer(I4P), intent(out) :: idx(:)   !< Sorted order.
   integer(I4P)              :: m, i, t  !< Work.

   m = size(idx, kind=I4P)
   idx = [(i, i = 1, m)]
   do i = m / 2_I4P, 1_I4P, -1_I4P
      call sift(i, m)
   enddo
   do i = m, 2_I4P, -1_I4P
      t = idx(1) ; idx(1) = idx(i) ; idx(i) = t
      call sift(1_I4P, i - 1_I4P)
   enddo
   contains
      subroutine sift(lo, hi)
      !< Sift down idx(lo) in the heap idx(1:hi).
      integer(I4P), intent(in) :: lo, hi !< Root, heap end.
      integer(I4P)             :: r, c, tt !< Work.

      r = lo
      do while (2_I4P * r <= hi)
         c = 2_I4P * r
         if (c < hi) then
            if (less(idx(c), idx(c + 1))) c = c + 1_I4P
         endif
         if (.not. less(idx(r), idx(c))) exit
         tt = idx(r) ; idx(r) = idx(c) ; idx(c) = tt ; r = c
      enddo
      endsubroutine sift

      logical function less(a, b)
      !< Lexicographic order of the keys a, b.
      integer(I4P), intent(in) :: a, b !< Keys.
      integer(I4P)             :: r    !< Row.

      less = .false.
      do r = 1, size(key, 1, kind=I4P)
         if (key(r, a) /= key(r, b)) then
            less = key(r, a) < key(r, b)
            return
         endif
      enddo
      endfunction less
   endsubroutine sort_points

   pure function cross3(a, b) result(c)
   !< Cross product.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< a x b.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   pure subroutine clip_flags(p, fl, n, ax, level, keep_above, bit)
   !< `clip_plane` (same kept points and crossings) carrying on-plane flags: a crossing gets the bit of
   !< the plane plus the flags its two end points share; kept points keep theirs.
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
      ik = kept(p(:, k)) ; ikk = kept(p(:, kk))
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
   contains
      pure logical function kept(q)
      !< Kept side.
      real(R8P), intent(in) :: q(3) !< Point.

      if (keep_above) then
         kept = q(ax) > level
      else
         kept = q(ax) <= level
      endif
      endfunction kept
   endsubroutine clip_flags

   pure subroutine tri_terms(p1, p2, p3, a)
   !< Divergence-theorem terms of the triangle p1 p2 p3 (local coordinates, outward orientation):
   !< a(1) += 1/3 p1 . A (volume), a(2:4) += 1/2 A_k int x_k^2 dA / area (moments), a(5) += area,
   !< a(6:8) += area * centroid, a(9:11) += A (vector area).
   real(R8P), intent(in)    :: p1(3), p2(3), p3(3) !< Vertices.
   real(R8P), intent(inout) :: a(11)               !< Accumulators.
   real(R8P)                :: va(3), ar, q(3)     !< Vector area, area, quadratic sums.

   va = 0.5_R8P * [(p2(2) - p1(2)) * (p3(3) - p1(3)) - (p2(3) - p1(3)) * (p3(2) - p1(2)), &
                   (p2(3) - p1(3)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(3) - p1(3)), &
                   (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))]
   ar = norm2(va)
   q = (p1 * p1 + p2 * p2 + p3 * p3 + p1 * p2 + p2 * p3 + p3 * p1) / 6._R8P ! int x_k^2 dA / area
   a(1) = a(1) + dot_product(p1, va) / 3._R8P
   a(2:4) = a(2:4) + 0.5_R8P * va * q
   a(5) = a(5) + ar
   a(6:8) = a(6:8) + ar * (p1 + p2 + p3) / 3._R8P
   a(9:11) = a(9:11) + va
   endsubroutine tri_terms

   pure subroutine clip_plane(p, n, ax, level, keep_above)
   !< Clip the polygon p(:, 1:n) by the plane x(ax) = level: keep x > level (keep_above) or x <= level.
   !< Crossing points are canonical (from the edge end points in lexicographic order, the axis
   !< coordinate set to level), so neighbouring cells get the same points.
   real(R8P),    intent(inout) :: p(3, 16)   !< Polygon.
   integer(I4P), intent(inout) :: n          !< Vertices.
   integer,      intent(in)    :: ax         !< Axis.
   real(R8P),    intent(in)    :: level      !< Plane.
   logical,      intent(in)    :: keep_above !< Side kept.
   real(R8P)                   :: o(3, 16)   !< Output.
   integer(I4P)                :: k, kk, no  !< Work.
   logical                     :: ik, ikk    !< Inside flags.

   if (n == 0) return
   no = 0
   do k = 1, n
      kk = mod(k, n) + 1
      ik = inside(p(:, k)) ; ikk = inside(p(:, kk))
      if (ik) then
         no = no + 1 ; o(:, no) = p(:, k)
      endif
      if (ik .neqv. ikk) then
         no = no + 1 ; o(:, no) = crossing(p(:, k), p(:, kk))
      endif
   enddo
   n = no ; p(:, 1:n) = o(:, 1:n)
   contains
      pure logical function inside(q)
      !< Kept side.
      real(R8P), intent(in) :: q(3) !< Point.

      if (keep_above) then
         inside = q(ax) > level
      else
         inside = q(ax) <= level
      endif
      endfunction inside

      pure function crossing(a, b) result(c)
      !< Canonical crossing of the segment ab with the plane.
      real(R8P), intent(in) :: a(3), b(3) !< End points.
      real(R8P)             :: c(3)       !< Crossing.
      real(R8P)             :: lo(3), hi(3) !< Ordered end points.

      if (a(1) < b(1) .or. (a(1) == b(1) .and. (a(2) < b(2) .or. (a(2) == b(2) .and. a(3) <= b(3))))) then
         lo = a ; hi = b
      else
         lo = b ; hi = a
      endif
      c = lo + (level - lo(ax)) / (hi(ax) - lo(ax)) * (hi - lo)
      c(ax) = level
      endfunction crossing
   endsubroutine clip_plane

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
endmodule fossil_cut_cells_3d
