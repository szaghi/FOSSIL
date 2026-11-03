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
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_curve_slice,        only : slice_surface
use fossil_facet_object,       only : facet_object
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
endtype cut_cells_3d_t

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
   endsubroutine cut_cells_3d

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
