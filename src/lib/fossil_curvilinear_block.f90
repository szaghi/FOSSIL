!< FOSSIL, curvilinear structured blocks: vertices, topology records and their detection (issue #71).
module fossil_curvilinear_block
!< FOSSIL, curvilinear structured blocks: vertices, topology records and their detection (issue #71).
!<
!< A block of ni x nj x nk vertices (2-D: ni x nj) stored as structure-of-arrays x(ni, nj, nk), y, z (2-D: u, v).
!< Cell (i, j, k) is spanned by the vertices i..i+1, j..j+1, k..k+1. Faces are numbered 1 to 6 (-i, +i, -j, +j, -k,
!< +k; 1 to 4 in 2-D). Design: docs/guide/advanced/curvilinear-blocks.md.
!<
!< Topology records (`topo_record_t`) live on the faces, over a vertex range: a named physical boundary
!< (TOPO_BOUNDARY), an abutting range (TOPO_ABUT: coincident vertex by vertex with a range of a face of this block or
!< another one), a periodic range (TOPO_PERIODIC: the same after the rigid motion x -> R x + t) or a singular face
!< (TOPO_SINGULAR: collapsed along direction `collapsed`). An abutting or periodic record maps the vertex p to
!< donor_begin + M (p - begin), M(:, d) = sign(T(d)) e_|T(d)| (the CGNS transform), off the face too: one layer outside
!< this block maps to the donor's first interior layer.
!<
!< `detect_topology` (`detect_topology_2d`) proposes the records of a set of blocks: boundary vertices sorted by x
!< and matched in a window of the tolerance (exact coincidence by default), singular faces first, then each uncovered
!< vertex whose image is unique, with unique images of its face neighbours, seeds a transform that grows to the
!< largest rectangle of the face on which every vertex maps onto a coincident vertex.
use penf, only : I4P, R8P
use fossil_block_geometry, only : hex_volume, hex_centroid, bilinear_face_vector, quad_area, quad_centroid, hex_validity, &
                                  quad_validity

implicit none
private
public :: curvilinear_block_t, curvilinear_block_2d_t, topo_record_t
public :: detect_topology, detect_topology_2d, record_map
public :: TOPO_BOUNDARY, TOPO_ABUT, TOPO_PERIODIC, TOPO_SINGULAR
public :: BLOCK_STATUS_OK, BLOCK_STATUS_BAD_INPUT

integer(I4P), parameter :: TOPO_BOUNDARY = 1_I4P !< Named physical boundary.
integer(I4P), parameter :: TOPO_ABUT     = 2_I4P !< Abutting range (seam, cut, interface).
integer(I4P), parameter :: TOPO_PERIODIC = 3_I4P !< Periodic range (rigid motion).
integer(I4P), parameter :: TOPO_SINGULAR = 4_I4P !< Singular face (collapsed along one direction).

integer(I4P), parameter :: BLOCK_STATUS_OK        = 0_I4P !< Success.
integer(I4P), parameter :: BLOCK_STATUS_BAD_INPUT = 1_I4P !< Bad sizes or non-finite coordinates.

type :: topo_record_t
   !< A topology record on a face of a block.
   integer(I4P)      :: kind = 0_I4P                  !< TOPO_*.
   integer(I4P)      :: face = 0_I4P                  !< Face 1..6 (1..4 in 2-D).
   integer(I4P)      :: lo(3) = 1_I4P                 !< Vertex range, low corner (normal index fixed).
   integer(I4P)      :: hi(3) = 1_I4P                 !< Vertex range, high corner.
   character(len=32) :: name = ''                     !< Name (boundary patches).
   integer(I4P)      :: donor_block = 0_I4P           !< Donor block (index in the set).
   integer(I4P)      :: donor_face = 0_I4P            !< Donor face.
   integer(I4P)      :: begin(3) = 1_I4P              !< A vertex of this range...
   integer(I4P)      :: donor_begin(3) = 1_I4P        !< ...and its image on the donor.
   integer(I4P)      :: transform(3) = [1, 2, 3]      !< CGNS transform.
   integer(I4P)      :: collapsed = 0_I4P             !< Singular: the collapsed direction.
   real(R8P)         :: rotation(3, 3) = reshape([1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], &
                                                 [3, 3]) !< Periodic: donor point = rotation x + translation.
   real(R8P)         :: translation(3) = 0._R8P       !< Periodic translation.
endtype topo_record_t

type :: curvilinear_block_t
   !< Curvilinear structured block in 3-D.
   integer(I4P)                     :: ni = 0_I4P, nj = 0_I4P, nk = 0_I4P !< Vertices per direction.
   character(len=:),    allocatable :: name                   !< Block name.
   real(R8P),           allocatable :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   integer(I4P)                     :: nrec = 0_I4P           !< Number of records.
   type(topo_record_t), allocatable :: rec(:)                 !< Topology records [nrec].
   contains
      procedure, pass(self) :: initialize => initialize_3d    !< Set the vertices.
      procedure, pass(self) :: destroy => destroy_3d          !< Reset.
      procedure, pass(self) :: add_record => add_record_3d    !< Append a record.
      procedure, pass(self) :: cells_number => cells_number_3d !< Number of cells.
      procedure, pass(self) :: corners => corners_3d          !< Corners of a cell.
      procedure, pass(self) :: geometry => geometry_3d        !< Exact cell geometry.
      procedure, pass(self) :: orientation => orientation_3d  !< Handedness (sign of the block volume).
      procedure, pass(self) :: validity => validity_3d        !< Certified cell validity.
endtype curvilinear_block_t

type :: curvilinear_block_2d_t
   !< Curvilinear structured block in 2-D.
   integer(I4P)                     :: ni = 0_I4P, nj = 0_I4P !< Vertices per direction.
   character(len=:),    allocatable :: name                   !< Block name.
   real(R8P),           allocatable :: u(:,:), v(:,:)         !< Vertices.
   integer(I4P)                     :: nrec = 0_I4P           !< Number of records.
   type(topo_record_t), allocatable :: rec(:)                 !< Topology records [nrec].
   contains
      procedure, pass(self) :: initialize => initialize_2d    !< Set the vertices.
      procedure, pass(self) :: destroy => destroy_2d          !< Reset.
      procedure, pass(self) :: add_record => add_record_2d    !< Append a record.
      procedure, pass(self) :: cells_number => cells_number_2d !< Number of cells.
      procedure, pass(self) :: corners => corners_2d          !< Corners of a cell.
      procedure, pass(self) :: geometry => geometry_2d        !< Exact cell geometry.
      procedure, pass(self) :: orientation => orientation_2d  !< Handedness (sign of the block area).
      procedure, pass(self) :: validity => validity_2d        !< Exact cell validity.
endtype curvilinear_block_2d_t

type :: vx_t
   !< Work copy of a block's vertices for the detection (2-D: nk = 1, z = 0).
   integer(I4P)           :: n(3) = 1_I4P
   real(R8P), allocatable :: c(:,:,:,:) !< Coordinates [3, ni, nj, nk].
endtype vx_t

type :: bv_t
   !< Sorted boundary vertices of a set of blocks.
   integer(I4P)              :: m = 0_I4P
   integer(I4P), allocatable :: id(:,:)  !< [block, i, j, k] per entry, sorted by coordinates.
   real(R8P),    allocatable :: xs(:)    !< x of the sorted entries.
endtype bv_t

contains
   ! ---------------------------------------------------------------------------------------------------------- blocks
   subroutine initialize_3d(self, x, y, z, name, status)
   !< Set the vertices (arrays of the same shape, at least 2 x 2 x 2, finite).
   class(curvilinear_block_t), intent(inout)         :: self    !< Block.
   real(R8P),                  intent(in)            :: x(:,:,:), y(:,:,:), z(:,:,:) !< Vertices.
   character(*),               intent(in),  optional :: name    !< Name.
   integer(I4P),               intent(out), optional :: status  !< BLOCK_STATUS_*.
   integer(I4P)                                      :: st      !< Status.

   call self%destroy
   st = BLOCK_STATUS_OK
   if (any(shape(x) < 2) .or. any(shape(y) /= shape(x)) .or. any(shape(z) /= shape(x))) then
      st = BLOCK_STATUS_BAD_INPUT
   elseif (.not. (all(abs(x) <= huge(1._R8P)) .and. all(abs(y) <= huge(1._R8P)) .and. all(abs(z) <= huge(1._R8P)))) then
      st = BLOCK_STATUS_BAD_INPUT
   endif
   if (st == BLOCK_STATUS_OK) then
      self%ni = size(x, 1, I4P) ; self%nj = size(x, 2, I4P) ; self%nk = size(x, 3, I4P)
      self%x = x ; self%y = y ; self%z = z
      allocate(self%rec(0))
      self%name = 'block' ; if (present(name)) self%name = name
   endif
   if (present(status)) status = st
   endsubroutine initialize_3d

   subroutine initialize_2d(self, u, v, name, status)
   !< Set the vertices (arrays of the same shape, at least 2 x 2, finite).
   class(curvilinear_block_2d_t), intent(inout)         :: self   !< Block.
   real(R8P),                     intent(in)            :: u(:,:), v(:,:) !< Vertices.
   character(*),                  intent(in),  optional :: name   !< Name.
   integer(I4P),                  intent(out), optional :: status !< BLOCK_STATUS_*.
   integer(I4P)                                         :: st     !< Status.

   call self%destroy
   st = BLOCK_STATUS_OK
   if (any(shape(u) < 2) .or. any(shape(v) /= shape(u))) then
      st = BLOCK_STATUS_BAD_INPUT
   elseif (.not. (all(abs(u) <= huge(1._R8P)) .and. all(abs(v) <= huge(1._R8P)))) then
      st = BLOCK_STATUS_BAD_INPUT
   endif
   if (st == BLOCK_STATUS_OK) then
      self%ni = size(u, 1, I4P) ; self%nj = size(u, 2, I4P)
      self%u = u ; self%v = v
      allocate(self%rec(0))
      self%name = 'block' ; if (present(name)) self%name = name
   endif
   if (present(status)) status = st
   endsubroutine initialize_2d

   elemental subroutine destroy_3d(self)
   !< Reset.
   class(curvilinear_block_t), intent(inout) :: self  !< Block.
   type(curvilinear_block_t)                 :: fresh !< Reset value.

   select type(self)
   type is(curvilinear_block_t)
      self = fresh
   endselect
   endsubroutine destroy_3d

   elemental subroutine destroy_2d(self)
   !< Reset.
   class(curvilinear_block_2d_t), intent(inout) :: self  !< Block.
   type(curvilinear_block_2d_t)                 :: fresh !< Reset value.

   select type(self)
   type is(curvilinear_block_2d_t)
      self = fresh
   endselect
   endsubroutine destroy_2d

   pure subroutine add_record_3d(self, rec)
   !< Append a record.
   class(curvilinear_block_t), intent(inout) :: self !< Block.
   type(topo_record_t),        intent(in)    :: rec  !< Record.

   if (.not. allocated(self%rec)) allocate(self%rec(0))
   self%rec = [self%rec, rec] ; self%nrec = size(self%rec, kind=I4P)
   endsubroutine add_record_3d

   pure subroutine add_record_2d(self, rec)
   !< Append a record.
   class(curvilinear_block_2d_t), intent(inout) :: self !< Block.
   type(topo_record_t),           intent(in)    :: rec  !< Record.

   if (.not. allocated(self%rec)) allocate(self%rec(0))
   self%rec = [self%rec, rec] ; self%nrec = size(self%rec, kind=I4P)
   endsubroutine add_record_2d

   elemental integer(I4P) function cells_number_3d(self) result(n)
   !< Number of cells.
   class(curvilinear_block_t), intent(in) :: self !< Block.

   n = max(self%ni - 1_I4P, 0_I4P) * max(self%nj - 1_I4P, 0_I4P) * max(self%nk - 1_I4P, 0_I4P)
   endfunction cells_number_3d

   elemental integer(I4P) function cells_number_2d(self) result(n)
   !< Number of cells.
   class(curvilinear_block_2d_t), intent(in) :: self !< Block.

   n = max(self%ni - 1_I4P, 0_I4P) * max(self%nj - 1_I4P, 0_I4P)
   endfunction cells_number_2d

   pure function record_map(rec, p) result(q)
   !< Image of the vertex index p (on or off the face) through an abutting or periodic record.
   type(topo_record_t), intent(in) :: rec  !< Record.
   integer(I4P),        intent(in) :: p(3) !< Vertex index.
   integer(I4P)                    :: q(3) !< Donor vertex index.
   integer(I4P)                    :: d    !< Direction.

   q = rec%donor_begin
   do d = 1, 3
      q(abs(rec%transform(d))) = q(abs(rec%transform(d))) + sign(1_I4P, rec%transform(d)) * (p(d) - rec%begin(d))
   enddo
   endfunction record_map

   ! -------------------------------------------------------------------------------------------------------- geometry
   pure function corners_3d(self, i, j, k) result(c)
   !< Corners of cell (i, j, k): corner (a, b, e) at 1 + a + 2 b + 4 e.
   class(curvilinear_block_t), intent(in) :: self    !< Block.
   integer(I4P),               intent(in) :: i, j, k !< Cell.
   real(R8P)                              :: c(3, 8) !< Corners.
   integer(I4P)                           :: a, b, e !< Corner bits.

   do e = 0, 1
      do b = 0, 1
         do a = 0, 1
            c(:, 1 + a + 2 * b + 4 * e) = [self%x(i + a, j + b, k + e), self%y(i + a, j + b, k + e), &
                                           self%z(i + a, j + b, k + e)]
         enddo
      enddo
   enddo
   endfunction corners_3d

   pure function corners_2d(self, i, j) result(c)
   !< Corners of cell (i, j): corner (a, b) at 1 + a + 2 b.
   class(curvilinear_block_2d_t), intent(in) :: self    !< Block.
   integer(I4P),                  intent(in) :: i, j    !< Cell.
   real(R8P)                                 :: c(2, 4) !< Corners.

   c(:, 1) = [self%u(i, j), self%v(i, j)] ; c(:, 2) = [self%u(i + 1, j), self%v(i + 1, j)]
   c(:, 3) = [self%u(i, j + 1), self%v(i, j + 1)] ; c(:, 4) = [self%u(i + 1, j + 1), self%v(i + 1, j + 1)]
   endfunction corners_2d

   pure function vtx(self, i, j, k) result(p)
   !< Vertex (i, j, k).
   class(curvilinear_block_t), intent(in) :: self    !< Block.
   integer(I4P),               intent(in) :: i, j, k !< Vertex.
   real(R8P)                              :: p(3)    !< Coordinates.

   p = [self%x(i, j, k), self%y(i, j, k), self%z(i, j, k)]
   endfunction vtx

   subroutine geometry_3d(self, volume, centroid, face_i, face_j, face_k)
   !< Exact cell geometry: volumes, centroids, and the area vectors of the i, j, k faces, oriented towards +i, +j, +k
   !< (outward on the high faces of a right-handed cell). Cell (i, j, k) is bounded by face_i(:, i..i+1, j, k), etc.
   class(curvilinear_block_t), intent(in)                         :: self            !< Block.
   real(R8P), allocatable,     intent(out), optional              :: volume(:,:,:)   !< Volumes [ni-1, nj-1, nk-1].
   real(R8P), allocatable,     intent(out), optional              :: centroid(:,:,:,:) !< Centroids [3, cells].
   real(R8P), allocatable,     intent(out), optional              :: face_i(:,:,:,:) !< i faces [3, ni, nj-1, nk-1].
   real(R8P), allocatable,     intent(out), optional              :: face_j(:,:,:,:) !< j faces [3, ni-1, nj, nk-1].
   real(R8P), allocatable,     intent(out), optional              :: face_k(:,:,:,:) !< k faces [3, ni-1, nj-1, nk].
   integer(I4P)                                                   :: i, j, k, ni, nj, nk !< Counters.

   ni = self%ni ; nj = self%nj ; nk = self%nk
   if (present(volume)) then
      allocate(volume(ni - 1, nj - 1, nk - 1))
      do k = 1, nk - 1
         do j = 1, nj - 1
            do i = 1, ni - 1
               volume(i, j, k) = hex_volume(self%corners(i, j, k))
            enddo
         enddo
      enddo
   endif
   if (present(centroid)) then
      allocate(centroid(3, ni - 1, nj - 1, nk - 1))
      do k = 1, nk - 1
         do j = 1, nj - 1
            do i = 1, ni - 1
               centroid(:, i, j, k) = hex_centroid(self%corners(i, j, k))
            enddo
         enddo
      enddo
   endif
   if (present(face_i)) then
      allocate(face_i(3, ni, nj - 1, nk - 1))
      do k = 1, nk - 1
         do j = 1, nj - 1
            do i = 1, ni
               face_i(:, i, j, k) = bilinear_face_vector(vtx(self, i, j, k), vtx(self, i, j + 1, k), &
                                                         vtx(self, i, j + 1, k + 1), vtx(self, i, j, k + 1))
            enddo
         enddo
      enddo
   endif
   if (present(face_j)) then
      allocate(face_j(3, ni - 1, nj, nk - 1))
      do k = 1, nk - 1
         do j = 1, nj
            do i = 1, ni - 1
               face_j(:, i, j, k) = bilinear_face_vector(vtx(self, i, j, k), vtx(self, i, j, k + 1), &
                                                         vtx(self, i + 1, j, k + 1), vtx(self, i + 1, j, k))
            enddo
         enddo
      enddo
   endif
   if (present(face_k)) then
      allocate(face_k(3, ni - 1, nj - 1, nk))
      do k = 1, nk
         do j = 1, nj - 1
            do i = 1, ni - 1
               face_k(:, i, j, k) = bilinear_face_vector(vtx(self, i, j, k), vtx(self, i + 1, j, k), &
                                                         vtx(self, i + 1, j + 1, k), vtx(self, i, j + 1, k))
            enddo
         enddo
      enddo
   endif
   endsubroutine geometry_3d

   subroutine geometry_2d(self, area, centroid, face_i, face_j)
   !< Exact cell geometry: areas, centroids, and the normals (length = edge length) of the i and j edges, oriented
   !< towards +i, +j.
   class(curvilinear_block_2d_t), intent(in)            :: self            !< Block.
   real(R8P), allocatable,        intent(out), optional :: area(:,:)       !< Areas [ni-1, nj-1].
   real(R8P), allocatable,        intent(out), optional :: centroid(:,:,:) !< Centroids [2, cells].
   real(R8P), allocatable,        intent(out), optional :: face_i(:,:,:)   !< i edges [2, ni, nj-1].
   real(R8P), allocatable,        intent(out), optional :: face_j(:,:,:)   !< j edges [2, ni-1, nj].
   integer(I4P)                                         :: i, j, ni, nj    !< Counters.

   ni = self%ni ; nj = self%nj
   if (present(area)) then
      allocate(area(ni - 1, nj - 1))
      do j = 1, nj - 1
         do i = 1, ni - 1
            area(i, j) = quad_area(self%corners(i, j))
         enddo
      enddo
   endif
   if (present(centroid)) then
      allocate(centroid(2, ni - 1, nj - 1))
      do j = 1, nj - 1
         do i = 1, ni - 1
            centroid(:, i, j) = quad_centroid(self%corners(i, j))
         enddo
      enddo
   endif
   if (present(face_i)) then
      allocate(face_i(2, ni, nj - 1))
      do j = 1, nj - 1
         do i = 1, ni
            face_i(:, i, j) = [self%v(i, j + 1) - self%v(i, j), -(self%u(i, j + 1) - self%u(i, j))]
         enddo
      enddo
   endif
   if (present(face_j)) then
      allocate(face_j(2, ni - 1, nj))
      do j = 1, nj
         do i = 1, ni - 1
            face_j(:, i, j) = [-(self%v(i + 1, j) - self%v(i, j)), self%u(i + 1, j) - self%u(i, j)]
         enddo
      enddo
   endif
   endsubroutine geometry_2d

   real(R8P) function orientation_3d(self) result(o)
   !< Handedness: +1 if the block volume (sum of the cell volumes) is >= 0, else -1.
   class(curvilinear_block_t), intent(in) :: self    !< Block.
   real(R8P)                              :: vol     !< Volume.
   integer(I4P)                           :: i, j, k !< Counters.

   vol = 0._R8P
   do k = 1, self%nk - 1
      do j = 1, self%nj - 1
         do i = 1, self%ni - 1
            vol = vol + hex_volume(self%corners(i, j, k))
         enddo
      enddo
   enddo
   o = merge(1._R8P, -1._R8P, vol >= 0._R8P)
   endfunction orientation_3d

   real(R8P) function orientation_2d(self) result(o)
   !< Handedness: +1 if the block area (sum of the cell areas) is >= 0, else -1.
   class(curvilinear_block_2d_t), intent(in) :: self !< Block.
   real(R8P)                                 :: ar   !< Area.
   integer(I4P)                              :: i, j !< Counters.

   ar = 0._R8P
   do j = 1, self%nj - 1
      do i = 1, self%ni - 1
         ar = ar + quad_area(self%corners(i, j))
      enddo
   enddo
   o = merge(1._R8P, -1._R8P, ar >= 0._R8P)
   endfunction orientation_2d

   subroutine validity_3d(self, status, min_jacobian, orientation)
   !< Certified validity of every cell (CELL_*) and the certified lower bound of det J times the block orientation.
   class(curvilinear_block_t), intent(in)            :: self             !< Block.
   integer(I4P), allocatable,  intent(out)           :: status(:,:,:)    !< CELL_* [ni-1, nj-1, nk-1].
   real(R8P),    allocatable,  intent(out), optional :: min_jacobian(:,:,:) !< Lower bound of orient det J.
   real(R8P),                  intent(in),  optional :: orientation      !< Orientation (default: the block's).
   real(R8P)                                         :: o, jm            !< Orientation, bound.
   integer(I4P)                                      :: i, j, k          !< Counters.

   o = 0._R8P ; if (present(orientation)) o = sign(1._R8P, orientation)
   if (o == 0._R8P) o = self%orientation()
   allocate(status(self%ni - 1, self%nj - 1, self%nk - 1))
   if (present(min_jacobian)) allocate(min_jacobian(self%ni - 1, self%nj - 1, self%nk - 1))
   do k = 1, self%nk - 1
      do j = 1, self%nj - 1
         do i = 1, self%ni - 1
            call hex_validity(self%corners(i, j, k), o, status(i, j, k), jm)
            if (present(min_jacobian)) min_jacobian(i, j, k) = jm
         enddo
      enddo
   enddo
   endsubroutine validity_3d

   subroutine validity_2d(self, status, min_jacobian, orientation)
   !< Exact validity of every cell (CELL_*) and the minimum of det J times the block orientation (at a corner).
   class(curvilinear_block_2d_t), intent(in)            :: self             !< Block.
   integer(I4P), allocatable,     intent(out)           :: status(:,:)      !< CELL_* [ni-1, nj-1].
   real(R8P),    allocatable,     intent(out), optional :: min_jacobian(:,:) !< min of orient det J.
   real(R8P),                     intent(in),  optional :: orientation      !< Orientation (default: the block's).
   real(R8P)                                            :: o, jm            !< Orientation, minimum.
   integer(I4P)                                         :: i, j             !< Counters.

   o = 0._R8P ; if (present(orientation)) o = sign(1._R8P, orientation)
   if (o == 0._R8P) o = self%orientation()
   allocate(status(self%ni - 1, self%nj - 1))
   if (present(min_jacobian)) allocate(min_jacobian(self%ni - 1, self%nj - 1))
   do j = 1, self%nj - 1
      do i = 1, self%ni - 1
         call quad_validity(self%corners(i, j), o, status(i, j), jm)
         if (present(min_jacobian)) min_jacobian(i, j) = jm
      enddo
   enddo
   endsubroutine validity_2d

   ! ------------------------------------------------------------------------------------------------------- detection
   subroutine detect_topology(blocks, tolerance, translations)
   !< Propose the topology records of a set of 3-D blocks (appended to each block's records).
   type(curvilinear_block_t), intent(inout)        :: blocks(:)         !< Blocks.
   real(R8P),                 intent(in), optional :: tolerance         !< Coincidence in the max norm (default exact).
   real(R8P),                 intent(in), optional :: translations(:,:) !< Periodic translations [3, nt].
   type(vx_t),          allocatable                :: vx(:)             !< Work copies.
   type(topo_record_t), allocatable                :: found(:)          !< Records found.
   integer(I4P),        allocatable                :: owner(:)          !< Their blocks.
   integer(I4P)                                    :: b, r              !< Counters.

   allocate(vx(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      vx(b)%n = [blocks(b)%ni, blocks(b)%nj, blocks(b)%nk]
      allocate(vx(b)%c(3, vx(b)%n(1), vx(b)%n(2), vx(b)%n(3)))
      vx(b)%c(1,:,:,:) = blocks(b)%x ; vx(b)%c(2,:,:,:) = blocks(b)%y ; vx(b)%c(3,:,:,:) = blocks(b)%z
   enddo
   call detect_core(vx, 3_I4P, tolerance, translations, found, owner)
   do r = 1, size(found, kind=I4P)
      call blocks(owner(r))%add_record(found(r))
   enddo
   endsubroutine detect_topology

   subroutine detect_topology_2d(blocks, tolerance, translations)
   !< Propose the topology records of a set of 2-D blocks (appended to each block's records).
   type(curvilinear_block_2d_t), intent(inout)        :: blocks(:)         !< Blocks.
   real(R8P),                    intent(in), optional :: tolerance         !< Coincidence in the max norm (default exact).
   real(R8P),                    intent(in), optional :: translations(:,:) !< Periodic translations [2, nt].
   type(vx_t),          allocatable                   :: vx(:)             !< Work copies.
   type(topo_record_t), allocatable                   :: found(:)          !< Records found.
   integer(I4P),        allocatable                   :: owner(:)          !< Their blocks.
   real(R8P),           allocatable                   :: t3(:,:)           !< Translations in 3-D.
   integer(I4P)                                       :: b, r              !< Counters.

   allocate(vx(size(blocks)))
   do b = 1, size(blocks, kind=I4P)
      vx(b)%n = [blocks(b)%ni, blocks(b)%nj, 1_I4P]
      allocate(vx(b)%c(3, vx(b)%n(1), vx(b)%n(2), 1))
      vx(b)%c(1,:,:,1) = blocks(b)%u ; vx(b)%c(2,:,:,1) = blocks(b)%v ; vx(b)%c(3,:,:,1) = 0._R8P
   enddo
   if (present(translations)) then
      allocate(t3(3, size(translations, 2))) ; t3 = 0._R8P ; t3(1:2,:) = translations(1:2,:)
      call detect_core(vx, 2_I4P, tolerance, t3, found, owner)
   else
      call detect_core(vx, 2_I4P, tolerance, found=found, owner=owner)
   endif
   do r = 1, size(found, kind=I4P)
      call blocks(owner(r))%add_record(found(r))
   enddo
   endsubroutine detect_topology_2d

   subroutine detect_core(vx, dm, tolerance, translations, found, owner)
   !< Detection on work copies: singular faces, then abutting ranges, then periodic ranges per translation (+ and -).
   type(vx_t),                       intent(in)           :: vx(:)             !< Blocks.
   integer(I4P),                     intent(in)           :: dm                !< Dimension.
   real(R8P),                        intent(in), optional :: tolerance         !< Coincidence tolerance.
   real(R8P),                        intent(in), optional :: translations(:,:) !< Translations [3, nt].
   type(topo_record_t), allocatable, intent(out)          :: found(:)          !< Records.
   integer(I4P),        allocatable, intent(out)          :: owner(:)          !< Their blocks.
   type(bv_t)                                             :: bv                !< Boundary vertices.
   type(topo_record_t)                                    :: rec               !< A record.
   logical,             allocatable                       :: covered(:,:,:)    !< Covered vertices of a face.
   real(R8P)                                              :: tol, shift(3)     !< Tolerance, translation.
   integer(I4P)                                           :: b, f, a, pass, nt !< Counters.

   tol = 0._R8P ; if (present(tolerance)) tol = max(tolerance, 0._R8P)
   nt = 0_I4P ; if (present(translations)) nt = size(translations, 2, I4P)
   allocate(found(0), owner(0))
   call boundary_vertices(vx, dm, bv)
   do b = 1, size(vx, kind=I4P)
      do f = 1, 2_I4P * dm
         call face_cover(vx(b)%n, f, covered)
         ! singular: every vertex line along an in-face direction collapses
         do a = 1, dm
            if (a == (f + 1_I4P) / 2_I4P) cycle
            if (face_collapses(vx(b), f, a, tol)) then
               rec = topo_record_t()
               rec%kind = TOPO_SINGULAR ; rec%face = f ; rec%collapsed = a
               call face_range(vx(b)%n, f, rec%lo, rec%hi)
               found = [found, rec] ; owner = [owner, b]
               covered = .true.
            endif
         enddo
         do pass = 0, 2_I4P * nt
            shift = 0._R8P
            if (pass > 0) shift = translations(:, (pass + 1_I4P) / 2_I4P) * merge(1._R8P, -1._R8P, mod(pass, 2_I4P) == 1)
            call face_ranges(vx, dm, bv, tol, b, f, shift, pass > 0, covered, found, owner)
         enddo
      enddo
   enddo
   endsubroutine detect_core

   subroutine face_ranges(vx, dm, bv, tol, b, f, shift, periodic, covered, found, owner)
   !< Abutting (periodic) ranges of face f of block b: seeds in index order, each grown to its largest rectangle.
   type(vx_t),                       intent(in)    :: vx(:)          !< Blocks.
   integer(I4P),                     intent(in)    :: dm, b, f       !< Dimension, block, face.
   type(bv_t),                       intent(in)    :: bv             !< Boundary vertices.
   real(R8P),                        intent(in)    :: tol, shift(3)  !< Tolerance, translation.
   logical,                          intent(in)    :: periodic       !< Periodic pass.
   logical,             allocatable, intent(inout) :: covered(:,:,:) !< Covered vertices (face bounds).
   type(topo_record_t), allocatable, intent(inout) :: found(:)       !< Records.
   integer(I4P),        allocatable, intent(inout) :: owner(:)       !< Their blocks.
   type(topo_record_t)                             :: rec            !< Record.
   integer(I4P)                                    :: lo(3), hi(3), p(3), i, j, k !< Work.
   logical                                         :: ok             !< Seed accepted.

   call face_range(vx(b)%n, f, lo, hi)
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            if (covered(i, j, k)) cycle
            p = [i, j, k]
            call seed_record(vx, dm, bv, tol, b, f, p, shift, rec, ok)
            if (.not. ok) cycle
            call grow_record(vx, dm, tol, b, f, shift, rec)
            if (count(rec%hi(1:dm) > rec%lo(1:dm)) < dm - 1_I4P) cycle ! a range without a cell face
            covered(rec%lo(1):rec%hi(1), rec%lo(2):rec%hi(2), rec%lo(3):rec%hi(3)) = .true.
            rec%kind = merge(TOPO_PERIODIC, TOPO_ABUT, periodic) ; rec%translation = shift
            found = [found, rec] ; owner = [owner, b]
         enddo
      enddo
   enddo
   endsubroutine face_ranges

   subroutine seed_record(vx, dm, bv, tol, b, f, p, shift, rec, ok)
   !< Transform seeded at vertex p of face f: unique image of p and of one face neighbour per in-face direction.
   type(vx_t),          intent(in)  :: vx(:)          !< Blocks.
   integer(I4P),        intent(in)  :: dm, b, f, p(3) !< Dimension, block, face, vertex.
   type(bv_t),          intent(in)  :: bv             !< Boundary vertices.
   real(R8P),           intent(in)  :: tol, shift(3)  !< Tolerance, translation.
   type(topo_record_t), intent(out) :: rec            !< Record (transform, begin, donor).
   logical,             intent(out) :: ok             !< Seed accepted.
   integer(I4P)                     :: t(4), tn(4), nd, a, s, c, dd(3), used(3), nb, side_a, side_b, e !< Work.
   logical                          :: uniq           !< Unique image.

   ok = .false. ; rec = topo_record_t()
   call unique_image(vx, bv, tol, b, p, shift, t, uniq)
   if (.not. uniq) return
   nd = (f + 1_I4P) / 2_I4P
   used = 0_I4P ; rec%transform = [1, 2, 3]
   do a = 1, dm
      if (a == nd) cycle
      dd = 0_I4P
      do s = 1, -1, -2
         e = p(a) + s
         if (e < 1 .or. e > vx(b)%n(a)) cycle
         dd(a) = s
         call unique_image(vx, bv, tol, b, p + dd, shift, tn, uniq)
         if (.not. uniq .or. tn(1) /= t(1)) cycle
         if (count(tn(2:4) /= t(2:4)) /= 1 .or. sum(abs(tn(2:4) - t(2:4))) /= 1) cycle
         do c = 1, 3
            if (tn(1 + c) /= t(1 + c)) exit
         enddo
         if (any(used == c)) cycle
         used(a) = c ; rec%transform(a) = c * s * (tn(1 + c) - t(1 + c))
         exit
      enddo
      if (used(a) == 0_I4P) return
   enddo
   ! donor normal: the direction no in-face direction maps to
   nb = 0_I4P
   do c = 1, dm
      if (all(used /= c)) nb = c
   enddo
   if (nb == 0_I4P) return
   if (t(1 + nb) == 1_I4P) then
      side_b = 1_I4P
   elseif (t(1 + nb) == vx(t(1))%n(nb)) then
      side_b = 2_I4P
   else
      return
   endif
   side_a = 2_I4P - mod(f, 2_I4P)
   rec%transform(nd) = nb * merge(1_I4P, -1_I4P, side_a /= side_b)
   if (dm == 2_I4P) rec%transform(3) = 3_I4P
   rec%face = f ; rec%donor_block = t(1) ; rec%donor_face = 2_I4P * (nb - 1_I4P) + side_b
   rec%begin = p ; rec%donor_begin = t(2:4) ; rec%lo = p ; rec%hi = p
   ok = .true.
   endsubroutine seed_record

   subroutine grow_record(vx, dm, tol, b, f, shift, rec)
   !< Grow the record's range, one line at a time on each side, while every vertex of the new line maps onto a
   !< coincident vertex (inside the donor, not onto itself).
   type(vx_t),          intent(in)    :: vx(:)         !< Blocks.
   integer(I4P),        intent(in)    :: dm, b, f      !< Dimension, block, face.
   real(R8P),           intent(in)    :: tol, shift(3) !< Tolerance, translation.
   type(topo_record_t), intent(inout) :: rec           !< Record.
   integer(I4P)                       :: a, s, nd, lo(3), hi(3) !< Work.
   logical                            :: grown         !< A side moved.

   nd = (f + 1_I4P) / 2_I4P
   grown = .true.
   do while (grown)
      grown = .false.
      do a = 1, dm
         if (a == nd) cycle
         do s = -1, 1, 2
            lo = rec%lo ; hi = rec%hi
            if (s < 0) then
               if (lo(a) == 1_I4P) cycle
               lo(a) = lo(a) - 1_I4P ; hi(a) = lo(a)
            else
               if (hi(a) == vx(b)%n(a)) cycle
               hi(a) = hi(a) + 1_I4P ; lo(a) = hi(a)
            endif
            if (line_maps(vx, tol, b, shift, rec, lo, hi)) then
               if (s < 0) then
                  rec%lo(a) = lo(a)
               else
                  rec%hi(a) = hi(a)
               endif
               grown = .true.
            endif
         enddo
      enddo
   enddo
   endsubroutine grow_record

   logical function line_maps(vx, tol, b, shift, rec, lo, hi) result(ok)
   !< Every vertex of the range [lo, hi] maps onto a coincident donor vertex, other than itself.
   type(vx_t),          intent(in) :: vx(:)          !< Blocks.
   real(R8P),           intent(in) :: tol, shift(3)  !< Tolerance, translation.
   integer(I4P),        intent(in) :: b, lo(3), hi(3) !< Block, range.
   type(topo_record_t), intent(in) :: rec            !< Record.
   integer(I4P)                    :: i, j, k, q(3)  !< Work.

   ok = .false.
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            q = record_map(rec, [i, j, k])
            if (any(q < 1) .or. any(q > vx(rec%donor_block)%n)) return
            if (rec%donor_block == b .and. all(q == [i, j, k])) return
            if (.not. coincide(vx(b)%c(:, i, j, k) + shift, vx(rec%donor_block)%c(:, q(1), q(2), q(3)), tol)) return
         enddo
      enddo
   enddo
   ok = .true.
   endfunction line_maps

   subroutine unique_image(vx, bv, tol, b, p, shift, t, uniq)
   !< The boundary vertices coinciding with vertex p of block b moved by shift, other than p itself: unique or not.
   type(vx_t),   intent(in)  :: vx(:)          !< Blocks.
   type(bv_t),   intent(in)  :: bv             !< Boundary vertices.
   real(R8P),    intent(in)  :: tol, shift(3)  !< Tolerance, translation.
   integer(I4P), intent(in)  :: b, p(3)        !< Block, vertex.
   integer(I4P), intent(out) :: t(4)           !< The image [block, i, j, k].
   logical,      intent(out) :: uniq           !< Exactly one image.
   real(R8P)                 :: c(3)           !< Moved point.
   integer(I4P)              :: e, first, n    !< Work.

   t = 0_I4P ; n = 0_I4P
   c = vx(b)%c(:, p(1), p(2), p(3)) + shift
   first = lower_bound(bv%xs, c(1) - tol)
   do e = first, bv%m
      if (bv%xs(e) > c(1) + tol) exit
      if (bv%id(1, e) == b .and. all(bv%id(2:4, e) == p)) cycle
      if (.not. coincide(c, vx(bv%id(1, e))%c(:, bv%id(2, e), bv%id(3, e), bv%id(4, e)), tol)) cycle
      if (n > 0 .and. all(bv%id(:, e) == t)) cycle
      n = n + 1_I4P ; t = bv%id(:, e)
      if (n > 1) exit
   enddo
   uniq = n == 1_I4P
   endsubroutine unique_image

   logical function face_collapses(v, f, a, tol) result(yes)
   !< Every vertex line of face f along direction a has coincident vertices.
   type(vx_t),   intent(in) :: v          !< Block.
   integer(I4P), intent(in) :: f, a       !< Face, direction.
   real(R8P),    intent(in) :: tol        !< Tolerance.
   integer(I4P)             :: lo(3), hi(3), i, j, k, q(3) !< Work.

   yes = .false.
   if (v%n(a) < 2) return
   call face_range(v%n, f, lo, hi)
   hi(a) = hi(a) - 1_I4P
   do k = lo(3), hi(3)
      do j = lo(2), hi(2)
         do i = lo(1), hi(1)
            q = [i, j, k] ; q(a) = q(a) + 1_I4P
            if (.not. coincide(v%c(:, i, j, k), v%c(:, q(1), q(2), q(3)), tol)) return
         enddo
      enddo
   enddo
   yes = .true.
   endfunction face_collapses

   pure subroutine face_range(n, f, lo, hi)
   !< Vertex range of face f.
   integer(I4P), intent(in)  :: n(3), f    !< Vertices per direction, face.
   integer(I4P), intent(out) :: lo(3), hi(3) !< Range.
   integer(I4P)              :: nd         !< Normal direction.

   nd = (f + 1_I4P) / 2_I4P
   lo = 1_I4P ; hi = n
   if (mod(f, 2_I4P) == 1_I4P) then
      hi(nd) = 1_I4P
   else
      lo(nd) = n(nd)
   endif
   endsubroutine face_range

   subroutine face_cover(n, f, covered)
   !< Covered flags of face f, all false.
   integer(I4P),         intent(in)  :: n(3), f      !< Vertices per direction, face.
   logical, allocatable, intent(out) :: covered(:,:,:) !< Flags over the whole block index box (face slice used).
   integer(I4P)                      :: lo(3), hi(3) !< Range.

   call face_range(n, f, lo, hi)
   allocate(covered(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3))) ; covered = .false.
   endsubroutine face_cover

   pure logical function coincide(p, q, tol) result(yes)
   !< Coincidence in the max norm (exact when tol = 0).
   real(R8P), intent(in) :: p(3), q(3), tol !< Points, tolerance.

   if (tol > 0._R8P) then
      yes = maxval(abs(p - q)) <= tol
   else
      yes = all(p == q)
   endif
   endfunction coincide

   subroutine boundary_vertices(vx, dm, bv)
   !< Every boundary vertex once (block, i, j, k), sorted by x (ties by y, z, then identity).
   type(vx_t),   intent(in)  :: vx(:)  !< Blocks.
   integer(I4P), intent(in)  :: dm     !< Dimension.
   type(bv_t),   intent(out) :: bv     !< Sorted boundary vertices.
   integer(I4P), allocatable :: ids(:,:), perm(:), tmp(:) !< Work.
   real(R8P),    allocatable :: key(:,:) !< Sort keys.
   integer(I4P)              :: b, i, j, k, m, nb !< Counters.
   logical                   :: on   !< On the boundary.

   nb = 0_I4P
   do b = 1, size(vx, kind=I4P)
      nb = nb + product(vx(b)%n) - product(max(vx(b)%n(1:dm) - 2_I4P, 0_I4P))
   enddo
   allocate(ids(4, nb), key(3, nb)) ; m = 0_I4P
   do b = 1, size(vx, kind=I4P)
      do k = 1, vx(b)%n(3)
         do j = 1, vx(b)%n(2)
            do i = 1, vx(b)%n(1)
               on = i == 1 .or. i == vx(b)%n(1) .or. j == 1 .or. j == vx(b)%n(2)
               if (dm == 3_I4P) on = on .or. k == 1 .or. k == vx(b)%n(3)
               if (.not. on) cycle
               m = m + 1_I4P ; ids(:, m) = [b, i, j, k] ; key(:, m) = vx(b)%c(:, i, j, k)
            enddo
         enddo
      enddo
   enddo
   allocate(perm(m), tmp(m)) ; perm = [(i, i = 1, m)]
   call merge_sort(key, ids, perm, tmp, 1_I4P, m)
   bv%m = m ; bv%id = ids(:, perm) ; bv%xs = key(1, perm)
   endsubroutine boundary_vertices

   recursive subroutine merge_sort(key, ids, perm, tmp, lo, hi)
   !< Stable merge sort of perm(lo:hi) by (key, ids).
   real(R8P),    intent(in)    :: key(:,:)  !< Keys [3, m].
   integer(I4P), intent(in)    :: ids(:,:)  !< Identities [4, m].
   integer(I4P), intent(inout) :: perm(:), tmp(:) !< Permutation, work.
   integer(I4P), intent(in)    :: lo, hi    !< Range.
   integer(I4P)                :: mid, a, b, o !< Work.

   if (hi <= lo) return
   mid = (lo + hi) / 2_I4P
   call merge_sort(key, ids, perm, tmp, lo, mid)
   call merge_sort(key, ids, perm, tmp, mid + 1_I4P, hi)
   a = lo ; b = mid + 1_I4P ; o = lo
   do while (a <= mid .and. b <= hi)
      if (precedes(key(:, perm(b)), ids(:, perm(b)), key(:, perm(a)), ids(:, perm(a)))) then
         tmp(o) = perm(b) ; b = b + 1_I4P
      else
         tmp(o) = perm(a) ; a = a + 1_I4P
      endif
      o = o + 1_I4P
   enddo
   do while (a <= mid)
      tmp(o) = perm(a) ; a = a + 1_I4P ; o = o + 1_I4P
   enddo
   do while (b <= hi)
      tmp(o) = perm(b) ; b = b + 1_I4P ; o = o + 1_I4P
   enddo
   perm(lo:hi) = tmp(lo:hi)
   endsubroutine merge_sort

   pure logical function precedes(ka, ia, kb, ib) result(yes)
   !< Lexicographic order of (key, identity).
   real(R8P),    intent(in) :: ka(3), kb(3) !< Keys.
   integer(I4P), intent(in) :: ia(4), ib(4) !< Identities.
   integer(I4P)             :: d            !< Component.

   do d = 1, 3
      if (ka(d) /= kb(d)) then
         yes = ka(d) < kb(d) ; return
      endif
   enddo
   do d = 1, 4
      if (ia(d) /= ib(d)) then
         yes = ia(d) < ib(d) ; return
      endif
   enddo
   yes = .false.
   endfunction precedes

   pure integer(I4P) function lower_bound(xs, v) result(e)
   !< First entry with xs(e) >= v (size + 1 if none).
   real(R8P), intent(in) :: xs(:) !< Sorted values.
   real(R8P), intent(in) :: v     !< Value.
   integer(I4P)          :: lo, hi, mid !< Bisection.

   lo = 1_I4P ; hi = size(xs, kind=I4P) + 1_I4P
   do while (lo < hi)
      mid = (lo + hi) / 2_I4P
      if (xs(mid) < v) then
         lo = mid + 1_I4P
      else
         hi = mid
      endif
   enddo
   e = lo
   endfunction lower_bound
endmodule fossil_curvilinear_block
