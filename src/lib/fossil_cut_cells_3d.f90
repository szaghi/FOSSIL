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
use fossil_curve_object,       only : curve_object, cut_cells_2d_t
use fossil_curve_slice,        only : slice_surface
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
   integer(I4P)         :: nx, ny, nz, i, j, k, st !< Work.
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
   endsubroutine cut_cells_3d
endmodule fossil_cut_cells_3d
