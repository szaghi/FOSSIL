!< FOSSIL cartesian block class definition.

module fossil_block_object
!< FOSSIL cartesian block class definition.

use fossil_block_aabb_object, only : aabb_object
use fossil, only : surface_stl_object
use penf, only : I4P, I8P, R8P, FR8P, MaxR8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use vtk_fortran, only : vtk_file

implicit none
private
public :: block_object

type :: block_object
   !< Cartesian block class definition.
   type(vector_R8P)               :: bmin, bmax                      !< Bounding box extents.
   integer(I4P)                   :: ni=0, nj=0, nk=0                !< Grid dimensions where distance is computed exactly.
   integer(I4P)                   :: gi(2)=0, gj(2)=0, gk(2)=0       !< Frame around the grid where distance is computed exactly.
   integer(I4P)                   :: ei(2)=0, ej(2)=0, ek(2)=0       !< Extrapolation of the grid where distance is extraplated.
   real(R8P)                      :: Dx=0._R8P, Dy=0._R8P, Dz=0._R8P !< Space steps.
   integer(I4P)                   :: refinement_levels=0             !< Total number of refinement levels used.
   integer(I4P)                   :: Naabb=0                         !< Number of AABB in each directions, 2**refinement_levels.
   type(aabb_object), allocatable :: aabb(:,:,:)                     !< AABB refinement blocks.
   type(vector_R8P),  allocatable :: nodes(:,:,:)                    !< Grid nodes.
   type(vector_R8P),  allocatable :: centers(:,:,:)                  !< Grid centers.
   real(R8P),         allocatable :: distances(:,:,:)                !< Distance of grid centers to STL surface.
   contains
      ! public methods
      procedure, pass(self) :: compute_cells_centers     !< Compute cells centers from nodes.
      procedure, pass(self) :: compute_distances         !< Compute distances of cells centers from Immersed Boundary.
      procedure, pass(self) :: export_aabb_tecplot_ascii !< Export AABB boxes geometry into Tecplot ascii file.
      procedure, pass(self) :: export_vtk_file           !< Export block data into VTK file format.
      procedure, pass(self) :: export_xall_files         !< Export block data into XALL files format.
      procedure, pass(self) :: extrapolate_distances     !< Extrapolate the body-close distances in the body-far mesh.
      procedure, pass(self) :: get_closest_cells_indexes !< Get the closest cells indexes in the mesh given a point.
      procedure, pass(self) :: initialize                !< Initialize block.
      ! operators
      generic :: assignment(=) => block_assign_block      !< Overload `=`.
      procedure, pass(lhs), private :: block_assign_block !< Operator `=`.
   endtype block_object

contains
   ! public methods
   pure subroutine compute_cells_centers(self)
   !< Compute cells centers from nodes.
   class(block_object), intent(inout) :: self    !< Block.
   integer(I4P)                       :: i, j, k !< Counter.

   do k=1 - self%gk(1), self%nk + self%gk(2)
      do j=1 - self%gj(1), self%nj + self%gj(2)
         do i=1 - self%gi(1), self%ni + self%gi(2)
            self%centers(i, j, k) = (self%nodes(i-1,j-1,k-1) + self%nodes(i  ,j-1,k-1) + &
                                     self%nodes(i-1,j  ,k-1) + self%nodes(i  ,j  ,k-1) + &
                                     self%nodes(i-1,j-1,k  ) + self%nodes(i  ,j-1,k  ) + &
                                     self%nodes(i-1,j  ,k  ) + self%nodes(i  ,j  ,k  ))/8._R8P
         enddo
      enddo
   enddo
   endsubroutine compute_cells_centers

   subroutine compute_distances(self, surface_stl, is_signed, sign_algorithm, invert_sign)
   !< Compute distances of cells centers from Immersed Boundary.
   class(block_object),      intent(inout)        :: self           !< Block.
   type(surface_stl_object), intent(in)           :: surface_stl    !< STL surface.
   logical,                  intent(in)           :: is_signed      !< Signed distance or not.
   character(*),             intent(in)           :: sign_algorithm !< Algorithm used for "point in polyhedron" test.
   logical,                  intent(in), optional :: invert_sign    !< Invert sign of distances.
   logical                                        :: invert_sign_   !< Invert sign of distances, local variable.
   integer(I4P)                                   :: i, j, k, c, f  !< Counter.
   integer(I4P), parameter                        :: Nc=10          !< Number of closest cells to facets centroid.
   integer(I4P)                                   :: cindexes(3,Nc) !< Indexes (i,j,k) of the Nc closest cells.
   real(R8P)                                      :: distance       !< Temporary distance buffer.
   integer(I8P)                                   :: timing(0:2)    !< Tic toc timing.

   invert_sign_ = .false. ; if (present(invert_sign)) invert_sign_ = invert_sign
   self%distances = MaxR8P
   call system_clock(timing(1))
   do f=1, surface_stl%facets_number
      call self%get_closest_cells_indexes(point=surface_stl%facet(f)%centroid, cindexes=cindexes)
      do c=1, Nc
         i = cindexes(1, c)
         j = cindexes(2, c)
         k = cindexes(3, c)
         call surface_stl%facet(f)%compute_distance(point=self%centers(i, j, k), distance=distance)
         if (distance < self%distances(i,j,k)) self%distances(i,j,k) = distance
         if (surface_stl%is_point_inside(point=self%centers(i, j, k), sign_algorithm=trim(sign_algorithm))) &
            self%distances(i, j, k) = - self%distances(i, j, k)
      enddo
   enddo
   call system_clock(timing(2), timing(0))
   print '(A, F8.3)', 'compute forces timing: ', real(timing(2) - timing(1)) / timing(0)
   if (invert_sign_) self%distances = -self%distances
   endsubroutine compute_distances

   subroutine export_aabb_tecplot_ascii(self, file_name)
   !< Export AABB boxes geometry into Tecplot ascii file.
   class(block_object), intent(in) :: self       !< Block.
   character(*),        intent(in) :: file_name  !< File name.
   integer(I4P)                    :: file_unit  !< File unit.
   type(vector_R8P)                :: vertex(8)  !< AABB vertices.
   integer(I4P)                    :: i, j, k, v !< Counter.

   if (self%Naabb > 0) then
      open(newunit=file_unit, file=trim(adjustl(file_name)))
      write(file_unit, '(A)') 'VARIABLES=x y z'
      do k=1, self%Naabb
         do j=1, self%Naabb
            do i=1, self%Naabb
               write(file_unit, '(A)') 'ZONE T="aabb-'//trim(str(i, .true.))//'_'// &
                                                        trim(str(j, .true.))//'_'// &
                                                        trim(str(k, .true.))//'_'//'", I=2, J=2, K=2'
               vertex = self%aabb(i,j,k)%vertex()
               do v=1, 8
                  write(file_unit, '(3('//FR8P//',1X))') vertex(v)%x, vertex(v)%y, vertex(v)%z
               enddo
            enddo
         enddo
      enddo
      close(file_unit)
   endif
   endsubroutine export_aabb_tecplot_ascii

   subroutine export_vtk_file(self, file_name)
   !< Export block data into VTK file format.
   class(block_object), intent(in) :: self           !< Block.
   character(*),        intent(in) :: file_name      !< File name.
   integer(I4P)                    :: ni, nj, nk, nn !< Grid dimensions.
   type(vtk_file)                  :: vtk            !< VTK file handler.
   integer(I4P)                    :: error          !< Status error.

   ni = size(self%distances, dim=1)
   nj = size(self%distances, dim=2)
   nk = size(self%distances, dim=3)
   nn = (ni+1) * (nj+1) * (nk+1)
   error = vtk%initialize(format='raw', filename=adjustl(trim(file_name)), mesh_topology='StructuredGrid', &
                          nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
   error = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
   error = vtk%xml_writer%write_geo(n=nn, x=self%nodes%x, y=self%nodes%y, z=self%nodes%z)
   error = vtk%xml_writer%write_dataarray(location='cell', action='open')
   error = vtk%xml_writer%write_dataarray(data_name='distance', x=self%distances, one_component=.true.)
   error = vtk%xml_writer%write_dataarray(location='cell', action='close')
   error = vtk%xml_writer%write_piece()
   error = vtk%finalize()
   endsubroutine export_vtk_file

   subroutine export_xall_files(self, basename)
   !< Export block data into XALL files format.
   class(block_object), intent(in) :: self             !< Block.
   character(*),        intent(in) :: basename         !< Base files name.
   integer(I4P)                    :: ni, nj, nk       !< Grid dimensions.
   integer(I4P)                    :: i(2), j(2), k(2) !< Grid bounds.
   integer(I4P)                    :: funit            !< File unit.

   ! number of cells subtracted ghost cells
   ni = size(self%distances, dim=1) - 4
   nj = size(self%distances, dim=2) - 4
   nk = size(self%distances, dim=3) - 4
   ! grid bounds
   i(1) = lbound(self%nodes, dim=1)
   i(2) = ubound(self%nodes, dim=1)
   j(1) = lbound(self%nodes, dim=2)
   j(2) = ubound(self%nodes, dim=2)
   k(1) = lbound(self%nodes, dim=3)
   k(2) = ubound(self%nodes, dim=3)
   ! grid
   open(newunit=funit, file=adjustl(trim(basename))//'.grd', form='UNFORMATTED', action='WRITE')
   write(funit) 1 ! number of blocks
   write(funit) ni, nj, nk
   write(funit) self%nodes(i(1)+2:i(2)-2,j(1)+2:j(2)-2,k(1)+2:k(2)-2)%x
   write(funit) self%nodes(i(1)+2:i(2)-2,j(1)+2:j(2)-2,k(1)+2:k(2)-2)%y
   write(funit) self%nodes(i(1)+2:i(2)-2,j(1)+2:j(2)-2,k(1)+2:k(2)-2)%z
   close(funit)

   ! distance
   open(newunit=funit, file=adjustl(trim(basename))//'.ib', form='UNFORMATTED', action='WRITE')
   write(funit) 0._R8P ! time
   write(funit) ni, nj, nk
   write(funit) self%distances
   close(funit)
   endsubroutine export_xall_files

   subroutine extrapolate_distances(self)
   !< Extrapolate the body-close distances in the body-far mesh.
   class(block_object), intent(inout) :: self                 !< Block.
   type(vector_R8P), allocatable      :: nodes_ext(:,:,:)     !< Grid nodes extrapolated.
   real(R8P),        allocatable      :: distances_ext(:,:,:) !< Distance of grid centers extrapolated.
   real(R8P)                          :: Dxyz_ext             !< Space steps for extending domain.
   integer(I4P)                       :: i, j, k              !< Counter.

   associate(ni=>self%ni, nj=>self%nj, nk=>self%nk, &
             ei=>self%ei, ej=>self%ej, ek=>self%ek, &
             Dx=>self%Dx, Dy=>self%Dy, Dz=>self%Dz, &
             nodes=>self%nodes, distances=>self%distances)
      if (ei(1)>=2.or.ei(2)>=2.or.ej(1)>=2.or.ej(2)>=2.or.ek(1)>=2.or.ek(2)>=2) then
         if (ei(1)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1)+ei(1), size(nodes,     dim=2), size(nodes,     dim=3)))
            allocate(distances_ext(size(distances, dim=1)+ei(1), size(distances, dim=2), size(distances, dim=3)))
            nodes_ext(1:size(nodes, dim=1), :, :) = nodes
            Dxyz_ext = Dx
            do i=size(nodes, dim=1)+1, size(nodes, dim=1)+ei(1)
               nodes_ext(i, :, :)%x = nodes_ext(i-1, :, :)%x + Dxyz_ext
               nodes_ext(i, :, :)%y = nodes_ext(i-1, :, :)%y
               nodes_ext(i, :, :)%z = nodes_ext(i-1, :, :)%z
               if (i<size(nodes, dim=1)+ei(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(1:size(distances, dim=1), :, :) = distances
            do i=size(distances, dim=1)+1, size(distances, dim=1)+ei(1)
               distances_ext(i, :, :) = distances_ext(i-1, :, :)
            enddo
            ni = ni + ei(1)
            ei(1) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif

         if (ei(2)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1)+ei(2), size(nodes,     dim=2), size(nodes,     dim=3)))
            allocate(distances_ext(size(distances, dim=1)+ei(2), size(distances, dim=2), size(distances, dim=3)))
            nodes_ext(ei(2)+1:, :, :) = nodes
            Dxyz_ext = Dx
            do i=ei(2), 1, -1
               nodes_ext(i, :, :)%x = nodes_ext(i+1, :, :)%x - Dxyz_ext
               nodes_ext(i, :, :)%y = nodes_ext(i+1, :, :)%y
               nodes_ext(i, :, :)%z = nodes_ext(i+1, :, :)%z
               if (i>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(ei(2)+1:, :, :) = distances
            do i=ei(2), 1, -1
               distances_ext(i, :, :) = distances_ext(i+1, :, :)
            enddo
            ni = ni + ei(2)
            ei(2) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif

         if (ej(1)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1), size(nodes,     dim=2)+ej(1), size(nodes,     dim=3)))
            allocate(distances_ext(size(distances, dim=1), size(distances, dim=2)+ej(1), size(distances, dim=3)))
            nodes_ext(:, 1:size(nodes, dim=2), :) = nodes
            Dxyz_ext = Dy
            do j=size(nodes, dim=2)+1, size(nodes, dim=2)+ej(1)
               nodes_ext(:, j, :)%x = nodes_ext(:, j-1, :)%x
               nodes_ext(:, j, :)%y = nodes_ext(:, j-1, :)%y + Dxyz_ext
               nodes_ext(:, j, :)%z = nodes_ext(:, j-1, :)%z
               if (j<size(nodes, dim=2)+ej(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(:, 1:size(distances, dim=2), :) = distances
            do j=size(distances, dim=2)+1, size(distances, dim=2)+ej(1)
               distances_ext(:, j, :) = distances_ext(:, j-1, :)
            enddo
            nj = nj + ej(1)
            ej(1) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif

         if (ej(2)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1), size(nodes,     dim=2)+ej(2), size(nodes,     dim=3)))
            allocate(distances_ext(size(distances, dim=1), size(distances, dim=2)+ej(2), size(distances, dim=3)))
            nodes_ext(:, ej(2)+1:, :) = nodes
            Dxyz_ext = Dy
            do j=ej(2), 1, -1
               nodes_ext(:, j, :)%x = nodes_ext(:, j+1, :)%x
               nodes_ext(:, j, :)%y = nodes_ext(:, j+1, :)%y - Dxyz_ext
               nodes_ext(:, j, :)%z = nodes_ext(:, j+1, :)%z
               if (j>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(:, ej(2)+1:, :) = distances
            do j=ej(2), 1, -1
               distances_ext(:, j, :) = distances_ext(:, j+1, :)
            enddo
            nj = nj + ej(2)
            ej(2) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif

         if (ek(1)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1), size(nodes,     dim=2), size(nodes,     dim=3)+ek(1)))
            allocate(distances_ext(size(distances, dim=1), size(distances, dim=2), size(distances, dim=3)+ek(1)))
            nodes_ext(:, :, 1:size(nodes, dim=3)) = nodes
            Dxyz_ext = Dz
            do k=size(nodes, dim=3)+1, size(nodes, dim=3)+ek(1)
               nodes_ext(:, :, k)%x = nodes_ext(:, :, k-1)%x
               nodes_ext(:, :, k)%y = nodes_ext(:, :, k-1)%y
               nodes_ext(:, :, k)%z = nodes_ext(:, :, k-1)%z + Dxyz_ext
               if (k<size(nodes, dim=3)+ek(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(:, :, 1:size(distances, dim=3)) = distances
            do k=size(distances, dim=3)+1, size(distances, dim=3)+ek(1)
               distances_ext(:, :, k) = distances_ext(:, :, k-1)
            enddo
            nk = nk + ek(1)
            ek(1) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif

         if (ek(2)>=2) then
            allocate(nodes_ext(    size(nodes,     dim=1), size(nodes,     dim=2), size(nodes,     dim=3)+ek(2)))
            allocate(distances_ext(size(distances, dim=1), size(distances, dim=2), size(distances, dim=3)+ek(2)))
            nodes_ext(:, :, ek(2)+1:) = nodes
            Dxyz_ext = Dz
            do k=ek(2), 1, -1
               nodes_ext(:, :, k)%x = nodes_ext(:, :, k+1)%x
               nodes_ext(:, :, k)%y = nodes_ext(:, :, k+1)%y
               nodes_ext(:, :, k)%z = nodes_ext(:, :, k+1)%z - Dxyz_ext
               if (k>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
            enddo
            distances_ext(:, :, ek(2)+1:) = distances
            do k=ek(2), 1, -1
               distances_ext(:, :, k) = distances_ext(:, :, k+1)
            enddo
            nk = nk + ek(2)
            ek(2) = 0
            call move_alloc(from=nodes_ext,     to=self%nodes)
            call move_alloc(from=distances_ext, to=self%distances)
         endif
        call self%compute_cells_centers
      endif
   endassociate
   endsubroutine extrapolate_distances

   subroutine get_closest_cells_indexes(self, point, cindexes)
   !< Get the closest cells indexes in the mesh given a point.
   class(block_object), intent(in)  :: self            !< Block
   type(vector_R8P),    intent(in)  :: point           !< Point coordinates.
   integer(I4P),        intent(out) :: cindexes(1:,1:) !< Indexes (i,j,k) of the Nc closest cells.
   endsubroutine get_closest_cells_indexes

   pure subroutine initialize(self, bmin, bmax, ni, nj, nk, gi, gj, gk, ei, ej, ek, refinement_levels)
   !< Initialize block.
   class(block_object), intent(inout)        :: self                !< Block
   type(vector_R8P),    intent(in), optional :: bmin, bmax          !< Bounding box extents.
   integer(I4P),        intent(in), optional :: ni, nj, nk          !< Grid dimensions where distance is computed exactly.
   integer(I4P),        intent(in), optional :: gi(2), gj(2), gk(2) !< Frame around the grid where distance is computed exactly.
   integer(I4P),        intent(in), optional :: ei(2), ej(2), ek(2) !< Extrapolation of the grid where distance is extraplated.
   integer(I4P),        intent(in), optional :: refinement_levels   !< Total number of refinement levels used.
   type(block_object)                        :: clean               !< Clean block.
   type(vector_R8P)                          :: Dbb                 !< Delta bounding box extents.
   integer(I4P)                              :: i, j, k             !< Counter.

   self = clean
   if (present(bmin).and.present(bmax).and.present(ni).and.present(ni).and.present(ni)) then
                                      self%bmin              = bmin
                                      self%bmax              = bmax
                                      self%ni                = ni
                                      self%nj                = nj
                                      self%nk                = nk
      if (present(gi))                self%gi                = gi
      if (present(gj))                self%gj                = gj
      if (present(gk))                self%gk                = gk
      if (present(ei))                self%ei                = ei
      if (present(ej))                self%ej                = ej
      if (present(ek))                self%ek                = ek
      if (present(refinement_levels)) self%refinement_levels = refinement_levels

      ! create body-close mesh where distances are computed exactly
      self%Dx = (self%bmax%x - self%bmin%x) / self%ni
      self%Dy = (self%bmax%y - self%bmin%y) / self%nj
      self%Dz = (self%bmax%z - self%bmin%z) / self%nk
      allocate(self%nodes(    0-self%gi(1):self%ni+self%gi(2), 0-self%gj(1):self%nj+self%gj(2), 0-self%gk(1):self%nk+self%gk(2)))
      allocate(self%centers(  1-self%gi(1):self%ni+self%gi(2), 1-self%gj(1):self%nj+self%gj(2), 1-self%gk(1):self%nk+self%gk(2)))
      allocate(self%distances(1-self%gi(1):self%ni+self%gi(2), 1-self%gj(1):self%nj+self%gj(2), 1-self%gk(1):self%nk+self%gk(2)))
      do k=0 - self%gk(1), self%nk + self%gk(2)
         do j=0 - self%gj(1), self%nj + self%gj(2)
            do i=0 - self%gi(1), self%ni + self%gi(2)
               self%nodes(i, j, k) = self%bmin + (i * self%Dx) * ex_R8P + (j * self%Dy) * ey_R8P + (k * self%Dz) * ez_R8P
            enddo
         enddo
      enddo
      call self%compute_cells_centers
      if (self%refinement_levels > 0) then
         self%Naabb = 2**self%refinement_levels
         Dbb = (self%bmax - self%bmin) / self%Naabb
         allocate(self%aabb(1:self%Naabb,1:self%Naabb,1:self%Naabb))
         do k=1, self%Naabb
            do j=1, self%Naabb
               do i=1, self%Naabb
                  self%aabb(i,j,k)%bmin%x = self%bmin%x + (i-1) * Dbb%x ; self%aabb(i,j,k)%bmax%x = self%bmin%x + i * Dbb%x
                  self%aabb(i,j,k)%bmin%y = self%bmin%y + (j-1) * Dbb%y ; self%aabb(i,j,k)%bmax%y = self%bmin%y + j * Dbb%y
                  self%aabb(i,j,k)%bmin%z = self%bmin%z + (k-1) * Dbb%z ; self%aabb(i,j,k)%bmax%z = self%bmin%z + k * Dbb%z
               enddo
            enddo
         enddo
      endif
   endif
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine block_assign_block(lhs, rhs)
   !< Operator `=`.
   class(block_object), intent(inout) :: lhs !< Left hand side.
   type(block_object),  intent(in)    :: rhs !< Right hand side.

   lhs%bmin              = rhs%bmin
   lhs%bmax              = rhs%bmax
   lhs%ni                = rhs%ni
   lhs%nj                = rhs%nj
   lhs%nk                = rhs%nk
   lhs%gi                = rhs%gi
   lhs%gj                = rhs%gj
   lhs%gk                = rhs%gk
   lhs%ei                = rhs%ei
   lhs%ej                = rhs%ej
   lhs%ek                = rhs%ek
   lhs%Dx                = rhs%Dx
   lhs%Dy                = rhs%Dy
   lhs%Dz                = rhs%Dz
   lhs%refinement_levels = rhs%refinement_levels
   lhs%Naabb             = rhs%Naabb
   if (allocated(rhs%aabb)) then
      lhs%aabb = rhs%aabb
   else
      if (allocated(lhs%aabb)) deallocate(lhs%aabb)
   endif
   if (allocated(rhs%nodes)) then
      lhs%nodes = rhs%nodes
   else
      if (allocated(lhs%nodes)) deallocate(lhs%nodes)
   endif
   if (allocated(rhs%centers)) then
      lhs%centers = rhs%centers
   else
      if (allocated(lhs%centers)) deallocate(lhs%centers)
   endif
   if (allocated(rhs%distances)) then
      lhs%distances = rhs%distances
   else
      if (allocated(lhs%distances)) deallocate(lhs%distances)
   endif
endsubroutine block_assign_block
endmodule fossil_block_object

