!< FOSSIL, generate Immersed Boundary distance function.

module fossil_block_object
!< Cartesian block class definition module.

use fossil, only : surface_stl_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use vtk_fortran, only : vtk_file

implicit none
private
public :: block_object

type block_object
   type(vector_R8P)              :: bmin, bmax                      !< Bounding box extents.
   integer(I4P)                  :: ni=0, nj=0, nk=0                !< Grid dimensions where distance is computed exactly.
   integer(I4P)                  :: gi(2)=0, gj(2)=0, gk(2)=0       !< Frame around the grid where distance is computed exactly.
   integer(I4P)                  :: ei(2)=0, ej(2)=0, ek(2)=0       !< Extrapolation of the grid where distance is extraplated.
   real(R8P)                     :: Dx=0._R8P, Dy=0._R8P, Dz=0._R8P !< Space steps.
   type(vector_R8P), allocatable :: nodes(:,:,:)                    !< Grid nodes.
   type(vector_R8P), allocatable :: centers(:,:,:)                  !< Grid centers.
   real(R8P),        allocatable :: distances(:,:,:)                !< Distance of grid centers to STL surface.
   contains
      ! public methods
      procedure, pass(self) :: initialize            !< Initialize block.
      procedure, pass(self) :: compute_cells_centers !< Compute cells centers from nodes.
      procedure, pass(self) :: compute_distances     !< Compute distances of cells centers from Immersed Boundary.
      procedure, pass(self) :: export_vtk_file       !< Export block data into VTK file format.
      procedure, pass(self) :: export_xall_files     !< Export block data into XALL files format.
      procedure, pass(self) :: extrapolate_distances !< Extrapolate the body-close distances in the body-far mesh.
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
   integer(I4P)                                   :: i, j, k        !< Counter.
   integer(I8P)                                   :: timing(0:2)    !< Tic toc timing.

   invert_sign_ = .false. ; if (present(invert_sign)) invert_sign_ = invert_sign
   call system_clock(timing(1))
   do k=1 - self%gk(1), self%nk + self%gk(2)
      do j=1 - self%gj(1), self%nj + self%gj(2)
         do i=1 - self%gi(1), self%ni + self%gi(2)
            self%distances(i, j, k) = surface_stl%distance(point=self%centers(i, j, k), is_signed=.true., &
                                                           is_signed=is_signed, sign_algorithm=trim(sign_algorithm))
         enddo
      enddo
   enddo
   call system_clock(timing(2), timing(0))
   print '(A, F8.3)', 'compute forces timing: ', real(timing(2) - timing(1)) / timing(0)
   if (invert_sign_) self%distances = -self%distances
   endsubroutine compute_distances

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

   pure subroutine initialize(self, bmin, bmax, ni, nj, nk, gi, gj, gk, ei, ej, ek)
   !< Initialize block.
   class(block_object), intent(inout)        :: self                !< Block
   type(vector_R8P),    intent(in), optional :: bmin, bmax          !< Bounding box extents.
   integer(I4P),        intent(in), optional :: ni, nj, nk          !< Grid dimensions where distance is computed exactly.
   integer(I4P),        intent(in), optional :: gi(2), gj(2), gk(2) !< Frame around the grid where distance is computed exactly.
   integer(I4P),        intent(in), optional :: ei(2), ej(2), ek(2) !< Extrapolation of the grid where distance is extraplated.
   type(block_object)                        :: clean               !< Clean block.
   integer(I4P)                              :: i, j, k             !< Counter.

   self = clean
   if (present(bmin).and.present(bmax).and.present(ni).and.present(ni).and.present(ni)) then
                       self%bmin = bmin
                       self%bmax = bmax
                       self%ni   = ni
                       self%nj   = nj
                       self%nk   = nk
      if (present(gi)) self%gi   = gi
      if (present(gj)) self%gj   = gj
      if (present(gk)) self%gk   = gk
      if (present(ei)) self%ei   = ei
      if (present(ej)) self%ej   = ej
      if (present(ek)) self%ek   = ek

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
   endif
   endsubroutine initialize

   ! operators
   ! =
   pure subroutine block_assign_block(lhs, rhs)
   !< Operator `=`.
   class(block_object), intent(inout) :: lhs !< Left hand side.
   type(block_object),  intent(in)    :: rhs !< Right hand side.

   lhs%bmin = rhs%bmin
   lhs%bmax = rhs%bmax
   lhs%ni   = rhs%ni
   lhs%nj   = rhs%nj
   lhs%nk   = rhs%nk
   lhs%gi   = rhs%gi
   lhs%gj   = rhs%gj
   lhs%gk   = rhs%gk
   lhs%ei   = rhs%ei
   lhs%ej   = rhs%ej
   lhs%ek   = rhs%ek
   lhs%Dx   = rhs%Dx
   lhs%Dy   = rhs%Dy
   lhs%Dz   = rhs%Dz
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

program fossil_ib_generator
!< FOSSIL, generate distances function from Immersed Boundary.

use fossil_block_object, only : block_object
use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use fossil_aabb_tree_object, only : aabb_tree_object
use vtk_fortran, only : vtk_file

implicit none

type(command_line_interface) :: cli                     !< Test command line interface.
type(file_stl_object)        :: file_stl                !< STL file.
type(surface_stl_object)     :: surface_stl             !< STL surface.
type(vector_R8P)             :: bmin, bmax              !< Bounding box extents.
character(999)               :: file_name_stl           !< Input STL file name.
character(999)               :: output_base_name        !< Output base name.
integer(I4P)                 :: refinement_levels       !< AABB refinement levels used.
logical                      :: is_signed               !< Signed distance or not.
character(999)               :: sign_algorithm          !< Algorithm used for "point in polyhedron" test.
logical                      :: unsigned                !< Compute unsigned distance.
logical                      :: save_aabb_tree_geometry !< Sentinel to save AABB geometry.
logical                      :: save_aabb_tree_stl      !< Sentinel to save AABB stl.
integer(I4P)                 :: ni, nj, nk              !< Grid dimensions, close to STL: where distance is computed exactly.
integer(I4P)                 :: gi(2), gj(2), gk(2)     !< Frame dimensions around the grid: where distance is computed exactly.
integer(I4P)                 :: ei(2), ej(2), ek(2)     !< Dimensions of extrapolation of the grid: where distance is extraplated.
type(block_object)           :: cblock                  !< Cartesian block.

! parse command line input e load STL file
call cli_parse
call file_stl%load_from_file(facet=surface_stl%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
call surface_stl%analize(aabb_refinement_levels=refinement_levels)
call surface_stl%sanitize
call surface_stl%analize(aabb_refinement_levels=refinement_levels)
print '(A)', surface_stl%statistics()
if (save_aabb_tree_geometry) call surface_stl%aabb%save_geometry_tecplot_ascii(file_name=trim(output_base_name)//'_aabb_tree.dat')
if (save_aabb_tree_stl) call file_stl%save_aabb_into_file(surface=surface_stl, base_file_name=trim(output_base_name), &
                                                          is_ascii=.false.)

if (.not.cli%is_passed(switch='--bmin')) then
   bmin = surface_stl%bmin
endif
if (.not.cli%is_passed(switch='--bmax')) then
   bmax = surface_stl%bmax
endif

call cblock%initialize(bmin=bmin, bmax=bmax, ni=ni, nj=nj, nk=nk, gi=gi, gj=gj, gk=gk, ei=ei, ej=ej, ek=ek)
call cblock%compute_distances(surface_stl=surface_stl, sign_algorithm=sign_algorithm, invert_sign=.true.)
call cblock%export_vtk_file(file_name=trim(output_base_name)//'.vts')
call cblock%export_xall_files(basename=trim(output_base_name))

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  real(R8P)    :: bbmin(3), bbmax(3) !< Bounding box extents, local variable.
  integer(I4P) :: error              !< Error trapping flag.

  call cli%init(progname='fossil_ib_generator',                              &
                authors='S. Zaghi',                                          &
                help='Usage: ',                                              &
                examples=["fossil_ib_generator --stl src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='--stl',                      &
               help='STL (input) file name',        &
               required=.false.,                    &
               def='src/tests/naca0012-binary.stl', &
               act='store')

  call cli%add(switch='--out',          &
               help='output base name', &
               required=.false.,        &
               def='ib',                &
               act='store')

  call cli%add(switch='--bmin',             &
               help='bounding box minimum', &
               required=.false.,            &
               nargs='3',                   &
               def='0.0 0.0 0.0',           &
               act='store')

  call cli%add(switch='--bmax',             &
               help='bounding box maximum', &
               required=.false.,            &
               nargs='3',                   &
               def='0.0 0.0 0.0',           &
               act='store')

  call cli%add(switch='--ni',                      &
               help='cells number in i direction', &
               required=.false.,                   &
               def='32',                           &
               act='store')

  call cli%add(switch='--nj',                      &
               help='cells number in j direction', &
               required=.false.,                   &
               def='32',                           &
               act='store')

  call cli%add(switch='--nk',                      &
               help='cells number in k direction', &
               required=.false.,                   &
               def='32',                           &
               act='store')

  call cli%add(switch='--gi',                            &
               help='ghost cells number in i direction', &
               required=.false.,                         &
               nargs='2',                                &
               def='2 2',                                &
               act='store')

  call cli%add(switch='--gj',                            &
               help='ghost cells number in j direction', &
               required=.false.,                         &
               nargs='2',                                &
               def='2 2',                                &
               act='store')

  call cli%add(switch='--gk',                            &
               help='ghost cells number in k direction', &
               required=.false.,                         &
               nargs='2',                                &
               def='2 2',                                &
               act='store')

  call cli%add(switch='--ei',                                                  &
               help='extrapolation cells number in i directions (front/back)', &
               required=.false.,                                               &
               nargs='2',                                                      &
               def='0 0',                                                      &
               act='store')

  call cli%add(switch='--ej',                                                  &
               help='extrapolation cells number in j directions (front/back)', &
               required=.false.,                                               &
               nargs='2',                                                      &
               def='0 0',                                                      &
               act='store')

  call cli%add(switch='--ek',                                                  &
               help='extrapolation cells number in k directions (front/back)', &
               required=.false.,                                               &
               nargs='2',                                                      &
               def='0 0',                                                      &
               act='store')

  call cli%add(switch='--ref_levels',         &
               help='AABB refinement levels', &
               required=.false.,              &
               def='2',                       &
               act='store')

  call cli%add(switch='--is_signed',          &
               help='signed distance or not', &
               required=.false.,              &
               def='.true.',                  &
               act='store_true')

  call cli%add(switch='--sign_algorithm',                         &
               help='algorithm used to compute sign of distance', &
               required=.false.,                                  &
               def='ray_intersections',                           &
               act='store')

  call cli%add(switch='--save_aabb_tree_geometry', &
               help='save AABB tree geometry',     &
               required=.false.,                   &
               def='.true.',                       &
               act='store_true')

  call cli%add(switch='--save_aabb_tree_stl', &
               help='save AABB tree STL',     &
               required=.false.,              &
               def='.true.',                  &
               act='store_true')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',                     val=file_name_stl,           error=error) ; if (error/=0) stop
  call cli%get(switch='--out',                     val=output_base_name,        error=error) ; if (error/=0) stop
  call cli%get(switch='--bmin',                    val=bbmin,                   error=error) ; if (error/=0) stop
  call cli%get(switch='--bmax',                    val=bbmax,                   error=error) ; if (error/=0) stop
  call cli%get(switch='--ni',                      val=ni,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--nj',                      val=nj,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--nk',                      val=nk,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gi',                      val=gi,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gj',                      val=gj,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--gk',                      val=gk,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--ei',                      val=ei,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--ej',                      val=ej,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--ek',                      val=ek,                      error=error) ; if (error/=0) stop
  call cli%get(switch='--ref_levels',              val=refinement_levels,       error=error) ; if (error/=0) stop
  call cli%get(switch='--is_signed',               val=is_signed,               error=error) ; if (error/=0) stop
  call cli%get(switch='--sign_algorithm',          val=sign_algorithm,          error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_geometry', val=save_aabb_tree_geometry, error=error) ; if (error/=0) stop
  call cli%get(switch='--save_aabb_tree_stl',      val=save_aabb_tree_stl,      error=error) ; if (error/=0) stop
  bmin%x = bbmin(1) ; bmin%y = bbmin(2) ; bmin%z = bbmin(3)
  bmax%x = bbmax(1) ; bmax%y = bbmax(2) ; bmax%z = bbmax(3)
  endsubroutine cli_parse
endprogram fossil_ib_generator
