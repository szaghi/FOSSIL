!< FOSSIL, generate Immersed Boundary distance function.

program fossil_ib_generator
!< FOSSIL, generate Immersed Boundary distance function.

use flap, only : command_line_interface
use fossil, only : file_stl_object, surface_stl_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : ex_R8P, ey_R8P, ez_R8P, vector_R8P
use fossil_aabb_tree_object, only : aabb_tree_object
use vtk_fortran, only : vtk_file

implicit none

type(command_line_interface)   :: cli                  !< Test command line interface.
type(file_stl_object)          :: file_stl             !< STL file.
type(surface_stl_object)       :: surface_stl          !< STL surface.
type(vector_R8P), allocatable  :: nodes(:,:,:)         !< Grid nodes.
type(vector_R8P), allocatable  :: centers(:,:,:)       !< Grid centers.
real(R8P),        allocatable  :: distance(:,:,:)      !< Distance of grid centers to STL surface.
type(vector_R8P), allocatable  :: nodes_ext(:,:,:)     !< Grid nodes extrapolated.
real(R8P),        allocatable  :: distance_ext(:,:,:)  !< Distance of grid centers extrapolated.
type(vector_R8P)               :: bmin, bmax           !< Bounding box extents.
character(999)                 :: file_name_stl        !< Input STL file name.
character(999)                 :: output_base_name     !< Output base name.
integer(I4P)                   :: refinement_levels    !< AABB refinement levels used.
character(999)                 :: sign_algorithm       !< Algorithm used for "point in polyhedron" test.
logical                        :: unsigned             !< Compute unsigned distance.
integer(I4P)                   :: ni, nj, nk           !< Grid dimensions, close to STL: where distance is computed exactly.
integer(I4P)                   :: gi, gj(2), gk        !< Frame dimensions around the grid: where distance is computed exactly.
integer(I4P)                   :: ei(2), ej(2), ek(2)  !< Dimensions of extrapolation of the grid: where distance is extraplated.
integer(I4P)                   :: i, j, k              !< Counter.
real(R8P)                      :: Dx, Dy, Dz           !< Space steps.
real(R8P)                      :: Dxyz_ext             !< Space steps for extending domain.
type(vtk_file)                 :: vtk                  !< VTK file handler.
integer(I4P)                   :: error                !< Status error.
integer(I8P)                   :: timing(0:4)          !< Tic toc timing.

! parse command line input e load STL file
call cli_parse
call file_stl%load_from_file(facet=surface_stl%facet, file_name=trim(adjustl(file_name_stl)), guess_format=.true.)
call surface_stl%analize(aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics before sanitization'
print '(A)', file_stl%statistics()
print '(A)', surface_stl%statistics()
call surface_stl%sanitize
call surface_stl%analize(aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics after sanitization'
print '(A)', surface_stl%statistics()

if (.not.cli%is_passed(switch='--bmin')) then
   bmin = surface_stl%bmin
endif
if (.not.cli%is_passed(switch='--bmax')) then
   bmax = surface_stl%bmax
endif

! create body-close mesh
Dx = (bmax%x - bmin%x) / ni
Dy = (bmax%y - bmin%y) / nj
Dz = (bmax%z - bmin%z) / nk

allocate(nodes(   0-gi:ni+gi, 0-gj(1):nj+gj(2), 0-gk:nk+gk))
allocate(centers( 1-gi:ni+gi, 1-gj(1):nj+gj(2), 1-gk:nk+gk))
allocate(distance(1-gi:ni+gi, 1-gj(1):nj+gj(2), 1-gk:nk+gk))

do k=0 - gk, nk + gk
   do j=0 - gj(1), nj + gj(2)
      do i=0 - gi, ni + gi
         nodes(i, j, k) = bmin + (i * Dx) * ex_R8P + (j * Dy) * ey_R8P + (k * Dz) * ez_R8P
      enddo
   enddo
enddo

do k=1 - gk, nk + gk
   do j=1 - gj(1), nj + gj(2)
      do i=1 - gi, ni + gi
         centers(i, j, k) = (nodes(i-1,j-1,k-1) + nodes(i  ,j-1,k-1) + nodes(i-1,j  ,k-1) + nodes(i  ,j  ,k-1) + &
                             nodes(i-1,j-1,k  ) + nodes(i  ,j-1,k  ) + nodes(i-1,j  ,k  ) + nodes(i  ,j  ,k  ))/8._R8P
      enddo
   enddo
enddo

! compute signed distance field in the body-close mesh
surface_stl%aabb%is_initialized = .true.
print '(A)', 'compute distances AABB'
call system_clock(timing(3))
do k=1 - gk, nk + gk
   do j=1 - gj(1), nj + gj(2)
      do i=1 - gi, ni + gi
         distance(i, j, k) = surface_stl%distance(point=centers(i, j, k), is_signed=.not.unsigned, &
                                                  sign_algorithm=trim(sign_algorithm))
      enddo
   enddo
enddo
call system_clock(timing(4), timing(0))
print '(A, F8.3)', 'AABB timing: ', real(timing(4) - timing(3))/ timing(0)
distance = -distance

! extrapolation of the body-close signed distance in the body-far mesh
if (ei(1)>=2.or.ei(2)>=2.or.ej(1)>=2.or.ej(2)>=2.or.ek(1)>=2.or.ek(2)>=2) then
   if (ei(1)>=2) then
      if (allocated(nodes_ext))    deallocate(nodes_ext)
      if (allocated(distance_ext)) deallocate(distance_ext)
      allocate(nodes_ext(   size(nodes,    dim=1)+ei(1), size(nodes,    dim=2), size(nodes,    dim=3)))
      allocate(distance_ext(size(distance, dim=1)+ei(1), size(distance, dim=2), size(distance, dim=3)))
      nodes_ext(1:size(nodes, dim=1), :, :) = nodes
      Dxyz_ext = Dx
      do i=size(nodes, dim=1)+1, size(nodes, dim=1)+ei(1)
         nodes_ext(i, :, :)%x = nodes_ext(i-1, :, :)%x + Dxyz_ext
         nodes_ext(i, :, :)%y = nodes_ext(i-1, :, :)%y
         nodes_ext(i, :, :)%z = nodes_ext(i-1, :, :)%z
         if (i<size(nodes, dim=1)+ei(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(1:size(distance, dim=1), :, :) = distance
      do i=size(distance, dim=1)+1, size(distance, dim=1)+ei(1)
         distance_ext(i, :, :) = distance_ext(i-1, :, :)
      enddo
   endif
   if (ei(2)>=2) then
      call move_alloc(from=nodes_ext, to=nodes)
      call move_alloc(from=distance_ext, to=distance)
      allocate(nodes_ext(   size(nodes,    dim=1)+ei(2), size(nodes,    dim=2), size(nodes,    dim=3)))
      allocate(distance_ext(size(distance, dim=1)+ei(2), size(distance, dim=2), size(distance, dim=3)))
      nodes_ext(ei(2)+1:, :, :) = nodes
      Dxyz_ext = Dx
      do i=ei(2), 1, -1
         nodes_ext(i, :, :)%x = nodes_ext(i+1, :, :)%x - Dxyz_ext
         nodes_ext(i, :, :)%y = nodes_ext(i+1, :, :)%y
         nodes_ext(i, :, :)%z = nodes_ext(i+1, :, :)%z
         if (i>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(ei(2)+1:, :, :) = distance
      do i=ei(2), 1, -1
         distance_ext(i, :, :) = distance_ext(i+1, :, :)
      enddo
   endif
   if (ej(1)>=2) then
      call move_alloc(from=nodes_ext, to=nodes)
      call move_alloc(from=distance_ext, to=distance)
      allocate(nodes_ext(   size(nodes,    dim=1), size(nodes,    dim=2)+ej(1), size(nodes,    dim=3)))
      allocate(distance_ext(size(distance, dim=1), size(distance, dim=2)+ej(1), size(distance, dim=3)))
      nodes_ext(:, 1:size(nodes, dim=2), :) = nodes
      Dxyz_ext = Dy
      do j=size(nodes, dim=2)+1, size(nodes, dim=2)+ej(1)
         nodes_ext(:, j, :)%x = nodes_ext(:, j-1, :)%x
         nodes_ext(:, j, :)%y = nodes_ext(:, j-1, :)%y + Dxyz_ext
         nodes_ext(:, j, :)%z = nodes_ext(:, j-1, :)%z
         if (j<size(nodes, dim=2)+ej(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(:, 1:size(distance, dim=2), :) = distance
      do j=size(distance, dim=2)+1, size(distance, dim=2)+ej(1)
         distance_ext(:, j, :) = distance_ext(:, j-1, :)
      enddo
   endif
   if (ej(2)>=2) then
      call move_alloc(from=nodes_ext, to=nodes)
      call move_alloc(from=distance_ext, to=distance)
      allocate(nodes_ext(   size(nodes,    dim=1), size(nodes,    dim=2)+ej(2), size(nodes,    dim=3)))
      allocate(distance_ext(size(distance, dim=1), size(distance, dim=2)+ej(2), size(distance, dim=3)))
      nodes_ext(:, ej(2)+1:, :) = nodes
      Dxyz_ext = Dy
      do j=ej(2), 1, -1
         nodes_ext(:, j, :)%x = nodes_ext(:, j+1, :)%x
         nodes_ext(:, j, :)%y = nodes_ext(:, j+1, :)%y - Dxyz_ext
         nodes_ext(:, j, :)%z = nodes_ext(:, j+1, :)%z
         if (j>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(:, ej(2)+1:, :) = distance
      do j=ej(2), 1, -1
         distance_ext(:, j, :) = distance_ext(:, j+1, :)
      enddo
   endif
   if (ek(1)>=2) then
      call move_alloc(from=nodes_ext, to=nodes)
      call move_alloc(from=distance_ext, to=distance)
      allocate(nodes_ext(   size(nodes,    dim=1), size(nodes,    dim=2), size(nodes,    dim=3)+ek(1)))
      allocate(distance_ext(size(distance, dim=1), size(distance, dim=2), size(distance, dim=3)+ek(1)))
      nodes_ext(:, :, 1:size(nodes, dim=3)) = nodes
      Dxyz_ext = Dz
      do k=size(nodes, dim=3)+1, size(nodes, dim=3)+ek(1)
         nodes_ext(:, :, k)%x = nodes_ext(:, :, k-1)%x
         nodes_ext(:, :, k)%y = nodes_ext(:, :, k-1)%y
         nodes_ext(:, :, k)%z = nodes_ext(:, :, k-1)%z + Dxyz_ext
         if (k<size(nodes, dim=3)+ek(1)-2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(:, :, 1:size(distance, dim=3)) = distance
      do k=size(distance, dim=3)+1, size(distance, dim=3)+ek(1)
         distance_ext(:, :, k) = distance_ext(:, :, k-1)
      enddo
   endif
   if (ek(2)>=2) then
      call move_alloc(from=nodes_ext, to=nodes)
      call move_alloc(from=distance_ext, to=distance)
      allocate(nodes_ext(   size(nodes,    dim=1), size(nodes,    dim=2), size(nodes,    dim=3)+ek(2)))
      allocate(distance_ext(size(distance, dim=1), size(distance, dim=2), size(distance, dim=3)+ek(2)))
      nodes_ext(:, :, ek(2)+1:) = nodes
      Dxyz_ext = Dz
      do k=ek(2), 1, -1
         nodes_ext(:, :, k)%x = nodes_ext(:, :, k+1)%x
         nodes_ext(:, :, k)%y = nodes_ext(:, :, k+1)%y
         nodes_ext(:, :, k)%z = nodes_ext(:, :, k+1)%z - Dxyz_ext
         if (k>2) Dxyz_ext = 1.1_R8P * Dxyz_ext
      enddo
      distance_ext(:, :, ek(2)+1:) = distance
      do k=ek(2), 1, -1
         distance_ext(:, :, k) = distance_ext(:, :, k+1)
      enddo
   endif
   print '(A)', 'save output'
   call export_vtk_file(nodes_=nodes_ext, distance_=distance_ext, file_name=trim(output_base_name)//'.vts')
   call export_xall_files(nodes_=nodes_ext, distance_=distance_ext, basename=trim(output_base_name))
else
   print '(A)', 'save output'
   call export_vtk_file(nodes_=nodes, distance_=distance, file_name=trim(output_base_name)//'.vts')
   call export_xall_files(nodes_=nodes, distance_=distance, basename=trim(output_base_name))
endif

contains
  subroutine cli_parse()
  !< Build and parse test cli.
  real(R8P)    :: bb(3) !< Bounding box extents, local variable.
  integer(I4P) :: error !< Error trapping flag.

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
               def='0.0 0.0 -0.07',         &
               act='store')

  call cli%add(switch='--bmax',             &
               help='bounding box maximum', &
               required=.false.,            &
               nargs='3',                   &
               def='1.0 0.22 0.15',         &
               act='store')

  call cli%add(switch='--ni',                      &
               help='cells number in i direction', &
               required=.false.,                   &
               def='256',                          &
               act='store')

  call cli%add(switch='--nj',                      &
               help='cells number in j direction', &
               required=.false.,                   &
               def='128',                          &
               act='store')

  call cli%add(switch='--nk',                      &
               help='cells number in k direction', &
               required=.false.,                   &
               def='128',                          &
               act='store')

  call cli%add(switch='--gi',                            &
               help='ghost cells number in i direction', &
               required=.false.,                         &
               def='32',                                 &
               act='store')

  call cli%add(switch='--gj',                            &
               help='ghost cells number in j direction', &
               required=.false.,                         &
               nargs='2',                                &
               def='32 0',                               &
               act='store')

  call cli%add(switch='--gk',                            &
               help='ghost cells number in k direction', &
               required=.false.,                         &
               def='32',                                 &
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

  call cli%add(switch='--sign_algorithm',                         &
               help='algorithm used to compute sign of distance', &
               required=.false.,                                  &
               def='ray_intersections',                           &
               act='store')

  call cli%add(switch='--unsigned',              &
               help='compute unsigned distance', &
               required=.false.,                 &
               def='.false.',                    &
               act='store_true')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get(switch='--stl',                     val=file_name_stl,     error=error) ; if (error/=0) stop
  call cli%get(switch='--out',                     val=output_base_name,  error=error) ; if (error/=0) stop
  call cli%get(switch='--bmin',                    val=bb,                error=error) ; if (error/=0) stop
  bmin%x = bb(1) ; bmin%y = bb(2) ;bmin%z = bb(3)
  call cli%get(switch='--bmax',                    val=bb,                error=error) ; if (error/=0) stop
  bmax%x = bb(1) ; bmax%y = bb(2) ;bmax%z = bb(3)
  call cli%get(switch='--ni',                      val=ni,                error=error) ; if (error/=0) stop
  call cli%get(switch='--nj',                      val=nj,                error=error) ; if (error/=0) stop
  call cli%get(switch='--nk',                      val=nk,                error=error) ; if (error/=0) stop
  call cli%get(switch='--gi',                      val=gi,                error=error) ; if (error/=0) stop
  call cli%get(switch='--gj',                      val=gj,                error=error) ; if (error/=0) stop
  call cli%get(switch='--gk',                      val=gk,                error=error) ; if (error/=0) stop
  call cli%get(switch='--ei',                      val=ei,                error=error) ; if (error/=0) stop
  call cli%get(switch='--ej',                      val=ej,                error=error) ; if (error/=0) stop
  call cli%get(switch='--ek',                      val=ek,                error=error) ; if (error/=0) stop
  call cli%get(switch='--ref_levels',              val=refinement_levels, error=error) ; if (error/=0) stop
  call cli%get(switch='--sign_algorithm',          val=sign_algorithm,    error=error) ; if (error/=0) stop
  call cli%get(switch='--unsigned',                val=unsigned,          error=error) ; if (error/=0) stop
  endsubroutine cli_parse

   subroutine export_vtk_file(nodes_, distance_, file_name)
   !< Export IB data into VTK file format.
   type(vector_R8P), intent(in) :: nodes_(:,:,:)    !< Grid nodes.
   real(R8P),        intent(in) :: distance_(:,:,:) !< Distance of grid centers.
   character(*),     intent(in) :: file_name        !< File name.
   integer(I4P)                 :: ni, nj, nk, nn   !< Grid dimensions.

   ni = size(distance_, dim=1)
   nj = size(distance_, dim=2)
   nk = size(distance_, dim=3)
   nn = (ni+1) * (nj+1) * (nk+1)
   error = vtk%initialize(format='raw', filename=adjustl(trim(file_name)), mesh_topology='StructuredGrid', &
                          nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
   error = vtk%xml_writer%write_piece(nx1=0, nx2=ni, ny1=0, ny2=nj, nz1=0, nz2=nk)
   error = vtk%xml_writer%write_geo(n=nn, x=nodes_%x, y=nodes_%y, z=nodes_%z)
   error = vtk%xml_writer%write_dataarray(location='cell', action='open')
   error = vtk%xml_writer%write_dataarray(data_name='distance', x=distance_, one_component=.true.)
   error = vtk%xml_writer%write_dataarray(location='cell', action='close')
   error = vtk%xml_writer%write_piece()
   error = vtk%finalize()
   endsubroutine export_vtk_file

   subroutine export_xall_files(nodes_, distance_, basename)
   !< Export IB data into XALL files format.
   type(vector_R8P), intent(in) :: nodes_(-2:,-2:,-2:)    !< Grid nodes.
   real(R8P),        intent(in) :: distance_(-1:,-1:,-1:) !< Distance of grid centers.
   character(*),     intent(in) :: basename               !< Base files name.
   integer(I4P)                 :: ni, nj, nk             !< Grid dimensions.
   integer(I4P)                 :: funit                  !< File unit.

   ! number of cells subtracted ghost cells
   ni = size(distance_, dim=1) - 4
   nj = size(distance_, dim=2) - 4
   nk = size(distance_, dim=3) - 4
   ! grid
   open(newunit=funit, file=adjustl(trim(basename))//'.grd', form='UNFORMATTED', action='WRITE')
   write(funit) 1 ! number of blocks
   write(funit) ni, nj, nk
   write(funit) nodes_(0:ni,0:nj,0:nk)%x
   write(funit) nodes_(0:ni,0:nj,0:nk)%y
   write(funit) nodes_(0:ni,0:nj,0:nk)%z
   close(funit)

   ! distance
   open(newunit=funit, file=adjustl(trim(basename))//'.ib', form='UNFORMATTED', action='WRITE')
   write(funit) 0._R8P ! time
   write(funit) ni, nj, nk
   write(funit) distance_
   close(funit)
   endsubroutine export_xall_files
endprogram fossil_ib_generator
