!< FOSSIL, generate Immersed Boundary distance function.

program fossil_ib_generator
!< FOSSIL, generate Immersed Boundary distance function.
!<
!< Signed distance of the cell centres of a Cartesian block to an STL surface, written as a
!< VTK structured grid (and, without extension cells, as Xall `.grd` / `.ib` files). The
!< block is the body-close grid (`--ni/--nj/--nk` cells over `--bmin`..`--bmax`), a frame of
!< `--gi/--gj/--gk` cells with the same spacing, and optional extension cells
!< (`--ei/--ej/--ek`) whose spacing grows by 1.1.
!<
!< The field comes from `distance_field_object` (issue #32): exact distances in a narrow band
!< around the surface (`--band`, default twice the local spacing), signs outside it by flood
!< fill, and the far field by fast sweeping (`--sweep_tol`), so every cell, extension cells
!< included, holds a distance. `--device` computes the band on the GPU (nvfortran build;
!< sign algorithm pseudo_normal or fast_winding_number). The written field is positive
!< inside the body (Xall convention).

use flap, only : command_line_interface
use fossil, only : surface_stl_object, sign_algorithm_from_string, STATUS_OK, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER
use fossil_distance_field, only : distance_field_object
use penf, only : I4P, I8P, R8P, str
use vecfor, only : vector_R8P
use vtk_fortran, only : vtk_file

implicit none

type(command_line_interface)   :: cli                  !< Test command line interface.
type(surface_stl_object)       :: surface_stl          !< STL surface.
type(distance_field_object)    :: field                !< Signed-distance field.
type(vector_R8P), allocatable  :: nodes(:,:,:)         !< Grid nodes.
real(R8P),        allocatable  :: distance(:,:,:)      !< Distance of grid centers to STL surface (positive inside).
real(R8P),        allocatable  :: xn(:), yn(:), zn(:)  !< Node coordinates per axis.
real(R8P),        allocatable  :: xc(:), yc(:), zc(:)  !< Centre coordinates per axis.
type(vector_R8P)               :: bmin, bmax           !< Bounding box extents.
character(999)                 :: file_name_stl        !< Input STL file name.
character(999)                 :: output_base_name     !< Output base name.
integer(I4P)                   :: refinement_levels    !< AABB refinement levels used.
character(999)                 :: sign_algorithm       !< Algorithm name as typed on the CLI.
integer(I4P)                   :: sign_algorithm_code  !< Algorithm code derived from sign_algorithm.
logical                        :: unsigned             !< Compute unsigned distance.
logical                        :: use_device           !< Narrow band on the GPU.
real(R8P)                      :: band                 !< Band half-width (<= 0: default).
real(R8P)                      :: sweep_tol            !< Fast-sweeping relative tolerance.
integer(I4P)                   :: ni, nj, nk           !< Grid dimensions, close to STL.
integer(I4P)                   :: gi, gj(2), gk        !< Frame dimensions around the grid.
integer(I4P)                   :: ei(2), ej(2), ek(2)  !< Extension cells (high side, low side).
integer(I4P)                   :: i, j, k              !< Counter.
integer(I4P)                   :: sweeps, fallback     !< Field counters.
type(vtk_file)                 :: vtk                  !< VTK file handler.
integer(I4P)                   :: error                !< Status error.
integer(I8P)                   :: timing(0:2)          !< Tic toc timing.

! parse command line input e load STL file
call cli_parse
call cli_print
call surface_stl%load_from_file(file_name=trim(adjustl(file_name_stl)), guess_format=.true., &
                                aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics before sanitization'
print '(A)', surface_stl%statistics()
call surface_stl%sanitize
call surface_stl%analyze(aabb_refinement_levels=refinement_levels)
print '(A)', 'STL statistics after sanitization'
print '(A)', surface_stl%statistics()

if (.not.cli%is_passed(switch='--bmin')) bmin = surface_stl%get_bmin()
if (.not.cli%is_passed(switch='--bmax')) bmax = surface_stl%get_bmax()
if (use_device .and. sign_algorithm_code /= SIGN_PSEUDO_NORMAL .and. sign_algorithm_code /= SIGN_FAST_WINDING_NUMBER) then
   write(*, '(A)') 'fossil_ib_generator: --device needs --sign_algorithm pseudo_normal or fast_winding_number'
   error stop 1
endif

! grid: per-axis nodes (body-close + frame + extension) and cell centres
xn = node_axis(bmin%x, bmax%x, ni, gi,    gi,    ei)
yn = node_axis(bmin%y, bmax%y, nj, gj(1), gj(2), ej)
zn = node_axis(bmin%z, bmax%z, nk, gk,    gk,    ek)
xc = centres(xn)
yc = centres(yn)
zc = centres(zn)

! signed-distance field: band, signs, fast sweeping
print '(A)', 'compute distance field'
call system_clock(timing(1))
call field%initialize(x=xc, y=yc, z=zc, status=error)
if (error /= STATUS_OK) error stop 'fossil_ib_generator: invalid grid'
if (band > 0._R8P) then
   call field%compute(surface=surface_stl, band=band, is_signed=.not.unsigned, sign_algorithm=sign_algorithm_code, &
                      use_device=use_device, tol=sweep_tol, sweeps=sweeps, fallback_cells=fallback, status=error)
else
   call field%compute(surface=surface_stl, is_signed=.not.unsigned, sign_algorithm=sign_algorithm_code, &
                      use_device=use_device, tol=sweep_tol, sweeps=sweeps, fallback_cells=fallback, status=error)
endif
if (error /= STATUS_OK) error stop 'fossil_ib_generator: distance field failed'
call system_clock(timing(2), timing(0))
print '(A,F8.3,A,I0,A,ES10.3,A,I0,A,I0)', 'distance field timing: ', real(timing(2) - timing(1)) / timing(0), &
      ' s  band cells ', field%get_band_cells(), '  band ', field%get_band_width(), '  sweeps ', sweeps, &
      '  fallback cells ', fallback

! output, positive inside (Xall convention)
allocate(distance(lbound(xc, 1):ubound(xc, 1), lbound(yc, 1):ubound(yc, 1), lbound(zc, 1):ubound(zc, 1)))
distance = -field%phi
allocate(nodes(lbound(xn, 1):ubound(xn, 1), lbound(yn, 1):ubound(yn, 1), lbound(zn, 1):ubound(zn, 1)))
do k = lbound(zn, 1), ubound(zn, 1)
   do j = lbound(yn, 1), ubound(yn, 1)
      do i = lbound(xn, 1), ubound(xn, 1)
         nodes(i, j, k)%x = xn(i) ; nodes(i, j, k)%y = yn(j) ; nodes(i, j, k)%z = zn(k)
      enddo
   enddo
enddo
print '(A)', 'save output'
call export_vtk_file(nodes_=nodes, distance_=distance, file_name=trim(output_base_name)//'.vts')
if (all([ei, ej, ek] < 2) .and. min(gi, gj(1), gj(2), gk) >= 2) &
   call export_xall_files(ni=ni,nj=nj,nk=nk,gi=gi,gj=gj,gk=gk,nodes_=nodes,distance_=distance,basename=trim(output_base_name))

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
               def='2',                                  &
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
               def='2',                                  &
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

  call cli%add(switch='--ref_levels',                                                     &
               help='AABB refinement levels (-1 = auto-tune from facet count)',           &
               required=.false.,                                                          &
               def='-1',                                                                  &
               act='store')

  call cli%add(switch='--sign_algorithm',                         &
               help='algorithm used to compute sign of distance: ray_intersections, solid_angle, '// &
                    'pseudo_normal, fast_winding_number (most robust on surfaces with holes)', &
               required=.false.,                                  &
               def='ray_intersections',                           &
               act='store')

  call cli%add(switch='--unsigned',              &
               help='compute unsigned distance', &
               required=.false.,                 &
               def='.false.',                    &
               act='store_true')

  call cli%add(switch='--band',                                                              &
               help='narrow-band half-width with exact distances (<= 0: twice the local spacing)', &
               required=.false.,                                                             &
               def='0.0',                                                                    &
               act='store')

  call cli%add(switch='--sweep_tol',                                               &
               help='fast-sweeping tolerance, relative to the domain diagonal',    &
               required=.false.,                                                   &
               def='1.e-12',                                                       &
               act='store')

  call cli%add(switch='--device',                                                         &
               help='narrow band on the GPU (nvfortran build; pseudo_normal or fast_winding_number)', &
               required=.false.,                                                          &
               def='.false.',                                                             &
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
  sign_algorithm_code = sign_algorithm_from_string(trim(adjustl(sign_algorithm)))
  call cli%get(switch='--unsigned',                val=unsigned,          error=error) ; if (error/=0) stop
  call cli%get(switch='--band',                    val=band,              error=error) ; if (error/=0) stop
  call cli%get(switch='--sweep_tol',               val=sweep_tol,         error=error) ; if (error/=0) stop
  call cli%get(switch='--device',                  val=use_device,        error=error) ; if (error/=0) stop
  endsubroutine cli_parse

  subroutine cli_print()
  !< Print cli options.
  real(R8P)    :: bb(3) !< Bounding box extents, local variable.
  integer(I4P) :: error !< Error trapping flag.

  print  '(A)' ,'CLI options'
  print  '(A)' ,'--stl             '//trim(    file_name_stl          )
  print  '(A)' ,'--out             '//trim(    output_base_name       )
  print  '(A)' ,'--bmin            '//trim(str([bmin%x,bmin%y,bmin%z]))
  print  '(A)' ,'--bmax            '//trim(str([bmax%x,bmax%y,bmax%z]))
  print  '(A)' ,'--ni              '//trim(str(ni                    ))
  print  '(A)' ,'--nj              '//trim(str(nj                    ))
  print  '(A)' ,'--nk              '//trim(str(nk                    ))
  print  '(A)' ,'--gi              '//trim(str(gi                    ))
  print  '(A)' ,'--gj              '//trim(str(gj                    ))
  print  '(A)' ,'--gk              '//trim(str(gk                    ))
  print  '(A)' ,'--ei              '//trim(str(ei                    ))
  print  '(A)' ,'--ej              '//trim(str(ej                    ))
  print  '(A)' ,'--ek              '//trim(str(ek                    ))
  print  '(A)' ,'--ref_levels      '//trim(str(refinement_levels     ))
  print  '(A)' ,'--sign_algorithm  '//trim(    sign_algorithm         )
  print  '(A)' ,'--unsigned        '//trim(str(unsigned              ))
  print  '(A)' ,'--band            '//trim(str(band                  ))
  print  '(A)' ,'--sweep_tol       '//trim(str(sweep_tol             ))
  print  '(A)' ,'--device          '//trim(str(use_device            ))
  endsubroutine cli_print

   function node_axis(lo, hi, n, g_lo, g_hi, e) result(c)
   !< Node coordinates of one axis, indexed from -g_lo - e_lo to n + g_hi + e_hi: n uniform cells
   !< over [lo, hi], frames of g_lo / g_hi cells with the same spacing, then extension cells on
   !< the high side (e(1)) and the low side (e(2)). The extension spacing starts at the uniform
   !< one and grows by 1.1, except for the last cells, and a side with e < 2 gets none: the
   !< spacing rules of the former generator, kept so grids are reproduced.
   real(R8P),    intent(in) :: lo, hi     !< Body-close extent.
   integer(I4P), intent(in) :: n          !< Body-close cells.
   integer(I4P), intent(in) :: g_lo, g_hi !< Frame cells (low, high).
   integer(I4P), intent(in) :: e(2)       !< Extension cells (high side, low side).
   real(R8P), allocatable   :: c(:)       !< Nodes.
   real(R8P)                :: dd, h      !< Uniform and extension spacing.
   integer(I4P)             :: e_hi, e_lo, m

   e_hi = 0 ; if (e(1) >= 2) e_hi = e(1)
   e_lo = 0 ; if (e(2) >= 2) e_lo = e(2)
   dd = (hi - lo) / n
   allocate(c(-g_lo - e_lo:n + g_hi + e_hi))
   do m = -g_lo, n + g_hi
      c(m) = lo + m * dd
   enddo
   h = dd
   do m = 1, e_hi
      c(n + g_hi + m) = c(n + g_hi + m - 1) + h
      if (m < e_hi - 2) h = 1.1_R8P * h
   enddo
   h = dd
   do m = 1, e_lo
      c(-g_lo - m) = c(-g_lo - m + 1) - h
      if (m < e_lo - 1) h = 1.1_R8P * h
   enddo
   endfunction node_axis

   function centres(c) result(x)
   !< Cell centres of an axis: midpoints of consecutive nodes, cell i between nodes i-1 and i.
   real(R8P), intent(in)  :: c(:) !< Nodes (assumed shape: indexed from 1 here).
   real(R8P), allocatable :: x(:) !< Centres.

   x = 0.5_R8P * (c(1:size(c) - 1) + c(2:size(c)))
   endfunction centres

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

   subroutine export_xall_files(ni, nj, nk, gi, gj, gk, nodes_, distance_, basename)
   !< Export IB data into XALL files format.
   integer(I4P),     intent(in) :: ni, nj, nk                      !< Grid dimensions.
   integer(I4P),     intent(in) :: gi, gj(2), gk                   !< Frame dimensions around the grid.
   type(vector_R8P), intent(in) :: nodes_(   0-gi:,0-gj(1):,0-gk:) !< Grid nodes.
   real(R8P),        intent(in) :: distance_(1-gi:,1-gj(1):,1-gk:) !< Distance of grid centers.
   character(*),     intent(in) :: basename                        !< Base files name.
   integer(I4P)                 :: funit                           !< File unit.

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
   write(funit) distance_(-1:ni+2,-1:nj+2,-1:nk+2)
   close(funit)
   endsubroutine export_xall_files
endprogram fossil_ib_generator
