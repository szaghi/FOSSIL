!< FOSSIL, generate Immersed Boundary distance function.

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

call cblock%initialize(bmin=bmin, bmax=bmax, ni=ni, nj=nj, nk=nk, gi=gi, gj=gj, gk=gk, ei=ei, ej=ej, ek=ek, &
                       refinement_levels=refinement_levels)
call cblock%export_aabb_tecplot_ascii(file_name='ib-cart-block-aabb.dat')
! call cblock%compute_distances(surface_stl=surface_stl, is_signed=is_signed, sign_algorithm=sign_algorithm, invert_sign=.true.)
call cblock%export_vtk_file(file_name=trim(output_base_name)//'.vts')
! call cblock%export_xall_files(basename=trim(output_base_name))

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
