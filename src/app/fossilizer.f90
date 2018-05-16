!< FOSSILer, analyze, modify and repair STL files.

program fossiler
!< FOSSILer, analyze, modify and repair STL files.

use flap, only : command_line_interface
use fossil, only : file_stl_object
use penf, only : I4P, R8P, str
use vecfor, only : vector_R8P

implicit none

type(command_line_interface)       :: cli                                   !< Command line interface.
type(file_stl_object), allocatable :: file_stl(:)                           !< STL files.
character(999),        allocatable :: file_name_stl(:)                      !< Input STL file name.
logical                            :: connect_nearby                        !< Connect nearby vertices belong to disconnect edges.
character(999)                     :: output                                !< Output file name.
integer(I4P)                       :: stl_numbers                           !< Number of STL files.
integer(I4P)                       :: s                                     !< Counter.
type(vector_R8P)                   :: mirror_normal                         !< Normal of mirroring plane.
type(vector_R8P)                   :: resize_factor                         !< Resize factor.
real(R8P)                          :: resize_x, resize_y, resize_z          !< Scalar resize factors.
type(vector_R8P)                   :: rotate_axis                           !< Axis of rotation.
real(R8P)                          :: rotate_angle                          !< Angle of rotation.
logical                            :: sanitize                              !< Sanitize STL.
type(vector_R8P)                   :: translate_delta                       !< Translate delta.
real(R8P)                          :: translate_x, translate_y, translate_z !< Translate scalar deltas.

call cli_parse
do s=1, stl_numbers
   call file_stl(s)%load_from_file(file_name=trim(adjustl(file_name_stl(s))), guess_format=.true.)
   print '(A)', file_stl(s)%statistics()
   if (connect_nearby) then
      print '(A)', 'connect nearby vertices belong to disconnected edges'
      call file_stl(s)%connect_nearby_vertices
      print '(A)', file_stl(s)%statistics()
   endif
   if (cli%is_passed(switch='--mirror_normal')) then
      print '(A)', 'mirror surface with respect mirroring plane normal: '//trim(str(mirror_normal%x))//', '//&
                                                                           trim(str(mirror_normal%y))//', '//&
                                                                           trim(str(mirror_normal%z))
      call file_stl(s)%mirror(normal=mirror_normal)
   endif
   if (cli%is_passed(switch='--resize_factor')) then
      print '(A)', 'resize surface by: '//trim(str(resize_factor%x))//', '//&
                                          trim(str(resize_factor%y))//', '//&
                                          trim(str(resize_factor%z))
      call file_stl(s)%resize(factor=resize_factor)
   endif
   if (cli%is_passed(switch='--resize_x')) then
      print '(A)', 'resize surface along x by: '//trim(str(resize_x))
      call file_stl(s)%resize(x=resize_x)
   endif
   if (cli%is_passed(switch='--resize_y')) then
      print '(A)', 'resize surface along y by: '//trim(str(resize_y))
      call file_stl(s)%resize(y=resize_y)
   endif
   if (cli%is_passed(switch='--resize_z')) then
      print '(A)', 'resize surface along z by: '//trim(str(resize_z))
      call file_stl(s)%resize(z=resize_z)
   endif
   if (cli%is_passed(switch='--rotate_axis').and.cli%is_passed(switch='--rotate_angle')) then
      print '(A)', 'rotate surface around: '//trim(str(rotate_axis%x))//', '//&
                                              trim(str(rotate_axis%y))//', '//&
                                              trim(str(rotate_axis%z))//' by '//trim(str(rotate_angle))//'[rad]'
      call file_stl(s)%rotate(axis=rotate_axis, angle=rotate_angle)
   endif
   if (cli%is_passed(switch='--translate_delta')) then
      print '(A)', 'translate surface by: '//trim(str(translate_delta%x))//', '//&
                                             trim(str(translate_delta%y))//', '//&
                                             trim(str(translate_delta%z))
      call file_stl(s)%translate(delta=translate_delta)
   endif
   if (cli%is_passed(switch='--translate_x')) then
      print '(A)', 'translate surface along x by: '//trim(str(translate_x))
      call file_stl(s)%translate(x=translate_x)
   endif
   if (cli%is_passed(switch='--translate_y')) then
      print '(A)', 'translate surface along y by: '//trim(str(translate_y))
      call file_stl(s)%translate(y=translate_y)
   endif
   if (cli%is_passed(switch='--translate_z')) then
      print '(A)', 'translate surface along z by: '//trim(str(translate_z))
      call file_stl(s)%translate(z=translate_z)
   endif
   if (sanitize) then
      print '(A)', 'sanitize STL'
      call file_stl(s)%sanitize
      print '(A)', file_stl(s)%statistics()
   endif
   if (cli%is_passed(switch='--output')) then
      print '(A)', 'save output in: '//trim(adjustl(output))
      call file_stl(s)%save_into_file(file_name=adjustl(output))
   endif
   print*
enddo

contains
  subroutine cli_parse()
  !< Build and parse command line interface.
  integer(I4P) :: error               !< Error trapping flag.
  real(R8P)    :: mirror_normal_(3)   !< Normal of mirroring plane.
  real(R8P)    :: resize_factor_(3)   !< Resize factor.
  real(R8P)    :: rotate_axis_(3)     !< Axis of rotation.
  real(R8P)    :: translate_delta_(3) !< Translate delta.

  call cli%init(progname='fossilizer',                           &
                authors='S. Zaghi',                              &
                help='Usage: ',                                  &
                examples=["fossilizer -i src/tests/dragon.stl"], &
                epilog=new_line('a')//"all done")

  call cli%add(switch='-i',                   &
               help='STL (input) file names', &
               nargs='+',                     &
               required=.true.,               &
               act='store')

  call cli%add(switch='--connect_nearby',                                   &
               help='connect_nearby vertices belong to disconnected edges', &
               required=.false.,                                            &
               def='.false.',                                               &
               act='store_true')

  call cli%add(switch='--mirror_normal',         &
               help='normal of mirroring plane', &
               required=.false.,                 &
               nargs='+',                        &
               def='1.0 0.0 0.0',                &
               act='store')

  call cli%add(switch='--output',           &
               help='output file name',     &
               required=.false.,            &
               def='fossilizer_output.stl', &
               act='store')

  call cli%add(switch='--resize_factor', &
               help='resize factor',     &
               required=.false.,         &
               nargs='+',                &
               def='2.0 2.0 2.0',        &
               act='store')

  call cli%add(switch='--resize_x',    &
               help='resize factor x', &
               required=.false.,       &
               def='2.0',              &
               act='store')

  call cli%add(switch='--resize_y',    &
               help='resize factor y', &
               required=.false.,       &
               def='2.0',              &
               act='store')

  call cli%add(switch='--resize_z',    &
               help='resize factor z', &
               required=.false.,       &
               def='2.0',              &
               act='store')

  call cli%add(switch='--rotate_axis',  &
               help='axis of rotation', &
               required=.false.,        &
               nargs='+',               &
               def='1.0 0.0 0.0',       &
               act='store')

  call cli%add(switch='--rotate_angle',  &
               help='angle of rotation', &
               required=.false.,         &
               def='1.57079633',         &
               act='store')

  call cli%add(switch='--sanitize', &
               help='sanitize STL', &
               required=.false.,    &
               def='.false.',       &
               act='store_true')

  call cli%add(switch='--translate_delta', &
               help='translate delta',     &
               required=.false.,           &
               nargs='+',                  &
               def='1.0 0.0 0.0',          &
               act='store')

  call cli%add(switch='--translate_x', &
               help='translate x',     &
               required=.false.,       &
               def='1.0',              &
               act='store')

  call cli%add(switch='--translate_y', &
               help='translate y',     &
               required=.false.,       &
               def='1.0',              &
               act='store')

  call cli%add(switch='--translate_z', &
               help='translate z',     &
               required=.false.,       &
               def='1.0',              &
               act='store')

  call cli%parse(error=error) ; if (error/=0) stop

  call cli%get_varying(switch='-i',                val=file_name_stl,    error=error) ; if (error/=0) stop
  call cli%get(        switch='--connect_nearby',  val=connect_nearby,   error=error) ; if (error/=0) stop
  call cli%get(        switch='--mirror_normal',   val=mirror_normal_,   error=error) ; if (error/=0) stop
  call cli%get(        switch='--output',          val=output,           error=error) ; if (error/=0) stop
  call cli%get(        switch='--resize_factor',   val=resize_factor_,   error=error) ; if (error/=0) stop
  call cli%get(        switch='--resize_x',        val=resize_x,         error=error) ; if (error/=0) stop
  call cli%get(        switch='--resize_y',        val=resize_y,         error=error) ; if (error/=0) stop
  call cli%get(        switch='--resize_z',        val=resize_z,         error=error) ; if (error/=0) stop
  call cli%get(        switch='--rotate_axis',     val=rotate_axis_,     error=error) ; if (error/=0) stop
  call cli%get(        switch='--rotate_angle',    val=rotate_angle,     error=error) ; if (error/=0) stop
  call cli%get(        switch='--sanitize',        val=sanitize,         error=error) ; if (error/=0) stop
  call cli%get(        switch='--translate_delta', val=translate_delta_, error=error) ; if (error/=0) stop
  call cli%get(        switch='--translate_x',     val=translate_x,      error=error) ; if (error/=0) stop
  call cli%get(        switch='--translate_y',     val=translate_y,      error=error) ; if (error/=0) stop
  call cli%get(        switch='--translate_z',     val=translate_z,      error=error) ; if (error/=0) stop

  stl_numbers = size(file_name_stl, dim=1)
  allocate(file_stl(1:stl_numbers))
  mirror_normal%x = mirror_normal_(1)
  mirror_normal%y = mirror_normal_(2)
  mirror_normal%z = mirror_normal_(3)
  resize_factor%x = resize_factor_(1)
  resize_factor%y = resize_factor_(2)
  resize_factor%z = resize_factor_(3)
  rotate_axis%x = rotate_axis_(1)
  rotate_axis%y = rotate_axis_(2)
  rotate_axis%z = rotate_axis_(3)
  translate_delta%x = translate_delta_(1)
  translate_delta%y = translate_delta_(2)
  translate_delta%z = translate_delta_(3)
  endsubroutine cli_parse
endprogram fossiler
