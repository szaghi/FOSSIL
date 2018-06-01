!< FOSSIL,  STL file class definition.

module fossil_file_stl_object
!< FOSSIL,  STL file class definition.

! use fossil_aabb_tree_object, only : aabb_tree_object
use fossil_facet_object, only : facet_object
use fossil_utils, only : FRLEN, is_inside_bb
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P
use vecfor, only : vector_R8P

implicit none
private
public :: file_stl_object

type :: file_stl_object
   !< FOSSIL STL file class.
   character(len=:), allocatable :: file_name       !< File name
   integer(I4P)                  :: file_unit=0     !< File unit.
   character(FRLEN)              :: header=''       !< File header.
   logical                       :: is_ascii=.true. !< Sentinel to check if file is ASCII.
   logical                       :: is_open=.false. !< Sentinel to check if file is open.
   contains
      ! public methods
      procedure, pass(self) :: close_file     !< Close file.
      procedure, pass(self) :: destroy        !< Destroy file.
      procedure, pass(self) :: initialize     !< Initialize file.
      procedure, pass(self) :: load_from_file !< Load from file.
      procedure, pass(self) :: open_file      !< Open file, once initialized.
      procedure, pass(self) :: save_into_file !< Save into file.
      procedure, pass(self) :: statistics     !< Return STL statistics.
      ! operators
      generic :: assignment(=) => file_stl_assign_file_stl       !< Overload `=`.
      procedure, pass(lhs),  private :: file_stl_assign_file_stl !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: load_facets_number_from_file !< Load facets number from file.
      procedure, pass(self), private :: load_header_from_file        !< Load header from file.
      procedure, pass(self), private :: save_header_into_file        !< Save header into file.
      procedure, pass(self), private :: save_trailer_into_file       !< Save trailer into file.
endtype file_stl_object

contains
   ! public methods
   subroutine close_file(self)
   !< Close file.
   class(file_stl_object), intent(inout) :: self       !< File STL.
   logical                               :: file_exist !< Sentinel to check if file exist.

   if (self%is_open) then
      close(unit=self%file_unit)
      self%file_unit = 0
      self%is_open = .false.
   endif
   endsubroutine close_file

   elemental subroutine destroy(self)
   !< Destroy file.
   class(file_stl_object), intent(inout) :: self  !< File STL.
   type(file_stl_object)                 :: fresh !< Fresh instance of file STL.

   self = fresh
   endsubroutine destroy

   elemental subroutine initialize(self, skip_destroy, file_name, is_ascii)
   !< Initialize file.
   class(file_stl_object), intent(inout)        :: self          !< File STL.
   logical,                intent(in), optional :: skip_destroy  !< Flag to skip destroy file.
   character(*),           intent(in), optional :: file_name     !< File name.
   logical,                intent(in), optional :: is_ascii      !< Sentinel to check if file is ASCII.
   logical                                      :: skip_destroy_ !< Flag to skip destroy file, local variable.

   skip_destroy_ = .false. ; if (present(skip_destroy)) skip_destroy_ = skip_destroy
   if (.not.skip_destroy_) call self%destroy
   if (present(file_name)) self%file_name = trim(adjustl(file_name))
   if (present(is_ascii)) self%is_ascii = is_ascii
   endsubroutine initialize

   subroutine load_from_file(self, facet, file_name, is_ascii, guess_format, clip_min, clip_max)
   !< Load from file.
   class(file_stl_object),          intent(inout)        :: self               !< File STL.
   type(facet_object), allocatable, intent(out)          :: facet(:)           !< Surface facets.
   character(*),                    intent(in), optional :: file_name          !< File name.
   logical,                         intent(in), optional :: is_ascii           !< Sentinel to check if file is ASCII.
   logical,                         intent(in), optional :: guess_format       !< Try to guess format directly from file.
   type(vector_R8P),                intent(in), optional :: clip_min, clip_max !< Clip bounding box extents.
   integer(I4P)                                          :: facets_number      !< Facets number.
   type(facet_object)                                    :: facet_clip         !< Buffer for clipped loading.
   integer(I4P)                                          :: f, ff              !< Counter.

   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='read', guess_format=guess_format)
   call self%load_facets_number_from_file(facets_number=facets_number)
   call self%load_header_from_file
   if (present(clip_min).and.present(clip_max)) then
      call facet_clip%set_io_methods(is_ascii=self%is_ascii)
      ! count the facets that are actually inside the clipping bounding box
      ff = 0
      do f=1, facets_number
         call facet_clip%load_from_file(file_unit=self%file_unit)
         if (is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(1)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(2)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(3))) then
            ff = ff + 1
         endif
      enddo
      call self%load_header_from_file
      allocate(facet(1:ff))
      call facet%set_io_methods(is_ascii=self%is_ascii)
      ff = 0
      do f=1, facets_number
         call facet_clip%load_from_file(file_unit=self%file_unit)
         if (is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(1)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(2)).and.&
             is_inside_bb(bmin=clip_min, bmax=clip_max, point=facet_clip%vertex(3))) then
            ff = ff + 1
            facet(ff) = facet_clip
            facet(ff)%id = ff
         endif
      enddo
   else
      allocate(facet(1:facets_number))
      call facet%set_io_methods(is_ascii=self%is_ascii)
      do f=1, facets_number
         call facet(f)%load_from_file(file_unit=self%file_unit)
         facet(f)%id = f
      enddo
   endif
   call self%close_file
   endsubroutine load_from_file

   subroutine open_file(self, file_action, guess_format)
   !< Open file, once initialized.
   class(file_stl_object), intent(inout)        :: self          !< File STL.
   character(*),           intent(in)           :: file_action   !< File action, "read" or "write".
   logical,                intent(in), optional :: guess_format  !< Sentinel to try to guess format directly from file.
   logical                                      :: guess_format_ !< Sentinel to try to guess format directly from file, local var.
   logical                                      :: file_exist    !< Sentinel to check if file exist.
   character(5)                                 :: ascii_header  !< Ascii header sentinel.

   if (allocated(self%file_name)) then
      select case(trim(adjustl(file_action)))
      case('read')
         guess_format_ = .false. ; if (present(guess_format)) guess_format_ = guess_format
         inquire(file=self%file_name, exist=file_exist)
         if (file_exist) then
            if (guess_format_) then
               open(newunit=self%file_unit, file=self%file_name, form='formatted')
               read(self%file_unit, '(A)') ascii_header
               close(self%file_unit)
               if (ascii_header=='solid') then
                  self%is_ascii = .true.
               else
                  self%is_ascii = .false.
               endif
            endif
            if (self%is_ascii) then
               open(newunit=self%file_unit, file=self%file_name,                  form='formatted')
            else
               open(newunit=self%file_unit, file=self%file_name, access='stream', form='unformatted')
            endif
            self%is_open = .true.
         else
            write(stderr, '(A)') 'error: file "'//self%file_name//'" does not exist, impossible to open file!'
         endif
      case('write')
         if (self%is_ascii) then
            open(newunit=self%file_unit, file=self%file_name,                  form='formatted')
         else
            open(newunit=self%file_unit, file=self%file_name, access='stream', form='unformatted')
         endif
         self%is_open = .true.
      case default
         write(stderr, '(A)') 'error: file action "'//trim(adjustl(file_action))//'" unknown!'
      endselect
   else
      write(stderr, '(A)') 'error: file name has not be initialized, impossible to open file!'
   endif
   endsubroutine open_file

   subroutine save_into_file(self, facet, file_name, is_ascii)
   !< Save into file.
   class(file_stl_object),    intent(inout)        :: self          !< File STL.
   type(facet_object),        intent(inout)        :: facet(:)      !< Surface facets.
   character(*),              intent(in), optional :: file_name     !< File name.
   logical,                   intent(in), optional :: is_ascii      !< Sentinel to check if file is ASCII.
   integer(I4P)                                    :: facets_number !< Facets number.
   integer(I4P)                                    :: f             !< Counter.

   facets_number = size(facet, dim=1)
   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='write')
   call self%save_header_into_file(facets_number=facets_number)
   call facet%set_io_methods(is_ascii=self%is_ascii)
   do f=1, facets_number
      call facet(f)%save_into_file(file_unit=self%file_unit)
   enddo
   call self%save_trailer_into_file
   call self%close_file
   endsubroutine save_into_file

   pure function statistics(self, prefix) result(stats)
   !< Return STL statistics.
   class(file_stl_object), intent(in)           :: self             !< File STL.
   character(*),           intent(in), optional :: prefix           !< Lines prefix.
   character(len=:), allocatable                :: stats            !< STL statistics.
   character(len=:), allocatable                :: prefix_          !< Lines prefix, local variable.
   character(1), parameter                      :: NL=new_line('a') !< Line terminator.

   prefix_ = '' ; if (present(prefix)) prefix_ = prefix
   stats = prefix_//self%header//NL
   if (allocated(self%file_name)) stats=stats//prefix_//'file name:   '//self%file_name//NL
   if (self%is_ascii) then
      stats=stats//prefix_//'file format: ascii'
   else
      stats=stats//prefix_//'file format: binary'
   endif
   endfunction statistics

   ! operators
   ! =
   pure subroutine file_stl_assign_file_stl(lhs, rhs)
   !< Operator `=`.
   class(file_stl_object), intent(inout) :: lhs !< Left hand side.
   type(file_stl_object),  intent(in)    :: rhs !< Right hand side.
   integer(I4P)                          :: f   !< Counter.

   if (allocated(lhs%file_name)) deallocate(lhs%file_name)
   if (allocated(rhs%file_name)) lhs%file_name = rhs%file_name
   lhs%file_unit = rhs%file_unit
   lhs%header = rhs%header
   endsubroutine file_stl_assign_file_stl

   ! private methods
   subroutine load_facets_number_from_file(self, facets_number)
   !< Load facets number from file.
   !<
   !< @note File is rewinded.
   class(file_stl_object), intent(inout) :: self          !< File STL.
   integer(I4P),           intent(out)   :: facets_number !< Facets number.
   character(FRLEN)                      :: facet_record  !< Facet record string buffer.

   facets_number = 0
   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         do
            read(self%file_unit, '(A)', end=10, err=10) facet_record
            if (index(string=facet_record, substring='facet normal') > 0) facets_number = facets_number + 1
         enddo
      else
         read(self%file_unit, end=10, err=10) facet_record
         read(self%file_unit, end=10, err=10) facets_number
      endif
      10 rewind(self%file_unit)
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load facets number from file!'
   endif
   endsubroutine load_facets_number_from_file

   subroutine load_header_from_file(self)
   !< Load header from file.
   class(file_stl_object), intent(inout) :: self          !< File STL.
   integer(I4P)                          :: facets_number !< Facets number.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         read(self%file_unit, '(A)') self%header
         self%header = trim(adjustl(self%header(index(self%header, 'solid')+6:)))
      else
         read(self%file_unit) self%header
         read(self%file_unit) facets_number ! actually read elsewhere
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine load_header_from_file

   subroutine save_header_into_file(self, facets_number)
   !< Save header into file.
   class(file_stl_object), intent(inout) :: self          !< File STL.
   integer(I4P),           intent(in)    :: facets_number !< Facets number.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         write(self%file_unit, '(A)') 'solid '//trim(self%header)
      else
         write(self%file_unit) self%header
         write(self%file_unit) facets_number
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine save_header_into_file

   subroutine save_trailer_into_file(self)
   !< Save trailer into file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      if (self%is_ascii) write(self%file_unit, '(A)') 'endsolid '//trim(self%header)
   else
      write(stderr, '(A)') 'error: file is not open, impossible to write trailer into file!'
   endif
   endsubroutine save_trailer_into_file
endmodule fossil_file_stl_object
