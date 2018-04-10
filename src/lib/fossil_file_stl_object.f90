!< FOSSIL,  STL file class definition.

module fossil_file_stl_object
!< FOSSIL,  STL file class definition.

use fossil_facet_object, only : facet_object, FRLEN
use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I4P

implicit none
private
public :: file_stl_object

type :: file_stl_object
   !< FOSSIL STL file class.
   character(len=:), allocatable   :: file_name       !< File name
   integer(I4P)                    :: file_unit=0     !< File unit.
   character(FRLEN)                :: header          !< File header.
   integer(I4P)                    :: facets_number=0 !< Facets number.
   type(facet_object), allocatable :: facet(:)        !< Facets.
   logical                         :: is_ascii=.true. !< Sentinel to check if file is ASCII.
   logical                         :: is_open=.false. !< Sentinel to check if file is open.
   contains
      ! public methods
      procedure, pass(self) :: close_file      !< Close file.
      procedure, pass(self) :: destroy         !< Destroy file.
      procedure, pass(self) :: initialize      !< Initialize file.
      procedure, pass(self) :: load_from_file  !< Load from file.
      procedure, pass(self) :: open_file       !< Open file, once initialized.
      procedure, pass(self) :: save_into_file  !< Save into file.
      ! operators
      generic :: assignment(=) => file_stl_assign_file_stl       !< Overload `=`.
      procedure, pass(lhs),  private :: file_stl_assign_file_stl !< Operator `=`.
      ! private methods
      procedure, pass(self), private :: allocate_facets              !< Allocate facets.
      procedure, pass(self), private :: load_facets_number_from_file !< Load facets number from file.
      procedure, pass(self), private :: load_header_from_file        !< Load header from file.
      procedure, pass(self), private :: save_header_into_file        !< Save header into file.
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

   subroutine load_from_file(self, file_name, is_ascii)
   !< Load from file.
   class(file_stl_object), intent(inout)        :: self      !< File STL.
   character(*),           intent(in), optional :: file_name !< File name.
   logical,                intent(in), optional :: is_ascii  !< Sentinel to check if file is ASCII.
   integer(I4P)                                 :: f         !< Counter.

   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='read')
   call self%load_facets_number_from_file
   call self%allocate_facets
   call self%load_header_from_file
   if (self%is_ascii) then
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_ascii(file_unit=self%file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%load_from_file_binary(file_unit=self%file_unit)
      enddo
   endif
   call self%close_file
   endsubroutine load_from_file

   subroutine open_file(self, file_action)
   !< Open file, once initialized.
   class(file_stl_object), intent(inout) :: self        !< File STL.
   character(*),           intent(in)    :: file_action !< File action, "read" or "write".
   logical                               :: file_exist  !< Sentinel to check if file exist.

   if (allocated(self%file_name)) then
      select case(trim(adjustl(file_action)))
      case('read')
         inquire(file=self%file_name, exist=file_exist)
         if (file_exist) then
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

   subroutine save_into_file(self, file_name, is_ascii)
   !< Save into file.
   class(file_stl_object), intent(inout)        :: self      !< File STL.
   character(*),           intent(in), optional :: file_name !< File name.
   logical,                intent(in), optional :: is_ascii  !< Sentinel to check if file is ASCII.
   integer(I4P)                                 :: f         !< Counter.

   call self%initialize(skip_destroy=.true., file_name=file_name, is_ascii=is_ascii)
   call self%open_file(file_action='write')
   call self%save_header_into_file
   if (self%is_ascii) then
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_ascii(file_unit=self%file_unit)
      enddo
   else
      do f=1, self%facets_number
         call self%facet(f)%save_into_file_binary(file_unit=self%file_unit)
      enddo
   endif
   call self%close_file
   endsubroutine save_into_file

   ! operators
   ! =
   pure subroutine file_stl_assign_file_stl(lhs, rhs)
   !< Operator `=`.
   class(file_stl_object), intent(inout) :: lhs !< Left hand side.
   type(file_stl_object),  intent(in)    :: rhs !< Right hand side.

   if (allocated(rhs%file_name)) then
      lhs%file_name = rhs%file_name
   else
      if (allocated(lhs%file_name))  deallocate(lhs%file_name)
   endif
   lhs%file_unit = rhs%file_unit
   lhs%header = rhs%header
   lhs%facets_number = rhs%facets_number
   if (allocated(rhs%facet)) then
      lhs%facet = rhs%facet
   else
      if (allocated(lhs%facet))  deallocate(lhs%facet)
   endif
   lhs%is_ascii = rhs%is_ascii
   lhs%is_open = rhs%is_open
   endsubroutine file_stl_assign_file_stl

   ! private methods
   elemental subroutine allocate_facets(self)
   !< Allocate facets.
   !<
   !< @note Facets previously allocated are lost.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%facets_number>0) then
      if (allocated(self%facet))  deallocate(self%facet)
      allocate(self%facet(1:self%facets_number))
   endif
   endsubroutine allocate_facets

   subroutine load_facets_number_from_file(self)
   !< Load facets number from file.
   !<
   !< @note File is rewinded.
   class(file_stl_object), intent(inout) :: self         !< File STL.
   character(FRLEN)                      :: facet_record !< Facet record string buffer.

   if (self%is_open) then
      self%facets_number = 0
      rewind(self%file_unit)
      if (self%is_ascii) then
         do
            read(self%file_unit, '(A)', end=10, err=10) facet_record
            if (index(string=facet_record, substring='facet normal') > 0) self%facets_number = self%facets_number + 1
         enddo
      else
         read(self%file_unit, end=10, err=10) facet_record
         read(self%file_unit, end=10, err=10) self%facets_number
      endif
      10 rewind(self%file_unit)
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load facets number from file!'
   endif
   endsubroutine load_facets_number_from_file

   subroutine load_header_from_file(self)
   !< Load header from file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         read(self%file_unit, '(A)') self%header
         self%header = trim(self%header(index(self%header, 'solid')+1:))
      else
         read(self%file_unit) self%header
         read(self%file_unit) self%facets_number
      endif
   else
      write(stderr, '(A)') 'error: file is not open, impossible to load header from file!'
   endif
   endsubroutine load_header_from_file

   subroutine save_header_into_file(self)
   !< Save header into file.
   class(file_stl_object), intent(inout) :: self !< File STL.

   if (self%is_open) then
      rewind(self%file_unit)
      if (self%is_ascii) then
         write(self%file_unit, '(A)') 'solid '//trim(self%header)
      else
         write(self%file_unit) self%header
         write(self%file_unit) self%facets_number
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
