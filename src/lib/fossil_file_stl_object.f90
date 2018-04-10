!< FOSSIL,  STL file class definition.

module fossil_file_stl_object
!< FOSSIL,  STL file class definition.

! use, intrinsic :: iso_fortran_env, only : stdout => output_unit

implicit none
private
public :: file_stl_object

type :: file_stl_object
  !< FOSSIL STL file class.
  private
  character(len=:), allocatable, public :: file_name             !< File name
  ! integer(I4P)                          :: Ns = 0                !< Number of sections.
  ! character(1)                          :: opt_sep = DEF_OPT_SEP !< Separator character of option name/value.
  ! type(section), allocatable            :: sections(:)           !< Sections.
  ! contains
    ! public methods
    ! generic               :: add          => add_section, &             !< Add a section.
    !                                          add_option,  &             !< Add an option to a section (scalar).
    !                                          add_a_option               !< Add an option to a section (array).
    ! procedure, pass(self) :: count_values                               !< Count option value(s).
    ! generic               :: del          => free_option_of_section, &  !< Remove (freeing) an option of a section.
    !                                          free_section               !< Remove (freeing) a section.
    ! procedure, pass(self) :: free                                       !< Free dynamic memory destroyng file data.
    ! generic               :: free_options => free_options_all,        & !< Free all options.
    !                                          free_options_of_section, & !< Free all options of a section.
    !                                          free_option_of_section     !< Free an option of a section.
    ! generic               :: get          => get_option, &              !< Get option value (scalar).
    !                                          get_a_option               !< Get option value (array).
    ! procedure, pass(self) :: get_items                                  !< Get list of pairs option name/value.
    ! procedure, pass(self) :: get_sections_list                          !< Get sections names list.
    ! procedure, pass(self) :: initialize                                 !< Initialize file.
    ! procedure, pass(self) :: has_option                                 !< Inquire the presence of an option.
    ! procedure, pass(self) :: has_section                                !< Inquire the presence of a section.
    ! generic               :: index        => index_section, &           !< Return the index of a section.
    !                                          index_option               !< Return the index of an option.
    ! procedure, pass(self) :: load                                       !< Load file data.
    ! generic               :: loop         => loop_options_section, &    !< Loop over options of a section.
    !                                          loop_options               !< Loop over all options.
    ! procedure, pass(self) :: print        => print_file_stl_object             !< Pretty printing data.
    ! procedure, pass(self) :: save         => save_file_stl_object              !< Save data.
    ! procedure, pass(self) :: section      => section_file_stl_object           !< Get section name once provided an index.
    ! ! operators overloading
    ! generic :: assignment(=) => assign_file_stl_object !< Procedure for section assignment overloading.
    ! ! private methods
    ! procedure, private, pass(self) :: add_a_option            !< Add an option to a section (array).
    ! procedure, private, pass(self) :: add_option              !< Add an option to a section (scalar).
    ! procedure, private, pass(self) :: add_section             !< Add a section.
    ! procedure, private, pass(self) :: free_options_all        !< Free all options of all sections.
    ! procedure, private, pass(self) :: free_options_of_section !< Free all options of a section.
    ! procedure, private, pass(self) :: free_option_of_section  !< Free an option of a section.
    ! procedure, private, pass(self) :: free_section            !< Free a section.
    ! procedure, private, pass(self) :: get_a_option            !< Get option value (array).
    ! procedure, private, pass(self) :: get_option              !< Get option value (scalar).
    ! procedure, private, pass(self) :: index_option            !< Return the index of an option.
    ! procedure, private, pass(self) :: index_section           !< Return the index of a section.
    ! procedure, private, pass(self) :: loop_options            !< Loop over all options.
    ! procedure, private, pass(self) :: loop_options_section    !< Loop over options of a section.
    ! procedure, private, pass(self) :: parse                   !< Parse file data.
    ! ! assignments
    ! procedure, private, pass(lhs) :: assign_file_stl_object !< Assignment overloading.
endtype file_stl_object

contains
  ! public methods

  ! private methods
endmodule fossil_file_stl_object
