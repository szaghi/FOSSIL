!< FOSSIL, facet class definition.

module fossil_facet_object
!< FOSSIL, facet class definition.

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I2P, I4P, R4P, str

implicit none
private
public :: facet_object
public :: FRLEN

integer(I4P), parameter :: FRLEN=80 !< Maximum length of facet record string.

type :: facet_object
   !< FOSSIL, facet class.
   real(R4P) :: normal(1:3)=[0._R4P,0._R4P,0._R4P]   !< Facet (outward) normal (versor).
   real(R4P) :: vertex_1(1:3)=[0._R4P,0._R4P,0._R4P] !< Facet vertex 1.
   real(R4P) :: vertex_2(1:3)=[0._R4P,0._R4P,0._R4P] !< Facet vertex 2.
   real(R4P) :: vertex_3(1:3)=[0._R4P,0._R4P,0._R4P] !< Facet vertex 3.
   contains
      ! public methods
      procedure, pass(self) :: initialize            !< Initialize facet.
      procedure, pass(self) :: load_from_file_ascii  !< Load facet from ASCII file.
      procedure, pass(self) :: load_from_file_binary !< Load facet from binary file.
      procedure, pass(self) :: save_into_file_ascii  !< Save facet into ASCII file.
      procedure, pass(self) :: save_into_file_binary !< Save facet into binary file.
      ! operators
      generic :: assignment(=) => facet_assign_facet !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: facet_assign_facet !< Operator `=`.
endtype facet_object

contains
   ! public methods
   elemental subroutine initialize(self)
   !< Initialize facet.
   class(facet_object), intent(inout) :: self  !< Facet.
   type(facet_object)                 :: fresh !< Fresh instance of facet.

   self = fresh
   endsubroutine initialize

   subroutine load_from_file_ascii(self, file_unit)
   !< Load facet from ASCII file.
   class(facet_object), intent(inout) :: self      !< Facet.
   integer(I4P),        intent(in)    :: file_unit !< File unit.

   call load_facet_record(prefix='facet normal', record=self%normal)
   read(file_unit, *) ! outer loop
   call load_facet_record(prefix='vertex', record=self%vertex_1)
   call load_facet_record(prefix='vertex', record=self%vertex_2)
   call load_facet_record(prefix='vertex', record=self%vertex_3)
   read(file_unit, *) ! endloop
   read(file_unit, *) ! endfacet
   contains
      subroutine load_facet_record(prefix, record)
      !< Load a facet *record*, namely normal or vertex data.
      character(*), intent(in)  :: prefix       !< Record prefix string.
      real(R4P),    intent(out) :: record(1:3)  !< Record data.
      character(FRLEN)          :: facet_record !< Facet record string buffer.
      integer(I4P)              :: i            !< Counter.

      read(file_unit, '(A)') facet_record
      i = index(string=facet_record, substring=prefix)
      if (i>0) then
         read(facet_record(i+len(prefix):), *) record(1), record(2), record(3)
      else
         write(stderr, '(A)') 'error: impossible to read "'//prefix//'" from file unit "'//trim(str(file_unit))//'"!'
      endif
      endsubroutine load_facet_record
   endsubroutine load_from_file_ascii

   subroutine load_from_file_binary(self, file_unit)
   !< Load facet from binary file.
   class(facet_object), intent(inout) :: self      !< Facet.
   integer(I4P),        intent(in)    :: file_unit !< File unit.
   integer(I2P)                       :: padding   !< Facet padding.

   read(file_unit) self%normal(1), self%normal(2), self%normal(3)
   read(file_unit) self%vertex_1(1), self%vertex_1(2), self%vertex_1(3)
   read(file_unit) self%vertex_2(1), self%vertex_2(2), self%vertex_2(3)
   read(file_unit) self%vertex_3(1), self%vertex_3(2), self%vertex_3(3)
   read(file_unit) padding
   endsubroutine load_from_file_binary

   subroutine save_into_file_ascii(self, file_unit)
   !< Save facet into ASCII file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit, '(A)') '  facet normal '//trim(str(self%normal(1)))//' '//&
                                              trim(str(self%normal(2)))//' '//&
                                              trim(str(self%normal(3)))
   write(file_unit, '(A)') '    outer loop'
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_1(1)))//' '//&
                                            trim(str(self%vertex_1(2)))//' '//&
                                            trim(str(self%vertex_1(3)))
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_2(1)))//' '//&
                                            trim(str(self%vertex_2(2)))//' '//&
                                            trim(str(self%vertex_2(3)))
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_3(1)))//' '//&
                                            trim(str(self%vertex_3(2)))//' '//&
                                            trim(str(self%vertex_3(3)))
   write(file_unit, '(A)') '    endloop'
   write(file_unit, '(A)') '  endfacet'
   endsubroutine save_into_file_ascii

   subroutine save_into_file_binary(self, file_unit)
   !< Save facet into binary file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit) self%normal(1), self%normal(2), self%normal(3)
   write(file_unit) self%vertex_1(1), self%vertex_1(2), self%vertex_1(3)
   write(file_unit) self%vertex_2(1), self%vertex_2(2), self%vertex_2(3)
   write(file_unit) self%vertex_3(1), self%vertex_3(2), self%vertex_3(3)
   write(file_unit) 0_I2P
   endsubroutine save_into_file_binary

   ! private methods
   ! `=` operator
   pure subroutine facet_assign_facet(lhs, rhs)
   !< Operator `=`.
   class(facet_object), intent(inout) :: lhs !< Left hand side.
   type(facet_object),  intent(in)    :: rhs !< Right hand side.

   lhs%normal = rhs%normal
   lhs%vertex_1 = rhs%vertex_1
   lhs%vertex_2 = rhs%vertex_2
   lhs%vertex_3 = rhs%vertex_3
   endsubroutine facet_assign_facet
endmodule fossil_facet_object
