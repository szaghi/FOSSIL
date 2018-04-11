!< FOSSIL, facet class definition.

module fossil_facet_object
!< FOSSIL, facet class definition.

use, intrinsic :: iso_fortran_env, only : stderr => error_unit
use penf, only : I2P, I4P, R4P, str, ZeroR4P
use vecfor, only : face_normal3_R4P, vector_R4P

implicit none
private
public :: facet_object
public :: FRLEN

integer(I4P), parameter :: FRLEN=80 !< Maximum length of facet record string.

type :: facet_object
   !< FOSSIL, facet class.
   type(vector_R4P) :: normal   !< Facet (outward) normal (versor).
   type(vector_R4P) :: vertex_1 !< Facet vertex 1.
   type(vector_R4P) :: vertex_2 !< Facet vertex 2.
   type(vector_R4P) :: vertex_3 !< Facet vertex 3.
   contains
      ! public methods
      procedure, pass(self) :: check_normal          !< Check normal consistency.
      procedure, pass(self) :: initialize            !< Initialize facet.
      procedure, pass(self) :: load_from_file_ascii  !< Load facet from ASCII file.
      procedure, pass(self) :: load_from_file_binary !< Load facet from binary file.
      procedure, pass(self) :: sanitize_normal       !< Sanitize normal, make normal consistent with vertices.
      procedure, pass(self) :: save_into_file_ascii  !< Save facet into ASCII file.
      procedure, pass(self) :: save_into_file_binary !< Save facet into binary file.
      ! operators
      generic :: assignment(=) => facet_assign_facet !< Overload `=`.
      ! private methods
      procedure, pass(lhs)  :: facet_assign_facet !< Operator `=`.
endtype facet_object

contains
   ! public methods
   elemental function check_normal(self) result(is_consistent)
   !< Check normal consistency.
   class(facet_object), intent(in) :: self          !< Facet.
   logical                         :: is_consistent !< Consistency check result.
   type(vector_R4P)                :: normal        !< Normal computed by means of vertices data.

   normal = face_normal3_R4P(pt1=self%vertex_1, pt2=self%vertex_2, pt3=self%vertex_3, norm='y')
   is_consistent = ((abs(normal%x - self%normal%x)<=2*ZeroR4P).and.&
                    (abs(normal%y - self%normal%y)<=2*ZeroR4P).and.&
                    (abs(normal%z - self%normal%z)<=2*ZeroR4P))
   endfunction check_normal

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
      character(*),     intent(in)  :: prefix       !< Record prefix string.
      type(vector_R4P), intent(out) :: record       !< Record data.
      character(FRLEN)              :: facet_record !< Facet record string buffer.
      integer(I4P)                  :: i            !< Counter.

      read(file_unit, '(A)') facet_record
      i = index(string=facet_record, substring=prefix)
      if (i>0) then
         read(facet_record(i+len(prefix):), *) record%x, record%y, record%z
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

   read(file_unit) self%normal%x, self%normal%y, self%normal%z
   read(file_unit) self%vertex_1%x, self%vertex_1%y, self%vertex_1%z
   read(file_unit) self%vertex_2%x, self%vertex_2%y, self%vertex_2%z
   read(file_unit) self%vertex_3%x, self%vertex_3%y, self%vertex_3%z
   read(file_unit) padding
   endsubroutine load_from_file_binary

   elemental subroutine sanitize_normal(self)
   !< Sanitize normal, make normal consistent with vertices.
   !<
   !<```fortran
   !< type(facet_object) :: facet
   !< facet%vertex_1 = -0.231369_R4P * ex_R4P + 0.0226865_R4P * ey_R4P + 1._R4P * ez_R4P
   !< facet%vertex_2 = -0.227740_R4P * ex_R4P + 0.0245457_R4P * ey_R4P + 0._R4P * ez_R4P
   !< facet%vertex_2 = -0.235254_R4P * ex_R4P + 0.0201881_R4P * ey_R4P + 0._R4P * ez_R4P
   !< call facet%sanitize_normal
   !< print "(3(F3.1,1X))", facet%normal%x, facet%normal%y, facet%normal%z
   !<```
   !=> -0.501673222 0.865057290 -2.12257713<<<
   class(facet_object), intent(inout) :: self !< Facet.

   self%normal = face_normal3_R4P(pt1=self%vertex_1, pt2=self%vertex_2, pt3=self%vertex_3, norm='y')
   endsubroutine sanitize_normal

   subroutine save_into_file_ascii(self, file_unit)
   !< Save facet into ASCII file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit, '(A)') '  facet normal '//trim(str(self%normal%x))//' '//&
                                              trim(str(self%normal%y))//' '//&
                                              trim(str(self%normal%z))
   write(file_unit, '(A)') '    outer loop'
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_1%x))//' '//&
                                            trim(str(self%vertex_1%y))//' '//&
                                            trim(str(self%vertex_1%z))
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_2%x))//' '//&
                                            trim(str(self%vertex_2%y))//' '//&
                                            trim(str(self%vertex_2%z))
   write(file_unit, '(A)') '      vertex '//trim(str(self%vertex_3%x))//' '//&
                                            trim(str(self%vertex_3%y))//' '//&
                                            trim(str(self%vertex_3%z))
   write(file_unit, '(A)') '    endloop'
   write(file_unit, '(A)') '  endfacet'
   endsubroutine save_into_file_ascii

   subroutine save_into_file_binary(self, file_unit)
   !< Save facet into binary file.
   class(facet_object), intent(in) :: self      !< Facet.
   integer(I4P),        intent(in) :: file_unit !< File unit.

   write(file_unit) self%normal%x, self%normal%y, self%normal%z
   write(file_unit) self%vertex_1%x, self%vertex_1%y, self%vertex_1%z
   write(file_unit) self%vertex_2%x, self%vertex_2%y, self%vertex_2%z
   write(file_unit) self%vertex_3%x, self%vertex_3%y, self%vertex_3%z
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
