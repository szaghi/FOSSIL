!< FOSSIL, list of IDs class definition.

module fossil_list_id_object
!< FOSSIL, list of IDs class definition.

use penf, only : I4P

implicit none
private
public :: list_id_object

type :: list_id_object
   !< FOSSIL, list of IDs class.
   integer(I4P)              :: ids_number=0 !< Number of IDs.
   integer(I4P), allocatable :: id(:)        !< IDs list.
   contains
      ! public methods
      procedure, pass(self) :: del        !< Delete ID from list.
      procedure, pass(self) :: destroy    !< Destroy list.
      procedure, pass(self) :: initialize !< Initialize list.
      procedure, pass(self) :: put        !< Put ID in list.
      ! operators
      generic :: assignment(=) => list_id_assign_list_id !< Overload `=`.
      ! private methods
      procedure, pass(lhs),  private :: list_id_assign_list_id !< Operator `=`.
endtype list_id_object

contains
   ! public methods
   elemental subroutine del(self, id)
   !< Delete ID from list.
   class(list_id_object), intent(inout) :: self      !< List.
   integer(I4P),          intent(in)    :: id        !< Given ID.
   integer(I4P), allocatable            :: id_tmp(:) !< Temporary list.
   integer(I4P)                         :: i, ii     !< Counter.

   if (self%ids_number>0) then
      if (any(self%id==id)) then
         allocate(id_tmp(1:self%ids_number-1))
         ii = 0
         do i=1, self%ids_number
            if (self%id(i) /= id) then
               ii = ii + 1
               id_tmp(ii) = self%id(i)
            endif
         enddo
         call move_alloc(from=id_tmp, to=self%id)
         self%ids_number = self%ids_number - 1
      endif
   endif
   endsubroutine del

   elemental subroutine destroy(self)
   !< Destroy list.
   class(list_id_object), intent(inout) :: self  !< List.
   type(list_id_object)                 :: fresh !< Fresh instance of list.

   self = fresh
   endsubroutine destroy

   pure subroutine initialize(self, id)
   !< Initialize list.
   class(list_id_object), intent(inout)        :: self  !< List.
   integer(I4P),          intent(in), optional :: id(:) !< IDs list.

   call self%destroy
   if (present(id)) then
      self%ids_number = size(id, dim=1)
      allocate(self%id, source=id)
   endif
   endsubroutine initialize

   elemental subroutine put(self, id)
   !< Put ID in list.
   !<
   !< ID is put in list only if it is not already present.
   class(list_id_object), intent(inout) :: self      !< List.
   integer(I4P),          intent(in)    :: id        !< Given ID.
   integer(I4P), allocatable            :: id_tmp(:) !< Temporary list.

   if (self%ids_number>0) then
      if (.not.any(self%id==id)) then
         allocate(id_tmp(1:self%ids_number+1))
         id_tmp(1:self%ids_number) = self%id(1:self%ids_number)
         id_tmp(self%ids_number+1) = id
         call move_alloc(from=id_tmp, to=self%id)
         self%ids_number = self%ids_number + 1
      endif
   else
      allocate(self%id(1))
      self%id(1) = id
      self%ids_number = 1
   endif
   endsubroutine put

   ! `=` operator
   pure subroutine list_id_assign_list_id(lhs, rhs)
   !< Operator `=`.
   class(list_id_object), intent(inout) :: lhs !< Left hand side.
   type(list_id_object),  intent(in)    :: rhs !< Right hand side.

   lhs%ids_number = rhs%ids_number
   if (allocated(lhs%id)) deallocate(lhs%id)
   if (allocated(rhs%id)) lhs%id = rhs%id
   endsubroutine list_id_assign_list_id
endmodule fossil_list_id_object
