!< FOSSIL test: named patches, I/O (issue #66, step 1).
!<
!< 1. A three-block ASCII STL (wall, inlet, outlet: the faces of the unit cube split 8 / 2 / 2): three patches with
!<    their names and counts, facets in file order; before #66 the reader stopped at the first `endsolid`.
!< 2. Round trip: written back (one block per patch) and read again, the names, the patch of every facet and the
!<    vertex bits are the same.
!< 3. `load_from_files` of the three blocks as three files equals the three-block file (names, counts, vertices).
!< 4. A single-solid ASCII file and a binary file: one patch named after the header; the written ASCII file is the
!<    single block it was before #66 (first line `solid <header>`, last line `endsolid <header>`).
!< 5. `set_patch` creates a named patch and moves facets into it; an empty patch writes no block.
program fossil_test_patches_io

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P

implicit none

character(*), parameter  :: D = 'fossil_test_patches_io'
character(3), parameter  :: SUFFIX(6) = ['-3 ', '-rt', '-a ', '-b ', '-c ', '-1 ']
type(surface_stl_object) :: s, r, m, c
type(facet_object), pointer :: f, g
integer(I4P)             :: k, u, n
logical                  :: ok1, ok2, ok3, ok4, ok5
character(256)           :: first, last, line

call write_cube_blocks(D//'-3.stl')
call s%load_from_file(file_name=D//'-3.stl', guess_format=.true.)
ok1 = s%get_facets_number() == 12 .and. s%get_patches_number() == 3 .and. s%patch_name(1) == 'wall' .and. &
      s%patch_name(2) == 'inlet' .and. s%patch_name(3) == 'outlet' .and. count_patch(s, 1) == 8 .and. &
      count_patch(s, 2) == 2 .and. count_patch(s, 3) == 2 .and. s%facet_patch(9) == 2 .and. s%facet_patch(12) == 3
print '(A,L1)', '1. three-block file: ', ok1

call s%save_into_file(file_name=D//'-rt.stl', is_ascii=.true.)
call r%load_from_file(file_name=D//'-rt.stl', guess_format=.true.)
ok2 = same(s, r)
print '(A,L1)', '2. round trip: ', ok2

call write_cube_blocks(D//'-a.stl', only=1) ; call write_cube_blocks(D//'-b.stl', only=2)
call write_cube_blocks(D//'-c.stl', only=3)
call m%load_from_files([character(len=64) :: D//'-a.stl', D//'-b.stl', D//'-c.stl'])
ok3 = same(s, m)
print '(A,L1)', '3. load_from_files: ', ok3

call c%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
ok4 = c%get_patches_number() == 1 .and. c%patch_name(1) == trim(c%get_header()) .and. all(patches(c) == 1)
call c%save_into_file(file_name=D//'-1.stl', is_ascii=.true.)
open(newunit=u, file=D//'-1.stl', action='read')
read(u, '(A)') first ; n = 1
do
   read(u, '(A)', end=10) line ; last = line ; n = n + 1
enddo
10 close(u)
ok4 = ok4 .and. trim(first) == 'solid '//trim(c%get_header()) .and. trim(last) == 'endsolid '//trim(c%get_header()) .and. &
      n == 2 + 7 * c%get_facets_number()
call r%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.) ! binary
ok4 = ok4 .and. r%get_patches_number() == 1 .and. all(patches(r) == 1)
print '(A,L1)', '4. single solid and binary unchanged: ', ok4

call c%set_patch([1_I4P, 2_I4P], 'inlet')
ok5 = c%get_patches_number() == 2 .and. c%patch_name(2) == 'inlet' .and. c%facet_patch(1) == 2 .and. &
      c%facet_patch(2) == 2 .and. c%facet_patch(3) == 1
call c%set_patch([3_I4P], 'inlet')
ok5 = ok5 .and. c%get_patches_number() == 2 .and. c%facet_patch(3) == 2
! an empty patch (here the header's, once every facet moved out of it) writes no block
call c%set_patch([(k, k = 1, c%get_facets_number())], 'all')
call c%save_into_file(file_name=D//'-1.stl', is_ascii=.true.)
call r%load_from_file(file_name=D//'-1.stl', guess_format=.true.)
ok5 = ok5 .and. r%get_patches_number() == 1 .and. r%patch_name(1) == 'all' .and. r%get_facets_number() == 12
print '(A,L1)', '5. set_patch: ', ok5

do k = 1, 6
   open(newunit=u, file=D//trim(SUFFIX(k))//'.stl')
   close(u, status='delete')
enddo
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4 .and. ok5

contains
   subroutine write_cube_blocks(name, only)
   !< The 12 facets of src/tests/cube.stl as blocks wall (8), inlet (2), outlet (2); or one block alone.
   character(*),           intent(in) :: name !< File.
   integer(I4P), optional, intent(in) :: only !< Write this block only.
   type(surface_stl_object)           :: cube
   type(facet_object), pointer        :: ff
   character(6), parameter            :: NAMES(3) = ['wall  ', 'inlet ', 'outlet']
   integer(I4P), parameter            :: FIRST(3) = [1, 9, 11], LAST(3) = [8, 10, 12]
   integer(I4P)                       :: uu, b, i, j

   call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   open(newunit=uu, file=name, action='write', status='replace')
   do b = 1, 3
      if (present(only)) then
         if (b /= only) cycle
      endif
      write(uu, '(A)') 'solid '//trim(NAMES(b))
      do i = FIRST(b), LAST(b)
         ff => cube%facet_at(i)
         write(uu, '(A,3ES26.17)') '  facet normal', ff%normal%x, ff%normal%y, ff%normal%z
         write(uu, '(A)') '    outer loop'
         do j = 1, 3
            write(uu, '(A,3ES26.17)') '      vertex', ff%vertex(j)%x, ff%vertex(j)%y, ff%vertex(j)%z
         enddo
         write(uu, '(A)') '    endloop'
         write(uu, '(A)') '  endfacet'
      enddo
      write(uu, '(A)') 'endsolid '//trim(NAMES(b))
   enddo
   close(uu)
   endsubroutine write_cube_blocks

   function count_patch(x, p) result(n)
   type(surface_stl_object), intent(in) :: x
   integer(I4P),             intent(in) :: p
   integer(I4P)                         :: n

   n = count(patches(x) == p)
   endfunction count_patch

   function patches(x) result(p)
   type(surface_stl_object), intent(in) :: x
   integer(I4P), allocatable            :: p(:)
   integer(I4P)                         :: i

   allocate(p(x%get_facets_number()))
   do i = 1, size(p)
      p(i) = x%facet_patch(i)
   enddo
   endfunction patches

   function same(x, y) result(eq)
   !< Same patch names, patch of every facet and vertex bits.
   type(surface_stl_object), intent(in) :: x, y
   logical                              :: eq
   integer(I4P)                         :: i, j

   eq = x%get_facets_number() == y%get_facets_number() .and. x%get_patches_number() == y%get_patches_number()
   if (.not. eq) return
   do i = 1, x%get_patches_number()
      eq = eq .and. x%patch_name(i) == y%patch_name(i)
   enddo
   do i = 1, x%get_facets_number()
      f => x%facet_at(i) ; g => y%facet_at(i)
      eq = eq .and. f%patch == g%patch
      do j = 1, 3
         eq = eq .and. f%vertex(j)%x == g%vertex(j)%x .and. f%vertex(j)%y == g%vertex(j)%y .and. &
                       f%vertex(j)%z == g%vertex(j)%z
      enddo
   enddo
   endfunction same
endprogram fossil_test_patches_io
