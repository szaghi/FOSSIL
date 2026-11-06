!< FOSSIL test: named patches through the pipeline (issue #66, step 2).
!<
!< Body: cube-02.stl (the unit cube, 192 facets), written as three blocks by the normal of each facet: z- face
!< `inlet`, z+ face `outlet`, the four side faces `wall` (areas 1, 1, 4).
!< 1. sanitize: a degenerate facet in `outlet` and a duplicate of a `wall` facet in `inlet` are removed; the area of
!<    every patch is unchanged (a duplicate in another patch: the first one wins).
!< 2. Rigid motions (rotate, translate, mirror) keep the area of every patch.
!< 3. clip keeps the patch of every facet, in the clipped part and in the remainder (which keeps the names).
!< 4. merge_solids with a translated copy whose patches are `wall` and `heated`: the tables merge by name (wall,
!<    inlet, outlet, heated), the areas add up per name.
!< 5. fill_holes: a hole inside the wall gets the patch `wall`; a hole on the edge between the wall and the inlet
!<    gets the patch `filled`; wall + filled = 4, inlet and outlet 1 each, total 6.
program fossil_test_patches_pipeline

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

character(*), parameter  :: D = 'fossil_test_patches_pipeline'
type(surface_stl_object) :: s, r, t, c
type(facet_object), pointer :: f
integer(I4P), allocatable :: block(:)
integer(I4P)             :: i, u, interior, edge
logical                  :: ok1, ok2, ok3, ok4, ok5
real(R8P)                :: a(4)

call blocks_of_cube(block, interior, edge)

! ---- 1. sanitize
call write_blocks(D//'-dirty.stl', block, extra=.true.)
call s%load_from_file(file_name=D//'-dirty.stl', guess_format=.true.)
call s%sanitize
ok1 = s%get_patches_number() == 3 .and. s%get_facets_number() == 192 .and. &
      close3(areas(s, 3), [4._R8P, 1._R8P, 1._R8P])
print '(A,L1,A,3F10.6)', '1. sanitize: ', ok1, '  areas ', areas(s, 3)

! ---- 2. rigid motions
r = s
call r%rotate(axis=vector_R8P(1._R8P, 2._R8P, 3._R8P), angle=0.7_R8P)
call r%translate(x=0.3_R8P, y=-2._R8P, z=1._R8P)
call r%mirror(normal=vector_R8P(0._R8P, 1._R8P, 0._R8P))
call r%analyze
ok2 = close3(areas(r, 3), [4._R8P, 1._R8P, 1._R8P])
print '(A,L1)', '2. rigid motions: ', ok2

! ---- 3. clip
r = s
call r%clip(bmin=vector_R8P(-1._R8P, -1._R8P, -1._R8P), bmax=vector_R8P(0.3_R8P, 2._R8P, 2._R8P), remainder=t)
ok3 = r%get_patches_number() == 3 .and. t%get_patches_number() == 3 .and. t%patch_name(2) == 'inlet' .and. &
      close3(areas(r, 3) + areas(t, 3), [4._R8P, 1._R8P, 1._R8P])
print '(A,L1)', '3. clip: ', ok3

! ---- 4. merge_solids
r = s
c = s
call c%translate(x=3._R8P) ; call c%analyze
call c%set_patch([(i, i = 1, c%get_facets_number())], 'heated')
call c%set_patch([1_I4P, 2_I4P, 3_I4P, 4_I4P], 'wall')
call r%merge_solids(c)
ok4 = r%get_patches_number() == 4 .and. r%patch_name(4) == 'heated' .and. &
      close3(areas(r, 3), [4._R8P + area_of(c, 1), 1._R8P, 1._R8P]) .and. abs(areas4(r) - (6._R8P - area_of(c, 1))) < 1e-12_R8P
print '(A,L1)', '4. merge_solids: ', ok4

! ---- 5. fill_holes
call write_blocks(D//'-holes.stl', block, skip=[interior, edge])
call r%load_from_file(file_name=D//'-holes.stl', guess_format=.true.)
call r%fill_holes(max_hole_size=0.5_R8P)
a = 0._R8P
do i = 1, r%get_facets_number()
   f => r%facet_at(i)
   a(f%patch) = a(f%patch) + f%area()
enddo
ok5 = r%get_patches_number() == 4 .and. r%patch_name(4) == 'filled' .and. &
      abs(a(1) + a(4) - 4._R8P) < 1e-12_R8P .and. a(4) > 0._R8P .and. abs(sum(a) - 6._R8P) < 1e-12_R8P .and. &
      abs(a(2) - 1._R8P) < 1e-12_R8P .and. abs(a(3) - 1._R8P) < 1e-12_R8P
print '(A,L1,A,4F10.6)', '5. fill_holes: ', ok5, '  areas ', a

do i = 1, 2
   open(newunit=u, file=D//trim(merge('-dirty', '-holes', i == 1))//'.stl')
   close(u, status="delete")
enddo
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4 .and. ok5

contains
   subroutine blocks_of_cube(b, fin, fedge)
   !< Block of every facet of cube-02 (1 wall, 2 inlet z-, 3 outlet z+); a wall facet strictly inside the face x = 0,
   !< and a wall facet of x = 0 with an edge on z = 0 (next to the inlet).
   integer(I4P), allocatable, intent(out) :: b(:)
   integer(I4P),              intent(out) :: fin, fedge
   type(surface_stl_object)               :: cube
   type(facet_object), pointer            :: ff
   integer(I4P)                           :: k, nz

   call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   allocate(b(cube%get_facets_number()))
   fin = 0 ; fedge = 0
   do k = 1, size(b)
      ff => cube%facet_at(k)
      b(k) = 1
      if (ff%normal%z < -0.5_R8P) b(k) = 2
      if (ff%normal%z > 0.5_R8P) b(k) = 3
      if (ff%normal%x < -0.5_R8P) then
         nz = count([ff%vertex(1)%z, ff%vertex(2)%z, ff%vertex(3)%z] == 0._R8P)
         if (fin == 0 .and. all([ff%vertex(1)%y, ff%vertex(2)%y, ff%vertex(3)%y, ff%vertex(1)%z, ff%vertex(2)%z, &
             ff%vertex(3)%z] > 0._R8P) .and. all([ff%vertex(1)%y, ff%vertex(2)%y, ff%vertex(3)%y, ff%vertex(1)%z, &
             ff%vertex(2)%z, ff%vertex(3)%z] < 1._R8P)) fin = k
         if (fedge == 0 .and. nz == 2) fedge = k
      endif
   enddo
   endsubroutine blocks_of_cube

   subroutine write_blocks(name, b, extra, skip)
   !< cube-02 as blocks wall, inlet, outlet (file order kept inside each block); `extra`: a degenerate facet in
   !< outlet and a duplicate of the first wall facet in inlet; `skip`: facets left out.
   character(*),           intent(in) :: name
   integer(I4P),           intent(in) :: b(:)
   logical,      optional, intent(in) :: extra
   integer(I4P), optional, intent(in) :: skip(:)
   type(surface_stl_object)           :: cube
   type(facet_object), pointer        :: ff
   character(6), parameter            :: NAMES(3) = ['wall  ', 'inlet ', 'outlet']
   integer(I4P)                       :: uu, k, j, first_wall

   call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   first_wall = findloc(b, 1, dim=1)
   open(newunit=uu, file=name, action='write', status='replace')
   do j = 1, 3
      write(uu, '(A)') 'solid '//trim(NAMES(j))
      do k = 1, size(b)
         if (b(k) /= j) cycle
         if (present(skip)) then
            if (any(skip == k)) cycle
         endif
         ff => cube%facet_at(k)
         call put(uu, ff%vertex(1), ff%vertex(2), ff%vertex(3))
      enddo
      if (present(extra)) then
         if (j == 2) then
            ff => cube%facet_at(first_wall)
            call put(uu, ff%vertex(1), ff%vertex(2), ff%vertex(3))
         endif
         if (j == 3) call put(uu, vector_R8P(0._R8P, 0._R8P, 1._R8P), vector_R8P(0.5_R8P, 0._R8P, 1._R8P), &
                              vector_R8P(1._R8P, 0._R8P, 1._R8P))
      endif
      write(uu, '(A)') 'endsolid '//trim(NAMES(j))
   enddo
   close(uu)
   endsubroutine write_blocks

   subroutine put(uu, v1, v2, v3)
   !< One ASCII facet.
   integer(I4P),     intent(in) :: uu
   type(vector_R8P), intent(in) :: v1, v2, v3

   write(uu, '(A)') '  facet normal 0 0 0'
   write(uu, '(A)') '    outer loop'
   write(uu, '(A,3ES26.17)') '      vertex', v1%x, v1%y, v1%z
   write(uu, '(A,3ES26.17)') '      vertex', v2%x, v2%y, v2%z
   write(uu, '(A,3ES26.17)') '      vertex', v3%x, v3%y, v3%z
   write(uu, '(A)') '    endloop'
   write(uu, '(A)') '  endfacet'
   endsubroutine put

   function areas(x, n) result(ar)
   type(surface_stl_object), intent(in) :: x
   integer(I4P),             intent(in) :: n
   real(R8P)                            :: ar(n)
   type(facet_object), pointer          :: ff
   integer(I4P)                         :: k

   ar = 0._R8P
   do k = 1, x%get_facets_number()
      ff => x%facet_at(k)
      if (ff%patch <= n) ar(ff%patch) = ar(ff%patch) + ff%area()
   enddo
   endfunction areas

   function areas4(x) result(a4)
   type(surface_stl_object), intent(in) :: x
   real(R8P)                            :: a4, ar(4)
   type(facet_object), pointer          :: ff
   integer(I4P)                         :: k

   ar = 0._R8P
   do k = 1, x%get_facets_number()
      ff => x%facet_at(k)
      ar(ff%patch) = ar(ff%patch) + ff%area()
   enddo
   a4 = ar(4)
   endfunction areas4

   function area_of(x, p) result(ap)
   type(surface_stl_object), intent(in) :: x
   integer(I4P),             intent(in) :: p
   real(R8P)                            :: ap
   type(facet_object), pointer          :: ff
   integer(I4P)                         :: k

   ap = 0._R8P
   do k = 1, x%get_facets_number()
      ff => x%facet_at(k)
      if (ff%patch == p) ap = ap + ff%area()
   enddo
   endfunction area_of

   pure function close3(x, y) result(ok)
   real(R8P), intent(in) :: x(3), y(3)
   logical               :: ok

   ok = all(abs(x - y) < 1e-12_R8P)
   endfunction close3
endprogram fossil_test_patches_pipeline
