program ck_patches
use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)    :: duct
type(facet_object), pointer :: f
integer(I4P), allocatable   :: inlet(:), outlet(:)
type(vector_R8P)            :: probes(3)
real(R8P)                   :: d(3)
integer(I4P)                :: i, patch(3)

! tag the faces of a box: z- is the inlet, z+ the outlet, the rest is wall
call duct%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(inlet(0), outlet(0))
do i = 1, duct%get_facets_number()
   f => duct%facet_at(i)
   if (f%normal%z < -0.5_R8P) inlet = [inlet, i]
   if (f%normal%z > 0.5_R8P) outlet = [outlet, i]
enddo
call duct%set_patch(facets=[(i, i = 1, duct%get_facets_number())], name='wall')
call duct%set_patch(facets=inlet, name='inlet')
call duct%set_patch(facets=outlet, name='outlet')
call duct%save_into_file(file_name='duct.stl', is_ascii=.true.) ! one solid block per patch

! read it back: one patch per solid block
call duct%load_from_file(file_name='duct.stl', guess_format=.true.)
do i = 1, duct%get_patches_number()
   print '(A,I0,3A,I0,A)', 'patch ', i, ': ', duct%patch_name(i), ', ', count(facets_of(i)), ' facets'
enddo

! the boundary condition at a point is the patch of its closest facet
probes = [vector_R8P(0.5_R8P, 0.5_R8P, -0.2_R8P), vector_R8P(0.5_R8P, 0.5_R8P, 1.3_R8P), &
          vector_R8P(1.4_R8P, 0.5_R8P, 0.5_R8P)]
call duct%query_many(points=probes, distances=d, patches=patch)
do i = 1, 3
   print '(A,3F5.1,2A)', 'point ', probes(i)%x, probes(i)%y, probes(i)%z, ' -> ', duct%patch_name(patch(i))
enddo
contains
   function facets_of(p) result(m)
   integer(I4P), intent(in) :: p
   logical, allocatable     :: m(:)
   integer(I4P)             :: k

   allocate(m(duct%get_facets_number()))
   do k = 1, size(m)
      m(k) = duct%facet_patch(k) == p
   enddo
   endfunction facets_of
endprogram ck_patches
