!run tut01-patches tut01_patches
!fig tut01-patches
!< Tutorial 1, third run: named patches. The dragon's facets are split into front, middle and rear by the x of their
!< centroids, saved as an ASCII STL with one solid block per patch, and read back.
program tut01_patches
use fossil_fig,          only : fig_open, fig_triangle, fig_close
!region name
use fossil,              only : surface_stl_object
use fossil_facet_object, only : facet_object
use penf,                only : I4P, R8P

implicit none
type(surface_stl_object)    :: surface, back
type(facet_object), pointer :: f
integer(I4P), allocatable   :: ids(:)
real(R8P),    allocatable   :: xc(:)
integer(I4P)                :: i, p, n, u

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)
n = surface%get_facets_number()
allocate(xc(n))
do i = 1, n
   f => surface%facet_at(i) ; xc(i) = f%centroid%x
enddo
ids = [(i, i = 1, n)]
call surface%set_patch(facets=ids, name='middle')                        ! every facet, then two subsets
call surface%set_patch(facets=pack(ids, xc < -4._R8P), name='front')
call surface%set_patch(facets=pack(ids, xc >  4._R8P), name='rear')
call surface%save_into_file(file_name='dragon-patches.stl', is_ascii=.true.) ! one solid block per patch
!endregion name
!region read
call back%load_from_file(file_name='dragon-patches.stl', guess_format=.true.)
print '(A,I0)', 'patches read back ', back%get_patches_number()
do p = 1, back%get_patches_number()
   print '(A,I0,A,A8,I6,A)', 'patch ', p, '  ', back%patch_name(p), &
                             count(back%facet_patch(ids) == p), ' facets'
enddo
!endregion read
call fig_open(u, 'tut01-patches', 'dragon-patches.stl read back: patches front (teal), middle (cream), rear (ember)', &
              -60._R8P, 22._R8P)
do i = 1, n
   f => back%facet_at(i)
   call fig_triangle(u, 'field', f%vertex(1), f%vertex(2), f%vertex(3), &
                     merge(-0.8_R8P, merge(0.8_R8P, 0._R8P, back%patch_name(back%facet_patch(i)) == 'rear'), &
                           back%patch_name(back%facet_patch(i)) == 'front'))
enddo
call fig_close(u)
endprogram tut01_patches
