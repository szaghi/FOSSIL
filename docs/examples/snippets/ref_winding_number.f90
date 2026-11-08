program ref_winding_number
use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object)                :: cube, holed, dragon
type(facet_object), pointer             :: facets(:)
type(facet_object), allocatable, target :: kept(:)
type(vector_R8P)                        :: q(5), bmin, bmax
character(8),                 parameter :: WHERE(5) = ['centre  ', 'face    ', 'edge    ', 'corner  ', 'outside ']
real(R8P)                               :: err(2)
integer(I4P)                            :: i, j, k
call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
q = [vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P), vector_R8P(0.5_R8P, 0.5_R8P, 0._R8P), vector_R8P(0.5_R8P, 0._R8P, 0._R8P), &
     vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(2._R8P, 2._R8P, 2._R8P)]
do i = 1, 5
   print '(A,A,F12.9)', WHERE(i), ' w = ', cube%winding_number(point=q(i))
enddo
facets => cube%facets_ref()                                   ! drop the 8 facets of the top face z = 1
allocate(kept(count(facets%normal%z < 0.5_R8P)))
j = 0
do i = 1, size(facets)
   if (facets(i)%normal%z > 0.5_R8P) cycle
   j = j + 1 ; kept(j) = facets(i)
enddo
call holed%adopt_facets(facets=kept)
print '(A,L1,A,F6.4,A,F6.4)', 'open box: watertight ', holed%is_watertight(), ', w at the centre ', &
      holed%winding_number(point=q(1)), ', below the opening ', holed%winding_number(point=vector_R8P(0.5_R8P, &
      0.5_R8P, 0.9_R8P))
call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = dragon%get_bmin() ; bmax = dragon%get_bmax()
err = 0._R8P
do k = 1, 5 ; do j = 1, 5 ; do i = 1, 5                       ! 125 points over the dragon's box
   q(1) = bmin + ((i - 0.5_R8P) / 5) * (bmax%x - bmin%x) * vector_R8P(1._R8P, 0._R8P, 0._R8P) + &
                 ((j - 0.5_R8P) / 5) * (bmax%y - bmin%y) * vector_R8P(0._R8P, 1._R8P, 0._R8P) + &
                 ((k - 0.5_R8P) / 5) * (bmax%z - bmin%z) * vector_R8P(0._R8P, 0._R8P, 1._R8P)
   err(1) = max(err(1), abs(dragon%winding_number(point=q(1)) - dragon%winding_number(point=q(1), beta=-1._R8P)))
   err(2) = max(err(2), abs(dragon%winding_number(point=q(1), beta=4._R8P) - &
                            dragon%winding_number(point=q(1), beta=-1._R8P)))
enddo ; enddo ; enddo
print '(A,ES8.2,A,ES8.2)', 'dragon, max |fast - exact|: beta 2 ', err(1), ', beta 4 ', err(2)
endprogram ref_winding_number
