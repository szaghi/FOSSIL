logical function naive_inside(q)
!< The naive sign: the normal of the closest facet, not the pseudo-normal of the closest feature.
type(vector_R8P), intent(in) :: q
type(facet_object), pointer  :: fp
type(vector_R8P)             :: c
real(R8P)                    :: d
integer(I4P)                 :: k

call wedge%compute_distance(point=q, distance=d, facet_index=k, closest_point=c)
fp => wedge%facet_at(k)
naive_inside = fp%normal%x * (q%x - c%x) + fp%normal%y * (q%y - c%y) + fp%normal%z * (q%z - c%z) < 0._R8P
endfunction naive_inside
