!run ck_ray_hit ck_ray_hit
program ck_ray_hit
use fossil, only : surface_stl_object, ray_hit_t
use penf,   only : R8P
use vecfor, only : vector_R8P
implicit none
type(surface_stl_object) :: dragon
type(ray_hit_t)          :: hit
logical                  :: has_hit

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%intersect_ray_first(ray_origin=vector_R8P(0._R8P, -20._R8P, 12._R8P), &
                                ray_direction=vector_R8P(0._R8P, 1._R8P, 0._R8P), hit=hit, has_hit=has_hit)
if (has_hit) then
   print '(A,I0,A,F8.4)', 'hit facet ', hit%facet_id, ' at t = ', hit%t
   print '(A,3F9.4)',     'hit point ', hit%point%x, hit%point%y, hit%point%z
else
   print '(A)', 'no hit'
endif
endprogram ck_ray_hit
