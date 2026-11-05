program ck_feature_edges
use fossil, only : surface_stl_object, FEATURE_EDGE_SHARP, FEATURE_EDGE_BOUNDARY, FEATURE_EDGE_NONMANIFOLD, STATUS_OK
use penf,   only : I4P, R8P
implicit none
type(surface_stl_object)  :: dragon
integer(I4P), allocatable :: kinds(:), corners(:), ps(:), pv(:)
integer(I4P)              :: status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0)', 'sharp ', count(kinds == FEATURE_EDGE_SHARP), ', boundary ', &
                          count(kinds == FEATURE_EDGE_BOUNDARY), ', non-manifold ', count(kinds == FEATURE_EDGE_NONMANIFOLD)
print '(A,I0,A,I0)', 'corners ', size(corners), ', polylines ', size(ps) - 1
call dragon%feature_edges(kinds=kinds, angle=60._R8P, include_boundary=.false., include_nonmanifold=.false.)
print '(A,I0)', 'sharp above 60 degrees only: ', size(kinds)
endprogram ck_feature_edges
