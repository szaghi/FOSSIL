call part%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv)
print '(A,I0,A,I0,A,I0)', 'part:     feature edges ', size(kinds), ', corners ', size(corners), &
                          ', polylines ', size(ps) - 1
call mesh%feature_edges(kinds=kinds, corners=corners, polyline_start=ps, polyline_vertices=pv, status=status)
if (status /= STATUS_OK) error stop 'feature_edges failed'
print '(A,I0,A,I0,A,I0,A,I0)', 'remeshed: feature edges ', size(kinds), ', corners ', size(corners), &
                               ', polylines ', size(ps) - 1, ', all sharp ', count(kinds == FEATURE_EDGE_SHARP)
n = count([(pv(ps(i)) == pv(ps(i + 1) - 1), i = 1, size(ps) - 1)])         ! first vertex = last: a loop
print '(A,I0,A,I0)', '  closed loops without a corner ', n, ', chains between corners ', size(ps) - 1 - n
