call cyl%feature_edges(edges=edges)
call cyl%build_features(status=status)
if (status /= STATUS_OK .or. .not. cyl%is_watertight()) error stop 'cylinder: build_features failed'
print '(A,I0,A,I0,A,F8.5)', 'cylinder: ', cyl%get_facets_number(), ' facets, feature edges ', size(edges, 2), &
                            ', volume ', cyl%get_volume()
do j = 1, N
   do i = 1, N
      p(i + N * (j - 1)) = vector_R8P(-1._R8P + H * (i - 0.5_R8P), 0._R8P, -0.5_R8P + H * (j - 0.5_R8P))
   enddo
enddo
call cyl%query_many(points=p, distances=d, feature_distances=fd, status=status)
if (status /= STATUS_OK) error stop 'query_many failed'
print '(A,F6.4,A,F6.4)', 'plane y = 0: min |phi| ', minval(abs(d)), ', min feature distance ', minval(fd)
print '(A,F6.4,A,F6.4)', 'at the side (x = 0.479, z = 0.479): |phi| ', abs(d(36 + N * 23)), &
                         ', feature distance ', fd(36 + N * 23)
