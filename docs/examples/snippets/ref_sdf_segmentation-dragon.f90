call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call dragon%segment_sdf(facet_labels=labels, sdf=sdf, status=status)    ! defaults: 4 clusters, 30 rays, 120 degrees
if (status /= SDF_STATUS_OK) error stop 'segment_sdf failed'
print '(A,I0,A,I0)', 'dragon: ', size(labels), ' facets, unassigned ', count(labels == SDF_LABEL_UNASSIGNED)
do k = 1, 4
   print '(A,I0,A,I5,A,F6.3,A,F6.3)', '  label ', k, ': ', count(labels == k), ' facets, sdf ', &
         minval(sdf, labels == k), ' to ', maxval(sdf, labels == k)
enddo
