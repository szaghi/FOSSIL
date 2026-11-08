call b%surface_markers(dx=0.1_R8P, centroids=mc, areas=ma, normals=mn, spacing_fraction=fraction, &
                       preserve_features=.true.)
print '(A,I0,A,F9.6,A,F5.3)', 'markers on one cube, dx 0.1: ', size(ma), ', area ', sum(ma), &
                              ', spacings in [0.8, 1.2] dx: ', fraction
