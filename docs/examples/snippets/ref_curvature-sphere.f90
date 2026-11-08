call implicit_body(sphere, torus=.false.)                ! unit sphere by marching cubes, then sanitize
call sphere%gaussian_curvature(K=K, status=status)       ! angle defect / barycentric area
if (status /= CURV_STATUS_OK) error stop 'gaussian_curvature failed'
call sphere%mean_curvature(H=H)                          ! |M^-1 L x| / 2, signed by the pseudo-normal
call sphere%jet_curvature(H=Hjet, K=Kjet)                ! cubic fit over the 2-ring
print '(A,I0,A,L1)', 'sphere: vertices ', size(K), ', watertight ', sphere%is_watertight()
print '(A)', '                   mean     2%    98%   (percentiles over the vertices)'
print '(A,3F7.3)', '  K angle defect', sum(K) / size(K), pct(K, 2), pct(K, 98)
print '(A,3F7.3)', '  H cotangent   ', sum(H) / size(H), pct(H, 2), pct(H, 98)
print '(A,3F7.3)', '  H jet         ', sum(Hjet) / size(Hjet), pct(Hjet, 2), pct(Hjet, 98)
