call c%wall_distance_many(points=p, distances=d, bodies=bod, cutoff=CUTOFF)
print '(A,I0,A,I0,A,I0,A,I0)', 'wall distance, cutoff 0.4: ', NX * NY, ' points, body 1 ', count(bod == 1), &
                               ', body 2 ', count(bod == 2), ', beyond the cutoff ', count(bod == 0)
call field%initialize(x=xs, y=ys)
call field%compute(curve=c, status=status)
if (status /= STATUS_OK) error stop 'field: compute failed'
call field%gradient_diagnostic(q=q, flagged=flagged)
print '(A,I0,A,I0,A,F6.3)', 'gradient diagnostic: ', flagged, ' of ', NX * NY, ' nodes above 0.05, max q ', maxval(q)
print '(A,I0,A,I0)', '  flagged in the band ', count(q(:, :, 1) > 0.05_R8P .and. field%in_band(:, :, 1)), &
                     ', in the swept far field ', count(q(:, :, 1) > 0.05_R8P .and. .not. field%in_band(:, :, 1))
