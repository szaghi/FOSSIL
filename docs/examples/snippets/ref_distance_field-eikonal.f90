allocate(q(n, n, n))
call field%gradient_diagnostic(q=q, flagged=flagged)
print '(A,I0,A,I0,A,I0)', 'eikonal q > 0.05 at ', flagged, ' of ', n**3, ' nodes; unknown neighbour ', count(q < 0._R8P)
print '(A,I0,A,I0)', '  of them in the band ', count(q > 0.05_R8P .and. field%in_band), ', inside the cube ', &
                     count(q > 0.05_R8P .and. field%phi < 0._R8P)
