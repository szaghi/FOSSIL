call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
print '(A,I0,A,I0)', 'fresh (solid to fluid) ', count(fresh), ', dead (fluid to solid) ', count(dead)
print '(A,I0,A,I0,A)', 'second query needed for ', checked, ' of ', size(centres), ' centres'
