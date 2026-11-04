call surface%distance_many(points=centres, distances=d_posed, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, pose=p1)    ! the tree stays as built
moved = surface
call moved%move(p1)                                                       ! a moved copy, tree refitted
call moved%distance_many(points=centres, distances=d_moved, is_signed=.true., is_square_root=.true., &
                         sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,L1)', 'posed = moved to 1e-14: ', maxval(abs(d_posed - d_moved)) < 1.e-14_R8P
call surface%state_changes(pose0=p0, pose1=p1, points=centres, fresh=fresh, dead=dead, &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER, checked=checked)
print '(A,I0,A,I0,A,I0,A,I0)', 'cells: fresh ', count(fresh), ', dead ', count(dead), ', second query ', checked, &
                               ' of ', size(centres)
