call surface%query_many(points=centres, distances=d_posed, closest_points=foot_posed, &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER, pose=p1)        ! the tree stays as built
moved = surface
call moved%move(p1)                                                        ! a moved copy, tree refitted
call moved%query_many(points=centres, distances=d_moved, closest_points=foot_moved, &
                      sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,ES9.2)', 'distances:      posed - moved, max ', maxval(abs(d_posed - d_moved))
print '(A,ES9.2)', 'closest points: posed - moved, max ', maxval(abs(foot_posed - foot_moved))
print '(A,I0,A,I0)', 'centres inside at pose1: ', count(d_posed < 0._R8P), ' of ', size(centres)
