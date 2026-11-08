i = 24                                                ! the -x edges of column 24 lie on the line of step 1
print '(A,F9.6)', 'solid length of the line x = -0.025 from the apertures: ', sum(1._R8P - cut%aperture(1, i, :)) * h
c = maxloc(cut%wet_length)                            ! the cell with the longest wall
print '(A,I0,A,I0,A,F7.4,A,4F7.4)', 'cell (', c(1), ',', c(2), '): alpha ', cut%alpha(c(1), c(2)), &
                                    ', apertures -x +x -y +y', cut%aperture(:, c(1), c(2))
print '(A,F7.4,A,2F8.4)', '  wall length / h ', cut%wet_length(c(1), c(2)) / h, ', normal integral / h ', &
                          cut%wet_normal(:, c(1), c(2)) / h
print '(A,2F8.4)', '  aperture differences (+x - -x, +y - -y)    ', &
                   cut%aperture(2, c(1), c(2)) - cut%aperture(1, c(1), c(2)), &
                   cut%aperture(4, c(1), c(2)) - cut%aperture(3, c(1), c(2))
