call section%load_from_file(file_name='dragon-section.dat', status=status)        ! chapter 7: the section z = 0
if (status /= STATUS_OK) error stop 'cannot load dragon-section.dat'
call section%cut_cells(x=x, y=x, cut=cut2)                                         ! its 2-D cut cells, same grid
print '(A,I0,A,ES9.2)', 'faces z = 0: cut ', count(cut%face_z(:, :, 13) > 0._R8P .and. cut%face_z(:, :, 13) < 1._R8P), &
                        ', largest |face_z - 2-D alpha| ', maxval(abs(cut%face_z(:, :, 13) - cut2%alpha))
