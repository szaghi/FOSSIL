p0 = pose(k) ; p1 = pose(k + 1)
call swept_cells_2d(plate, p0, p1, x, y, sw, map=map, status=status)
if (status /= STATUS_OK) error stop 'swept_cells_2d: map refused'
write(caption, '(A,I0,A,I0,A,I0,A,I0,A,ES8.1)') 'step ', k + 1, ': ', map%nt, ' tracks, ', &
      map%multi_track_cells, ' cells with several, ', map%untrusted_cells, ' untrusted, max residual ', &
      maxval(abs(map%residual))
