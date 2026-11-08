p0 = pose(k) ; p1 = pose(k + 1)
call swept_cells_2d(square, p0, p1, x, y, sw, status=status)
if (status /= STATUS_OK .or. .not. sw%exact) error stop 'swept_cells_2d: no exact result'
worst = max(worst, maxval(abs(sw%residual)))
