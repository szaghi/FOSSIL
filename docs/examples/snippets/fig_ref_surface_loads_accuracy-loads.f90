call load_samples_2d(wall, samples=ls, order=2_I4P)                 ! FD: the load samples of order 2
p = [(pressure(ls%points(:, i)), i = 1, ls%n)]                      ! p of the true wall at each sample
call integrate_loads(ls, f, m, pressure=p)
e(k, 1) = abs(f(1)) / (RHO * U**2 * R) ; e(k, 2) = abs(f(2) - LIFT) / abs(LIFT)
h = 1.6_R8P / (n / 2)                                               ! FV: n / 2 cells per axis
x = [(ZC(1) - 0.8_R8P + 0.31_R8P * h + i * h, i = 0, n / 2)]
y = [(ZC(2) - 0.8_R8P + 0.17_R8P * h + i * h, i = 0, n / 2)]
call wall%cut_cells(x, y, cut, boundaries=.true.)
call cut_walls_2d(cut, walls)                                       ! the wall segments of the cut pieces
p = [(pressure(walls%centroid(:, i)), i = 1, walls%n)]              ! p of the true wall at each midpoint
call cut_loads_2d(walls, f, m, pressure=p)
e(k, 3) = abs(f(1)) / (RHO * U**2 * R) ; e(k, 4) = abs(f(2) - LIFT) / abs(LIFT)
