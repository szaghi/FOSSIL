call blocks(1)%build_search                                  ! validity, orientation, the cell-box tree
n = 200
allocate(pu(n), pv(n))
do i = 1, n                                                  ! points on a spiral around the body
   s = real(i - 1, R8P) / (n - 1)
   pu(i) = (0.9_R8P + 0.9_R8P * s) * cos(14 * PI * s) ; pv(i) = (0.3_R8P + 1.1_R8P * s) * sin(14 * PI * s)
enddo
call blocks(1)%locate_many(pu, pv, cell, xi, st)             ! robust: the cell-box tree, exact containment
print '(A,I0,A,I0,A,I0)', 'located ', count(st == LOC_INSIDE), ', outside ', count(st == LOC_OUTSIDE), ' of ', n
call blocks(1)%locate_many(pu, pv, cell_walk, xi_walk, st_walk, walk=.true., fallbacks=fallbacks)
print '(A,L1,A,I0)', 'walk equal to the robust path: ', all(cell_walk == cell) .and. all(xi_walk == xi) .and. &
      all(st_walk == st), ', fallbacks ', fallbacks
