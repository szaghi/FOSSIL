do n = 1, 400 ! points spread over the block, by a fixed sequence
   t = 2 * PI * modulo(0.618034_R8P * n, 1._R8P) ; r = 1.05_R8P + 0.9_R8P * modulo(0.4142136_R8P * n, 1._R8P)
   px(n) = r * cos(t) ; py(n) = r * sin(t) ; pz(n) = 0.1_R8P + 0.8_R8P * modulo(0.7320508_R8P * n, 1._R8P)
enddo
call grid(1)%build_search
call grid(1)%locate_many(px, py, pz, cell, xi, st)
call grid(1)%locate_many(px, py, pz, cell_walk, xi_walk, st_walk, walk=.true., jumps=jumps, fallbacks=fallbacks)
print '(A,I0,A,L1,A,I0,A,I0)', 'inside ', count(st == LOC_INSIDE), ', walk equal bit for bit ', &
      all(cell_walk == cell) .and. all(xi_walk == xi) .and. all(st_walk == st), ', jumps ', jumps, ', fallbacks ', fallbacks
