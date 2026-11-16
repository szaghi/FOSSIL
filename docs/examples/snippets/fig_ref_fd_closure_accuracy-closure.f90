do kk = 1, M                                                ! u at the probes: the solver's interpolation
   ui(kk) = interp(probe(:, kk), cell(:, kk))
enddo
wall = field(pw)                                            ! Dirichlet datum of the wall at the foot
call lagrange([0._R8P, t%lengths], [wall, ui], eta, q, l0)  ! along the normal, at the ghost's distance
ed(k, regime) = max(ed(k, regime), abs(q - field(t%points(:, g))))
dwall = dot_product(grad(pw), pw / R)                       ! Neumann datum: d_n u along the closure's normal
q = neumann(t%lengths, ui, dwall, eta)
en(k, regime) = max(en(k, regime), abs(q - field(t%points(:, g))))
