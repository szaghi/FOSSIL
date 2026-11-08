call torus%jet_curvature(H=Hjet, K=Kjet)
k1 = Hjet + sqrt(max(0._R8P, Hjet**2 - Kjet))            ! principal curvatures from H and K, k1 >= k2
k2 = Hjet - sqrt(max(0._R8P, Hjet**2 - Kjet))
call torus_exact(torus, exact)                           ! (rho - R) / (r rho) at every vertex, rho = |(x, y)|
k1 = k1 - 1._R8P / RT ; k2 = k2 - exact                 ! errors
print '(A,2F7.4)', 'torus:  k1 - 1/r,   rms and max ', sqrt(sum(k1**2) / size(k1)), maxval(abs(k1))
print '(A,2F7.4)', '        k2 - exact, rms and max ', sqrt(sum(k2**2) / size(k2)), maxval(abs(k2))
