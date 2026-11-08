call vertex_area(sphere, area)                           ! M_ii of cotangent_laplacian: the area each vertex owns
print '(A,ES8.1)', 'sphere: sum K_i A_i / (2 pi) - 2 ', sum(K * area) / (2 * PI) - 2._R8P ! Euler characteristic 2
call implicit_body(torus, torus=.true.)                  ! a torus, R = 1, r = 0.45: genus 1
call torus%gaussian_curvature(K=K)
call vertex_area(torus, area)
print '(A,ES8.1)', 'torus:  sum K_i A_i / (2 pi)     ', sum(K * area) / (2 * PI)          ! Euler characteristic 0
