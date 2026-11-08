call grid(1)%quality(sk, orth, asp, warp, mj)
print '(A,F6.3,A,F6.3,A,F6.2,A,F6.2,A,F6.3)', 'skewness max ', maxval(sk), ', orthogonality min ', minval(orth), &
      ', aspect max ', maxval(asp), ', warp max ', maxval(warp), ' deg, min J / mean ', minval(mj)
