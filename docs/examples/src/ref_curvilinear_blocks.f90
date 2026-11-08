!run ref_curvilinear_blocks ref_curvilinear_blocks
program ref_curvilinear_blocks
use fossil, only : curvilinear_block_t, detect_topology, TOPO_ABUT, TOPO_PERIODIC, CELL_VALID, LOC_INSIDE
use penf,   only : I4P, R8P
implicit none
real(R8P), parameter      :: PI = acos(-1._R8P)
character(9), parameter   :: KIND_NAME(4) = ['boundary ', 'abut     ', 'periodic ', 'singular ']
type(curvilinear_block_t) :: grid(1)
real(R8P)                 :: x(33, 9, 9), y(33, 9, 9), z(33, 9, 9), r, t, px(400), py(400), pz(400), f, e(2)
real(R8P), allocatable    :: vol(:,:,:), fi(:,:,:,:), fj(:,:,:,:), fk(:,:,:,:), xi(:,:), xi_walk(:,:), w(:,:)
real(R8P), allocatable    :: sk(:,:,:), orth(:,:,:), asp(:,:,:), warp(:,:,:), mj(:,:,:)
integer(I4P), allocatable :: valid(:,:,:), cell(:,:), st(:), cell_walk(:,:), st_walk(:), nodes(:,:,:), sst(:)
integer(I4P)              :: i, j, k, l, n, q, jumps, fallbacks
real(R8P)                 :: gcl

!region block
do k = 1, 9 ! a wavy annulus 1 <= r <= 2, 32 x 8 cells around (clockwise) and out, 8 cells over 0 <= z <= 1
   do j = 1, 9
      do i = 1, 33
         t = -2 * PI * (i - 1) / 32 ; r = 1._R8P + (j - 1) / 8._R8P
         x(i, j, k) = r * cos(t) ; y(i, j, k) = r * sin(t)
         z(i, j, k) = (k - 1) / 8._R8P + 0.05_R8P * sin(2 * t) * sin(PI * (j - 1) / 8) * sin(PI * (k - 1) / 8)
      enddo
   enddo
enddo
x(33,:,:) = x(1,:,:) ; y(33,:,:) = y(1,:,:) ; z(33,:,:) = z(1,:,:)              ! the seam, exact
call grid(1)%initialize(x, y, z)
call detect_topology(grid, translations=reshape([0._R8P, 0._R8P, 1._R8P], [3, 1]))
do n = 1, grid(1)%nrec
   print '(A,I0,3A,I0,A,I0,A,3I3)', 'record ', n, ': ', KIND_NAME(grid(1)%rec(n)%kind), ' face ', &
         grid(1)%rec(n)%face, ' -> face ', grid(1)%rec(n)%donor_face, ', transform', grid(1)%rec(n)%transform
enddo
!endregion block

!region geometry
call grid(1)%geometry(volume=vol, face_i=fi, face_j=fj, face_k=fk)
gcl = 0._R8P ! the six outward face vectors of each cell sum to zero
do k = 1, 8
   do j = 1, 8
      do i = 1, 32
         gcl = max(gcl, maxval(abs(fi(:, i + 1, j, k) - fi(:, i, j, k) + fj(:, i, j + 1, k) - fj(:, i, j, k) + &
                                   fk(:, i, j, k + 1) - fk(:, i, j, k))))
      enddo
   enddo
enddo
print '(A,F12.9,A,ES9.2)', 'volume ', sum(vol), ', largest face-vector sum of a cell ', gcl
call grid(1)%validity(valid)
print '(A,I0,A,I0)', 'cells certified valid ', count(valid == CELL_VALID), ' of ', size(valid)
!endregion geometry

!region locate
do n = 1, 400 ! points spread over the block, by a fixed sequence
   t = 2 * PI * modulo(0.618034_R8P * n, 1._R8P) ; r = 1.05_R8P + 0.9_R8P * modulo(0.4142136_R8P * n, 1._R8P)
   px(n) = r * cos(t) ; py(n) = r * sin(t) ; pz(n) = 0.1_R8P + 0.8_R8P * modulo(0.7320508_R8P * n, 1._R8P)
enddo
call grid(1)%build_search
call grid(1)%locate_many(px, py, pz, cell, xi, st)
call grid(1)%locate_many(px, py, pz, cell_walk, xi_walk, st_walk, walk=.true., jumps=jumps, fallbacks=fallbacks)
print '(A,I0,A,L1,A,I0,A,I0)', 'inside ', count(st == LOC_INSIDE), ', walk equal bit for bit ', &
      all(cell_walk == cell) .and. all(xi_walk == xi) .and. all(st_walk == st), ', jumps ', jumps, ', fallbacks ', fallbacks
!endregion locate

!region stencils
do q = 2, 4, 2
   call grid(1)%stencils(cell, xi, q, nodes, w, sst)
   e(q / 2) = 0._R8P
   do n = 1, 400 ! interpolate f = sin(x) cos(y) cos(2 pi z), periodic in z as the grid, from the vertices
      f = 0._R8P
      do k = 1, q**3
         i = nodes(1, k, n) ; j = nodes(2, k, n) ; l = nodes(3, k, n)
         f = f + w(k, n) * sin(x(i, j, l)) * cos(y(i, j, l)) * cos(2 * PI * z(i, j, l))
      enddo
      e(q / 2) = max(e(q / 2), abs(f - sin(px(n)) * cos(py(n)) * cos(2 * PI * pz(n))))
   enddo
enddo
print '(A,2ES10.2)', 'largest interpolation error, widths 2 and 4 ', e
!endregion stencils

!region quality
call grid(1)%quality(sk, orth, asp, warp, mj)
print '(A,F6.3,A,F6.3,A,F6.2,A,F6.2,A,F6.3)', 'skewness max ', maxval(sk), ', orthogonality min ', minval(orth), &
      ', aspect max ', maxval(asp), ', warp max ', maxval(warp), ' deg, min J / mean ', minval(mj)
!endregion quality
endprogram ref_curvilinear_blocks
