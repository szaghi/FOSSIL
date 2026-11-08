!< Figure of the reference page "Curvilinear structured blocks": the records detect_topology_2d proposes on two blocks.
!< Left, a C-grid around an ellipse with its wake to x = 3: two TOPO_ABUT records on face 3 (j = 1), the lower wake
!< i = 1..13 onto the upper one i = 61..49, reversed (teal and ember, drawn just apart). Right, a disc grid: the seam
!< (faces 1 and 2, gold) and the centre, where face 3 collapses (TOPO_SINGULAR, ember).
program fig_ref_curvilinear_topology
use fossil,     only : curvilinear_block_2d_t, detect_topology_2d, TOPO_ABUT, TOPO_SINGULAR
use fossil_fig, only : fig_open, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P
implicit none
real(R8P), parameter         :: PI = acos(-1._R8P)
character(9), parameter      :: KIND_NAME(4) = ['boundary ', 'abut     ', 'periodic ', 'singular ']
type(curvilinear_block_2d_t) :: c(1), d(1)
real(R8P)                    :: u(61, 7), v(61, 7), x(33, 9), y(33, 9), w0(2, 61), nrm(2, 61), phi, r
integer(I4P)                 :: i, j, k, w
character(60)                :: text

do i = 1, 61 ! the C-grid: 13 vertices per wake side, 36 around the ellipse (chord 1, thickness 0.2), 7 layers
   if (i <= 13) then
      w0(:, i) = [3._R8P - 2._R8P * (i - 1) / 12, 0._R8P] ; nrm(:, i) = [0._R8P, -1._R8P]
   elseif (i >= 49) then
      w0(:, i) = [3._R8P - 2._R8P * (61 - i) / 12, 0._R8P] ; nrm(:, i) = [0._R8P, 1._R8P]
   else
      phi = 2 * PI * (i - 13) / 36
      w0(:, i) = [0.5_R8P + 0.5_R8P * cos(phi), -0.1_R8P * sin(phi)]
      nrm(:, i) = [0.1_R8P * cos(phi), -0.5_R8P * sin(phi)] ; nrm(:, i) = nrm(:, i) / norm2(nrm(:, i))
   endif
enddo
do k = 1, 40 ! smooth the marching directions along i, so that the layers do not cross at the trailing edge
   nrm(:, 2:60) = (nrm(:, 1:59) + 2 * nrm(:, 2:60) + nrm(:, 3:61)) / 4
   do i = 2, 60
      nrm(:, i) = nrm(:, i) / norm2(nrm(:, i))
   enddo
enddo
do i = 1, 61
   do j = 1, 7
      u(i, j) = w0(1, i) + 0.15_R8P * (j - 1) * nrm(1, i) ; v(i, j) = w0(2, i) + 0.15_R8P * (j - 1) * nrm(2, i)
   enddo
enddo
call c(1)%initialize(u, v) ; call detect_topology_2d(c)
do j = 1, 9 ! the disc: radius 0.9 about (5.2, 0), 32 cells around, the centre collapsed
   r = 0.9_R8P * (j - 1) / 8
   do i = 1, 33
      x(i, j) = 5.2_R8P + r * cos(2 * PI * (i - 1) / 32) ; y(i, j) = r * sin(2 * PI * (i - 1) / 32)
   enddo
enddo
x(33,:) = x(1,:) ; y(33,:) = y(1,:)
call d(1)%initialize(x, y) ; call detect_topology_2d(d)
do k = 1, c(1)%nrec
   print '(A,I0,3A,I0,A,I0,A,I0,A,2I3,A,I0,A,2I3,A,2I3)', 'C-grid record ', k, ': ', &
         trim(KIND_NAME(c(1)%rec(k)%kind)), ', face ', c(1)%rec(k)%face, ', i ', c(1)%rec(k)%lo(1), '..', &
         c(1)%rec(k)%hi(1), ', begin', c(1)%rec(k)%begin(1:2), ' -> face ', c(1)%rec(k)%donor_face, ', begin', &
         c(1)%rec(k)%donor_begin(1:2), ', transform', c(1)%rec(k)%transform(1:2)
enddo
do k = 1, d(1)%nrec
   print '(A,I0,3A,I0,A,I0)', 'disc record ', k, ': ', trim(KIND_NAME(d(1)%rec(k)%kind)), ', face ', &
         d(1)%rec(k)%face, ', donor face ', d(1)%rec(k)%donor_face
enddo

call fig_open(w, 'ref-curvilinear-topology', 'records: a C-grid wake cut (teal, ember), a seam (gold) and a pole (ember)', &
              -90._R8P, 90._R8P, height=330_I4P)
call lines(u, v) ; call lines(x, y)
call fig_segments(w, 'probe', [(q(u(i, 1), v(i, 1) - 0.02_R8P), i = 1, 12)], [(q(u(i + 1, 1), v(i + 1, 1) - 0.02_R8P), &
                  i = 1, 12)])
call fig_points(w, 'probe', [(q(u(i, 1), v(i, 1) - 0.02_R8P), i = 1, 13)])
call fig_segments(w, 'hole', [(q(u(i, 1), v(i, 1) + 0.02_R8P), i = 49, 60)], [(q(u(i + 1, 1), v(i + 1, 1) + 0.02_R8P), &
                  i = 49, 60)])
call fig_points(w, 'hole', [(q(u(i, 1), v(i, 1) + 0.02_R8P), i = 49, 61)])
call fig_segments(w, 'accent', [(q(x(1, j), y(1, j)), j = 1, 8)], [(q(x(1, j + 1), y(1, j + 1)), j = 1, 8)])
call fig_points(w, 'hole', [q(x(1, 1), y(1, 1))])
write(text, '(A,I0,A,I0,A)') 'i = 1..13 <-> i = ', c(1)%rec(1)%donor_begin(1), '..', &
                             c(1)%rec(1)%donor_begin(1) + c(1)%rec(1)%transform(1) * 12, ', reversed'
call fig_label(w, 'probe', q(1.1_R8P, -1.2_R8P), trim(text))
call fig_label(w, 'accent', q(4.5_R8P, -1.2_R8P), 'seam: faces 1 and 2')
call fig_label(w, 'hole', q(4.6_R8P, 1.1_R8P), 'face 3 singular')
call fig_close(w)
contains
   subroutine lines(a, b)
   !< The grid lines of a 2-D block.
   real(R8P), intent(in)         :: a(:,:), b(:,:) !< Vertices.
   type(vector_R8P), allocatable :: s(:), e(:)
   integer(I4P)                  :: i, j, m

   allocate(s(2 * size(a)), e(2 * size(a))) ; m = 0
   do j = 1, size(a, 2)
      do i = 1, size(a, 1) - 1
         m = m + 1 ; s(m) = q(a(i, j), b(i, j)) ; e(m) = q(a(i + 1, j), b(i + 1, j))
      enddo
   enddo
   do j = 1, size(a, 2) - 1
      do i = 1, size(a, 1)
         m = m + 1 ; s(m) = q(a(i, j), b(i, j)) ; e(m) = q(a(i, j + 1), b(i, j + 1))
      enddo
   enddo
   call fig_segments(w, 'box', s(1:m), e(1:m))
   endsubroutine lines

   function q(a, b) result(pt)
   !< A point of the plane.
   real(R8P), intent(in) :: a, b !< Coordinates.
   type(vector_R8P)      :: pt   !< Point.

   pt = vector_R8P(a, b, 0._R8P)
   endfunction q
endprogram fig_ref_curvilinear_topology
