!< Animation data: a NACA 0012 section pitching and translating through a Cartesian grid. Per frame: the fluid
!< fraction of every cut cell (curve%cut_cells) and the signed distance at every grid node (distance_many with the
!< pose of the frame). scripts/anim2svg.py draws it (docs/examples/images/airfoil-cut-cells.svg, airfoil-level-sets.svg).
program anim_airfoil
use fossil, only : curve_object, cut_cells_2d_t
use fossil_rigid_pose, only : rigid_pose_2d_t, pose_2d
use penf, only : I4P, R8P

implicit none
integer(I4P), parameter :: NX = 48, NY = 24, NF = 30, NS = 60 ! cells, frames, points per side of the section
real(R8P),    parameter :: X0 = -0.5_R8P, X1 = 2.5_R8P, Y0 = -0.75_R8P, Y1 = 0.75_R8P, PI = acos(-1._R8P)
type(curve_object)      :: body, moved
type(cut_cells_2d_t)    :: cut
type(rigid_pose_2d_t)   :: pose
real(R8P)               :: x(NX + 1), y(NY + 1), u(2 * NS), v(2 * NS), nodes(2, (NX + 1) * (NY + 1))
real(R8P)               :: d((NX + 1) * (NY + 1)), t, s, th, a, xc
integer(I4P)            :: i, j, k, f, unit

do i = 1, NX + 1
   x(i) = X0 + (X1 - X0) * (i - 1) / NX
enddo
do j = 1, NY + 1
   y(j) = Y0 + (Y1 - Y0) * (j - 1) / NY
enddo
k = 0
do j = 1, NY + 1
   do i = 1, NX + 1
      k = k + 1 ; nodes(:, k) = [x(i), y(j)]
   enddo
enddo
! NACA 0012, chord 1, quarter chord at the origin, counter-clockwise (the solid on the left of every segment):
! the upper side from the trailing edge to the leading edge, then the lower side back
do k = 1, NS
   s = real(k - 1, R8P) / NS ; xc = 0.5_R8P * (1._R8P + cos(PI * s))            ! 1 -> 0 (cosine spacing)
   th = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.1260_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1036_R8P * xc**4)
   u(k) = xc - 0.25_R8P ; v(k) = th
   xc = 0.5_R8P * (1._R8P - cos(PI * s))                                        ! 0 -> 1
   th = 0.6_R8P * (0.2969_R8P * sqrt(xc) - 0.1260_R8P * xc - 0.3516_R8P * xc**2 + 0.2843_R8P * xc**3 - 0.1036_R8P * xc**4)
   u(NS + k) = xc - 0.25_R8P ; v(NS + k) = -th
enddo
call body%adopt_loops(u=u, v=v, loop_sizes=[2 * NS])

open(newunit=unit, file='airfoil.dat', status='replace', action='write')
write(unit, '(4I6)') NX, NY, NF, 2 * NS
write(unit, '(4F10.4)') X0, X1, Y0, Y1
do f = 0, NF - 1
   t = real(f, R8P) / NF
   a = -12._R8P * PI / 180._R8P * sin(2._R8P * PI * t)                          ! pitch, one period per loop
   pose = pose_2d(angle=a, translation=[0.35_R8P + 1.2_R8P * (0.5_R8P - 0.5_R8P * cos(2._R8P * PI * t)), &
                                       0.12_R8P * sin(2._R8P * PI * t)])
   moved = body
   call moved%move(pose)
   call moved%cut_cells(x=x, y=y, cut=cut)
   call body%distance_many(points=nodes, distances=d, pose=pose)
   write(unit, '(A,I4)') 'frame', f
   write(unit, '(10(1X,F7.4))') ((cut%alpha(i, j), i = 1, NX), j = 1, NY)
   write(unit, '(10(1X,F8.4))') d
   write(unit, '(10(1X,F8.4))') (pose%r(1, 1) * u(k) + pose%r(1, 2) * v(k) + pose%t(1), &
                            pose%r(2, 1) * u(k) + pose%r(2, 2) * v(k) + pose%t(2), k = 1, 2 * NS)
enddo
close(unit)
print '(A,I0,A,I0,A,I0,A)', 'airfoil.dat: ', NF, ' frames of a ', NX, ' x ', NY, ' grid'
endprogram anim_airfoil
