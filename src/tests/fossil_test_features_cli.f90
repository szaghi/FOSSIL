!< FOSSIL test: fossilizer --features and save_features_into_file (issue #61, roadmap #58 WP3, step 5).
!<
!< `fossilizer -i src/tests/cube-01.stl --features <file>`: the printed summary (24 sharp, 0 boundary, 0 non-manifold,
!< 8 corners, 12 polylines) and the legacy VTK file read back (20 points: the 8 corners and the 12 edge midpoints;
!< 8 VERTICES; 24 LINES; 32 cell scalars, kind 0 on the corners and FEATURE_EDGE_SHARP = 1 on the edges).
!< Runs exe/fossilizer from the repository root, like every FOSSIL test.
program fossil_test_features_cli

use penf, only : I4P

implicit none

character(*), parameter :: OUT = 'fossil_test_features_cli.txt'
character(*), parameter :: VTK = 'fossil_test_features_cli.vtk'
character(256)          :: line
integer(I4P)            :: u, ios, exitstat, sharp, bnd, nm, corners, polylines, npoints, nverts, nlines, ncell, k, ones
logical                 :: ok

call execute_command_line('exe/fossilizer -i src/tests/cube-01.stl --features '//VTK//' > '//OUT, exitstat=exitstat)
sharp = -1 ; corners = -1
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, '  sharp ') == 1) read(line, *) line(1:1), sharp, line(2:2), bnd, line(3:3), nm
   if (index(line, '  corners ') == 1) read(line, *) line(1:1), corners, line(2:2), polylines
enddo
close(u, status='delete')
npoints = -1 ; nverts = -1 ; nlines = -1 ; ncell = -1 ; ones = 0
open(newunit=u, file=VTK, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, 'POINTS ') == 1) read(line(8:), *) npoints
   if (index(line, 'VERTICES ') == 1) read(line(10:), *) nverts
   if (index(line, 'LINES ') == 1) read(line(7:), *) nlines
   if (index(line, 'CELL_DATA ') == 1) then
      read(line(11:), *) ncell
      read(u, '(A)') line ; read(u, '(A)') line
      do k = 1, ncell
         read(u, *) ios
         if (ios == 1) ones = ones + 1
      enddo
   endif
enddo
close(u, status='delete')
ok = exitstat == 0 .and. sharp == 24 .and. bnd == 0 .and. nm == 0 .and. corners == 8 .and. polylines == 12 .and. &
     npoints == 20 .and. nverts == 8 .and. nlines == 24 .and. ncell == 32 .and. ones == 24
print '(A,7I5,A,L1)', 'fossilizer --features: sharp, corners, polylines, points, verts, lines, cells', sharp, corners, &
      polylines, npoints, nverts, nlines, ncell, '  ok ', ok
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1
endprogram fossil_test_features_cli
