!< FOSSIL test: fossilizer --thickness and save_thickness_into_file (issue #62, roadmap #58 WP4, step 5).
!<
!< `fossilizer -i src/tests/cube-01.stl --thickness <file> --thickness_method medial --h 1` against the library on the
!< same file (thickness(method=THICKNESS_MEDIAL), thin_regions(h=1)): the printed thinnest wall (16 digits), "no gap" (the cube is
!< convex), the count of facets under h, the regions and the largest area; the legacy VTK file read back (the vertex
!< pool as POINTS, one triangle per facet in POLYGONS, the thickness bit for bit (17 digits), gap -1, the regions).
!< With --h 0.2 (below every medial value of the cube) nothing is reported. Runs exe/fossilizer from the repository
!< root, like every FOSSIL test.
program fossil_test_thickness_cli

use fossil, only : surface_stl_object, THICKNESS_MEDIAL
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,   only : I4P, R8P

implicit none

character(*), parameter   :: OUT = 'fossil_test_thickness_cli.txt'
character(*), parameter   :: VTK = 'fossil_test_thickness_cli.vtk'
type(surface_stl_object)  :: s
type(vertex_pool_object), pointer :: pool
real(R8P), allocatable    :: th(:), gp(:), area(:)
integer(I4P), allocatable :: region(:)
character(256)            :: line, name
real(R8P)                 :: thinnest, larea, val
integer(I4P)              :: u, ios, exitstat, nfac, nreg, npoints, npoly, ncell, k, good, under
logical                   :: ok, nogap

call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call s%thickness(thickness=th, gap=gp, method=THICKNESS_MEDIAL)
call s%thin_regions(thickness=th, gap=gp, h=1._R8P, factor=1._R8P, region=region, region_area=area)
call execute_command_line('exe/fossilizer -i src/tests/cube-01.stl --thickness '//VTK//' --thickness_method medial '// &
                          '--h 1 > '//OUT, exitstat=exitstat)
thinnest = -1._R8P ; larea = -1._R8P ; nfac = -1 ; nreg = -1 ; nogap = .false.
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, '  thinnest wall ') == 1) read(line(16:), *) thinnest
   if (index(line, '  narrowest gap  none') == 1) nogap = .true.
   if (index(line, '  under ') == 1) read(line(index(line, 'facets') + 6:), *) nfac, name, nreg
   if (index(line, '  largest region area ') == 1) read(line(23:), *) larea
enddo
close(u, status='delete')
ok = exitstat == 0 .and. abs(thinnest - minval(th)) <= 1e-14_R8P .and. nogap .and. all(gp == huge(1._R8P)) .and. &
     nfac == count(region > 0) .and. nreg == size(area) .and. nreg > 0
if (ok) ok = abs(larea - area(1)) <= 1e-14_R8P * area(1)
print '(A,ES12.4,A,L1,A,I0,A,I0,A,ES12.4,A,L1)', '1. summary: thinnest ', thinnest, ', no gap ', nogap, ', facets ', &
                                               nfac, ', regions ', nreg, ', area ', larea, '  ok ', ok
npoints = -1 ; npoly = -1 ; ncell = -1 ; good = 0
open(newunit=u, file=VTK, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, 'POINTS ') == 1) read(line(8:), *) npoints
   if (index(line, 'POLYGONS ') == 1) read(line(10:), *) npoly
   if (index(line, 'CELL_DATA ') == 1) read(line(11:), *) ncell
   if (index(line, 'SCALARS ') == 1) then
      name = line(9:)
      read(u, '(A)') line
      do k = 1, ncell
         read(u, *) val
         if (index(name, 'thickness') == 1 .and. val == th(k)) good = good + 1
         if (index(name, 'gap') == 1 .and. val == -1._R8P) good = good + 1
         if (index(name, 'region') == 1 .and. nint(val) == region(k)) good = good + 1
      enddo
   endif
enddo
close(u, status='delete')
pool => s%get_vertex_pool()
ok = ok .and. npoints == pool%vertex_count() .and. npoly == s%get_facets_number() .and. &
     ncell == s%get_facets_number() .and. good == 3 * ncell
print '(A,I0,A,I0,A,I0,A,L1)', '2. VTK: points ', npoints, ', polygons ', npoly, ', scalar values equal ', good, '  ok ', ok
call execute_command_line('exe/fossilizer -i src/tests/cube-01.stl --thickness '//VTK//' --thickness_method medial '// &
                          '--h 0.2 > '//OUT, exitstat=exitstat)
under = -1
open(newunit=u, file=OUT, status='old')
do
   read(u, '(A)', iostat=ios) line
   if (ios /= 0) exit
   if (index(line, '  under ') == 1) read(line(index(line, 'facets') + 6:), *) under
enddo
close(u, status='delete')
open(newunit=u, file=VTK, status='old') ; close(u, status='delete')
ok = ok .and. under == 0 .and. minval(th) > 0.2_R8P
print '(A,I0,A,L1)', '3. --h 0.2: facets under ', under, '  ok ', ok
print '(A,L1)', 'Are all tests passed? ', ok
endprogram fossil_test_thickness_cli
