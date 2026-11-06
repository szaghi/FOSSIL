!< FOSSIL test: named patches on 2-D curves (issue #66, step 7).
!<
!< 1. The loop format: `# patch <name>` sets the patch of the segments that follow (segment i = vertex i -> next);
!<    a unit square with three `wall` sides and one `inlet` side reads back with its patches, and a round trip keeps
!<    them. A file without directives is one unnamed patch and is written back without directives.
!< 2. repair: two open chains, wound clockwise, joined, closed and reversed: every segment keeps its patch (checked by
!<    where it lies: the side y = 0 is `inlet`).
!< 3. query_many(patches=), query_many_device(patches=), probe_many and image_points: the patch of the closest segment.
!< 4. slice_surface: each segment takes the patch of the facet it is cut from (the unit cube with its x faces `inlet`,
!<    cut at z = 0.5: the segments on x = 0 and x = 1 are `inlet`).
program fossil_test_patches_curve

use fossil, only : surface_stl_object, curve_object, facet_object, slice_surface
use penf,   only : I4P, I8P, R8P

implicit none

character(*), parameter     :: FN = 'fossil_test_patches_curve'
type(curve_object)          :: c, r
type(surface_stl_object)    :: cube
type(facet_object), pointer :: f
integer(I4P)                :: u, i, k, nq, pi
integer(I4P), allocatable   :: seg(:), pat(:), pat2(:), pat3(:), pat4(:), ids(:)
real(R8P),    allocatable   :: q(:,:), d(:), probes(:,:), images(:,:), lens(:)
logical,      allocatable   :: vis(:)
character(256)              :: line
logical                     :: ok1, ok2, ok3, ok4, any_directive
real(R8P)                   :: a(2), b(2)

! ---- 1. format
open(newunit=u, file=FN//'.dat', status='replace')
write(u, '(A)') '# patch wall'
write(u, '(A)') '0 0'
write(u, '(A)') '# patch wall'
write(u, '(A)') '1 0'
write(u, '(A)') '1 1'
write(u, '(A)') '# patch inlet'
write(u, '(A)') '0 1'
close(u)
call c%load_from_file(FN//'.dat')
ok1 = c%get_patches_number() == 2 .and. c%patch_name(1) == 'wall' .and. c%patch_name(2) == 'inlet' .and. &
      all([(c%segment_patch(i), i = 1, 4)] == [1, 1, 1, 2])
call c%save_into_file(FN//'-rt.dat')
call r%load_from_file(FN//'-rt.dat')
ok1 = ok1 .and. r%get_patches_number() == 2 .and. all([(r%segment_patch(i), i = 1, 4)] == [1, 1, 1, 2])
open(newunit=u, file=FN//'.dat', status='replace')
write(u, '(A)') '0 0' ; write(u, '(A)') '1 0' ; write(u, '(A)') '0 1'
close(u)
call r%load_from_file(FN//'.dat')
call r%save_into_file(FN//'-rt.dat')
any_directive = .false.
open(newunit=u, file=FN//'-rt.dat', action='read')
do
   read(u, '(A)', end=10) line
   if (index(line, '# patch') > 0) any_directive = .true.
enddo
10 close(u)
ok1 = ok1 .and. r%get_patches_number() == 1 .and. r%patch_name(1) == '' .and. .not. any_directive
print '(A,L1)', '1. format: ', ok1

! ---- 2. repair: chains (1,1)->(1,0)->(0,0) [inlet on y = 0] and (0,0)->(0,1)->(1,1) [wall], clockwise
open(newunit=u, file=FN//'.dat', status='replace')
write(u, '(A)') '# open'
write(u, '(A)') '# patch wall'
write(u, '(A)') '1 1'
write(u, '(A)') '# patch inlet'
write(u, '(A)') '1 0'
write(u, '(A)') '0 0'
write(u, '(A)') ''
write(u, '(A)') '# open'
write(u, '(A)') '# patch wall'
write(u, '(A)') '0 0'
write(u, '(A)') '0 1'
write(u, '(A)') '1 1'
close(u)
call c%load_from_file(FN//'.dat')
call c%repair(tolerance=1.e-12_R8P)
ok2 = c%get_loops_number() == 1 .and. c%is_loop_closed(1) .and. c%get_area() > 0._R8P
do i = 1, c%get_vertices_number()
   k = c%get_next(i)
   if (k == 0) cycle
   a = c%get_vertex(i) ; b = c%get_vertex(k)
   ! the segment on y = 0 is inlet, the others wall
   if (a(2) == 0._R8P .and. b(2) == 0._R8P) then
      ok2 = ok2 .and. c%patch_name(c%segment_patch(i)) == 'inlet'
   else
      ok2 = ok2 .and. c%patch_name(c%segment_patch(i)) == 'wall'
   endif
enddo
print '(A,L1)', '2. repair: ', ok2

! ---- 3. queries
nq = 200
allocate(q(2, nq), d(nq), seg(nq), pat(nq), pat2(nq), pat3(nq), pat4(nq), probes(2, nq), images(2, nq), lens(nq), vis(nq))
do i = 1, nq
   q(:, i) = [mod(37._R8P * i, 1.6_R8P) - 0.3_R8P, mod(53._R8P * i, 1.6_R8P) - 0.3_R8P]
enddo
call c%query_many(points=q, distances=d, segments=seg, patches=pat)
call c%enter_device
call c%query_many_device(points=q, distances=d, patches=pat2)
call c%exit_device
lens = 0.05_R8P
call c%probe_many(points=q, lengths=lens, probe_points=probes, visible=vis, patches=pat3)
call c%image_points(points=q, images=images, patches=pat4)
ok3 = all(pat == c%segment_patch(seg)) .and. all(pat2 == pat) .and. all(pat3 == pat) .and. all(pat4 == pat)
print '(A,L1)', '3. queries: ', ok3

! ---- 4. slice_surface
call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
allocate(ids(0))
do i = 1, cube%get_facets_number()
   f => cube%facet_at(i)
   if (abs(f%normal%x) > 0.5_R8P) ids = [ids, i]
enddo
call cube%set_patch(ids, 'inlet')
call slice_surface(surface=cube, normal_axis=3_I4P, level=0.5_R8P, curve=r)
ok4 = r%get_vertices_number() > 0
pi = 0
do i = 1, r%get_patches_number()
   if (r%patch_name(i) == 'inlet') pi = i
enddo
do i = 1, r%get_vertices_number()
   k = r%get_next(i)
   if (k == 0) cycle
   a = r%get_vertex(i) ; b = r%get_vertex(k)
   if (a(1) == b(1) .and. (a(1) == 0._R8P .or. a(1) == 1._R8P)) then
      ok4 = ok4 .and. r%segment_patch(i) == pi
   else
      ok4 = ok4 .and. r%segment_patch(i) /= pi
   endif
enddo
print '(A,L1)', '4. slice_surface: ', ok4

open(newunit=u, file=FN//'.dat') ; close(u, status='delete')
open(newunit=u, file=FN//'-rt.dat') ; close(u, status='delete')
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4
endprogram fossil_test_patches_curve
