!< FOSSIL, test the Liepa minimum-weight triangulation of boundary loops (issue #31, step 2).
!<
!< For each hole (cut as in fossil_test_boundary_loops): the patch has n - 2 triangles; every loop
!< edge vid(k) -> vid(k+1) appears in exactly one patch triangle in that direction (orientation
!< consistent with the mesh); every patch normal points out of the unit cube; on planar holes the
!< patch area equals the hole area and all normals equal the face normal; and the mesh plus the
!< patch is watertight with volume 1.
!<   a. cube minus 1 facet (3-loop); b. cube-02 minus 1 cell (planar 4-loop);
!<   c. cube-02 minus 2x2 cells (planar 8-loop); d. a hole across the edge between the z = 0 and
!<      y = 0 faces (non-planar 6-loop).

program fossil_test_triangulate_loop

use fossil,              only : surface_stl_object, facet_object
use fossil_hole_filling, only : boundary_loop_t, extract_boundary_loops, triangulate_loop
use penf,                only : I4P, R8P

implicit none

logical :: ok(4)
real(R8P), parameter :: NONE(6,0) = reshape([real(R8P) ::], [6, 0])

ok(1) = check('src/tests/cube.stl',    'a cube minus 1 facet  ', .true.,  NONE, -1._R8P)
ok(2) = check('src/tests/cube-02.stl', 'b planar 1 cell       ', .false., &
              reshape([0._R8P, 0.25_R8P, 0._R8P, 0.25_R8P, 0._R8P, 0._R8P], [6, 1]), 0.0625_R8P)
ok(3) = check('src/tests/cube-02.stl', 'c planar 2x2 cells    ', .false., &
              reshape([0._R8P, 0.5_R8P, 0._R8P, 0.5_R8P, 0._R8P, 0._R8P], [6, 1]), 0.25_R8P)
ok(4) = check('src/tests/cube-02.stl', 'd non-planar, 2 faces ', .false., &
              reshape([0._R8P, 0.25_R8P, 0._R8P, 0.25_R8P, 0._R8P, 0._R8P, &
                       0._R8P, 0.25_R8P, 0._R8P, 0._R8P, 0._R8P, 0.25_R8P], [6, 2]), -1._R8P)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(file_name, label, remove_first, boxes, hole_area) result(passed)
   character(*), intent(in)           :: file_name, label
   logical,      intent(in)           :: remove_first
   real(R8P),    intent(in)           :: boxes(:,:)  !< Centroid boxes (xmin,xmax,ymin,ymax,zmin,zmax) to remove.
   real(R8P),    intent(in)           :: hole_area   !< Expected planar hole area (< 0: not planar).
   logical                            :: passed
   type(surface_stl_object)           :: src, holed, filled
   type(facet_object), allocatable    :: keep(:), patch(:), all_f(:)
   type(facet_object), pointer        :: fp, facets(:)
   type(boundary_loop_t), allocatable :: loops(:)
   integer(I4P)                       :: i, k, c, n_keep, n_open, status, uses
   logical                            :: drop, orient_ok, outward_ok, planar_ok
   real(R8P)                          :: area
   real(R8P), parameter               :: E = 1.e-9_R8P

   call src%load_from_file(file_name=file_name, guess_format=.true.)
   call src%sanitize
   allocate(keep(src%get_facets_number())) ; n_keep = 0
   do i = 1, src%get_facets_number()
      fp => src%facet_at(i)
      drop = (remove_first .and. i == 1)
      do c = 1, size(boxes, dim=2)
         if (fp%centroid%x > boxes(1,c) - E .and. fp%centroid%x < boxes(2,c) + E .and. &
             fp%centroid%y > boxes(3,c) - E .and. fp%centroid%y < boxes(4,c) + E .and. &
             fp%centroid%z > boxes(5,c) - E .and. fp%centroid%z < boxes(6,c) + E) drop = .true.
      enddo
      if (drop) cycle
      n_keep = n_keep + 1 ; keep(n_keep) = fp
   enddo
   keep = keep(1:n_keep)
   call holed%adopt_facets(facets=keep)
   facets => holed%facets_ref()
   call extract_boundary_loops(facet=facets, loops=loops, n_open=n_open)
   passed = (size(loops) == 1 .and. n_open == 0)
   if (.not. passed) then
      print '(A,A,I0,A,I0)', label, ': unexpected loops=', size(loops), ' n_open=', n_open
      return
   endif
   call triangulate_loop(loop=loops(1), facet=facets, patch=patch, status=status)

   orient_ok = .true.
   do k = 1, size(loops(1)%vid)
      uses = 0
      do i = 1, size(patch)
         if (has_edge(patch(i), loops(1)%vid(k), loops(1)%vid(mod(k, size(loops(1)%vid)) + 1))) uses = uses + 1
      enddo
      if (uses /= 1) orient_ok = .false.
   enddo
   outward_ok = .true. ; planar_ok = .true. ; area = 0._R8P
   do i = 1, size(patch)
      if (patch(i)%normal%x * (patch(i)%centroid%x - 0.5_R8P) + patch(i)%normal%y * (patch(i)%centroid%y - 0.5_R8P) + &
          patch(i)%normal%z * (patch(i)%centroid%z - 0.5_R8P) <= 0._R8P) outward_ok = .false.
      if (hole_area > 0._R8P .and. abs(patch(i)%normal%z + 1._R8P) > 1.e-12_R8P) planar_ok = .false.
      area = area + patch(i)%area()
   enddo
   if (hole_area > 0._R8P) planar_ok = planar_ok .and. abs(area - hole_area) < 1.e-12_R8P

   allocate(all_f(size(facets) + size(patch)))
   all_f(1:size(facets)) = facets ; all_f(size(facets) + 1:) = patch
   call filled%adopt_facets(facets=all_f)

   passed = (size(patch) == size(loops(1)%vid) - 2 .and. status == 0) .and. orient_ok .and. outward_ok .and. &
            planar_ok .and. filled%is_watertight() .and. abs(filled%get_volume() - 1._R8P) < 1.e-12_R8P
   print '(A,A,I0,A,I0,A,L1,A,L1,A,L1,A,L1,A,F8.5,A,L1)', label, ': loop=', size(loops(1)%vid), ' patch=', size(patch), &
         ' orient=', orient_ok, ' outward=', outward_ok, ' planar=', planar_ok, ' watertight=', filled%is_watertight(), &
         ' volume=', filled%get_volume(), ' ok=', passed
   endfunction check

   pure function has_edge(f, a, b) result(yes)
   type(facet_object), intent(in) :: f
   integer(I4P),       intent(in) :: a, b
   logical                        :: yes
   integer(I4P)                   :: e
   yes = .false.
   do e = 1, 3
      if (f%vertex_id(e) == a .and. f%vertex_id(mod(e, 3) + 1) == b) yes = .true.
   enddo
   endfunction has_edge
endprogram fossil_test_triangulate_loop
