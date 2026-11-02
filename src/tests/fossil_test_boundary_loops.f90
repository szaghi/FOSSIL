!< FOSSIL, test boundary-loop extraction for hole filling (issue #31, step 1).
!<
!< Holes are cut out of the unit cube (12 facets) and of cube-02 (each face a 4x4 grid of 0.25
!< cells, two triangles per cell) by removing facets whose centroid lies in given cells of the
!< z = 0 face. Expected loops (count, edge counts), n_open = 0 everywhere, and orientation: every
!< loop edge vid(k) -> vid(k+1) is the reverse of an edge of its owner facet.
!<   a. watertight cube              -> 0 loops
!<   b. cube minus 1 facet           -> 1 loop of 3
!<   c. cube-02 minus 1 cell         -> 1 loop of 4
!<   d. cube-02 minus 2x2 cells      -> 1 loop of 8
!<   e. cube-02 minus 2 far cells    -> 2 loops of 4
!<   f. pinch: 2 cells sharing a corner -> 2 loops of 4 (the shared vertex on both)
!<   g. open sheet (only the z = 0 face of cube-02) -> 1 loop of 16

program fossil_test_boundary_loops

use fossil,              only : surface_stl_object, facet_object
use fossil_hole_filling, only : boundary_loop_t, extract_boundary_loops
use penf,                only : I4P, R8P

implicit none

logical :: ok(7)
real(R8P), parameter :: NONE(4,0) = reshape([real(R8P) ::], [4, 0])

ok(1) = check('src/tests/cube.stl',    'a watertight cube      ', remove_first=.false., cells=NONE, keep_only_z0=.false., &
              expected=[integer(I4P) ::])
ok(2) = check('src/tests/cube.stl',    'b cube minus 1 facet   ', remove_first=.true.,  cells=NONE, keep_only_z0=.false., &
              expected=[3_I4P])
ok(3) = check('src/tests/cube-02.stl', 'c cube-02 minus 1 cell ', remove_first=.false., &
              cells=reshape([0._R8P, 0.25_R8P, 0._R8P, 0.25_R8P], [4, 1]), keep_only_z0=.false., expected=[4_I4P])
ok(4) = check('src/tests/cube-02.stl', 'd cube-02 minus 2x2    ', remove_first=.false., &
              cells=reshape([0._R8P, 0.5_R8P, 0._R8P, 0.5_R8P], [4, 1]), keep_only_z0=.false., expected=[8_I4P])
ok(5) = check('src/tests/cube-02.stl', 'e two separate holes   ', remove_first=.false., &
              cells=reshape([0._R8P, 0.25_R8P, 0._R8P, 0.25_R8P, 0.5_R8P, 0.75_R8P, 0.5_R8P, 0.75_R8P], [4, 2]), &
              keep_only_z0=.false., expected=[4_I4P, 4_I4P])
ok(6) = check('src/tests/cube-02.stl', 'f pinch (shared corner)', remove_first=.false., &
              cells=reshape([0._R8P, 0.25_R8P, 0._R8P, 0.25_R8P, 0.25_R8P, 0.5_R8P, 0.25_R8P, 0.5_R8P], [4, 2]), &
              keep_only_z0=.false., expected=[4_I4P, 4_I4P])
ok(7) = check('src/tests/cube-02.stl', 'g open sheet           ', remove_first=.false., cells=NONE, keep_only_z0=.true., &
              expected=[16_I4P])

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check(file_name, label, remove_first, cells, keep_only_z0, expected) result(passed)
   character(*), intent(in)           :: file_name, label
   logical,      intent(in)           :: remove_first, keep_only_z0
   real(R8P),    intent(in)           :: cells(:,:)       !< (xmin, xmax, ymin, ymax) per cell to remove on z = 0.
   integer(I4P), intent(in)           :: expected(:)      !< Expected loop sizes (any order).
   logical                            :: passed
   type(surface_stl_object)           :: src, holed
   type(facet_object), allocatable    :: keep(:)
   type(facet_object), pointer        :: fp, facets(:)
   type(boundary_loop_t), allocatable :: loops(:)
   integer(I4P)                       :: i, k, c, n_keep, n_open, sizes(64)
   logical                            :: drop, orient_ok, shared_ok

   call src%load_from_file(file_name=file_name, guess_format=.true.)
   call src%sanitize
   allocate(keep(src%get_facets_number())) ; n_keep = 0
   do i = 1, src%get_facets_number()
      fp => src%facet_at(i)
      drop = (remove_first .and. i == 1)
      if (keep_only_z0 .and. abs(fp%centroid%z) > 1.e-12_R8P) drop = .true.
      do c = 1, size(cells, dim=2)
         if (abs(fp%centroid%z) < 1.e-12_R8P .and. fp%centroid%x > cells(1,c) .and. fp%centroid%x < cells(2,c) .and. &
             fp%centroid%y > cells(3,c) .and. fp%centroid%y < cells(4,c)) drop = .true.
      enddo
      if (drop) cycle
      n_keep = n_keep + 1 ; keep(n_keep) = fp
   enddo
   keep = keep(1:n_keep)
   call holed%adopt_facets(facets=keep)
   facets => holed%facets_ref()
   call extract_boundary_loops(facet=facets, loops=loops, n_open=n_open)

   orient_ok = .true.
   do k = 1, size(loops)
      do i = 1, size(loops(k)%vid)
         if (.not. has_edge(facets(loops(k)%owner(i)), loops(k)%vid(mod(i, size(loops(k)%vid)) + 1), loops(k)%vid(i))) &
            orient_ok = .false.
      enddo
   enddo
   shared_ok = .true.
   if (label(1:1) == 'f' .and. size(loops) == 2) shared_ok = any([(any(loops(2)%vid == loops(1)%vid(i)), i=1,size(loops(1)%vid))])

   sizes = 0
   do k = 1, size(loops)
      sizes(k) = size(loops(k)%vid)
   enddo
   passed = (size(loops) == size(expected)) .and. (n_open == 0) .and. orient_ok .and. shared_ok
   if (passed .and. size(expected) > 0) passed = all(sort_small(sizes(1:size(loops))) == sort_small(expected))
   print '(A,A,I0,A,8I4)', label, ': loops=', size(loops), ' sizes=', sizes(1:min(8, size(loops)))
   print '(A,I0,A,L1,A,L1,A,L1)', '                           n_open=', n_open, ' orientation=', orient_ok, &
         ' shared vertex=', shared_ok, ' ok=', passed
   endfunction check

   pure function has_edge(f, a, b) result(yes)
   !< True if facet f has the directed edge a -> b (pool ids) in its winding.
   type(facet_object), intent(in) :: f
   integer(I4P),       intent(in) :: a, b
   logical                        :: yes
   integer(I4P)                   :: e
   yes = .false.
   do e = 1, 3
      if (f%vertex_id(e) == a .and. f%vertex_id(mod(e, 3) + 1) == b) yes = .true.
   enddo
   endfunction has_edge

   pure function sort_small(v) result(s)
   integer(I4P), intent(in) :: v(:)
   integer(I4P)             :: s(size(v)), i, j, t
   s = v
   do i = 2, size(s)
      t = s(i) ; j = i - 1
      do while (j >= 1)
         if (s(j) <= t) exit
         s(j + 1) = s(j) ; j = j - 1
      enddo
      s(j + 1) = t
   enddo
   endfunction sort_small
endprogram fossil_test_boundary_loops
