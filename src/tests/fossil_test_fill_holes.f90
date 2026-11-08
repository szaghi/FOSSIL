!< FOSSIL, test `surface%fill_holes` (issue #31, step 3).
!<
!< Holes on the unit cube are large relative to the surface (a missing facet spans 82% of the bbox
!< diagonal), so cases a, b, f pass max_hole_size explicitly; the default 0.1 is meant for real
!< meshes, where it keeps genuine open boundaries open (case d).
!<   a. cube minus 1 facet: 1 hole filled; watertight, volume 1.
!<   b. cube-02 minus 2x2 cells and one far cell: 2 holes filled; watertight, volume 1.
!<   c. dragon (it replaced dragon-fine, issue #45): it has boundary loops (else the case passes
!<      trivially), filled + skipped = boundary loops before, and the loops left afterwards are
!<      exactly the skipped ones.
!<   d. open sheet (the z = 0 face of cube-02): its single loop is as large as the surface, so the
!<      default max_hole_size keeps it open; the facet count is unchanged.
!<   e. max_hole_edges below the loop size: the hole is skipped.
!<   g. issue #88: open_edges counts true boundary edges (one incident facet) only, nonmanifold_edges the edges
!<      with 3+ facets met along the boundary. The dragon after sanitize (duplicates gone): 0 boundary edges left
!<      after filling, 54 non-manifold edges (tutorial 2's measurement); before #88 open_edges counted the facet edges
!<      on non-manifold edges as boundary edges. The cube minus 1 facet keeps 0 / 0.
!<   f. signed distance on the filled cube-02 equals the intact cube-02 at grid points away
!<      from the hole (1e-12).

program fossil_test_fill_holes

use fossil,              only : surface_stl_object, facet_object
use fossil_hole_filling, only : boundary_loop_t, extract_boundary_loops
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none

logical :: ok(7)
real(R8P), parameter :: CELL_2X2(6,2) = reshape([0._R8P, 0.5_R8P, 0._R8P, 0.5_R8P, 0._R8P, 0._R8P, &
                                                 0.75_R8P, 1._R8P, 0.75_R8P, 1._R8P, 0._R8P, 0._R8P], [6, 2])
type(surface_stl_object) :: s, ref, dragon
integer(I4P)             :: filled, skipped, open_edges, loops_before, loops_after, n0, nonmanifold

! a.
call holed(s, 'src/tests/cube.stl', remove_first=.true., boxes=reshape([real(R8P) ::], [6, 0]), keep_only_z0=.false.)
call s%fill_holes(max_hole_size=1._R8P, holes_filled=filled, holes_skipped=skipped, open_edges=open_edges)
ok(1) = filled == 1 .and. skipped == 0 .and. open_edges == 0 .and. s%is_watertight() .and. &
        abs(s%get_volume() - 1._R8P) < 1.e-12_R8P
print '(A,3I4,A,L1,A,F9.6,A,L1)', 'a cube minus 1 facet: filled/skipped/open', filled, skipped, open_edges, &
      ' watertight=', s%is_watertight(), ' volume=', s%get_volume(), ' ok=', ok(1)
! b.
call holed(s, 'src/tests/cube-02.stl', remove_first=.false., boxes=CELL_2X2, keep_only_z0=.false.)
call s%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled, holes_skipped=skipped)
ok(2) = filled == 2 .and. skipped == 0 .and. s%is_watertight() .and. abs(s%get_volume() - 1._R8P) < 1.e-12_R8P
print '(A,2I4,A,L1,A,F9.6,A,L1)', 'b cube-02 two holes: filled/skipped', filled, skipped, ' watertight=', s%is_watertight(), &
      ' volume=', s%get_volume(), ' ok=', ok(2)
! c.
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call s%sanitize
dragon = s                                                                    ! for case g
loops_before = count_loops(s)
call s%fill_holes(holes_filled=filled, holes_skipped=skipped, open_edges=open_edges)
loops_after = count_loops(s)
ok(3) = (loops_before > 0) .and. (filled + skipped == loops_before) .and. (loops_after == skipped)
print '(A,I0,A,I0,A,I0,A,I0,A,I0,A,L1,A,L1)', 'c dragon: loops before=', loops_before, ' filled=', filled, &
      ' skipped=', skipped, ' open edges=', open_edges, ' loops after=', loops_after, ' watertight=', s%is_watertight(), &
      ' ok=', ok(3)
! d.
call holed(s, 'src/tests/cube-02.stl', remove_first=.false., boxes=reshape([real(R8P) ::], [6, 0]), keep_only_z0=.true.)
n0 = s%get_facets_number()
call s%fill_holes(holes_filled=filled, holes_skipped=skipped)
ok(4) = filled == 0 .and. skipped == 1 .and. s%get_facets_number() == n0
print '(A,2I4,A,L1)', 'd open sheet: filled/skipped', filled, skipped, ' ok=', ok(4)
! e.
call holed(s, 'src/tests/cube-02.stl', remove_first=.false., boxes=CELL_2X2(:, 1:1), keep_only_z0=.false.)
call s%fill_holes(max_hole_edges=7_I4P, max_hole_size=0.5_R8P, holes_filled=filled, holes_skipped=skipped)
ok(5) = filled == 0 .and. skipped == 1 .and. .not. s%is_watertight()
print '(A,2I4,A,L1)', 'e max_hole_edges=7 on an 8-loop: filled/skipped', filled, skipped, ' ok=', ok(5)
! f.
call holed(s, 'src/tests/cube-02.stl', remove_first=.false., boxes=CELL_2X2(:, 1:1), keep_only_z0=.false.)
call s%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled)
call ref%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
call ref%sanitize
ok(6) = filled == 1 .and. s%is_watertight() .and. same_far_distances(s, ref)
print '(A,L1)', 'f signed distance away from the patch unchanged: ', ok(6)

! g.
s = dragon
call s%fill_holes(holes_filled=filled, holes_skipped=skipped, open_edges=open_edges, nonmanifold_edges=nonmanifold)
ok(7) = open_edges == 0 .and. nonmanifold == 54 .and. nonmanifold == s%get_non_manifold_edges_number()
print '(A,I0,A,I0,A,I0,A,I0,A,L1)', 'g dragon: filled=', filled, ' open edges=', open_edges, ' non-manifold=', &
      nonmanifold, ' (surface ', s%get_non_manifold_edges_number(), ') ok=', ok(7)
call holed(s, 'src/tests/cube.stl', remove_first=.true., boxes=reshape([real(R8P) ::], [6, 0]), keep_only_z0=.false.)
call s%fill_holes(max_hole_size=1._R8P, open_edges=open_edges, nonmanifold_edges=nonmanifold)
ok(7) = ok(7) .and. open_edges == 0 .and. nonmanifold == 0
print '(A,I0,A,I0,A,L1)', 'g cube minus 1 facet: open edges=', open_edges, ' non-manifold=', nonmanifold, ' ok=', ok(7)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   subroutine holed(s, file_name, remove_first, boxes, keep_only_z0)
   !< Load, sanitize, remove facets (first one, or centroids in boxes, or all but the z = 0 face).
   type(surface_stl_object), intent(inout) :: s
   character(*),             intent(in)    :: file_name
   logical,                  intent(in)    :: remove_first, keep_only_z0
   real(R8P),                intent(in)    :: boxes(:,:)
   type(surface_stl_object)                :: src
   type(facet_object), allocatable         :: keep(:)
   type(facet_object), pointer             :: fp
   integer(I4P)                            :: i, c, n
   logical                                 :: drop
   real(R8P), parameter                    :: E = 1.e-9_R8P

   call src%load_from_file(file_name=file_name, guess_format=.true.)
   call src%sanitize
   allocate(keep(src%get_facets_number())) ; n = 0
   do i = 1, src%get_facets_number()
      fp => src%facet_at(i)
      drop = (remove_first .and. i == 1)
      if (keep_only_z0 .and. abs(fp%centroid%z) > 1.e-12_R8P) drop = .true.
      do c = 1, size(boxes, dim=2)
         if (fp%centroid%x > boxes(1,c) - E .and. fp%centroid%x < boxes(2,c) + E .and. &
             fp%centroid%y > boxes(3,c) - E .and. fp%centroid%y < boxes(4,c) + E .and. &
             fp%centroid%z > boxes(5,c) - E .and. fp%centroid%z < boxes(6,c) + E) drop = .true.
      enddo
      if (drop) cycle
      n = n + 1 ; keep(n) = fp
   enddo
   keep = keep(1:n)
   call s%adopt_facets(facets=keep)
   endsubroutine holed

   function count_loops(s) result(n)
   type(surface_stl_object), intent(in) :: s
   integer(I4P)                         :: n, n_open
   type(boundary_loop_t), allocatable   :: loops(:)
   type(facet_object), pointer          :: facets(:)
   facets => s%facets_ref()
   call extract_boundary_loops(facet=facets, loops=loops, n_open=n_open)
   n = size(loops)
   endfunction count_loops

   function same_far_distances(a, b) result(same)
   !< Signed distances equal at grid points farther than 0.6 from the hole centre (0.25, 0.25, 0).
   type(surface_stl_object), intent(in) :: a, b
   logical                              :: same
   type(vector_R8P)                     :: p
   integer(I4P)                         :: i, j, k
   same = .true.
   do k = 0, 8 ; do j = 0, 8 ; do i = 0, 8
      p = vector_R8P(-0.2_R8P + 0.175_R8P * i, -0.2_R8P + 0.175_R8P * j, -0.2_R8P + 0.175_R8P * k)
      if ((p%x - 0.25_R8P)**2 + (p%y - 0.25_R8P)**2 + p%z**2 < 0.36_R8P) cycle
      if (abs(a%distance(point=p, is_signed=.true., is_square_root=.true.) - &
              b%distance(point=p, is_signed=.true., is_square_root=.true.)) > 1.e-12_R8P) same = .false.
   enddo ; enddo ; enddo
   endfunction same_far_distances
endprogram fossil_test_fill_holes
