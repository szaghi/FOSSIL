!< FOSSIL test: isotropic_remesh coarsens (issue #79).
!<
!< Before #79 the collapse pass rejected almost every collapse: its 1-ring check counted the two apexes of the edge's
!< own facets, which neighbour both endpoints on every closed mesh, so a larger target length never reduced the facet
!< count. The check is now the link condition (the endpoints share no neighbour but the two apexes).
!< 1. cube-01 at L = 1, features on: the edges and corners are locked, the 6 face centres collapse onto them (36
!<    facets, 48 before #79), area and volume stay 6 and 1.
!< 2. The dragon at L = 1.5, 3 iterations: far fewer than its 6588 facets, no new boundary or non-manifold edge.
!< 3. A marching-cubes sphere (25^3, as on isotropic-remesh.md) at L = 0.1, one iteration per call: the number of edges
!<    shorter than 4L/5 falls from call to call (before #79 it grew: 696, 2017, 2382).
!< 4. The same sphere at L = 0.4, 3 iterations: coarser (2312 -> 298 facets), watertight, manifold, certified one-sided
!<    deviation from the input 4.14e-2 <= h <= 4.22e-2 at tolerance 1e-3 when written (bound 4.5e-2). Edges near 0.4 on a
!<    unit sphere: the chords alone put the facet centres about L^2/6 = 2.7e-2 inside.
program fossil_test_remesh_coarsen

use fossil, only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, REM_STATUS_OK, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)        :: s, input
type(facet_object), allocatable :: facets(:)
real(R8P)                       :: phi(25, 25, 25), lower, upper
integer(I4P)                    :: i, j, k, status, nbe0, nnm0, short(0:3)
logical                         :: ok(4)
real(R8P), parameter            :: HSPHERE = 4.5e-2_R8P !< Deviation bound of case 4 (measured 4.14e-2 .. 4.22e-2).

ok = .false.

! 1. cube-01, features on
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call s%isotropic_remesh(target_length=1._R8P, iterations=3_I4P, preserve_features=.true., status=status)
print '(A,I0,A,F12.9,A,F12.9,A,L1)', '1. cube L=1, features on: facets ', s%get_facets_number(), ', area ', s%get_area(), &
      ', volume ', s%get_volume(), ', watertight ', s%is_watertight()
ok(1) = status == REM_STATUS_OK .and. s%get_facets_number() == 36_I4P .and. s%is_watertight() .and. s%is_manifold() .and. &
        abs(s%get_area() - 6._R8P) < 1.e-12_R8P .and. abs(s%get_volume() - 1._R8P) < 1.e-12_R8P

! 2. dragon at L = 1.5
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call edge_classes(s, nbe0, nnm0)
call s%isotropic_remesh(target_length=1.5_R8P, iterations=3_I4P, preserve_features=.false., status=status)
call edge_classes(s, i, j)
print '(A,I0,A,I0,A,I0,A,I0,A,I0,A)', '2. dragon L=1.5: facets 6588 -> ', s%get_facets_number(), ', boundary edges ', nbe0, &
      ' -> ', i, ', non-manifold edges ', nnm0, ' -> ', j
ok(2) = status == REM_STATUS_OK .and. s%get_facets_number() < 6588_I4P / 2 .and. i <= nbe0 .and. j <= nnm0

! 3. sphere at L = 0.1, short edges per call
do k = 1, 25 ; do j = 1, 25 ; do i = 1, 25
   phi(i, j, k) = norm2(-1.5_R8P + 0.125_R8P * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call input%adopt_facets(facets=facets)
call input%sanitize(fill_holes=.true.)
s = input
short(0) = short_edges(s, 0.08_R8P)
do k = 1, 3
   call s%isotropic_remesh(target_length=0.1_R8P, iterations=1_I4P, status=status)
   short(k) = short_edges(s, 0.08_R8P)
enddo
print '(A,4I6,A,I0,A,L1)', '3. sphere L=0.1, edges < 4L/5 per call:', short, ', facets ', s%get_facets_number(), &
      ', watertight ', s%is_watertight()
ok(3) = short(2) < short(1) .and. short(3) < short(2) .and. short(3) < short(0) .and. s%is_watertight() .and. s%is_manifold()

! 4. sphere at L = 0.4: certified deviation from the input
s = input
call s%isotropic_remesh(target_length=0.4_R8P, iterations=3_I4P, status=status)
call s%deviation(input, lower=lower, upper=upper, tolerance=1.e-3_R8P, status=k)
print '(A,I0,A,I0,A,L1,A,L1,A,2ES12.4)', '4. sphere L=0.4: facets ', input%get_facets_number(), ' -> ', &
      s%get_facets_number(), ', watertight ', s%is_watertight(), ', manifold ', s%is_manifold(), ', deviation ', lower, upper
ok(4) = status == REM_STATUS_OK .and. k == DEV_STATUS_OK .and. s%get_facets_number() < input%get_facets_number() / 4 .and. &
        s%is_watertight() .and. s%is_manifold() .and. upper <= HSPHERE

print '(A,4L2)', 'per-case results: ', ok
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains

   function short_edges(surf, thr) result(n)
   !< Number of edges shorter than `thr` (each interior edge seen from its two facets).
   type(surface_stl_object), intent(in) :: surf
   real(R8P),                intent(in) :: thr
   integer(I4P)                         :: n
   type(facet_object), pointer          :: f
   integer(I4P)                         :: m, e

   n = 0_I4P
   do m = 1, surf%get_facets_number()
      f => surf%facet_at(m)
      do e = 1, 3
         if (norm2([f%vertex(e)%x - f%vertex(mod(e, 3) + 1)%x, f%vertex(e)%y - f%vertex(mod(e, 3) + 1)%y, &
                    f%vertex(e)%z - f%vertex(mod(e, 3) + 1)%z]) < thr) n = n + 1
      enddo
   enddo
   n = n / 2
   endfunction short_edges

   subroutine edge_classes(surf, nboundary, nnonmanifold)
   !< Boundary edges (one facet) and non-manifold edges (three or more), by vertex ids.
   type(surface_stl_object), intent(in)  :: surf
   integer(I4P),             intent(out) :: nboundary, nnonmanifold
   type(facet_object), pointer           :: f
   integer(I4P), allocatable             :: key(:, :)
   integer(I4P)                          :: m, e, h, a, b, run

   allocate(key(2, 3 * surf%get_facets_number()))
   h = 0
   do m = 1, surf%get_facets_number()
      f => surf%facet_at(m)
      do e = 1, 3
         a = f%vertex_id(e) ; b = f%vertex_id(mod(e, 3) + 1)
         h = h + 1 ; key(:, h) = [min(a, b), max(a, b)]
      enddo
   enddo
   call sort_keys(key)
   nboundary = 0 ; nnonmanifold = 0 ; run = 1
   do m = 2, h + 1
      if (m <= h) then
         if (all(key(:, m) == key(:, m - 1))) then
            run = run + 1 ; cycle
         endif
      endif
      if (run == 1) nboundary = nboundary + 1
      if (run >= 3) nnonmanifold = nnonmanifold + 1
      run = 1
   enddo
   endsubroutine edge_classes

   subroutine sort_keys(key)
   !< Sort (a, b) pairs lexicographically (shell sort: test helper).
   integer(I4P), intent(inout) :: key(:, :)
   integer(I4P)                :: gap, i, j, t(2)

   gap = size(key, dim=2) / 2
   do while (gap > 0)
      do i = gap + 1, size(key, dim=2)
         t = key(:, i) ; j = i
         do while (j > gap)
            if (key(1, j - gap) < t(1) .or. (key(1, j - gap) == t(1) .and. key(2, j - gap) <= t(2))) exit
            key(:, j) = key(:, j - gap) ; j = j - gap
         enddo
         key(:, j) = t
      enddo
      gap = gap / 2
   enddo
   endsubroutine sort_keys

endprogram fossil_test_remesh_coarsen
