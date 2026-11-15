!< FOSSIL test: 3-D cut cells of a facet within ulps of a grid plane (issue #107).
program fossil_test_cut_ulp_plane
!< FOSSIL test: 3-D cut cells of a facet within ulps of a grid plane (issue #107).
!<
!< A hexahedron [0.1, 0.9]^2 x [zk, zk + 0.3] (or below zk) on the 4^3 grid of spacing 0.25 over [0, 1]^3 (levels exact),
!< zk = 0.5 a grid level. The face on zk is two triangles (diagonal v1 v3) whose four vertices sit at zk + d_m ulp(zk):
!< some exactly on the plane, the others a few ulps off it (or 1000 ulps: a shallow genuine crossing). The face crosses
!< the x and y grid planes 0.25, 0.5, 0.75. The other axes are reached by a cyclic permutation of the coordinates; the
!< shifted cases add 1.3 to every grid level and body coordinate (levels and vertices no longer dyadic: the #106 case).
!< Every cell's solid volume (1 - alpha) h^3 is compared with its overlap with the flat box, to 1e-13 h^3 plus the
!< offsets' own share (max |d| ulp(zk) h^2: <= 4.4e-16 h^3 at one ulp); the total with the polyhedron's volume
!< (divergence theorem about the centre, error ~1e-16) to 1e-13. Paths: `cut_cells_3d`, with links, and
!< `cut_cells_3d_device` (the host twin without OpenACC).
!<
!< The defect (fixed by #107 step 1): the cut cells clipped every facet by x, then y, then z planes, and a crossing made
!< by an x (or y) clip interpolates the facet's z, which rounds onto zk for a facet within ulps of it; the z clip then
!< put that part of the facet in the cell below, while the plane section (`facet_segment_dev`, the facet's own
!< vertices) puts it above: the cell below got a wall polygon on its top face that no face term closes, its volume off
!< by up to 1/3 of the cell (1000 ulps: 2.7e-4). A nearly flat facet is now clipped along its flat axis first
!< (`flat_axis`, `clip_cell`). The cases with vertices one ulp below, exactly on the plane and on x planes were right
!< before too (controls). Any failure stops the program.

use fossil,           only : surface_stl_object, facet_object, cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device
use fossil_cut_graph, only : cut_links_t
use penf,             only : I4P, R8P
use vecfor,           only : vector_R8P

implicit none

integer(I4P), parameter :: UP2(4) = [0, 0, 1, 1], ALT(4) = [0, 1, 0, 1], STR(4) = [-1, 0, 1, 0] !< Face offsets [ulp].
integer(I4P), parameter :: STR1K(4) = [-1000, 0, 1000, 0], DN2(4) = [0, 0, -1, -1], ON(4) = 0      !< Face offsets [ulp].
logical                 :: ok !< All cases right.

ok = .true.
ok = check('z, [0, 0, 1, 1] ulp, host', 3_I4P, UP2, 1._R8P, 0_I4P, 0._R8P) .and. ok
ok = check('z, [0, 0, 1, 1] ulp, host + links', 3_I4P, UP2, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('z, [0, 0, 1, 1] ulp, device path', 3_I4P, UP2, 1._R8P, 2_I4P, 0._R8P) .and. ok
ok = check('z, [0, 1, 0, 1] ulp, host + links', 3_I4P, ALT, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('z, [-1, 0, 1, 0] ulp (straddle), host', 3_I4P, STR, 1._R8P, 0_I4P, 0._R8P) .and. ok
ok = check('z, [-1000, 0, 1000, 0] ulp (straddle), host + links', 3_I4P, STR1K, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('z, [-1000, 0, 1000, 0] ulp (straddle), device path', 3_I4P, STR1K, 1._R8P, 2_I4P, 0._R8P) .and. ok
ok = check('y, [-1000, 0, 1000, 0] ulp (straddle), host + links', 2_I4P, STR1K, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('z, [0, 0, 1, 1] ulp, body below zk', 3_I4P, UP2, -1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('y, [0, 0, 1, 1] ulp, host + links', 2_I4P, UP2, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('y, [-1, 0, 1, 0] ulp (straddle), host + links', 2_I4P, STR, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('z, [0, 0, 1, 1] ulp, shifted 1.3, host + links', 3_I4P, UP2, 1._R8P, 1_I4P, 1.3_R8P) .and. ok
ok = check('z, [-1, 0, 1, 0] ulp (straddle), shifted 1.3, device path', 3_I4P, STR, 1._R8P, 2_I4P, 1.3_R8P) .and. ok
ok = check('y, [0, 1, 0, 1] ulp, shifted 1.3, host + links', 2_I4P, ALT, 1._R8P, 1_I4P, 1.3_R8P) .and. ok
ok = check('z, [0, 0, 1, 1] ulp, body below zk, shifted 1.3', 3_I4P, UP2, -1._R8P, 1_I4P, 1.3_R8P) .and. ok
ok = check('control z, [0, 0, -1, -1] ulp (below)', 3_I4P, DN2, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('control z, [0, 0, 0, 0] (on the plane)', 3_I4P, ON, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('control x, [0, 0, 1, 1] ulp', 1_I4P, UP2, 1._R8P, 1_I4P, 0._R8P) .and. ok
ok = check('control x, [-1, 0, 1, 0] ulp (straddle)', 1_I4P, STR, 1._R8P, 2_I4P, 0._R8P) .and. ok
ok = check('control x, [0, 0, 1, 1] ulp, shifted 1.3', 1_I4P, UP2, 1._R8P, 1_I4P, 1.3_R8P) .and. ok
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   function check(name, ax, d, side, path, shift) result(good)
   !< Cut the hexahedron whose face near zk has vertex offsets d (ulps) and check every cell and the total.
   character(*), intent(in) :: name !< Case.
   integer(I4P), intent(in) :: ax   !< Axis normal to the face near the plane (1, 2, 3).
   integer(I4P), intent(in) :: d(4) !< Offsets of the face's vertices [ulp(zk)].
   real(R8P),    intent(in) :: side !< +1: the body is above zk, -1: below.
   integer(I4P), intent(in) :: path !< 0: cut_cells_3d, 1: with links, 2: cut_cells_3d_device.
   real(R8P),    intent(in) :: shift !< Added to every level and coordinate (0: none).
   logical                  :: good !< Every cell and the total right.
   type(facet_object), allocatable :: f(:)      !< Facets.
   type(surface_stl_object)        :: s         !< Body.
   type(cut_cells_3d_t)            :: cut       !< Cut cells.
   type(cut_links_t)               :: lk        !< Links.
   real(R8P)    :: v(3, 8), tri(3, 3, 12), p(3), g(5), zk, h, vex, vsum, ea, vc, ov, b0(3), b1(3), c0(3), c1(3) !< Work.
   real(R8P)    :: lo, hi               !< Body extent in x, y.
   real(R8P)    :: tol                  !< Cell tolerance (fraction of the cell).
   integer(I4P) :: perm(3), i, j, kk, m, q, a, bad !< Work.

   g = [(0.25_R8P * (i - 1) + shift, i = 1, 5)]
   zk = g(3) ; h = 0.3_R8P * side ; lo = 0.1_R8P + shift ; hi = 0.9_R8P + shift
   v(:, 1) = [lo, lo, zk + d(1) * spacing(zk)] ; v(:, 2) = [hi, lo, zk + d(2) * spacing(zk)]
   v(:, 3) = [hi, hi, zk + d(3) * spacing(zk)] ; v(:, 4) = [lo, hi, zk + d(4) * spacing(zk)]
   do q = 1, 4
      v(:, 4 + q) = [v(1, q), v(2, q), zk + h]
   enddo
   tol = 1.e-13_R8P + maxval(abs(d)) * spacing(zk) / (g(2) - g(1))
   call quad(1_I4P, v(:, 1), v(:, 4), v(:, 3), v(:, 2), tri)  ! the face near zk (-z for side = +1)
   call quad(3_I4P, v(:, 5), v(:, 6), v(:, 7), v(:, 8), tri)  ! the far face
   call quad(5_I4P, v(:, 1), v(:, 2), v(:, 6), v(:, 5), tri)  ! sides
   call quad(7_I4P, v(:, 2), v(:, 3), v(:, 7), v(:, 6), tri)
   call quad(9_I4P, v(:, 3), v(:, 4), v(:, 8), v(:, 7), tri)
   call quad(11_I4P, v(:, 4), v(:, 1), v(:, 5), v(:, 8), tri)
   if (side < 0._R8P) then ! the body below zk: reverse every facet
      do m = 1, 12
         p = tri(:, 2, m) ; tri(:, 2, m) = tri(:, 3, m) ; tri(:, 3, m) = p
      enddo
   endif
   ! cyclic permutation (orientation kept): coordinate a of the surface is coordinate perm(a) of the body above
   select case (ax)
   case (1)
      perm = [3_I4P, 1_I4P, 2_I4P]
   case (2)
      perm = [2_I4P, 3_I4P, 1_I4P]
   case default
      perm = [1_I4P, 2_I4P, 3_I4P]
   endselect
   allocate(f(12))
   vex = 0._R8P
   do m = 1, 12
      do q = 1, 3
         f(m)%vertex(q) = vector_R8P(tri(perm(1), q, m), tri(perm(2), q, m), tri(perm(3), q, m))
      enddo
      vex = vex + vol6(tri(:, 1, m) - g(3), tri(:, 2, m) - g(3), tri(:, 3, m) - g(3)) / 6._R8P
   enddo
   vex = abs(vex)
   call s%adopt_facets(facets=f)
   select case (path)
   case (1)
      call cut_cells_3d(surface=s, x=g, y=g, z=g, cut=cut, links=lk)
   case (2)
      call cut_cells_3d_device(surface=s, x=g, y=g, z=g, cut=cut)
   case default
      call cut_cells_3d(surface=s, x=g, y=g, z=g, cut=cut)
   endselect
   b0 = [lo, lo, min(zk, zk + h)] ; b1 = [hi, hi, max(zk, zk + h)]
   vsum = 0._R8P ; ea = 0._R8P ; bad = 0_I4P
   do kk = 1, 4
      do j = 1, 4
         do i = 1, 4
            c0 = [g(i), g(j), g(kk)] ; c1 = [g(i + 1), g(j + 1), g(kk + 1)] ; vc = product(c1 - c0)
            ov = 1._R8P
            do a = 1, 3
               ov = ov * max(0._R8P, min(c1(a), b1(perm(a))) - max(c0(a), b0(perm(a))))
            enddo
            vsum = vsum + (1._R8P - cut%alpha(i, j, kk)) * vc
            ea = max(ea, abs((1._R8P - cut%alpha(i, j, kk)) * vc - ov) / vc)
            if (abs((1._R8P - cut%alpha(i, j, kk)) * vc - ov) > tol * vc) bad = bad + 1_I4P
         enddo
      enddo
   enddo
   good = all(cut%closed) .and. bad == 0_I4P .and. abs(vsum - vex) <= 1.e-13_R8P * vex
   print '(A,A,ES10.2,A,ES10.2,A,I0,A,I0,A,L1)', name, ': volume err ', (vsum - vex) / vex, ', max cell err ', ea, &
         ', cells wrong ', bad, ', untrusted ', count(cut%pieces == -1_I4P), ', ok ', good
   endfunction check

   subroutine quad(m0, a, b, c, dd, tri)
   !< The quad a b c dd as the triangles (a, b, c), (a, c, dd) at tri(:, :, m0:m0+1).
   integer(I4P), intent(in)    :: m0                    !< First triangle.
   real(R8P),    intent(in)    :: a(3), b(3), c(3), dd(3) !< Vertices.
   real(R8P),    intent(inout) :: tri(3, 3, 12)          !< Triangles.

   tri(:, 1, m0) = a ; tri(:, 2, m0) = b ; tri(:, 3, m0) = c
   tri(:, 1, m0 + 1) = a ; tri(:, 2, m0 + 1) = c ; tri(:, 3, m0 + 1) = dd
   endsubroutine quad

   pure real(R8P) function vol6(a, b, c)
   !< Six times the signed volume of the tetrahedron (0, a, b, c).
   real(R8P), intent(in) :: a(3), b(3), c(3) !< Vertices.

   vol6 = a(1) * (b(2) * c(3) - b(3) * c(2)) - a(2) * (b(1) * c(3) - b(3) * c(1)) + a(3) * (b(1) * c(2) - b(2) * c(1))
   endfunction vol6
endprogram fossil_test_cut_ulp_plane
