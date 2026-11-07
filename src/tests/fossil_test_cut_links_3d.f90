!< FOSSIL test: the links between the 3-D cut-cell pieces (issue #70, step 3).
program fossil_test_cut_links_3d
!< FOSSIL test: the links between the 3-D cut-cell pieces (issue #70, step 3).
!<
!< Bodies: the 5-degree wedge tilted off every axis (thin cells with two pieces), the cone, the L-prism on a generic and
!< on a grid-aligned grid, the thin plate (two pieces per cell), a square tube whose hole is a fluid island of a face
!< (no fluid edge interval: matched by the face area left to the piece), and 8 lattice tetrahedra (vertices on the grid
!< nodes).
!< 1. per piece and face, the links sum to the piece's face fraction times the face area (relative 1e-13), and a face
!<    on the block boundary has no link;
!< 2. symmetry: p -> q over a face exists exactly when q -> p does, with the same area;
!< 3. every fluid region of a face is matched (unmatched = 0);
!< 4. the thin plate: a piece below the plate never links to a piece above it.

use fossil,                   only : surface_stl_object, cut_links_t, pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_facet_object,      only : facet_object
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable       :: b(:)
type(surface_stl_object)        :: s
type(cut_cells_3d_t)            :: cut
type(cut_links_t)               :: lk
type(facet_object), allocatable :: f(:)
real(R8P), allocatable          :: x(:), y(:), z(:)
real(R8P)                       :: emax(2), v(3, 4), vol
integer(I4P)                    :: i, k, nbad, nlinks, nunm, nmulti
integer(I8P)                    :: seed
logical                         :: ok

emax = 0._R8P ; nbad = 0 ; nlinks = 0 ; nunm = 0 ; nmulti = 0
call bench_bodies(b)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call check3(s, x, y, z)
call to_surface(b(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call check3(s, x, y, z)
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [-0.2_R8P, 0.45_R8P, 1.3_R8P]
call check3(s, x, y, z)
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call check3(s, x, y, z) ! aligned
call to_surface(b(5)%soup, s) ! plate [0, 1]^2 x [0, 1e-3]
x = [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P, 1.3_R8P] ; y = [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P] ; z = [-0.5_R8P, 0.5_R8P, 0.8_R8P]
call check3(s, x, y, z)
call plate3
call to_surface(tube([0._R8P, 3._R8P], [1._R8P, 2._R8P], [0._R8P, 3._R8P]), s) ! the hole is an island of the face z = 2.5
call check3(s, [-0.5_R8P, 2.5_R8P, 3.5_R8P], [-0.5_R8P, 2.5_R8P, 3.5_R8P], [0.5_R8P, 2.5_R8P, 4._R8P])
if (cut%pieces(1, 1, 1) /= 2 .or. size(lk%piece) == 0) nbad = nbad + 1
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)]
do k = 1, 8
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(s, v)
   call check3(s, x, x, x)
enddo
print '(A,I0,A,I0,A,I0,A,I0)', '3-D links ', nlinks, '  multi-piece cells ', nmulti, '  unmatched ', nunm, '  bad ', nbad
print '(A,2ES10.2)', 'max link-sum / symmetry errors: ', emax
ok = nbad == 0 .and. nunm == 0 .and. nlinks > 0 .and. nmulti > 20 .and. all(emax <= 1.e-13_R8P)
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check3(sf, xx, yy, zz)
   !< Checks 1-3 on the block xx x yy x zz.
   type(surface_stl_object), intent(in) :: sf
   real(R8P),                intent(in) :: xx(:), yy(:), zz(:)
   integer(I4P)                         :: p, fc, kk, q, k2, ic(3), n(3), ax
   real(R8P)                            :: sm, af(3), h(3)
   logical                              :: border, found

   call cut_cells_3d(sf, xx, yy, zz, cut, links=lk)
   n = [cut%nx, cut%ny, cut%nz]
   nlinks = nlinks + size(lk%piece) ; nunm = nunm + lk%unmatched
   nmulti = nmulti + count(cut%pieces > 1)
   do p = 1, cut%np
      ic = cut%piece_cell(:, p)
      h = [xx(ic(1) + 1) - xx(ic(1)), yy(ic(2) + 1) - yy(ic(2)), zz(ic(3) + 1) - zz(ic(3))]
      af = [h(2) * h(3), h(3) * h(1), h(1) * h(2)]
      do fc = 1, 6
         ax = (fc + 1) / 2
         sm = 0._R8P
         do kk = lk%start(p), lk%start(p + 1) - 1
            if (lk%face(kk) == fc) sm = sm + lk%area(kk)
         enddo
         border = (mod(fc, 2) == 1 .and. ic(ax) == 1) .or. (mod(fc, 2) == 0 .and. ic(ax) == n(ax))
         if (border) then
            if (sm /= 0._R8P) nbad = nbad + 1
         elseif (.not. next_flagged(ic, fc)) then
            emax(1) = max(emax(1), abs(sm - cut%piece_face(fc, p) * af(ax)) / af(ax))
         endif
      enddo
      do kk = lk%start(p), lk%start(p + 1) - 1
         q = lk%piece(kk) ; found = .false.
         do k2 = lk%start(q), lk%start(q + 1) - 1
            if (lk%piece(k2) == p .and. lk%face(k2) == lk%face(kk) + merge(1, -1, mod(lk%face(kk), 2) == 1)) then
               found = .true. ; emax(2) = max(emax(2), abs(lk%area(k2) - lk%area(kk)) / af((lk%face(kk) + 1) / 2))
            endif
         enddo
         if (.not. found) nbad = nbad + 1
      enddo
   enddo
   endsubroutine check3

   logical function next_flagged(ic, fc)
   !< The cell across face fc of cell ic is untrusted (no pieces to link to).
   integer(I4P), intent(in) :: ic(3), fc
   integer(I4P)             :: jc(3)

   jc = ic ; jc((fc + 1) / 2) = jc((fc + 1) / 2) + merge(-1, 1, mod(fc, 2) == 1)
   next_flagged = cut%pieces(jc(1), jc(2), jc(3)) < 0
   endfunction next_flagged

   subroutine plate3
   !< Check 4 on the last block (the plate): pieces below z = 5e-4 never link to pieces above.
   integer(I4P) :: p, kk, q

   do p = 1, cut%np
      if (cut%pieces(cut%piece_cell(1, p), cut%piece_cell(2, p), cut%piece_cell(3, p)) /= 2) cycle
      do kk = lk%start(p), lk%start(p + 1) - 1
         q = lk%piece(kk)
         if (cut%pieces(cut%piece_cell(1, q), cut%piece_cell(2, q), cut%piece_cell(3, q)) /= 2) cycle
         if ((cut%piece_centroid(3, p) - 5.e-4_R8P) * (cut%piece_centroid(3, q) - 5.e-4_R8P) < 0._R8P) nbad = nbad + 1
      enddo
   enddo
   endsubroutine plate3

   subroutine tetra(sf, vv)
   !< Tetrahedron vv (positive orientation) with outward facets.
   type(surface_stl_object), intent(inout) :: sf
   real(R8P),                intent(in)    :: vv(3, 4)
   integer(I4P), parameter                 :: T(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4])
   integer(I4P)                            :: q, j

   if (allocated(f)) deallocate(f)
   allocate(f(4))
   do q = 1, 4
      do j = 1, 3
         f(q)%vertex(j) = vector_R8P(vv(1, T(j, q)), vv(2, T(j, q)), vv(3, T(j, q)))
      enddo
   enddo
   call sf%adopt_facets(facets=f)
   endsubroutine tetra

   function tube(o, h, zr) result(sp)
   !< Square tube: outer [o(1), o(2)]^2, hole [h(1), h(2)]^2, z in zr.
   real(R8P), intent(in) :: o(2), h(2), zr(2) !< Outer and hole spans, height.
   type(soup_t)          :: sp                !< Soup.
   real(R8P)             :: oc(2, 4), hc(2, 4), c(2), mid(2) !< Corners, centre, edge midpoint.
   integer(I4P)          :: k, k1, m          !< Work.

   oc = reshape([o(1), o(1), o(2), o(1), o(2), o(2), o(1), o(2)], [2, 4])
   hc = reshape([h(1), h(1), h(2), h(1), h(2), h(2), h(1), h(2)], [2, 4])
   c = 0.5_R8P * [o(1) + o(2), o(1) + o(2)]
   do k = 1, 4
      k1 = mod(k, 4) + 1
      mid = 0.5_R8P * (oc(:, k) + oc(:, k1))
      call add_quad(sp, [oc(:, k), zr(1)], [oc(:, k1), zr(1)], [oc(:, k1), zr(2)], [oc(:, k), zr(2)], [mid - c, 0._R8P])
      mid = 0.5_R8P * (hc(:, k) + hc(:, k1))
      call add_quad(sp, [hc(:, k), zr(1)], [hc(:, k1), zr(1)], [hc(:, k1), zr(2)], [hc(:, k), zr(2)], [c - mid, 0._R8P])
      do m = 1, 2
         call add_quad(sp, [oc(:, k), zr(m)], [oc(:, k1), zr(m)], [hc(:, k1), zr(m)], [hc(:, k), zr(m)], &
                       [0._R8P, 0._R8P, merge(-1._R8P, 1._R8P, m == 1)])
      enddo
   enddo
   endfunction tube

   subroutine add_quad(sp, p1, p2, p3, p4, nout)
   !< Append the quad p1 p2 p3 p4 (planar) as two triangles, oriented along nout.
   type(soup_t), intent(inout) :: sp                 !< Soup.
   real(R8P),    intent(in)    :: p1(3), p2(3), p3(3), p4(3), nout(3) !< Corners, outward direction.
   real(R8P)                   :: t(3, 3, 2), nq(3)  !< Triangles, normal.

   nq = [(p2(2) - p1(2)) * (p3(3) - p1(3)) - (p2(3) - p1(3)) * (p3(2) - p1(2)), &
         (p2(3) - p1(3)) * (p3(1) - p1(1)) - (p2(1) - p1(1)) * (p3(3) - p1(3)), &
         (p2(1) - p1(1)) * (p3(2) - p1(2)) - (p2(2) - p1(2)) * (p3(1) - p1(1))]
   if (dot_product(nq, nout) > 0._R8P) then
      t(:, 1, 1) = p1 ; t(:, 2, 1) = p2 ; t(:, 3, 1) = p3 ; t(:, 1, 2) = p1 ; t(:, 2, 2) = p3 ; t(:, 3, 2) = p4
   else
      t(:, 1, 1) = p1 ; t(:, 2, 1) = p3 ; t(:, 3, 1) = p2 ; t(:, 1, 2) = p1 ; t(:, 2, 2) = p4 ; t(:, 3, 2) = p3
   endif
   if (.not. allocated(sp%v)) allocate(sp%v(3, 3, 0))
   sp%v = reshape([reshape(sp%v, [9 * size(sp%v, 3)]), reshape(t, [18])], [3, 3, size(sp%v, 3) + 2])
   endsubroutine add_quad

   pure function cross(a, c) result(w)
   real(R8P), intent(in) :: a(3), c(3)
   real(R8P)             :: w(3)

   w = [a(2) * c(3) - a(3) * c(2), a(3) * c(1) - a(1) * c(3), a(1) * c(2) - a(2) * c(1)]
   endfunction cross

   function rnd() result(r)
   real(R8P) :: r

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_links_3d
