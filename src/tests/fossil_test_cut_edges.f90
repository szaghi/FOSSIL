!< FOSSIL test: the fluid fractions of the 3-D cell edges and the node flags (issue #70, step 4).
program fossil_test_cut_edges
!< FOSSIL test: the fluid fractions of the 3-D cell edges and the node flags (issue #70, step 4).
!<
!< `cut_cells_3d(..., edges=.true.)` gives edge_x/y/z (fluid fraction of every cell edge) and node_in (exact inside flag of the perturbed
!< grid nodes), from the exact crossings of every grid line. They are checked against an independent path: the 2-D cut
!< cells of the plane sections (`slice_surface` + `curve%cut_cells`), whose apertures are the fluid fractions of the
!< same edges and whose node flags are the same nodes. Bodies: the 5-degree wedge tilted off every axis, the L-prism
!< on a generic and on a grid-aligned grid, cube-01 on a grid through its faces, and 6 lattice tetrahedra.
!< 1. every edge fraction equals the aperture of that edge in both sections that contain it (1e-14);
!< 2. every node flag equals the node flag of the three sections through it, except at a node on the surface (distance
!<    below 1e-13), where the sections decide on rounded crossing points (cut-cells limitations: 7 such nodes, all on
!<    lattice tetrahedra; the 3-D flags are exact predicates on the input triangles);
!< 3. fractions in [0, 1], and an edge with fraction 1 between two inside nodes never occurs.

use fossil,                   only : surface_stl_object, curve_object, cut_cells_2d_t, slice_surface, pose_from_axis_angle
use fossil_cut_cells_3d,      only : cut_cells_3d_t, cut_cells_3d
use fossil_facet_object,      only : facet_object
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable       :: b(:)
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:)
real(R8P), allocatable          :: x(:), y(:), z(:)
real(R8P)                       :: emax, v(3, 4), vol
integer(I4P)                    :: i, k, nnode, nbad, nedge, nsurf
integer(I8P)                    :: seed
logical                         :: ok

emax = 0._R8P ; nnode = 0 ; nbad = 0 ; nedge = 0 ; nsurf = 0
call bench_bodies(b)
call to_surface(b(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call check(s, x, y, z)
call to_surface(b(8)%soup, s) ! L-prism
x = [(-0.3_R8P + 0.37_R8P * i, i = 0, 7)] ; y = x ; z = [-0.2_R8P, 0.45_R8P, 1.3_R8P]
call check(s, x, y, z)
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call check(s, x, y, z) ! aligned
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.5_R8P + 0.25_R8P * i, i = 0, 8)]
call check(s, x, x, x)
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)]
do k = 1, 6
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(s, v)
   call check(s, x, x, x)
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'edges checked ', nedge, '  nodes checked ', nnode, '  mismatches ', nbad, &
      '  node flags differing on the surface ', nsurf
print '(A,ES10.2)', 'max |edge fraction - section aperture|: ', emax
ok = nbad == 0 .and. emax <= 1.e-14_R8P .and. nedge > 0
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine check(sf, xx, yy, zz)
   !< Checks 1-3 on the block xx x yy x zz.
   type(surface_stl_object), intent(in) :: sf
   real(R8P),                intent(in) :: xx(:), yy(:), zz(:)
   type(cut_cells_3d_t)                 :: cut
   type(curve_object)                   :: sec
   type(cut_cells_2d_t)                 :: c2
   integer(I4P)                         :: ia, ja, ka, nx, ny, nz

   call cut_cells_3d(sf, xx, yy, zz, cut, edges=.true.)
   nx = cut%nx ; ny = cut%ny ; nz = cut%nz
   if (any(cut%edge_x < 0._R8P) .or. any(cut%edge_x > 1._R8P) .or. any(cut%edge_y < 0._R8P) .or. &
       any(cut%edge_y > 1._R8P) .or. any(cut%edge_z < 0._R8P) .or. any(cut%edge_z > 1._R8P)) nbad = nbad + 1
   do ia = 1, nx + 1 ! planes x: (u, v) = (y, z), cells (ja, ka)
      call slice_surface(surface=sf, normal_axis=1_I4P, level=xx(ia), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=yy, y=zz, cut=c2)
      do ka = 1, nz
         do ja = 1, ny
            call cmp(c2, ja, ka, [cut%edge_z(ia, ja, ka), cut%edge_z(ia, ja + 1, ka), cut%edge_y(ia, ja, ka), &
                     cut%edge_y(ia, ja, ka + 1)])
         enddo
      enddo
      call nodes(c2%node_in, cut%node_in(ia, :, :), 1_I4P, xx(ia), yy, zz)
   enddo
   do ja = 1, ny + 1 ! planes y: (u, v) = (z, x), cells (ka, ia)
      call slice_surface(surface=sf, normal_axis=2_I4P, level=yy(ja), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=zz, y=xx, cut=c2)
      do ia = 1, nx
         do ka = 1, nz
            call cmp(c2, ka, ia, [cut%edge_x(ia, ja, ka), cut%edge_x(ia, ja, ka + 1), cut%edge_z(ia, ja, ka), &
                     cut%edge_z(ia + 1, ja, ka)])
         enddo
      enddo
      call nodes(c2%node_in, transpose(cut%node_in(:, ja, :)), 2_I4P, yy(ja), zz, xx)
   enddo
   do ka = 1, nz + 1 ! planes z: (u, v) = (x, y), cells (ia, ja)
      call slice_surface(surface=sf, normal_axis=3_I4P, level=zz(ka), curve=sec, on_plane_low=.true.)
      call sec%cut_cells(x=xx, y=yy, cut=c2)
      do ja = 1, ny
         do ia = 1, nx
            call cmp(c2, ia, ja, [cut%edge_y(ia, ja, ka), cut%edge_y(ia + 1, ja, ka), cut%edge_x(ia, ja, ka), &
                     cut%edge_x(ia, ja + 1, ka)])
         enddo
      enddo
      call nodes(c2%node_in, cut%node_in(:, :, ka), 3_I4P, zz(ka), xx, yy)
   enddo
   do ka = 1, nz + 1 ! check 3 along x
      do ja = 1, ny + 1
         do ia = 1, nx
            if (cut%edge_x(ia, ja, ka) == 1._R8P .and. cut%node_in(ia, ja, ka) .and. cut%node_in(ia + 1, ja, ka)) &
               nbad = nbad + 1
         enddo
      enddo
   enddo
   endsubroutine check

   subroutine cmp(c2, iu, iv, e3)
   !< The four apertures of section cell (iu, iv) against the 3-D edge fractions e3 (-u, +u, -v, +v).
   type(cut_cells_2d_t), intent(in) :: c2
   integer(I4P),         intent(in) :: iu, iv
   real(R8P),            intent(in) :: e3(4)

   if (.not. c2%closed(iu, iv)) return
   nedge = nedge + 4
   emax = max(emax, maxval(abs(c2%aperture(:, iu, iv) - e3)))
   endsubroutine cmp

   subroutine nodes(n2, n3, ax, lev, uu, vv)
   !< Node flags of a section (plane ax at lev, in-plane coordinates uu, vv) against the 3-D ones on its plane; a
   !< difference is allowed only at a node on the surface.
   logical,      intent(in) :: n2(:,:), n3(:,:)
   integer(I4P), intent(in) :: ax
   real(R8P),    intent(in) :: lev, uu(:), vv(:)
   integer(I4P)             :: iu, iv
   real(R8P)                :: p(3)

   nnode = nnode + size(n2)
   do iv = 1, size(n2, 2)
      do iu = 1, size(n2, 1)
         if (n2(iu, iv) .eqv. n3(iu, iv)) cycle
         p(ax) = lev ; p(mod(ax, 3) + 1) = uu(iu) ; p(mod(ax + 1, 3) + 1) = vv(iv)
         if (s%distance(point=vector_R8P(p(1), p(2), p(3)), is_signed=.false., is_square_root=.true.) < 1.e-13_R8P) then
            nsurf = nsurf + 1
         else
            nbad = nbad + 1
         endif
      enddo
   enddo
   endsubroutine nodes

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
endprogram fossil_test_cut_edges
