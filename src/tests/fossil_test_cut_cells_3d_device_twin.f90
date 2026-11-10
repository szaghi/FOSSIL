!< FOSSIL test: the device path of the 3-D cut cells equals the host path bit for bit (issue #100 step 3c).
program fossil_test_cut_cells_3d_device_twin
!< FOSSIL test: the device path of the 3-D cut cells equals the host path bit for bit (issue #100 step 3c).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, offset and
!< prefix-sum bugs of the device stages, not device codegen (the distance bench does that on the GPU, `cut cells 3-D
!< nodes device vs serial mismatches`).
!< 1. Node states and cell-edge fractions: `cut_cells_3d_nodes_device` (one thread per grid line: crossing count,
!<    offsets, sorted crossings, exact placement) against `cut_cells_3d(..., edges=.true.)`, node_in and edge_x/_y/_z
!<    bitwise, with and without the edges, on a tilted wedge, a cone, an L-prism on an aligned grid, lattice
!<    tetrahedra (vertices and edges on grid lines and nodes), the unit cube `cube-01.stl` on a grid through its faces
!<    and that cube without one facet (an open surface); the open dragon `dragon.stl` on a coarse grid.
!< 2. The cut cells: `cut_cells_3d_device` (sections, node states, cell totals, classification, piercings, the count
!<    and fill passes of the pieces on workspace slices, prefix sums) against `cut_cells_3d` without links, every output
!<    bitwise (`cut_cells_3d_mismatches`), with and without the edges, copied back and resident (then `update_host`),
!<    the resident cut overwritten by the next call (the writers release it), on the same bodies.
use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, &
                                     pose_from_axis_angle, facet_object
use fossil_cut_cells_compare, only : cut_cells_3d_mismatches
use fossil_cut_cells_3d,      only : cut_cells_3d_nodes_device
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t),     allocatable  :: bb(:)
type(surface_stl_object)       :: s
type(facet_object), allocatable :: fc(:)
real(R8P)                      :: v(3, 4), vol
real(R8P),        allocatable  :: x(:), y(:), z(:)
integer(I4P)                   :: k, total, i
integer(I8P)                   :: seed
logical                        :: ok

ok = .true. ; total = 0_I4P
call bench_bodies(bb)
call to_surface(bb(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call body('tilted wedge')
call to_surface(bb(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call body('cone')
call to_surface(bb(8)%soup, s) ! L-prism, aligned grid
x = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P, 2._R8P, 2.5_R8P] ; y = x ; z = [-0.5_R8P, 0._R8P, 0.5_R8P, 1._R8P, 1.5_R8P]
call body('L-prism, aligned grid')
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)] ; y = x ; z = x
do k = 1, 3
   do
      v = reshape([(real(int(9 * rnd()), R8P) / 8._R8P, i = 1, 12)], [3, 4])
      vol = dot_product(v(:, 2) - v(:, 1), cross3(v(:, 3) - v(:, 1), v(:, 4) - v(:, 1)))
      if (vol /= 0._R8P) exit
   enddo
   if (vol < 0._R8P) v(:, [3, 4]) = v(:, [4, 3])
   call tetra(v)
   call body('lattice tetrahedron')
enddo
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)] ; y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)] ; z = x
call body('cube-01, grid through its faces')
allocate(fc(s%get_facets_number() - 1))
do k = 1, size(fc)
   fc(k) = s%facet_at(k + 1)
enddo
call s%adopt_facets(facets=fc)
call body('cube-01 without one facet (open)')
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
x = [(-11.3_R8P + 1.7_R8P * i, i = 0, 12)] ; y = [(-4.6_R8P + 1.0_R8P * i, i = 0, 10)] ; z = [(4.9_R8P + 1.3_R8P * i, i = 0, 12)]
call body('dragon (open)')
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine body(label)
   !< The device stages on s against `cut_cells_3d` on the grid x, y, z.
   character(*), intent(in) :: label !< Case.
   type(cut_cells_3d_t)     :: ch, ce, cd !< Host cut cells (with edges, without), device path.
   logical,   allocatable   :: nd(:,:,:)                     !< Device path: node flags.
   real(R8P), allocatable   :: ex(:,:,:), ey(:,:,:), ez(:,:,:) !< Device path: edge fractions.
   integer(I4P)             :: m, st, st2, nm, nu            !< Work.

   m = 0_I4P
   call cut_cells_3d(s, x, y, z, ch, status=st, edges=.true.)
   call cut_cells_3d_nodes_device(s, x, y, z, nd, ex, ey, ez, edges=.true., status=st2)
   m = m + abs(st - st2) + count(ch%node_in .neqv. nd, kind=I4P) + rbits(ch%edge_x, ex) + rbits(ch%edge_y, ey) + &
       rbits(ch%edge_z, ez)
   call cut_cells_3d_nodes_device(s, x, y, z, nd, ex, ey, ez, status=st2)
   m = m + abs(st - st2) + count(ch%node_in .neqv. nd, kind=I4P)
   if (allocated(ex)) m = m + 1_I4P ! no edges asked: none returned
   ! the cut cells
   call cut_cells_3d(s, x, y, z, ce, status=st)
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, edges=.true.)
   m = m + abs(st - st2) + cut_cells_3d_mismatches(ch, cd)
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, resident=.true.)
   call cd%update_host
   m = m + abs(st - st2) + cut_cells_3d_mismatches(ce, cd)
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, edges=.true., resident=.true.) ! overwrites a resident cut
   call cd%update_host
   m = m + abs(st - st2) + cut_cells_3d_mismatches(ch, cd)
   call cd%exit_device
   call cut_cells_3d_device(s, x, y, z, cd, status=st2)
   m = m + abs(st - st2) + cut_cells_3d_mismatches(ce, cd)
   nm = count(ch%pieces > 1) ; nu = count(ch%pieces < 0)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', label//': nodes ', size(nd), ', cells ', size(ch%alpha), ' (', nm, &
         ' multi-piece, ', nu, ' untrusted), device path mismatches ', m
   endsubroutine body

   function rbits(a, b) result(n)
   !< Elements of a and b with different bits (or shapes).
   real(R8P), intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.
   integer(I4P)          :: n                  !< Mismatches.

   if (any(shape(a) /= shape(b))) then
      n = size(a, kind=I4P) + 1_I4P ; return
   endif
   n = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   endfunction rbits

   subroutine tetra(vt)
   !< s = the tetrahedron vt (positively oriented).
   real(R8P), intent(in) :: vt(3, 4) !< Vertices.
   type(facet_object), allocatable :: f(:) !< Facets.
   integer(I4P), parameter :: T(3, 4) = reshape([1, 3, 2, 1, 2, 4, 1, 4, 3, 2, 3, 4], [3, 4]) !< Outward faces.
   integer(I4P)          :: q, r     !< Counters.

   allocate(f(4))
   do q = 1, 4
      do r = 1, 3
         f(q)%vertex(r) = vector_R8P(vt(1, T(r, q)), vt(2, T(r, q)), vt(3, T(r, q)))
      enddo
   enddo
   call s%adopt_facets(facets=f)
   endsubroutine tetra

   pure function cross3(a, b) result(c)
   !< a x b.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< Product.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross3

   function rnd() result(r)
   !< Deterministic uniform in (0, 1) (Park-Miller, the generator of `fossil_test_cut_links_3d`).
   real(R8P) :: r !< Value.

   seed = mod(seed * 48271_I8P, 2147483647_I8P)
   r = real(seed, R8P) / 2147483647._R8P
   endfunction rnd
endprogram fossil_test_cut_cells_3d_device_twin
