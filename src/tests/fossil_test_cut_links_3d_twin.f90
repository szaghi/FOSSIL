!< FOSSIL test: the 3-D piece links on the device equal the host's bit for bit (issue #100 step 5).
program fossil_test_cut_links_3d_twin
!< FOSSIL test: the 3-D piece links on the device equal the host's bit for bit (issue #100 step 5).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, offset and
!< prefix-sum bugs of the device stages, not device codegen (the distance bench does that on the GPU, `cut links 3-D
!< device vs serial mismatches`).
!< `cut_cells_3d_device(..., links=)` against `cut_cells_3d(..., links=)`: the cut (`cut_cells_3d_mismatches`) and the
!< links (`cut_links_mismatches`: np, unmatched, start, piece, face, area bitwise), copied back and resident (then
!< `update_host`; the next call releases them) without the edges, resident with them (the 2 s budget); and the stage
!< `cut_links_3d_device` on link inputs copied back to the host (the copy-in path). A copy of resident links is not
!< resident: its `exit_device` leaves the original mapped. Cases: the tilted 5-degree wedge (thin cells with two pieces),
!< the thin plate (two pieces per cell), the cone, two lattice tetrahedra (vertices on grid nodes) and the unit cube
!< without one facet (an open surface: untrusted cells).
use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, cut_links_t, &
                                     pose_from_axis_angle, facet_object
use fossil_cut_cells_compare, only : cut_cells_3d_mismatches, cut_links_mismatches
use fossil_cut_cells_3d,      only : cut_link_inputs_3d_t, cut_cells_3d_link_inputs_device, cut_links_3d_device
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t),     allocatable  :: bb(:)
type(surface_stl_object)       :: s
type(facet_object), allocatable :: fc(:)
real(R8P)                      :: v(3, 4), vol
real(R8P),        allocatable  :: x(:), y(:), z(:)
integer(I4P)                   :: k, total, i, nlinks, nmulti
integer(I8P)                   :: seed
logical                        :: ok

ok = .true. ; total = 0_I4P ; nlinks = 0_I4P ; nmulti = 0_I4P
call bench_bodies(bb)
call to_surface(bb(2)%soup, s) ! wedge 5, tilted
call s%move(pose_from_axis_angle([0.3_R8P, -0.5_R8P, 0.81_R8P], 0.61_R8P))
x = [(-0.53_R8P + 0.071_R8P * i, i = 0, 22)] ; y = [(-0.67_R8P + 0.063_R8P * i, i = 0, 24)]
z = [(-0.59_R8P + 0.077_R8P * i, i = 0, 18)]
call body('tilted wedge')
call to_surface(bb(5)%soup, s) ! plate [0, 1]^2 x [0, 1e-3]
x = [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P, 1.3_R8P] ; y = [0.1_R8P, 0.4_R8P, 0.7_R8P, 0.9_R8P] ; z = [-0.5_R8P, 0.5_R8P, 0.8_R8P]
call body('thin plate')
call to_surface(bb(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call body('cone')
seed = 70_I8P
x = [(-0.125_R8P + 0.125_R8P * (i - 1), i = 1, 11)] ; y = x ; z = x
do k = 1, 2
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
allocate(fc(s%get_facets_number() - 1))
do k = 1, size(fc)
   fc(k) = s%facet_at(k + 1)
enddo
call s%adopt_facets(facets=fc)
x = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)] ; y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)] ; z = x
call body('cube-01 without one facet (open)')
if (nlinks == 0_I4P .or. nmulti == 0_I4P) ok = .false. ! the cases must exercise the links and multi-piece faces
print '(A,I0,A,I0,A,I0)', 'total mismatches: ', total, ', links: ', nlinks, ', multi-piece cells: ', nmulti
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine body(label)
   !< The links of the device path on s against the host's on the grid x, y, z.
   character(*), intent(in)   :: label    !< Case.
   type(cut_cells_3d_t)       :: ch, cd   !< Host, device cut.
   type(cut_links_t)          :: lh, ld, lb !< Host, device links, a copy of resident ones.
   type(cut_link_inputs_3d_t) :: li       !< Link inputs (copied back).
   integer(I4P)               :: m, st, st2, r, e !< Work.

   m = 0_I4P
   do e = 1, 2
      call cut_cells_3d(s, x, y, z, ch, status=st, edges=e == 2, links=lh)
      do r = e, 2 ! with the edges resident only
         call cut_cells_3d_device(s, x, y, z, cd, status=st2, edges=e == 2, resident=r == 2, links=ld)
         if (r == 2) then
            if (.not. (ld%on_device .and. cd%on_device)) m = m + 1_I4P
            lb = ld ; call lb%exit_device ! a copy: the original stays mapped
            if (lb%on_device .or. .not. ld%on_device) m = m + 1_I4P
            call cd%update_host ; call ld%update_host
         elseif (ld%on_device .or. cd%on_device) then
            m = m + 1_I4P
         endif
         m = m + abs(st - st2) + cut_cells_3d_mismatches(ch, cd) + cut_links_mismatches(lh, ld)
      enddo
   enddo
   call ld%exit_device ; call cd%exit_device
   if (ld%on_device) m = m + 1_I4P
   ! the stage on host inputs (copied in for the call)
   call cut_cells_3d_link_inputs_device(s, x, y, z, cd, li, status=st2)
   call cut_links_3d_device(cd, li, x, y, z, ld)
   m = m + abs(st - st2) + cut_links_mismatches(lh, ld)
   nlinks = nlinks + size(lh%piece, kind=I4P) ; nmulti = nmulti + count(ch%pieces > 1, kind=I4P)
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', label//': cells ', size(ch%alpha), ' (', count(ch%pieces > 1), ' multi-piece), '// &
         'links ', size(lh%piece), ', unmatched ', lh%unmatched, ', device path mismatches ', m
   endsubroutine body

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
endprogram fossil_test_cut_links_3d_twin
