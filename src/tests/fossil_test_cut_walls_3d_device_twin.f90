!< FOSSIL test: the labelled wall terms of the device 3-D cut cells equal the host's bit for bit (#109 WP2 step 05).
program fossil_test_cut_walls_3d_device_twin
!< FOSSIL test: the labelled wall terms of the device 3-D cut cells equal the host's bit for bit (#109 WP2 step 05).
!<
!< Built without OpenACC (this CI build) the device loops run on the host: the test catches the driver, offset and
!< prefix-sum bugs of the wall stage (count on the workspace slices or by the clip, `scan_dev`, fill at the offsets),
!< not device codegen (the distance bench does that on the GPU, `cut walls 3-D device vs serial mismatches`).
!< `cut_cells_3d_device(..., walls=)` against `cut_cells_3d(..., walls=)`, every field of `cut_walls_3d_t` bitwise
!< (`cut_walls_mismatches`: n, cell_start, piece, facet, area, wet_area, centroid) and the cut (`cut_cells_3d_mismatches`):
!< copied back and resident (then `update_host`), a resident result overwritten by the next call (the writers release
!< it), with links (copied back and resident), the alpha-only mode (the terms reset), the host path over resident terms.
!< Bodies: a marching-cubes sphere, the faceted cone, `cube-01.stl` on a grid through its faces, a tilted thin plate
!< (multi-piece cells: the terms read from the workspace slice the fill pass left), the open dragon (untrusted cells,
!< piece 0) and a grid that misses the body (no term: empty arrays).
use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_cells_3d_device, cut_walls_3d_t, &
                                     cut_links_t, pose_from_axis_angle
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_cut_cells_compare, only : cut_cells_3d_mismatches, cut_walls_mismatches, cut_links_mismatches
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(body_t), allocatable :: bb(:)
type(surface_stl_object)  :: s
real(R8P),    allocatable :: x(:), y(:), z(:)
integer(I4P)              :: total, i, nmulti
logical                   :: ok

ok = .true. ; total = 0_I4P ; nmulti = 0_I4P
call bench_bodies(bb)
call mc_body(BODY_SPHERE, 12_I4P, s) ; call box_grid(10_I4P)
call body('marching-cubes sphere')
call to_surface(bb(7)%soup, s) ! cone
x = [(-1.13_R8P + 0.19_R8P * i, i = 0, 12)] ; y = x ; z = [(-0.21_R8P + 0.23_R8P * i, i = 0, 18)]
call body('cone')
call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
x = [(-0.25_R8P + 0.125_R8P * i, i = 0, 12)] ; y = [(-0.3_R8P + 0.15_R8P * i, i = 0, 10)] ; z = x
call body('cube-01, grid through its faces')
x = [(2._R8P + 0.25_R8P * i, i = 0, 4)] ; y = x ; z = x
call body('cube-01, grid off the body')
call to_surface(bb(5)%soup, s) ! plate 1e-3, tilted
call s%move(pose_from_axis_angle([0.2_R8P, 0.1_R8P, 1._R8P], 0.3_R8P))
x = [(-0.53_R8P + 0.182_R8P * i, i = 0, 11)] ; y = x ; z = [(-0.0105_R8P + 0.0014_R8P * i, i = 0, 15)]
call body('thin plate, tilted', multi=.true.)
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
x = [(-11.3_R8P + 2.21_R8P * i, i = 0, 6)] ; y = [(-4.6_R8P + 1.3_R8P * i, i = 0, 6)] ; z = [(4.9_R8P + 1.69_R8P * i, i = 0, 6)]
call body('dragon (open)')
ok = ok .and. nmulti > 0_I4P
print '(A,I0)', 'total mismatches: ', total
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine body(label, multi)
   !< The device wall terms on s against `cut_cells_3d(..., walls=)` on the grid x, y, z.
   character(*), intent(in)           :: label !< Case.
   logical,      intent(in), optional :: multi !< The case must have multi-piece cells with terms.
   type(cut_cells_3d_t)               :: ch, cd !< Host, device cut.
   type(cut_walls_3d_t)               :: wh, wd !< Host, device wall terms.
   type(cut_links_t)                  :: lh, ld !< Host, device links.
   integer(I4P)                       :: m, st, st2, r, nm !< Work.

   m = 0_I4P
   call cut_cells_3d(s, x, y, z, ch, status=st, walls=wh)
   do r = 1, 2 ! copied back, resident
      call cut_cells_3d_device(s, x, y, z, cd, status=st2, resident=r == 2, walls=wd)
      if (r == 2) then
         if (.not. (wd%on_device .and. cd%on_device)) m = m + 1_I4P
         call cd%update_host ; call wd%update_host
      elseif (wd%on_device) then
         m = m + 1_I4P
      endif
      m = m + abs(st - st2) + cut_cells_3d_mismatches(ch, cd) + cut_walls_mismatches(wh, wd)
   enddo
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, walls=wd) ! over resident terms (released first)
   m = m + abs(st - st2) + cut_walls_mismatches(wh, wd)
   if (wd%on_device) m = m + 1_I4P
   ! with links: the same terms, copied back and resident
   call cut_cells_3d(s, x, y, z, ch, status=st, links=lh)
   do r = 1, 2
      call cut_cells_3d_device(s, x, y, z, cd, status=st2, resident=r == 2, links=ld, walls=wd)
      if (r == 2) then
         if (.not. (wd%on_device .and. ld%on_device .and. cd%on_device)) m = m + 1_I4P
         call cd%update_host ; call ld%update_host ; call wd%update_host
      elseif (wd%on_device .or. cd%on_device) then
         m = m + 1_I4P
      endif
      m = m + abs(st - st2) + cut_cells_3d_mismatches(ch, cd) + cut_links_mismatches(lh, ld) + cut_walls_mismatches(wh, wd)
   enddo
   call ld%exit_device ; call cd%exit_device
   ! alpha alone: the terms are reset (released first)
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, alpha_only=.true., walls=wd)
   if (wd%n /= 0_I4P .or. allocated(wd%cell_start) .or. allocated(wd%area) .or. wd%on_device) m = m + 1_I4P
   ! the host path over resident terms (released first)
   call cut_cells_3d_device(s, x, y, z, cd, status=st2, resident=.true., walls=wd)
   call cd%exit_device
   call cut_cells_3d(s, x, y, z, cd, status=st2, walls=wd)
   m = m + abs(st - st2) + cut_walls_mismatches(wh, wd)
   if (wd%on_device) m = m + 1_I4P
   nm = count_multi(wh, ch) ! ch: the cut with links, the same cut
   if (present(multi)) then
      if (multi .and. nm == 0_I4P) m = m + 1_I4P
   endif
   nmulti = nmulti + nm
   total = total + m ; ok = ok .and. m == 0_I4P .and. st == 0_I4P
   print '(A,I0,A,I0,A,I0,A,I0,A,I0)', label//': cells ', size(ch%alpha), ', wall terms ', wh%n, ' (', nm, &
         ' in multi-piece cells, ', count(wh%piece == 0), ' without piece), device path mismatches ', m
   endsubroutine body

   function count_multi(w, c) result(n)
   !< Terms of w in the multi-piece cells of the cut c.
   type(cut_walls_3d_t), intent(in) :: w !< Wall terms.
   type(cut_cells_3d_t), intent(in) :: c !< Cut.
   integer(I4P)                     :: n !< Terms.
   integer(I4P)                     :: i, j, k, l !< Work.

   n = 0_I4P
   do k = 1, c%nz
      do j = 1, c%ny
         do i = 1, c%nx
            if (c%pieces(i, j, k) < 2_I4P) cycle
            l = i + c%nx * ((j - 1_I4P) + c%ny * (k - 1_I4P))
            n = n + w%cell_start(l + 1) - w%cell_start(l)
         enddo
      enddo
   enddo
   endfunction count_multi

   subroutine box_grid(m)
   !< m^3 uniform cells over the box of s grown by 10 %.
   integer(I4P), intent(in) :: m      !< Cells per axis.
   type(vector_R8P)         :: b0, b1 !< Box.
   integer(I4P)             :: q      !< Counter.

   b0 = s%get_bmin() ; b1 = s%get_bmax()
   x = [(b0%x - 0.1_R8P * (b1%x - b0%x) + 1.2_R8P * (b1%x - b0%x) * q / m, q = 0, m)]
   y = [(b0%y - 0.1_R8P * (b1%y - b0%y) + 1.2_R8P * (b1%y - b0%y) * q / m, q = 0, m)]
   z = [(b0%z - 0.1_R8P * (b1%z - b0%z) + 1.2_R8P * (b1%z - b0%z) * q / m, q = 0, m)]
   endsubroutine box_grid
endprogram fossil_test_cut_walls_3d_device_twin
