!< FOSSIL test: the 3-D cut cells of several posed bodies (#109 WP3 step 03).
program fossil_test_cut_bodies
!< FOSSIL test: the 3-D cut cells of several posed bodies (#109 WP3 step 03).
!<
!< `cut_cells_3d_bodies` / `merge_posed_bodies` (`fossil_cut_bodies`):
!< 1. two posed marching-cubes spheres: the cut (with links, walls, edges; then moments and quadrature) equals the cut of
!<    the copies moved by `move` and merged by `merge_solids` bit for bit; the merged facets equal those of that copy;
!<    the facet map and the wall labels agree; the solid volume of the cut cells is V1 + V2; the device path (built
!<    without OpenACC: the device loops on the host) equals the host path bit for bit, copied back and resident, and
!<    gives the same labels; moments with `device` are refused;
!< 2. the two spheres in one cell: the cell's piece is bounded by both walls and lists both bodies;
!< 3. overlaps refused (STATUS_INVALID_INPUT, the cut reset): crossing spheres, a sphere nested in a big cube (both
!<    orders), coincident spheres, two cubes overlapping on coplanar faces of the same orientation, two cubes crossing;
!< 4. touching allowed: two cubes along a face (aligned, shifted), an edge (aligned, shifted), a corner, on grids whose
!<    planes miss and contain the contact: every cell trusted, the solid volume 2 to round-off;
!< 5. the merge alone: patches carried by name, bad poses refused.
use fossil,                   only : surface_stl_object, cut_cells_3d_t, cut_cells_3d, cut_walls_3d_t, cut_links_t, &
                                     rigid_pose_t, pose_from_axis_angle, piece_bodies_t, merge_posed_bodies, &
                                     cut_cells_3d_bodies, STATUS_OK, STATUS_INVALID_INPUT
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_cut_cells_compare, only : cut_cells_3d_mismatches, cut_walls_mismatches, cut_links_mismatches, &
                                     piece_bodies_mismatches
use fossil_facet_object,      only : facet_object
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

type(surface_stl_object) :: sph, cube, big
type(surface_stl_object) :: spheres(2), cubes(2), big_sph(2), sph_big(2) !< Body pairs (no array constructors of a
                                                                         !< finalizable type).
type(rigid_pose_t)       :: p2(2)
real(R8P), allocatable   :: x(:), y(:), z(:)
integer(I4P)             :: i
logical                  :: ok

ok = .true.
call mc_body(BODY_SPHERE, 10_I4P, sph)
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
big = cube
call big%resize(factor=vector_R8P(3._R8P, 3._R8P, 3._R8P))
call big%analyze
spheres(1) = sph ; spheres(2) = sph ; cubes(1) = cube ; cubes(2) = cube
big_sph(1) = big ; big_sph(2) = sph ; sph_big(1) = sph ; sph_big(2) = big
p2(1) = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [-0.62_R8P, 0.05_R8P, 0.1_R8P])
p2(2) = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], -0.7_R8P, [0.61_R8P, -0.03_R8P, 0.02_R8P])
x = [(-1.31_R8P + 0.2_R8P * i, i = 0, 13)] ; y = [(-0.73_R8P + 0.2_R8P * i, i = 0, 7)] ; z = y
call twin_case
call one_cell_case
call overlap_case
call touch_case
call merge_case
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine twin_case
   !< Case 1: equal to the moved + merged copy, labels, volume, device twin.
   type(surface_stl_object) :: a, b, mg
   type(cut_cells_3d_t)     :: ch, cb, cd
   type(cut_walls_3d_t)     :: wh, wb, wd
   type(cut_links_t)        :: lh, lb, ld
   type(piece_bodies_t)     :: pb, pd
   integer(I4P), allocatable :: fb(:), fs(:)
   integer(I4P)             :: m, st, nc, k, f, r, n1
   real(R8P)                :: vs, v12

   m = 0_I4P
   a = sph ; call a%move(p2(1)) ; b = sph ; call b%move(p2(2)) ; call a%merge_solids(b)
   call cut_cells_3d(a, x, y, z, ch, status=st, links=lh, edges=.true., walls=wh)
   call cut_cells_3d_bodies(spheres, p2, x, y, z, cb, status=st, links=lb, walls=wb, edges=.true., piece_bodies=pb, &
                            merged=mg, facet_body=fb, facet_sub_body=fs, contacts=nc)
   if (st /= STATUS_OK .or. nc /= 0_I4P) m = m + 1_I4P
   m = m + cut_cells_3d_mismatches(ch, cb) + cut_links_mismatches(lh, lb) + cut_walls_mismatches(wh, wb)
   m = m + facet_mismatches(a, mg)
   n1 = sph%get_facets_number()
   if (size(fb) /= 2 * n1 .or. size(fs) /= 2 * n1) then
      m = m + 1_I4P
   else
      if (any(fb(1:n1) /= 1_I4P) .or. any(fb(n1 + 1:) /= 2_I4P) .or. any(fs /= 1_I4P)) m = m + 1_I4P
   endif
   ! labels: each wall term carries its facet's body; each piece with a term lists bodies, sorted, distinct
   if (pb%np /= cb%np .or. pb%nw /= wb%n) m = m + 1_I4P
   do k = 1, wb%n
      f = wb%facet(k)
      if (pb%wall_body(k) /= fb(f) .or. pb%wall_sub_body(k) /= fs(f)) m = m + 1_I4P
   enddo
   do k = 1, pb%np
      if (pb%start(k + 1) < pb%start(k)) m = m + 1_I4P
      do f = pb%start(k) + 1, pb%start(k + 1) - 1
         if (pb%body(f) <= pb%body(f - 1)) m = m + 1_I4P
      enddo
   enddo
   if (count(pb%start(2:) > pb%start(:pb%np)) /= n_with_walls(wb, pb%np)) m = m + 1_I4P
   ! volume: the solid of the cut cells is V1 + V2
   vs = solid_volume(cb, x, y, z) ; v12 = 2._R8P * sph%get_volume()
   if (abs(vs - v12) > 1.e-13_R8P) m = m + 1_I4P
   print '(A,I0,A,I0,A,I0,A,ES9.2,A,I0)', 'two spheres: cells ', size(cb%alpha), ', pieces ', cb%np, ', wall terms ', wb%n, &
         ', |V_solid - V1 - V2| ', abs(vs - v12), ', mismatches ', m
   ok = ok .and. m == 0_I4P
   ! moments and quadrature (host only)
   m = 0_I4P
   call cut_cells_3d(a, x, y, z, ch, status=st, moments_order=2_I4P, quadrature=.true.)
   call cut_cells_3d_bodies(spheres, p2, x, y, z, cb, status=st, moments_order=2_I4P, quadrature=.true.)
   if (st /= STATUS_OK) m = m + 1_I4P
   m = m + bits_mismatches(reshape(ch%piece_moments, [size(ch%piece_moments)]), &
                           reshape(cb%piece_moments, [size(cb%piece_moments)]))
   m = m + bits_mismatches(ch%quad_weight, cb%quad_weight) + &
       bits_mismatches(reshape(ch%alpha, [size(ch%alpha)]), reshape(cb%alpha, [size(cb%alpha)]))
   call cut_cells_3d_bodies(spheres, p2, x, y, z, cd, status=st, moments_order=2_I4P, device=.true.)
   if (st /= STATUS_INVALID_INPUT .or. cd%nx /= 0_I4P) m = m + 1_I4P
   print '(A,I0)', 'two spheres, moments and quadrature: mismatches ', m
   ok = ok .and. m == 0_I4P
   ! device twin: copied back and resident
   m = 0_I4P
   call cut_cells_3d_bodies(spheres, p2, x, y, z, cb, status=st, links=lb, walls=wb, piece_bodies=pb)
   do r = 1, 2
      call cut_cells_3d_bodies(spheres, p2, x, y, z, cd, status=st, links=ld, walls=wd, piece_bodies=pd, device=.true., &
                               resident=r == 2)
      if (r == 2) then
         call cd%update_host ; call ld%update_host ; call wd%update_host
      endif
      m = m + abs(st) + cut_cells_3d_mismatches(cb, cd) + cut_links_mismatches(lb, ld) + cut_walls_mismatches(wb, wd) + &
          piece_bodies_mismatches(pb, pd)
   enddo
   call wd%exit_device ; call ld%exit_device ; call cd%exit_device
   call cut_cells_3d_bodies(spheres, p2, x, y, z, cd, status=st, piece_bodies=pd, device=.true., resident=.true.)
   call cd%update_host
   m = m + abs(st) + piece_bodies_mismatches(pb, pd)
   call cd%exit_device
   print '(A,I0)', 'two spheres, device path: mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine twin_case

   subroutine one_cell_case
   !< Case 2: both spheres in one cell, one piece bounded by both walls.
   type(cut_cells_3d_t) :: cb
   type(piece_bodies_t) :: pb
   real(R8P)            :: gx(3), gy(3), gz(3)
   integer(I4P)         :: m, st, p

   m = 0_I4P
   gx = [-1.4_R8P, 1.4_R8P, 1.7_R8P] ; gy = [-0.8_R8P, 0.8_R8P, 1.1_R8P] ; gz = gy
   call cut_cells_3d_bodies(spheres, p2, gx, gy, gz, cb, status=st, piece_bodies=pb)
   if (st /= STATUS_OK .or. cb%pieces(1, 1, 1) /= 1_I4P) then
      m = m + 1_I4P
   else
      p = cb%first_piece(1, 1, 1)
      if (pb%start(p + 1) - pb%start(p) /= 2_I4P) then
         m = m + 1_I4P
      else
         if (any(pb%body(pb%start(p):pb%start(p) + 1) /= [1_I4P, 2_I4P])) m = m + 1_I4P
      endif
   endif
   if (abs(solid_volume(cb, gx, gy, gz) - 2._R8P * sph%get_volume()) > 1.e-13_R8P) m = m + 1_I4P
   print '(A,I0,A,I0)', 'two spheres in one cell: bodies of its piece ', pb%start(2) - pb%start(1), ', mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine one_cell_case

   subroutine overlap_case
   !< Case 3: overlaps refused.
   type(rigid_pose_t) :: q(2)
   integer(I4P)       :: m

   m = 0_I4P
   q(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [-0.3_R8P, 0._R8P, 0._R8P])
   q(2) = pose_from_axis_angle([0._R8P, 0._R8P, 1._R8P], 0.3_R8P, [0.3_R8P, 0.01_R8P, 0._R8P])
   m = m + refused('crossing spheres', spheres, q)
   q(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [-1.5_R8P, -1.5_R8P, -1.5_R8P])
   q(2) = pose_from_axis_angle([1._R8P, 1._R8P, 0._R8P], 0.2_R8P, [0.1_R8P, 0._R8P, 0._R8P])
   m = m + refused('sphere in a cube', big_sph, q)
   m = m + refused('sphere in a cube (reversed)', sph_big, q(2:1:-1))
   q(1) = p2(1) ; q(2) = p2(1)
   m = m + refused('coincident spheres', spheres, q)
   q(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P)
   q(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0._R8P, 0.5_R8P, 0._R8P])
   m = m + refused('cubes on coplanar faces', cubes, q)
   q(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [0.5_R8P, 0.3_R8P, 0.2_R8P])
   m = m + refused('crossing cubes', cubes, q)
   ok = ok .and. m == 0_I4P
   endsubroutine overlap_case

   function refused(label, s, q) result(m)
   !< 1 unless the cut of s posed by q is refused (and the cut reset).
   character(*),             intent(in) :: label !< Case.
   type(surface_stl_object), intent(in) :: s(:)  !< Bodies.
   type(rigid_pose_t),       intent(in) :: q(:)  !< Poses.
   integer(I4P)                         :: m     !< Mismatches.
   type(cut_cells_3d_t)                 :: cb    !< Cut.
   integer(I4P)                         :: st, nc !< Status, contacts.

   call cut_cells_3d_bodies(s, q, x, y, z, cb, status=st, contacts=nc)
   m = 0_I4P ; if (st /= STATUS_INVALID_INPUT .or. cb%nx /= 0_I4P .or. allocated(cb%alpha)) m = 1_I4P
   print '(A,A,I0,A,I0)', label, ': refused (contacts ', nc, '), mismatches ', m
   endfunction refused

   subroutine touch_case
   !< Case 4: touching cubes allowed, every cell trusted, volume 2.
   real(R8P), parameter :: T(3, 5) = reshape([1._R8P, 0._R8P, 0._R8P, 1._R8P, 0.37_R8P, 0.21_R8P, 1._R8P, 1._R8P, 0._R8P, &
                                              1._R8P, 1._R8P, 0.4_R8P, 1._R8P, 1._R8P, 1._R8P], [3, 5])
   character(12), parameter :: NAMES(5) = [character(12) :: 'face', 'face shifted', 'edge', 'edge shifted', 'corner']
   real(R8P),    parameter :: H(2) = [0.13_R8P, 0.125_R8P]
   type(rigid_pose_t)      :: q(2)
   type(cut_cells_3d_t)    :: cb
   type(piece_bodies_t)    :: pb
   real(R8P), allocatable  :: gx(:)
   real(R8P)               :: vs
   integer(I4P)            :: c, g, m, st, nc, nun, nboth, p

   do c = 1, 5
      do g = 1, 2
         m = 0_I4P
         q(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P)
         q(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, T(:, c))
         gx = [(-0.27_R8P + H(g) * i, i = 0, nint(2.6_R8P / H(g)))]
         call cut_cells_3d_bodies(cubes, q, gx, gx, gx, cb, status=st, piece_bodies=pb, contacts=nc)
         if (st /= STATUS_OK) then
            m = m + 1_I4P ; nun = -1_I4P ; vs = 0._R8P ; nboth = 0_I4P
         else
            nun = count(cb%pieces < 0_I4P .or. .not. cb%closed)
            vs = solid_volume(cb, gx, gx, gx)
            if (nun /= 0_I4P .or. abs(vs - 2._R8P) > 1.e-13_R8P) m = m + 1_I4P
            nboth = 0_I4P
            do p = 1, pb%np
               if (pb%start(p + 1) - pb%start(p) == 2_I4P) nboth = nboth + 1_I4P
            enddo
            if (nboth == 0_I4P) m = m + 1_I4P ! some piece is bounded by both cubes near the contact
         endif
         print '(A,A,A,F6.3,A,I0,A,I0,A,I0,A,ES9.2,A,I0)', 'touching cubes, ', NAMES(c), ', h ', H(g), ': contacts ', nc, &
               ', untrusted ', nun, ', pieces of both ', nboth, ', |V_solid - 2| ', abs(vs - 2._R8P), ', mismatches ', m
         ok = ok .and. m == 0_I4P
      enddo
   enddo
   endsubroutine touch_case

   subroutine merge_case
   !< Case 5: patches carried, bad poses refused.
   type(surface_stl_object)  :: mg, two, pair(2)
   integer(I4P), allocatable :: fb(:), fs(:)
   type(rigid_pose_t)        :: q(2)
   integer(I4P)              :: m, st, f, l
   logical                   :: ov

   m = 0_I4P
   two = cube
   call two%set_patch([1_I4P, 2_I4P, 3_I4P], 'lid')
   pair(1) = cube ; pair(2) = two
   q(1) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P)
   q(2) = pose_from_axis_angle([1._R8P, 0._R8P, 0._R8P], 0._R8P, [3._R8P, 0._R8P, 0._R8P])
   call merge_posed_bodies(pair, q, mg, fb, fs, status=st, overlap=ov)
   if (st /= STATUS_OK .or. ov .or. mg%get_patches_number() /= 2_I4P) m = m + 1_I4P
   do f = 1, mg%get_facets_number()
      l = f - merge(0_I4P, cube%get_facets_number(), fb(f) == 1_I4P)
      if (fb(f) == 1_I4P) then
         if (mg%patch_name(mg%facet_patch(f)) /= cube%patch_name(cube%facet_patch(l))) m = m + 1_I4P
      else
         if (mg%patch_name(mg%facet_patch(f)) /= two%patch_name(two%facet_patch(l))) m = m + 1_I4P
      endif
   enddo
   call merge_posed_bodies(pair, q(1:1), mg, fb, fs, status=st)
   if (st /= STATUS_INVALID_INPUT .or. mg%get_facets_number() /= 0_I4P .or. size(fb) /= 0) m = m + 1_I4P
   print '(A,I0)', 'merge: patches by name, bad poses refused: mismatches ', m
   ok = ok .and. m == 0_I4P
   endsubroutine merge_case

   function facet_mismatches(a, b) result(m)
   !< Facets of a and b whose vertices differ in a bit (or a count mismatch).
   type(surface_stl_object), intent(in) :: a, b !< Surfaces.
   integer(I4P)                         :: m    !< Mismatches.
   type(facet_object), pointer          :: fa, fb !< Facets.
   integer(I4P)                         :: f, k   !< Work.

   m = 0_I4P
   if (a%get_facets_number() /= b%get_facets_number()) then
      m = 1_I4P ; return
   endif
   do f = 1, a%get_facets_number()
      fa => a%facet_at(f) ; fb => b%facet_at(f)
      do k = 1, 3
         if (any(transfer([fa%vertex(k)%x, fa%vertex(k)%y, fa%vertex(k)%z], 0_I8P, 3) /= &
                 transfer([fb%vertex(k)%x, fb%vertex(k)%y, fb%vertex(k)%z], 0_I8P, 3))) m = m + 1_I4P
      enddo
   enddo
   endfunction facet_mismatches

   function bits_mismatches(a, b) result(m)
   !< Elements of a and b differing in a bit (or a size mismatch).
   real(R8P), intent(in) :: a(:), b(:) !< Arrays.
   integer(I4P)          :: m          !< Mismatches.

   if (size(a) /= size(b)) then
      m = max(1_I4P, size(a, kind=I4P), size(b, kind=I4P)) ; return
   endif
   m = count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b)), kind=I4P)
   endfunction bits_mismatches

   function n_with_walls(w, np) result(n)
   !< Pieces with at least one wall term.
   type(cut_walls_3d_t), intent(in) :: w  !< Wall terms.
   integer(I4P),         intent(in) :: np !< Pieces.
   integer(I4P)                     :: n  !< Pieces with terms.
   logical                          :: has(np) !< Piece has a term.
   integer(I4P)                     :: k  !< Work.

   has = .false.
   do k = 1, w%n
      if (w%piece(k) >= 1_I4P .and. w%piece(k) <= np) has(w%piece(k)) = .true.
   enddo
   n = count(has, kind=I4P)
   endfunction n_with_walls

   function solid_volume(c, gx, gy, gz) result(v)
   !< Solid volume of the cut cells: sum of (1 - alpha) |cell| over the cells with solid.
   type(cut_cells_3d_t), intent(in) :: c              !< Cut.
   real(R8P),            intent(in) :: gx(:), gy(:), gz(:) !< Cell edges.
   real(R8P)                        :: v              !< Solid volume.
   integer(I4P)                     :: i, j, k        !< Work.

   v = 0._R8P
   do k = 1, c%nz
      do j = 1, c%ny
         do i = 1, c%nx
            if (c%alpha(i, j, k) < 1._R8P) v = v + (1._R8P - c%alpha(i, j, k)) * (gx(i + 1) - gx(i)) * (gy(j + 1) - gy(j)) * &
                                                (gz(k + 1) - gz(k))
         enddo
      enddo
   enddo
   endfunction solid_volume
endprogram fossil_test_cut_bodies
