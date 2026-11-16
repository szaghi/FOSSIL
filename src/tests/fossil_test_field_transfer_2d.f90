!< FOSSIL test: curve fields at points on segments and back to the vertices (#109 WP2 step 04).
program fossil_test_field_transfer_2d
!< FOSSIL test: curve fields at points on segments and back to the vertices (#109 WP2 step 04).
!<
!< The 2-D twin of `fossil_test_field_transfer_3d`: a 64-gon circle (unposed and at a pose), a disc with a hole, an open
!< polyline. A linear vertex field F(x) = (g . x + c0, h . x - 1) reproduced on the segments and off them (projected), at
!< the vertices (bits), the segment field copied, bad ids zero; the feet and ILW crossings of a posed closure table, the
!< load samples of orders 1-3; samples_to_vertices of a constant (bits) and of the linear field (orders 2-3, the round
!< trip); every device twin bit for bit (resident table and samples too); refusals.
use fossil, only : curve_object, load_samples_t, load_samples_2d, load_samples_reset, ib_closure_t, fd_closure_2d, &
                   ib_closure_reset, transfer_field_2d, transfer_to_closure, transfer_to_samples, samples_to_vertices, &
                   rigid_pose_2d_t, pose_2d, to_world, to_body, STATUS_OK, STATUS_INVALID_INPUT
use penf,   only : I4P, I8P, R8P

implicit none

real(R8P), parameter  :: PI = acos(-1._R8P), EPS = epsilon(1._R8P), TOL = 4._R8P
real(R8P), parameter  :: G(2) = [0.4_R8P, -0.9_R8P], H(2) = [-0.3_R8P, 0.7_R8P], C0 = 1.5_R8P
type(curve_object)    :: circle, holed, open
type(rigid_pose_2d_t) :: pose
real(R8P), allocatable :: u(:), v(:)
real(R8P)             :: worst, worst_off, worst_cl, worst_s, worst_v, err1
integer(I4P)          :: i
logical               :: ok

ok = .true. ; worst = 0._R8P ; worst_off = 0._R8P ; worst_cl = 0._R8P ; worst_s = 0._R8P ; worst_v = 0._R8P
err1 = 0._R8P
u = [(0.5_R8P + 0.4_R8P * cos(2 * PI * i / 64), i = 0, 63)] ; v = [(0.5_R8P + 0.4_R8P * sin(2 * PI * i / 64), i = 0, 63)]
call circle%adopt_loops(u=u, v=v, loop_sizes=[64_I4P])
u = [u, [(0.5_R8P + 0.15_R8P * cos(-2 * PI * i / 32), i = 0, 31)]]
v = [v, [(0.6_R8P + 0.15_R8P * sin(-2 * PI * i / 32), i = 0, 31)]]
call holed%adopt_loops(u=u, v=v, loop_sizes=[64_I4P, 32_I4P])
u = [(0.1_R8P * i, i = 0, 10)] ; v = [(0.2_R8P + 0.05_R8P * i + 0.02_R8P * (i * i), i = 0, 10)]
call open%adopt_loops(u=u, v=v, loop_sizes=[11_I4P], closed=[.false.])
pose = pose_2d(0.35_R8P, [0.4_R8P, -0.7_R8P])
call circle%enter_device ; call holed%enter_device ; call open%enter_device ! the resident paths need them mapped
call on_segments(circle, 'circle', .false.) ; call on_segments(circle, 'circle posed', .true.)
call on_segments(holed, 'disc with a hole', .false.) ; call on_segments(open, 'open polyline', .false.)
print '(A,F6.2,A,F6.2)', 'linear field on the segments, largest error / (eps S): ', worst, ';  off the segments ', worst_off
call closure_case(holed)
call samples_case(circle, 'circle') ; call samples_case(holed, 'disc with a hole') ; call samples_case(open, 'open polyline')
print '(A,F6.2,A,F6.2)', 'closure feet and ILW, largest error / (eps S): ', worst_cl, ';  load samples ', worst_s
call back_case(circle, 'circle') ; call back_case(holed, 'disc with a hole') ; call back_case(open, 'open polyline')
print '(A,F7.2,A,ES9.2)', 'samples_to_vertices, linear field, largest error / (eps S): ', worst_v, &
      ';  order 1 (segment means), largest error / S: ', err1
call refusals
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) stop 1

contains
   function field(x) result(f)
   !< The linear field F(x) [2].
   real(R8P), intent(in) :: x(2) !< Point.
   real(R8P)             :: f(2) !< Value.

   f(1) = (G(1) * x(1) + G(2) * x(2)) + C0
   f(2) = (H(1) * x(1) + H(2) * x(2)) - 1._R8P
   endfunction field

   real(R8P) function scale(x) result(s)
   !< The sum of the absolute terms of F at x.
   real(R8P), intent(in) :: x(2) !< Point.

   s = sum(abs(G) * abs(x)) + sum(abs(H) * abs(x)) + C0 + 1._R8P
   endfunction scale

   subroutine vertex_fields(c, vf, sf)
   !< The linear field at the vertices [2, nv] and a segment field [2, nv].
   type(curve_object),     intent(in)  :: c       !< Curve.
   real(R8P), allocatable, intent(out) :: vf(:,:) !< Vertex field.
   real(R8P), allocatable, intent(out) :: sf(:,:) !< Segment field.
   integer(I4P) :: k, nv                          !< Counters.

   nv = c%get_vertices_number()
   allocate(vf(2, nv), sf(2, nv))
   do k = 1, nv
      vf(:, k) = field(c%get_vertex(k)) ; sf(1, k) = real(k, R8P) ; sf(2, k) = -2._R8P * k + 0.5_R8P
   enddo
   endsubroutine vertex_fields

   subroutine on_segments(c, name, posed)
   !< Points on every segment (four parameters), off it along its normal, at its ends, and bad ids.
   type(curve_object), intent(in) :: c     !< Curve.
   character(*),       intent(in) :: name  !< Label.
   logical,            intent(in) :: posed !< At the pose.
   real(R8P), parameter :: SP(4) = [0.5_R8P, 0.13_R8P, 0.77_R8P, 0.3_R8P]
   real(R8P), allocatable :: vf(:,:), sf(:,:), pb(:,:), pw(:,:), po(:,:), val(:,:), vd(:,:), pv(:,:)
   integer(I4P), allocatable :: fa(:), fv(:)
   real(R8P)    :: a(2), b(2), nrm(2), r
   integer(I4P) :: nv, e, j, k, i, st, sd, np
   logical      :: good, bits, vgood

   call vertex_fields(c, vf, sf)
   nv = c%get_vertices_number()
   allocate(pb(2, 4 * nv), pw(2, 4 * nv), po(2, 4 * nv), fa(4 * nv), val(2, 4 * nv), vd(2, 4 * nv))
   allocate(pv(2, 2 * nv + 2), fv(2 * nv + 2))
   i = 0 ; np = 0
   do e = 1, nv
      j = c%get_next(e)
      a = c%get_vertex(e) ; b = a ; if (j > 0) b = c%get_vertex(j)
      nrm = [b(2) - a(2), a(1) - b(1)] ; if (j > 0) nrm = nrm / sqrt(sum(nrm * nrm))
      do k = 1, 4
         i = i + 1 ; fa(i) = e
         pb(:, i) = a + SP(k) * (b - a)
         pw(:, i) = pb(:, i) ; if (posed) pw(:, i) = to_world(pose, pb(:, i))
         po(:, i) = pw(:, i) + 0.02_R8P * merge(to_world(pose, nrm) - pose%t, nrm, posed)
      enddo
      pv(:, 2 * e - 1) = a ; pv(:, 2 * e) = b ; fv(2 * e - 1) = e ; fv(2 * e) = e
   enddo
   pv(:, 2 * nv + 1:) = 0._R8P ; fv(2 * nv + 1) = 0 ; fv(2 * nv + 2) = nv + 1
   if (posed) then
      call transfer_field_2d(c, fa, pw, val, status=st, vertex_field=vf, pose=pose)
      call transfer_field_2d(c, fa, pw, vd, status=sd, vertex_field=vf, pose=pose, device=.true.)
   else
      call transfer_field_2d(c, fa, pw, val, status=st, vertex_field=vf)
      call transfer_field_2d(c, fa, pw, vd, status=sd, vertex_field=vf, device=.true.)
   endif
   good = st == STATUS_OK .and. sd == STATUS_OK ; bits = same2(val, vd)
   do i = 1, 4 * nv
      r = maxval(abs(val(:, i) - field(pb(:, i)))) / (EPS * scale(pb(:, i)))
      worst = max(worst, r) ; good = good .and. r <= TOL
   enddo
   call check(name//': linear vertex field on the segments (TOL eps S)', good)
   if (posed) then
      call transfer_field_2d(c, fa, po, val, status=st, vertex_field=vf, pose=pose)
      call transfer_field_2d(c, fa, po, vd, status=sd, vertex_field=vf, pose=pose, device=.true.)
   else
      call transfer_field_2d(c, fa, po, val, status=st, vertex_field=vf)
      call transfer_field_2d(c, fa, po, vd, status=sd, vertex_field=vf, device=.true.)
   endif
   good = st == STATUS_OK .and. sd == STATUS_OK ; bits = bits .and. same2(val, vd)
   do i = 1, 4 * nv
      r = maxval(abs(val(:, i) - field(pb(:, i)))) / (EPS * scale(pb(:, i)))
      worst_off = max(worst_off, r) ; good = good .and. r <= 4 * TOL
   enddo
   call check(name//': points off the segments take the value of their projection', good)
   call transfer_field_2d(c, fa, pw, val, status=st, segment_field=sf)
   call transfer_field_2d(c, fa, pw, vd, status=sd, segment_field=sf, device=.true.)
   bits = bits .and. same2(val, vd)
   call check(name//': segment field copied (bits)', st == STATUS_OK .and. same2(val, sf(:, fa)))
   if (.not. posed) then
      deallocate(val, vd) ; allocate(val(2, size(fv)), vd(2, size(fv)))
      call transfer_field_2d(c, fv, pv, val, status=st, vertex_field=vf)
      call transfer_field_2d(c, fv, pv, vd, status=sd, vertex_field=vf, device=.true.)
      bits = bits .and. same2(val, vd)
      vgood = st == STATUS_OK
      do e = 1, nv
         j = c%get_next(e) ; if (j < 1) j = e
         vgood = vgood .and. same(val(:, 2 * e - 1), vf(:, e)) .and. same(val(:, 2 * e), vf(:, j))
      enddo
      call check(name//': a point at a segment end takes the vertex value (bits)', vgood)
      call check(name//': segment ids 0 and nv + 1 give zeros', all(val(:, 2 * nv + 1:) == 0._R8P))
   endif
   call check(name//': device twin bit for bit', bits)
   endsubroutine on_segments

   subroutine closure_case(c)
   !< The feet and ILW crossings of a posed closure table with ILW (host, device, resident).
   type(curve_object), intent(inout) :: c !< Curve.
   type(ib_closure_t)     :: cl, cr
   real(R8P), allocatable :: vf(:,:), sf(:,:), x(:), y(:), val(:,:), vd(:,:), iv(:,:,:), id(:,:,:), fv(:,:), fi(:,:,:)
   real(R8P)    :: q(2), r
   integer(I4P) :: k, g, sl, st, nx, ncr
   logical      :: good, bits, zgood

   call vertex_fields(c, vf, sf)
   nx = 48
   allocate(x(nx), y(nx))
   do k = 1, nx
      x(k) = 0.1_R8P + 1.2_R8P * (k - 0.5_R8P) / nx ; y(k) = -0.7_R8P + 1.2_R8P * (k - 0.5_R8P) / nx
   enddo
   call fd_closure_2d(c, x, y, [2_I4P, 2_I4P], cl, status=st, pose=pose, ilw=.true.)
   call check('closure: table built with ghosts', st == STATUS_OK .and. cl%n > 0)
   allocate(val(2, cl%n), vd(2, cl%n), iv(2, 4, cl%n), id(2, 4, cl%n), fv(2, cl%n), fi(2, 4, cl%n))
   call transfer_to_closure(c, cl, val, status=st, vertex_field=vf, ilw_values=iv, pose=pose)
   good = st == STATUS_OK ; zgood = .true. ; ncr = 0
   do g = 1, cl%n
      q = to_body(pose, cl%foot(:, g))
      r = maxval(abs(val(:, g) - field(q))) / (EPS * scale(q))
      worst_cl = max(worst_cl, r) ; good = good .and. r <= TOL
      do sl = 1, 4
         if (cl%ilw_element(sl, g) < 1) then
            zgood = zgood .and. all(iv(:, sl, g) == 0._R8P)
         else
            ncr = ncr + 1
            q = cl%points(:, g) ; q((sl + 1) / 2) = cl%ilw_t(sl, g) ; q = to_body(pose, q)
            r = maxval(abs(iv(:, sl, g) - field(q))) / (EPS * scale(q))
            worst_cl = max(worst_cl, r) ; good = good .and. r <= TOL
         endif
      enddo
   enddo
   call check('closure: linear field at the feet and ILW crossings (TOL eps S)', good .and. ncr > 0)
   call check('closure: ILW slots without a crossing give zeros', zgood)
   call transfer_to_closure(c, cl, fv, status=st, segment_field=sf, ilw_values=fi, pose=pose)
   good = st == STATUS_OK .and. same2(fv, sf(:, cl%element))
   do g = 1, cl%n
      do sl = 1, 4
         if (cl%ilw_element(sl, g) > 0) good = good .and. same(fi(:, sl, g), sf(:, cl%ilw_element(sl, g)))
      enddo
   enddo
   call check('closure: segment field copied at the feet and crossings (bits)', good)
   call transfer_to_closure(c, cl, vd, vertex_field=vf, ilw_values=id, pose=pose, device=.true.)
   bits = same2(val, vd) .and. same3(iv, id)
   call fd_closure_2d(c, x, y, [2_I4P, 2_I4P], cr, pose=pose, ilw=.true., resident=.true.)
   call transfer_to_closure(c, cr, vd, vertex_field=vf, ilw_values=id, pose=pose, device=.true.)
   bits = bits .and. same2(val, vd) .and. same3(iv, id)
   call transfer_to_closure(c, cr, vd, segment_field=sf, ilw_values=id, pose=pose, device=.true.)
   bits = bits .and. same2(fv, vd) .and. same3(fi, id)
   call check('closure: device twin bit for bit (host and resident table)', bits)
   call ib_closure_reset(cr)
   endsubroutine closure_case

   subroutine samples_case(c, name)
   !< The load samples of orders 1-3 (host, device, resident samples).
   type(curve_object), intent(in) :: c    !< Curve.
   character(*),       intent(in) :: name !< Label.
   type(load_samples_t)   :: ls, lr
   real(R8P), allocatable :: vf(:,:), sf(:,:), val(:,:), vd(:,:)
   real(R8P)    :: r
   integer(I4P) :: o, k, st
   logical      :: good, bits, fgood

   call vertex_fields(c, vf, sf)
   good = .true. ; bits = .true. ; fgood = .true.
   do o = 1, 3
      call load_samples_2d(c, samples=ls, order=o)
      call load_samples_2d(c, samples=lr, order=o, resident=.true.)
      if (allocated(val)) deallocate(val, vd)
      allocate(val(2, ls%n), vd(2, ls%n))
      call transfer_to_samples(c, ls, val, status=st, vertex_field=vf)
      good = good .and. st == STATUS_OK
      do k = 1, ls%n
         r = maxval(abs(val(:, k) - field(ls%points(:, k)))) / (EPS * scale(ls%points(:, k)))
         worst_s = max(worst_s, r) ; good = good .and. r <= TOL
      enddo
      call transfer_to_samples(c, ls, vd, vertex_field=vf, device=.true.) ; bits = bits .and. same2(val, vd)
      call transfer_to_samples(c, lr, vd, vertex_field=vf, device=.true.) ; bits = bits .and. same2(val, vd)
      call transfer_to_samples(c, ls, val, status=st, segment_field=sf)
      fgood = fgood .and. st == STATUS_OK .and. same2(val, sf(:, ls%element))
      call transfer_to_samples(c, lr, vd, segment_field=sf, device=.true.) ; bits = bits .and. same2(val, vd)
      call load_samples_reset(lr)
   enddo
   call check(name//': linear vertex field at the load samples, orders 1-3 (TOL eps S)', good)
   call check(name//': segment field at the load samples (bits)', fgood)
   call check(name//': load samples device twin bit for bit (host and resident samples)', bits)
   endsubroutine samples_case

   subroutine back_case(c, name)
   !< samples_to_vertices: constants (bits), the linear field from order 2, the round trip; device twin.
   type(curve_object), intent(in) :: c    !< Curve.
   character(*),       intent(in) :: name !< Label.
   type(load_samples_t)   :: ls, lr
   real(R8P), allocatable :: vf(:,:), sf(:,:), val(:,:), fl(:,:), fd(:,:), fr(:,:)
   real(R8P)    :: r
   integer(I4P) :: o, k, st, l
   logical      :: cgood, lgood, tgood, bits

   call vertex_fields(c, vf, sf)
   allocate(fl(2, size(vf, 2)), fd(2, size(vf, 2)), fr(2, size(vf, 2)))
   cgood = .true. ; lgood = .true. ; tgood = .true. ; bits = .true.
   do o = 1, 3
      call load_samples_2d(c, samples=ls, order=o)
      call load_samples_2d(c, samples=lr, order=o, resident=.true.)
      if (allocated(val)) deallocate(val)
      allocate(val(2, ls%n))
      val(1, :) = 3.25_R8P ; val(2, :) = -1.1_R8P
      do l = 1, 2
         call samples_to_vertices(c, ls, val, fl, status=st, linear=(l == 1))
         cgood = cgood .and. st == STATUS_OK .and. all(fl(1, :) == 3.25_R8P) .and. all(fl(2, :) == -1.1_R8P)
      enddo
      do k = 1, ls%n
         val(:, k) = field(ls%points(:, k))
      enddo
      call samples_to_vertices(c, ls, val, fl, status=st)
      call samples_to_vertices(c, ls, val, fd, device=.true.) ; bits = bits .and. same2(fl, fd)
      call samples_to_vertices(c, lr, val, fd, device=.true.) ; bits = bits .and. same2(fl, fd)
      call samples_to_vertices(c, lr, val, fd, linear=.false., device=.true.)
      call samples_to_vertices(c, ls, val, fr, linear=.false.) ; bits = bits .and. same2(fr, fd)
      do k = 1, size(vf, 2)
         r = maxval(abs(fl(:, k) - vf(:, k))) / (EPS * scale(c%get_vertex(k)))
         if (o == 1) then
            err1 = max(err1, maxval(abs(fl(:, k) - vf(:, k))) / scale(c%get_vertex(k)))
         else
            worst_v = max(worst_v, r) ; lgood = lgood .and. r <= 4 * TOL
         endif
      enddo
      if (o > 1) then
         call transfer_to_samples(c, ls, val, vertex_field=vf)
         call samples_to_vertices(c, ls, val, fl)
         do k = 1, size(vf, 2)
            tgood = tgood .and. maxval(abs(fl(:, k) - vf(:, k))) <= 4 * TOL * EPS * scale(c%get_vertex(k))
         enddo
      endif
      call load_samples_reset(lr)
   enddo
   call check(name//': samples_to_vertices of a constant is the constant (bits)', cgood)
   call check(name//': samples_to_vertices of a linear field exact at the vertices, orders 2-3 (4 TOL eps S)', lgood)
   call check(name//': round trip vertices -> samples -> vertices of a linear field (4 TOL eps S)', tgood)
   call check(name//': samples_to_vertices device twin bit for bit', bits)
   endsubroutine back_case

   subroutine refusals
   !< Bad inputs: refused, outputs zero.
   type(load_samples_t)   :: ls
   real(R8P), allocatable :: vf(:,:), sf(:,:), val(:,:), fl(:,:)
   real(R8P)    :: pts(2, 2)
   integer(I4P) :: st, fa(2)
   logical      :: good

   call vertex_fields(circle, vf, sf)
   pts = 0.1_R8P ; fa = [1, 2] ; allocate(val(2, 2))
   call transfer_field_2d(circle, fa, pts, val, status=st)
   good = st == STATUS_INVALID_INPUT .and. all(val == 0._R8P)
   call transfer_field_2d(circle, fa, pts, val, status=st, vertex_field=vf, segment_field=sf)
   good = good .and. st == STATUS_INVALID_INPUT
   call transfer_field_2d(circle, fa, pts, val, status=st, vertex_field=vf(:, 1:3))
   good = good .and. st == STATUS_INVALID_INPUT
   call transfer_field_2d(circle, fa, pts(:, 1:1), val, status=st, vertex_field=vf)
   good = good .and. st == STATUS_INVALID_INPUT
   call check('refused: no field, both, sizes of the field, points', good)
   call load_samples_2d(open, samples=ls)
   deallocate(val) ; allocate(val(2, ls%n), fl(2, size(vf, 2)))
   call transfer_to_samples(circle, ls, val, status=st, vertex_field=vf)
   good = st == STATUS_INVALID_INPUT
   call samples_to_vertices(circle, ls, val, fl, status=st)
   good = good .and. st == STATUS_INVALID_INPUT .and. all(fl == 0._R8P)
   call check('refused: samples of another curve', good)
   endsubroutine refusals

   logical function same(x, y)
   !< Equal bits.
   real(R8P), intent(in) :: x(:), y(:) !< Arrays.

   same = size(x) == size(y)
   if (same) same = all(transfer(x, 0_I8P, size(x)) == transfer(y, 0_I8P, size(y)))
   endfunction same

   logical function same2(x, y)
   !< Equal bits and shapes.
   real(R8P), intent(in) :: x(:,:), y(:,:) !< Arrays.

   same2 = all(shape(x) == shape(y))
   if (same2) same2 = same(reshape(x, [size(x)]), reshape(y, [size(y)]))
   endfunction same2

   logical function same3(x, y)
   !< Equal bits and shapes.
   real(R8P), intent(in) :: x(:,:,:), y(:,:,:) !< Arrays.

   same3 = all(shape(x) == shape(y))
   if (same3) same3 = same(reshape(x, [size(x)]), reshape(y, [size(y)]))
   endfunction same3

   subroutine check(label, cond)
   !< Print and accumulate one check.
   character(*), intent(in) :: label !< Label.
   logical,      intent(in) :: cond  !< Passed.

   print '(A,L1,2A)', '  [', cond, '] ', label
   ok = ok .and. cond
   endsubroutine check
endprogram fossil_test_field_transfer_2d
