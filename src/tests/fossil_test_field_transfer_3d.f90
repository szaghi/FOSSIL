!< FOSSIL test: wall fields at points on facets and back to the vertices (#109 WP2 step 04).
program fossil_test_field_transfer_3d
!< FOSSIL test: wall fields at points on facets and back to the vertices (#109 WP2 step 04).
!<
!< 1. A linear vertex field F(x) = (g . x + c0, h . x - 1) is reproduced at any point of a facet (barycentric points on
!<    every facet of an icosphere, a tetrahedron, cube-01; unposed and at a pose) and at points off the facet that project
!<    into it, to TOL eps times the sum of the absolute terms (largest ratio printed).
!< 2. A facet field is copied (bits); a point equal to a facet vertex takes that vertex's value (bits); ids outside
!<    1..nf give zeros.
!< 3. `transfer_to_closure` at the feet and the ILW crossings of a posed closure table, `transfer_to_samples` at the load
!<    samples (orders 1-3): the linear field to round-off, the facet field by bits, slots without a crossing zero.
!< 4. `samples_to_vertices`: a constant is the constant (bits; linear fit and mean); a linear field is exact at every
!<    vertex to round-off from order 2 (the round trip vertices -> samples -> vertices too); order 1 printed.
!< 5. Twin: `device=.true.` (built without OpenACC its loops run on the host) equals the host path bit for bit, a resident
!<    closure table and resident samples too. 6. Refusals.
use fossil,                  only : surface_stl_object, load_samples_t, load_samples_3d, load_samples_reset, ib_closure_t, &
                                    fd_closure_3d, ib_closure_reset, transfer_field_3d, transfer_to_closure, &
                                    transfer_to_samples, samples_to_vertices, rigid_pose_t, pose_from_axis_angle, &
                                    to_world, to_body, STATUS_OK, STATUS_INVALID_INPUT, curve_object, fd_closure_2d
use fossil_facet_object,     only : facet_object
use fossil_thickness_bodies, only : make_sphere
use fossil_vertex_pool_object, only : vertex_pool_object
use penf,                    only : I4P, I8P, R8P
use vecfor,                  only : vector_R8P

implicit none

real(R8P), parameter            :: EPS = epsilon(1._R8P), TOL = 4._R8P
real(R8P), parameter            :: G(3) = [0.3_R8P, -0.7_R8P, 1.1_R8P], H(3) = [-0.5_R8P, 0.2_R8P, 0.9_R8P], C0 = 2._R8P
type(surface_stl_object)        :: cube, sphere, tet
type(facet_object), allocatable :: fc(:)
type(rigid_pose_t)              :: pose
real(R8P)                       :: worst, worst_off, worst_cl, worst_s, worst_v, err1
logical                         :: ok

ok = .true. ; worst = 0._R8P ; worst_off = 0._R8P ; worst_cl = 0._R8P ; worst_s = 0._R8P ; worst_v = 0._R8P
err1 = 0._R8P
call cube%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
call make_sphere(3_I4P, fc) ; call sphere%adopt_facets(facets=fc)
call make_tet([0._R8P, 0._R8P, 0._R8P], [1.3_R8P, 0.1_R8P, 0.2_R8P], [0.2_R8P, 0.9_R8P, -0.1_R8P], &
              [0.4_R8P, 0.3_R8P, 1.1_R8P], fc)
call tet%adopt_facets(facets=fc)
pose = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.4_R8P, [0.7_R8P, -0.3_R8P, 1.2_R8P])
call sphere%enter_device ; call tet%enter_device ; call cube%enter_device ! the resident paths need them mapped
! 1, 2: points on and off the facets
call on_facets(sphere, 'icosphere', .false.) ; call on_facets(sphere, 'icosphere posed', .true.)
call on_facets(tet, 'tetrahedron', .false.) ; call on_facets(tet, 'tetrahedron posed', .true.)
call on_facets(cube, 'cube-01', .false.)
print '(A,F6.2,A,F6.2)', 'linear field on the facets, largest error / (eps S): ', worst, ';  off the facets ', worst_off
! 3. closure and samples
call closure_case(sphere)
call samples_case(sphere, 'icosphere') ; call samples_case(tet, 'tetrahedron') ; call samples_case(cube, 'cube-01')
print '(A,F6.2,A,F6.2)', 'closure feet and ILW, largest error / (eps S): ', worst_cl, ';  load samples ', worst_s
! 4. samples to vertices
call back_case(sphere, 'icosphere') ; call back_case(tet, 'tetrahedron') ; call back_case(cube, 'cube-01')
print '(A,F7.2,A,ES9.2)', 'samples_to_vertices, linear field, largest error / (eps S): ', worst_v, &
      ';  order 1 (facet means), largest error / S: ', err1
! 6. refusals
call refusals
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) stop 1

contains
   function field(x) result(f)
   !< The linear field F(x) [2].
   real(R8P), intent(in) :: x(3) !< Point.
   real(R8P)             :: f(2) !< Value.

   f(1) = (G(1) * x(1) + G(2) * x(2)) + G(3) * x(3) + C0
   f(2) = (H(1) * x(1) + H(2) * x(2)) + H(3) * x(3) - 1._R8P
   endfunction field

   real(R8P) function scale(x) result(s)
   !< The sum of the absolute terms of F at x (the error scale).
   real(R8P), intent(in) :: x(3) !< Point.

   s = sum(abs(G) * abs(x)) + sum(abs(H) * abs(x)) + C0 + 1._R8P
   endfunction scale

   subroutine vertex_fields(s, vf, ff)
   !< The linear field at the surface's vertices [2, nv] and a facet field [2, nf] (f, -2 f + 0.5).
   type(surface_stl_object), intent(in), target :: s       !< Surface.
   real(R8P), allocatable,   intent(out)        :: vf(:,:) !< Vertex field.
   real(R8P), allocatable,   intent(out)        :: ff(:,:) !< Facet field.
   type(vertex_pool_object), pointer :: pool !< Vertex pool.
   type(vector_R8P) :: c                     !< A vertex.
   integer(I4P)     :: v, f                  !< Counters.

   pool => s%get_vertex_pool()
   allocate(vf(2, pool%vertex_count()), ff(2, s%get_facets_number()))
   do v = 1, pool%vertex_count()
      c = pool%coord(v) ; vf(:, v) = field([c%x, c%y, c%z])
   enddo
   do f = 1, s%get_facets_number()
      ff(1, f) = real(f, R8P) ; ff(2, f) = -2._R8P * f + 0.5_R8P
   enddo
   endsubroutine vertex_fields

   subroutine on_facets(s, name, posed)
   !< Points on every facet (six barycentric combinations), off it along its normal, at its vertices, and bad ids.
   type(surface_stl_object), intent(in), target :: s     !< Surface.
   character(*),             intent(in)         :: name  !< Label.
   logical,                  intent(in)         :: posed !< At the pose.
   real(R8P), parameter :: LB(3, 6) = reshape([1._R8P / 3, 1._R8P / 3, 1._R8P / 3, 0.6_R8P, 0.3_R8P, 0.1_R8P, &
                                               0.15_R8P, 0.7_R8P, 0.15_R8P, 0.2_R8P, 0.2_R8P, 0.6_R8P, &
                                               0.5_R8P, 0.5_R8P, 0._R8P, 0._R8P, 0.25_R8P, 0.75_R8P], [3, 6])
   type(vertex_pool_object), pointer :: pool
   real(R8P), allocatable :: vf(:,:), ff(:,:), pb(:,:), pw(:,:), val(:,:), vd(:,:), po(:,:), pv(:,:)
   integer(I4P), allocatable :: fa(:), fv(:)
   real(R8P)    :: v(3, 3), nrm(3), r, ro
   integer(I4P) :: nf, f, k, i, st, sd
   logical      :: good, bits, vgood, zgood

   call vertex_fields(s, vf, ff)
   pool => s%get_vertex_pool()
   nf = s%get_facets_number()
   allocate(pb(3, 6 * nf), pw(3, 6 * nf), po(3, 6 * nf), fa(6 * nf), val(2, 6 * nf), vd(2, 6 * nf))
   allocate(pv(3, 3 * nf + 2), fv(3 * nf + 2))
   i = 0
   do f = 1, nf
      v = reshape(s%tri_vertices(:, f), [3, 3])
      nrm = cross(v(:, 2) - v(:, 1), v(:, 3) - v(:, 1)) ; nrm = nrm / sqrt(sum(nrm * nrm))
      do k = 1, 6
         i = i + 1 ; fa(i) = f
         pb(:, i) = (LB(1, k) * v(:, 1) + LB(2, k) * v(:, 2)) + LB(3, k) * v(:, 3)
         pw(:, i) = pb(:, i) ; if (posed) pw(:, i) = to_world(pose, pb(:, i))
         po(:, i) = pw(:, i) + 0.03_R8P * merge(to_world(pose, nrm) - pose%t, nrm, posed)
      enddo
      do k = 1, 3
         pv(:, 3 * (f - 1) + k) = v(:, k) ; fv(3 * (f - 1) + k) = f
      enddo
   enddo
   pv(:, 3 * nf + 1:) = 0._R8P ; fv(3 * nf + 1) = 0 ; fv(3 * nf + 2) = nf + 1
   ! on the facets
   if (posed) then
      call transfer_field_3d(s, fa, pw, val, status=st, vertex_field=vf, pose=pose)
      call transfer_field_3d(s, fa, pw, vd, status=sd, vertex_field=vf, pose=pose, device=.true.)
   else
      call transfer_field_3d(s, fa, pw, val, status=st, vertex_field=vf)
      call transfer_field_3d(s, fa, pw, vd, status=sd, vertex_field=vf, device=.true.)
   endif
   good = st == STATUS_OK .and. sd == STATUS_OK ; bits = same2(val, vd)
   do i = 1, 6 * nf
      r = maxval(abs(val(:, i) - field(pb(:, i)))) / (EPS * scale(pb(:, i)))
      worst = max(worst, r) ; good = good .and. r <= TOL
   enddo
   call check(name//': linear vertex field on the facets (TOL eps S)', good)
   ! off the facets (projected): the value of the projection
   if (posed) then
      call transfer_field_3d(s, fa, po, val, status=st, vertex_field=vf, pose=pose)
      call transfer_field_3d(s, fa, po, vd, status=sd, vertex_field=vf, pose=pose, device=.true.)
   else
      call transfer_field_3d(s, fa, po, val, status=st, vertex_field=vf)
      call transfer_field_3d(s, fa, po, vd, status=sd, vertex_field=vf, device=.true.)
   endif
   good = st == STATUS_OK .and. sd == STATUS_OK ; bits = bits .and. same2(val, vd)
   do i = 1, 6 * nf
      if (minval(LB(:, mod(i - 1, 6) + 1)) == 0._R8P) cycle ! on an edge: the projection may leave the facet's plane
      ro = maxval(abs(val(:, i) - field(pb(:, i)))) / (EPS * scale(pb(:, i)))
      worst_off = max(worst_off, ro) ; good = good .and. ro <= 4 * TOL
   enddo
   call check(name//': points off the facets take the value of their projection', good)
   ! facet field: copied
   call transfer_field_3d(s, fa, pw, val, status=st, facet_field=ff)
   call transfer_field_3d(s, fa, pw, vd, status=sd, facet_field=ff, device=.true.)
   bits = bits .and. same2(val, vd)
   call check(name//': facet field copied (bits)', st == STATUS_OK .and. same2(val, ff(:, fa)))
   ! vertices and bad ids (unposed points)
   if (.not. posed) then
      deallocate(val, vd) ; allocate(val(2, size(fv)), vd(2, size(fv)))
      call transfer_field_3d(s, fv, pv, val, status=st, vertex_field=vf)
      call transfer_field_3d(s, fv, pv, vd, status=sd, vertex_field=vf, device=.true.)
      bits = bits .and. same2(val, vd)
      vgood = st == STATUS_OK
      do f = 1, nf
         do k = 1, 3
            vgood = vgood .and. same(val(:, 3 * (f - 1) + k), vf(:, pool%facet_vid(f, k)))
         enddo
      enddo
      zgood = all(val(:, 3 * nf + 1:) == 0._R8P)
      call check(name//': a point at a vertex takes the vertex value (bits)', vgood)
      call check(name//': facet ids 0 and nf + 1 give zeros', zgood)
   endif
   call check(name//': device twin bit for bit', bits)
   endsubroutine on_facets

   subroutine closure_case(s)
   !< The feet and ILW crossings of a posed closure table with ILW (host, device, resident).
   type(surface_stl_object), intent(inout), target :: s !< Surface.
   type(ib_closure_t)     :: cl, cr
   real(R8P), allocatable :: vf(:,:), ff(:,:), x(:), y(:), z(:), val(:,:), vd(:,:), iv(:,:,:), id(:,:,:), fv(:,:), &
                             fi(:,:,:)
   real(R8P)    :: q(3), r
   integer(I4P) :: i, g, sl, st, nx, ncr
   logical      :: good, bits, zgood

   call vertex_fields(s, vf, ff)
   nx = 20
   allocate(x(nx), y(nx), z(nx))
   do i = 1, nx
      x(i) = -0.6_R8P + 2.4_R8P * (i - 0.5_R8P) / nx ; y(i) = -1.6_R8P + 2.4_R8P * (i - 0.5_R8P) / nx
      z(i) = 0.0_R8P + 2.4_R8P * (i - 0.5_R8P) / nx
   enddo
   call fd_closure_3d(s, x, y, z, [2_I4P, 2_I4P, 2_I4P], cl, status=st, pose=pose, ilw=.true.)
   call check('closure: table built with ghosts', st == STATUS_OK .and. cl%n > 0)
   allocate(val(2, cl%n), vd(2, cl%n), iv(2, 6, cl%n), id(2, 6, cl%n), fv(2, cl%n), fi(2, 6, cl%n))
   call transfer_to_closure(s, cl, val, status=st, vertex_field=vf, ilw_values=iv, pose=pose)
   good = st == STATUS_OK ; zgood = .true. ; ncr = 0
   do g = 1, cl%n
      q = to_body(pose, cl%foot(:, g))
      r = maxval(abs(val(:, g) - field(q))) / (EPS * scale(q))
      worst_cl = max(worst_cl, r) ; good = good .and. r <= TOL
      do sl = 1, 6
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
   call transfer_to_closure(s, cl, fv, status=st, facet_field=ff, ilw_values=fi, pose=pose)
   good = st == STATUS_OK .and. same2(fv, ff(:, cl%element))
   do g = 1, cl%n
      do sl = 1, 6
         if (cl%ilw_element(sl, g) > 0) good = good .and. same(fi(:, sl, g), ff(:, cl%ilw_element(sl, g)))
      enddo
   enddo
   call check('closure: facet field copied at the feet and crossings (bits)', good)
   ! twin: device on the host table, device on a resident table
   call transfer_to_closure(s, cl, vd, vertex_field=vf, ilw_values=id, pose=pose, device=.true.)
   bits = same2(val, vd) .and. same3(iv, id)
   call fd_closure_3d(s, x, y, z, [2_I4P, 2_I4P, 2_I4P], cr, pose=pose, ilw=.true., resident=.true.)
   call transfer_to_closure(s, cr, vd, vertex_field=vf, ilw_values=id, pose=pose, device=.true.)
   bits = bits .and. same2(val, vd) .and. same3(iv, id)
   call transfer_to_closure(s, cr, vd, facet_field=ff, ilw_values=id, pose=pose, device=.true.)
   bits = bits .and. same2(fv, vd) .and. same3(fi, id)
   call check('closure: device twin bit for bit (host and resident table)', bits)
   call ib_closure_reset(cr)
   endsubroutine closure_case

   subroutine samples_case(s, name)
   !< The load samples of orders 1-3 (host, device, resident samples).
   type(surface_stl_object), intent(in), target :: s    !< Surface.
   character(*),             intent(in)         :: name !< Label.
   type(load_samples_t)   :: ls, lr
   real(R8P), allocatable :: vf(:,:), ff(:,:), val(:,:), vd(:,:)
   real(R8P)    :: r
   integer(I4P) :: o, i, st
   logical      :: good, bits, fgood

   call vertex_fields(s, vf, ff)
   good = .true. ; bits = .true. ; fgood = .true.
   do o = 1, 3
      call load_samples_3d(s, samples=ls, order=o)
      call load_samples_3d(s, samples=lr, order=o, resident=.true.)
      if (allocated(val)) deallocate(val, vd)
      allocate(val(2, ls%n), vd(2, ls%n))
      call transfer_to_samples(s, ls, val, status=st, vertex_field=vf)
      good = good .and. st == STATUS_OK
      do i = 1, ls%n
         if (ls%weight(i) == 0._R8P) cycle
         r = maxval(abs(val(:, i) - field(ls%points(:, i)))) / (EPS * scale(ls%points(:, i)))
         worst_s = max(worst_s, r) ; good = good .and. r <= TOL
      enddo
      call transfer_to_samples(s, ls, vd, vertex_field=vf, device=.true.) ; bits = bits .and. same2(val, vd)
      call transfer_to_samples(s, lr, vd, vertex_field=vf, device=.true.) ; bits = bits .and. same2(val, vd)
      call transfer_to_samples(s, ls, val, status=st, facet_field=ff)
      fgood = fgood .and. st == STATUS_OK .and. same2(val, ff(:, ls%element))
      call transfer_to_samples(s, lr, vd, facet_field=ff, device=.true.) ; bits = bits .and. same2(val, vd)
      call load_samples_reset(lr)
   enddo
   call check(name//': linear vertex field at the load samples, orders 1-3 (TOL eps S)', good)
   call check(name//': facet field at the load samples (bits)', fgood)
   call check(name//': load samples device twin bit for bit (host and resident samples)', bits)
   endsubroutine samples_case

   subroutine back_case(s, name)
   !< samples_to_vertices: constants (bits), the linear field from order 2, the round trip; device twin.
   type(surface_stl_object), intent(in), target :: s    !< Surface.
   character(*),             intent(in)         :: name !< Label.
   type(vertex_pool_object), pointer :: pool
   type(load_samples_t)   :: ls, lr
   real(R8P), allocatable :: vf(:,:), ff(:,:), val(:,:), fl(:,:), fd(:,:), fr(:,:)
   type(vector_R8P) :: c
   real(R8P)    :: r, xs(3)
   integer(I4P) :: o, i, v, st, l
   logical      :: cgood, lgood, tgood, bits

   call vertex_fields(s, vf, ff)
   pool => s%get_vertex_pool()
   allocate(fl(2, size(vf, 2)), fd(2, size(vf, 2)), fr(2, size(vf, 2)))
   cgood = .true. ; lgood = .true. ; tgood = .true. ; bits = .true.
   do o = 1, 3
      call load_samples_3d(s, samples=ls, order=o)
      call load_samples_3d(s, samples=lr, order=o, resident=.true.)
      if (allocated(val)) deallocate(val)
      allocate(val(2, ls%n))
      ! constant
      val(1, :) = 3.25_R8P ; val(2, :) = -1.1_R8P
      do l = 1, 2
         call samples_to_vertices(s, ls, val, fl, status=st, linear=(l == 1))
         cgood = cgood .and. st == STATUS_OK .and. all(fl(1, :) == 3.25_R8P) .and. all(fl(2, :) == -1.1_R8P)
      enddo
      ! linear
      do i = 1, ls%n
         val(:, i) = field(ls%points(:, i))
      enddo
      call samples_to_vertices(s, ls, val, fl, status=st)
      call samples_to_vertices(s, ls, val, fd, device=.true.) ; bits = bits .and. same2(fl, fd)
      call samples_to_vertices(s, lr, val, fd, device=.true.) ; bits = bits .and. same2(fl, fd)
      call samples_to_vertices(s, lr, val, fd, linear=.false., device=.true.)
      call samples_to_vertices(s, ls, val, fr, linear=.false.) ; bits = bits .and. same2(fr, fd)
      do v = 1, size(vf, 2)
         c = pool%coord(v) ; xs = [c%x, c%y, c%z]
         r = maxval(abs(fl(:, v) - vf(:, v))) / (EPS * scale(xs))
         if (o == 1) then
            err1 = max(err1, maxval(abs(fl(:, v) - vf(:, v))) / scale(xs))
         else
            worst_v = max(worst_v, r) ; lgood = lgood .and. r <= 4 * TOL
         endif
      enddo
      ! round trip: vertex field -> samples -> vertices
      if (o > 1) then
         call transfer_to_samples(s, ls, val, vertex_field=vf)
         call samples_to_vertices(s, ls, val, fl)
         do v = 1, size(vf, 2)
            c = pool%coord(v) ; xs = [c%x, c%y, c%z]
            tgood = tgood .and. maxval(abs(fl(:, v) - vf(:, v))) <= 4 * TOL * EPS * scale(xs)
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
   type(ib_closure_t)     :: c2
   type(load_samples_t)   :: ls
   type(curve_object)     :: cur
   real(R8P), allocatable :: vf(:,:), ff(:,:), val(:,:), fl(:,:)
   real(R8P)    :: pts(3, 2)
   integer(I4P) :: st, fa(2)
   logical      :: good

   call vertex_fields(sphere, vf, ff)
   pts = 0.1_R8P ; fa = [1, 2] ; allocate(val(2, 2))
   call transfer_field_3d(sphere, fa, pts, val, status=st)
   good = st == STATUS_INVALID_INPUT .and. all(val == 0._R8P)
   call transfer_field_3d(sphere, fa, pts, val, status=st, vertex_field=vf, facet_field=ff)
   good = good .and. st == STATUS_INVALID_INPUT
   call transfer_field_3d(sphere, fa, pts, val, status=st, vertex_field=vf(:, 1:3))
   good = good .and. st == STATUS_INVALID_INPUT
   call transfer_field_3d(sphere, fa, pts(:, 1:1), val, status=st, vertex_field=vf)
   good = good .and. st == STATUS_INVALID_INPUT
   call transfer_field_3d(sphere, fa, pts, val(1:1, :), status=st, vertex_field=vf)
   good = good .and. st == STATUS_INVALID_INPUT
   call check('refused: no field, both, sizes of the field, points, values', good)
   call cur%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P], v=[0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P])
   call fd_closure_2d(cur, [(-0.5_R8P + 0.1_R8P * st, st = 0, 20)], [(-0.5_R8P + 0.1_R8P * st, st = 0, 20)], &
                      [2_I4P, 2_I4P], c2)
   call transfer_to_closure(sphere, c2, val, status=st, vertex_field=vf)
   good = st == STATUS_INVALID_INPUT
   call load_samples_3d(tet, samples=ls)
   deallocate(val) ; allocate(val(2, ls%n), fl(2, size(vf, 2)))
   call transfer_to_samples(sphere, ls, val, status=st, vertex_field=vf)
   good = good .and. st == STATUS_INVALID_INPUT
   call samples_to_vertices(sphere, ls, val, fl, status=st)
   good = good .and. st == STATUS_INVALID_INPUT .and. all(fl == 0._R8P)
   call check('refused: a 2-D table, samples of another surface', good)
   endsubroutine refusals

   subroutine make_tet(a, b, cc, d, ff)
   !< The tetrahedron a b cc d, facets outward.
   real(R8P),                       intent(in)  :: a(3), b(3), cc(3), d(3) !< Vertices.
   type(facet_object), allocatable, intent(out) :: ff(:)                   !< Facets.

   allocate(ff(4))
   call put(ff(1), a, cc, b) ; call put(ff(2), a, b, d) ; call put(ff(3), b, cc, d) ; call put(ff(4), cc, a, d)
   endsubroutine make_tet

   subroutine put(fo, a, b, cc)
   !< One facet.
   type(facet_object), intent(inout) :: fo                !< Facet.
   real(R8P),          intent(in)    :: a(3), b(3), cc(3) !< Vertices.

   fo%vertex(1) = vector_R8P(a(1), a(2), a(3)) ; fo%vertex(2) = vector_R8P(b(1), b(2), b(3))
   fo%vertex(3) = vector_R8P(cc(1), cc(2), cc(3)) ; call fo%compute_metrix
   endsubroutine put

   pure function cross(a, b) result(c)
   !< a x b.
   real(R8P), intent(in) :: a(3), b(3) !< Vectors.
   real(R8P)             :: c(3)       !< Cross product.

   c = [a(2) * b(3) - a(3) * b(2), a(3) * b(1) - a(1) * b(3), a(1) * b(2) - a(2) * b(1)]
   endfunction cross

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
endprogram fossil_test_field_transfer_3d
