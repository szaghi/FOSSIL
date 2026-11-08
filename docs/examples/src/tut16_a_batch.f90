!run tut16a tut16_a_batch
!fig tut16-slice
!fig tut16-signs
!fig tut16-slices
!< Tutorial 16, step 1: the distance batch on the device (OpenACC), checked against the host batch and against the
!< same points sent one slice at a time; what the device refuses. Then the figures of the step, from the same arrays:
!< the plane nearest z = 0, the points where the two device signs disagree, and an animation of four slices computed
!< one host call at a time against one device batch. Built with gfortran the !$acc directives are comments: the
!< device kernels run on the host, so the figures show the device code path, not a GPU run.
program tut16_a_batch
!region resident
use fossil, only : surface_stl_object, AABB_TREE_OCTREE, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, &
                   SIGN_RAY_INTERSECTIONS, STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
integer(I4P), parameter       :: n = 24                    ! points per direction
type(surface_stl_object)      :: surface, octree
type(vector_R8P), allocatable :: points(:)
real(R8P),        allocatable :: d_pn(:), d_fwn(:), d_host(:)
integer(I4P)                  :: i, j, k, p, s, status

call surface%load_from_file(file_name='dragon-placed.stl', guess_format=.true.)  ! SAH BVH: the device tree
call surface%enter_device                                   ! copy the flat BVH arrays to the device, once
call surface%enter_device                                   ! a second call is a no-op
print '(A,L1)', 'on the device: ', surface%is_on_device
!endregion resident
!region batch
allocate(points(n**3), d_pn(n**3), d_fwn(n**3), d_host(n**3))
p = 0
do k = 1, n ; do j = 1, n ; do i = 1, n                     ! x fastest, z slowest: slice k is contiguous
   p = p + 1
   points(p) = (-0.6_R8P + 1.2_R8P * (i - 1) / (n - 1)) * ex_R8P + (-0.6_R8P + 1.2_R8P * (j - 1) / (n - 1)) * ey_R8P + &
               (-0.6_R8P + 1.2_R8P * (k - 1) / (n - 1)) * ez_R8P
enddo ; enddo ; enddo
call surface%distance_many_device(points=points, distances=d_pn, is_signed=.true., is_square_root=.true., &
                                  sign_algorithm=SIGN_PSEUDO_NORMAL, status=status)
call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_PSEUDO_NORMAL)
print '(A,I0,A,I0)', 'pseudo-normal : device vs host mismatches ', count(d_pn /= d_host), ' of ', size(points)
call surface%distance_many_device(points=points, distances=d_fwn, is_signed=.true., is_square_root=.true., &
                                  sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
call surface%distance_many(points=points, distances=d_host, is_signed=.true., is_square_root=.true., &
                           sign_algorithm=SIGN_FAST_WINDING_NUMBER)
print '(A,I0,A,I0)', 'winding number: device vs host mismatches ', count(d_fwn /= d_host), ' of ', size(points)
print '(A,I0,A,I0)', 'the two signs disagree at ', count(sign(1._R8P, d_pn) /= sign(1._R8P, d_fwn)), &
                     ' points; |d| differs at ', count(abs(d_pn) /= abs(d_fwn))
!endregion batch
!region slices
do k = 1, n                                                 ! the same points, one host call per slice z = const
   s = (k - 1) * n**2
   call surface%distance_many(points=points(s + 1:s + n**2), distances=d_host(s + 1:s + n**2), is_signed=.true., &
                              is_square_root=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
enddo
print '(A,I0,A,I0)', '24 host calls of one slice vs one device batch: mismatches ', count(d_host /= d_fwn), &
                     ' of ', size(points)
!endregion slices
!region refused
call surface%distance_many_device(points=points, distances=d_host, is_signed=.true., &
                                  sign_algorithm=SIGN_RAY_INTERSECTIONS, status=status)  ! not ported: refused
print '(A,L1)', 'ray-intersection sign on the device refused: ', status /= STATUS_OK
call octree%load_from_file(file_name='dragon-placed.stl', guess_format=.true., aabb_tree_kind=AABB_TREE_OCTREE)
call octree%enter_device                                    ! the device traverses the SAH BVH only: refused
print '(A,L1)', 'octree on the device: ', octree%is_on_device
call surface%exit_device                                    ! release the device copies, before destroy or analyze
print '(A,L1)', 'after exit_device: ', surface%is_on_device
!endregion refused
call draw_figures
contains
   subroutine draw_figures
   !< The figures of the step: the plane nearest z = 0 of the winding-number batch, the points where the two device
   !< signs disagree, and four slices of the batch computed by four host calls, then by one device batch (cells shaded
   !< by the distance, saturated at 0.12).
   use fossil_fig, only : fig_open, fig_surface, fig_quad, fig_segments, fig_points, fig_label, fig_anim, fig_frame, &
                          fig_close
   integer(I4P), parameter       :: ks(4) = [9_I4P, 11_I4P, 13_I4P, 15_I4P] ! the slices of the animation
   type(vector_R8P)              :: e(4), c, lo, hi
   type(vector_R8P), allocatable :: a(:), b(:), q(:), r(:)
   real(R8P)                     :: dx, zc
   character(100)                :: caption
   integer(I4P)                  :: u, kc, f, t

   dx = 1.2_R8P / (n - 1)
   kc = n / 2                                               ! the plane nearest z = 0
   zc = -0.6_R8P + dx * (kc - 1)
   s = (kc - 1) * n**2
   write(caption, '(A,F6.3,A,I0,A,I0,A)') 'distance_many_device, z = ', zc, ': ', &
         count(d_fwn(s + 1:s + n**2) /= d_host(s + 1:s + n**2)), ' of ', n**2, ' differ from the host'
   call fig_open(u, 'tut16-slice', trim(caption), -90._R8P, 90._R8P)
   do p = s + 1, s + n**2                                   ! a cell per point, shaded by its distance (to 0.12)
      call corners(points(p), 0._R8P * ex_R8P, dx, e)
      call fig_quad(u, 'field', e(1), e(2), e(3), e(4), value=d_fwn(p) / 0.12_R8P)
   enddo
   call section_segments(zc, 0._R8P * ex_R8P, a, b)
   call fig_segments(u, 'accent', a, b)
   call fig_close(u)

   q = pack(points, d_pn < 0._R8P .and. d_fwn > 0._R8P)      ! the two device signs disagree
   r = pack(points, d_pn > 0._R8P .and. d_fwn < 0._R8P)
   write(caption, '(A,I0,A,I0,A)') 'inside by the pseudo-normal only: ', size(q), &
                                    ' (ember); by the winding number only: ', size(r), ' (teal)'
   call fig_open(u, 'tut16-signs', trim(caption), -60._R8P, 22._R8P)
   call fig_surface(u, 'ghost', surface)
   if (size(q) > 0) call fig_points(u, 'outside', q)
   if (size(r) > 0) call fig_points(u, 'inside', r)
   call fig_close(u)

   call fig_open(u, 'tut16-slices', 'four slices of the batch: four host calls, then one device batch', &
                 -90._R8P, 90._R8P, height=340_I4P)
   lo = surface%get_bmin() - 0.08_R8P * (ex_R8P + ey_R8P)   ! the cells about the dragon only
   hi = surface%get_bmax() + 0.08_R8P * (ex_R8P + ey_R8P)
   call fig_anim(u, 1.2_R8P)
   do t = 1, 4                                              ! the sections of the four slices, in every frame
      call section_segments(-0.6_R8P + dx * (ks(t) - 1), panel(t), a, b)
      call fig_segments(u, 'accent', a, b)
      c = panel(t) + lo%x * ex_R8P + (hi%y + 0.04_R8P) * ey_R8P
      write(caption, '(A,F6.3)') 'z = ', -0.6_R8P + dx * (ks(t) - 1)
      call fig_label(u, 'box', c, trim(caption))
   enddo
   do f = 1, 5                                              ! frames 1-4: host call f; frame 5: the device batch
      if (f < 5) then
         write(caption, '(A,I0,A,I0,A)') 'host: distance_many on slice ', f, ' of 4 (', n**2, ' points)'
      else
         write(caption, '(A,I0,A,I0,A)') 'device: one distance_many_device of ', 4 * n**2, ' points, ', &
               count(d_host /= d_fwn), ' differ from the host calls'
      endif
      call fig_frame(u, trim(caption), merge(1.2_R8P, 2.4_R8P, f < 5))
      do t = 1, 4
         if (f < 5 .and. t /= f) cycle
         s = (ks(t) - 1) * n**2
         do p = s + 1, s + n**2
            if (points(p)%x < lo%x .or. points(p)%x > hi%x .or. points(p)%y < lo%y .or. points(p)%y > hi%y) cycle
            call corners(points(p), panel(t), dx, e)
            if (f < 5) then
               call fig_quad(u, 'field', e(1), e(2), e(3), e(4), value=d_host(p) / 0.12_R8P)
            else
               call fig_quad(u, 'field', e(1), e(2), e(3), e(4), value=d_fwn(p) / 0.12_R8P)
            endif
         enddo
      enddo
   enddo
   call fig_close(u)
   endsubroutine draw_figures

   function panel(t) result(o)
   !< Offset of panel t of the montage, the four panels in a row.
   integer(I4P), intent(in) :: t !< Panel.
   type(vector_R8P)         :: o !< Offset.
   type(vector_R8P)         :: w !< Panel size and gap.

   w = surface%get_bmax() - surface%get_bmin() + 0.3_R8P * (ex_R8P + ey_R8P)
   o = (w%x * (t - 1)) * ex_R8P
   endfunction panel

   subroutine corners(pt, o, dx, e)
   !< The corners of the cell of size dx about the point pt in its plane z = const, flattened to z = 0, shifted by o.
   type(vector_R8P), intent(in)  :: pt, o !< Point, offset.
   real(R8P),        intent(in)  :: dx    !< Cell size.
   type(vector_R8P), intent(out) :: e(4)  !< Corners.

   e = pt + o ; e%z = 0._R8P
   e(1)%x = e(1)%x - 0.5_R8P * dx ; e(1)%y = e(1)%y - 0.5_R8P * dx
   e(2)%x = e(2)%x + 0.5_R8P * dx ; e(2)%y = e(2)%y - 0.5_R8P * dx
   e(3)%x = e(3)%x + 0.5_R8P * dx ; e(3)%y = e(3)%y + 0.5_R8P * dx
   e(4)%x = e(4)%x - 0.5_R8P * dx ; e(4)%y = e(4)%y + 0.5_R8P * dx
   endsubroutine corners

   subroutine section_segments(z, o, a, b)
   !< The dragon cut by the plane z, as segments drawn just above z = 0, shifted by o.
   use fossil, only : curve_object, slice_surface
   real(R8P),                     intent(in)  :: z    !< Plane.
   type(vector_R8P),              intent(in)  :: o    !< Offset.
   type(vector_R8P), allocatable, intent(out) :: a(:), b(:) !< Segments.
   type(curve_object)                         :: cut
   real(R8P)                                  :: uv(2)
   integer(I4P)                               :: v, m, st

   call slice_surface(surface, 3_I4P, z, cut, status=st)
   allocate(a(cut%get_vertices_number()), b(cut%get_vertices_number())) ; m = 0
   do v = 1, cut%get_vertices_number()
      if (cut%get_next(v) == 0_I4P) cycle
      m = m + 1
      uv = cut%get_vertex(v) ; a(m) = uv(1) * ex_R8P + uv(2) * ey_R8P + 1.e-3_R8P * ez_R8P + o
      uv = cut%get_vertex(cut%get_next(v)) ; b(m) = uv(1) * ex_R8P + uv(2) * ey_R8P + 1.e-3_R8P * ez_R8P + o
   enddo
   a = a(1:m) ; b = b(1:m)
   endsubroutine section_segments
endprogram tut16_a_batch
