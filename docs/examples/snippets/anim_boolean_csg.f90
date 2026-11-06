!< Animation of the reference page "Boolean operations" (issue #63): the unit cube drilled by a 3 x 3 grid of
!< capped cylinders (marching cubes of a signed field, closed by sanitize), one difference per frame. Each frame
!< checks its step against an independent boolean: the volume removed must be V(cube n C), the cylinder's part
!< inside the original cube (the cylinders are disjoint); the caption gives that residual relative to V(cube).
program anim_boolean_csg
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, BOOL_INTERSECT, BOOL_DIFFERENCE, &
                       BOOL_STATUS_OK, MC_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: cube, r, c, ic
character(96)            :: caption
real(R8P)                :: v_prev, res, worst
integer(I4P)             :: u, i, j, st

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
r = cube
call fig_open(u, 'anim-boolean-csg', '', -55._R8P, 30._R8P)
call fig_anim(u, 0.9_R8P)
call fig_frame(u, caption='the unit cube: nine differences follow, each checked by an independent boolean')
call fig_surface(u, 'body', r)
worst = 0._R8P
do j = 1, 3 ; do i = 1, 3
   call cylinder(0.2_R8P + 0.3_R8P * (i - 1), 0.2_R8P + 0.3_R8P * (j - 1), c)
   ic = cube ; call ic%boolean(other=c, op=BOOL_INTERSECT, status=st) ; if (st /= BOOL_STATUS_OK) error stop 'check'
   v_prev = r%get_volume()
   call r%boolean(other=c, op=BOOL_DIFFERENCE, status=st) ; if (st /= BOOL_STATUS_OK) error stop 'difference failed'
   res = abs(r%get_volume() - (v_prev - ic%get_volume())) / cube%get_volume() ; worst = max(worst, res)
   write(caption, '(A,I0,A,I0,A,ES7.1,A,L1)') 'step ', i + 3 * (j - 1), ': ', r%get_facets_number(), &
                                              ' facets, volume residual ', res, ', watertight ', r%is_watertight()
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption))
   call fig_surface(u, 'body', r)
enddo ; enddo
print '(A,ES7.1)', 'worst step residual ', worst
call fig_close(u)

contains
   subroutine cylinder(x0, y0, s)
   !< A capped cylinder of radius 0.08 along z through the cube, centred at (x0, y0): the zero level set of
   !< max(distance from the axis - 0.08, |z - 0.5| - 0.6) on 8^3 samples of its box.
   real(R8P),                intent(in)    :: x0, y0 !< Axis position.
   type(surface_stl_object), intent(inout) :: s      !< Cylinder.
   type(surface_stl_object)                :: fresh  !< Reset.
   type(facet_object), allocatable         :: f(:)   !< Facets.
   real(R8P)                               :: phi(8, 8, 8), lo(3), hi(3), x(3)
   integer(I4P)                            :: a, b, k, st

   lo = [x0 - 0.1_R8P, y0 - 0.1_R8P, -0.15_R8P] ; hi = [x0 + 0.1_R8P, y0 + 0.1_R8P, 1.15_R8P]
   do k = 1, 8 ; do b = 1, 8 ; do a = 1, 8
      x = lo + (hi - lo) * [a - 1, b - 1, k - 1] / 7._R8P
      phi(a, b, k) = max(sqrt((x(1) - x0) ** 2 + (x(2) - y0) ** 2) - 0.08_R8P, abs(x(3) - 0.5_R8P) - 0.6_R8P)
   enddo ; enddo ; enddo
   call extract_isosurface(values=phi, bmin=vector_R8P(lo(1), lo(2), lo(3)), bmax=vector_R8P(hi(1), hi(2), hi(3)), &
                           iso=0._R8P, surface=f, status=st)
   if (st /= MC_STATUS_OK) error stop 'marching cubes failed'
   s = fresh
   call s%adopt_facets(facets=f)
   call s%sanitize(fill_holes=.true.)
   if (.not. s%is_watertight()) error stop 'cylinder not closed'
   endsubroutine cylinder
endprogram anim_boolean_csg
