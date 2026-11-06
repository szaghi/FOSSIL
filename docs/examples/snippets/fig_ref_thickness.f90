!< Figures of the reference page "Thickness and gaps", ray method (default 60-degree cone, 31 rays):
!< - a tapered plate, 4 x 1, its thickness growing linearly from 0.05 to 0.45 along x, coloured by the ray thickness of
!<   each facet (teal thin, ember thick, saturating at 0.5);
!< - a flat plate under a tilted one (one surface, two shells), the gap between them growing from 0.15 to 0.35, coloured
!<   by the ray gap where something is in front, plain where nothing is (gap huge).
program fig_ref_thickness
use fossil,              only : surface_stl_object, STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_triangle, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P), parameter         :: NX = 24, NY = 6
type(surface_stl_object)        :: s
type(facet_object), allocatable :: f(:), g(:), h(:)
real(R8P), allocatable          :: th(:), gp(:)
real(R8P)                       :: cen(3)
integer(I4P)                    :: status, invalid, m, i

! ---- tapered plate
call plate(0.05_R8P, 0.10_R8P, f, 0._R8P, 0._R8P)                   ! bottom z = -t/2, top z = t/2, t = 0.05 + 0.1 x
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(thickness=th, invalid=invalid, status=status)
if (status /= STATUS_OK .or. invalid /= 0) error stop 'thickness failed'
print '(A,I0,A,F6.4,A,F6.4,A,L1)', 'tapered plate: ', s%get_facets_number(), ' facets, ray thickness ', minval(th), &
                                   ' to ', maxval(th), ', watertight ', s%is_watertight()
call draw('ref-thickness-ray', th, 'ray thickness, plate tapering 0.05 to 0.45 (teal 0, ember 0.5)')

! ---- two plates
call plate(0.30_R8P, 0._R8P, f, 0._R8P, 0._R8P)                     ! flat: z in [-0.15, 0.15]
call plate(0.30_R8P, 0._R8P, g, 0.45_R8P, 0.05_R8P)                 ! tilted: lifted by 0.45 + 0.05 x
allocate(h(size(f) + size(g)))                                    ! one surface, two shells
do i = 1, size(f)
   h(i) = f(i)
enddo
do i = 1, size(g)
   h(size(f) + i) = g(i)
enddo
call move_alloc(h, f)
call s%adopt_facets(facets=f) ; call s%analyze
call s%thickness(gap=gp, invalid=invalid, status=status)
if (status /= STATUS_OK .or. invalid /= 0) error stop 'gap failed'
print '(A,F6.4,A,F6.4,A,I0,A)', 'two plates: ray gap ', minval(gp), ' to ', maxval(gp, gp < huge(1._R8P)), ', ', &
                                count(gp == huge(1._R8P)), ' facets with nothing in front'
call draw('ref-thickness-gap', gp, 'ray gap, flat plate under a tilted one (teal 0, ember 0.5; plain: none)')
contains
   subroutine plate(t0, t1, ff, lift0, lift1)
   !< The plate [0, 4] x [0, 1], z in [-t/2, t/2] + lift with t = t0 + t1 x, lift = lift0 + lift1 x; NX x NY cells on
   !< the top and bottom faces, the side walls split to match (watertight), facets oriented outward (the body is
   !< convex: away from its centre).
   real(R8P),                       intent(in)  :: t0, t1, lift0, lift1
   type(facet_object), allocatable, intent(out) :: ff(:)
   real(R8P)                                    :: x(0:NX), y(0:NY), zb(0:NX), zt(0:NX)
   integer(I4P)                                 :: i, j

   x = [(4._R8P * i / NX, i = 0, NX)] ; y = [(1._R8P * j / NY, j = 0, NY)]
   zb = lift0 + lift1 * x - 0.5_R8P * (t0 + t1 * x) ; zt = lift0 + lift1 * x + 0.5_R8P * (t0 + t1 * x)
   cen = [2._R8P, 0.5_R8P, 0.5_R8P * (zb(NX / 2) + zt(NX / 2))]
   allocate(ff(4 * NX * NY + 4 * NX + 4 * NY))
   m = 0
   do j = 0, NY - 1
      do i = 0, NX - 1
         call quad(ff, [x(i), y(j), zt(i)], [x(i + 1), y(j), zt(i + 1)], [x(i + 1), y(j + 1), zt(i + 1)], [x(i), y(j + 1), zt(i)])
         call quad(ff, [x(i), y(j), zb(i)], [x(i + 1), y(j), zb(i + 1)], [x(i + 1), y(j + 1), zb(i + 1)], [x(i), y(j + 1), zb(i)])
      enddo
   enddo
   do i = 0, NX - 1
      call quad(ff, [x(i), y(0), zb(i)], [x(i + 1), y(0), zb(i + 1)], [x(i + 1), y(0), zt(i + 1)], [x(i), y(0), zt(i)])
      call quad(ff, [x(i), y(NY), zb(i)], [x(i + 1), y(NY), zb(i + 1)], [x(i + 1), y(NY), zt(i + 1)], [x(i), y(NY), zt(i)])
   enddo
   do j = 0, NY - 1
      call quad(ff, [x(0), y(j), zb(0)], [x(0), y(j + 1), zb(0)], [x(0), y(j + 1), zt(0)], [x(0), y(j), zt(0)])
      call quad(ff, [x(NX), y(j), zb(NX)], [x(NX), y(j + 1), zb(NX)], [x(NX), y(j + 1), zt(NX)], [x(NX), y(j), zt(NX)])
   enddo
   endsubroutine plate

   subroutine quad(ff, p1, p2, p3, p4)
   !< Two facets of the quad p1 p2 p3 p4, appended at ff(m + 1:m + 2).
   type(facet_object), intent(inout) :: ff(:)
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3), p4(3)

   m = m + 1 ; call tri(ff(m), p1, p2, p3)
   m = m + 1 ; call tri(ff(m), p1, p3, p4)
   endsubroutine quad

   subroutine tri(ft, p1, p2, p3)
   !< A facet oriented away from the plate's centre `cen`.
   type(facet_object), intent(inout) :: ft
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)
   real(R8P)                         :: e1(3), e2(3), n(3)

   e1 = p2 - p1 ; e2 = p3 - p1
   n = [e1(2) * e2(3) - e1(3) * e2(2), e1(3) * e2(1) - e1(1) * e2(3), e1(1) * e2(2) - e1(2) * e2(1)]
   ft%vertex(1) = vector_R8P(p1(1), p1(2), p1(3))
   if (dot_product(n, (p1 + p2 + p3) / 3._R8P - cen) >= 0._R8P) then
      ft%vertex(2) = vector_R8P(p2(1), p2(2), p2(3)) ; ft%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   else
      ft%vertex(2) = vector_R8P(p3(1), p3(2), p3(3)) ; ft%vertex(3) = vector_R8P(p2(1), p2(2), p2(3))
   endif
   call ft%compute_metrix
   endsubroutine tri

   subroutine draw(name, v, caption)
   !< The surface, each facet coloured by v from 0 (teal) to 0.5 (ember) and above (field), plain where v is huge.
   character(*), intent(in)    :: name, caption
   real(R8P),    intent(in)    :: v(:)
   type(facet_object), pointer :: fp
   integer(I4P)                :: u, i

   call fig_open(u, name, caption, -60._R8P, 30._R8P)
   do i = 1, s%get_facets_number()
      fp => s%facet_at(i)
      if (v(i) < huge(1._R8P)) then
         call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), value=2._R8P * min(1._R8P, v(i) / 0.5_R8P) - 1._R8P)
      else
         call fig_triangle(u, 'body', fp%vertex(1), fp%vertex(2), fp%vertex(3))
      endif
   enddo
   call fig_close(u)
   endsubroutine draw
endprogram fig_ref_thickness
