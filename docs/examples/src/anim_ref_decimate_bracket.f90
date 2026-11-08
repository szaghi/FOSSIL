!run anim-ref-decimate-bracket anim_ref_decimate_bracket
!fig anim-ref-decimate-bracket
!< Animation of the reference page "Mesh decimation", section "Patches": the bracket of tutorial 12 with 128-sided
!< prisms (3 patches, 1868 facets) decimated to 1000, 600 and 300 facets, then with max_error = 1e-2, then with the rims
!< locked (rim_collapses=.false.): each frame a fresh decimate of the input, coloured by patch (base sepia, boss teal,
!< bore ember); the caption gives the status, the facets and the volume change.
program anim_ref_decimate_bracket
use fossil,     only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P),    parameter     :: PI = acos(-1._R8P)
integer(I4P), parameter     :: TARGETS(3) = [1000_I4P, 600_I4P, 300_I4P]
type(surface_stl_object)    :: part, boss, bore, dec
type(facet_object), pointer :: g
character(100)              :: caption
integer(I4P)                :: u, k, status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, 128_I4P, 'boss')
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, 128_I4P, 'bore')
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
call fig_open(u, 'anim-ref-decimate-bracket', '', -60._R8P, 30._R8P)
call fig_anim(u, 1.6_R8P)
write(caption, '(A,I0,A)') 'the bracket, 128-sided prisms: ', part%get_facets_number(), ' facets, 3 patches'
print '(A)', trim(caption)
call fig_frame(u, trim(caption))
call by_patch(part)
do k = 1, size(TARGETS)
   dec = part
   call dec%decimate(target_facets=TARGETS(k), status=status)
   write(caption, '(A,I0)') 'target_facets=', TARGETS(k)
   call frame(trim(caption))
enddo
dec = part
call dec%decimate(target_facets=4_I4P, max_error=1.e-2_R8P, status=status)
call frame('max_error=1e-2')
dec = part
call dec%decimate(target_facets=600_I4P, rim_collapses=.false., status=status)
call frame('rim_collapses=.false., 600')
call fig_close(u)
contains
   subroutine frame(what)
   !< One frame of dec: status, facets, volume change.
   character(*), intent(in) :: what

   write(caption, '(2A,I0,A,I0,A,ES9.2)') what, ': status ', status, ', ', dec%get_facets_number(), &
         ' facets, volume change ', dec%get_volume() / part%get_volume() - 1._R8P
   print '(A,A,L1,L2,I2)', trim(caption), '; watertight, manifold, shells ', dec%is_watertight(), dec%is_manifold(), &
         dec%get_shells_number()
   call fig_frame(u, caption=trim(caption), seconds=merge(3._R8P, 1.6_R8P, index(what, 'rim_collapses') > 0))
   call by_patch(dec)
   endsubroutine frame

   subroutine by_patch(s)
   !< The facets of s coloured by patch: base sepia, boss teal, bore ember.
   type(surface_stl_object), intent(in) :: s
   integer(I4P)                         :: j

   do j = 1, s%get_facets_number()
      g => s%facet_at(j)
      select case (s%patch_name(g%patch))
      case ('boss') ; call fig_triangle(u, 'probe', g%vertex(1), g%vertex(2), g%vertex(3))
      case ('bore') ; call fig_triangle(u, 'patch', g%vertex(1), g%vertex(2), g%vertex(3))
      case default  ; call fig_triangle(u, 'body', g%vertex(1), g%vertex(2), g%vertex(3))
      endselect
   enddo
   endsubroutine by_patch

   subroutine box(s, lo, hi, name)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)
   character(*),             intent(in)    :: name
   integer(I4P)                            :: k

   call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call s%resize(factor=(hi(1) - lo(1)) * ex_R8P + (hi(2) - lo(2)) * ey_R8P + (hi(3) - lo(3)) * ez_R8P)
   call s%translate(delta=lo(1) * ex_R8P + lo(2) * ey_R8P + lo(3) * ez_R8P)
   call s%analyze
   call s%set_patch(facets=[(k, k = 1, s%get_facets_number())], name=name)
   endsubroutine box

   subroutine cylinder(s, c, r, z0, z1, m, name)
   !< A regular m-gon of radius r about c, extruded over [z0, z1], caps fanned, outward, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: c(2), r, z0, z1
   integer(I4P),             intent(in)    :: m
   character(*),             intent(in)    :: name
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: x(m), y(m)
   integer(I4P)                            :: k, j, n

   x = [(c(1) + r * cos(2 * PI * (k - 1) / m), k = 1, m)] ; y = [(c(2) + r * sin(2 * PI * (k - 1) / m), k = 1, m)]
   allocate(f(2 * (m - 2) + 2 * m))
   n = 0
   do k = 2, m - 1
      n = n + 1 ; call tri(f(n), [x(1), y(1), z1], [x(k), y(k), z1], [x(k+1), y(k+1), z1])
      n = n + 1 ; call tri(f(n), [x(1), y(1), z0], [x(k+1), y(k+1), z0], [x(k), y(k), z0])
   enddo
   do k = 1, m
      j = mod(k, m) + 1
      n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z0], [x(j), y(j), z1])
      n = n + 1 ; call tri(f(n), [x(k), y(k), z0], [x(j), y(j), z1], [x(k), y(k), z1])
   enddo
   call s%destroy
   call s%adopt_facets(facets=f)
   call s%set_patch(facets=[(k, k = 1, s%get_facets_number())], name=name)
   endsubroutine cylinder

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_metrix
   endsubroutine tri
endprogram anim_ref_decimate_bracket
