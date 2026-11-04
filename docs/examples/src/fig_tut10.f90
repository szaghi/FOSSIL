!run fig10 fig_tut10
!fig tut10-cells2d
!fig tut10-cells3d
!< Figures of tutorial 10. The 2-D cut cells of the dragon section, each cell shaded by its fluid fraction alpha (solid
!< cells cream, fluid cells sepia), the section (the wall) on top, the cells with two fluid pieces outlined. The 3-D
!< cells of the dragon below z = 0 (the plane of the 2-D section) that are not all fluid, as cubes shaded by alpha
!< (flagged cells, pieces = -1, in ember), and the dragon above z = 0 as a ghost (clip keeps the facets that reach the
!< box, so a few go just below the plane). The figure programs run before the tutorial ones, so the placed dragon
!< (chapters 2-3) and its section (chapter 7) are rebuilt here with the same calls.
program fig_tut10
use fossil,     only : surface_stl_object, curve_object, cut_cells_2d_t, cut_cells_3d_t, cut_cells_3d, slice_surface
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_quad, fig_segments, fig_box, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ez_R8P

implicit none
real(R8P), parameter          :: pi = 4._R8P * atan(1._R8P)
integer(I4P), parameter       :: KC = 12_I4P       ! the 3-D cells drawn: layers 1-12, below z = 0
real(R8P), parameter          :: LIFT = 1.e-3_R8P  ! the wall and outlines drawn just above the cells (seen from +z)
type(surface_stl_object)      :: surface, top
type(curve_object)            :: section
type(cut_cells_2d_t)          :: cut2
type(cut_cells_3d_t)          :: cut3
type(vector_R8P)              :: bmin, bmax, centre
type(vector_R8P), allocatable :: a(:), b(:)
real(R8P)                     :: x(49), scale, uv(2), r
integer(I4P)                  :: u, i, j, k, n, nv, open_loops, status

call surface%load_from_file(file_name='dragon.stl', guess_format=.true.)          ! chapter 2: repair
call surface%sanitize(fill_holes=.true.)
call surface%save_into_file(file_name='fig10-dragon-clean.stl')
call surface%load_from_file(file_name='fig10-dragon-clean.stl', guess_format=.true.) ! chapter 3: place
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=pi / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
call surface%save_into_file(file_name='fig10-dragon-placed.stl')
call surface%load_from_file(file_name='fig10-dragon-placed.stl', guess_format=.true.)
call slice_surface(surface, 3_I4P, 0._R8P, section, merge_collinear=.true., open_loops=open_loops, status=status)
call section%save_into_file(file_name='fig10-dragon-section.dat')                   ! chapter 7: the section
call section%load_from_file(file_name='fig10-dragon-section.dat', status=status)

x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]                                          ! chapter 10: the cut cells
call section%cut_cells(x=x, y=x, cut=cut2)

call fig_open(u, 'tut10-cells2d', 'cut_cells: fluid fraction (cream 0, sepia 1), the wall, cells with 2 pieces', &
              0._R8P, 90._R8P)
bmin = section_box(1) ; bmax = section_box(2)
centre = 0.5_R8P * (bmin + bmax)
r = max(0.5_R8P * (bmax%x - bmin%x), 0.3_R8P * (bmax%y - bmin%y)) + 0.03_R8P      ! x across, y along the card
call fig_zoom(u, centre, r)
do j = 1, 48 ; do i = 1, 48                                                          ! the cells in the card, by alpha
   if (abs(0.5_R8P * (x(i) + x(i + 1)) - centre%x) > r + 0.025_R8P) cycle
   if (abs(0.5_R8P * (x(j) + x(j + 1)) - centre%y) > 1.7_R8P * r + 0.025_R8P) cycle
   call fig_quad(u, 'frac', vec(x(i), x(j)), vec(x(i + 1), x(j)), vec(x(i + 1), x(j + 1)), vec(x(i), x(j + 1)), &
                 value=cut2%alpha(i, j))
enddo ; enddo
nv = section%get_vertices_number()                                                   ! the wall
allocate(a(nv), b(nv)) ; n = 0
do i = 1, nv
   if (section%get_next(i) == 0_I4P) cycle
   n = n + 1
   uv = section%get_vertex(i) ; a(n) = vec(uv(1), uv(2), LIFT)
   uv = section%get_vertex(section%get_next(i)) ; b(n) = vec(uv(1), uv(2), LIFT)
enddo
call fig_segments(u, 'hole', a(1:n), b(1:n))
do j = 1, 48 ; do i = 1, 48                                                          ! cells with two fluid pieces
   if (cut2%pieces(i, j) < 2) cycle
   call fig_box(u, 'accent', vec(x(i), x(j), LIFT), vec(x(i + 1), x(j + 1), LIFT))
   call fig_label(u, 'accent', vec(x(i) - 0.035_R8P, x(j), LIFT), '2 pieces')
enddo ; enddo
call fig_close(u)

call cut_cells_3d(surface=surface, x=x(1:49:2), y=x(1:49:2), z=x(1:49:2), cut=cut3)   ! 24^3 cells, h = 0.05
top = surface                                                                       ! the dragon above z = 0
call top%clip(bmin=vec(-1._R8P, -1._R8P, 0._R8P), bmax=vec(1._R8P, 1._R8P, 1._R8P))
call fig_open(u, 'tut10-cells3d', 'cut_cells_3d below z = 0: fluid fraction, solid cells cream, flagged in ember', &
              -60._R8P, 50._R8P)
call fig_surface(u, 'ghost', top)
do k = 1, KC ; do j = 1, 24 ; do i = 1, 24                                          ! the cells below z = 0
   if (.not. is_drawn(i, j, k)) cycle
   if (.not. is_drawn(i - 1, j, k)) call cube_face(i, j, k, 1)
   if (.not. is_drawn(i + 1, j, k)) call cube_face(i, j, k, 2)
   if (.not. is_drawn(i, j - 1, k)) call cube_face(i, j, k, 3)
   if (.not. is_drawn(i, j + 1, k)) call cube_face(i, j, k, 4)
   if (.not. is_drawn(i, j, k - 1)) call cube_face(i, j, k, 5)
   if (.not. is_drawn(i, j, k + 1)) call cube_face(i, j, k, 6)
enddo ; enddo ; enddo
call fig_close(u)
contains
   pure function is_drawn(i, j, k) result(yes)
   !< The 3-D cell (i, j, k) is below z = 0 and not all fluid: cut, solid or flagged (pieces = -1, alpha = -1).
   integer(I4P), intent(in) :: i, j, k !< Cell.
   logical                  :: yes     !< Drawn.

   yes = .false.
   if (min(i, j, k) < 1 .or. max(i, j) > 24 .or. k > KC) return
   yes = cut3%alpha(i, j, k) < 1._R8P
   endfunction is_drawn

   subroutine cube_face(i, j, k, f)
   !< Face f (-x, +x, -y, +y, -z, +z) of the 3-D cell (i, j, k), by alpha, or in ember if the cell is flagged.
   integer(I4P), intent(in) :: i, j, k, f !< Cell and face.
   real(R8P)                :: lo(3), hi(3), q(3, 4)
   integer(I4P)             :: a, b, m

   lo = [x(2 * i - 1), x(2 * j - 1), x(2 * k - 1)] ; hi = [x(2 * i + 1), x(2 * j + 1), x(2 * k + 1)]
   m = (f + 1) / 2 ; a = mod(m, 3) + 1 ; b = mod(m + 1, 3) + 1                     ! normal axis m, face axes a, b
   q(m, :) = merge(lo(m), hi(m), mod(f, 2) == 1)
   q(a, :) = [lo(a), hi(a), hi(a), lo(a)] ; q(b, :) = [lo(b), lo(b), hi(b), hi(b)]
   if (cut3%pieces(i, j, k) < 0) then
      call fig_quad(u, 'patch', vec(q(1, 1), q(2, 1), q(3, 1)), vec(q(1, 2), q(2, 2), q(3, 2)), &
                    vec(q(1, 3), q(2, 3), q(3, 3)), vec(q(1, 4), q(2, 4), q(3, 4)))
   else
      call fig_quad(u, 'frac', vec(q(1, 1), q(2, 1), q(3, 1)), vec(q(1, 2), q(2, 2), q(3, 2)), &
                    vec(q(1, 3), q(2, 3), q(3, 3)), vec(q(1, 4), q(2, 4), q(3, 4)), value=cut3%alpha(i, j, k))
   endif
   endsubroutine cube_face

   function section_box(which) result(v)
   !< Corner of the box of the section (1: min, 2: max).
   integer(I4P), intent(in) :: which !< Corner.
   type(vector_R8P)         :: v     !< Point.
   real(R8P)                :: c(2)  !< In-plane coordinates.

   if (which == 1) then
      c = section%get_bmin()
   else
      c = section%get_bmax()
   endif
   v = vec(c(1), c(2))
   endfunction section_box

   function vec(p, q, z) result(v)
   !< The point (p, q, z), z = 0 by default.
   real(R8P), intent(in)           :: p, q !< In-plane coordinates.
   real(R8P), intent(in), optional :: z    !< Height.
   type(vector_R8P)                :: v    !< Point.

   v%x = p ; v%y = q ; v%z = 0._R8P
   if (present(z)) v%z = z
   endfunction vec
endprogram fig_tut10
