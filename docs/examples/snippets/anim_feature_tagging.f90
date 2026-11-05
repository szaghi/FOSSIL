!< Animation of the reference page "Feature edges": the machined part of the figure (a block, a boss, a drilled hole),
!< build_features at the default 30 degrees, then box_feature_many on blocks of 8, 16 and 32 cells per side over its
!< bounding box: the cells touching a feature edge (wireframe, gold) concentrate along the edges and corners as the grid
!< refines, the part drawn translucent.
program anim_feature_tagging
use fossil,                    only : surface_stl_object, facet_object, STATUS_OK
use fossil_boolean,            only : BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK
use fossil_fig,                only : fig_open, fig_anim, fig_frame, fig_surface, fig_box, fig_close
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none
real(R8P), parameter     :: PI = acos(-1._R8P)
integer(I4P), parameter  :: LEVELS(3) = [8_I4P, 16_I4P, 32_I4P]
type(surface_stl_object) :: part, boss, hole
type(vector_R8P)         :: bmin, bmax
real(R8P), allocatable   :: lo(:,:), hi(:,:)
logical,   allocatable   :: touched(:)
real(R8P)                :: o(3), h(3)
character(96)            :: caption
integer(I4P)             :: u, l, n, i, j, k, c, status

call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P])
call box(boss, [1.4_R8P, 0.2_R8P, 0._R8P], [1.9_R8P, 0.8_R8P, 0.9_R8P])
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call cylinder(hole, [0.6_R8P, 0.5_R8P], 0.25_R8P, -0.1_R8P, 0.6_R8P, 24_I4P)
call part%boolean(hole, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
call part%build_features(status=status)
if (status /= STATUS_OK) error stop 'build_features failed'
bmin = part%get_bmin() ; bmax = part%get_bmax()
o = [bmin%x, bmin%y, bmin%z] - 0.05_R8P
call fig_open(u, 'anim-feature-tagging', '', -55._R8P, 30._R8P)
call fig_anim(u, 1.6_R8P)
do l = 1, size(LEVELS)
   n = LEVELS(l)
   h = ([bmax%x, bmax%y, bmax%z] + 0.05_R8P - o) / n
   if (allocated(lo)) deallocate(lo, hi, touched)
   allocate(lo(3, n**3), hi(3, n**3), touched(n**3))
   c = 0
   do k = 1, n ; do j = 1, n ; do i = 1, n
      c = c + 1 ; lo(:, c) = o + h * [i - 1, j - 1, k - 1] ; hi(:, c) = lo(:, c) + h
   enddo ; enddo ; enddo
   call part%box_feature_many(lo, hi, touched, status=status)
   write(caption, '(I0,A,I0,A,I0,A)') n, '^3 cells: ', count(touched), ' touch a feature edge (', &
                                      nint(100._R8P * count(touched) / n**3), ' %)'
   print '(A)', trim(caption)
   call fig_frame(u, caption=trim(caption), seconds=merge(2.5_R8P, 1.6_R8P, l == size(LEVELS)))
   call fig_surface(u, 'ghost', part)
   do c = 1, n**3
      if (touched(c)) call fig_box(u, 'accent', vector_R8P(lo(1, c), lo(2, c), lo(3, c)), &
                                   vector_R8P(hi(1, c), hi(2, c), hi(3, c)))
   enddo
enddo
call fig_close(u)
contains
   subroutine box(s, lo, hi)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)

   call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call s%resize(factor=(hi(1) - lo(1)) * ex_R8P + (hi(2) - lo(2)) * ey_R8P + (hi(3) - lo(3)) * ez_R8P)
   call s%translate(delta=lo(1) * ex_R8P + lo(2) * ey_R8P + lo(3) * ez_R8P)
   call s%analyze
   endsubroutine box

   subroutine cylinder(s, centre, r, z0, z1, m)
   !< A regular m-gon of radius r about centre, extruded over [z0, z1], caps fanned, outward.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: centre(2), r, z0, z1
   integer(I4P),             intent(in)    :: m
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: x(m), y(m)
   integer(I4P)                            :: k, j, n

   x = [(centre(1) + r * cos(2 * PI * (k - 1) / m), k = 1, m)] ; y = [(centre(2) + r * sin(2 * PI * (k - 1) / m), k = 1, m)]
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
   endsubroutine cylinder

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1) = vector_R8P(p1(1), p1(2), p1(3)) ; f%vertex(2) = vector_R8P(p2(1), p2(2), p2(3))
   f%vertex(3) = vector_R8P(p3(1), p3(2), p3(3))
   call f%compute_metrix
   endsubroutine tri
endprogram anim_feature_tagging
