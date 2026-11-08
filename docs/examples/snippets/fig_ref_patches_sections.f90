!< Figures of the reference page "Named patches", one per section, every facet in the colour of its patch (patch 1
!< stone, 2 teal, 3 ember, 4 gold, 5 brown). "Reading and writing": the dragon split by x into head, middle and tail
!< (set_patch), saved as ASCII and loaded back: one solid block per patch, the names read from the file. "Through the
!< pipeline": a box remeshed to edge 0.15 with an inlet (z = 0), an outlet (z = 1) and walls, two inlet facets removed,
!< then fill_holes: the hole inside the inlet takes the inlet's patch, the one on the inlet-wall edge the patch
!< "filled". "Through booleans": a block with a lid minus a sphere: the dimple carries the sphere's patch. "Through
!< remeshing": a sphere in a north and a south patch, remeshed and decimated: the zigzag patch boundary stays where it
!< is. "2-D curves": the union of two split boxes (the page's example) cut at z = 0.25 and 0.75, every segment of the
!< sections in the colour of its patch.
program fig_ref_patches_sections
use fossil,     only : surface_stl_object, facet_object, curve_object, slice_surface, extract_isosurface, MC_STATUS_OK, &
                       BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, REM_STATUS_OK, DEC_STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_triangle, fig_surface, fig_segment, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
character(6), parameter         :: STYLE(5) = ['body  ', 'probe ', 'patch ', 'accent', 'cell  ']
type(surface_stl_object)        :: s, t, w
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(curve_object)              :: section
type(vector_R8P)                :: right
integer(I4P),     allocatable   :: ids(:)
real(R8P)                       :: phi(25, 25, 25)
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, n, status, filled, inner, edge, z

! Reading and writing
call s%load_from_file(file_name='dragon.stl', guess_format=.true.)
call tag(s, 'head', 1) ; call tag(s, 'middle', 2) ; call tag(s, 'tail', 3)
call s%save_into_file(file_name='dragon-patches.stl', is_ascii=.true.)
call t%load_from_file(file_name='dragon-patches.stl', guess_format=.true.)
write(caption, '(A,I0,A,I0,A)') 'dragon-patches.stl read back: ', t%get_patches_number(), ' solid blocks, ', &
      t%get_facets_number(), ' facets'
print '(A)', trim(caption)
call fig_open(u, 'ref-patches-files', trim(caption), -60._R8P, 22._R8P)
call draw(t, vector_R8P(0._R8P, 0._R8P, 0._R8P))
call labels(t, 3._R8P)
call fig_close(u)

! Through the pipeline: fill_holes
call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call s%isotropic_remesh(target_length=0.15_R8P, iterations=3_I4P, preserve_features=.true., status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
call by_normal(s)
inner = nearest(s, vector_R8P(0.5_R8P, 0.5_R8P, 0._R8P))
edge = nearest(s, vector_R8P(0.5_R8P, 0.02_R8P, 0._R8P))
allocate(facets(s%get_facets_number() - 2))
j = 0
do i = 1, s%get_facets_number()
   if (i == inner .or. i == edge) cycle
   f => s%facet_at(i) ; j = j + 1 ; facets(j) = f
enddo
call s%adopt_facets(facets=facets)
n = s%get_facets_number()
call s%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled)
write(caption, '(A,I0,A)') 'fill_holes: ', filled, ' holes; the new facets take the inlet or the patch "filled"'
print '(A)', trim(caption)
do i = n + 1, s%get_facets_number()
   print '(A,I0,A,A)', '  new facet ', i, ': patch ', s%patch_name(s%facet_patch(i))
enddo
call fig_open(u, 'ref-patches-fill', trim(caption), -60._R8P, -40._R8P)
call fig_zoom(u, vector_R8P(0.5_R8P, 0.3_R8P, 0._R8P), 0.42_R8P)
call draw(s, vector_R8P(0._R8P, 0._R8P, 0._R8P))
do i = n + 1, s%get_facets_number()                                  ! the new facets outlined
   f => s%facet_at(i)
   call fig_segment(u, 'accent', f%vertex(1), f%vertex(2)) ; call fig_segment(u, 'accent', f%vertex(2), f%vertex(3))
   call fig_segment(u, 'accent', f%vertex(3), f%vertex(1))
   call fig_label(u, 'accent', f%centroid, 'new: '//s%patch_name(s%facet_patch(i)))
enddo
call fig_close(u)

! Through booleans
call s%load_from_file(file_name='cube-01.stl', guess_format=.true.)
allocate(ids(0))
do i = 1, s%get_facets_number()
   f => s%facet_at(i)
   if (f%normal%z > 0.5_R8P) ids = [ids, i]
enddo
call s%set_patch(facets=[(i, i = 1, s%get_facets_number())], name='block')
call s%set_patch(facets=ids, name='lid')
do k = 1, 25 ; do j = 1, 25 ; do i = 1, 25
   phi(i, j, k) = norm2([0.05_R8P, 0.05_R8P, 0.55_R8P] + 0.9_R8P / 24 * ([i, j, k] - 1) - [0.5_R8P, 0.5_R8P, 1._R8P]) - &
                  0.4_R8P
enddo ; enddo ; enddo
call t%destroy
call ball(t, vector_R8P(0.05_R8P, 0.05_R8P, 0.55_R8P), vector_R8P(0.95_R8P, 0.95_R8P, 1.45_R8P))
call t%set_patch(facets=[(i, i = 1, t%get_facets_number())], name='sphere')
call s%boolean(other=t, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
write(caption, '(A,I0,A)') 'block minus sphere: ', s%get_facets_number(), ' facets, the dimple in patch "sphere"'
print '(A)', trim(caption)
call fig_open(u, 'ref-patches-difference', trim(caption), -60._R8P, 40._R8P)
call draw(s, vector_R8P(0._R8P, 0._R8P, 0._R8P))
call labels(s, 0._R8P)
call fig_close(u)

! Through remeshing
do k = 1, 25 ; do j = 1, 25 ; do i = 1, 25
   phi(i, j, k) = norm2(-1.5_R8P + 3._R8P / 24 * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call s%destroy
call ball(s, vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P))
deallocate(ids) ; allocate(ids(0))
do i = 1, s%get_facets_number()
   f => s%facet_at(i)
   if (f%centroid%z > 0.1_R8P * f%centroid%x) ids = [ids, i]
enddo
call s%set_patch(facets=[(i, i = 1, s%get_facets_number())], name='south')
call s%set_patch(facets=ids, name='north')
t = s
call t%isotropic_remesh(target_length=0.08_R8P, iterations=3_I4P, status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
w = s
call w%decimate(target_facets=s%get_facets_number() / 6, status=status)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
write(caption, '(A,I0,A,I0,A,I0,A)') 'input ', s%get_facets_number(), ', isotropic_remesh ', t%get_facets_number(), &
      ', decimate ', w%get_facets_number(), ' facets: the boundary is kept'
print '(A)', trim(caption)
right = vector_R8P(2.4_R8P * 0.866_R8P, 2.4_R8P * 0.5_R8P, 0._R8P)
call fig_open(u, 'ref-patches-remesh', trim(caption), -60._R8P, 15._R8P, 760_I4P, 300_I4P)
call draw(s, vector_R8P(0._R8P, 0._R8P, 0._R8P))
call draw(t, right)
call draw(w, 2._R8P * right)
call fig_close(u)

! 2-D curves
call split_box(s, 'low', 'high', vector_R8P(0._R8P, 0._R8P, 0._R8P))
call split_box(t, 'b_low', 'b_high', vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P))
call s%boolean(other=t, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
call fig_open(u, 'ref-patches-curves', 'sections z = 0.25 and 0.75 of the union: segments coloured by patch', &
              -60._R8P, 30._R8P)
call fig_surface(u, 'ghost', s)
do z = 1, 2
   call slice_surface(surface=s, normal_axis=3_I4P, level=0.5_R8P * z - 0.25_R8P, curve=section)
   print '(A,F4.2,A,I0,A,I0,A)', 'section z = ', 0.5_R8P * z - 0.25_R8P, ': ', section%get_vertices_number(), &
         ' segments, ', section%get_patches_number(), ' patches'
   do i = 1, section%get_vertices_number()
      j = section%get_next(i)
      if (j == 0) cycle
      call fig_segment(u, trim(STYLE(section%segment_patch(i))), point(section%get_vertex(i), 0.5_R8P * z - 0.25_R8P), &
                       point(section%get_vertex(j), 0.5_R8P * z - 0.25_R8P))
   enddo
   do k = 1, section%get_patches_number()
      do i = 1, section%get_vertices_number()
         if (section%segment_patch(i) /= k .or. section%get_next(i) == 0) cycle
         call fig_label(u, trim(STYLE(k)), point(section%get_vertex(i), 0.5_R8P * z - 0.25_R8P), &
                        section%patch_name(k))
         exit
      enddo
   enddo
enddo
call fig_close(u)
contains
   subroutine tag(body, name, part)
   !< Put the facets of the dragon in a third of its x range (part 1, 2, 3) in the patch `name`.
   type(surface_stl_object), intent(inout) :: body !< Surface.
   character(*),             intent(in)    :: name !< Patch name.
   integer(I4P),             intent(in)    :: part !< Third.
   type(facet_object), pointer             :: g
   integer(I4P), allocatable               :: list(:)
   type(vector_R8P)                        :: c0, c1
   real(R8P)                               :: x0, x1
   integer(I4P)                            :: m

   c0 = body%get_bmin() ; c1 = body%get_bmax() ; x0 = c0%x ; x1 = c1%x
   allocate(list(0))
   do m = 1, body%get_facets_number()
      g => body%facet_at(m)
      if (g%centroid%x >= x0 + (x1 - x0) * (part - 1) / 3 .and. g%centroid%x < x0 + (x1 - x0) * part / 3 + &
          merge(1._R8P, 0._R8P, part == 3)) list = [list, m]
   enddo
   call body%set_patch(facets=list, name=name)
   endsubroutine tag

   subroutine by_normal(body)
   !< Patches by outward normal: inlet (z-), outlet (z+), wall (the rest).
   type(surface_stl_object), intent(inout) :: body !< Surface.
   type(facet_object), pointer             :: g
   integer(I4P), allocatable               :: inlet(:), outlet(:)
   integer(I4P)                            :: m

   allocate(inlet(0), outlet(0))
   do m = 1, body%get_facets_number()
      g => body%facet_at(m)
      if (g%normal%z < -0.5_R8P) inlet = [inlet, m]
      if (g%normal%z > 0.5_R8P) outlet = [outlet, m]
   enddo
   call body%set_patch(facets=[(m, m = 1, body%get_facets_number())], name='wall')
   call body%set_patch(facets=inlet, name='inlet')
   call body%set_patch(facets=outlet, name='outlet')
   endsubroutine by_normal

   function nearest(body, p) result(id)
   !< The facet whose centroid is nearest to p.
   type(surface_stl_object), intent(in) :: body !< Surface.
   type(vector_R8P),         intent(in) :: p    !< Point.
   integer(I4P)                         :: id   !< Facet.
   type(facet_object), pointer          :: g
   real(R8P)                            :: best, dd
   integer(I4P)                         :: m

   best = huge(1._R8P) ; id = 0
   do m = 1, body%get_facets_number()
      g => body%facet_at(m)
      dd = norm2([g%centroid%x - p%x, g%centroid%y - p%y, g%centroid%z - p%z])
      if (dd < best) then
         best = dd ; id = m
      endif
   enddo
   endfunction nearest

   subroutine ball(body, bmin, bmax)
   !< The marching-cubes surface of phi on the box [bmin, bmax], sanitized.
   type(surface_stl_object), intent(inout) :: body       !< Surface.
   type(vector_R8P),         intent(in)    :: bmin, bmax !< Box.
   type(facet_object), allocatable         :: tri(:)
   integer(I4P)                            :: st

   call extract_isosurface(values=phi, bmin=bmin, bmax=bmax, iso=0._R8P, surface=tri, status=st)
   if (st /= MC_STATUS_OK) error stop 'marching cubes failed'
   call body%adopt_facets(facets=tri)
   call body%sanitize
   endsubroutine ball

   subroutine draw(body, delta)
   !< The facets of body moved by delta, in the colour of their patch.
   type(surface_stl_object), intent(in) :: body  !< Surface.
   type(vector_R8P),         intent(in) :: delta !< Shift.
   type(facet_object), pointer          :: g
   integer(I4P)                         :: m

   do m = 1, body%get_facets_number()
      g => body%facet_at(m)
      call fig_triangle(u, trim(STYLE(min(5, body%facet_patch(m)))), g%vertex(1) + delta, g%vertex(2) + delta, &
                        g%vertex(3) + delta)
   enddo
   endsubroutine draw

   subroutine labels(body, lift)
   !< The name and facet count of every non-empty patch, at the centroid of its facets lifted by `lift` along z.
   type(surface_stl_object), intent(in) :: body !< Surface.
   real(R8P),                intent(in) :: lift !< Lift of the label.
   type(facet_object), pointer          :: g
   type(vector_R8P)                     :: c
   character(64)                        :: text
   integer(I4P)                         :: m, q, cnt

   do q = 1, body%get_patches_number()
      c = vector_R8P(0._R8P, 0._R8P, 0._R8P) ; cnt = 0
      do m = 1, body%get_facets_number()
         if (body%facet_patch(m) /= q) cycle
         g => body%facet_at(m)
         c = c + g%centroid ; cnt = cnt + 1
      enddo
      if (cnt == 0) cycle
      write(text, '(A,A,I0,A)') body%patch_name(q), ' (', cnt, ')'
      print '(A,A)', '  patch ', trim(text)
      call fig_label(u, trim(STYLE(min(5, q))), c / real(cnt, R8P) + vector_R8P(0._R8P, 0._R8P, lift), trim(text))
   enddo
   endsubroutine labels

   function point(uv, level) result(p)
   !< The 3-D point of a section of the plane z = level.
   real(R8P), intent(in) :: uv(2) !< (x, y).
   real(R8P), intent(in) :: level !< z.
   type(vector_R8P)      :: p     !< Point.

   p = vector_R8P(uv(1), uv(2), level)
   endfunction point

   subroutine split_box(box, lower, upper, shift)
   !< The unit box moved by shift, its facets below its mid-height in patch `lower`, the others in `upper`.
   type(surface_stl_object), intent(inout) :: box          !< Box.
   character(*),             intent(in)    :: lower, upper !< Patch names.
   type(vector_R8P),         intent(in)    :: shift        !< Shift.
   type(facet_object), pointer             :: g
   integer(I4P), allocatable               :: low(:), high(:)
   integer(I4P)                            :: m

   call box%load_from_file(file_name='cube-01.stl', guess_format=.true.)
   call box%translate(delta=shift) ; call box%analyze
   allocate(low(0), high(0))
   do m = 1, box%get_facets_number()
      g => box%facet_at(m)
      if (g%centroid%z < shift%z + 0.5_R8P) then ; low = [low, m] ; else ; high = [high, m] ; endif
   enddo
   call box%set_patch(facets=low, name=lower)
   call box%set_patch(facets=high, name=upper)
   endsubroutine split_box
endprogram fig_ref_patches_sections
