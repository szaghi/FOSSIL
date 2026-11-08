!run anim-ref-usage-sanitize anim_ref_usage_sanitize
!fig anim-ref-usage-sanitize
!< Animation of the usage reference, "Sanitize (the full repair pipeline)": the passes of sanitize called one by one, in
!< its order, on the unit cube with one defect of each class: a degenerate facet (three collinear vertices on the edge
!< y = 0, z = 1), a facet shifted by 1e-4 (a crack along its edges), a duplicate facet, two facets wound inwards and a
!< missing facet (a hole). Each frame draws the boundary edges (red), the facet normals (teal outward, ember inward),
!< the duplicate (ember facet) and the degenerate facet (gold points), with the counters of the pass in the caption.
program anim_ref_usage_sanitize
use fossil,     only : surface_stl_object, facet_object
use fossil_fig, only : fig_open, fig_surface, fig_boundary, fig_segments, fig_triangle, fig_points, fig_label, &
                       fig_anim, fig_frame, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P, ex_R8P

implicit none
type(surface_stl_object)        :: cube
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P), parameter     :: MID = vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)
type(vector_R8P)                :: anchor(4)
character(120)                  :: caption
integer(I4P)                    :: u, i, n, m, k, crack, hole, dup, flip(2), flipped, filled, skipped

call cube%load_from_file(file_name='cube-01.stl', guess_format=.true.)
n = cube%get_facets_number()
crack = 0 ; hole = 0 ; dup = 0 ; flip = 0 ; k = 0
do i = 1, n                                                   ! pick the facets of the defects on the visible faces
   f => cube%facet_at(i)
   if (crack == 0 .and. f%normal%x > 0.5_R8P) crack = i
   if (f%normal%y < -0.5_R8P) then
      if (hole == 0) then
         hole = i
      elseif (k < 2 .and. f%centroid%z > 0.5_R8P) then
         k = k + 1 ; flip(k) = i
      endif
   endif
   if (dup == 0 .and. f%normal%z > 0.5_R8P .and. f%centroid%x > 0.5_R8P) dup = i
enddo
f => cube%facet_at(crack) ; anchor(1) = f%centroid
f => cube%facet_at(hole) ; anchor(2) = f%centroid
f => cube%facet_at(flip(1)) ; anchor(3) = f%centroid
f => cube%facet_at(dup) ; anchor(4) = f%centroid
allocate(facets(n + 1))                                       ! n - 1 kept, the duplicate, the degenerate facet
m = 0
do i = 1, n
   if (i == hole) cycle                                       ! the hole
   f => cube%facet_at(i)
   m = m + 1 ; facets(m) = f
   if (i == crack) then                                       ! the crack
      do k = 1, 3
         facets(m)%vertex(k) = facets(m)%vertex(k) + 1.e-4_R8P * ex_R8P
      enddo
   elseif (i == flip(1) .or. i == flip(2)) then               ! wound inwards
      facets(m)%vertex(2:3) = facets(m)%vertex([3, 2])
   endif
   if (i == dup) then                                         ! the duplicate
      m = m + 1 ; facets(m) = f
   endif
enddo
facets(n + 1) = facets(1)                                     ! the degenerate facet: three collinear vertices
facets(size(facets))%vertex(1) = vector_R8P(0._R8P, 0._R8P, 1._R8P)
facets(size(facets))%vertex(2) = vector_R8P(0.5_R8P, 0._R8P, 1._R8P)
facets(size(facets))%vertex(3) = vector_R8P(1._R8P, 0._R8P, 1._R8P)
do i = 1, size(facets)
   call facets(i)%compute_normal
enddo
call cube%adopt_facets(facets=facets)

call fig_open(u, 'anim-ref-usage-sanitize', '', -55._R8P, 24._R8P)
call fig_anim(u, 2.2_R8P)
write(caption, '(A,I0,A,I0,A)') 'as loaded: ', cube%get_facets_number(), ' facets, ', open_edges(), &
                                ' boundary edges, watertight F'
call frame(trim(caption), labels=anchor)
call cube%remove_degenerate_facets
call cube%analyze
write(caption, '(A,I0,A)') '1. remove_degenerate_facets: ', cube%get_degenerate_facets_removed(), &
                           ' removed (the gold points)'
call frame(trim(caption))
call cube%connect_nearby_vertices
call cube%analyze
write(caption, '(A,I0,A)') '2-3. connect_nearby_vertices, analyze: crack closed, ', open_edges(), &
                           ' boundary edges left'
call frame(trim(caption))
call cube%remove_duplicate_facets
call cube%analyze
write(caption, '(A,I0,A,I0,A)') '4-5. remove_duplicate_facets, analyze: ', cube%get_duplicate_facets_removed(), &
                           ' removed, ', open_edges(), ' boundary edges (the hole)'
call frame(trim(caption))
call cube%sanitize_normals(flipped_count=flipped)
call cube%analyze
write(caption, '(A,I0,A)') '6. sanitize_normals: ', flipped, ' facets flipped, every normal outward'
call frame(trim(caption))
call cube%fill_holes(max_hole_size=0.5_R8P, holes_filled=filled, holes_skipped=skipped)
write(caption, '(A,I0,A,L1,A,F6.4)') '7. fill_holes: ', filled, &
                                     ' filled, watertight ', cube%is_watertight(), ', volume ', cube%get_volume()
call frame(trim(caption), 3.5_R8P)
call fig_close(u)
contains
   function open_edges() result(m)
   !< The number of boundary edges (facet edges without a neighbour).
   integer(I4P) :: m !< Boundary edges.
   integer(I4P) :: j

   m = 0
   do j = 1, cube%get_facets_number()
      f => cube%facet_at(j)
      m = m + count(f%fcon_edge == 0_I4P)
   enddo
   endfunction open_edges

   subroutine frame(text, seconds, labels)
   !< One frame: the cube, its boundary edges, its normals, the duplicate and the degenerate facets.
   character(*),     intent(in)           :: text      !< Caption.
   real(R8P),        intent(in), optional :: seconds   !< Duration.
   type(vector_R8P), intent(in), optional :: labels(4) !< Anchors of the defect labels (first frame).
   type(vector_R8P), allocatable          :: a(:), b(:), dg(:)
   type(facet_object), pointer            :: g
   integer(I4P)                           :: j, l

   print '(A)', text
   call fig_frame(u, caption=text, seconds=seconds)
   call fig_surface(u, 'body', cube)
   call fig_boundary(u, 'hole', cube)
   allocate(a(0), b(0), dg(0))
   do j = 1, cube%get_facets_number()
      f => cube%facet_at(j)
      if (f%area() < 1.e-12_R8P) then                         ! degenerate
         dg = [dg, f%vertex]
         cycle
      endif
      if ((f%normal .dot. (f%centroid - MID)) > 0._R8P) then
         a = [a, f%centroid] ; b = [b, f%centroid + 0.15_R8P * f%normal]
      else                                                    ! wound inwards: the facet in ember
         call fig_triangle(u, 'outside', f%vertex(1) - 0.005_R8P * f%normal, f%vertex(2) - 0.005_R8P * f%normal, &
                           f%vertex(3) - 0.005_R8P * f%normal)
      endif
      do l = j + 1, cube%get_facets_number()                 ! a duplicate: same centroid and area
         g => cube%facet_at(l)
         if (norm2([g%centroid%x - f%centroid%x, g%centroid%y - f%centroid%y, g%centroid%z - f%centroid%z]) < &
             1.e-12_R8P .and. abs(g%area() - f%area()) < 1.e-12_R8P) &
            call fig_triangle(u, 'cell', f%vertex(1) + 0.01_R8P * f%normal, f%vertex(2) + 0.01_R8P * f%normal, &
                              f%vertex(3) + 0.01_R8P * f%normal)
      enddo
   enddo
   if (size(a) > 0) call fig_segments(u, 'probe', a, b)
   if (size(dg) > 0) call fig_points(u, 'accent', dg)
   if (present(labels)) then
      call fig_label(u, 'accent', dg(1), 'degenerate')
      call fig_label(u, 'hole', labels(1), 'crack')
      call fig_label(u, 'hole', labels(2), 'hole')
      call fig_label(u, 'accent', labels(3), 'wound inwards')
      call fig_label(u, 'accent', labels(4), 'duplicate')
   endif
   endsubroutine frame
endprogram anim_ref_usage_sanitize
