!run figref-booleans-ops fig_ref_booleans_ops
!fig ref-booleans-ops
!fig ref-booleans-patches
!< Figures of the reference page "Boolean operations": A = the unit cube, B = a sphere of radius 0.6 about the corner
!< (1, 1, 1) (marching cubes on 22^3 samples). First the four operations side by side, each with its volume; the
!< symmetric difference cut open by the plane x + y = 1.6, so that the void between its two parts shows. Then the
!< arrangement and the selection: the sub-triangles of A and of B, split along the cut curve (red) and tagged inside or
!< outside the other body. They are read back from the results: A's sub-triangles inside B are those of A ∩ B that lie
!< on A, A's outside ones those of A \ B on A; B's from A ∩ B and A ∪ B. B is drawn moved away along (1, 1, 1).
program fig_ref_booleans_ops
use fossil,     only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, BOOL_UNION, BOOL_INTERSECT, &
                       BOOL_DIFFERENCE, BOOL_SYMDIFF, BOOL_STATUS_OK
use fossil_fig, only : fig_open, fig_triangle, fig_boundary, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter         :: NG = 22
character(10), parameter        :: OPNAME(4) = ['UNION     ', 'INTERSECT ', 'DIFFERENCE', 'SYMDIFF   ']
integer(I4P), parameter         :: OPS(4) = [BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF]
type(surface_stl_object)        :: a, b, r(4), part
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: right, shift, apart
real(R8P)                       :: phi(NG, NG, NG)
character(120)                  :: caption
integer(I4P)                    :: u, i, j, k, status, nin, nout, bin, bout

call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
do k = 1, NG ; do j = 1, NG ; do i = 1, NG
   phi(i, j, k) = norm2(0.3_R8P + 1.4_R8P / (NG - 1) * ([i, j, k] - 1) - 1._R8P) - 0.6_R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(0.3_R8P, 0.3_R8P, 0.3_R8P), &
                        bmax=vector_R8P(1.7_R8P, 1.7_R8P, 1.7_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call b%adopt_facets(facets=facets)
call b%sanitize
do k = 1, 4
   r(k) = a
   call r(k)%boolean(other=b, op=OPS(k), status=status)
   if (status /= BOOL_STATUS_OK) error stop 'boolean failed'
   print '(A,A,I5,A,F8.5)', trim(OPNAME(k)), ': facets ', r(k)%get_facets_number(), ', volume ', r(k)%get_volume()
enddo

apart = vector_R8P(0.35_R8P, 0.35_R8P, 0.35_R8P)
right = vector_R8P(0.866_R8P, 0.5_R8P, 0._R8P)                     ! the screen's right for azimuth -60
call fig_open(u, 'ref-booleans-ops', &
              'A = unit cube, B = sphere of radius 0.6 about (1, 1, 1); SYMDIFF: B \ A (ember) moved apart', &
              -60._R8P, 25._R8P, 760_I4P, 300_I4P)
do k = 1, 4
   shift = (2.0_R8P * (k - 1)) * right
   do i = 1, r(k)%get_facets_number()
      f => r(k)%facet_at(i)
      if (k == 4) then                                              ! B \ A, the facets not in A \ B, moved apart
         if (.not. in_difference(f)) then
            call fig_triangle(u, 'patch', f%vertex(1) + shift + apart, f%vertex(2) + shift + apart, &
                              f%vertex(3) + shift + apart)
            cycle
         endif
      endif
      call fig_triangle(u, 'body', f%vertex(1) + shift, f%vertex(2) + shift, f%vertex(3) + shift)
   enddo
   write(caption, '(A,A,F6.4)') trim(OPNAME(k)), ' ', r(k)%get_volume()
   call fig_label(u, 'accent', vector_R8P(0._R8P, 0._R8P, -0.25_R8P) + shift, trim(caption))
enddo
call fig_close(u)

apart = vector_R8P(0.55_R8P, 0.55_R8P, 0.55_R8P)
if (allocated(facets)) deallocate(facets)               ! adopt_facets took the sphere's array
nin = 0 ; bin = 0
do i = 1, r(2)%get_facets_number()
   f => r(2)%facet_at(i)
   if (on_a(f)) then ; nin = nin + 1 ; else ; bin = bin + 1 ; endif
enddo
nout = 0 ; bout = 0
do i = 1, r(3)%get_facets_number()
   f => r(3)%facet_at(i)
   if (on_a(f)) nout = nout + 1
enddo
do i = 1, r(1)%get_facets_number()
   f => r(1)%facet_at(i)
   if (.not. on_a(f)) bout = bout + 1
enddo
write(caption, '(A,I0,A,I0,A,I0,A,I0,A)') 'A: ', nin, ' inside B (ember), ', nout, ' outside; B: ', bin, &
      ' inside A (teal), ', bout, ' outside'
print '(A)', trim(caption)
call fig_open(u, 'ref-booleans-patches', trim(caption), -60._R8P, 25._R8P, 640_I4P, 460_I4P)
allocate(facets(nin))
j = 0
do i = 1, r(2)%get_facets_number()                                  ! A ∩ B: A's inside, B's inside
   f => r(2)%facet_at(i)
   if (on_a(f)) then
      call fig_triangle(u, 'patch', f%vertex(1), f%vertex(2), f%vertex(3))
      j = j + 1 ; facets(j) = f
   else
      call fig_triangle(u, 'probe', f%vertex(1) + apart, f%vertex(2) + apart, f%vertex(3) + apart)
   endif
enddo
call part%adopt_facets(facets=facets)
call fig_boundary(u, 'hole', part)                                  ! the cut curve on A
do i = 1, r(3)%get_facets_number()                                  ! A \ B: A's outside
   f => r(3)%facet_at(i)
   if (on_a(f)) call fig_triangle(u, 'body', f%vertex(1), f%vertex(2), f%vertex(3))
enddo
do i = 1, r(1)%get_facets_number()                                  ! A ∪ B: B's outside
   f => r(1)%facet_at(i)
   if (.not. on_a(f)) call fig_triangle(u, 'ghost', f%vertex(1) + apart, f%vertex(2) + apart, f%vertex(3) + apart)
enddo
call fig_label(u, 'hole', vector_R8P(1._R8P, 0.45_R8P, 1._R8P), 'cut curve')
call fig_close(u)
contains
   function in_difference(g) result(yes)
   !< True if the facet of the symmetric difference is a facet of A \ B (same centroid).
   type(facet_object), intent(in) :: g   !< Facet.
   logical                        :: yes !< In A \ B.
   type(facet_object), pointer    :: h
   integer(I4P)                   :: m

   yes = .false.
   do m = 1, r(3)%get_facets_number()
      h => r(3)%facet_at(m)
      if (abs(h%centroid%x - g%centroid%x) + abs(h%centroid%y - g%centroid%y) + abs(h%centroid%z - g%centroid%z) &
          < 1.e-12_R8P) then
         yes = .true. ; return
      endif
   enddo
   endfunction in_difference

   pure function on_a(g) result(yes)
   !< True if the facet lies on the surface of the unit cube: its three vertices on one face plane.
   type(facet_object), intent(in) :: g   !< Facet.
   logical                        :: yes !< On A.
   real(R8P)                      :: c(3, 3)
   integer(I4P)                   :: m

   c(:, 1) = [g%vertex(1)%x, g%vertex(2)%x, g%vertex(3)%x]
   c(:, 2) = [g%vertex(1)%y, g%vertex(2)%y, g%vertex(3)%y]
   c(:, 3) = [g%vertex(1)%z, g%vertex(2)%z, g%vertex(3)%z]
   yes = .false.
   do m = 1, 3
      if (all(abs(c(:, m)) < 1.e-12_R8P) .or. all(abs(c(:, m) - 1._R8P) < 1.e-12_R8P)) yes = .true.
   enddo
   endfunction on_a
endprogram fig_ref_booleans_ops
