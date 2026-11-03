!< FOSSIL test: refinement of hole patches (issue #31, step 5).
program fossil_test_hole_refine
!< FOSSIL test: refinement of hole patches (issue #31, step 5).
!<
!< cube-02 (unit cube, 4 x 4 grid squares per face, two triangles each) with a 3 x 3 block of
!< squares removed from the face z = 1: one planar hole, 12 loop edges of length 0.25, area
!< 0.5625. Filled without and with `refine`:
!<  - both results watertight, volume 1 to 1e-12;
!<  - the refined patch has more triangles than the plain one, covers the hole exactly (patch
!<    area 0.5625 to 1e-12) and stays in the plane z = 1 (bit for bit);
!<  - its edges are near the loop's density: longest patch edge below 2 x 0.25, mean within
!<    [0.5, 1.5] x 0.25.
!< Tiny: 192 facets.

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P

implicit none

integer(I4P) :: n_plain, n_refined
logical      :: ok

ok = .true.
call run(refine=.false., n_patch=n_plain)
call run(refine=.true., n_patch=n_refined)
ok = ok .and. n_refined > n_plain
print '(A,I0,A,I0)', 'patch triangles: plain ', n_plain, '  refined ', n_refined
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine run(refine, n_patch)
   logical,      intent(in)        :: refine
   integer(I4P), intent(out)       :: n_patch
   type(surface_stl_object)        :: cube, holed
   type(facet_object), allocatable :: facets(:)
   type(facet_object), pointer     :: f
   integer(I4P)                    :: i, n_keep, n_orig, filled, e
   real(R8P)                       :: area, lmax, lsum
   logical                         :: in_plane, case_ok

   call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   allocate(facets(cube%get_facets_number()))
   n_keep = 0
   do i = 1, cube%get_facets_number()
      f => cube%facet_at(i)
      ! drop the facets of z = 1 inside the 3 x 3 block [0, 0.75]^2
      if (f%centroid%z > 1._R8P - 1.e-9_R8P .and. f%centroid%x < 0.75_R8P .and. f%centroid%y < 0.75_R8P) cycle
      n_keep = n_keep + 1
      facets(n_keep) = f
      facets(n_keep)%vertex_id = 0_I4P
   enddo
   facets = facets(1:n_keep)
   n_orig = n_keep
   call holed%adopt_facets(facets=facets)
   call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, refine=refine)
   n_patch = holed%get_facets_number() - n_orig
   area = 0._R8P ; lmax = 0._R8P ; lsum = 0._R8P ; in_plane = .true.
   do i = n_orig + 1, holed%get_facets_number()
      f => holed%facet_at(i)
      area = area + f%area()
      do e = 1, 3
         in_plane = in_plane .and. f%vertex(e)%z == 1._R8P
         lmax = max(lmax, edge(f, e))
         lsum = lsum + edge(f, e)
      enddo
   enddo
   case_ok = filled == 1 .and. holed%is_watertight() .and. abs(holed%get_volume() - 1._R8P) < 1.e-12_R8P .and. &
             abs(area - 0.5625_R8P) < 1.e-12_R8P .and. in_plane
   if (refine) case_ok = case_ok .and. lmax < 0.5_R8P .and. lsum / (3 * n_patch) >= 0.125_R8P .and. &
                         lsum / (3 * n_patch) <= 0.375_R8P
   print '(A,L1,A,I0,A,I0,A,L1,A,F14.12,A,F14.12,A,L1,A,F6.3,A,F6.3,A,L1)', 'refine ', refine, ': removed ', &
         cube%get_facets_number() - n_orig, '  patch ', n_patch, '  watertight ', holed%is_watertight(), '  volume ', &
         holed%get_volume(), '  patch area ', area, '  in plane ', in_plane, '  edge max ', lmax, '  mean ', &
         lsum / (3 * n_patch), '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine run

   pure function edge(f, e) result(l)
   type(facet_object), intent(in) :: f
   integer(I4P),       intent(in) :: e
   real(R8P)                      :: l
   integer(I4P)                   :: e2

   e2 = mod(e, 3) + 1
   l = sqrt((f%vertex(e2)%x - f%vertex(e)%x)**2 + (f%vertex(e2)%y - f%vertex(e)%y)**2 + (f%vertex(e2)%z - f%vertex(e)%z)**2)
   endfunction edge
endprogram fossil_test_hole_refine
