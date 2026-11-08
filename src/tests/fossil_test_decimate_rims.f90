!< FOSSIL test: decimate coarsens patch rims (issue #81 item 3, plan #91 step 2).
!<
!< Before, two patch-boundary vertices never merged (issue #66), so a body whose vertices are all on rims could not be
!< decimated at all. Now two rim vertices merge along their own rim edge, onto one of its endpoints (the rim coarsens as
!< a curve, its vertices a subset of the input rim's); a corner, where three or more patches meet, never moves.
!<   1. cube-01 with one patch per face (8 corners where 3 patches meet, 12 rims) to 12 facets: the 12-facet cube, its
!<      volume 1, the 8 corners kept, every facet still on its face. Before: NO_PROGRESS at 36 facets (rims locked, only
!<      the face centres went).
!<   2. A 64-gon prism, caps and side in three patches (every vertex on a rim) to 64 facets: OK, watertight, one shell,
!<      the three patches kept, every output vertex an input vertex (bit for bit), cap vertices still on their cap, the
!<      volume within a measured bound. Before: NO_PROGRESS at 252 facets.
!<   3. The same prism with `rim_collapses=.false.`: the old rule (NO_PROGRESS, nothing collapses).
!<   4. `max_error`: the prism to 8 facets with max_error = 0.02 stops with NO_PROGRESS, and the certified Hausdorff
!<      distance to the input is <= max_error (1.46e-2 at 102 facets, certified to 5e-3); the cube of case 1 with
!<      max_error = 1e-9 (flat faces, every collapse of case 1 costs 0) still reaches 12. The open NACA 0012 wing (188
!<      facets) with max_error = 1e-3 of its diagonal (1.42e-3): Hausdorff <= max_error (1.25e-3 at 101 facets,
!<      certified to 4e-4 of the diagonal; 1.03e-3 to 2.5e-4); without the constraint planes along its open edges,
!<      which max_error adds, the distance was 1.65 (a tip vertex slid along the span).

program fossil_test_decimate_rims

use fossil, only : surface_stl_object, facet_object, DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P),    parameter         :: PI = acos(-1._R8P)
integer(I4P), parameter         :: M  = 64_I4P                     !< Prism sides.
type(surface_stl_object)        :: cube, prism, wing, dec
type(vector_R8P)                :: diagonal
integer(I4P)                    :: status, dev_status
real(R8P)                       :: v0, upper
logical                         :: are_tests_passed(4)

are_tests_passed = .false.

! ----- 1. cube-01, one patch per face, to 12 facets.
call face_patches(cube)
dec = cube
call dec%decimate(target_facets=12_I4P, status=status)
call report('cube-01, 6 patches, to 12', dec, status)
print '(A,ES10.2,A,L1,A,L1,A,I0)', '  volume - 1 ', dec%get_volume() - 1._R8P, ', corners kept ', corners_kept(dec), &
      ', facets on their faces ', on_faces(dec), ', patches ', dec%get_patches_number()
are_tests_passed(1) = status == DEC_STATUS_OK .and. dec%get_facets_number() == 12_I4P .and. valid_closed(dec) .and. &
                      abs(dec%get_volume() - 1._R8P) <= 1.e-14_R8P .and. corners_kept(dec) .and. on_faces(dec) .and. &
                      dec%get_patches_number() == cube%get_patches_number()

! ----- 2. The prism to 64 facets.
call make_prism(prism)
v0 = prism%get_volume()
dec = prism
call dec%decimate(target_facets=64_I4P, status=status)
call report('prism to 64', dec, status)
print '(A,ES10.2)', '  relative volume change ', (dec%get_volume() - v0) / v0
are_tests_passed(2) = status == DEC_STATUS_OK .and. dec%get_facets_number() <= 64_I4P .and. valid_closed(dec) .and. &
                      vertices_subset(dec, prism) .and. caps_ok(dec) .and. &
                      abs(dec%get_volume() - v0) <= 0.05_R8P * v0

! ----- 3. The prism with the rims locked.
dec = prism
call dec%decimate(target_facets=64_I4P, rim_collapses=.false., status=status)
call report('prism to 64, rims locked', dec, status)
are_tests_passed(3) = status == DEC_STATUS_NO_PROGRESS .and. dec%get_facets_number() == prism%get_facets_number()

! ----- 4. max_error.
dec = prism
call dec%decimate(target_facets=8_I4P, max_error=0.02_R8P, status=status)
call report('prism to 8, max_error 0.02', dec, status)
call dec%hausdorff(prism, upper=upper, tolerance=5.e-3_R8P, status=dev_status)
print '(A,ES10.3)', '  certified Hausdorff distance to the input <= ', upper
are_tests_passed(4) = status == DEC_STATUS_NO_PROGRESS .and. dev_status == DEV_STATUS_OK .and. upper <= 0.02_R8P .and. &
                      valid_closed(dec) .and. vertices_subset(dec, prism) .and. caps_ok(dec)
dec = cube
call dec%decimate(target_facets=12_I4P, max_error=1.e-9_R8P, status=status)
call report('cube-01, 6 patches, to 12, max_error 1e-9', dec, status)
are_tests_passed(4) = are_tests_passed(4) .and. status == DEC_STATUS_OK .and. dec%get_facets_number() == 12_I4P
call wing%load_from_file(file_name='src/tests/naca0012-ascii.stl', guess_format=.true.)
diagonal = wing%get_bmax() - wing%get_bmin()
dec = wing
call dec%decimate(target_facets=4_I4P, max_error=1.e-3_R8P * diagonal%normL2(), status=status)
call report('NACA 0012 wing, max_error 1e-3 of the diagonal', dec, status)
call dec%hausdorff(wing, upper=upper, tolerance=4.e-4_R8P * diagonal%normL2(), status=dev_status)
print '(A,ES10.3,A,ES10.3)', '  certified Hausdorff distance to the input <= ', upper, ', max_error ', &
      1.e-3_R8P * diagonal%normL2()
are_tests_passed(4) = are_tests_passed(4) .and. status == DEC_STATUS_NO_PROGRESS .and. dev_status == DEV_STATUS_OK .and. &
                      upper <= 1.e-3_R8P * diagonal%normL2()

print '(A,4L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   subroutine face_patches(s)
   !< cube-01 with patches x0, x1, y0, y1, z0, z1 (the face each facet lies on).
   type(surface_stl_object), intent(inout) :: s
   character(2), parameter                 :: NAMES(6) = ['x0', 'x1', 'y0', 'y1', 'z0', 'z1']
   integer(I4P)                            :: lists(48, 6), n(6), k, a, side

   call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
   n = 0
   do k = 1, s%get_facets_number()
      associate(f => s%facet_at(k))
         do a = 1, 3
            if (abs(f%normal%x * merge(1, 0, a == 1) + f%normal%y * merge(1, 0, a == 2) + &
                    f%normal%z * merge(1, 0, a == 3)) > 0.5_R8P) exit
         enddo
         side = 2 * a - 1
         if (f%normal%x + f%normal%y + f%normal%z > 0._R8P) side = side + 1
      endassociate
      n(side) = n(side) + 1 ; lists(n(side), side) = k
   enddo
   do side = 1, 6
      call s%set_patch(facets=lists(1:n(side), side), name=NAMES(side))
   enddo
   endsubroutine face_patches

   subroutine make_prism(s)
   !< A regular M-gon of radius 1 extruded over [0, 1], caps fanned from vertex 1; patches top, bottom, side.
   type(surface_stl_object), intent(inout) :: s
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: x(M), y(M)
   integer(I4P)                            :: k, j, n

   allocate(f(4 * M - 4))
   x = [(cos(2 * PI * (k - 1) / M), k = 1, M)] ; y = [(sin(2 * PI * (k - 1) / M), k = 1, M)]
   n = 0
   do k = 2, M - 1
      n = n + 1 ; call tri(f(n), [x(1), y(1), 1._R8P], [x(k), y(k), 1._R8P], [x(k+1), y(k+1), 1._R8P])
      n = n + 1 ; call tri(f(n), [x(1), y(1), 0._R8P], [x(k+1), y(k+1), 0._R8P], [x(k), y(k), 0._R8P])
   enddo
   do k = 1, M
      j = mod(k, M) + 1
      n = n + 1 ; call tri(f(n), [x(k), y(k), 0._R8P], [x(j), y(j), 0._R8P], [x(j), y(j), 1._R8P])
      n = n + 1 ; call tri(f(n), [x(k), y(k), 0._R8P], [x(j), y(j), 1._R8P], [x(k), y(k), 1._R8P])
   enddo
   call s%adopt_facets(facets=f)
   call s%set_patch(facets=[(2 * k - 1, k = 1, M - 2)], name='top')
   call s%set_patch(facets=[(2 * k, k = 1, M - 2)], name='bottom')
   call s%set_patch(facets=[(k, k = 2 * M - 3, 4 * M - 4)], name='side')
   endsubroutine make_prism

   subroutine tri(f, p1, p2, p3)
   !< A facet from three vertices.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p1(3), p2(3), p3(3)

   f%vertex(1)%x = p1(1) ; f%vertex(1)%y = p1(2) ; f%vertex(1)%z = p1(3)
   f%vertex(2)%x = p2(1) ; f%vertex(2)%y = p2(2) ; f%vertex(2)%z = p2(3)
   f%vertex(3)%x = p3(1) ; f%vertex(3)%y = p3(2) ; f%vertex(3)%z = p3(3)
   call f%compute_metrix
   endsubroutine tri

   function valid_closed(s) result(yes)
   !< Watertight, manifold, one shell.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes

   yes = s%is_watertight() .and. s%is_manifold() .and. s%get_shells_number() == 1_I4P
   endfunction valid_closed

   function corners_kept(s) result(yes)
   !< The 8 corners of the unit cube are vertices of s.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes
   integer(I4P)                         :: c, k, j
   real(R8P)                            :: p(3)
   logical                              :: found

   yes = .true.
   do c = 0, 7
      p = real([mod(c, 2), mod(c / 2, 2), c / 4], R8P)
      found = .false.
      do k = 1, s%get_facets_number()
         associate(f => s%facet_at(k))
            do j = 1, 3
               found = found .or. (f%vertex(j)%x == p(1) .and. f%vertex(j)%y == p(2) .and. f%vertex(j)%z == p(3))
            enddo
         endassociate
      enddo
      yes = yes .and. found
   enddo
   endfunction corners_kept

   function on_faces(s) result(yes)
   !< Every facet's three vertices lie on the face of its patch (x0: x = 0, x1: x = 1, ...).
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes
   integer(I4P)                         :: k, j, a
   real(R8P)                            :: c, p(3)
   character(:), allocatable            :: name

   yes = .true.
   do k = 1, s%get_facets_number()
      associate(f => s%facet_at(k))
         name = s%patch_name(f%patch)
         a = index('xyz', name(1:1)) ; c = merge(1._R8P, 0._R8P, name(2:2) == '1')
         do j = 1, 3
            p = [f%vertex(j)%x, f%vertex(j)%y, f%vertex(j)%z]
            yes = yes .and. p(a) == c
         enddo
      endassociate
   enddo
   endfunction on_faces

   function vertices_subset(s, s0) result(yes)
   !< Every vertex of s is a vertex of s0, bit for bit.
   type(surface_stl_object), intent(in) :: s, s0
   logical                              :: yes
   integer(I4P)                         :: k, j, k0, j0
   logical                              :: found

   yes = .true.
   do k = 1, s%get_facets_number()
      associate(f => s%facet_at(k))
         do j = 1, 3
            found = .false.
            do k0 = 1, s0%get_facets_number()
               associate(g => s0%facet_at(k0))
                  do j0 = 1, 3
                     found = found .or. (f%vertex(j)%x == g%vertex(j0)%x .and. f%vertex(j)%y == g%vertex(j0)%y .and. &
                                         f%vertex(j)%z == g%vertex(j0)%z)
                  enddo
               endassociate
               if (found) exit
            enddo
            yes = yes .and. found
         enddo
      endassociate
   enddo
   endfunction vertices_subset

   function caps_ok(s) result(yes)
   !< Facets of top have z = 1, of bottom z = 0, at all three vertices; the three patches are all there.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes
   integer(I4P)                         :: k, j, nf(3)
   character(:), allocatable            :: name

   yes = .true. ; nf = 0
   do k = 1, s%get_facets_number()
      associate(f => s%facet_at(k))
         name = s%patch_name(f%patch)
         do j = 1, 3
            if (name == 'top') yes = yes .and. f%vertex(j)%z == 1._R8P
            if (name == 'bottom') yes = yes .and. f%vertex(j)%z == 0._R8P
         enddo
         if (name == 'top') nf(1) = nf(1) + 1
         if (name == 'bottom') nf(2) = nf(2) + 1
         if (name == 'side') nf(3) = nf(3) + 1
      endassociate
   enddo
   yes = yes .and. all(nf > 0)
   endfunction caps_ok

   subroutine report(what, s, status)
   !< One line per case.
   character(*),             intent(in) :: what
   type(surface_stl_object), intent(in) :: s
   integer(I4P),             intent(in) :: status

   print '(2A,I0,A,I0,A,F16.12,A,I0,A,L1,A,L1)', what, ': status ', status, ', facets ', s%get_facets_number(), &
         ', volume ', s%get_volume(), ', shells ', s%get_shells_number(), ', watertight ', s%is_watertight(), &
         ', manifold ', s%is_manifold()
   endsubroutine report
endprogram fossil_test_decimate_rims
