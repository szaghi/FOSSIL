!< FOSSIL validation: decimate on a body built by booleans (issue #81, plan #91 step 2).
!<
!< The bracket of tutorial 12 with 128-sided prisms: cube-01 resized to a 2 x 1 x 0.5 base, a boss prism united, a
!< bore prism subtracted (exact booleans, 3 named patches, 1868 facets). Each call on a fresh copy of the input.
!<   1. To 1000 and to 600 facets: OK, valid (watertight, manifold, one shell), the patches kept (same table, every
!<      patch with facets), every output rim vertex an input rim vertex (bit for bit), |dV/V| <= 1e-4 (measured 7.0e-6
!<      and 2.4e-5). Before (rims locked, no ceiling): NO_PROGRESS at 1582 facets for both, dV/V = -36.6 %: the 3 rims
!<      of 128 vertices could not coarsen, the cheap collapses ran out and the loop went on with costly ones (box
!<      corners sheared along the box edges) until the heap was empty.
!<   2. max_error = 3e-3 and 1e-2 (target 4): NO_PROGRESS, valid, the certified Hausdorff distance to the input (to
!<      max_error / 4) <= max_error (measured 1.8e-3 and 6.3e-3 at 458 and 294 facets).
!<   3. rim_collapses=.false. to 600: the old rule; NO_PROGRESS, valid (printed: 1580 facets, dV/V = -41.8 %, the cost of
!<      running on with no ceiling).

program fossil_val_decimate_bracket

use fossil, only : surface_stl_object, facet_object, BOOL_UNION, BOOL_DIFFERENCE, BOOL_STATUS_OK, &
                   DEC_STATUS_OK, DEC_STATUS_NO_PROGRESS, DEV_STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ey_R8P, ez_R8P

implicit none

real(R8P),    parameter  :: PI        = acos(-1._R8P)
integer(I4P), parameter  :: M         = 128_I4P                        !< Sides of the prisms.
integer(I4P), parameter  :: TARGET(2) = [1000_I4P, 600_I4P]
real(R8P),    parameter  :: ERRS(2)   = [3.e-3_R8P, 1.e-2_R8P]
type(surface_stl_object) :: part, boss, bore, dec
real(R8P),   allocatable :: rim0(:, :)                                 !< Rim vertices of the input.
real(R8P)                :: v0, dv, h
integer(I4P)             :: status, dev_status, t
logical                  :: are_tests_passed(5)

are_tests_passed = .false.
call box(part, [0._R8P, 0._R8P, 0._R8P], [2._R8P, 1._R8P, 0.5_R8P], 'base')
call cylinder(boss, [1.5_R8P, 0.5_R8P], 0.35_R8P, 0.25_R8P, 0.9_R8P, M, 'boss')
call cylinder(bore, [0.37_R8P, 0.5_R8P], 0.31_R8P, -0.1_R8P, 0.6_R8P, M, 'bore')
call part%boolean(boss, op=BOOL_UNION, status=status)
if (status /= BOOL_STATUS_OK) error stop 'union failed'
call part%boolean(bore, op=BOOL_DIFFERENCE, status=status)
if (status /= BOOL_STATUS_OK) error stop 'difference failed'
v0 = part%get_volume()
call rim_vertices(part, rim0)
print '(A,I0,A,F10.6,A,L1,A,I0)', 'bracket: ', part%get_facets_number(), ' facets, volume ', v0, ', watertight ', &
      part%is_watertight(), ', rim vertices ', size(rim0, dim=2)

do t = 1, size(TARGET)
   dec = part
   call dec%decimate(target_facets=TARGET(t), status=status)
   call report('decimate to ', TARGET(t), 0._R8P)
   are_tests_passed(t) = status == DEC_STATUS_OK .and. dec%get_facets_number() <= TARGET(t) .and. valid(dec) .and. &
                         abs(dv) <= 1.e-4_R8P .and. rims_kept(dec)
enddo
do t = 1, size(ERRS)
   dec = part
   call dec%decimate(target_facets=4_I4P, max_error=ERRS(t), status=status)
   call dec%hausdorff(part, upper=h, tolerance=ERRS(t) / 4, max_triangles=40000000_I4P, status=dev_status)
   call report('max_error ', 4_I4P, ERRS(t))
   are_tests_passed(2 + t) = status == DEC_STATUS_NO_PROGRESS .and. valid(dec) .and. dev_status == DEV_STATUS_OK .and. &
                             h <= ERRS(t) .and. rims_kept(dec)
enddo
dec = part
call dec%decimate(target_facets=600_I4P, rim_collapses=.false., status=status)
call report('rims locked, to ', 600_I4P, 0._R8P)
are_tests_passed(5) = status == DEC_STATUS_NO_PROGRESS .and. valid(dec)

print '(A,5L2)', 'per-case results: ', are_tests_passed
print '(A,L1)',  'Are all tests passed? ', all(are_tests_passed)
if (.not. all(are_tests_passed)) error stop 1

contains
   subroutine report(what, t, e)
   !< One line per case (h printed for the max_error cases).
   character(*), intent(in) :: what
   integer(I4P), intent(in) :: t
   real(R8P),    intent(in) :: e

   dv = (dec%get_volume() - v0) / v0
   if (e > 0._R8P) then
      print '(A,ES8.1,A,I0,A,I0,A,ES10.2,A,ES10.3,A,L1)', what, e, ': status ', status, ', facets ', &
            dec%get_facets_number(), ', dV/V ', dv, ', certified H <= ', h, ', valid ', valid(dec)
   else
      print '(A,I0,A,I0,A,I0,A,ES10.2,A,L1)', what, t, ': status ', status, ', facets ', dec%get_facets_number(), &
            ', dV/V ', dv, ', valid ', valid(dec)
   endif
   endsubroutine report

   function valid(s) result(yes)
   !< Watertight, manifold, one shell.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes

   yes = s%is_watertight() .and. s%is_manifold() .and. s%get_shells_number() == 1_I4P
   endfunction valid

   subroutine rim_vertices(s, v)
   !< The vertices of s shared by facets of two or more patches (coordinates, unique).
   type(surface_stl_object), intent(in)  :: s
   real(R8P), allocatable,   intent(out) :: v(:, :)
   real(R8P), allocatable                :: p(:, :)
   integer(I4P), allocatable             :: pat(:)
   logical, allocatable                  :: rim(:)
   integer(I4P)                          :: k, j, i, n

   allocate(p(3, 3 * s%get_facets_number()), pat(3 * s%get_facets_number()), rim(3 * s%get_facets_number()))
   n = 0
   do k = 1, s%get_facets_number()
      associate(f => s%facet_at(k))
         do j = 1, 3
            do i = 1, n
               if (p(1, i) == f%vertex(j)%x .and. p(2, i) == f%vertex(j)%y .and. p(3, i) == f%vertex(j)%z) exit
            enddo
            if (i > n) then
               n = n + 1 ; p(:, n) = [f%vertex(j)%x, f%vertex(j)%y, f%vertex(j)%z] ; pat(n) = f%patch
               rim(n) = .false.
            elseif (pat(i) /= f%patch) then
               rim(i) = .true.
            endif
         enddo
      endassociate
   enddo
   v = p(:, pack([(i, i = 1, n)], rim(1:n)))
   endsubroutine rim_vertices

   function rims_kept(s) result(yes)
   !< The patches of the input are all there (same table, each with facets) and every rim vertex of s is one of the input.
   type(surface_stl_object), intent(in) :: s
   logical                              :: yes
   real(R8P), allocatable               :: v(:, :)
   integer(I4P)                         :: i, j, k

   yes = s%get_patches_number() == part%get_patches_number()
   do k = 1, s%get_patches_number() ! a patch with facets in the input has facets in s
      if (count([(part%facet_patch(i) == k, i = 1, part%get_facets_number())]) == 0) cycle
      yes = yes .and. count([(s%facet_patch(i) == k, i = 1, s%get_facets_number())]) > 0
   enddo
   call rim_vertices(s, v)
   do i = 1, size(v, dim=2)
      yes = yes .and. any([(all(v(:, i) == rim0(:, j)), j = 1, size(rim0, dim=2))])
   enddo
   endfunction rims_kept

   subroutine box(s, lo, hi, name)
   !< The box [lo, hi] from cube-01 ([0, 1]^3), resized and moved, all in one patch.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: lo(3), hi(3)
   character(*),             intent(in)    :: name
   integer(I4P)                            :: k

   call s%load_from_file(file_name='src/tests/cube-01.stl', guess_format=.true.)
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
endprogram fossil_val_decimate_bracket
