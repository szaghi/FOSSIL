!< FOSSIL test: body-body distance at the edges of the per-call skips of issue #94.
!<
!< `body_distance` skips a facet's row of a leaf pair when its box is farther than the pruning bound from the union of
!< the other leaf's facet boxes, and the exact `pair_cut` when two triangle boxes are strictly apart and both triangles
!< are regular. Every case compares `distance_to` against `body_distance_brute` with == on distance, facets, witness
!< points and contact, without poses, with an exact pose (a quarter turn by a permutation matrix: contacts stay exact)
!< and with rotations:
!< - two 8 x 8 unit grids (128 facets, several leaves) 1 apart: exact ties between many pairs; also shifted by
!<   (0.5, 0.25);
!< - the grids coplanar and overlapping, and coplanar sharing an edge line (boxes touching);
!< - a tetrahedron with its apex exactly on a cube face, and 2^-40 off it;
!< - an exactly degenerate (collinear) facet above the cube, boxes apart along z only: `pair_cut` is not exact on it
!<   (it reports a contact), so the skip must not apply there.
!< Breaking any of the three skip conditions (row skip with >=, disjointness with >=, no regularity guard) fails it.
program fossil_test_body_distance_ties

use fossil,               only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle
use fossil_body_distance, only : body_distance_brute
use penf,                 only : I4P, R8P
use vecfor,               only : vector_R8P

implicit none

type(surface_stl_object) :: a, b
type(rigid_pose_t)       :: quarter, turn_a, turn_b
integer(I4P)             :: nbad, ncase, nties
logical                  :: ok

quarter%r = reshape([0._R8P, 1._R8P, 0._R8P, -1._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ! exact
quarter%t = [3._R8P, -2._R8P, 5._R8P]
turn_a = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.7_R8P, [0.1_R8P, 0.2_R8P, 0.3_R8P])
turn_b = pose_from_axis_angle([-2._R8P, 1._R8P, 0.5_R8P], 1.3_R8P, [0.3_R8P, -0.2_R8P, 1.1_R8P])
nbad = 0 ; ncase = 0 ; nties = 0
call grid(a, [0._R8P, 0._R8P, 0._R8P]) ; call grid(b, [0._R8P, 0._R8P, 1._R8P])     ! 1 apart: exact ties
call compare(a, b, 1._R8P, .false.) ; call compare_posed(a, b)
call grid(b, [0.5_R8P, 0.25_R8P, 1._R8P])                                          ! shifted ties
call compare(a, b, 1._R8P, .false.) ; call compare_posed(a, b)
call grid(b, [0._R8P, 0._R8P, 0._R8P])                                             ! coplanar, overlapping
call compare(a, b, 0._R8P, .true.) ; call compare_posed(a, b)
call grid(b, [8._R8P, 0._R8P, 0._R8P])                                             ! coplanar, an edge line shared
call compare(a, b, 0._R8P, .true.) ; call compare_posed(a, b)
call cube(a)
call tetra(b, [1._R8P, 0.5_R8P, 0.5_R8P])                                          ! apex on the face x = 1
call compare(a, b, 0._R8P, .true.) ; call compare_posed(a, b)
call tetra(b, [1._R8P + 2._R8P**(-40), 0.5_R8P, 0.5_R8P])                         ! 2^-40 off the face
call compare(a, b, 2._R8P**(-40), .false.) ; call compare_posed(a, b)
call needle(b)                                                                     ! a degenerate facet
call compare(a, b, -1._R8P, .false.) ; call compare_posed(a, b)
ok = nbad == 0 .and. ncase == 42 .and. nties == 3
print '(A,L1,A,I0,A,I0,A,I0)', 'skip edges vs brute force: ', ok, '  cases ', ncase, '  mismatches ', nbad, &
                                '  tied geometries ', nties
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine compare(s, t, dexp, hexp, pose_s, pose_t)
   !< distance_to(s, t) against body_distance_brute, bit for bit; without poses, with dexp >= 0 also the expected
   !< distance and contact, and (dexp > 0) the geometry counted if several facet pairs tie at the minimum.
   type(surface_stl_object), intent(in)           :: s, t
   real(R8P),                intent(in)           :: dexp
   logical,                  intent(in)           :: hexp
   type(rigid_pose_t),       intent(in), optional :: pose_s, pose_t
   type(facet_object), allocatable                :: fs(:), ft(:)
   real(R8P)                                      :: dd, ps(3), pt(3), e2, qs(3), qt(3)
   integer(I4P)                                   :: is_, it_, js, jt
   logical                                        :: h, m, posed

   ncase = ncase + 1 ; posed = present(pose_s) .or. present(pose_t)
   call s%distance_to(t, dd, point_self=ps, point_other=pt, facet_self=is_, facet_other=it_, intersecting=h, &
                      pose_self=pose_s, pose_other=pose_t)
   call facets_of(s, fs) ; call facets_of(t, ft)
   call body_distance_brute(fs, ft, e2, js, jt, qs, qt, m, pose_s, pose_t)
   if (dd /= sqrt(e2) .or. is_ /= js .or. it_ /= jt .or. any(ps /= qs) .or. any(pt /= qt) .or. (h .neqv. m)) then
      nbad = nbad + 1
      print '(A,I0,2ES24.16,4I6,2L2)', '   mismatch in case ', ncase, dd, sqrt(e2), is_, js, it_, jt, h, m
   endif
   if (.not. posed) then
      if (dexp >= 0._R8P .and. (dd /= dexp .or. (h .neqv. hexp))) then
         nbad = nbad + 1
         print '(A,I0,2ES24.16,2L2)', '   unexpected in case ', ncase, dd, dexp, h, hexp
      endif
      if (dexp > 0._R8P) then ! the cases apart: count the geometries where several pairs tie at the minimum
         if (ties(fs, ft, e2)) nties = nties + 1
      endif
   endif
   endsubroutine compare

   subroutine compare_posed(s, t)
   !< `compare` with the exact quarter turn on both bodies, the rotations on both, and each body alone posed.
   type(surface_stl_object), intent(in) :: s, t

   call compare(s, t, -1._R8P, .false., quarter, quarter)
   call compare(s, t, -1._R8P, .false., turn_a, turn_b)
   call compare(s, t, -1._R8P, .false., pose_s=turn_a)
   call compare(s, t, -1._R8P, .false., pose_t=turn_b)
   call compare(s, t, -1._R8P, .false., quarter, turn_b)
   endsubroutine compare_posed

   function ties(fs, ft, e2) result(yes)
   !< More than one facet pair at the minimum e2 (unposed; brute force over the pairs).
   type(facet_object), intent(in) :: fs(:), ft(:)
   real(R8P),          intent(in) :: e2
   logical                        :: yes
   real(R8P)                      :: f2, q1(3), q2(3)
   integer(I4P)                   :: i, j, n, i1, j1
   logical                        :: m

   n = 0
   do i = 1, size(fs)
      do j = 1, size(ft)
         call body_distance_brute(fs(i:i), ft(j:j), f2, i1, j1, q1, q2, m)
         if (f2 == e2) n = n + 1
      enddo
   enddo
   yes = n > 1
   endfunction ties

   subroutine grid(s, o)
   !< An 8 x 8 grid of unit squares (two facets each) in the plane z = o(3), corner (o(1), o(2)).
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: o(3)
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   integer(I4P)                            :: i, j, k

   s = fresh
   allocate(f(128)) ; k = 0
   do j = 0, 7
      do i = 0, 7
         k = k + 1
         call tri(f(k), o + [real(i, R8P), real(j, R8P), 0._R8P], o + [real(i + 1, R8P), real(j, R8P), 0._R8P], &
                  o + [real(i + 1, R8P), real(j + 1, R8P), 0._R8P])
         k = k + 1
         call tri(f(k), o + [real(i, R8P), real(j, R8P), 0._R8P], o + [real(i + 1, R8P), real(j + 1, R8P), 0._R8P], &
                  o + [real(i, R8P), real(j + 1, R8P), 0._R8P])
      enddo
   enddo
   call s%adopt_facets(f)
   endsubroutine grid

   subroutine cube(s)
   !< The unit cube [0, 1]^3.
   type(surface_stl_object), intent(inout) :: s
   type(surface_stl_object)                :: fresh

   s = fresh
   call s%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   endsubroutine cube

   subroutine tetra(s, apex)
   !< A tetrahedron with its apex at `apex` and its base in the plane x = apex(1) + 1.
   type(surface_stl_object), intent(inout) :: s
   real(R8P),                intent(in)    :: apex(3)
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   real(R8P)                               :: b1(3), b2(3), b3(3)

   s = fresh
   b1 = [apex(1) + 1._R8P, 0._R8P, 0._R8P] ; b2 = [apex(1) + 1._R8P, 1._R8P, 0._R8P]
   b3 = [apex(1) + 1._R8P, 0.5_R8P, 1._R8P]
   allocate(f(4))
   call tri(f(1), apex, b2, b1) ; call tri(f(2), apex, b3, b2) ; call tri(f(3), apex, b1, b3) ; call tri(f(4), b1, b2, b3)
   call s%adopt_facets(f)
   endsubroutine tetra

   subroutine needle(s)
   !< An exactly collinear (degenerate) facet above the unit cube (boxes apart along z, overlapping in x and y, a vertex
   !< in the plane x = 0 of two cube facets) and a triangle at z = 1.25. `pair_cut` is not exact on a degenerate facet,
   !< so its classification must not be skipped there (whatever it says, the result must stay that of brute force).
   type(surface_stl_object), intent(inout) :: s
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)

   s = fresh
   allocate(f(2))
   call tri(f(1), [0._R8P, 0._R8P, 1.5_R8P], [0.5_R8P, 0.5_R8P, 1.75_R8P], [1._R8P, 1._R8P, 2._R8P]) ! collinear
   call tri(f(2), [0._R8P, 0._R8P, 1.25_R8P], [1._R8P, 0._R8P, 1.25_R8P], [0.5_R8P, 1._R8P, 1.25_R8P])
   call s%adopt_facets(f)
   endsubroutine needle

   subroutine tri(f, p, q, r)
   !< Facet with vertices p, q, r.
   type(facet_object), intent(inout) :: f
   real(R8P),          intent(in)    :: p(3), q(3), r(3)

   f%vertex(1) = vector_R8P(p(1), p(2), p(3)) ; f%vertex(2) = vector_R8P(q(1), q(2), q(3))
   f%vertex(3) = vector_R8P(r(1), r(2), r(3))
   call f%compute_normal
   endsubroutine tri

   subroutine facets_of(s, f)
   !< A copy of the facets of s.
   type(surface_stl_object),        intent(in)  :: s
   type(facet_object), allocatable, intent(out) :: f(:)
   type(facet_object), pointer                  :: p
   integer(I4P)                                 :: i

   allocate(f(s%get_facets_number()))
   do i = 1, size(f)
      p => s%facet_at(i) ; f(i) = p
   enddo
   endsubroutine facets_of
endprogram fossil_test_body_distance_ties
