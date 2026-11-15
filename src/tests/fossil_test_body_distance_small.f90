!< FOSSIL test: body-body distance on small bodies, its cost and its bits (issue #92 WP5 d).
!<
!< `body_distance` evaluates the pair of least box gap of each leaf pair first, and once the best distance is 0 skips
!< every pair after the best in (i, j) order. Both only change the work, never the result: every case compares
!< `distance_to` against `body_distance_brute` with == on distance, facets, witness points and contact.
!< 1. Two 12-facet cubes (one leaf each), the setup of the bench line `distance_to.cubes`: 2000 random poses of the
!<    first, 53 of them in contact. The pairs evaluated by the triangle-triangle distance (`pairs=`) are bounded: on
!<    average at most 30 per call apart (27.3 measured; 37.1 before issue #92 WP5 d) and at most 5 per call in contact
!<    (3.4; 35.5 before); this part fails on the code before it.
!< 2. Several leaves: an 8 x 8 grid of unit squares (128 facets) against the cube, SAH BVH and octree, 200 random
!<    poses of the cube about the grid's centre, 141 in contact.
!< 3. Exact ties at contact: the cube and its copy touching on a face (moved by exactly 1), sharing an edge, a vertex,
!<    and overlapping (moved by 1/2), each also under an exact quarter turn of both: many pairs at distance 0, the
!<    smallest (facet_a, facet_b) wins.
!< 4. Degenerate facets (issue #96): a body with a collinear and a coincident-vertex facet next to two regular ones,
!<    200 random poses about the cube, 132 in contact.
program fossil_test_body_distance_small

use fossil,               only : surface_stl_object, facet_object, rigid_pose_t, pose_from_axis_angle, &
                                 AABB_TREE_OCTREE, AABB_TREE_SAH_BVH
use fossil_body_distance, only : body_distance_brute
use penf,                 only : I4P, I8P, R8P
use vecfor,               only : vector_R8P

implicit none

type(surface_stl_object) :: a, b
type(rigid_pose_t)       :: quarter, ps
integer(I4P)             :: nbad, ncase, ncont, k, np, nq, nk, napart, ncontact, sum_apart, sum_contact, kind_
integer(I8P)             :: seed
real(R8P)                :: w(6)
logical                  :: ok1, ok2, ok3, ok4, hit
character(*), parameter  :: CUBE_FILE = 'src/tests/cube.stl'

nbad = 0 ; ncase = 0 ; ncont = 0
! 1. two cubes, the bench setup
call a%load_from_file(file_name=CUBE_FILE, guess_format=.true.)
call b%load_from_file(file_name=CUBE_FILE, guess_format=.true.)
call b%translate(x=1.8_R8P, y=0.3_R8P, z=0.1_R8P) ; call b%analyze
seed = 20261109_I8P ; napart = 0 ; ncontact = 0 ; sum_apart = 0 ; sum_contact = 0
do k = 1, 2000
   call prng(seed, w)
   ps = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [0.5_R8P * w(5), 0.2_R8P * w(6), 0._R8P])
   call compare(a, b, np, hit, pose_s=ps)
   if (hit) then
      ncontact = ncontact + 1 ; sum_contact = sum_contact + np
   else
      napart = napart + 1 ; sum_apart = sum_apart + np
   endif
enddo
ok1 = nbad == 0 .and. ncontact > 0 .and. sum_apart <= 30 * napart .and. sum_contact <= 5 * ncontact
print '(A,L1,A,I0,A,I0,A,F6.2,A,F6.2)', 'cubes vs brute force: ', ok1, '  contacts ', ncontact, '  mismatches ', nbad, &
   '  pairs per call apart ', real(sum_apart, R8P) / max(napart, 1), '  in contact ', &
   real(sum_contact, R8P) / max(ncontact, 1)
nq = nbad
! 2. several leaves, both tree kinds
do kind_ = 1, 2
   if (kind_ == 1) then
      call grid(b, AABB_TREE_SAH_BVH)
   else
      call grid(b, AABB_TREE_OCTREE)
   endif
   seed = 20261115_I8P + kind_
   do k = 1, 100
      call prng(seed, w)
      ps = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [3._R8P + 2._R8P * w(5), 3._R8P + 2._R8P * w(6), &
                                                                   1.6_R8P * w(4) - 0.6_R8P])
      call compare(a, b, np, hit, pose_s=ps)
   enddo
enddo
ok2 = nbad == nq .and. ncont > ncontact
print '(A,L1,A,I0,A,I0)', 'grid (several leaves) vs brute force: ', ok2, '  contacts ', ncont - ncontact, &
   '  mismatches ', nbad - nq
nq = nbad ; np = ncont
! 3. exact ties at contact
quarter%r = reshape([0._R8P, 1._R8P, 0._R8P, -1._R8P, 0._R8P, 0._R8P, 0._R8P, 0._R8P, 1._R8P], [3, 3]) ! exact
quarter%t = [3._R8P, -2._R8P, 5._R8P]
call contact_case([1._R8P, 0._R8P, 0._R8P])   ! a face
call contact_case([1._R8P, 1._R8P, 0._R8P])   ! an edge
call contact_case([1._R8P, 1._R8P, 1._R8P])   ! a vertex
call contact_case([0.5_R8P, 0._R8P, 0._R8P])  ! overlapping
ok3 = nbad == nq .and. ncont - np == 8
print '(A,L1,A,I0,A,I0)', 'exact ties at contact vs brute force: ', ok3, '  contacts ', ncont - np, '  mismatches ', &
   nbad - nq
nq = nbad ; np = ncont
! 4. degenerate facets
nk = nbad
call needles(b)
seed = 20261116_I8P
do k = 1, 200
   call prng(seed, w)
   ps = pose_from_axis_angle(w(1:3) - 0.5_R8P, 6.28_R8P * w(4), [0.5_R8P + 1.2_R8P * w(5), 0.5_R8P + 1.2_R8P * w(6), &
                                                                1.2_R8P * w(4)])
   call compare(a, b, nq, hit, pose_t=ps)
enddo
ok4 = nbad == nk .and. ncont > np
print '(A,L1,A,I0,A,I0)', 'degenerate facets vs brute force: ', ok4, '  contacts ', ncont - np, '  mismatches ', nbad - nk
print '(A,L1)', 'Are all tests passed? ', ok1 .and. ok2 .and. ok3 .and. ok4

contains
   subroutine compare(s, t, pairs, h, pose_s, pose_t)
   !< distance_to(s, t) against body_distance_brute, bit for bit; the pairs distance_to evaluated and its contact.
   type(surface_stl_object), intent(in)           :: s, t
   integer(I4P),             intent(out)          :: pairs
   logical,                  intent(out)          :: h
   type(rigid_pose_t),       intent(in), optional :: pose_s, pose_t
   type(facet_object), allocatable                :: fs(:), ft(:)
   real(R8P)                                      :: dd, ps_(3), pt(3), e2, qs(3), qt(3)
   integer(I4P)                                   :: is_, it_, js, jt
   logical                                        :: m

   ncase = ncase + 1
   call s%distance_to(t, dd, point_self=ps_, point_other=pt, facet_self=is_, facet_other=it_, intersecting=h, &
                      pose_self=pose_s, pose_other=pose_t, pairs=pairs)
   call facets_of(s, fs) ; call facets_of(t, ft)
   call body_distance_brute(fs, ft, e2, js, jt, qs, qt, m, pose_s, pose_t)
   if (dd /= sqrt(e2) .or. is_ /= js .or. it_ /= jt .or. any(ps_ /= qs) .or. any(pt /= qt) .or. (h .neqv. m)) then
      nbad = nbad + 1
      print '(A,I0,2ES24.16,4I6,2L2)', '   mismatch in case ', ncase, dd, sqrt(e2), is_, js, it_, jt, h, m
   endif
   if (h) ncont = ncont + 1
   endsubroutine compare

   subroutine contact_case(shift)
   !< The cube against its copy moved by `shift` (exact), as they are and both under the exact quarter turn.
   real(R8P), intent(in) :: shift(3)
   integer(I4P)          :: n
   logical               :: h

   call b%load_from_file(file_name=CUBE_FILE, guess_format=.true.)
   call b%translate(x=shift(1), y=shift(2), z=shift(3)) ; call b%analyze
   call compare(a, b, n, h)
   call compare(a, b, n, h, quarter, quarter)
   endsubroutine contact_case

   subroutine grid(s, tree_kind)
   !< An 8 x 8 grid of unit squares (two facets each) in the plane z = 0, corner at the origin.
   type(surface_stl_object), intent(inout) :: s
   integer(I4P),             intent(in)    :: tree_kind
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)
   integer(I4P)                            :: i, j, k

   s = fresh
   allocate(f(128)) ; k = 0
   do j = 0, 7
      do i = 0, 7
         k = k + 1
         call tri(f(k), [real(i, R8P), real(j, R8P), 0._R8P], [real(i + 1, R8P), real(j, R8P), 0._R8P], &
                  [real(i + 1, R8P), real(j + 1, R8P), 0._R8P])
         k = k + 1
         call tri(f(k), [real(i, R8P), real(j, R8P), 0._R8P], [real(i + 1, R8P), real(j + 1, R8P), 0._R8P], &
                  [real(i, R8P), real(j + 1, R8P), 0._R8P])
      enddo
   enddo
   call s%adopt_facets(f, aabb_tree_kind=tree_kind)
   endsubroutine grid

   subroutine needles(s)
   !< Two regular facets and two degenerate ones (collinear; two coincident vertices) about the origin, side ~1.
   type(surface_stl_object), intent(inout) :: s
   type(surface_stl_object)                :: fresh
   type(facet_object), allocatable         :: f(:)

   s = fresh
   allocate(f(4))
   call tri(f(1), [-0.5_R8P, -0.5_R8P, 0._R8P], [0.5_R8P, -0.5_R8P, 0._R8P], [0._R8P, 0.5_R8P, 0._R8P])
   call tri(f(2), [-0.5_R8P, 0._R8P, -0.5_R8P], [0.5_R8P, 0._R8P, -0.5_R8P], [0._R8P, 0._R8P, 0.5_R8P])
   call tri(f(3), [-0.5_R8P, -0.5_R8P, -0.5_R8P], [0._R8P, 0._R8P, 0._R8P], [0.5_R8P, 0.5_R8P, 0.5_R8P]) ! collinear
   call tri(f(4), [0._R8P, -0.5_R8P, 0.25_R8P], [0._R8P, -0.5_R8P, 0.25_R8P], [0._R8P, 0.5_R8P, -0.25_R8P]) ! coincident
   call s%adopt_facets(f)
   endsubroutine needles

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

   pure subroutine prng(seed, r)
   !< Park-Miller draws in (0, 1), fixed seed (the bench's generator).
   integer(I8P), intent(inout) :: seed
   real(R8P),    intent(out)   :: r(:)
   integer(I4P)                :: i

   do i = 1, size(r)
      seed = mod(seed * 48271_I8P, 2147483647_I8P)
      r(i) = real(seed, R8P) / 2147483647._R8P
   enddo
   endsubroutine prng
endprogram fossil_test_body_distance_small
