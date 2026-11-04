!< FOSSIL test: rigid move with tree refit (issue #41 step 3).
program fossil_test_move
!< FOSSIL test: rigid move with tree refit (issue #41 step 3).
!<
!< `surface%move(pose)` (refit) against the same mesh moved and analyzed from scratch (rebuild): the
!< cube, the cone and the L-prism of the sign benchmark and dragon.stl (it replaced the bunny, issue #45),
!< under a rotation with a large translation:
!< - every facet vertex, volume, area, centroid and box bit for bit;
!< - signed distances (pseudo-normal and fast winding number), inside flags, box classification codes,
!<   line crossings and cut-cell fractions bit for bit; closest points and normals bit for bit at face
!<   feet, within 4 ulps at edge and vertex feet (the facets tied there parametrize the point
!<   differently, and which one is reported depends on the tree);
!< - the refit is faster than the rebuild (printed ratio, dragon);
!< - two moves p then q equal one move q o p to round-off;
!< - an octree surface is rebuilt by move and agrees too.
!< `curve%move(pose)` against the curve moved and rebuilt: a NACA 0012 and an annulus, distances,
!< rich queries, winding numbers, line crossings and cut cells bit for bit.

use fossil,                   only : surface_stl_object, curve_object, facet_object, cut_cells_2d_t, cut_cells_3d_t, &
                                     cut_cells_3d, rigid_pose_t, rigid_pose_2d_t, pose_from_axis_angle, pose_2d, &
                                     pose_compose, to_world, SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER, AABB_TREE_OCTREE, &
                                     FEATURE_FACE
use fossil_sign_bench_bodies, only : body_t, soup_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none

real(R8P), parameter        :: PI = acos(-1._R8P)
type(body_t), allocatable   :: b(:)
type(soup_t)                :: drg
type(surface_stl_object)    :: s
type(curve_object)          :: c
type(rigid_pose_t)          :: p, q
type(facet_object), pointer :: fp
real(R8P)                   :: u(200), v(200)
integer(I4P)                :: k, f, m
logical                     :: ok(4)

call random_seed(size=k)
call random_seed(put=[(777 + m, m = 1, k)])
call bench_bodies(b)
call s%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
allocate(drg%v(3, 3, s%get_facets_number()))
do f = 1, s%get_facets_number()
   fp => s%facet_at(f)
   do m = 1, 3
      drg%v(:, m, f) = [fp%vertex(m)%x, fp%vertex(m)%y, fp%vertex(m)%z]
   enddo
enddo
p = pose_from_axis_angle([1._R8P, -2._R8P, 0.5_R8P], 1.1_R8P, [20._R8P, 3._R8P, -7._R8P])
q = pose_from_axis_angle([0._R8P, 1._R8P, 1._R8P], -0.6_R8P, [-1._R8P, 2._R8P, 0.5_R8P])
ok(1) = same3(b(1)%soup, p, 'cube', .false.) .and. same3(b(7)%soup, p, 'cone', .false.) .and. &
        same3(b(8)%soup, p, 'L-prism', .false.) .and. same3(drg, p, 'dragon', .false.) .and. &
        same3(b(1)%soup, p, 'cube on an octree', .true.)
ok(2) = two_moves(b(8)%soup, p, q)
do k = 0, 99
   u(k + 1) = 0.5_R8P * (1._R8P - cos(PI * (100 - k) / 100)) ; v(k + 1) = thickness(u(k + 1))
enddo
u(101) = 0._R8P ; v(101) = 0._R8P
do k = 1, 99
   u(101 + k) = 0.5_R8P * (1._R8P - cos(PI * k / 100)) ; v(101 + k) = -thickness(u(101 + k))
enddo
ok(3) = same2(u, v, [200_I4P], pose_2d(0.9_R8P, [5._R8P, -2._R8P]), 'NACA 0012')
ok(4) = same2([0._R8P, 6._R8P, 6._R8P, 0._R8P, 1._R8P, 1._R8P, 5._R8P, 5._R8P], &
              [0._R8P, 0._R8P, 6._R8P, 6._R8P, 1._R8P, 5._R8P, 5._R8P, 1._R8P], [4_I4P, 4_I4P], &
              pose_2d(-0.4_R8P, [30._R8P, 7._R8P]), 'annulus')

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function same3(sp, pose, label, octree) result(same)
   !< move (refit) vs moved-and-rebuilt mesh: geometry and queries bit for bit.
   type(soup_t),       intent(in) :: sp     !< Body-frame mesh.
   type(rigid_pose_t), intent(in) :: pose   !< Pose.
   character(*),       intent(in) :: label  !< Label.
   logical,            intent(in) :: octree !< Use an octree.
   logical                        :: same   !< Result.
   integer, parameter       :: N = 300 !< Points.
   type(soup_t)             :: mv      !< Moved soup.
   type(surface_stl_object) :: sm, sr  !< Moved by move, rebuilt.
   type(vector_R8P)         :: qp(N)   !< Points.
   type(cut_cells_3d_t)     :: cm, cr  !< Cut cells.
   type(facet_object), pointer :: f1, f2 !< Facets.
   real(R8P)    :: lo(3), hi(3), r(3), w(3), d1(N), d2(N), c1(3, N), c2(3, N), n1(3, N), n2(3, N), bl(3, N), bh(3, N), e(7) !< Work.
   real(R8P),    allocatable :: t1(:), t2(:) !< Crossings.
   integer(I4P), allocatable :: g1(:), g2(:) !< Crossing facets.
   logical,      allocatable :: e1(:), e2(:) !< Crossing flags.
   integer(I4P) :: i, ff, mm, nb, k1(N), k2(N), n_1, n_2, ax, ft(N) !< Work.
   integer(I8P) :: c0, cx, rate !< Clock.
   real(R8P)    :: t_refit, t_build !< Timings.
   logical      :: amb(N) !< Ambiguous feet.
   type(vector_R8P) :: vb1, vb2 !< Boxes.

   mv = sp
   do ff = 1, size(sp%v, 3)
      do mm = 1, 3
         mv%v(:, mm, ff) = to_world(pose, sp%v(:, mm, ff))
      enddo
   enddo
   call to_surface(sp, sm)
   if (octree) call sm%analyze(aabb_tree_kind=AABB_TREE_OCTREE)
   call system_clock(c0, rate) ; call sm%move(pose) ; call system_clock(cx) ; t_refit = real(cx - c0, R8P) / rate
   call system_clock(c0) ; call to_surface(mv, sr) ; call system_clock(cx) ; t_build = real(cx - c0, R8P) / rate
   if (octree) call sr%analyze(aabb_tree_kind=AABB_TREE_OCTREE)
   nb = 0
   do ff = 1, sm%get_facets_number()
      f1 => sm%facet_at(ff) ; f2 => sr%facet_at(ff)
      do mm = 1, 3
         if (f1%vertex(mm)%x /= f2%vertex(mm)%x .or. f1%vertex(mm)%y /= f2%vertex(mm)%y .or. &
             f1%vertex(mm)%z /= f2%vertex(mm)%z) nb = nb + 1
      enddo
   enddo
   vb1 = sm%get_bmin() ; vb2 = sr%get_bmin()
   if (vb1%x /= vb2%x .or. vb1%y /= vb2%y .or. vb1%z /= vb2%z) nb = nb + 1
   if (sm%get_volume() /= sr%get_volume() .or. sm%get_area() /= sr%get_area()) nb = nb + 1
   vb1 = sm%get_centroid() ; vb2 = sr%get_centroid()
   if (vb1%x /= vb2%x .or. vb1%y /= vb2%y .or. vb1%z /= vb2%z) nb = nb + 1
   lo = [minval(mv%v(1, :, :)), minval(mv%v(2, :, :)), minval(mv%v(3, :, :))]
   hi = [maxval(mv%v(1, :, :)), maxval(mv%v(2, :, :)), maxval(mv%v(3, :, :))]
   do i = 1, N
      call random_number(r) ; w = lo - 0.2_R8P * (hi - lo) + 1.4_R8P * r * (hi - lo) ; qp(i) = vector_R8P(w(1), w(2), w(3))
      call random_number(r) ; bl(:, i) = w ; bh(:, i) = w + 0.1_R8P * r * (hi - lo)
   enddo
   call sm%distance_many(points=qp, distances=d1, is_signed=.true.)
   call sr%distance_many(points=qp, distances=d2, is_signed=.true.)
   nb = nb + count(d1 /= d2)
   call sm%distance_many(points=qp, distances=d1, is_signed=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   call sr%distance_many(points=qp, distances=d2, is_signed=.true., sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   nb = nb + count(d1 /= d2)
   call sm%query_many(points=qp, distances=d1, closest_points=c1, normals=n1)
   call sr%query_many(points=qp, distances=d2, closest_points=c2, normals=n2, features=ft, ambiguous=amb, &
                      rel_tolerance=1.e-12_R8P)
   do i = 1, N
      if (d1(i) /= d2(i)) nb = nb + 1
      if (amb(i)) cycle
      if (ft(i) == FEATURE_FACE) then ! a unique facet: bit for bit
         if (any(c1(:, i) /= c2(:, i)) .or. any(n1(:, i) /= n2(:, i))) nb = nb + 1
      else ! an edge or vertex foot: the tied facets parametrize it differently, to the last bits
         if (any(abs(c1(:, i) - c2(:, i)) > 4._R8P * spacing(maxval(abs(c2(:, i))))) .or. &
             any(abs(n1(:, i) - n2(:, i)) > 4._R8P * epsilon(1._R8P))) nb = nb + 1
      endif
   enddo
   call sm%box_classify_many(lo=bl, hi=bh, codes=k1, radius=0.05_R8P) ; call sr%box_classify_many(lo=bl, hi=bh, codes=k2, &
                                                                                                  radius=0.05_R8P)
   nb = nb + count(k1 /= k2)
   do ax = 1, 3
      do i = 1, 20
         call sm%line_crossings(axis=int(ax, I4P), a=qp(i)%y, b=qp(i)%z, t0=-1.e3_R8P, t1=1.e3_R8P, t=t1, facets=g1, &
                                entering=e1, n=n_1)
         call sr%line_crossings(axis=int(ax, I4P), a=qp(i)%y, b=qp(i)%z, t0=-1.e3_R8P, t1=1.e3_R8P, t=t2, facets=g2, &
                                entering=e2, n=n_2)
         if (n_1 /= n_2) then
            nb = nb + 1
         elseif (n_1 > 0) then
            if (any(t1 /= t2) .or. any(g1 /= g2) .or. any(e1 .neqv. e2)) nb = nb + 1
         endif
      enddo
   enddo
   call cut_cells_3d(surface=sm, x=lin(lo(1), hi(1)), y=lin(lo(2), hi(2)), z=lin(lo(3), hi(3)), cut=cm)
   call cut_cells_3d(surface=sr, x=lin(lo(1), hi(1)), y=lin(lo(2), hi(2)), z=lin(lo(3), hi(3)), cut=cr)
   nb = nb + count(cm%alpha /= cr%alpha)
   same = nb == 0
   e = 0._R8P
   print '(A,A,I0,A,F7.3,A,L1)', label, ': mismatches vs rebuilt ', nb, '  refit/rebuild time ', &
         t_refit / max(t_build, 1.e-9_R8P), '  ok ', same
   endfunction same3

   function lin(a, bb) result(e)
   !< 6 cells over [a, b] widened by 10%.
   real(R8P), intent(in) :: a, bb !< Range.
   real(R8P)             :: e(7)  !< Edges.
   integer               :: k     !< Counter.

   do k = 1, 7
      e(k) = a - 0.1_R8P * (bb - a) + 1.2_R8P * (bb - a) * (k - 1) / 6._R8P
   enddo
   endfunction lin

   function two_moves(sp, p1, p2) result(good)
   !< move(p1) then move(p2) = move(p2 o p1) to round-off.
   type(soup_t),       intent(in) :: sp     !< Mesh.
   type(rigid_pose_t), intent(in) :: p1, p2 !< Poses.
   logical                        :: good   !< Result.
   type(surface_stl_object)       :: sa, sb !< Surfaces.
   type(facet_object), pointer    :: fa, fb !< Facets.
   real(R8P)                      :: e      !< Max difference.
   integer(I4P)                   :: ff, mm !< Counters.

   call to_surface(sp, sa) ; call to_surface(sp, sb)
   call sa%move(p1) ; call sa%move(p2)
   call sb%move(pose_compose(p2, p1))
   e = 0._R8P
   do ff = 1, sa%get_facets_number()
      fa => sa%facet_at(ff) ; fb => sb%facet_at(ff)
      do mm = 1, 3
         e = max(e, abs(fa%vertex(mm)%x - fb%vertex(mm)%x), abs(fa%vertex(mm)%y - fb%vertex(mm)%y), &
                 abs(fa%vertex(mm)%z - fb%vertex(mm)%z))
      enddo
   enddo
   good = e < 1.e-13_R8P * 25._R8P
   print '(A,ES9.1,A,L1)', 'two moves vs one composed move: max vertex difference ', e, '  ok ', good
   endfunction two_moves

   function same2(uu, vv, sizes, pose, label) result(same)
   !< curve move vs moved-and-rebuilt curve: queries bit for bit.
   real(R8P),             intent(in) :: uu(:), vv(:) !< Vertices.
   integer(I4P),          intent(in) :: sizes(:)     !< Loop sizes.
   type(rigid_pose_2d_t), intent(in) :: pose         !< Pose.
   character(*),          intent(in) :: label        !< Label.
   logical                           :: same         !< Result.
   integer, parameter   :: N = 300 !< Points.
   type(curve_object)   :: cm, cr  !< Moved, rebuilt.
   type(cut_cells_2d_t) :: km, kr  !< Cut cells.
   real(R8P) :: mu(size(uu)), mvv(size(uu)), w(2), lo(2), hi(2), r(2), qq(2, N), d1(N), d2(N), c1(2, N), c2(2, N), e(9) !< Work.
   real(R8P),    allocatable :: t1(:), t2(:) !< Crossings.
   integer(I4P), allocatable :: g1(:), g2(:) !< Segments.
   logical,      allocatable :: e1(:), e2(:) !< Flags.
   logical      :: amb(N) !< Ambiguous.
   integer(I4P) :: i, nb, n_1, n_2, k !< Work.

   do i = 1, size(uu)
      w = to_world(pose, [uu(i), vv(i)]) ; mu(i) = w(1) ; mvv(i) = w(2)
   enddo
   call cm%adopt_loops(u=uu, v=vv, loop_sizes=sizes) ; call cm%move(pose)
   call cr%adopt_loops(u=mu, v=mvv, loop_sizes=sizes)
   lo = cr%get_bmin() ; hi = cr%get_bmax()
   nb = 0
   if (any(cm%get_bmin() /= lo) .or. any(cm%get_bmax() /= hi) .or. cm%get_area() /= cr%get_area()) nb = nb + 1
   do i = 1, N
      call random_number(r) ; qq(:, i) = lo - 0.2_R8P * (hi - lo) + 1.4_R8P * r * (hi - lo)
   enddo
   call cm%distance_many(points=qq, distances=d1) ; call cr%distance_many(points=qq, distances=d2)
   nb = nb + count(d1 /= d2)
   call cm%distance_many(points=qq, distances=d1, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   call cr%distance_many(points=qq, distances=d2, sign_algorithm=SIGN_FAST_WINDING_NUMBER)
   nb = nb + count(d1 /= d2)
   call cm%query_many(points=qq, distances=d1, closest_points=c1)
   call cr%query_many(points=qq, distances=d2, closest_points=c2, ambiguous=amb, rel_tolerance=1.e-12_R8P)
   do i = 1, N
      if (d1(i) /= d2(i)) nb = nb + 1
      if (.not. amb(i) .and. any(c1(:, i) /= c2(:, i))) nb = nb + 1
   enddo
   do i = 1, 30
      call cm%line_crossings(axis=1_I4P, a=qq(2, i), t0=-1.e3_R8P, t1=1.e3_R8P, t=t1, segments=g1, entering=e1, n=n_1)
      call cr%line_crossings(axis=1_I4P, a=qq(2, i), t0=-1.e3_R8P, t1=1.e3_R8P, t=t2, segments=g2, entering=e2, n=n_2)
      if (n_1 /= n_2) then
         nb = nb + 1
      elseif (n_1 > 0) then
         if (any(t1 /= t2) .or. any(g1 /= g2) .or. any(e1 .neqv. e2)) nb = nb + 1
      endif
   enddo
   e = [(lo(1) - 0.1_R8P * (hi(1) - lo(1)) + 1.2_R8P * (hi(1) - lo(1)) * (k - 1) / 8._R8P, k = 1, 9)]
   call cm%cut_cells(x=e, y=[(lo(2) - 0.1_R8P * (hi(2) - lo(2)) + 1.2_R8P * (hi(2) - lo(2)) * (k - 1) / 8._R8P, k = 1, 9)], cut=km)
   call cr%cut_cells(x=e, y=[(lo(2) - 0.1_R8P * (hi(2) - lo(2)) + 1.2_R8P * (hi(2) - lo(2)) * (k - 1) / 8._R8P, k = 1, 9)], cut=kr)
   nb = nb + count(km%alpha /= kr%alpha)
   same = nb == 0
   print '(A,A,I0,A,L1)', label, ': mismatches vs rebuilt ', nb, '  ok ', same
   endfunction same2

   pure function thickness(xx) result(y)
   !< NACA 0012 half thickness (closed trailing edge).
   real(R8P), intent(in) :: xx !< Chord position.
   real(R8P)             :: y  !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(xx) - 0.126_R8P * xx - 0.3516_R8P * xx**2 + 0.2843_R8P * xx**3 - 0.1036_R8P * xx**4)
   endfunction thickness
endprogram fossil_test_move
