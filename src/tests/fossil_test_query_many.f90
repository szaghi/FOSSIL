!< FOSSIL test: batched rich queries on surfaces (issue #37 step 2).
program fossil_test_query_many
!< FOSSIL test: batched rich queries on surfaces (issue #37 step 2).
!<
!< - Unit cube (cube-02, grid vertices every 0.25), analytic: a facet-interior point (distance 1, closest point on the face, normal
!<   +z, FEATURE_FACE), an edge point (FEATURE_EDGE, normal (1,1,0)/sqrt2), a vertex point
!<   (FEATURE_VERTEX, normal (1,1,1)/sqrt3); signed negative inside.
!< - cube-02 and dragon.stl (it replaced the bunny, issue #45), a 12^3 grid around each: every
!<   output of query_many equals the scalar compute_distance (signed, pseudo-normal, square root;
!<   facet, edge/vertex index, closest point, normal) bit for bit, and bodies equal facet_body of
!<   the facet.
!< - Fast-winding-number sign: query_many equals compute_distance(sign_algorithm=FWN) bit for bit.
!< - A size mismatch is refused with STATUS_INVALID_INPUT.
!< - query_many_device (issue #37 step 3) equals query_many bit for bit, every output, both signs,
!<   on cube-02 and dragon. Built without OpenACC (tests-gnu) the device loop runs on the host: this
!<   checks the device path's logic; the GPU check is the manual bench (pending). Before
!<   enter_device it is refused.

use fossil, only : surface_stl_object, FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX, SIGN_PSEUDO_NORMAL, &
                   SIGN_FAST_WINDING_NUMBER, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, dragon
type(vector_R8P)         :: p(3)
real(R8P)                :: d(3), cp(3, 3), nrm(3, 3)
integer(I4P)             :: ft(3), fi(3), st
logical                  :: ok(7)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
p(1) = vector_R8P(0.55_R8P, 0.4_R8P, 2._R8P) ! inside a facet (cube-02 has grid vertices every 0.25)
p(2) = vector_R8P(2._R8P, 2._R8P, 0.6_R8P)  ! mid-edge, between grid vertices
p(3) = vector_R8P(2._R8P, 2._R8P, 2._R8P)
call cube%query_many(points=p, distances=d, closest_points=cp, features=ft, feature_indexes=fi, normals=nrm)
ok(1) = abs(d(1) - 1._R8P) < 1.e-15_R8P .and. all(abs(cp(:, 1) - [0.55_R8P, 0.4_R8P, 1._R8P]) < 1.e-15_R8P) .and. &
        ft(1) == FEATURE_FACE .and. &
        all(abs(nrm(:, 1) - [0._R8P, 0._R8P, 1._R8P]) < 1.e-15_R8P) .and. &
        abs(d(2) - sqrt(2._R8P)) < 1.e-15_R8P .and. all(abs(cp(:, 2) - [1._R8P, 1._R8P, 0.6_R8P]) < 1.e-15_R8P) .and. &
        ft(2) == FEATURE_EDGE .and. &
        all(abs(nrm(:, 2) - [1._R8P, 1._R8P, 0._R8P] / sqrt(2._R8P)) < 1.e-12_R8P) .and. &
        abs(d(3) - sqrt(3._R8P)) < 1.e-15_R8P .and. all(cp(:, 3) == [1._R8P, 1._R8P, 1._R8P]) .and. ft(3) == FEATURE_VERTEX .and. &
        all(abs(nrm(:, 3) - [1._R8P, 1._R8P, 1._R8P] / sqrt(3._R8P)) < 1.e-12_R8P) .and. all(fi(2:3) >= 1) .and. fi(1) == 0
call cube%query_many(points=[vector_R8P(0.5_R8P, 0.5_R8P, 0.25_R8P)], distances=d(1:1))
ok(1) = ok(1) .and. abs(d(1) + 0.25_R8P) < 1.e-15_R8P
print '(A,L1)', 'cube, face/edge/vertex analytic: ', ok(1)

ok(2) = compare(cube, 'cube-02', SIGN_PSEUDO_NORMAL)
call dragon%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
ok(3) = compare(dragon, 'dragon', SIGN_PSEUDO_NORMAL)
ok(4) = compare(dragon, 'dragon, FWN sign', SIGN_FAST_WINDING_NUMBER)

call cube%query_many(points=p, distances=d, closest_points=cp(:, 1:2), status=st)
ok(5) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'size mismatch refused: ', ok(5)

call cube%query_many_device(points=p, distances=d, status=st)
ok(6) = st == STATUS_INVALID_INPUT
ok(6) = ok(6) .and. device_path(cube, 'cube-02') .and. device_path(dragon, 'dragon')
call dragon%query_many_device(points=p, distances=d, sign_algorithm=3_I4P - 1_I4P, status=st)
ok(7) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'device path refusals (not entered, unsupported sign): ', ok(6) .and. ok(7)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function device_path(s, label) result(good)
   !< query_many_device == query_many on a 12^3 grid, every output, both signs.
   type(surface_stl_object), intent(inout) :: s     !< Surface.
   character(*),             intent(in)    :: label !< Label.
   logical                                 :: good  !< Result.
   integer(I4P), parameter                 :: M = 12, N = M**3 !< Grid.
   type(vector_R8P)                        :: q(N), lo, hi !< Points, box.
   real(R8P)    :: d1(N), d2(N), c1(3, N), c2(3, N), n1(3, N), n2(3, N) !< Outputs.
   integer(I4P) :: f1(N), f2(N), t1(N), t2(N), i1(N), i2(N), b1(N), b2(N), i, j, k, l, a, mism !< Outputs, work.
   integer(I4P), parameter :: ALGO(2) = [SIGN_PSEUDO_NORMAL, SIGN_FAST_WINDING_NUMBER] !< Signs.

   lo = s%get_bmin() ; hi = s%get_bmax()
   l = 0
   do k = 1, M
      do j = 1, M
         do i = 1, M
            l = l + 1
            q(l)%x = lo%x + (hi%x - lo%x) * (1.3_R8P * (i - 1) / (M - 1) - 0.15_R8P)
            q(l)%y = lo%y + (hi%y - lo%y) * (1.3_R8P * (j - 1) / (M - 1) - 0.15_R8P)
            q(l)%z = lo%z + (hi%z - lo%z) * (1.3_R8P * (k - 1) / (M - 1) - 0.15_R8P)
         enddo
      enddo
   enddo
   call s%enter_device
   mism = 0
   do a = 1, 2
      call s%query_many(points=q, distances=d1, closest_points=c1, facets=f1, features=t1, feature_indexes=i1, normals=n1, &
                        bodies=b1, sign_algorithm=ALGO(a))
      call s%query_many_device(points=q, distances=d2, closest_points=c2, facets=f2, features=t2, feature_indexes=i2, &
                               normals=n2, bodies=b2, sign_algorithm=ALGO(a))
      do l = 1, N
         if (d1(l) /= d2(l) .or. any(c1(:, l) /= c2(:, l)) .or. any(n1(:, l) /= n2(:, l)) .or. f1(l) /= f2(l) .or. &
             t1(l) /= t2(l) .or. i1(l) /= i2(l) .or. b1(l) /= b2(l)) mism = mism + 1
      enddo
   enddo
   call s%exit_device
   good = mism == 0
   print '(A,A,I0,A,I0)', label, ': device path vs query_many mismatches ', mism, ' of ', 2 * N
   endfunction device_path

   function compare(s, label, algo) result(good)
   !< query_many against compute_distance on a 12^3 grid around the surface, every output.
   type(surface_stl_object), intent(in) :: s     !< Surface.
   character(*),             intent(in) :: label !< Label.
   integer(I4P),             intent(in) :: algo  !< Sign algorithm.
   logical                              :: good  !< Result.
   integer(I4P), parameter              :: M = 12, N = M**3 !< Grid.
   type(vector_R8P)                     :: q(N), lo, hi, c, nn !< Points, box, scalar outputs.
   real(R8P)                            :: dq(N), cq(3, N), nq(3, N), ds !< Batch outputs.
   integer(I4P)                         :: fq(N), ftq(N), fiq(N), bq(N), fs, es, vs, i, j, k, l, mism !< Work.

   lo = s%get_bmin() ; hi = s%get_bmax()
   l = 0
   do k = 1, M
      do j = 1, M
         do i = 1, M
            l = l + 1
            q(l)%x = lo%x + (hi%x - lo%x) * (1.3_R8P * (i - 1) / (M - 1) - 0.15_R8P)
            q(l)%y = lo%y + (hi%y - lo%y) * (1.3_R8P * (j - 1) / (M - 1) - 0.15_R8P)
            q(l)%z = lo%z + (hi%z - lo%z) * (1.3_R8P * (k - 1) / (M - 1) - 0.15_R8P)
         enddo
      enddo
   enddo
   call s%query_many(points=q, distances=dq, closest_points=cq, facets=fq, features=ftq, feature_indexes=fiq, normals=nq, &
                     bodies=bq, sign_algorithm=algo)
   mism = 0
   do l = 1, N
      call s%compute_distance(point=q(l), distance=ds, is_signed=.true., sign_algorithm=algo, is_square_root=.true., &
                              facet_index=fs, edge_index=es, vertex_index=vs, closest_point=c, normal=nn)
      if (dq(l) /= ds .or. fq(l) /= fs .or. any(cq(:, l) /= [c%x, c%y, c%z]) .or. any(nq(:, l) /= [nn%x, nn%y, nn%z]) .or. &
          bq(l) /= s%facet_body(fs) .or. fiq(l) /= max(es, vs)) mism = mism + 1
      if (es > 0 .and. ftq(l) /= FEATURE_EDGE) mism = mism + 1
      if (vs > 0 .and. ftq(l) /= FEATURE_VERTEX) mism = mism + 1
      if (es == 0 .and. vs == 0 .and. ftq(l) /= FEATURE_FACE) mism = mism + 1
   enddo
   good = mism == 0
   print '(A,A,I0,A,I0)', label, ': query_many vs compute_distance mismatches ', mism, ' of ', N
   endfunction compare
endprogram fossil_test_query_many
