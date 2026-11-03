!< FOSSIL test: batched rich queries on surfaces (issue #37 step 2).
program fossil_test_query_many
!< FOSSIL test: batched rich queries on surfaces (issue #37 step 2).
!<
!< - Unit cube (cube-02, grid vertices every 0.25), analytic: a facet-interior point (distance 1, closest point on the face, normal
!<   +z, FEATURE_FACE), an edge point (FEATURE_EDGE, normal (1,1,0)/sqrt2), a vertex point
!<   (FEATURE_VERTEX, normal (1,1,1)/sqrt3); signed negative inside.
!< - cube-02 and bunny, a 12^3 grid around each: every output of query_many equals the scalar
!<   compute_distance (signed, pseudo-normal, square root; facet, edge/vertex index, closest point,
!<   normal) bit for bit, and bodies equal facet_body of the facet.
!< - Fast-winding-number sign: query_many equals compute_distance(sign_algorithm=FWN) bit for bit.
!< - A size mismatch is refused with STATUS_INVALID_INPUT.

use fossil, only : surface_stl_object, FEATURE_FACE, FEATURE_EDGE, FEATURE_VERTEX, SIGN_PSEUDO_NORMAL, &
                   SIGN_FAST_WINDING_NUMBER, STATUS_INVALID_INPUT
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object) :: cube, bunny
type(vector_R8P)         :: p(3)
real(R8P)                :: d(3), cp(3, 3), nrm(3, 3)
integer(I4P)             :: ft(3), fi(3), st
logical                  :: ok(5)

call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
p(1) = vector_R8P(0.55_R8P, 0.4_R8P, 2._R8P) ! inside a facet (cube-02 has grid vertices every 0.25)
p(2) = vector_R8P(2._R8P, 2._R8P, 0.6_R8P)  ! mid-edge, between grid vertices
p(3) = vector_R8P(2._R8P, 2._R8P, 2._R8P)
call cube%query_many(points=p, distances=d, closest_points=cp, features=ft, feature_indexes=fi, normals=nrm)
ok(1) = abs(d(1) - 1._R8P) < 1.e-15_R8P .and. all(abs(cp(:, 1) - [0.55_R8P, 0.4_R8P, 1._R8P]) < 1.e-15_R8P) .and. ft(1) == FEATURE_FACE .and. &
        all(abs(nrm(:, 1) - [0._R8P, 0._R8P, 1._R8P]) < 1.e-15_R8P) .and. &
        abs(d(2) - sqrt(2._R8P)) < 1.e-15_R8P .and. all(abs(cp(:, 2) - [1._R8P, 1._R8P, 0.6_R8P]) < 1.e-15_R8P) .and. ft(2) == FEATURE_EDGE .and. &
        all(abs(nrm(:, 2) - [1._R8P, 1._R8P, 0._R8P] / sqrt(2._R8P)) < 1.e-12_R8P) .and. &
        abs(d(3) - sqrt(3._R8P)) < 1.e-15_R8P .and. all(cp(:, 3) == [1._R8P, 1._R8P, 1._R8P]) .and. ft(3) == FEATURE_VERTEX .and. &
        all(abs(nrm(:, 3) - [1._R8P, 1._R8P, 1._R8P] / sqrt(3._R8P)) < 1.e-12_R8P) .and. all(fi(2:3) >= 1) .and. fi(1) == 0
call cube%query_many(points=[vector_R8P(0.5_R8P, 0.5_R8P, 0.25_R8P)], distances=d(1:1))
ok(1) = ok(1) .and. abs(d(1) + 0.25_R8P) < 1.e-15_R8P
print '(A,L1)', 'cube, face/edge/vertex analytic: ', ok(1)

ok(2) = compare(cube, 'cube-02', SIGN_PSEUDO_NORMAL)
call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
ok(3) = compare(bunny, 'bunny', SIGN_PSEUDO_NORMAL)
ok(4) = compare(bunny, 'bunny, FWN sign', SIGN_FAST_WINDING_NUMBER)

call cube%query_many(points=p, distances=d, closest_points=cp(:, 1:2), status=st)
ok(5) = st == STATUS_INVALID_INPUT
print '(A,L1)', 'size mismatch refused: ', ok(5)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
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
