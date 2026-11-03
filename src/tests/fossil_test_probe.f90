!< FOSSIL test: probe and image points with visibility (issue #37 step 6).
program fossil_test_probe
!< FOSSIL test: probe and image points with visibility (issue #37 step 6).
!<
!< 3-D: cube A = [0,1]^3 (cube-02) and cube B = A shifted by 1.2 along x (a gap of 0.2), merged.
!< - From (0.9, 0.55, 0.4) (foot (1, 0.55, 0.4), normal +x): probes of 0.1 and -0.3 are visible, a
!<   probe of 0.5 crosses B and is not; the image point (1.1, 0.55, 0.4) is visible; the image of
!<   (0.7, 0.55, 0.4) is (1.3, ...), inside B, not visible.
!< - Placement rules at dx = 0.05 from the same point: sqrt(2) dx (Wang), sqrt(3) dx (image line in
!<   3-D, Funada & Imamura), fixed 0.08 (van Noordt) visible; a second layer at
!<   0.1 + 2 sqrt(3) dx = 0.273 (Constant et al.) reaches past x = 1.2, not visible.
!< 2-D: the same with squares [0,1]^2 and [1.2,2.2] x [0,1], exact visibility.
!< - A 256-gon of radius 0.5: the image of an inside point equals 2 foot - x.

use fossil, only : surface_stl_object, facet_object, curve_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter            :: PI = acos(-1._R8P), DX = 0.05_R8P
type(surface_stl_object)        :: cube, two
type(curve_object)              :: c
type(facet_object), allocatable :: f(:), g(:)
type(facet_object), pointer     :: fp
real(R8P)                       :: pp(3, 7), img(3, 2), ft(3, 7), pp2(2, 7), img2(2, 2), u(256), v(256), fo(2, 1)
logical                         :: vis(7), vis2(7), ivis(2), ok(4)
real(R8P)                       :: lengths(7)
integer(I4P)                    :: i, n, k

lengths = [0.1_R8P, -0.3_R8P, 0.5_R8P, sqrt(2._R8P) * DX, sqrt(3._R8P) * DX, 0.08_R8P, 0.1_R8P + 2._R8P * sqrt(3._R8P) * DX]
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
n = cube%get_facets_number()
allocate(f(n), g(n))
do i = 1, n
   fp => cube%facet_at(i)
   f(i) = fp ; f(i)%vertex_id = 0_I4P
   g(i) = fp ; g(i)%vertex_id = 0_I4P
   do k = 1, 3
      g(i)%vertex(k)%x = g(i)%vertex(k)%x + 1.2_R8P
   enddo
   call g(i)%compute_metrix
enddo
f = [f, g]
call two%adopt_facets(facets=f)
call two%probe_many(points=[(vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P), i = 1, 7)], lengths=lengths, probe_points=pp, &
                    visible=vis, foot_points=ft)
ok(1) = all(abs(ft(:, 1) - [1._R8P, 0.55_R8P, 0.4_R8P]) < 1.e-15_R8P) .and. &
        all(abs(pp(:, 1) - [1.1_R8P, 0.55_R8P, 0.4_R8P]) < 1.e-15_R8P) .and. &
        all(vis .eqv. [.true., .true., .false., .true., .true., .true., .false.])
call two%image_points(points=[vector_R8P(0.9_R8P, 0.55_R8P, 0.4_R8P), vector_R8P(0.7_R8P, 0.55_R8P, 0.4_R8P)], images=img, &
                      visible=ivis)
ok(1) = ok(1) .and. all(abs(img(:, 1) - [1.1_R8P, 0.55_R8P, 0.4_R8P]) < 1.e-15_R8P) .and. &
        all(abs(img(:, 2) - [1.3_R8P, 0.55_R8P, 0.4_R8P]) < 1.e-15_R8P) .and. ivis(1) .and. .not. ivis(2)
print '(A,L1,A,7L2)', '3-D two cubes, probes and images: ', ok(1), '  visible:', vis

call c%adopt_loops(u=[0._R8P, 1._R8P, 1._R8P, 0._R8P,  1.2_R8P, 2.2_R8P, 2.2_R8P, 1.2_R8P], &
                   v=[0._R8P, 0._R8P, 1._R8P, 1._R8P,  0._R8P, 0._R8P, 1._R8P, 1._R8P], loop_sizes=[4_I4P, 4_I4P])
call c%probe_many(points=reshape([(0.9_R8P, 0.4_R8P, i = 1, 7)], [2, 7]), lengths=lengths, probe_points=pp2, visible=vis2)
ok(2) = all(vis2 .eqv. [.true., .true., .false., .true., .true., .true., .false.]) .and. &
        all(abs(pp2(:, 1) - [1.1_R8P, 0.4_R8P]) < 1.e-15_R8P)
call c%image_points(points=reshape([0.9_R8P, 0.4_R8P, 0.7_R8P, 0.4_R8P], [2, 2]), images=img2, visible=ivis)
ok(2) = ok(2) .and. all(abs(img2(:, 2) - [1.3_R8P, 0.4_R8P]) < 1.e-15_R8P) .and. ivis(1) .and. .not. ivis(2)
print '(A,L1,A,7L2)', '2-D two squares, probes and images: ', ok(2), '  visible:', vis2

do k = 1, 256
   u(k) = 0.5_R8P * cos(2._R8P * PI * (k - 1) / 256) ; v(k) = 0.5_R8P * sin(2._R8P * PI * (k - 1) / 256)
enddo
call c%adopt_loops(u=u, v=v, loop_sizes=[256_I4P])
call c%image_points(points=reshape([0.1_R8P, 0.2_R8P], [2, 1]), images=img2(:, 1:1), visible=ivis(1:1), foot_points=fo)
ok(3) = all(abs(img2(:, 1) - (2._R8P * fo(:, 1) - [0.1_R8P, 0.2_R8P])) < 1.e-15_R8P) .and. ivis(1) .and. &
        abs(hypot(fo(1, 1), fo(2, 1)) - 0.5_R8P) < 1.e-4_R8P
print '(A,L1)', '2-D 256-gon image: ', ok(3)

call c%probe_many(points=reshape([0.1_R8P, 0.2_R8P], [2, 1]), lengths=[0.1_R8P, 0.2_R8P], probe_points=pp2(:, 1:1), &
                  visible=vis2(1:1), status=k)
ok(4) = k /= 0
print '(A,L1)', 'size mismatch refused: ', ok(4)

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_test_probe
