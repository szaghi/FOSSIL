!< FOSSIL, assert reflections keep an outward-oriented surface outward.
!<
!< `mirror` and `resize` with an odd number of negative factors are reflections: they
!< reverse the handedness of every facet's vertex order. Before the fix the facets kept
!< their vertex order, so a reflected outward surface came out inside-out (unit cube:
!< volume -1 after analyze, centre reported outside). Now two vertices are swapped when
!< the map reverses orientation. Cases on the sanitized unit cube, each followed by
!< `analyze`: volume +V, signed distance at the image of the centre negative, inside:
!<   mirror(normal), mirror(matrix), resize(-1,1,1), resize(-2,-2,-2) (odd, reflection),
!<   resize(-1,-1,1) (even, a rotation: must stay outward too),
!<   rotate(matrix=M) with an improper M (a reflection) and with a proper 90-degree rotation.
!< Then dragon (it replaced dragon-fine, issue #45) mirrored through x = 0 against the original at mirrored points: same
!< sign everywhere, |d| equal to a relative 1e-12.

program fossil_test_reflection_orientation

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ex_R8P, ez_R8P, mirror_matrix_R8P, rotation_matrix_R8P

implicit none

real(R8P), parameter :: REL_TOL = 1.e-12_R8P !< Relative |d| tolerance for the dragon comparison.
logical              :: ok(8)                !< Per-case results.

ok(1) = check_cube(kind=1_I4P, centre=vector_R8P(-0.5_R8P,  0.5_R8P, 0.5_R8P), volume=1._R8P,  label='mirror(normal=ex)      ')
ok(2) = check_cube(kind=2_I4P, centre=vector_R8P(-0.5_R8P,  0.5_R8P, 0.5_R8P), volume=1._R8P,  label='mirror(matrix)         ')
ok(3) = check_cube(kind=3_I4P, centre=vector_R8P(-0.5_R8P,  0.5_R8P, 0.5_R8P), volume=1._R8P,  label='resize(-1, 1, 1)       ')
ok(4) = check_cube(kind=4_I4P, centre=vector_R8P(-1._R8P,  -1._R8P, -1._R8P),  volume=8._R8P,  label='resize(-2,-2,-2)       ')
ok(5) = check_cube(kind=5_I4P, centre=vector_R8P(-0.5_R8P, -0.5_R8P, 0.5_R8P), volume=1._R8P,  label='resize(-1,-1, 1) (even)')
ok(6) = check_cube(kind=6_I4P, centre=vector_R8P(-0.5_R8P,  0.5_R8P, 0.5_R8P), volume=1._R8P,  label='rotate(matrix=improper)')
ok(7) = check_cube(kind=7_I4P, centre=vector_R8P(-0.5_R8P,  0.5_R8P, 0.5_R8P), volume=1._R8P,  label='rotate(matrix=proper)  ')
ok(8) = check_dragon()

print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1

contains
   function check_cube(kind, centre, volume, label) result(passed)
   !< Apply one reflection (or rotation) to the sanitized unit cube and check orientation.
   integer(I4P),     intent(in) :: kind    !< Transform selector.
   type(vector_R8P), intent(in) :: centre  !< Image of the cube centre.
   real(R8P),        intent(in) :: volume  !< Expected volume.
   character(*),     intent(in) :: label   !< Case label.
   logical                      :: passed  !< Check result.
   type(surface_stl_object)     :: cube    !< Unit cube.
   real(R8P)                    :: d       !< Signed distance at the centre image.
   logical                      :: inside  !< Inside test at the centre image.

   call cube%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   call cube%sanitize
   select case (kind)
   case (1) ; call cube%mirror(normal=ex_R8P)
   case (2) ; call cube%mirror(matrix=mirror_matrix_R8P(normal=ex_R8P))
   case (3) ; call cube%resize(x=-1._R8P)
   case (4) ; call cube%resize(factor=vector_R8P(-2._R8P, -2._R8P, -2._R8P))
   case (5) ; call cube%resize(factor=vector_R8P(-1._R8P, -1._R8P, 1._R8P))
   case (6) ; call cube%rotate(matrix=mirror_matrix_R8P(normal=ex_R8P))
   case (7) ; call cube%rotate(matrix=rotation_matrix_R8P(axis=ez_R8P, angle=acos(-1._R8P)/2._R8P))
   end select
   call cube%analyze
   d      = cube%distance(point=centre, is_signed=.true., is_square_root=.true.)
   inside = cube%is_point_inside(point=centre)
   passed = (abs(cube%get_volume() - volume) < 1.e-12_R8P * volume) .and. (d < 0._R8P) .and. inside
   print '(A,A,F8.4,A,F8.4,A,L1,A,L1)', label, ': volume=', cube%get_volume(), ' d(centre)=', d, &
                                        ' inside=', inside, ' ok=', passed
   endfunction check_cube

   function check_dragon() result(passed)
   !< Mirror dragon through x = 0 and compare with the original at mirrored points.
   logical                  :: passed            !< Check result.
   type(surface_stl_object) :: original, mirrored !< Surfaces.
   type(vector_R8P)         :: bmin, bmax, h, p, q !< Grid geometry, query point and its mirror image.
   real(R8P)                :: d_orig, d_mirr    !< Signed distances.
   integer(I4P), parameter  :: N = 12_I4P        !< Grid points per direction.
   integer(I4P)             :: i, j, k           !< Counters.
   integer(I4P)             :: n_sign, n_mag     !< Mismatch counters.
   type(facet_object), allocatable :: fs(:)      !< Copy of the sanitized facets.
   type(facet_object), pointer     :: fp         !< Facet.

   call original%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
   call original%sanitize
   allocate(fs(original%get_facets_number())) ! the sanitized dragon, copied (cheaper than a second load + sanitize)
   do i = 1, size(fs)
      fp => original%facet_at(i)
      fs(i) = fp
   enddo
   call mirrored%adopt_facets(facets=fs)
   call mirrored%mirror(normal=ex_R8P)
   call mirrored%analyze

   bmin = original%get_bmin()
   bmax = original%get_bmax()
   h    = bmax - bmin
   bmin = bmin - 0.1_R8P * h
   h    = 1.2_R8P * h / real(N - 1, R8P)
   n_sign = 0 ; n_mag = 0
   do k = 0, N - 1
      do j = 0, N - 1
         do i = 0, N - 1
            p%x = bmin%x + real(i, R8P) * h%x
            p%y = bmin%y + real(j, R8P) * h%y
            p%z = bmin%z + real(k, R8P) * h%z
            q   = vector_R8P(-p%x, p%y, p%z)
            d_orig = original%distance(point=p, is_signed=.true., is_square_root=.true.)
            d_mirr = mirrored%distance(point=q, is_signed=.true., is_square_root=.true.)
            if (sign(1._R8P, d_orig) /= sign(1._R8P, d_mirr)) n_sign = n_sign + 1
            if (abs(abs(d_orig) - abs(d_mirr)) > REL_TOL * max(abs(d_orig), tiny(1._R8P))) n_mag = n_mag + 1
         enddo
      enddo
   enddo
   passed = (n_sign == 0) .and. (n_mag == 0) .and. (mirrored%get_volume() > 0._R8P)
   print '(A,I0,A,I0,A,I0,A,ES12.5,A,L1)', 'dragon mirrored vs original on ', N**3, ' points: sign mismatches=', &
         n_sign, ' |d| mismatches=', n_mag, ' volume=', mirrored%get_volume(), ' ok=', passed
   endfunction check_dragon
endprogram fossil_test_reflection_orientation
