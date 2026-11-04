!< FOSSIL, assert sanitize leaves signed-distance queries correct on wrongly oriented input.
!<
!< `sanitize` fixes the winding in `sanitize_normals` (BFS flips + a global flip when the
!< volume is negative). The pseudo-normals, `pnormal_payload` and the AABB distance
!< payload are built by `analyze` from the winding at that time, so after a flip they
!< describe the old orientation. Before the fix, an inside-out cube put through
!< `sanitize` had volume +1 but reported its centre as outside (signed distance +0.5,
!< `is_point_inside` false). Cases:
!<   1. Fully inside-out cube (global flip).
!<   2. Cube with only some facets reversed (BFS flips, no global flip needed).
!<   3. Inside-out dragon (it replaced dragon-fine, issue #45) vs the outward original on a point grid: same sign at
!<      every point, same |d| to a relative 1e-12 (the flipped winding stores the
!<      vertices in another order, so the last bits of |d| may move).

program fossil_test_sanitize_inside_out

use fossil, only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

real(R8P), parameter :: REL_TOL = 1.e-12_R8P !< Relative |d| tolerance for the dragon comparison.
logical              :: ok_full, ok_partial, ok_dragon

ok_full    = check_cube(n_reversed=-1_I4P, label='cube, fully inside-out  ')
ok_partial = check_cube(n_reversed=5_I4P,  label='cube, 5 facets reversed ')
ok_dragon  = check_dragon()

print '(A,L1)', 'Are all tests passed? ', (ok_full .and. ok_partial .and. ok_dragon)
if (.not. (ok_full .and. ok_partial .and. ok_dragon)) error stop 1

contains
   function check_cube(n_reversed, label) result(ok)
   !< Reverse the first `n_reversed` facets of the unit cube (all if negative), analyze,
   !< sanitize, then check volume, signed distance and inside test at the centre.
   integer(I4P), intent(in)  :: n_reversed !< Facets to reverse (negative = all).
   character(*), intent(in)  :: label      !< Case label.
   logical                   :: ok         !< Check result.
   type(surface_stl_object)  :: surface    !< Cube.
   type(vector_R8P)          :: centre     !< Cube centre.
   real(R8P)                 :: d          !< Signed distance at the centre.
   logical                   :: inside     !< Inside test at the centre.

   call surface%load_from_file(file_name='src/tests/cube.stl', guess_format=.true.)
   if (n_reversed < 0_I4P) then
      call surface%reverse_normals
   else
      call reverse_first(surface, n_reversed)
   endif
   call surface%analyze           ! derived data built for the wrong winding, as after load_from_file
   call surface%sanitize
   centre = vector_R8P(0.5_R8P, 0.5_R8P, 0.5_R8P)
   d      = surface%distance(point=centre, is_signed=.true., is_square_root=.true., sign_algorithm=SIGN_PSEUDO_NORMAL)
   inside = surface%is_point_inside(point=centre, sign_algorithm=SIGN_PSEUDO_NORMAL)
   ok = (abs(surface%get_volume() - 1._R8P) < 1.e-12_R8P) .and. (abs(d + 0.5_R8P) < 1.e-12_R8P) .and. inside
   print '(A,A,F8.4,A,F8.4,A,L1,A,L1)', label, ': volume=', surface%get_volume(), ' d(centre)=', d, &
                                        ' inside=', inside, ' ok=', ok
   endfunction check_cube

   subroutine reverse_first(surface, n)
   !< Reverse the winding of facets 1..n through the public facet accessor.
   type(surface_stl_object), intent(inout) :: surface !< Surface.
   integer(I4P),             intent(in)    :: n       !< Facets to reverse.
   integer(I4P)                            :: i       !< Counter.

   do i = 1_I4P, min(n, surface%get_facets_number())
      call reverse_one(surface, i)
   enddo
   endsubroutine reverse_first

   subroutine reverse_one(surface, i)
   !< Reverse facet i.
   use fossil, only : facet_object
   type(surface_stl_object), intent(inout) :: surface !< Surface.
   integer(I4P),             intent(in)    :: i       !< Facet index.
   type(facet_object), pointer             :: f       !< Facet.

   f => surface%facet_at(i)
   call f%reverse_normal
   endsubroutine reverse_one

   function check_dragon() result(ok)
   !< Outward dragon vs an inside-out copy put through sanitize, on a 12^3 grid
   !< spanning the bounding box grown by 10%.
   logical                  :: ok              !< Check result.
   type(surface_stl_object) :: good, bad       !< Outward reference and inside-out copy.
   type(vector_R8P)         :: bmin, bmax, h, p !< Grid geometry and query point.
   real(R8P)                :: d_good, d_bad   !< Signed distances.
   integer(I4P), parameter  :: N = 12_I4P      !< Grid points per direction.
   integer(I4P)             :: i, j, k         !< Counters.
   integer(I4P)             :: n_sign, n_mag, n_inside !< Mismatch counters.
   type(facet_object), allocatable :: fs(:)    !< Copy of the sanitized facets.
   type(facet_object), pointer     :: fp       !< Facet.

   call good%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
   call good%sanitize
   allocate(fs(good%get_facets_number())) ! the sanitized dragon, copied (cheaper than a second load + sanitize)
   do i = 1, size(fs)
      fp => good%facet_at(i)
      fs(i) = fp
   enddo
   call bad%adopt_facets(facets=fs)
   call bad%reverse_normals
   call bad%analyze               ! fully analyzed inside-out dragon
   call bad%sanitize

   bmin = good%get_bmin()
   bmax = good%get_bmax()
   h    = bmax - bmin
   bmin = bmin - 0.1_R8P * h
   h    = 1.2_R8P * h / real(N - 1, R8P)
   n_sign = 0 ; n_mag = 0 ; n_inside = 0
   do k = 0, N - 1
      do j = 0, N - 1
         do i = 0, N - 1
            p%x = bmin%x + real(i, R8P) * h%x
            p%y = bmin%y + real(j, R8P) * h%y
            p%z = bmin%z + real(k, R8P) * h%z
            d_good = good%distance(point=p, is_signed=.true., is_square_root=.true.)
            d_bad  = bad%distance(point=p, is_signed=.true., is_square_root=.true.)
            if (sign(1._R8P, d_good) /= sign(1._R8P, d_bad)) n_sign = n_sign + 1
            if (abs(abs(d_good) - abs(d_bad)) > REL_TOL * max(abs(d_good), tiny(1._R8P))) n_mag = n_mag + 1
            if (good%is_point_inside(point=p) .neqv. bad%is_point_inside(point=p)) n_inside = n_inside + 1
         enddo
      enddo
   enddo
   ok = (n_sign == 0) .and. (n_mag == 0) .and. (n_inside == 0)
   print '(A,I0,A,I0,A,I0,A,I0,A,L1)', 'dragon inside-out vs outward on ', N**3, ' points: sign mismatches=', n_sign, &
                                       ' |d| mismatches=', n_mag, ' inside mismatches=', n_inside, ' ok=', ok
   endfunction check_dragon
endprogram fossil_test_sanitize_inside_out
