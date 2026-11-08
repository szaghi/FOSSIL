!< FOSSIL validation: exact 2-D swept geometry under rotation on dragon sections (issue #86).
!<
!< The dragon is repaired and placed as in tutorial chapters 2 and 3 (`sanitize(fill_holes=.true.)`, box centre to the
!< origin, a quarter turn about z, longest side 1), sliced at z = -0.35, -0.20, -0.10, -0.05 (`merge_collinear`), and
!< each section turned by pi/18 about the origin on the 48 x 48 grid x = y = -0.6 + 0.025 k. Every section must take
!< the exact rotation path (`sw%exact`) with the rigid-field residual at round-off (<= 1e-14 in every closed cell),
!< the projected walls summing to zero and the swept areas summing to zero. Before the fix (base 9335a01) the sections
!< at z = -0.20 (and -0.35 in the docs build) gave `exact = F` with residual 7e-18: a piece's Gauss values stopped at
!< the round-off floor of its evaluations, above the strict per-piece target.
program fossil_val_swept_rotation_2d

use fossil, only : surface_stl_object, curve_object, pose_2d, slice_surface, swept_cells_2d_t, swept_cells_2d, &
                   STATUS_OK
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P, ez_R8P

implicit none

real(R8P), parameter     :: PI = 4._R8P * atan(1._R8P)
type(surface_stl_object) :: surface
type(curve_object)       :: section
type(swept_cells_2d_t)   :: sw
type(vector_R8P)         :: bmin, bmax
real(R8P)                :: x(49), zs(4), r, w(2), scale
integer(I4P)             :: i, k, st, sst
logical                  :: ok(4)

call surface%load_from_file(file_name='src/tests/dragon.stl', guess_format=.true.)
call surface%sanitize(fill_holes=.true.)
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
call surface%translate(delta=-0.5_R8P * (bmin + bmax))
call surface%rotate(axis=ez_R8P, angle=PI / 2._R8P)
call surface%analyze
bmin = surface%get_bmin() ; bmax = surface%get_bmax()
scale = 1._R8P / maxval([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
call surface%resize(x=scale, y=scale, z=scale)
call surface%analyze
x = [(-0.6_R8P + 0.025_R8P * i, i = 0, 48)]
zs = [-0.35_R8P, -0.20_R8P, -0.10_R8P, -0.05_R8P]
do k = 1, 4
   call slice_surface(surface, 3_I4P, zs(k), section, merge_collinear=.true., status=sst)
   call swept_cells_2d(section, pose_2d(angle=0._R8P), pose_2d(angle=PI / 18._R8P), x, x, sw, status=st)
   r = maxval(abs(sw%residual), mask=sw%alpha0 >= 0._R8P .and. sw%alpha1 >= 0._R8P)
   w = [sum(sw%wall(1, :, :)), sum(sw%wall(2, :, :))]
   ok(k) = sst == STATUS_OK .and. st == STATUS_OK .and. sw%exact .and. r <= 1e-14_R8P .and. &
           all(abs(w) <= 1e-14_R8P) .and. abs(sum(sw%swept)) <= 1e-14_R8P
   print '(A,F6.2,A,I0,A,L1,A,ES9.2,A,2ES9.1,A,ES9.1,A,I0,A,L1)', 'z = ', zs(k), ': ', &
         section%get_vertices_number(), ' vertices, exact ', sw%exact, ', residual ', r, ', walls ', w, &
         ', sum of swept ', sum(sw%swept), ', evaluations/face ', sw%evaluations, ', ok ', ok(k)
enddo
print '(A,L1)', 'Are all tests passed? ', all(ok)
if (.not. all(ok)) error stop 1
endprogram fossil_val_swept_rotation_2d
