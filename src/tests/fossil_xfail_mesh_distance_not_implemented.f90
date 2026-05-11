!< FOSSIL, expected-failure test: compute_mesh_distance raises `not implemented`.
!<
!< Regression guard for audit S5 + H. The routine has a stub body and now `error stop`s
!< to prevent silent garbage results. If a future change re-implements it, this test
!< must be deleted or converted to a positive-path test. Until then, calling it from
!< the public API is expected to terminate with a non-zero exit code.

program fossil_xfail_mesh_distance_not_implemented

use fossil, only : surface_stl_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)      :: surface           !< Empty surface; dispatch test only.
type(vector_R8P), allocatable :: mesh(:, :, :)     !< Trivial 1x1x1 mesh.
real(R8P),        allocatable :: distance(:, :, :) !< Distance output.

allocate(mesh(1, 1, 1))
allocate(distance(1, 1, 1))

call surface%compute_mesh_distance(mesh=mesh, distance=distance)

! Unreachable.
print '(A)', 'FAIL: compute_mesh_distance did not error_stop'
endprogram fossil_xfail_mesh_distance_not_implemented
