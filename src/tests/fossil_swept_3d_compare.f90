!< FOSSIL test helper: bitwise comparison of two 3-D swept geometries (issue #100 step 4).
module fossil_swept_3d_compare
!< FOSSIL test helper: bitwise comparison of two 3-D swept geometries (issue #100 step 4).
!<
!< `swept_cells_3d_mismatches(a, b)` counts the differing values of every output of `swept_cells_3d_t`: the cell
!< counts, alpha0, alpha1, swept, face_x, face_y, face_z, wall, residual (bit patterns), cleared, covered, exact,
!< error_estimate (bits), evaluations, self_intersecting; an array allocated in one and not the other, or of another
!< shape, counts its size (at least 1). Used by `fossil_test_swept_3d_device_twin` (the CI host twin) and the distance
!< bench (`swept cells 3-D device vs serial mismatches`).
use fossil, only : swept_cells_3d_t
use penf,   only : I4P, I8P, R8P

implicit none
private
public :: swept_cells_3d_mismatches

contains
   function swept_cells_3d_mismatches(a, b) result(n)
   !< Differing values of every output of a and b (see the module notes).
   type(swept_cells_3d_t), intent(in) :: a, b !< Swept geometries.
   integer(I4P)                       :: n    !< Mismatches.

   n = 0_I4P
   if (a%nx /= b%nx) n = n + 1_I4P
   if (a%ny /= b%ny) n = n + 1_I4P
   if (a%nz /= b%nz) n = n + 1_I4P
   n = n + real3(a%alpha0, b%alpha0) + real3(a%alpha1, b%alpha1) + real3(a%swept, b%swept)
   n = n + real3(a%face_x, b%face_x) + real3(a%face_y, b%face_y) + real3(a%face_z, b%face_z)
   n = n + real3(a%residual, b%residual) + logic3(a%cleared, b%cleared) + logic3(a%covered, b%covered)
   if (allocated(a%wall) .neqv. allocated(b%wall)) then
      n = n + 1_I4P
   elseif (allocated(a%wall)) then
      if (any(shape(a%wall) /= shape(b%wall))) then
         n = n + max(1_I4P, int(size(a%wall), I4P))
      else
         n = n + int(count(transfer(a%wall, 0_I8P, size(a%wall)) /= transfer(b%wall, 0_I8P, size(b%wall))), I4P)
      endif
   endif
   if (a%exact .neqv. b%exact) n = n + 1_I4P
   if (transfer(a%error_estimate, 0_I8P) /= transfer(b%error_estimate, 0_I8P)) n = n + 1_I4P
   if (a%evaluations /= b%evaluations) n = n + 1_I4P
   if (a%self_intersecting .neqv. b%self_intersecting) n = n + 1_I4P
   endfunction swept_cells_3d_mismatches

   function real3(a, b) result(n)
   !< Differing bit patterns of two rank-3 real arrays (allocation or shape mismatch: the size, at least 1).
   real(R8P), allocatable, intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.
   integer(I4P)                       :: n                  !< Mismatches.

   n = 0_I4P
   if (allocated(a) .neqv. allocated(b)) then
      n = 1_I4P
   elseif (allocated(a)) then
      if (any(shape(a) /= shape(b))) then
         n = max(1_I4P, int(size(a), I4P))
      else
         n = int(count(transfer(a, 0_I8P, size(a)) /= transfer(b, 0_I8P, size(b))), I4P)
      endif
   endif
   endfunction real3

   function logic3(a, b) result(n)
   !< Differing values of two rank-3 logical arrays (allocation or shape mismatch: the size, at least 1).
   logical, allocatable, intent(in) :: a(:,:,:), b(:,:,:) !< Arrays.
   integer(I4P)                     :: n                  !< Mismatches.

   n = 0_I4P
   if (allocated(a) .neqv. allocated(b)) then
      n = 1_I4P
   elseif (allocated(a)) then
      if (any(shape(a) /= shape(b))) then
         n = max(1_I4P, int(size(a), I4P))
      else
         n = int(count(a .neqv. b), I4P)
      endif
   endif
   endfunction logic3
endmodule fossil_swept_3d_compare
