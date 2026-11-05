!< FOSSIL, mass properties of closed bodies (issue #59, roadmap #58 WP1).
module fossil_mass_properties
!< FOSSIL, mass properties of closed bodies (issue #59, roadmap #58 WP1).
!<
!< Exact polyhedral (3-D) and polygonal (2-D) integrals: the only error is floating-point round-off.
!<
!< 3-D: each facet (a, b, c) and a reference point o form a signed tetrahedron. With a, b, c relative to o and
!< D = a . (b x c): volume sum D/6, first moments sum D/24 (a + b + c), second moments
!< int x_i x_j = sum D/120 [sum_p v_pi v_pj + (sum_p v_pi)(sum_p v_pj)] (covariance form, Tonon 2004; the same
!< integrals as Mirtich 1996). 2-D: Green's theorem per segment (u_i, v_i) -> (u_j, v_j), c = u_i v_j - u_j v_i:
!< area sum c/2, first moments sum c (u_i + u_j)/6, int u^2 = sum c (u_i^2 + u_i u_j + u_j^2)/12,
!< int u v = sum c (2 u_i v_i + u_i v_j + u_j v_i + 2 u_j v_j)/24.
!<
!< The callers pass coordinates relative to a reference point near the body (the bounding-box centre), so the sums
!< do not cancel far from the origin. The sums are compensated (Neumaier) and run in a fixed order: O(eps) error,
!< independent of the number of terms and of the thread count. The routines are pure and take plain arrays.
use fossil_rigid_pose, only : rigid_pose_t
use penf, only : I4P, R8P

implicit none
private
public :: moments_accumulator_t
public :: accumulate_tetrahedron
public :: moments_to_properties_3d
public :: symmetric3_eigen
public :: inertia_to_world, principal_frame

type :: moments_accumulator_t
   !< Compensated sums of the raw moments about the reference point.
   !< 3-D: [V, Mx, My, Mz, Cxx, Cyy, Czz, Cxy, Cxz, Cyz]; 2-D: [A, Mu, Mv, Cuu, Cvv, Cuv, 0, 0, 0, 0].
   real(R8P) :: s(10) = 0._R8P !< Running sums.
   real(R8P) :: c(10) = 0._R8P !< Neumaier compensations.
endtype moments_accumulator_t

contains
   pure subroutine add_term(acc, k, x)
   !< Neumaier compensated addition of `x` to sum `k`.
   type(moments_accumulator_t), intent(inout) :: acc !< Accumulator.
   integer(I4P),                intent(in)    :: k   !< Sum index.
   real(R8P),                   intent(in)    :: x   !< Term.
   real(R8P)                                  :: t   !< New sum.

   t = acc%s(k) + x
   if (abs(acc%s(k)) >= abs(x)) then
      acc%c(k) = acc%c(k) + ((acc%s(k) - t) + x)
   else
      acc%c(k) = acc%c(k) + ((x - t) + acc%s(k))
   endif
   acc%s(k) = t
   endsubroutine add_term

   pure function accumulator_total(acc) result(m)
   !< The compensated totals.
   type(moments_accumulator_t), intent(in) :: acc   !< Accumulator.
   real(R8P)                               :: m(10) !< Totals.

   m = acc%s + acc%c
   endfunction accumulator_total

   pure subroutine accumulate_tetrahedron(acc, a, b, c)
   !< Add the signed tetrahedron (o, a, b, c), vertices relative to the reference point o.
   type(moments_accumulator_t), intent(inout) :: acc     !< Accumulator.
   real(R8P),                   intent(in)    :: a(3)    !< First vertex.
   real(R8P),                   intent(in)    :: b(3)    !< Second vertex.
   real(R8P),                   intent(in)    :: c(3)    !< Third vertex.
   real(R8P)                                  :: d       !< Triple product a . (b x c).
   real(R8P)                                  :: s(3)    !< a + b + c.
   real(R8P)                                  :: q       !< D / 120.
   integer(I4P)                               :: i, j, k !< Counters.
   integer(I4P), parameter                    :: II(3) = [1_I4P, 1_I4P, 2_I4P] !< Off-diagonal rows.
   integer(I4P), parameter                    :: JJ(3) = [2_I4P, 3_I4P, 3_I4P] !< Off-diagonal columns.

   d = a(1) * (b(2) * c(3) - b(3) * c(2)) + a(2) * (b(3) * c(1) - b(1) * c(3)) + a(3) * (b(1) * c(2) - b(2) * c(1))
   s = a + b + c
   q = d / 120._R8P
   call add_term(acc, 1_I4P, d / 6._R8P)
   do i = 1, 3
      call add_term(acc, 1_I4P + i, d / 24._R8P * s(i))
      call add_term(acc, 4_I4P + i, q * (a(i) * a(i) + b(i) * b(i) + c(i) * c(i) + s(i) * s(i)))
   enddo
   do k = 1, 3
      i = II(k) ; j = JJ(k)
      call add_term(acc, 7_I4P + k, q * (a(i) * a(j) + b(i) * b(j) + c(i) * c(j) + s(i) * s(j)))
   enddo
   endsubroutine accumulate_tetrahedron

   pure subroutine moments_to_properties_3d(acc, reference, density, about, volume, mass, center, inertia)
   !< Volume, mass, centre of mass and inertia tensor (about the centre, or about `about`) from the raw moments.
   type(moments_accumulator_t), intent(in)           :: acc          !< Raw moments about `reference`.
   real(R8P),                   intent(in)           :: reference(3) !< Reference point o.
   real(R8P),                   intent(in)           :: density      !< Uniform density.
   real(R8P),                   intent(in), optional :: about(3)     !< Point of the tensor (default: the centre).
   real(R8P),                   intent(out)          :: volume       !< Volume.
   real(R8P),                   intent(out)          :: mass         !< Mass.
   real(R8P),                   intent(out)          :: center(3)    !< Centre of mass.
   real(R8P),                   intent(out)          :: inertia(3,3) !< Inertia tensor.
   real(R8P)                                         :: m(10)        !< Totals.
   real(R8P)                                         :: d(3)         !< Centre relative to o.
   real(R8P)                                         :: e(3)         !< Centre relative to `about`.
   real(R8P)                                         :: cm(3,3)      !< Second moments about the centre (or `about`).
   integer(I4P)                                      :: i            !< Counter.

   m = accumulator_total(acc)
   volume = m(1)
   mass = density * volume
   d = m(2:4) / volume
   center = reference + d
   cm(1,1) = m(5) - volume * d(1) * d(1)
   cm(2,2) = m(6) - volume * d(2) * d(2)
   cm(3,3) = m(7) - volume * d(3) * d(3)
   cm(1,2) = m(8) - volume * d(1) * d(2) ; cm(2,1) = cm(1,2)
   cm(1,3) = m(9) - volume * d(1) * d(3) ; cm(3,1) = cm(1,3)
   cm(2,3) = m(10) - volume * d(2) * d(3) ; cm(3,2) = cm(2,3)
   if (present(about)) then
      e = center - about
      cm = cm + volume * spread(e, 2, 3) * spread(e, 1, 3)
   endif
   inertia = -density * cm
   do i = 1, 3
      inertia(i,i) = density * (cm(1,1) + cm(2,2) + cm(3,3) - cm(i,i))
   enddo
   endsubroutine moments_to_properties_3d

   pure subroutine symmetric3_eigen(a, w, v)
   !< Eigen-decomposition of a symmetric 3x3 matrix by cyclic Jacobi rotations: a v(:,k) = w(k) v(:,k).
   !< `w` ascending; `v` orthonormal and right-handed (det = +1). Sign convention, so the result is deterministic:
   !< the largest-magnitude component of columns 1 and 2 is positive (the first one on ties), column 3 = col1 x col2.
   !< With repeated eigenvalues the axes are not unique: any orthonormal set is a valid answer.
   real(R8P), intent(in)  :: a(3,3)    !< Symmetric matrix.
   real(R8P), intent(out) :: w(3)      !< Eigenvalues, ascending.
   real(R8P), intent(out) :: v(3,3)    !< Eigenvectors (columns).
   real(R8P)              :: b(3,3)    !< Working copy.
   real(R8P)              :: off       !< Off-diagonal norm.
   real(R8P)              :: scale     !< Matrix norm.
   real(R8P)              :: theta     !< Rotation parameter.
   real(R8P)              :: t, c, s   !< Tangent, cosine, sine.
   real(R8P)              :: bip, biq  !< Row entries.
   real(R8P)              :: vip, viq  !< Eigenvector entries.
   real(R8P)              :: tmp(3)    !< Swap buffer.
   real(R8P)              :: wt        !< Swap buffer.
   integer(I4P)           :: sweep     !< Sweep counter.
   integer(I4P)           :: p, q, i   !< Counters.
   integer(I4P)           :: k, kmax   !< Counters.

   b = 0.5_R8P * (a + transpose(a))
   v = 0._R8P ; v(1,1) = 1._R8P ; v(2,2) = 1._R8P ; v(3,3) = 1._R8P
   scale = sqrt(sum(b * b))
   do sweep = 1, 50
      off = sqrt(b(1,2) ** 2 + b(1,3) ** 2 + b(2,3) ** 2)
      if (off <= epsilon(1._R8P) * 1.e-3_R8P * scale .or. off == 0._R8P) exit
      do p = 1, 2
         do q = p + 1, 3
            if (b(p,q) == 0._R8P) cycle
            theta = (b(q,q) - b(p,p)) / (2._R8P * b(p,q))
            t = sign(1._R8P, theta) / (abs(theta) + sqrt(theta * theta + 1._R8P))
            c = 1._R8P / sqrt(t * t + 1._R8P)
            s = t * c
            do i = 1, 3
               bip = b(i,p) ; biq = b(i,q)
               b(i,p) = c * bip - s * biq
               b(i,q) = s * bip + c * biq
            enddo
            do i = 1, 3
               bip = b(p,i) ; biq = b(q,i)
               b(p,i) = c * bip - s * biq
               b(q,i) = s * bip + c * biq
            enddo
            do i = 1, 3
               vip = v(i,p) ; viq = v(i,q)
               v(i,p) = c * vip - s * viq
               v(i,q) = s * vip + c * viq
            enddo
         enddo
      enddo
   enddo
   w = [b(1,1), b(2,2), b(3,3)]
   ! Sort ascending (selection sort on three entries).
   do p = 1, 2
      k = p
      do q = p + 1, 3
         if (w(q) < w(k)) k = q
      enddo
      if (k /= p) then
         wt = w(p) ; w(p) = w(k) ; w(k) = wt
         tmp = v(:,p) ; v(:,p) = v(:,k) ; v(:,k) = tmp
      endif
   enddo
   do p = 1, 2
      kmax = 1
      do k = 2, 3
         if (abs(v(k,p)) > abs(v(kmax,p))) kmax = k
      enddo
      if (v(kmax,p) < 0._R8P) v(:,p) = -v(:,p)
   enddo
   v(1,3) = v(2,1) * v(3,2) - v(3,1) * v(2,2)
   v(2,3) = v(3,1) * v(1,2) - v(1,1) * v(3,2)
   v(3,3) = v(1,1) * v(2,2) - v(2,1) * v(1,2)
   endsubroutine symmetric3_eigen

   pure function inertia_to_world(inertia, pose) result(inertia_w)
   !< Inertia tensor of a body placed at `pose` (x_world = R x_body + t): R I R^T. The tensor is about the centre of mass
   !< (or any body point); it moves with that point, so the translation does not enter.
   real(R8P),          intent(in) :: inertia(3,3)   !< Tensor in the body frame.
   type(rigid_pose_t), intent(in) :: pose           !< Pose.
   real(R8P)                      :: inertia_w(3,3) !< Tensor in the world frame.

   inertia_w = matmul(pose%r, matmul(inertia, transpose(pose%r)))
   endfunction inertia_to_world

   pure function principal_frame(center, principal_axes) result(pose)
   !< The pose of the principal frame: x = R y + t maps principal coordinates y to the frame of `center`, with
   !< R = the principal axes (columns) and t = the centre of mass. In that frame the inertia tensor is diagonal.
   real(R8P),          intent(in) :: center(3)           !< Centre of mass.
   real(R8P),          intent(in) :: principal_axes(3,3) !< Principal axes (columns, right-handed).
   type(rigid_pose_t)             :: pose                !< Principal frame.

   pose%r = principal_axes
   pose%t = center
   endfunction principal_frame
endmodule fossil_mass_properties
