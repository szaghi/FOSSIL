!< Shape operator at a foot point (roadmap #109 WP1 step 03): plain-real routines shared by the host and the device.
module fossil_shape_operator
!< Shape operator at a foot point (roadmap #109 WP1 step 03): plain-real routines shared by the host and the device.
!<
!< The per-vertex principal curvatures k1 >= k2 and directions d1, d2 (`surface%build_curvature`, from the jet fit of
!< `fossil_curvature%build_jet_curvature`) define at every vertex the tangent tensor S_v = k1 d1 d1^T + k2 d2 d2^T. At a
!< foot point on facet f, `foot_principal_oac` interpolates the three S_v with the foot's barycentric weights (the rule of
!< `query_many(curvatures=)`: clamped to the triangle, along the longest edge on a near-collinear facet), restricts the
!< result to the tangent plane of the foot normal and eigen-decomposes the 2 x 2 restriction in closed form
!< (`eig_sym2_oac`, square roots only: no transcendental function, whose device and host results differ by an ulp).
!< Sign: positive where the surface is convex seen from the side the normals point to (a sphere of radius R: k1 = k2 =
!< 1/R); the normal curvature in a tangent direction t is k1 (t.d1)^2 + k2 (t.d2)^2 (`normal_curvature`).
!<
!< Every routine is `!$acc routine seq`, takes plain reals and integers (no vector_R8P, no class) and runs unchanged on
!< the host and on the device, so `query_many` and `query_many_device` agree bit for bit (with the -gpu=nofma /
!< -ffp-contract=off flags of the fobos modes).
use penf, only : I4P, R8P

implicit none
private
public :: eig_sym2_oac
public :: tangent_frame_oac
public :: foot_principal_oac
public :: curve_foot_curvature_oac
public :: normal_curvature

contains
   pure subroutine eig_sym2_oac(p, q, r, k1, k2, c, s)
   !< Eigen-decomposition of the symmetric 2 x 2 matrix [[p, q], [q, r]] in closed form: k1 >= k2 and the unit eigenvector
   !< (c, s) of k1 (the one of k2 is (-s, c)). The eigenvector is taken from the column without cancellation: (sq + h, q)
   !< if h = (p - r) / 2 >= 0, else (q, sq - h), sq = sqrt(h^2 + q^2). An umbilic (h = q = 0, k1 = k2) has every
   !< direction as eigenvector: it returns (1, 0), the first basis vector, deterministically. Square roots and the four
   !< operations only, so the device and the host give the same bits.
   !$acc routine seq
   real(R8P), intent(in)  :: p, q, r !< Matrix entries.
   real(R8P), intent(out) :: k1, k2  !< Eigenvalues, k1 >= k2.
   real(R8P), intent(out) :: c, s    !< Unit eigenvector of k1.
   real(R8P)              :: h, m, sq, a, b, nrm !< Work.

   h = 0.5_R8P * (p - r) ; m = 0.5_R8P * (p + r)
   sq = sqrt(h * h + q * q)
   k1 = m + sq ; k2 = m - sq
   if (h >= 0._R8P) then
      a = sq + h ; b = q
   else
      a = q ; b = sq - h
   endif
   nrm = sqrt(a * a + b * b)
   if (nrm > 0._R8P) then
      c = a / nrm ; s = b / nrm
   else ! umbilic (or NaN input): the first basis vector
      c = 1._R8P ; s = 0._R8P
   endif
   endsubroutine eig_sym2_oac

   pure subroutine tangent_frame_oac(nx, ny, nz, t1x, t1y, t1z, t2x, t2y, t2z)
   !< Orthonormal tangent frame (t1, t2) of a unit normal n with (t1, t2, n) right-handed: t1 = e_x x n normalised if
   !< |n_x| < 0.9, else e_y x n normalised (the frame of the jet fit); t2 = n x t1.
   !$acc routine seq
   real(R8P), intent(in)  :: nx, ny, nz           !< Unit normal.
   real(R8P), intent(out) :: t1x, t1y, t1z        !< First tangent.
   real(R8P), intent(out) :: t2x, t2y, t2z        !< Second tangent.
   real(R8P)              :: l                    !< Length.

   if (abs(nx) < 0.9_R8P) then
      t1x = 0._R8P ; t1y = -nz ; t1z = ny
   else
      t1x = nz ; t1y = 0._R8P ; t1z = -nx
   endif
   l = sqrt(t1x * t1x + t1y * t1y + t1z * t1z)
   t1x = t1x / l ; t1y = t1y / l ; t1z = t1z / l
   t2x = ny * t1z - nz * t1y ; t2y = nz * t1x - nx * t1z ; t2z = nx * t1y - ny * t1x
   endsubroutine tangent_frame_oac

   pure subroutine foot_principal_oac(tv, cvid, ck, cd, f, px, py, pz, nx, ny, nz, k1, k2, d1x, d1y, d1z, d2x, d2y, d2z)
   !< Principal curvatures k1 >= k2 and directions d1, d2 at a foot point p on facet f with normal n (the pseudo-normal of
   !< the closest feature): the per-vertex tensors S_v = k1_v d1_v d1_v^T + k2_v d2_v d2_v^T of the facet's three
   !< vertices weighted by the barycentric coordinates of p (clamped to the triangle and renormalised; on a near-collinear
   !< facet, Gram determinant <= 1e-12 d11 d22, linear along the longest edge at p's projection, the mean of the three if
   !< they coincide: the rule of `query_many(curvatures=)`), restricted to the tangent plane of n in the frame of
   !< `tangent_frame_oac`, then `eig_sym2_oac`: d1 = c t1 + s t2, d2 = n x d1 (an umbilic gives d1 = t1). An edge or a
   !< vertex foot takes the same rule (its weights have zeros). A vertex value that is NaN (a failed jet fit) makes the
   !< result NaN, as `curvatures=`. No facet (f outside 1..size(tv, 2)) or a zero normal: every output 0.
   !$acc routine seq
   real(R8P),    intent(in)  :: tv(:,:)              !< Facet vertices [9, nf].
   integer(I4P), intent(in)  :: cvid(:,:)            !< Vertex index of each facet corner into ck, cd [3, nf].
   real(R8P),    intent(in)  :: ck(:,:)              !< Per-vertex principal curvatures [2, nv].
   real(R8P),    intent(in)  :: cd(:,:,:)            !< Per-vertex principal directions [3, 2, nv].
   integer(I4P), intent(in)  :: f                    !< Facet (0: none).
   real(R8P),    intent(in)  :: px, py, pz           !< Foot point.
   real(R8P),    intent(in)  :: nx, ny, nz           !< Foot normal.
   real(R8P),    intent(out) :: k1, k2               !< Principal curvatures, k1 >= k2.
   real(R8P),    intent(out) :: d1x, d1y, d1z        !< Direction of k1.
   real(R8P),    intent(out) :: d2x, d2y, d2z        !< Direction of k2.
   real(R8P) :: e1x, e1y, e1z, e2x, e2y, e2z, wx, wy, wz !< Edges, foot offset.
   real(R8P) :: d11, d12, d22, d1w, d2w, den, b1, b2, b3, bs !< Gram system, weights.
   real(R8P) :: l1, l2, l3, t, lm                    !< Squared edge lengths, edge parameter, largest.
   real(R8P) :: sxx, syy, szz, sxy, sxz, syz          !< Interpolated tensor.
   real(R8P) :: ux, uy, uz, vx, vy, vz, mx, my, mz, ln !< Unit normal, frame, work.
   real(R8P) :: w, ka, kb, ax, ay, az, bx, by, bz     !< Vertex weight, values, directions.
   real(R8P) :: m11, m12, m22, c, s                   !< Restriction, eigenvector.
   integer(I4P) :: a, j, ib                           !< Longest edge, corner, index.

   k1 = 0._R8P ; k2 = 0._R8P
   d1x = 0._R8P ; d1y = 0._R8P ; d1z = 0._R8P ; d2x = 0._R8P ; d2y = 0._R8P ; d2z = 0._R8P
   if (f <= 0_I4P .or. f > size(tv, 2)) return
   ln = sqrt(nx * nx + ny * ny + nz * nz)
   if (.not. (ln > 0._R8P)) return
   mx = nx / ln ; my = ny / ln ; mz = nz / ln
   ! barycentric weights of the foot
   e1x = tv(4, f) - tv(1, f) ; e1y = tv(5, f) - tv(2, f) ; e1z = tv(6, f) - tv(3, f)
   e2x = tv(7, f) - tv(1, f) ; e2y = tv(8, f) - tv(2, f) ; e2z = tv(9, f) - tv(3, f)
   wx = px - tv(1, f) ; wy = py - tv(2, f) ; wz = pz - tv(3, f)
   d11 = e1x * e1x + e1y * e1y + e1z * e1z ; d12 = e1x * e2x + e1y * e2y + e1z * e2z
   d22 = e2x * e2x + e2y * e2y + e2z * e2z
   d1w = e1x * wx + e1y * wy + e1z * wz ; d2w = e2x * wx + e2y * wy + e2z * wz
   den = d11 * d22 - d12 * d12
   if (den > 1.e-12_R8P * d11 * d22) then
      b2 = max((d22 * d1w - d12 * d2w) / den, 0._R8P)
      b3 = max((d11 * d2w - d12 * d1w) / den, 0._R8P)
      b1 = max(1._R8P - b2 - b3, 0._R8P)
      bs = b1 + b2 + b3
      b1 = b1 / bs ; b2 = b2 / bs ; b3 = b3 / bs
   else ! near-collinear facet: linear along the longest edge (a -> a + 1)
      l1 = d11
      l2 = (tv(7, f) - tv(4, f))**2 + (tv(8, f) - tv(5, f))**2 + (tv(9, f) - tv(6, f))**2
      l3 = d22
      a = 1_I4P ; lm = l1
      if (l2 > lm) then
         a = 2_I4P ; lm = l2
      endif
      if (l3 > lm) then
         a = 3_I4P ; lm = l3
      endif
      if (lm <= 0._R8P) then
         b1 = 1._R8P / 3._R8P ; b2 = b1 ; b3 = b1
      else
         j = mod(a, 3_I4P) + 1_I4P ! far end of the edge a -> j
         t = ((px - tv(3 * a - 2, f)) * (tv(3 * j - 2, f) - tv(3 * a - 2, f)) + &
              (py - tv(3 * a - 1, f)) * (tv(3 * j - 1, f) - tv(3 * a - 1, f)) + &
              (pz - tv(3 * a, f)) * (tv(3 * j, f) - tv(3 * a, f))) / lm
         t = min(max(t, 0._R8P), 1._R8P)
         b1 = 0._R8P ; b2 = 0._R8P ; b3 = 0._R8P
         if (a == 1_I4P) then
            b1 = 1._R8P - t ; b2 = t
         elseif (a == 2_I4P) then
            b2 = 1._R8P - t ; b3 = t
         else
            b3 = 1._R8P - t ; b1 = t
         endif
      endif
   endif
   ! interpolated tensor
   sxx = 0._R8P ; syy = 0._R8P ; szz = 0._R8P ; sxy = 0._R8P ; sxz = 0._R8P ; syz = 0._R8P
   do j = 1_I4P, 3_I4P
      w = b1
      if (j == 2_I4P) w = b2
      if (j == 3_I4P) w = b3
      ib = cvid(j, f)
      ka = w * ck(1, ib) ; kb = w * ck(2, ib)
      ax = cd(1, 1, ib) ; ay = cd(2, 1, ib) ; az = cd(3, 1, ib)
      bx = cd(1, 2, ib) ; by = cd(2, 2, ib) ; bz = cd(3, 2, ib)
      sxx = sxx + (ka * ax * ax + kb * bx * bx) ; syy = syy + (ka * ay * ay + kb * by * by)
      szz = szz + (ka * az * az + kb * bz * bz) ; sxy = sxy + (ka * ax * ay + kb * bx * by)
      sxz = sxz + (ka * ax * az + kb * bx * bz) ; syz = syz + (ka * ay * az + kb * by * bz)
   enddo
   ! restriction to the tangent plane of the foot normal
   call tangent_frame_oac(mx, my, mz, ux, uy, uz, vx, vy, vz)
   m11 = ux * (sxx * ux + sxy * uy + sxz * uz) + uy * (sxy * ux + syy * uy + syz * uz) + uz * (sxz * ux + syz * uy + szz * uz)
   m12 = ux * (sxx * vx + sxy * vy + sxz * vz) + uy * (sxy * vx + syy * vy + syz * vz) + uz * (sxz * vx + syz * vy + szz * vz)
   m22 = vx * (sxx * vx + sxy * vy + sxz * vz) + vy * (sxy * vx + syy * vy + syz * vz) + vz * (sxz * vx + syz * vy + szz * vz)
   call eig_sym2_oac(m11, m12, m22, k1, k2, c, s)
   d1x = c * ux + s * vx ; d1y = c * uy + s * vy ; d1z = c * uz + s * vz
   d2x = c * vx - s * ux ; d2y = c * vy - s * uy ; d2z = c * vz - s * uz
   endsubroutine foot_principal_oac

   pure function curve_foot_curvature_oac(u, v, next, ck, sg, pu, pv) result(kappa)
   !< Curvature at a foot point (pu, pv) on segment sg (its first vertex) of a curve: the two per-vertex curvatures
   !< (`curve%build_curvature`) interpolated linearly in the foot's parameter along the segment, clamped to [0, 1]. 0 for
   !< no segment.
   !$acc routine seq
   real(R8P),    intent(in) :: u(:), v(:) !< Vertices [nv].
   integer(I4P), intent(in) :: next(:)    !< Next vertex [nv] (0: none).
   real(R8P),    intent(in) :: ck(:)      !< Per-vertex curvature [nv].
   integer(I4P), intent(in) :: sg         !< Segment (0: none).
   real(R8P),    intent(in) :: pu, pv     !< Foot point.
   real(R8P)                :: kappa      !< Curvature.
   real(R8P)                :: du, dv, l2, t !< Segment, parameter.
   integer(I4P)             :: j          !< End vertex.

   kappa = 0._R8P
   if (sg <= 0_I4P .or. sg > size(u)) return
   j = next(sg)
   if (j <= 0_I4P) return
   du = u(j) - u(sg) ; dv = v(j) - v(sg)
   l2 = du * du + dv * dv
   t = 0._R8P
   if (l2 > 0._R8P) t = min(max(((pu - u(sg)) * du + (pv - v(sg)) * dv) / l2, 0._R8P), 1._R8P)
   kappa = (1._R8P - t) * ck(sg) + t * ck(j)
   endfunction curve_foot_curvature_oac

   pure function normal_curvature(k, d, t) result(kn)
   !< Normal curvature in a tangent direction t from the principal curvatures k = (k1, k2) and directions d(:, 1), d(:, 2)
   !< (as `query_many(principal_curvatures=, principal_directions=)` returns them): Euler's formula
   !< kn = (k1 (t.d1)^2 + k2 (t.d2)^2) / (t.t); t need not be unit, its normal component is ignored. A sphere of radius R
   !< gives 1/R in every direction, so the curvature condition d_n v_t = v_t / R of a wall closure reads
   !< d_n v_t = kn v_t with t the direction of v_t. 0 for t = 0. For a curve pass k(1:1) and d(:, 1:1) of size 2: kn =
   !< k1 (t.d1)^2 / (t.t).
   !$acc routine seq
   real(R8P), intent(in) :: k(:)   !< Principal curvatures [2] (a curve: [1]).
   real(R8P), intent(in) :: d(:,:) !< Principal directions [3, 2] (a curve: [2, 1]).
   real(R8P), intent(in) :: t(:)   !< Tangent direction [3] (a curve: [2]).
   real(R8P)             :: kn     !< Normal curvature.
   real(R8P)             :: tt, a  !< |t|^2, projection.
   integer(I4P)          :: i, l   !< Counters.

   kn = 0._R8P
   tt = 0._R8P
   do l = 1_I4P, size(t, kind=I4P)
      tt = tt + t(l) * t(l)
   enddo
   if (.not. (tt > 0._R8P)) return
   do i = 1_I4P, size(k, kind=I4P)
      a = 0._R8P
      do l = 1_I4P, size(t, kind=I4P)
         a = a + t(l) * d(l, i)
      enddo
      kn = kn + k(i) * a * a
   enddo
   kn = kn / tt
   endfunction normal_curvature
endmodule fossil_shape_operator
