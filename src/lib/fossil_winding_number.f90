!< FOSSIL, generalized / fast winding number on a triangulated surface.

module fossil_winding_number
!< FOSSIL, generalized / fast winding number on a triangulated surface.
!<
!< Implements the generalized winding number of Jacobson, Kavan & Sorkine-Hornung
!< (SIGGRAPH 2013) with optional hierarchical evaluation following Barill, Dickson,
!< Schmidt, Levin & Jacobson (SIGGRAPH 2018).
!<
!< For a triangulated surface S the winding number of a query point q is
!<
!<     w(q) = (1/4 pi) * sum_{f in S} Omega_f(q)
!<
!< where Omega_f(q) is the signed solid angle subtended by triangle f at q
!< (Van Oosterom-Strackee). For a closed, outward-oriented surface w(q) = 1
!< strictly inside, 0 strictly outside, and intermediate on the boundary or for
!< open / non-watertight meshes — which is precisely why this is preferred over
!< the parity / Baerentzen-Aanaes sign tests on dirty STL.
!<
!< The hierarchical (fast) variant (Barill et al. 2018) replaces a distant subtree by a Taylor
!< expansion of its solid angle about the node centre c (r' = c - q):
!<
!<     Omega_node(q) ~= D . r' / |r'|^3  -  (3 r'^T C r' - |r'|^2 tr C) / |r'|^5
!<
!< with the dipole D = sum_f area_f n_f and the second moment C = sum_f area_f n_f (x_f - c)^T,
!< admissible when |q - c| > beta * (bbox diagonal of the node). D, C, c and the diagonal live in
!< the persistent per-node cache `surface%wn_payload` (`wn_node_t`), built once by `analyze`
!< (issue #29).
!<
!< Two defects of the former implementation are fixed here (issue #29):
!<  - the far-field term had the wrong sign (it used q - c instead of c - q): the error was ~2 on
!<    bunny at beta = 2, i.e. the result was meaningless wherever the far field was used;
!<  - the cache was rebuilt over the whole tree on every call (O(N) per query) and, on octrees,
!<    counted every facet once per tree level (the octree stores a facet at each level).
!<
!< Measured |w_fast - w_exact| with the second-order expansion: beta = 2: bunny 8.7e-3,
!< dragon-fine 9.8e-4; beta = 4: 1.7e-3 and 7.2e-5. Ample for an inside/outside decision at
!< w = 0.5; pass beta <= 0 for the exact O(N) sum.

use fossil_aabb_tree_object, only : aabb_tree_object, AABB_TREE_SAH_BVH, wn_node_t, bvh_winding_number_oac
use fossil_facet_object,     only : facet_object
use fossil_utils,            only : PI
use penf,                    only : I4P, R8P
use vecfor,                  only : vector_R8P

implicit none
private
public :: winding_number

real(R8P),    parameter :: FOUR_PI      = 4._R8P * PI
real(R8P),    parameter :: DEFAULT_BETA = 2._R8P
integer(I4P), parameter :: MAX_CHILDREN = 8_I4P   !< enumerate_children buffer (octree fan-out).
integer(I4P), parameter :: MAX_STACK    = 128_I4P !< Octree traversal stack (7 * levels + 1 entries needed).

contains
   function winding_number(facet, tree, wn, point, beta) result(w)
   !< Generalized winding number of `point` with respect to the surface (`facet`, `tree`, `wn`).
   !<
   !< SAH BVH: the device-callable `bvh_winding_number_oac` (the same routine the GPU path uses).
   !< Octree: `octree_solid_angle`, the same expansion over up to 8 children. Without a tree or
   !< cache, or with beta <= 0, the exact per-facet sum.
   type(facet_object),     intent(in)           :: facet(:) !< Facets array (caller's surface).
   type(aabb_tree_object), intent(in)           :: tree     !< AABB tree built over `facet`.
   type(wn_node_t),        intent(in)           :: wn(0:)   !< Per-node cache built over `tree`.
   type(vector_R8P),       intent(in)           :: point    !< Query point.
   real(R8P),              intent(in), optional :: beta     !< Admissibility ratio (default 2.0; <= 0: exact).
   real(R8P)                                    :: w        !< Winding number.
   real(R8P)                                    :: beta_    !< Local copy of beta.

   beta_ = DEFAULT_BETA ; if (present(beta)) beta_ = beta
   if (tree%get_nodes_number() <= 0 .or. size(wn) /= tree%get_nodes_number() .or. beta_ <= 0._R8P) then
      w = exact_solid_angle_sum(facet=facet, point=point) / FOUR_PI
   elseif (tree%get_tree_kind() == AABB_TREE_SAH_BVH) then
      w = bvh_winding_number_oac(tree%node, tree%payload, wn, point%x, point%y, point%z, beta_ * beta_) / FOUR_PI
   else
      w = octree_solid_angle(facet=facet, tree=tree, wn=wn, point=point, beta2=beta_ * beta_) / FOUR_PI
   endif
   endfunction winding_number

   pure function exact_solid_angle_sum(facet, point) result(omega)
   !< Brute-force sum of signed solid angles over every facet: the fallback when no tree or cache
   !< is available, and the ground truth for the hierarchical paths.
   type(facet_object), intent(in) :: facet(:) !< Facets array.
   type(vector_R8P),   intent(in) :: point    !< Query point.
   real(R8P)                      :: omega    !< Sum of solid angles.
   integer(I4P)                   :: f        !< Facet counter.

   omega = 0._R8P
   do f = 1, size(facet, dim=1)
      omega = omega + facet(f)%solid_angle(point=point)
   enddo
   endfunction exact_solid_angle_sum

   function octree_solid_angle(facet, tree, wn, point, beta2) result(omega)
   !< Fast solid-angle sum over an octree: the expansion of `bvh_winding_number_oac` with up to 8
   !< children per node. Facets are summed at leaves only (the octree repeats each facet at every
   !< level, matching the leaf-only aggregation of the cache). CPU only.
   type(facet_object),     intent(in) :: facet(:)  !< Facets array.
   type(aabb_tree_object), intent(in) :: tree      !< Octree.
   type(wn_node_t),        intent(in) :: wn(0:)    !< Per-node cache.
   type(vector_R8P),       intent(in) :: point     !< Query point.
   real(R8P),              intent(in) :: beta2     !< Squared admissibility ratio (> 0).
   real(R8P)                          :: omega     !< Total solid angle.
   integer(I4P)                       :: stack(MAX_STACK), sp, cur, c, k, fid
   integer(I4P)                       :: child_idx(MAX_CHILDREN), nchild
   real(R8P)                          :: rx, ry, rz, r2, r1, rcr, trc

   omega    = 0._R8P
   stack(1) = 0_I4P
   sp       = 1_I4P
   do while (sp > 0_I4P)
      cur = stack(sp)
      sp  = sp - 1_I4P
      if (.not. wn(cur)%ready) cycle
      if (wn(cur)%diag2 > 0._R8P) then
         rx = point%x - wn(cur)%cx ; ry = point%y - wn(cur)%cy ; rz = point%z - wn(cur)%cz
         r2 = rx * rx + ry * ry + rz * rz
         if (r2 > beta2 * wn(cur)%diag2) then
            r1  = sqrt(r2)
            rcr = rx * (wn(cur)%c2(1,1) * rx + wn(cur)%c2(1,2) * ry + wn(cur)%c2(1,3) * rz) + &
                  ry * (wn(cur)%c2(2,1) * rx + wn(cur)%c2(2,2) * ry + wn(cur)%c2(2,3) * rz) + &
                  rz * (wn(cur)%c2(3,1) * rx + wn(cur)%c2(3,2) * ry + wn(cur)%c2(3,3) * rz)
            trc = wn(cur)%c2(1,1) + wn(cur)%c2(2,2) + wn(cur)%c2(3,3)
            omega = omega - (wn(cur)%dx * rx + wn(cur)%dy * ry + wn(cur)%dz * rz) / (r2 * r1) &
                          - (3._R8P * rcr - r2 * trc) / (r2 * r2 * r1)
            cycle
         endif
      endif
      call tree%enumerate_children(n=cur, out_idx=child_idx, nchild=nchild)
      if (nchild == 0) then
         if (tree%node(cur)%has_facets()) then
            do k = 1, tree%node(cur)%facet_id_count()
               fid = tree%node(cur)%facet_id_at(k)
               if (fid < 1 .or. fid > size(facet, dim=1)) cycle
               omega = omega + facet(fid)%solid_angle(point=point)
            enddo
         endif
      else
         do c = 1, nchild
            sp = sp + 1_I4P ; stack(sp) = child_idx(c)
         enddo
      endif
   enddo
   endfunction octree_solid_angle
endmodule fossil_winding_number
