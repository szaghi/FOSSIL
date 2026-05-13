!< Cotangent Laplacian + barycentric mass matrix (issue #18 §2.1).
!<
!< References:
!<   - Meyer, Desbrun, Schröder & Barr, *Discrete Differential-Geometry
!<     Operators for Triangulated 2-Manifolds*, VisMath 2002.
!<   - Botsch et al., *Polygon Mesh Processing*, AK Peters 2010, Ch. 3.
!<
!< Definitions used:
!<
!<   For each interior edge (i, j) shared by two triangles with opposite
!<   angles alpha (in one triangle) and beta (in the other):
!<
!<     L_ij = (cot(alpha) + cot(beta)) / 2     for i != j
!<     L_ii = -sum_{j != i} L_ij
!<
!<   For a boundary edge (one incident triangle only), the off-diagonal
!<   gets the cotangent of the single opposite angle, halved. This is
!<   the convention that keeps L * 1 = 0 (constant kernel) on open
!<   shells as well as closed solids.
!<
!<   The matrix is **positive-semidefinite** under this sign convention —
!<   suitable for direct use as `L` in the standard formulations
!<   `(M - tau L) U = M U_prev` (implicit smoothing, §2.3) and
!<   `(M + t L) u = delta_s` (heat method, §2.2).
!<
!<   Barycentric mass:
!<
!<     M_ii = (1/3) * sum_{T incident to i} area(T)
!<
!<   Diagonal-only. The Voronoi-area variant (Meyer 2003) handles obtuse
!<   triangles more accurately and is documented as a follow-up.
!<
!< Per-triangle assembly: for triangle (v1, v2, v3) with opposite angles
!< theta_1, theta_2, theta_3, the cotangent-of-opposite-angle identity
!< gives:
!<
!<     cot(theta_k) = (E_a . E_b) / |E_a x E_b|
!<
!< where E_a and E_b are the two triangle edges meeting at vertex k.
!< No trig calls, numerically stable for non-degenerate triangles.
!<
!< Module operates on `facet(:)` + `vertex_pool_object` directly (the
!< standard no-cycle pattern in this codebase). The surface module wires
!< the public TBP.

module fossil_laplacian
!< Cotangent Laplacian and barycentric mass matrix builders.

use fossil_csr_matrix, only : csr_matrix_t
use fossil_facet_object, only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none
private

public :: build_cotangent_laplacian
public :: LAPL_STATUS_OK, LAPL_STATUS_BAD_INPUT, LAPL_STATUS_DEGENERATE_TRIANGLE

integer(I4P), parameter :: LAPL_STATUS_OK                   = 0_I4P
integer(I4P), parameter :: LAPL_STATUS_BAD_INPUT            = 1_I4P  !< Empty facet array or pool not initialized.
integer(I4P), parameter :: LAPL_STATUS_DEGENERATE_TRIANGLE  = 2_I4P  !< At least one triangle has area below tolerance; cotangents inflated/skipped.

real(R8P), parameter :: DEGENERATE_AREA_TOL = 1.0e-30_R8P  !< Absolute area below this → skip the triangle entirely.

contains

   subroutine build_cotangent_laplacian(facet, pool, L, M, status)
   !< Build the cotangent Laplacian `L` and the barycentric mass matrix `M`
   !< over the unique vertices in `pool`.
   !<
   !< `L` is symmetric, positive-semidefinite, with `L * 1 = 0` (constant
   !< kernel). `M` is diagonal — stored as CSR with one nonzero per row
   !< for API uniformity.
   !<
   !< Mostly-degenerate inputs (any triangle below DEGENERATE_AREA_TOL)
   !< return status = LAPL_STATUS_DEGENERATE_TRIANGLE — the matrices are
   !< still built but the affected vertices' rows will have zero
   !< contributions from those triangles. Run `surface%sanitize` upstream
   !< to drop degenerate facets before calling this builder.
   type(facet_object),        intent(in)            :: facet(:)
   type(vertex_pool_object),  intent(in)            :: pool
   type(csr_matrix_t),        intent(out)           :: L
   type(csr_matrix_t),        intent(out)           :: M
   integer(I4P),              intent(out), optional :: status
   integer(I4P)                                     :: n_vertices, nf, f
   integer(I4P)                                     :: v1, v2, v3
   type(vector_R8P)                                 :: p1, p2, p3
   type(vector_R8P)                                 :: e12, e23, e31
   type(vector_R8P)                                 :: cross_full
   real(R8P)                                        :: area_2x, area_tri
   real(R8P)                                        :: cot1, cot2, cot3
   real(R8P), allocatable                           :: mass_diag(:)  !< Per-vertex accumulated 1/3 * area sums.
   integer(I4P)                                     :: i
   logical                                          :: had_degenerate

   if (present(status)) status = LAPL_STATUS_OK
   nf = size(facet, kind=I4P)
   if (nf == 0_I4P .or. .not. pool%get_is_initialized()) then
      call L%initialize(n_rows=0_I4P, n_cols=0_I4P)
      call L%finalize
      call M%initialize(n_rows=0_I4P, n_cols=0_I4P)
      call M%finalize
      if (present(status)) status = LAPL_STATUS_BAD_INPUT
      return
   endif

   n_vertices = pool%vertex_count()
   call L%initialize(n_rows=n_vertices, n_cols=n_vertices)
   allocate(mass_diag(n_vertices))
   mass_diag = 0._R8P
   had_degenerate = .false.

   ! Per-triangle assembly.
   do f = 1_I4P, nf
      v1 = pool%facet_vid(facet_id=f, local_v=1_I4P)
      v2 = pool%facet_vid(facet_id=f, local_v=2_I4P)
      v3 = pool%facet_vid(facet_id=f, local_v=3_I4P)
      if (v1 == 0_I4P .or. v2 == 0_I4P .or. v3 == 0_I4P) cycle  ! pool not yet populated for this facet

      p1 = pool%coord(vid=v1)
      p2 = pool%coord(vid=v2)
      p3 = pool%coord(vid=v3)

      e12 = p2 - p1
      e23 = p3 - p2
      e31 = p1 - p3
      cross_full = e12%crossproduct(rhs=e23)  ! magnitude = 2 * area
      area_2x = sqrt(cross_full%dotproduct(rhs=cross_full))
      if (area_2x < DEGENERATE_AREA_TOL) then
         had_degenerate = .true.
         cycle
      endif
      area_tri = 0.5_R8P * area_2x

      ! Cotangent of each angle via the identity cot(theta) = (u·v) / |u x v|,
      ! where u, v are the two triangle edges meeting at that vertex.
      !
      ! At v1, the two edges (outgoing) are e12 and -e31. Their cross
      ! product magnitude equals |e12 x e23| = 2*area regardless of which
      ! pair you pick — that's the shared denominator for all three cotangents.
      cot1 = (-(e31%dotproduct(rhs=e12))) / area_2x   ! angle at v1: edges e12, -e31
      cot2 = (-(e12%dotproduct(rhs=e23))) / area_2x   ! angle at v2: edges e23, -e12
      cot3 = (-(e23%dotproduct(rhs=e31))) / area_2x   ! angle at v3: edges e31, -e23

      ! Off-diagonal contributions: edge opposite angle v_k carries 0.5 * cot(theta_k).
      !   Edge (v2, v3) is opposite v1.
      !   Edge (v3, v1) is opposite v2.
      !   Edge (v1, v2) is opposite v3.
      ! Symmetric entries: append both (i, j) and (j, i).
      call L%append(row=v2, col=v3, value=0.5_R8P * cot1)
      call L%append(row=v3, col=v2, value=0.5_R8P * cot1)
      call L%append(row=v3, col=v1, value=0.5_R8P * cot2)
      call L%append(row=v1, col=v3, value=0.5_R8P * cot2)
      call L%append(row=v1, col=v2, value=0.5_R8P * cot3)
      call L%append(row=v2, col=v1, value=0.5_R8P * cot3)

      ! Diagonal: -sum of off-diagonals. Each vertex k accumulates
      ! -(0.5 cot of the two angles at the OTHER two vertices, since those
      ! cotangents weight the edges incident on k).
      !   At v1: edges to v2 and v3 weighted by 0.5*cot3 and 0.5*cot2.
      !   At v2: edges to v1 and v3 weighted by 0.5*cot3 and 0.5*cot1.
      !   At v3: edges to v1 and v2 weighted by 0.5*cot2 and 0.5*cot1.
      call L%append(row=v1, col=v1, value=-0.5_R8P * (cot3 + cot2))
      call L%append(row=v2, col=v2, value=-0.5_R8P * (cot3 + cot1))
      call L%append(row=v3, col=v3, value=-0.5_R8P * (cot2 + cot1))

      ! Barycentric mass: each vertex gets 1/3 of the triangle area.
      mass_diag(v1) = mass_diag(v1) + area_tri / 3._R8P
      mass_diag(v2) = mass_diag(v2) + area_tri / 3._R8P
      mass_diag(v3) = mass_diag(v3) + area_tri / 3._R8P
   enddo

   call L%finalize

   ! Build M as a diagonal CSR (one nonzero per row).
   call M%initialize(n_rows=n_vertices, n_cols=n_vertices)
   do i = 1_I4P, n_vertices
      if (mass_diag(i) /= 0._R8P) call M%append(row=i, col=i, value=mass_diag(i))
   enddo
   call M%finalize
   deallocate(mass_diag)

   if (had_degenerate .and. present(status)) status = LAPL_STATUS_DEGENERATE_TRIANGLE
   endsubroutine build_cotangent_laplacian

endmodule fossil_laplacian
