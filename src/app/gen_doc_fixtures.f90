!< Documentation-figure fixture generator.
!<
!< Produces the STL and VTU files consumed by the documentation
!< image-rendering scripts to populate the "Visual reference" sections
!< of the smoothing, curvature, and cotangent-Laplacian advanced
!< pages.
!<
!< Outputs (in `docs/pictures/advanced/_fixtures/`):
!<
!<   smoothing-bunny-noisy.stl
!<     bunny.stl perturbed with deterministic per-vertex Gaussian
!<     noise. Renders as the "input" panel of the smoothing visual.
!<
!<   smoothing-bunny-taubin.stl
!<     The noisy bunny after 5 Taubin (lambda=0.5, mu=-0.53) pairs.
!<     Renders as the "after smoothing" panel.
!<
!<   curvature-sphere.vtu / curvature-sphere.stl
!<     Marching-cubes sphere (resolution 32^3) plus point-scalar fields
!<     `K` (Gaussian curvature) and `H` (signed mean curvature). The
!<     `.stl` is for reference / cross-checking; the `.vtu` carries the
!<     scalar fields the heatmap renderer reads.
!<
!<   laplacian-sphere.vtu
!<     Same sphere geometry with a point-scalar field `Lx = (L * x)`
!<     where `x` is the x-coordinate field. This is a deliberately
!<     non-trivial field (constant input `x=1` would give zero, which is
!<     the unit-test invariant — boring as a figure). `L * x` highlights
!<     where the discrete Laplacian "feels" the curvature of the
!<     sphere; the sign pattern reads cleanly as a heatmap.
!<
!< Determinism: every random offset uses a fixed seed and a simple
!< xorshift PRNG embedded in this file. The fixtures are byte-identical
!< across runs so the committed PNGs do not churn on rebuilds.

program gen_doc_fixtures

use fossil, only : surface_stl_object, extract_isosurface,                       &
                   csr_matrix_t,                                                  &
                   CURV_STATUS_OK, LAPL_STATUS_OK,                                &
                   SMOOTH_METHOD_TAUBIN
use fossil_facet_object, only : facet_object
use fossil_vertex_pool_object, only : vertex_pool_object
use penf, only : I4P, R8P, str
use vecfor, only : vector_R8P
use vtk_fortran, only : vtk_file

implicit none

character(len=*), parameter :: OUTDIR = 'docs/pictures/advanced/_fixtures/'

call ensure_outdir()
call gen_smoothing_pair()
call gen_curvature_fields()
call gen_laplacian_field()

print '(A)', 'gen_doc_fixtures: done.'

contains

   subroutine ensure_outdir()
   !< `mkdir -p` equivalent — VTKFortran's initializer does not create
   !< the directory itself, so we delegate to the shell.
   integer :: ios
   call execute_command_line('mkdir -p ' // OUTDIR, exitstat=ios)
   if (ios /= 0) error stop 'gen_doc_fixtures: could not create output directory'
   endsubroutine ensure_outdir

   subroutine gen_smoothing_pair()
   !< Load bunny.stl, perturb every pool vertex with a deterministic
   !< Gaussian offset of magnitude ~ 1% of the bounding-box diagonal,
   !< write that as `smoothing-bunny-noisy.stl`, then run Taubin and
   !< write `smoothing-bunny-taubin.stl`.
   type(surface_stl_object)            :: bunny
   type(surface_stl_object)            :: bunny_noisy
   type(facet_object), allocatable     :: facets(:)
   type(facet_object), pointer         :: src(:)
   type(vertex_pool_object), pointer   :: pool
   real(R8P)                           :: diag, sigma
   real(R8P)                           :: ox, oy, oz
   integer(I4P)                        :: nv, nf, vi, f, lv, status
   integer(I4P)                        :: seed
   type(vector_R8P), allocatable       :: offset(:)
   type(vector_R8P)                    :: bmin, bmax, v

   call bunny%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true., status=status)
   bmin = bunny%get_bmin(); bmax = bunny%get_bmax()
   v = bmax - bmin
   diag = sqrt(v%x**2 + v%y**2 + v%z**2)
   sigma = 0.01_R8P * diag

   pool => bunny%get_vertex_pool()
   nv = pool%vertex_count()
   allocate(offset(nv))

   seed = 1234567_I4P
   do vi = 1_I4P, nv
      ox = sigma * box_muller(seed)
      oy = sigma * box_muller(seed)
      oz = sigma * box_muller(seed)
      offset(vi) = vector_R8P(ox, oy, oz)
   enddo

   ! Deep-copy the facet array via the public pointer accessor, then
   ! perturb each cached vertex by the offset of its pool vertex id.
   ! Per-element copy (not whole-array) keeps gfortran from emitting a
   ! benign "type-bound defined assignment not done" notice.
   nf = bunny%get_facets_number()
   src => bunny%facets_ref()
   allocate(facets(nf))
   do f = 1_I4P, nf
      facets(f) = src(f)
      do lv = 1_I4P, 3_I4P
         vi = pool%facet_vid(facet_id=f, local_v=lv)
         if (vi == 0_I4P) cycle
         facets(f)%vertex(lv) = facets(f)%vertex(lv) + offset(vi)
      enddo
   enddo
   deallocate(offset)

   call bunny_noisy%adopt_facets(facets=facets)
   call bunny_noisy%save_into_file(file_name=OUTDIR // 'smoothing-bunny-noisy.stl', is_ascii=.true.)

   ! Smooth in place via Taubin (defaults). Save.
   call bunny_noisy%smooth(method=SMOOTH_METHOD_TAUBIN, status=status)
   call bunny_noisy%save_into_file(file_name=OUTDIR // 'smoothing-bunny-taubin.stl', is_ascii=.true.)

   print '(A,I0,A)', 'smoothing pair: bunny ', bunny%get_facets_number(), ' facets written'
   endsubroutine gen_smoothing_pair

   subroutine gen_curvature_fields()
   !< Build a unit-sphere surface and emit a VTU carrying point scalars
   !< `K` (Gaussian) and `H` (signed mean curvature). The matching
   !< `.stl` is written too for reference / cross-checking.
   !<
   !< The sphere is built by marching cubes then **isotropically
   !< remeshed** before curvature is computed. This is deliberate and
   !< matches the guidance on the curvature documentation page: raw MC
   !< output has many tessellation-degenerate vertices where the
   !< discrete angle-defect / mean-curvature formulae spike by orders
   !< of magnitude, and a figure dominated by those artefacts would
   !< hide the true near-uniform `K ~ 1`, `H ~ 1` signal of a sphere.
   !< Remeshing produces the well-conditioned triangulation that
   !< discrete curvature is designed to run on.
   type(surface_stl_object)         :: sphere
   real(R8P), allocatable           :: K(:), H(:)
   integer(I4P)                     :: status

   call build_remeshed_sphere(sphere)
   call sphere%save_into_file(file_name=OUTDIR // 'curvature-sphere.stl', is_ascii=.true.)

   call sphere%gaussian_curvature(K=K, status=status)
   if (status /= CURV_STATUS_OK) error stop 'gen_doc_fixtures: gaussian_curvature failed'
   call sphere%mean_curvature(H=H, status=status)
   if (status /= CURV_STATUS_OK) error stop 'gen_doc_fixtures: mean_curvature failed'

   call write_surface_vtu_two_scalars(surface=sphere,                              &
                                      file_name=OUTDIR // 'curvature-sphere.vtu',  &
                                      scalar1_name='K', scalar1=K,                 &
                                      scalar2_name='H', scalar2=H)
   print '(A,I0,A)', 'curvature: ', sphere%get_facets_number(), ' facets written'
   endsubroutine gen_curvature_fields

   subroutine gen_laplacian_field()
   !< Same remeshed sphere as gen_curvature_fields. Build (L, M) via
   !< §2.1, form Lx = L * (x-coordinate field), write VTU. The
   !< x-coordinate is a harmonic-on-the-plane field; on the curved
   !< sphere it is NOT harmonic, so `L * x` is non-zero and exhibits a
   !< clean dipole pattern aligned with +/- x — visually informative
   !< as a heatmap.
   type(surface_stl_object)          :: sphere
   type(vertex_pool_object), pointer :: pool
   type(csr_matrix_t)                :: L, M
   real(R8P), allocatable            :: x_field(:), Lx(:)
   type(vector_R8P)                  :: v
   integer(I4P)                      :: vi, n, status

   call build_remeshed_sphere(sphere)

   call sphere%cotangent_laplacian(L=L, M=M, status=status)
   if (status /= LAPL_STATUS_OK) error stop 'gen_doc_fixtures: cotangent_laplacian failed'

   pool => sphere%get_vertex_pool()
   n = pool%vertex_count()
   allocate(x_field(n), Lx(n))
   do vi = 1_I4P, n
      v = pool%coord(vid=vi)
      x_field(vi) = v%x
   enddo
   call L%multiply_vector(x=x_field, y=Lx)

   call write_surface_vtu_two_scalars(surface=sphere,                              &
                                      file_name=OUTDIR // 'laplacian-sphere.vtu',  &
                                      scalar1_name='x_field', scalar1=x_field,     &
                                      scalar2_name='Lx',       scalar2=Lx)
   print '(A,I0,A,I0,A)', 'laplacian: ', sphere%get_facets_number(),               &
                          ' facets, ', n, ' vertices written'
   endsubroutine gen_laplacian_field

   subroutine build_remeshed_sphere(sphere)
   !< Build a unit sphere by marching cubes, then isotropically remesh
   !< it into a well-conditioned triangulation suitable for discrete
   !< differential geometry. Shared by the curvature and Laplacian
   !< fixtures so both render the same geometry.
   !<
   !< The raw 32^3 MC output has many sliver triangles and near-
   !< degenerate vertices; 3 isotropic-remesh iterations at the median
   !< edge length regularise valences and edge lengths while the
   !< Botsch-Kobbelt projection pass keeps the vertices on the sphere.
   type(surface_stl_object), intent(out) :: sphere
   type(facet_object), allocatable       :: sphere_facets(:)
   real(R8P), allocatable                :: values(:,:,:)
   type(vector_R8P)                      :: bmin, bmax
   integer(I4P), parameter               :: N_GRID = 32_I4P
   integer(I4P)                          :: i, j, kk, status
   real(R8P)                             :: x, y, z, dx

   allocate(values(N_GRID, N_GRID, N_GRID))
   dx = 4._R8P / real(N_GRID - 1, R8P)
   do kk = 1_I4P, N_GRID
      z = -2._R8P + (kk - 1) * dx
      do j = 1_I4P, N_GRID
         y = -2._R8P + (j - 1) * dx
         do i = 1_I4P, N_GRID
            x = -2._R8P + (i - 1) * dx
            values(i, j, kk) = sqrt(x**2 + y**2 + z**2) - 1._R8P
         enddo
      enddo
   enddo
   bmin = vector_R8P(-2._R8P, -2._R8P, -2._R8P)
   bmax = vector_R8P( 2._R8P,  2._R8P,  2._R8P)
   call extract_isosurface(values=values, bmin=bmin, bmax=bmax, iso=0._R8P, &
                           surface=sphere_facets, status=status)
   call sphere%adopt_facets(facets=sphere_facets)
   deallocate(values)

   ! target_length <= 0 -> use the median input edge length.
   call sphere%isotropic_remesh(target_length=-1._R8P, iterations=3_I4P, &
                                preserve_features=.false., status=status)
   endsubroutine build_remeshed_sphere

   subroutine write_surface_vtu_two_scalars(surface, file_name, &
                                            scalar1_name, scalar1, &
                                            scalar2_name, scalar2)
   !< Emit an UnstructuredGrid VTU with two point-scalar fields.
   !< Connectivity comes from the surface's vertex pool (one triangle
   !< per facet, three pool vids per triangle).
   type(surface_stl_object), intent(in) :: surface
   character(len=*),         intent(in) :: file_name
   character(len=*),         intent(in) :: scalar1_name
   real(R8P),                intent(in) :: scalar1(:)
   character(len=*),         intent(in) :: scalar2_name
   real(R8P),                intent(in) :: scalar2(:)
   type(vertex_pool_object), pointer    :: pool
   type(vtk_file)                       :: vtu
   real(R8P),    allocatable            :: px(:), py(:), pz(:)
   integer(I4P), allocatable            :: connectivity(:), offsets(:)
   integer(1),   allocatable            :: cell_types(:)
   integer(I4P)                         :: n, nf, vi, f, error
   type(vector_R8P)                     :: v

   pool => surface%get_vertex_pool()
   n  = pool%vertex_count()
   nf = surface%get_facets_number()

   allocate(px(n), py(n), pz(n))
   do vi = 1_I4P, n
      v = pool%coord(vid=vi)
      px(vi) = v%x; py(vi) = v%y; pz(vi) = v%z
   enddo

   allocate(connectivity(3_I4P * nf))
   allocate(offsets(nf))
   allocate(cell_types(nf))
   do f = 1_I4P, nf
      ! VTKFortran connectivity is 0-based.
      connectivity(3_I4P*(f - 1_I4P) + 1_I4P) = pool%facet_vid(facet_id=f, local_v=1_I4P) - 1_I4P
      connectivity(3_I4P*(f - 1_I4P) + 2_I4P) = pool%facet_vid(facet_id=f, local_v=2_I4P) - 1_I4P
      connectivity(3_I4P*(f - 1_I4P) + 3_I4P) = pool%facet_vid(facet_id=f, local_v=3_I4P) - 1_I4P
      offsets(f) = 3_I4P * f
      cell_types(f) = int(5, 1)   ! VTK_TRIANGLE
   enddo

   error = vtu%initialize(format='raw', filename=trim(file_name), &
                          mesh_topology='UnstructuredGrid')
   error = vtu%xml_writer%write_piece(np=n, nc=nf)
   error = vtu%xml_writer%write_geo(np=n, nc=nf, x=px, y=py, z=pz)
   error = vtu%xml_writer%write_connectivity(nc=nf, connectivity=connectivity, &
                                             offset=offsets, cell_type=cell_types)
   error = vtu%xml_writer%write_dataarray(location='node', action='open')
   error = vtu%xml_writer%write_dataarray(data_name=trim(scalar1_name), x=scalar1)
   error = vtu%xml_writer%write_dataarray(data_name=trim(scalar2_name), x=scalar2)
   error = vtu%xml_writer%write_dataarray(location='node', action='close')
   error = vtu%xml_writer%write_piece()
   error = vtu%finalize()

   deallocate(px, py, pz, connectivity, offsets, cell_types)
   endsubroutine write_surface_vtu_two_scalars

   function box_muller(seed) result(z)
   !< Standard Gauss(0, 1) via Box-Muller from an embedded xorshift32.
   !< Deterministic given the input seed; advances the seed by two
   !< draws per call.
   integer(I4P), intent(inout) :: seed
   real(R8P)                    :: z
   real(R8P)                    :: u1, u2

   u1 = xorshift_uniform(seed)
   u2 = xorshift_uniform(seed)
   if (u1 < 1.e-300_R8P) u1 = 1.e-300_R8P
   z = sqrt(-2._R8P * log(u1)) * cos(2._R8P * acos(-1._R8P) * u2)
   endfunction box_muller

   function xorshift_uniform(seed) result(u)
   !< xorshift32 -> [0, 1) uniform. Advances `seed` in place.
   integer(I4P), intent(inout) :: seed
   real(R8P)                    :: u
   integer(I4P)                 :: s

   s = seed
   s = ieor(s, ishft(s,  13))
   s = ieor(s, ishft(s, -17))
   s = ieor(s, ishft(s,   5))
   seed = s
   u = real(iand(s, 2147483647_I4P), R8P) / 2147483648._R8P
   endfunction xorshift_uniform

endprogram gen_doc_fixtures
