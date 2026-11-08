!< Figures of the reference page "Surface deviation", What it does, Method and On the device:
!< - the two one-sided deviations of a small cube A (side 0.2) inside a unit cube B, both translucent: h(A, B) is A's
!<   farthest vertex to B, the gap; h(B, A) is B's farthest vertex to A, a far corner; each drawn as a segment from
!<   the vertex (per_vertex of deviation) to its closest point on the other surface;
!< - the samples of the lower bound on an octahedron A inscribed in a unit sphere B (marching cubes): every vertex, every
!<   edge midpoint and the centroids of the 16 sub-triangles of each facet (samples_per_area = 16 / facet area), coloured
!<   by their distance to B (query_many: teal below half the largest, ember above); the farthest sample, the lower
!<   bound, in gold;
!< - the device: the dragon turned by 2 degrees about its centre, against the dragon, with device=.true. (the dragon on
!<   the device; with gfortran the device loops run on the host), each facet coloured by its vertices' mean distance
!<   (teal 0, cream 0.05, ember 0.1 and above); the caption
!<   compares the bounds with the host run.
program fig_ref_surface_deviation_method
use fossil,                    only : surface_stl_object, facet_object, extract_isosurface, MC_STATUS_OK, DEV_STATUS_OK
use fossil_vertex_pool_object, only : vertex_pool_object
use fossil_fig,                only : fig_open, fig_surface, fig_segment, fig_points, fig_label, fig_triangle, fig_close
use penf,                      only : I4P, R8P
use vecfor,                    only : vector_R8P

implicit none
integer(I4P), parameter           :: M = 31, NS = 4
type(surface_stl_object)          :: a, b, sphere, octa, rotated
type(facet_object), allocatable   :: f(:)
type(facet_object), pointer       :: fp
type(vertex_pool_object), pointer :: pool
type(vector_R8P)                  :: x, c, ov(6), smp(8 * (3 + NS * NS)), cen
type(vector_R8P), allocatable     :: verts(:)
real(R8P),        allocatable     :: pv(:), pv2(:), ds(:)
real(R8P)                         :: phi(M, M, M), h, lab, uab, lba, uba, lower, low_h, up_h, low_d, up_d, dmax
integer(I4P)                      :: i, j, k, n, u, status, ia, ib, imax, oct(3, 8)
character(100)                    :: caption
character(40)                     :: text

! ---- what it does: one-sided deviations of a small cube inside a large one
call b%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call a%load_from_file(file_name='cube-01.stl', guess_format=.true.)
call a%resize(x=0.2_R8P, y=0.2_R8P, z=0.2_R8P)
call a%translate(delta=vector_R8P(0.05_R8P, 0.4_R8P, 0.4_R8P))
call a%analyze
call a%deviation(b, lower=lab, upper=uab, per_vertex=pv)
call b%deviation(a, lower=lba, upper=uba, per_vertex=pv2)
print '(A,2F8.4,A,2F8.4)', 'h(A, B) in ', lab, uab, ';  h(B, A) in ', lba, uba
call fig_open(u, 'ref-deviation-one-sided', 'h(A, B) = '//trim(num(uab))//' (gold), h(B, A) = '//trim(num(uba))// &
              ' (teal): not symmetric', -60._R8P, 25._R8P)
call fig_surface(u, 'ghost', b)
call fig_surface(u, 'ghost', a)
pool => a%get_vertex_pool() ; ia = maxloc(pv, dim=1) ; x = pool%coord(ia)
c = foot(b, x)
call fig_segment(u, 'accent', x, c) ; call fig_points(u, 'foot', [x, c])
call fig_label(u, 'accent', 0.5_R8P * (x + c), 'h(A, B)')
pool => b%get_vertex_pool() ; ib = maxloc(pv2, dim=1) ; x = pool%coord(ib)
c = foot(a, x)
call fig_segment(u, 'probe', x, c) ; call fig_points(u, 'probe', [x, c])
call fig_label(u, 'probe', 0.5_R8P * (x + c), 'h(B, A)')
call fig_label(u, 'body', a%get_bmax(), 'A')
call fig_label(u, 'ghost', b%get_bmax(), 'B')
call fig_close(u)

! ---- method: the samples on an octahedron inscribed in a sphere
h = 3._R8P / real(M - 1, R8P)
do k = 1, M ; do j = 1, M ; do i = 1, M
   phi(i, j, k) = norm2(-1.5_R8P + h * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=f, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
call sphere%adopt_facets(facets=f)
ov = [vector_R8P(1._R8P, 0._R8P, 0._R8P), vector_R8P(-1._R8P, 0._R8P, 0._R8P), vector_R8P(0._R8P, 1._R8P, 0._R8P), &
      vector_R8P(0._R8P, -1._R8P, 0._R8P), vector_R8P(0._R8P, 0._R8P, 1._R8P), vector_R8P(0._R8P, 0._R8P, -1._R8P)]
oct = reshape([1, 3, 5,  3, 2, 5,  2, 4, 5,  4, 1, 5,  3, 1, 6,  2, 3, 6,  4, 2, 6,  1, 4, 6], [3, 8])
if (allocated(f)) deallocate(f)
allocate(f(8))
do i = 1, 8
   f(i)%vertex(1) = ov(oct(1, i)) ; f(i)%vertex(2) = ov(oct(2, i)) ; f(i)%vertex(3) = ov(oct(3, i))
   call f(i)%compute_metrix
enddo
call octa%adopt_facets(facets=f) ; call octa%analyze
fp => octa%facet_at(1)
call octa%deviation(sphere, lower=lower, samples_per_area=real(NS * NS, R8P) / fp%area(), samples=n)
k = 0                                                               ! the same samples, here: mid-edges, sub-centroids
do i = 1, 8
   fp => octa%facet_at(i)
   k = k + 1 ; smp(k) = 0.5_R8P * (fp%vertex(1) + fp%vertex(2))
   k = k + 1 ; smp(k) = 0.5_R8P * (fp%vertex(2) + fp%vertex(3))
   k = k + 1 ; smp(k) = 0.5_R8P * (fp%vertex(3) + fp%vertex(1))
   do j = 0, NS - 1
      do n = 0, NS - 1 - j
         k = k + 1 ; smp(k) = lattice(fp, 3._R8P * j + 1._R8P, 3._R8P * n + 1._R8P)
         if (j + n <= NS - 2) then
            k = k + 1 ; smp(k) = lattice(fp, 3._R8P * j + 2._R8P, 3._R8P * n + 2._R8P)
         endif
      enddo
   enddo
enddo
verts = [ov, smp(1:k)]
allocate(ds(size(verts)))
call sphere%query_many(points=verts, distances=ds, is_signed=.false.)
imax = maxloc(ds, dim=1) ; dmax = ds(imax)
print '(A,I0,A,F8.5,A,F8.5)', 'octahedron: ', size(verts), ' samples, largest distance ', dmax, &
                              ', deviation lower ', lower
x = verts(imax)                                                     ! look at the farthest sample
call fig_open(u, 'ref-deviation-samples', trim(num(real(size(verts), R8P), '(I0)'))//' samples, coloured by d; '// &
              'the largest, '//trim(num(dmax))//' = lower bound (gold)', &
              atan2(x%y, x%x) * 180._R8P / acos(-1._R8P) - 20._R8P, asin(x%z / sqrt(x%x**2 + x%y**2 + x%z**2)) * &
              180._R8P / acos(-1._R8P) + 10._R8P)
call fig_surface(u, 'body', octa)
call fig_points(u, 'probe', pack(verts, ds < 0.5_R8P * dmax))
call fig_points(u, 'outside', pack(verts, ds >= 0.5_R8P * dmax .and. ds < dmax))
call fig_points(u, 'foot', [verts(imax)])
call fig_close(u)

! ---- on the device: the dragon turned by 2 degrees, against the dragon
call a%load_from_file(file_name='dragon.stl', guess_format=.true.)
rotated = a
cen = 0.5_R8P * (a%get_bmin() + a%get_bmax())
call rotated%rotate(axis=vector_R8P(0._R8P, 0._R8P, 1._R8P), angle=2._R8P * acos(-1._R8P) / 180._R8P, center=cen)
call rotated%analyze
call rotated%deviation(a, lower=low_h, upper=up_h, per_vertex=pv)
call a%enter_device
call rotated%deviation(a, lower=low_d, upper=up_d, per_vertex=pv2, device=.true., status=status)
call a%exit_device
if (status /= DEV_STATUS_OK) error stop 'device deviation failed'
print '(A,2F9.5,A,2F9.5)', 'turned dragon: host ', low_h, up_h, ', device ', low_d, up_d
print '(A,L1,A,I0)', 'bounds bit-identical ', low_h == low_d .and. up_h == up_d, ', per-vertex mismatches ', &
                     count(pv /= pv2)
write(caption, '(A,F6.4,A,F6.4,A,I0,A)') 'device: ', low_d, ' <= h <= ', up_d, ', equal to the host; ', &
                                         count(pv /= pv2), ' vertex values differ'
call fig_open(u, 'ref-deviation-device', trim(caption), -40._R8P, 20._R8P)
pool => rotated%get_vertex_pool()
do i = 1, rotated%get_facets_number()
   fp => rotated%facet_at(i)
   call fig_triangle(u, 'field', fp%vertex(1), fp%vertex(2), fp%vertex(3), &
                     value=min(1._R8P, 2._R8P * sum(pv2(fp%vertex_id)) / (3._R8P * 0.1_R8P) - 1._R8P))
enddo
call fig_close(u)
contains
   function foot(s, p) result(q)
   !< The closest point of s to p.
   type(surface_stl_object), intent(in) :: s
   type(vector_R8P),         intent(in) :: p
   type(vector_R8P)                     :: q
   real(R8P)                            :: d

   call s%compute_distance(point=p, distance=d, closest_point=q)
   endfunction foot

   function lattice(ft, s, t) result(q)
   !< The point a + (b - a) s / (3 n) + (c - a) t / (3 n) of facet ft (n = NS): a sub-triangle centroid.
   type(facet_object), intent(in) :: ft
   real(R8P),          intent(in) :: s, t
   type(vector_R8P)               :: q

   q = ft%vertex(1) + (s / (3._R8P * NS)) * (ft%vertex(2) - ft%vertex(1)) + (t / (3._R8P * NS)) * (ft%vertex(3) - ft%vertex(1))
   endfunction lattice

   function num(v, fmt) result(txt)
   !< A number as text.
   real(R8P),    intent(in)           :: v
   character(*), intent(in), optional :: fmt
   character(24)                      :: txt

   if (present(fmt)) then
      write(txt, fmt) nint(v)
   else
      write(txt, '(F6.4)') v
   endif
   txt = adjustl(txt)
   endfunction num
endprogram fig_ref_surface_deviation_method
