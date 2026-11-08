!run figref-device-distance-batch fig_ref_device_distance_batch
!fig ref-device-distance-batch
!fig ref-device-distance-regions
!< Figures of the reference page "GPU distance queries", What it does and Pipeline. Built with gfortran, the device
!< loops run on the host: the same code and the same bits as on the GPU.
!< - one batch: the 20^3 cell centres of a block around the open dragon, signed by the winding number in one
!<   distance_many_device call; the centres of the cut cells (|d| within half a cell diagonal) drawn, teal on the inside
!<   of the surface, ember outside, the dragon translucent; the caption counts the bits that differ from distance_many;
!< - the sign step: query_many_device on a plane of 8 x 8 points over a corner of the unit cube, each joined to its
!<   closest point, coloured by the Voronoi region of that point (face teal, edge gold, vertex ember), and the
!<   pseudo-normal of the edge and vertex regions drawn at their closest points (red; a face's is its normal).
program fig_ref_device_distance_batch
use fossil,     only : surface_stl_object, SIGN_FAST_WINDING_NUMBER, SIGN_PSEUDO_NORMAL, FEATURE_FACE, FEATURE_EDGE, &
                       FEATURE_VERTEX, STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_segments, fig_points, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
integer(I4P), parameter       :: NC = 20, NB = 8
type(surface_stl_object)      :: body
type(vector_R8P)              :: bmin, bmax, h, q(NB * NB), foot(NB * NB), tip(NB * NB)
type(vector_R8P), allocatable :: centre(:)
real(R8P),        allocatable :: d(:), d_host(:)
real(R8P)                     :: half_diag, dq(NB * NB), cp(3, NB * NB), nr(3, NB * NB)
integer(I4P)                  :: feat(NB * NB), i, j, k, c, u, status
logical,          allocatable :: cut(:)
character(100)                :: caption

! ---- one batch over a block around the open dragon
call body%load_from_file(file_name='dragon.stl', guess_format=.true.)
bmin = body%get_bmin() ; bmax = body%get_bmax() ; h = bmax - bmin
bmin = bmin - 0.1_R8P * h ; bmax = bmax + 0.1_R8P * h
h = (bmax - bmin) / real(NC, R8P)
half_diag = 0.5_R8P * sqrt(h%x**2 + h%y**2 + h%z**2)
allocate(centre(NC**3), d(NC**3), d_host(NC**3))
c = 0
do k = 1, NC ; do j = 1, NC ; do i = 1, NC
   c = c + 1
   centre(c) = vector_R8P(bmin%x + (i - 0.5_R8P) * h%x, bmin%y + (j - 0.5_R8P) * h%y, bmin%z + (k - 0.5_R8P) * h%z)
enddo ; enddo ; enddo
call body%enter_device
call body%distance_many_device(points=centre, distances=d, is_signed=.true., is_square_root=.true., &
                               sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
call body%exit_device
if (status /= STATUS_OK) error stop 'distance_many_device failed'
call body%distance_many(points=centre, distances=d_host, is_signed=.true., is_square_root=.true., &
                        sign_algorithm=SIGN_FAST_WINDING_NUMBER)
cut = abs(d) <= half_diag
print '(A,I0,A,I0,A,I0,A,I0)', 'one call: ', NC**3, ' points, cut ', count(cut), ' (inside ', &
                               count(cut .and. d < 0._R8P), '), differing from distance_many ', count(d /= d_host)
write(caption, '(I0,A,I0,A,I0,A)') NC**3, ' points in one call; ', count(cut), ' cut cells shown; ', &
                                   count(d /= d_host), ' differ from distance_many'
call fig_open(u, 'ref-device-distance-batch', trim(caption), -60._R8P, 20._R8P)
call fig_surface(u, 'ghost', body)
call fig_points(u, 'inside', pack(centre, cut .and. d < 0._R8P))
call fig_points(u, 'outside', pack(centre, cut .and. d >= 0._R8P))
call fig_close(u)

! ---- the sign step: Voronoi regions and pseudo-normals around a corner of the unit cube
call body%load_from_file(file_name='cube-01.stl', guess_format=.true.)
c = 0
do j = 1, NB                                                       ! a plane of points above the corner (1, 1, 1)
   do i = 1, NB
      c = c + 1
      q(c) = vector_R8P(0.55_R8P + 0.9_R8P * (i - 1) / (NB - 1), 0.55_R8P + 0.9_R8P * (j - 1) / (NB - 1), 1.25_R8P)
   enddo
enddo
call body%enter_device
call body%query_many_device(points=q, distances=dq, closest_points=cp, features=feat, normals=nr, &
                            sign_algorithm=SIGN_PSEUDO_NORMAL, status=status)
call body%exit_device
if (status /= STATUS_OK) error stop 'query_many_device failed'
do c = 1, NB * NB
   foot(c) = vector_R8P(cp(1, c), cp(2, c), cp(3, c))
   tip(c) = foot(c) + 0.3_R8P * vector_R8P(nr(1, c), nr(2, c), nr(3, c))
enddo
print '(A,I0,A,I0,A,I0,A,I0)', 'query_many_device: ', NB * NB, ' points, closest in a face ', &
      count(feat == FEATURE_FACE), ', an edge ', count(feat == FEATURE_EDGE), ', the vertex ', &
      count(feat == FEATURE_VERTEX)
write(caption, '(A,I0,A,I0,A,I0,A)') 'closest point in a face: ', count(feat == FEATURE_FACE), ', edge: ', &
      count(feat == FEATURE_EDGE), ', vertex: ', count(feat == FEATURE_VERTEX), '; pseudo-normals red'
call fig_open(u, 'ref-device-distance-regions', trim(caption), -35._R8P, 30._R8P)
call fig_zoom(u, vector_R8P(0.9_R8P, 0.9_R8P, 0.95_R8P), 0.6_R8P)
call fig_surface(u, 'body', body)
call fig_segments(u, 'probe', pack(q, feat == FEATURE_FACE), pack(foot, feat == FEATURE_FACE))
call fig_segments(u, 'accent', pack(q, feat == FEATURE_EDGE), pack(foot, feat == FEATURE_EDGE))
call fig_segments(u, 'outside', pack(q, feat == FEATURE_VERTEX), pack(foot, feat == FEATURE_VERTEX))
call fig_segments(u, 'hole', pack(foot, feat /= FEATURE_FACE), pack(tip, feat /= FEATURE_FACE))
call fig_points(u, 'probe', q)
call fig_close(u)
endprogram fig_ref_device_distance_batch
