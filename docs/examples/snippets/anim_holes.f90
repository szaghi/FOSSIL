!< Animation of tutorial 2: fill_holes closing holes one after the other. A sphere (marching cubes on |x| - 1, as the
!< curvature recipe builds it; the raw extraction has slits, so sanitize(fill_holes) repairs it, then it is decimated
!< to 2200 facets) gets four round holes of different sizes (the facets whose centroid is within r_k of a point of the
!< surface). fill_holes(refine, fair) then runs with a max_hole_size that grows from frame to frame, set between the
!< loop sizes (printed), so each frame closes the next larger hole: the patches in ember, the open loops outlined.
program anim_holes
use fossil,              only : surface_stl_object, extract_isosurface, MC_STATUS_OK, DEC_STATUS_OK
use fossil_facet_object, only : facet_object
use fossil_fig,          only : fig_open, fig_anim, fig_frame, fig_surface, fig_boundary, fig_close
use penf,                only : I4P, R8P
use vecfor,              only : vector_R8P

implicit none
integer(I4P),       parameter   :: M = 31, NH = 4                          ! samples per direction, holes
real(R8P),          parameter   :: DEG = 57.29577951308232_R8P, AZ = -60._R8P, EL = 25._R8P
real(R8P),          parameter   :: RADIUS(NH) = [0.13_R8P, 0.24_R8P, 0.36_R8P, 0.50_R8P]    ! hole radii
real(R8P),          parameter   :: SIZES(NH) = [0.15_R8P, 0.25_R8P, 0.35_R8P, 0.45_R8P]     ! max_hole_size per frame
real(R8P),          parameter   :: OFFSET(2, NH) = reshape([0.62_R8P, -0.55_R8P, -0.62_R8P, -0.55_R8P, &
                                                            0.62_R8P, 0.55_R8P, -0.62_R8P, 0.55_R8P], [2, NH])
type(surface_stl_object)        :: sphere
type(facet_object), allocatable :: facets(:), kept(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: centre(NH), view, right, up
real(R8P)                       :: phi(M, M, M), h, a, e, sizes_cut(NH)
logical                         :: cut
character(80)                   :: caption
integer(I4P)                    :: u, i, j, k, n, n0, filled, skipped, status

h = 3._R8P / (M - 1)                                                        ! the sphere of radius 1 on [-1.5, 1.5]^3
do k = 1, M ; do j = 1, M ; do i = 1, M
   phi(i, j, k) = norm2(-1.5_R8P + h * ([i, j, k] - 1)) - 1._R8P
enddo ; enddo ; enddo
call extract_isosurface(values=phi, bmin=vector_R8P(-1.5_R8P, -1.5_R8P, -1.5_R8P), &
                        bmax=vector_R8P(1.5_R8P, 1.5_R8P, 1.5_R8P), iso=0._R8P, surface=facets, status=status)
if (status /= MC_STATUS_OK) error stop 'marching cubes failed'
a = AZ / DEG ; e = EL / DEG                                                 ! the hole centres: on the side facing the
view = vector_R8P(cos(e) * cos(a), cos(e) * sin(a), sin(e))                 ! camera, two by two around its axis
right = vector_R8P(-sin(a), cos(a), 0._R8P)
up = view .cross. right
do k = 1, NH
   centre(k) = view + OFFSET(1, k) * right + OFFSET(2, k) * up
   centre(k) = centre(k) / norm2([centre(k)%x, centre(k)%y, centre(k)%z])
enddo
call sphere%adopt_facets(facets=facets)
call sphere%sanitize(fill_holes=.true.)                                     ! the raw extraction has slits
call sphere%decimate(target_facets=2200_I4P, status=status)                 ! (and a light drawing)
if (status /= DEC_STATUS_OK) error stop 'decimate failed'
print '(A,I0,A,L1)', 'sphere: ', sphere%get_facets_number(), ' facets, watertight ', sphere%is_watertight()
n = 0
allocate(kept(sphere%get_facets_number()))
do i = 1, sphere%get_facets_number()
   f => sphere%facet_at(i)
   cut = .false.
   do k = 1, NH
      cut = cut .or. norm2([f%centroid%x - centre(k)%x, f%centroid%y - centre(k)%y, f%centroid%z - centre(k)%z]) &
                     <= RADIUS(k)
   enddo
   if (cut) cycle
   n = n + 1 ; kept(n) = f
enddo
if (allocated(facets)) deallocate(facets)
allocate(facets(n))
do i = 1, n
   facets(i) = kept(i)
enddo
call sphere%adopt_facets(facets=facets)
n0 = sphere%get_facets_number()
do k = 1, NH                                                                ! each loop's bbox diagonal / the sphere's
   sizes_cut(k) = loop_size(centre(k), RADIUS(k) + 0.15_R8P)
   print '(A,I0,A,F6.3)', 'hole ', k, ': loop size ', sizes_cut(k)
enddo

call fig_open(u, 'anim-holes', '', AZ, EL)
call fig_anim(u, 1.2_R8P)
write(caption, '(I0,A,4F5.2)') NH, ' holes, loop size / sphere size:', sizes_cut
call fig_frame(u, caption=trim(caption))
call fig_surface(u, 'body', sphere)
call fig_boundary(u, 'hole', sphere)
do k = 1, NH
   n = sphere%get_facets_number()
   call sphere%fill_holes(max_hole_size=SIZES(k), refine=.true., fair=.true., holes_filled=filled, &
                          holes_skipped=skipped)
   print '(A,F5.2,A,I0,A,I0,A,I0,A)', 'max_hole_size ', SIZES(k), ': filled ', filled, ', left open ', skipped, &
                                      ', patch of ', sphere%get_facets_number() - n, ' facets'
   write(caption, '(A,F4.2,A,I0,A,I0,A)') 'fill_holes(max_hole_size=', SIZES(k), '): ', filled, ' filled, ', &
                                         skipped, ' left open'
   if (k < NH) then
      call fig_frame(u, caption=trim(caption))
   else
      call fig_frame(u, caption=trim(caption), seconds=2.5_R8P)
   endif
   call fig_surface(u, 'body', sphere, last=n0)
   call fig_surface(u, 'patch', sphere, first=n0 + 1_I4P)
   call fig_boundary(u, 'hole', sphere)
enddo
call fig_close(u)
print '(A,L1)', 'watertight: ', sphere%is_watertight()
contains
   function loop_size(c, reach) result(ratio)
   !< Bounding-box diagonal of the boundary vertices within `reach` of c, over the sphere's: what max_hole_size tests.
   type(vector_R8P), intent(in) :: c     !< Hole centre.
   real(R8P),        intent(in) :: reach !< Search radius.
   real(R8P)                    :: ratio !< Ratio.
   type(vector_R8P)             :: bmin, bmax
   real(R8P)                    :: lo(3), hi(3), p(3)
   integer(I4P)                 :: ii, ee

   lo = huge(1._R8P) ; hi = -huge(1._R8P)
   do ii = 1, sphere%get_facets_number()
      f => sphere%facet_at(ii)
      do ee = 1, 3
         if (f%fcon_edge(ee) /= 0_I4P) cycle
         p = [f%vertex(ee)%x, f%vertex(ee)%y, f%vertex(ee)%z]
         if (norm2(p - [c%x, c%y, c%z]) > reach) cycle
         lo = min(lo, p) ; hi = max(hi, p)
      enddo
   enddo
   bmin = sphere%get_bmin() ; bmax = sphere%get_bmax()
   ratio = norm2(hi - lo) / norm2([bmax%x - bmin%x, bmax%y - bmin%y, bmax%z - bmin%z])
   endfunction loop_size
endprogram anim_holes
