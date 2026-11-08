!< Figures of the reference page "Sign algorithms at sharp features", What it measures, Method and Device, on an
!< L-shaped prism (profile (0,0) (2,0) (2,1) (1,1) (1,2) (0,2), height 1: a reentrant edge and two saddle vertices):
!< - what is measured: points x + delta d at every vertex, edge midpoint and face centre of the prism, d the 8 diagonal
!<   directions (never in a face plane), delta = 0.04; coloured by their side by construction (teal inside, ember
!<   outside); the caption counts the points each sign algorithm puts on the wrong side;
!< - the construction at the reentrant edge, in the section z = 0.5: directions in the inside sector (270 degrees) and
!<   the outside sector (90), at least a sixth of the opening away from both faces (the excluded wedges drawn faint), and
!<   delta halving from 0.2 to 0.003 along each; coloured by the side the pseudo-normal gives;
!< - the device: the prism with one facet removed (a hole, in red), the points of the first figure signed by the fast
!<   winding number through query_many_device, against query_many on the host.
program fig_ref_sign_benchmark_points
use fossil,     only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                       SIGN_FAST_WINDING_NUMBER, STATUS_OK
use fossil_fig, only : fig_open, fig_zoom, fig_surface, fig_boundary, fig_segments, fig_points, fig_label, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
real(R8P),     parameter        :: PI = acos(-1._R8P), DELTA = 0.04_R8P
integer(I4P),  parameter        :: NV = 6, ALGO(4) = [SIGN_PSEUDO_NORMAL, SIGN_RAY_INTERSECTIONS, SIGN_SOLID_ANGLE, &
                                                      SIGN_FAST_WINDING_NUMBER]
type(surface_stl_object)        :: prism, holed
type(facet_object), allocatable :: f(:)
type(vector_R8P)                :: prof(NV), feat(64), q(512), sq(80), wb(4), centre
real(R8P)                       :: dq(512), dh(512), ang
logical                         :: truth(512), side(512), naive(512)
integer(I4P)                    :: wrong(5), i, j, k, m, n, nf, nq, ns, u, status
character(110)                  :: caption

prof = [vector_R8P(0._R8P, 0._R8P, 0._R8P), vector_R8P(2._R8P, 0._R8P, 0._R8P), vector_R8P(2._R8P, 1._R8P, 0._R8P), &
        vector_R8P(1._R8P, 1._R8P, 0._R8P), vector_R8P(1._R8P, 2._R8P, 0._R8P), vector_R8P(0._R8P, 2._R8P, 0._R8P)]
call make_prism(.false.)
call prism%adopt_facets(facets=f) ; call prism%analyze

! ---- what it measures: points at every feature, 8 diagonal directions, side known by construction
nf = 0
do k = 0, 1                                                          ! the 12 vertices
   do i = 1, NV
      nf = nf + 1 ; feat(nf) = prof(i) + real(k, R8P) * up()
   enddo
enddo
do i = 1, NV                                                         ! the 18 edge midpoints
   j = mod(i, NV) + 1
   nf = nf + 1 ; feat(nf) = 0.5_R8P * (prof(i) + prof(j))
   nf = nf + 1 ; feat(nf) = 0.5_R8P * (prof(i) + prof(j)) + up()
   nf = nf + 1 ; feat(nf) = prof(i) + 0.5_R8P * up()
enddo
do i = 1, NV                                                         ! the 8 face centres: 6 walls, 2 caps
   j = mod(i, NV) + 1
   nf = nf + 1 ; feat(nf) = 0.5_R8P * (prof(i) + prof(j)) + 0.5_R8P * up()
enddo
nf = nf + 1 ; feat(nf) = vector_R8P(0.5_R8P, 0.5_R8P, 0._R8P)
nf = nf + 1 ; feat(nf) = vector_R8P(0.5_R8P, 0.5_R8P, 1._R8P)
nq = 0
do i = 1, nf
   do k = 0, 7
      nq = nq + 1
      q(nq) = feat(i) + (DELTA / sqrt(3._R8P)) * vector_R8P(merge(1._R8P, -1._R8P, btest(k, 0)), &
                                       merge(1._R8P, -1._R8P, btest(k, 1)), merge(1._R8P, -1._R8P, btest(k, 2)))
      truth(nq) = in_l(q(nq))
   enddo
enddo
do m = 1, 4
   do i = 1, nq
      side(i) = prism%is_point_inside(point=q(i), sign_algorithm=ALGO(m))
   enddo
   wrong(m) = count(side(1:nq) .neqv. truth(1:nq))
enddo
do i = 1, nq
   naive(i) = naive_inside(prism, q(i))
enddo
wrong(5) = count(naive(1:nq) .neqv. truth(1:nq))
print '(A,I0,A,I0,A,I0,A)', 'L-prism: ', nf, ' feature points x 8 directions = ', nq, ' points, ', count(truth(1:nq)), &
                            ' inside by construction'
print '(A,5I5)', 'wrong: pseudo-normal, ray, solid angle, winding number, naive ', wrong
write(caption, '(I0,A,4(I0,A),I0)') nq, ' points; wrong: pseudo-normal ', wrong(1), ', ray ', wrong(2), &
                                    ', solid angle ', wrong(3), ', winding ', wrong(4), ', naive ', wrong(5)
call fig_open(u, 'ref-sign-benchmark-points', trim(caption), -125._R8P, 30._R8P)
call fig_surface(u, 'ghost', prism)
call fig_points(u, 'inside', pack(q(1:nq), truth(1:nq)))
call fig_points(u, 'outside', pack(q(1:nq), .not. truth(1:nq)))
call fig_close(u)

! ---- the construction at the reentrant edge (1, 1), section z = 0.5
centre = vector_R8P(1._R8P, 1._R8P, 0.5_R8P)
ns = 0
do i = 1, 10                                                         ! 3 outside directions, 7 inside
   if (i <= 3) then
      ang = (15._R8P + 30._R8P * (i - 1)) * PI / 180._R8P             ! outside: 15 .. 75 degrees
   else
      ang = (135._R8P + 30._R8P * (i - 4)) * PI / 180._R8P            ! inside: 135 .. 315 degrees
   endif
   do k = 0, 6
      ns = ns + 1
      sq(ns) = centre + (0.2_R8P * 0.5_R8P**k) * vector_R8P(cos(ang), sin(ang), 0._R8P)
   enddo
enddo
do i = 1, ns
   side(i) = prism%is_point_inside(point=sq(i), sign_algorithm=SIGN_PSEUDO_NORMAL)
   truth(i) = in_l(sq(i))
enddo
print '(A,I0,A,I0,A,I0)', 'reentrant edge: ', ns, ' points, ', count(truth(1:ns)), ' inside, pseudo-normal wrong ', &
                          count(side(1:ns) .neqv. truth(1:ns))
write(caption, '(A,I0,A,I0)') 'reentrant edge, z = 0.5: ', ns, ' points, delta 0.2 to 0.003; pseudo-normal wrong ', &
                                count(side(1:ns) .neqv. truth(1:ns))
wb = [centre + 0.25_R8P * dirxy(-45._R8P), centre + 0.25_R8P * dirxy(15._R8P), &
      centre + 0.25_R8P * dirxy(75._R8P), centre + 0.25_R8P * dirxy(135._R8P)]
call fig_open(u, 'ref-sign-benchmark-construction', trim(caption), -90._R8P, 90._R8P)
call fig_zoom(u, centre, 0.3_R8P)
call fig_segments(u, 'box', [centre + vector_R8P(0._R8P, 0.4_R8P, 0._R8P), centre], &           ! the section's faces
                  [centre, centre + vector_R8P(0.4_R8P, 0._R8P, 0._R8P)])
call fig_segments(u, 'box', [centre, centre, centre, centre], wb)
call fig_points(u, 'inside', pack(sq(1:ns), side(1:ns)))
call fig_points(u, 'outside', pack(sq(1:ns), .not. side(1:ns)))
call fig_label(u, 'outside', centre + 0.2_R8P * dirxy(45._R8P), 'outside sector, 90 deg')
call fig_label(u, 'inside', centre + 0.2_R8P * dirxy(225._R8P), 'inside sector, 270 deg')
call fig_close(u)

! ---- the device: the prism with a hole, winding-number signs on the device and on the host
call make_prism(.true.)
call holed%adopt_facets(facets=f) ; call holed%analyze
call holed%enter_device
call holed%query_many_device(points=q(1:nq), distances=dq(1:nq), sign_algorithm=SIGN_FAST_WINDING_NUMBER, status=status)
call holed%exit_device
if (status /= STATUS_OK) error stop 'query_many_device failed'
call holed%query_many(points=q(1:nq), distances=dh(1:nq), sign_algorithm=SIGN_FAST_WINDING_NUMBER)
side(1:nq) = dq(1:nq) < 0._R8P
print '(A,I0,A,I0,A,I0)', 'holed prism, winding number: device vs host mismatches ', count(dq(1:nq) /= dh(1:nq)), &
                          ' of ', nq, ', wrong side ', count(side(1:nq) .neqv. truth(1:nq))
write(caption, '(A,I0,A,I0,A,I0,A)') 'hole in red; device vs host: ', count(dq(1:nq) /= dh(1:nq)), ' of ', nq, &
                                     ' differ; ', count(side(1:nq) .neqv. truth(1:nq)), ' on the wrong side'
call fig_open(u, 'ref-sign-benchmark-device', trim(caption), -125._R8P, 30._R8P)
call fig_surface(u, 'ghost', holed)
call fig_boundary(u, 'hole', holed)
call fig_points(u, 'inside', pack(q(1:nq), side(1:nq)))
call fig_points(u, 'outside', pack(q(1:nq), .not. side(1:nq)))
call fig_close(u)
contains
   function up() result(v)
   !< The extrusion vector.
   type(vector_R8P) :: v

   v = vector_R8P(0._R8P, 0._R8P, 1._R8P)
   endfunction up

   function dirxy(deg) result(v)
   !< The unit vector at `deg` degrees in the xy plane.
   real(R8P), intent(in) :: deg
   type(vector_R8P)      :: v

   v = vector_R8P(cos(deg * PI / 180._R8P), sin(deg * PI / 180._R8P), 0._R8P)
   endfunction dirxy

   logical function in_l(p)
   !< Strictly inside the L-prism.
   type(vector_R8P), intent(in) :: p

   in_l = p%z > 0._R8P .and. p%z < 1._R8P .and. p%x > 0._R8P .and. p%y > 0._R8P .and. &
          ((p%x < 2._R8P .and. p%y < 1._R8P) .or. (p%x < 1._R8P .and. p%y < 2._R8P))
   endfunction in_l

   logical function naive_inside(s, p)
   !< The naive sign: the normal of the closest facet.
   type(surface_stl_object), intent(in) :: s
   type(vector_R8P),         intent(in) :: p
   type(facet_object), pointer          :: fp
   type(vector_R8P)                     :: c
   real(R8P)                            :: d
   integer(I4P)                         :: kf

   call s%compute_distance(point=p, distance=d, facet_index=kf, closest_point=c)
   fp => s%facet_at(kf)
   naive_inside = fp%normal%x * (p%x - c%x) + fp%normal%y * (p%y - c%y) + fp%normal%z * (p%z - c%z) < 0._R8P
   endfunction naive_inside

   subroutine make_prism(hole)
   !< The L-prism, outward: caps fanned from the profile's first vertex (the L is star-shaped from it), two facets per
   !< wall; `hole`: without the first wall facet.
   logical, intent(in) :: hole
   integer(I4P)        :: a, b

   if (allocated(f)) deallocate(f)
   allocate(f(4 * (NV - 2) + 2 * NV - merge(1, 0, hole)))
   m = 0
   do a = 2, NV - 1
      call tri(prof(1), prof(a + 1), prof(a))                        ! bottom, normal -z
      call tri(prof(1) + up(), prof(a) + up(), prof(a + 1) + up())   ! top, normal +z
   enddo
   do a = 1, NV
      b = mod(a, NV) + 1
      if (.not. (hole .and. a == 1)) call tri(prof(a), prof(b), prof(b) + up())
      call tri(prof(a), prof(b) + up(), prof(a) + up())
   enddo
   endsubroutine make_prism

   subroutine tri(p1, p2, p3)
   !< The facet p1 p2 p3 (counter-clockwise seen from outside), appended to f.
   type(vector_R8P), intent(in) :: p1, p2, p3

   m = m + 1
   f(m)%vertex(1) = p1 ; f(m)%vertex(2) = p2 ; f(m)%vertex(3) = p3
   call f(m)%compute_metrix
   endsubroutine tri
endprogram fig_ref_sign_benchmark_points
