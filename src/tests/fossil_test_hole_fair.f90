!< FOSSIL test: fairing of hole patches on a curved surface (issue #31, step 5).
program fossil_test_hole_fair
!< FOSSIL test: fairing of hole patches on a curved surface (issue #31, step 5).
!<
!< A unit icosphere built here (icosahedron subdivided twice: 320 facets, tiny) with the cap
!< z > 0.75 removed: one curved hole. Filled three ways (plain, refine, fair):
!<  - every result is watertight (the loop vertices are kept bit for bit, or the weld would fail);
!<  - the faired patch follows the sphere: max |r - 1| at the patch triangles' centroids is lower
!<    than for the plain (flat) patch;
!<  - the faired patch has a lower maximum dihedral angle (between the normals of facets sharing
!<    an edge, over the patch's facets) than the plain patch: it continues the surface with
!<    matching tangents;
!<  - the enclosed volume gets closer to the sphere's polyhedron volume than the plain fill.
!< And a planar hole along two sharp edges: cube-02 with a 3 x 3 block of squares removed from
!< the face z = 1 at the corner x = y = 0. The loop runs along the cube's edges, where the loop
!< edges belong to the side faces; fairing must keep the 90-degree creases (feature-aware ring),
!< so the faired patch stays in the plane z = 1 and the volume stays 1 (to 1e-12).
!< Issue #60: the certified deviation of the faired sphere from the full icosphere is below that of the plain fill.

use fossil, only : surface_stl_object, facet_object
use penf,   only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(facet_object), allocatable :: sphere(:)
real(R8P)                       :: rdev(3), dih(3), vol(3), vol_full
real(R8P)                       :: hl(3), hu(3) !< Deviation of the filled sphere from the full one (issue #60).
integer(I4P)                    :: n_patch(3)
logical                         :: ok, wt(3)
type(surface_stl_object)        :: full

ok = .true.
call icosphere(levels=2_I4P, f=sphere)
call full%adopt_facets(facets=sphere)
vol_full = full%get_volume()
call run(mode=1_I4P, rdev=rdev(1), dih=dih(1), vol=vol(1), wt=wt(1), n_patch=n_patch(1), hl=hl(1), hu=hu(1))
call run(mode=2_I4P, rdev=rdev(2), dih=dih(2), vol=vol(2), wt=wt(2), n_patch=n_patch(2), hl=hl(2), hu=hu(2))
call run(mode=3_I4P, rdev=rdev(3), dih=dih(3), vol=vol(3), wt=wt(3), n_patch=n_patch(3), hl=hl(3), hu=hu(3))
print '(A,F10.6)', 'closed icosphere volume: ', vol_full
print '(A)', '            patch  watertight  max|r-1|   max dihedral (deg)  volume'
print '(A,I6,L8,F14.6,F14.3,F16.6)', 'plain   ', n_patch(1), wt(1), rdev(1), dih(1), vol(1)
print '(A,I6,L8,F14.6,F14.3,F16.6)', 'refine  ', n_patch(2), wt(2), rdev(2), dih(2), vol(2)
print '(A,I6,L8,F14.6,F14.3,F16.6)', 'fair    ', n_patch(3), wt(3), rdev(3), dih(3), vol(3)
ok = all(wt) .and. rdev(3) < rdev(1) .and. dih(3) < dih(1) .and. abs(vol(3) - vol_full) < abs(vol(1) - vol_full)
! Issue #60: the one-sided deviation of the filled sphere from the full icosphere, certified (default tolerance). The
! faired patch follows the sphere: its upper bound is below the plain patch's lower bound.
print '(A,2ES11.3,A,2ES11.3,A,2ES11.3)', 'deviation from the full sphere: plain ', hl(1), hu(1), '  refine ', hl(2), hu(2), &
                                         '  fair ', hl(3), hu(3)
ok = ok .and. hu(3) < hl(1)
call crease_case
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine crease_case
   !< Planar hole along two cube edges, faired: must stay flat.
   type(surface_stl_object)        :: cube, holed
   type(facet_object), allocatable :: f(:)
   type(facet_object), pointer     :: p
   integer(I4P)                    :: i, n, filled, n_orig
   logical                         :: flat, case_ok

   call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)
   allocate(f(cube%get_facets_number())) ; n = 0
   do i = 1, cube%get_facets_number()
      p => cube%facet_at(i)
      if (p%centroid%z > 1._R8P - 1.e-9_R8P .and. p%centroid%x < 0.75_R8P .and. p%centroid%y < 0.75_R8P) cycle
      n = n + 1 ; f(n) = p ; f(n)%vertex_id = 0_I4P
   enddo
   f = f(1:n) ; n_orig = n
   call holed%adopt_facets(facets=f)
   call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, fair=.true.)
   flat = .true.
   do i = n_orig + 1, holed%get_facets_number()
      p => holed%facet_at(i)
      flat = flat .and. all([p%vertex(1)%z, p%vertex(2)%z, p%vertex(3)%z] == 1._R8P)
   enddo
   case_ok = filled == 1 .and. holed%is_watertight() .and. abs(holed%get_volume() - 1._R8P) < 1.e-12_R8P .and. flat
   print '(A,I0,A,L1,A,F16.12,A,L1)', 'crease hole, faired: patch ', holed%get_facets_number() - n_orig, '  flat ', flat, &
         '  volume ', holed%get_volume(), '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine crease_case

   subroutine run(mode, rdev, dih, vol, wt, n_patch, hl, hu)
   integer(I4P), intent(in)        :: mode
   real(R8P),    intent(out)       :: rdev, dih, vol
   real(R8P),    intent(out)       :: hl, hu
   logical,      intent(out)       :: wt
   integer(I4P), intent(out)       :: n_patch
   type(surface_stl_object)        :: holed
   type(facet_object), allocatable :: f(:)
   type(facet_object), pointer     :: p, q
   integer(I4P)                    :: i, e, n_orig, filled, nb
   real(R8P)                       :: c

   call icosphere(levels=2_I4P, f=f)
   f = pack(f, f(:)%centroid%z <= 0.75_R8P)
   n_orig = size(f)
   call holed%adopt_facets(facets=f)
   select case (mode)
   case (1) ; call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled)
   case (2) ; call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, refine=.true.)
   case (3) ; call holed%fill_holes(max_hole_size=1._R8P, holes_filled=filled, fair=.true.)
   endselect
   n_patch = holed%get_facets_number() - n_orig
   call holed%deviation(full, lower=hl, upper=hu)
   wt  = holed%is_watertight() .and. filled == 1
   vol = holed%get_volume()
   rdev = 0._R8P ; dih = 0._R8P
   do i = n_orig + 1, holed%get_facets_number()
      p => holed%facet_at(i)
      c = sqrt(p%centroid%x**2 + p%centroid%y**2 + p%centroid%z**2)
      rdev = max(rdev, abs(c - 1._R8P))
      do e = 1, 3
         nb = p%fcon_edge(e)
         if (nb <= 0) cycle
         q => holed%facet_at(nb)
         dih = max(dih, acos(max(-1._R8P, min(1._R8P, p%normal%x * q%normal%x + p%normal%y * q%normal%y + &
                                                      p%normal%z * q%normal%z))) * 180._R8P / acos(-1._R8P))
      enddo
   enddo
   endsubroutine run

   subroutine icosphere(levels, f)
   !< Unit icosphere: icosahedron subdivided `levels` times, midpoints projected on the sphere
   !< ((a + b) / |a + b|, commutative, so a midpoint shared by two triangles has the same bits).
   integer(I4P),                    intent(in)  :: levels
   type(facet_object), allocatable, intent(out) :: f(:)
   type(vector_R8P), allocatable :: tri(:,:), nxt(:,:)
   type(vector_R8P)              :: v(12), m12, m23, m31
   integer(I4P), parameter       :: IDX(3, 20) = reshape([1, 12, 6,  1, 6, 2,  1, 2, 8,  1, 8, 11,  1, 11, 12, &
                                                        2, 6, 10,  6, 12, 5,  12, 11, 3,  11, 8, 7,  8, 2, 9,  &
                                                        4, 10, 5,  4, 5, 3,  4, 3, 7,  4, 7, 9,  4, 9, 10, &
                                                        5, 10, 6,  3, 5, 12,  7, 3, 11,  9, 7, 8,  10, 9, 2], [3, 20])
   real(R8P)                     :: g
   integer(I4P)                  :: k, l, n

   g = (1._R8P + sqrt(5._R8P)) / 2._R8P
   v = [vec(-1._R8P, g, 0._R8P), vec(1._R8P, g, 0._R8P), vec(-1._R8P, -g, 0._R8P), vec(1._R8P, -g, 0._R8P), &
        vec(0._R8P, -1._R8P, g), vec(0._R8P, 1._R8P, g), vec(0._R8P, -1._R8P, -g), vec(0._R8P, 1._R8P, -g), &
        vec(g, 0._R8P, -1._R8P), vec(g, 0._R8P, 1._R8P), vec(-g, 0._R8P, -1._R8P), vec(-g, 0._R8P, 1._R8P)]
   do k = 1, 12
      v(k) = unit(v(k))
   enddo
   allocate(tri(3, 20))
   do k = 1, 20
      tri(:, k) = v(IDX(:, k))
   enddo
   do l = 1, levels
      n = size(tri, dim=2)
      allocate(nxt(3, 4 * n))
      do k = 1, n
         m12 = unit(tri(1, k) + tri(2, k)) ; m23 = unit(tri(2, k) + tri(3, k)) ; m31 = unit(tri(3, k) + tri(1, k))
         nxt(:, 4 * k - 3) = [tri(1, k), m12, m31]
         nxt(:, 4 * k - 2) = [m12, tri(2, k), m23]
         nxt(:, 4 * k - 1) = [m31, m23, tri(3, k)]
         nxt(:, 4 * k)     = [m12, m23, m31]
      enddo
      call move_alloc(from=nxt, to=tri)
   enddo
   n = size(tri, dim=2)
   allocate(f(n))
   do k = 1, n
      f(k)%vertex = tri(:, k)
      call f(k)%compute_metrix
   enddo
   endsubroutine icosphere

   pure function vec(x, y, z) result(r)
   real(R8P), intent(in) :: x, y, z
   type(vector_R8P)      :: r

   r%x = x ; r%y = y ; r%z = z
   endfunction vec

   pure function unit(a) result(r)
   type(vector_R8P), intent(in) :: a
   type(vector_R8P)             :: r
   real(R8P)                    :: n

   n = sqrt(a%x**2 + a%y**2 + a%z**2)
   r%x = a%x / n ; r%y = a%y / n ; r%z = a%z / n
   endfunction unit
endprogram fossil_test_hole_fair
