!< FOSSIL test: signs outside the band of `distance_field_object` (issue #32, step 3).
program fossil_test_distance_field_sign
!< FOSSIL test: signs outside the band of `distance_field_object` (issue #32, step 3).
!<
!< The sign of every interior cell must equal the sign of a per-cell pseudo-normal
!< `surface%distance`, on:
!<  - cube-02 (unit cube) on a 32^3 grid around it;
!<  - a hollow cube: cube-02 plus an inward-facing inner cube [0.25, 0.75]^3, so the cavity
!<    cells are outside the solid;
!<  - the hollow cube with a band thinner than the grid spacing: the fallback must run
!<    (fallback cells > 0) and the signs must still match;
!<  - grids that do not see the surface (all outside, all inside): one point query each.

use fossil, only : surface_stl_object, facet_object, SIGN_PSEUDO_NORMAL, STATUS_OK
use fossil_distance_field, only : distance_field_object
use penf, only : I4P, R8P
use vecfor, only : vector_R8P

implicit none

type(surface_stl_object)        :: cube    !< Unit cube.
type(surface_stl_object)        :: hollow  !< Unit cube with a cavity.
type(facet_object), allocatable :: facets(:)
type(facet_object), pointer     :: f
type(vector_R8P)                :: v
real(R8P),          allocatable :: xg(:)   !< Grid axis.
logical                         :: ok
integer(I4P)                    :: i, nf, n

ok = .true.
call cube%load_from_file(file_name='src/tests/cube-02.stl', guess_format=.true.)

! hollow cube: outer cube-02 plus the inner cube [0.25, 0.75]^3 with reversed winding
nf = cube%get_facets_number()
allocate(facets(2 * nf))
do i = 1, nf
   f => cube%facet_at(i)
   facets(i) = f
   facets(i)%vertex_id = 0_I4P
   facets(nf + i) = f
   facets(nf + i)%vertex_id = 0_I4P
   facets(nf + i)%vertex(1) = 0.25_R8P + 0.5_R8P * f%vertex(1)
   v                        = 0.25_R8P + 0.5_R8P * f%vertex(2)
   facets(nf + i)%vertex(2) = 0.25_R8P + 0.5_R8P * f%vertex(3)
   facets(nf + i)%vertex(3) = v
   call facets(nf + i)%compute_normal
enddo
call hollow%adopt_facets(facets=facets)
print '(A,F10.6,A,L1)', 'hollow cube volume (expect 0.875): ', hollow%get_volume(), '  watertight: ', hollow%is_watertight()
ok = ok .and. abs(hollow%get_volume() - 0.875_R8P) < 1.e-12_R8P

n = 32
allocate(xg(n))
do i = 1, n
   xg(i) = -0.3_R8P + 1.6_R8P * real(i - 1, R8P) / real(n - 1, R8P)
enddo
call check_case(label='cube-02', surface=cube, x=xg, expect_fallback=.false.)
call check_case(label='hollow cube', surface=hollow, x=xg, expect_fallback=.false.)
call check_case(label='hollow cube, thin band', surface=hollow, x=xg, band=0.3_R8P * (xg(2) - xg(1)), &
                expect_fallback=.true.)
call check_case(label='grid outside the surface', surface=cube, x=3._R8P + xg, expect_fallback=.true.)
call check_case(label='grid inside the solid', surface=cube, x=0.4_R8P + 0.1_R8P * (xg + 0.3_R8P) / 1.6_R8P, &
                band=0.05_R8P, expect_fallback=.true.)

print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   subroutine check_case(label, surface, x, expect_fallback, band)
   !< Band + signs on an x = y = z grid; compare every interior sign with the per-cell PN sign.
   character(*),             intent(in)           :: label           !< Case label.
   type(surface_stl_object), intent(inout)        :: surface         !< Surface.
   real(R8P),                intent(in)           :: x(:)            !< Axis.
   logical,                  intent(in)           :: expect_fallback !< Fallback expected.
   real(R8P),                intent(in), optional :: band            !< Explicit band.
   type(distance_field_object)                    :: field           !< Field.
   type(vector_R8P)                               :: p               !< Point.
   real(R8P)                                      :: d               !< Reference distance.
   integer(I4P)                                   :: a, b, c, m, st1, st2, st3, n_fb, n_bad, n_in
   logical                                        :: case_ok

   m = size(x, kind=I4P)
   call field%initialize(x=x, y=x, z=x, status=st1)
   call field%compute_band(surface=surface, band=band, sign_algorithm=SIGN_PSEUDO_NORMAL, status=st2)
   call field%compute_signs(surface=surface, fallback_cells=n_fb, status=st3)
   n_bad = 0 ; n_in = 0
   do c = 1, m
      do b = 1, m
         do a = 1, m
            p%x = x(a) ; p%y = x(b) ; p%z = x(c)
            d = surface%distance(point=p, is_signed=.true., sign_algorithm=SIGN_PSEUDO_NORMAL, is_square_root=.true.)
            if ((d < 0._R8P) .neqv. (field%phi(a, b, c) < 0._R8P)) n_bad = n_bad + 1
            if (field%phi(a, b, c) < 0._R8P) n_in = n_in + 1
         enddo
      enddo
   enddo
   case_ok = st1 == STATUS_OK .and. st2 == STATUS_OK .and. st3 == STATUS_OK .and. n_bad == 0 .and. &
             ((n_fb > 0) .eqv. expect_fallback)
   print '(A,A,I0,A,I0,A,I0,A,I0,A,L1)', label, ': band cells ', field%get_band_cells(), '  inside ', n_in, &
         '  fallback ', n_fb, '  sign mismatches ', n_bad, '  -> ', case_ok
   ok = ok .and. case_ok
   endsubroutine check_case
endprogram fossil_test_distance_field_sign
