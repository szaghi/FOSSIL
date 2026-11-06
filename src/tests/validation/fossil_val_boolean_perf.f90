!< FOSSIL validation: boolean run time against facet count (issue #63, roadmap #58 WP5, step 5).
!<
!< Three families, each boolean timed (best of 3) at -O2:
!< 1. a marching-cubes sphere (16^3 .. 96^3) minus a small capped cylinder through it: a large body, short cuts (the
!<    step of a CSG chain);
!< 2. the same sphere minus its copy moved by (0.3, 0.2, 0.1): two large bodies, a long intersection curve;
!< 3. the alpha-wrapped bunny against the alpha-wrapped dragon (alpha d/50, offset d/500, unit diagonal), the four ops.
!< Pass: every status OK and every result watertight; the times are printed (the scaling is in booleans.md).
program fossil_val_boolean_perf

use fossil,                only : surface_stl_object, BOOL_UNION, BOOL_INTERSECT, BOOL_DIFFERENCE, BOOL_SYMDIFF, &
                                  BOOL_STATUS_OK
use fossil_boolean_bodies, only : BODY_SPHERE, BODY_CYLINDER, mc_body
use penf,                  only : I4P, I8P, R8P
use vecfor,                only : vector_R8P

implicit none

integer(I4P), parameter  :: NS(6) = [16_I4P, 24_I4P, 32_I4P, 48_I4P, 64_I4P, 96_I4P]
character(9), parameter  :: OPNAME(4) = ['union    ', 'intersect', 'diff     ', 'symdiff  ']
type(surface_stl_object) :: a, b, c
real(R8P)                :: t
integer(I4P)             :: i, k, nf
logical                  :: ok

ok = .true.
call mc_body(BODY_CYLINDER, 12_I4P, c)
call c%resize(factor=vector_R8P(0.3_R8P, 0.3_R8P, 1.5_R8P), respect_centroid=.true.)
call c%translate(x=0.1_R8P) ; call c%analyze
print '(A)', 'sphere - small cylinder    facets(A+B)   time [s]   us/facet'
do i = 1, size(NS)
   call mc_body(BODY_SPHERE, NS(i), a)
   call timed(a, c, BOOL_DIFFERENCE, t, nf)
   print '(A,I3,A,I10,F11.4,F11.2)', '   ', NS(i), '^3                ', nf, t, 1e6_R8P * t / nf
enddo
print '(A)', 'sphere - moved copy        facets(A+B)   time [s]   us/facet'
do i = 1, size(NS)
   call mc_body(BODY_SPHERE, NS(i), a) ; b = a
   call b%translate(x=0.3_R8P, y=0.2_R8P, z=0.1_R8P) ; call b%analyze
   call timed(a, b, BOOL_DIFFERENCE, t, nf)
   print '(A,I3,A,I10,F11.4,F11.2)', '   ', NS(i), '^3                ', nf, t, 1e6_R8P * t / nf
enddo
call wrapped('src/tests/bunny.stl', a)
call wrapped('src/tests/dragon.stl', b)
call b%rotate(axis=vector_R8P(0.2_R8P, 1._R8P, -0.3_R8P), angle=0.7_R8P) ; call b%translate(x=0.1_R8P) ; call b%analyze
print '(A)', 'wrapped bunny x dragon     facets(A+B)   time [s]   us/facet'
do k = 1, 4
   call timed(a, b, k, t, nf)
   print '(3A,I10,F11.4,F11.2)', '   ', OPNAME(k), '              ', nf, t, 1e6_R8P * t / nf
enddo
print '(A,L1)', 'Are all tests passed? ', ok

contains
   subroutine timed(a, b, op, t, nf)
   !< Best of three runs of op(a, b); status and watertightness folded into ok.
   type(surface_stl_object), intent(in)  :: a, b !< Inputs.
   integer(I4P),             intent(in)  :: op   !< BOOL_*.
   real(R8P),                intent(out) :: t    !< Best time [s].
   integer(I4P),             intent(out) :: nf   !< Input facets.
   type(surface_stl_object)              :: r    !< Result.
   integer(I8P)                          :: t0, t1, rate
   integer(I4P)                          :: rep, st

   t = huge(1._R8P) ; nf = a%get_facets_number() + b%get_facets_number()
   do rep = 1, 3
      r = a
      call system_clock(t0, rate)
      call r%boolean(other=b, op=op, status=st)
      call system_clock(t1)
      t = min(t, real(t1 - t0, R8P) / rate)
   enddo
   ok = ok .and. st == BOOL_STATUS_OK
   if (op /= BOOL_SYMDIFF) ok = ok .and. r%is_watertight()
   endsubroutine timed

   subroutine wrapped(file_name, s)
   !< The alpha wrap of a mesh, scaled to a unit diagonal and centred at the origin.
   character(*),             intent(in)    :: file_name !< STL file.
   type(surface_stl_object), intent(inout) :: s         !< Wrapped surface.
   type(surface_stl_object)                :: raw       !< Input mesh.
   type(vector_R8P)                        :: lo, hi    !< Box.
   real(R8P)                               :: d         !< Diagonal.
   integer(I4P)                            :: st        !< Status.

   call raw%load_from_file(file_name=file_name, guess_format=.true.)
   lo = raw%get_bmin() ; hi = raw%get_bmax() ; d = norm2([hi%x - lo%x, hi%y - lo%y, hi%z - lo%z])
   call raw%alpha_wrap(alpha=d / 50._R8P, offset=d / 500._R8P, wrapped=s, status=st)
   if (.not. s%is_watertight()) error stop 'alpha wrap not closed'
   lo = s%get_bmin() ; hi = s%get_bmax()
   call s%translate(x=-0.5_R8P * (lo%x + hi%x), y=-0.5_R8P * (lo%y + hi%y), z=-0.5_R8P * (lo%z + hi%z))
   call s%resize(factor=vector_R8P(1._R8P / d, 1._R8P / d, 1._R8P / d))
   call s%analyze
   endsubroutine wrapped
endprogram fossil_val_boolean_perf
