!< FOSSIL validation: ray against medial thickness and gap on the Stanford bunny (issue #62, roadmap #58 WP4, step 2).
!<
!< Thickness and gap by both methods (defaults: 60-degree cone, 31 rays; shrinking ball with a 30-degree separation angle,
!< or the angle given as the first argument); prints the facets with a value, the quantiles of ray / medial over the
!< facets where both are finite, and the times. The two measure different things: the ray method the chord along the
!< normal, the medial method the local feature size (on a blobby body, curvature bounds the inscribed ball well below the
!< chord). Pass: no facet invalid for the medial method, and the medial thickness never above the axial chord (rays=1),
!< its starting bound. The bunny is not watertight (401 facets have too many escaping rays).
program fossil_val_thickness

use fossil, only : surface_stl_object, THICKNESS_RAY, THICKNESS_MEDIAL
use penf,   only : I4P, I8P, R8P

implicit none

type(surface_stl_object) :: s
real(R8P), allocatable   :: tr(:), gr(:), tm(:), gm(:), ta(:), q(:)
integer(I4P)             :: ir, im
real(R8P)                :: sepa
character(16)            :: arg
integer(I8P)             :: c0, c1, c2, rate
logical                  :: ok

call s%load_from_file(file_name='src/tests/bunny.stl', guess_format=.true.)
sepa = 30._R8P ; if (command_argument_count() > 0) then ; call get_command_argument(1, arg) ; read(arg, *) sepa ; endif
call system_clock(c0, rate)
call s%thickness(thickness=tr, gap=gr, method=THICKNESS_RAY, invalid=ir)
call system_clock(c1)
call s%thickness(thickness=tm, gap=gm, method=THICKNESS_MEDIAL, invalid=im, separation_angle=sepa)
call system_clock(c2)
call s%thickness(thickness=ta, rays=1_I4P)
print '(A,I0,A,L1)', 'bunny: ', s%get_facets_number(), ' facets, watertight ', s%is_watertight()
print '(A,F7.2,A,F7.2,A)', 'time ray ', real(c1 - c0, R8P) / rate, ' s, medial ', real(c2 - c1, R8P) / rate, ' s'
print '(A,I0,A,I0)', 'invalid facets: ray ', ir, ', medial ', im
print '(A,I0,A,I0)', 'facets with a gap: ray ', count(gr < huge(1._R8P)), ', medial ', count(gm < huge(1._R8P))
q = quantiles(pack(tr / tm, tr < huge(1._R8P) .and. tm < huge(1._R8P)))
print '(A,5F8.3)', 'thickness ray / medial, quantiles 5 25 50 75 95 %: ', q
ok = im == 0 .and. all(tm <= ta * (1._R8P + 1e-12_R8P) .or. ta == huge(1._R8P))
print '(A,L1)', 'medial <= axial chord on every facet: ', all(tm <= ta * (1._R8P + 1e-12_R8P) .or. ta == huge(1._R8P))
if (count(gr < huge(1._R8P) .and. gm < huge(1._R8P)) > 0) then
   print '(A,5F8.3)', 'gap ray / medial, quantiles 5 25 50 75 95 %:       ', &
         quantiles(pack(gr / gm, gr < huge(1._R8P) .and. gm < huge(1._R8P)))
endif
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) error stop 1

contains
   function quantiles(x) result(qq)
   !< The 5, 25, 50, 75, 95 % quantiles of x (heap-sorted copy, nearest rank).
   real(R8P), intent(in) :: x(:)
   real(R8P)             :: qq(5)
   real(R8P)             :: w(size(x)), key
   integer(I4P)          :: i, n

   n = size(x) ; qq = 0._R8P
   if (n == 0) return
   w = x
   do i = n / 2, 1, -1                                 ! heap sort (n ~ 1e5)
      call sift(w, i, n)
   enddo
   do i = n, 2, -1
      key = w(1) ; w(1) = w(i) ; w(i) = key
      call sift(w, 1, i - 1)
   enddo
   qq = w(max(1, min(n, nint([0.05_R8P, 0.25_R8P, 0.5_R8P, 0.75_R8P, 0.95_R8P] * n))))
   endfunction quantiles

      subroutine sift(w, lo, hi)
      !< Sift w(lo) down the max-heap w(lo:hi).
      real(R8P),    intent(inout) :: w(:)
      integer(I4P), intent(in) :: lo, hi
      integer(I4P)             :: r, c
      real(R8P)                :: t

      r = lo
      do while (2 * r <= hi)
         c = 2 * r
         if (c < hi) then
            if (w(c + 1) > w(c)) c = c + 1
         endif
         if (w(r) >= w(c)) exit
         t = w(r) ; w(r) = w(c) ; w(c) = t ; r = c
      enddo
      endsubroutine sift
endprogram fossil_val_thickness
