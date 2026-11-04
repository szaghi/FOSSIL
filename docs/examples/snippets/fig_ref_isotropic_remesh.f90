!< Figures of the reference page "Isotropic remeshing": the dragon as loaded, then after 3 iterations of
!< isotropic_remesh at the default target (the median input edge), seen from the same camera; the captions give the
!< facet count and the median and standard deviation of the facet edge lengths (each interior edge counted twice).
program fig_ref_isotropic_remesh
use fossil,     only : surface_stl_object, facet_object, REM_STATUS_OK
use fossil_fig, only : fig_open, fig_surface, fig_close
use penf,       only : I4P, R8P
use vecfor,     only : vector_R8P

implicit none
type(surface_stl_object) :: dragon
integer(I4P)             :: status

call dragon%load_from_file(file_name='dragon.stl', guess_format=.true.)
call scene('ref-remesh-dragon-input', 'dragon.stl')
call dragon%isotropic_remesh(iterations=3_I4P, status=status)
if (status /= REM_STATUS_OK) error stop 'remesh failed'
call scene('ref-remesh-dragon-output', 'isotropic_remesh, 3 iterations')
contains
   subroutine scene(name, title)
   !< The dragon with its edge-length statistics in the caption.
   character(*), intent(in)    :: name, title !< Figure name, caption head.
   type(facet_object), pointer :: f
   type(vector_R8P)            :: d
   real(R8P), allocatable      :: e(:)
   real(R8P)                   :: mean, sd
   character(100)              :: caption
   integer(I4P)                :: u, i, n

   n = dragon%get_facets_number()
   allocate(e(3 * n))
   do i = 1, n
      f => dragon%facet_at(i)
      d = f%vertex(2) - f%vertex(1) ; e(3 * i - 2) = d%normL2()
      d = f%vertex(3) - f%vertex(2) ; e(3 * i - 1) = d%normL2()
      d = f%vertex(1) - f%vertex(3) ; e(3 * i)     = d%normL2()
   enddo
   mean = sum(e) / size(e)
   sd = sqrt(sum((e - mean)**2) / size(e))
   write(caption, '(A,A,I0,A,F5.3,A,F5.3)') title, ': ', n, ' facets, edge median ', median(e), ', sd ', sd
   call fig_open(u, name, trim(caption), -60._R8P, 22._R8P)
   call fig_surface(u, 'body', dragon)
   call fig_close(u)
   endsubroutine scene

   function median(a) result(m)
   !< Median of `a` (partial selection sort up to the middle; a few thousand values).
   real(R8P), intent(in)  :: a(:) !< Values.
   real(R8P)              :: m    !< Median.
   real(R8P), allocatable :: b(:)
   real(R8P)              :: t
   integer(I4P)           :: i, k

   b = a
   do i = 1, size(b) / 2 + 1
      k = i - 1 + minloc(b(i:), dim=1)
      t = b(i) ; b(i) = b(k) ; b(k) = t
   enddo
   m = b(size(b) / 2 + 1)
   endfunction median
endprogram fig_ref_isotropic_remesh
