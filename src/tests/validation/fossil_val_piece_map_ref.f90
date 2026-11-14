!< FOSSIL validation: reference hashes of the piece map (issue #105 step 1a).
!<
!< The step-0 matrix of #105: the ellipse (1024 segments) and the NACA 0012 (1022) on the device-chain grid at 256^2 and
!< 1024^2; the marching-cubes sphere (24^3 samples) and the centred cone of the sign bench on the box of
!< `fossil_val_piece_map_3d` at 32^3 and 64^3. pose0 is rotated by 0.1 rad. Motions (the step-0 probe's): a translation
!< (exact, `translation`), and a rotation about (0.2, 1, 0.3) (2-D: the body centre) of 0.75 delta h / r_max plus a
!< translation of 0.25 delta h, integrated exactly (`rotation`), by Gauss with nq = 4 (`gauss4`) and by Gauss with
!< `refine` (`gauss4+refine`); delta = 2 and 0.5 cells per step (the largest vertex displacement). 64 cases.
!<
!< One line per case: its name, the map's counters (key=value) and a 64-bit hash (16 hex digits) of the bits of every
!< array field of `piece_map_2d_t` / `piece_map_3d_t` (allocation flag, shape, then the elements: integers as 64-bit
!< words, reals by `transfer`; the scheme of `fossil_val_perf`'s hash). Equal lines mean equal maps, field for field.
!<
!< Usage (from the repository root):
!<   fossil_val_piece_map_ref ref FILE [quick] [case=N ...] [note=TEXT]   write the lines to FILE (refused if it exists)
!<   fossil_val_piece_map_ref check FILE [quick] [case=N ...]              recompute, fail on any difference (printed
!<                                                                         field by field) or a case missing in FILE
!<   fossil_val_piece_map_ref                                              = check perf/map-ref-f855194.txt quick
!< `quick` runs the 256^2 and 32^3 cases only (32 cases), `case=N` the N-th case of the matrix (repeatable). The
!< reference `src/tests/validation/perf/map-ref-f855194.txt` was generated from the tree of f855194 (before any #105
!< change); it is never overwritten: a deliberate change of the map writes a new named file.
module fossil_val_piece_map_ref_kit
!< Cases, hash and lines of the piece-map reference.
use fossil,                   only : surface_stl_object, curve_object, rigid_pose_t, rigid_pose_2d_t, &
                                     pose_from_axis_angle, pose_2d, swept_cells_2d_t, swept_cells_3d_t, swept_cells_2d, &
                                     swept_cells_3d, piece_map_2d_t, piece_map_3d_t
use fossil_boolean_bodies,    only : BODY_SPHERE, mc_body
use fossil_sign_bench_bodies, only : body_t, bench_bodies, to_surface
use penf,                     only : I4P, I8P, R8P
use vecfor,                   only : vector_R8P

implicit none
private
public :: NCASES, MAXLEN, case_name, case_quick, case_line, setup_bodies, split_tokens

integer(I4P), parameter :: NCASES = 64_I4P      !< Cases of the matrix.
integer(I4P), parameter :: MAXLEN = 2048_I4P    !< Line length.
integer(I8P), parameter :: P1 = 2147483647_I8P, P2 = 2147483629_I8P !< Primes < 2^31 (`fossil_val_perf`).
integer(I8P), parameter :: B1 = 1000003_I8P, B2 = 65599_I8P          !< Multipliers.
real(R8P),    parameter :: PI = acos(-1._R8P)   !< Pi.
character(16), parameter :: PNAME(4) = [character(16) :: 'translation', 'rotation', 'gauss4', 'gauss4+refine'] !< Paths.
character(8),  parameter :: BNAME(2, 2) = reshape([character(8) :: 'ellipse', 'naca0012', 'sphere', 'cone'], [2, 2])
                                                !< Bodies (body, dimension - 1).

type :: hash_t
   !< Two polynomial hashes mod primes below 2^31 over the 16-bit chunks of 64-bit words (no integer overflow).
   integer(I8P) :: a = 1_I8P !< Hash mod P1.
   integer(I8P) :: b = 7_I8P !< Hash mod P2.
endtype hash_t

interface field
   module procedure field_i1, field_i2, field_i3, field_r1, field_r2
endinterface field

type(curve_object)       :: curves(2)   !< Ellipse, NACA 0012.
type(surface_stl_object) :: surfaces(2) !< Sphere, cone (centred).
real(R8P)                :: half(2)     !< Half side of the 3-D box of each body.

contains
   subroutine case_params(ic, dim, ib, n, path, delta)
   !< Parameters of case ic (order: dimension, body, size, path, delta).
   integer(I4P), intent(in)  :: ic    !< Case (1..NCASES).
   integer(I4P), intent(out) :: dim   !< 2 or 3.
   integer(I4P), intent(out) :: ib    !< Body 1 or 2.
   integer(I4P), intent(out) :: n     !< Cells per direction.
   integer(I4P), intent(out) :: path  !< 1 translation, 2 rotation, 3 Gauss, 4 Gauss + refine.
   real(R8P),    intent(out) :: delta !< Cells per step.
   integer(I4P)              :: k, is !< Work.

   k = ic - 1_I4P
   delta = merge(2._R8P, 0.5_R8P, mod(k, 2_I4P) == 0_I4P) ; k = k / 2_I4P
   path = mod(k, 4_I4P) + 1_I4P ; k = k / 4_I4P
   is = mod(k, 2_I4P) ; k = k / 2_I4P
   ib = mod(k, 2_I4P) + 1_I4P ; k = k / 2_I4P
   dim = k + 2_I4P
   if (dim == 2_I4P) then
      n = merge(256_I4P, 1024_I4P, is == 0_I4P)
   else
      n = merge(32_I4P, 64_I4P, is == 0_I4P)
   endif
   endsubroutine case_params

   function case_name(ic) result(s)
   !< Name of case ic.
   integer(I4P), intent(in)  :: ic    !< Case.
   character(:), allocatable :: s     !< Name.
   integer(I4P)              :: dim, ib, n, path !< Parameters.
   real(R8P)                 :: delta !< Cells per step.
   character(8)              :: d     !< Delta.

   call case_params(ic, dim, ib, n, path, delta)
   write(d, '(F3.1)') delta
   s = istr(dim) // 'd-' // trim(BNAME(ib, dim - 1)) // '-n' // istr(n) // '-' // &
       trim(PNAME(path)) // '-d' // trim(d)
   endfunction case_name

   function case_quick(ic) result(q)
   !< True for the cases of the small grids (256^2, 32^3).
   integer(I4P), intent(in) :: ic    !< Case.
   logical                  :: q     !< Quick.
   integer(I4P)             :: dim, ib, n, path !< Parameters.
   real(R8P)                :: delta !< Cells per step.

   call case_params(ic, dim, ib, n, path, delta)
   q = n == 256_I4P .or. n == 32_I4P
   endfunction case_quick

   subroutine setup_bodies
   !< The four bodies (the step-0 probe's).
   real(R8P), allocatable    :: u(:), v(:) !< Loop.
   type(body_t), allocatable :: bb(:)      !< Sign-bench bodies.
   type(vector_R8P)          :: b0, b1     !< Box.
   integer(I4P)              :: k          !< Counter.

   allocate(u(1024), v(1024))
   do k = 1_I4P, 1024_I4P
      u(k) = 0.5_R8P + 0.4_R8P * cos(2._R8P * PI * (k - 1) / 1024)
      v(k) = 0.5_R8P + 0.25_R8P * sin(2._R8P * PI * (k - 1) / 1024)
   enddo
   call curves(1)%adopt_loops(u=u, v=v, loop_sizes=[1024_I4P])
   deallocate(u, v) ; allocate(u(1022), v(1022))
   call naca0012(511_I4P, u, v)
   u = 0.15_R8P + 0.7_R8P * u ; v = 0.5_R8P + 0.7_R8P * v
   call curves(2)%adopt_loops(u=u, v=v, loop_sizes=[1022_I4P])
   call mc_body(BODY_SPHERE, 24_I4P, surfaces(1))
   call bench_bodies(bb) ; call to_surface(bb(7)%soup, surfaces(2))
   b0 = surfaces(2)%get_bmin() ; b1 = surfaces(2)%get_bmax()
   call surfaces(2)%translate(delta=-0.5_R8P * (b0 + b1)) ; call surfaces(2)%analyze
   half = [0.75_R8P, 1.1_R8P]
   endsubroutine setup_bodies

   subroutine case_line(ic, line, seconds)
   !< Build the map of case ic; its reference line and the time of the call.
   integer(I4P),              intent(in)  :: ic      !< Case.
   character(:), allocatable, intent(out) :: line    !< Name, counters, hashes.
   real(R8P),                 intent(out) :: seconds !< Time of the swept call with the map.
   integer(I4P)                           :: dim, ib, n, path, i !< Parameters, counter.
   real(R8P)                              :: delta, h, rmax, d3(3), d2(2), cen(2), dth !< Motion.
   real(R8P), allocatable                 :: x(:), y(:) !< Cell edges.
   type(rigid_pose_t)                     :: p0, p1, st !< 3-D poses.
   type(rigid_pose_2d_t)                  :: q0, q1     !< 2-D poses.
   type(swept_cells_3d_t)                 :: w3         !< 3-D swept geometry.
   type(swept_cells_2d_t)                 :: w2         !< 2-D swept geometry.
   type(piece_map_3d_t)                   :: m3         !< 3-D map.
   type(piece_map_2d_t)                   :: m2         !< 2-D map.
   type(vector_R8P)                       :: v0, v1     !< Box.
   integer(I8P)                           :: c0, c1, cr !< Clock.

   call case_params(ic, dim, ib, n, path, delta)
   if (dim == 3_I4P) then
      h = 2._R8P * half(ib) / n
      x = [(-half(ib) + h * i, i = 0, n)]
      v0 = surfaces(ib)%get_bmin() ; v1 = surfaces(ib)%get_bmax()
      rmax = norm2(max(abs([v0%x, v0%y, v0%z]), abs([v1%x, v1%y, v1%z])))
      p0 = pose_from_axis_angle([1._R8P, 2._R8P, 3._R8P], 0.1_R8P)
      if (path == 1_I4P) then
         d3 = [0.7_R8P, 0.4_R8P, -0.3_R8P] ; d3 = delta * h * d3 / norm2(d3)
         p1 = p0 ; p1%t = p0%t + d3
      else
         d3 = [0.3_R8P, 0.2_R8P, -0.1_R8P] ; d3 = 0.25_R8P * delta * h * d3 / norm2(d3)
         st = pose_from_axis_angle([0.2_R8P, 1._R8P, 0.3_R8P], 0.75_R8P * delta * h / rmax)
         p1%r = matmul(st%r, p0%r) ; p1%t = matmul(st%r, p0%t) + d3
      endif
      call system_clock(c0, cr)
      call swept_cells_3d(surfaces(ib), p0, p1, x, x, x, w3, exact=path <= 2_I4P, refine=path == 4_I4P, map=m3)
      call system_clock(c1)
      line = case_name(ic) // ' nx=' // istr(m3%nx) // ' ny=' // istr(m3%ny) // ' nz=' // &
             istr(m3%nz) // counters(m3%nt, m3%np0, m3%np1, m3%exact, m3%n_status, m3%untrusted_cells, &
             m3%multi_track_cells, m3%evaluated_cells, m3%evaluations, m3%skipped_nodes, m3%contact_cells, &
             m3%narrow_cells) // &
             field('tracks', m3%tracks) // field('first_track', m3%first_track) // field('track_cell', m3%track_cell) // &
             field('status', m3%status) // field('old_start', m3%old_start) // field('old_piece', m3%old_piece) // &
             field('new_start', m3%new_start) // field('new_piece', m3%new_piece) // &
             field('old_track', m3%old_track) // field('new_track', m3%new_track) // field('volume0', m3%volume0) // &
             field('volume1', m3%volume1) // field('faces', m3%faces) // field('outflow', m3%outflow) // &
             field('wall', m3%wall) // field('residual', m3%residual)
   else
      call grid_edges([0.0031_R8P, 0.0031_R8P], [0.9969_R8P, 0.9969_R8P], n, x, y)
      h = x(2) - x(1)
      rmax = merge(0.4_R8P, 0.35_R8P, ib == 1_I4P)
      cen = [0.5_R8P, 0.5_R8P]
      q0 = pose_2d(0.1_R8P) ; q0%t = cen - matmul(q0%r, cen)
      if (path == 1_I4P) then
         d2 = [0.7_R8P, 0.4_R8P] ; d2 = delta * h * d2 / norm2(d2)
         q1 = q0 ; q1%t = q0%t + d2
      else
         dth = 0.75_R8P * delta * h / rmax
         d2 = [0.3_R8P, 0.2_R8P] ; d2 = 0.25_R8P * delta * h * d2 / norm2(d2)
         q1 = pose_2d(0.1_R8P + dth) ; q1%t = cen - matmul(q1%r, cen) + d2
      endif
      call system_clock(c0, cr)
      call swept_cells_2d(curves(ib), q0, q1, x, y, w2, exact=path <= 2_I4P, refine=path == 4_I4P, map=m2)
      call system_clock(c1)
      line = case_name(ic) // ' nx=' // istr(m2%nx) // ' ny=' // istr(m2%ny) // &
             counters(m2%nt, m2%np0, m2%np1, m2%exact, m2%n_status, m2%untrusted_cells, &
             m2%multi_track_cells, m2%evaluated_cells, m2%evaluations, m2%skipped_nodes, m2%contact_cells, &
             m2%narrow_cells) // &
             field('tracks', m2%tracks) // field('first_track', m2%first_track) // field('track_cell', m2%track_cell) // &
             field('status', m2%status) // field('old_start', m2%old_start) // field('old_piece', m2%old_piece) // &
             field('new_start', m2%new_start) // field('new_piece', m2%new_piece) // &
             field('old_track', m2%old_track) // field('new_track', m2%new_track) // field('volume0', m2%volume0) // &
             field('volume1', m2%volume1) // field('faces', m2%faces) // field('outflow', m2%outflow) // &
             field('wall', m2%wall) // field('residual', m2%residual)
   endif
   seconds = real(c1 - c0, R8P) / real(cr, R8P)
   endsubroutine case_line

   function counters(nt, np0, np1, exact, n_status, untrusted, multi, evaluated, evaluations, skipped, contact, &
                     narrow) result(s)
   !< The scalar fields as key=value tokens.
   integer(I4P), intent(in)  :: nt, np0, np1        !< Tracks, pieces.
   logical,      intent(in)  :: exact               !< Converged.
   integer(I4P), intent(in)  :: n_status(-1:7)      !< Tracks per status.
   integer(I4P), intent(in)  :: untrusted, multi, evaluated, evaluations, skipped, contact, narrow !< Counters.
   character(:), allocatable :: s                   !< Tokens.
   integer(I4P)              :: k                   !< Counter.

   s = ' nt=' // istr(nt) // ' np0=' // istr(np0) // ' np1=' // istr(np1) // ' exact=' // &
       merge('T', 'F', exact) // ' n_status='
   do k = -1_I4P, 7_I4P
      s = s // istr(n_status(k)) // merge(',', ' ', k < 7_I4P)
   enddo
   s = s // 'untrusted_cells=' // istr(untrusted) // ' multi_track_cells=' // istr(multi) // &
       ' evaluated_cells=' // istr(evaluated) // ' evaluations=' // istr(evaluations) // &
       ' skipped_nodes=' // istr(skipped) // ' contact_cells=' // istr(contact) // ' narrow_cells=' // &
       istr(narrow)
   endfunction counters

   function istr(n) result(s)
   !< An integer without blanks.
   integer(I4P), intent(in)  :: n   !< Value.
   character(:), allocatable :: s   !< Digits.
   character(16)             :: buf !< Work.

   write(buf, '(I0)') n
   s = trim(buf)
   endfunction istr

   pure subroutine hash_word(h, w)
   !< Add a 64-bit word: its four 16-bit chunks, low first (`fossil_val_perf`).
   type(hash_t), intent(inout) :: h !< Hash.
   integer(I8P), intent(in)    :: w !< Word.
   integer(I8P)                :: c !< Chunk.
   integer(I4P)                :: k !< Counter.

   do k = 0, 3
      c = ibits(w, 16 * k, 16)
      h%a = mod(h%a * B1 + c + 1_I8P, P1) ; h%b = mod(h%b * B2 + c + 1_I8P, P2)
   enddo
   endsubroutine hash_word

   function token(name, alloc, shp, w) result(s)
   !< ' name=HEX': the hash of the allocation flag, the shape and the words.
   character(*), intent(in)  :: name   !< Field.
   logical,      intent(in)  :: alloc  !< Allocated.
   integer(I4P), intent(in)  :: shp(:) !< Shape (if allocated).
   integer(I8P), intent(in)  :: w(:)   !< Elements as words (if allocated).
   character(:), allocatable :: s      !< Token.
   type(hash_t)              :: h      !< Hash.
   character(16)             :: hx     !< Digits.
   integer(I4P)              :: k      !< Counter.

   call hash_word(h, merge(1_I8P, 0_I8P, alloc))
   if (alloc) then
      do k = 1_I4P, size(shp, kind=I4P)
         call hash_word(h, int(shp(k), I8P))
      enddo
      do k = 1_I4P, size(w, kind=I4P)
         call hash_word(h, w(k))
      enddo
   endif
   write(hx, '(Z8.8,Z8.8)') h%a, h%b
   s = ' ' // name // '=' // hx
   endfunction token

   function field_i1(name, a) result(s)
   !< Token of an integer rank-1 field.
   character(*),              intent(in) :: name !< Field.
   integer(I4P), allocatable, intent(in) :: a(:) !< Values.
   character(:), allocatable             :: s    !< Token.

   if (allocated(a)) then
      s = token(name, .true., shape(a, kind=I4P), int(a, I8P))
   else
      s = token(name, .false., [0_I4P], [0_I8P])
   endif
   endfunction field_i1

   function field_i2(name, a) result(s)
   !< Token of an integer rank-2 field.
   character(*),              intent(in) :: name   !< Field.
   integer(I4P), allocatable, intent(in) :: a(:,:) !< Values.
   character(:), allocatable             :: s      !< Token.

   if (allocated(a)) then
      s = token(name, .true., shape(a, kind=I4P), int(reshape(a, [size(a)]), I8P))
   else
      s = token(name, .false., [0_I4P], [0_I8P])
   endif
   endfunction field_i2

   function field_i3(name, a) result(s)
   !< Token of an integer rank-3 field.
   character(*),              intent(in) :: name     !< Field.
   integer(I4P), allocatable, intent(in) :: a(:,:,:) !< Values.
   character(:), allocatable             :: s        !< Token.

   if (allocated(a)) then
      s = token(name, .true., shape(a, kind=I4P), int(reshape(a, [size(a)]), I8P))
   else
      s = token(name, .false., [0_I4P], [0_I8P])
   endif
   endfunction field_i3

   function field_r1(name, a) result(s)
   !< Token of a real rank-1 field (bits).
   character(*),           intent(in) :: name !< Field.
   real(R8P), allocatable, intent(in) :: a(:) !< Values.
   character(:), allocatable          :: s    !< Token.

   if (allocated(a)) then
      s = token(name, .true., shape(a, kind=I4P), transfer(a, 0_I8P, size(a)))
   else
      s = token(name, .false., [0_I4P], [0_I8P])
   endif
   endfunction field_r1

   function field_r2(name, a) result(s)
   !< Token of a real rank-2 field (bits).
   character(*),           intent(in) :: name   !< Field.
   real(R8P), allocatable, intent(in) :: a(:,:) !< Values.
   character(:), allocatable          :: s      !< Token.

   if (allocated(a)) then
      s = token(name, .true., shape(a, kind=I4P), transfer(a, 0_I8P, size(a)))
   else
      s = token(name, .false., [0_I4P], [0_I8P])
   endif
   endfunction field_r2

   subroutine split_tokens(line, tok, nt)
   !< The blank-separated tokens of a line.
   character(*),               intent(in)  :: line   !< Line.
   character(128), allocatable, intent(out) :: tok(:) !< Tokens.
   integer(I4P),               intent(out) :: nt     !< Their number.
   integer(I4P)                            :: i, j, l !< Work.

   allocate(tok(64)) ; nt = 0_I4P
   l = len_trim(line, kind=I4P) ; i = 1_I4P
   do while (i <= l)
      if (line(i:i) == ' ') then
         i = i + 1_I4P ; cycle
      endif
      j = i
      do while (j < l)
         if (line(j + 1:j + 1) == ' ') exit
         j = j + 1_I4P
      enddo
      nt = nt + 1_I4P
      if (nt <= size(tok, kind=I4P)) tok(nt) = line(i:j)
      i = j + 1_I4P
   enddo
   nt = min(nt, size(tok, kind=I4P))
   endsubroutine split_tokens

   subroutine grid_edges(lo, hi, n, x, y)
   !< The device-chain grid: [lo, hi] widened by 2 % per side, n cells per direction.
   real(R8P),              intent(in)  :: lo(2), hi(2) !< Body box.
   integer(I4P),           intent(in)  :: n            !< Cells.
   real(R8P), allocatable, intent(out) :: x(:), y(:)   !< Edges.
   real(R8P)                           :: a(2), b(2)   !< Grid box.
   integer(I4P)                        :: k            !< Counter.

   a = lo - 0.02_R8P * (hi - lo) ; b = hi + 0.02_R8P * (hi - lo)
   allocate(x(n + 1), y(n + 1))
   do k = 1_I4P, n + 1_I4P
      x(k) = a(1) + (b(1) - a(1)) * (k - 1) / real(n, R8P) ; y(k) = a(2) + (b(2) - a(2)) * (k - 1) / real(n, R8P)
   enddo
   endsubroutine grid_edges

   subroutine naca0012(m, x, y)
   !< NACA 0012 with 2 m points, cosine spacing, sharp trailing edge.
   integer(I4P), intent(in)  :: m      !< Points per side.
   real(R8P),    intent(out) :: x(:), y(:) !< Points.
   real(R8P)                 :: xi     !< Chord station.
   integer(I4P)              :: kk, j  !< Counters.

   j = 0_I4P
   do kk = m, 0_I4P, -1_I4P
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1_I4P ; x(j) = xi ; y(j) = ht(xi)
   enddo
   do kk = 1_I4P, m - 1_I4P
      xi = 0.5_R8P * (1._R8P - cos(PI * kk / m))
      j = j + 1_I4P ; x(j) = xi ; y(j) = -ht(xi)
   enddo
   endsubroutine naca0012

   pure function ht(x) result(y)
   !< Half thickness of the NACA 0012 at chord station x.
   real(R8P), intent(in) :: x !< Station.
   real(R8P)             :: y !< Half thickness.

   y = 0.6_R8P * (0.2969_R8P * sqrt(x) - 0.1260_R8P * x - 0.3516_R8P * x**2 + 0.2843_R8P * x**3 - 0.1036_R8P * x**4)
   endfunction ht
endmodule fossil_val_piece_map_ref_kit

program fossil_val_piece_map_ref
!< Write or check the piece-map reference (see the file notes).
use fossil_val_piece_map_ref_kit, only : NCASES, MAXLEN, case_name, case_quick, case_line, setup_bodies, split_tokens
use penf,                         only : I4P, R8P

implicit none

character(*), parameter     :: DEFAULT_FILE = 'src/tests/validation/perf/map-ref-f855194.txt' !< Committed reference.
character(MAXLEN)           :: arg, mode, file, note        !< Arguments.
character(MAXLEN), allocatable :: ref(:)                    !< Reference lines.
character(:), allocatable   :: line                         !< Computed line.
character(128), allocatable :: ta(:), tb(:)                 !< Tokens (reference, computed).
logical                     :: sel(NCASES), quick, exists, found, ok !< Selection, flags.
integer(I4P)                :: ia, ic, k, nref, u, ios, na, nb, ndiff, nfail, nrun !< Work.
real(R8P)                   :: sec, total                   !< Times.

mode = 'check' ; file = DEFAULT_FILE ; note = '' ; quick = .false. ; sel = .false.
if (command_argument_count() == 0) then
   quick = .true.
else
   call get_command_argument(1, mode) ; call get_command_argument(2, file)
   do ia = 3, command_argument_count()
      call get_command_argument(ia, arg)
      if (trim(arg) == 'quick') then
         quick = .true.
      elseif (arg(1:5) == 'case=') then
         read(arg(6:), *) ic ; if (ic >= 1_I4P .and. ic <= NCASES) sel(ic) = .true.
      elseif (arg(1:5) == 'note=') then
         note = arg(6:)
      else
         print '(A,A)', 'unknown argument ', trim(arg) ; stop 2
      endif
   enddo
endif
if (.not. any(sel)) sel = .true.
if (quick) then
   do ic = 1_I4P, NCASES
      sel(ic) = sel(ic) .and. case_quick(ic)
   enddo
endif
if (trim(mode) /= 'ref' .and. trim(mode) /= 'check') then
   print '(A)', 'usage: fossil_val_piece_map_ref ref|check FILE [quick] [case=N ...] [note=TEXT]' ; stop 2
endif
inquire(file=trim(file), exist=exists)
if (trim(mode) == 'ref') then
   if (exists) then
      print '(A,A,A)', 'refused: ', trim(file), ' exists (references are never overwritten)' ; stop 1
   endif
   open(newunit=u, file=trim(file), status='new', action='write')
   write(u, '(A)') '# piece map reference (fossil_val_piece_map_ref, issue #105 step 1a): case, counters, field hashes'
   if (len_trim(note) > 0) write(u, '(A)') '# ' // trim(note)
else
   if (.not. exists) then
      print '(A,A)', 'no reference ', trim(file) ; stop 1
   endif
   allocate(ref(NCASES + 16)) ; nref = 0_I4P
   open(newunit=u, file=trim(file), status='old', action='read')
   do
      read(u, '(A)', iostat=ios) arg
      if (ios /= 0) exit
      if (arg(1:1) == '#' .or. len_trim(arg) == 0) cycle
      nref = nref + 1_I4P
      if (nref > size(ref, kind=I4P)) then
         print '(A)', 'reference with too many lines' ; stop 1
      endif
      ref(nref) = arg
   enddo
   close(u)
endif
call setup_bodies
ok = .true. ; nfail = 0_I4P ; nrun = 0_I4P ; total = 0._R8P
do ic = 1_I4P, NCASES
   if (.not. sel(ic)) cycle
   call case_line(ic, line, sec)
   nrun = nrun + 1_I4P ; total = total + sec
   if (len(line) > MAXLEN) then
      print '(A)', 'line too long' ; stop 1
   endif
   if (trim(mode) == 'ref') then
      write(u, '(A)') line ; flush(u)
      print '(A,A,A,F9.2,A)', '  ', case_name(ic), ' written (', sec, ' s)'
      cycle
   endif
   found = .false.
   do k = 1_I4P, nref
      if (index(ref(k), case_name(ic) // ' ') == 1) then
         found = .true. ; exit
      endif
   enddo
   if (.not. found) then
      print '(A,A,A,F9.2,A)', '  ', case_name(ic), ': MISSING in the reference (', sec, ' s)'
      ok = .false. ; nfail = nfail + 1_I4P ; cycle
   endif
   call split_tokens(ref(k), ta, na) ; call split_tokens(line, tb, nb)
   ndiff = 0_I4P
   do ia = 1_I4P, max(na, nb)
      if (ia > na .or. ia > nb) then
         ndiff = ndiff + 1_I4P ; cycle
      endif
      if (ta(ia) /= tb(ia)) then
         ndiff = ndiff + 1_I4P
         print '(A,A,A,A)', '    reference ', trim(ta(ia)), '  computed ', trim(tb(ia))
      endif
   enddo
   if (ndiff == 0_I4P) then
      print '(A,A,A,F9.2,A)', '  ', case_name(ic), ': same (', sec, ' s)'
   else
      print '(A,A,A,I0,A,F9.2,A)', '  ', case_name(ic), ': DIFFERS in ', ndiff, ' fields (', sec, ' s)'
      ok = .false. ; nfail = nfail + 1_I4P
   endif
enddo
if (trim(mode) == 'ref') close(u)
print '(A,I0,A,I0,A,F10.1,A)', 'cases ', nrun, ', differing or missing ', nfail, ', swept calls ', total, ' s'
print '(A,L1)', 'Are all tests passed? ', ok
if (.not. ok) stop 1
endprogram fossil_val_piece_map_ref
