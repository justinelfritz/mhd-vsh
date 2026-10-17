MODULE GRID_RADIAL
!> Radial grid: finite-difference nodes, either a spherical shell
!> (r_min..r_max, genuine boundary conditions at both ends) or a full
!> sphere (r=0..r_max, half-step-staggered so no node sits exactly at the
!> origin, where `l(l+1)/r**2` is singular for l>=1).
!>
!> @note RADIAL_GEOMETRIC_OUTER (2026-10-16, per user request): clusters
!>   nodes toward R_MAX (shell grids only) via a geometric cell-width
!>   sequence, so a Hall-active near-surface layer gets genuinely more
!>   resolution instead of the same uniform dr as the rest of the domain.
!>   Motivated directly by a measured gap: on the production driver's own
!>   N_R=40 uniform grid, Rmag at the outermost node reached ~67 for the
!>   l=1 1e12G force-free dipole IC, while the SAME physical quantity
!>   evaluated on a much finer (387-point, log(P)-clustered) TOV-profile
!>   grid reached ~8e5 at the true surface -- a ~4-order-of-magnitude
!>   gap purely from the uniform grid's last node sitting measurably
!>   short of where the physics actually peaks (see
!>   results/ic_rmag_compare/). RADIAL_OPERATORS::BUILD_RADIAL_OPERATORS
!>   already builds its FD stencils via FD_WEIGHTS_GENERIC, which works
!>   for arbitrary (not just uniform) node spacing, so no change is
!>   needed there -- only this module needed a way to actually produce a
!>   non-uniform R(:) array.
USE KINDS, ONLY: dp, i4
IMPLICIT NONE
PRIVATE
PUBLIC :: RADIAL_GRID_T, BUILD_RADIAL_GRID, BUILD_RADIAL_GRID_FROM_ARRAY, &
          BUILD_RADIAL_GRID_FROM_FILE, RADIAL_UNIFORM, RADIAL_GEOMETRIC_OUTER, RADIAL_EXPLICIT

INTEGER(KIND=i4), PARAMETER :: RADIAL_UNIFORM         = 0
!> Geometric clustering toward R_MAX. STRETCH_PARAM is the cell-width
!> ratio q in (0,1): consecutive cell widths shrink by a factor q moving
!> outward (cell_i = cell_0 * q**i, i=0..N-2), so SMALLER q means
!> stronger clustering at the surface. q close to 1 reduces to
!> (nearly) uniform spacing.
INTEGER(KIND=i4), PARAMETER :: RADIAL_GEOMETRIC_OUTER = 1
!> Set by BUILD_RADIAL_GRID_FROM_ARRAY only -- record-keeping value for
!> RGRID%STRETCH_TYPE, not dispatched on inside BUILD_RADIAL_GRID itself.
INTEGER(KIND=i4), PARAMETER :: RADIAL_EXPLICIT        = 2

TYPE :: RADIAL_GRID_T
  INTEGER(KIND=i4) :: N            = 0
  LOGICAL          :: FULL_SPHERE  = .FALSE.
  REAL(KIND=dp)    :: R_MIN        = 0.0_dp
  REAL(KIND=dp)    :: R_MAX        = 0.0_dp
  INTEGER(KIND=i4) :: STRETCH_TYPE = RADIAL_UNIFORM
  REAL(KIND=dp)    :: STRETCH_PARAM = 0.0_dp
  INTEGER(KIND=i4) :: FD_ORDER     = 4
  REAL(KIND=dp), ALLOCATABLE :: R(:)
END TYPE RADIAL_GRID_T

CONTAINS

!> @param RGRID Output grid.
!> @param N Number of radial nodes.
!> @param R_MIN Inner boundary (ignored, forced to 0, when FULL_SPHERE).
!> @param R_MAX Outer boundary.
!> @param FULL_SPHERE .TRUE. for a full sphere (half-step-staggered grid,
!>   origin regularity handled by RADIAL_OPERATORS), .FALSE. for a shell
!>   (uniform grid spanning [R_MIN,R_MAX] with genuine BCs at both ends).
!> @param FD_ORDER Optional nominal interior finite-difference order,
!>   default 4.
!> @param STRETCH_TYPE Optional, default RADIAL_UNIFORM. RADIAL_GEOMETRIC_OUTER
!>   requires .NOT. FULL_SPHERE (hard STOP otherwise -- the full-sphere
!>   half-step-staggered construction is a different, unrelated scheme,
!>   not something this clustering generalizes) and STRETCH_PARAM present.
!> @param STRETCH_PARAM Optional; required (hard STOP if absent) when
!>   STRETCH_TYPE=RADIAL_GEOMETRIC_OUTER. See RADIAL_GEOMETRIC_OUTER's
!>   own declaration comment for its meaning.
SUBROUTINE BUILD_RADIAL_GRID(RGRID, N, R_MIN, R_MAX, FULL_SPHERE, FD_ORDER, &
    STRETCH_TYPE, STRETCH_PARAM)
  TYPE(RADIAL_GRID_T), INTENT(OUT) :: RGRID
  INTEGER(KIND=i4),    INTENT(IN)  :: N
  REAL(KIND=dp),       INTENT(IN)  :: R_MIN, R_MAX
  LOGICAL,              INTENT(IN)  :: FULL_SPHERE
  INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: FD_ORDER
  INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: STRETCH_TYPE
  REAL(KIND=dp),    INTENT(IN), OPTIONAL :: STRETCH_PARAM
  INTEGER(KIND=i4) :: I, REQ_STRETCH
  REAL(KIND=dp) :: DR, Q, CELL0, WIDTH, S

  REQ_STRETCH = RADIAL_UNIFORM
  IF (PRESENT(STRETCH_TYPE)) REQ_STRETCH = STRETCH_TYPE

  RGRID%N            = N
  RGRID%FULL_SPHERE   = FULL_SPHERE
  RGRID%STRETCH_TYPE  = REQ_STRETCH
  RGRID%STRETCH_PARAM = 0.0_dp
  IF (PRESENT(STRETCH_PARAM)) RGRID%STRETCH_PARAM = STRETCH_PARAM
  RGRID%FD_ORDER      = 4
  IF (PRESENT(FD_ORDER)) RGRID%FD_ORDER = FD_ORDER

  ALLOCATE(RGRID%R(N))

  IF (FULL_SPHERE) THEN
    IF (REQ_STRETCH /= RADIAL_UNIFORM) THEN
      WRITE(*,'(A)') 'BUILD_RADIAL_GRID: STRETCH_TYPE/=RADIAL_UNIFORM is not supported with FULL_SPHERE'
      STOP 1
    END IF
    RGRID%R_MIN = 0.0_dp
    RGRID%R_MAX = R_MAX
    DR = R_MAX / (REAL(N, KIND=dp) - 0.5_dp)
    DO I = 1, N
      RGRID%R(I) = (REAL(I, KIND=dp) - 0.5_dp) * DR
    END DO
  ELSE
    RGRID%R_MIN = R_MIN
    RGRID%R_MAX = R_MAX
    SELECT CASE (REQ_STRETCH)
    CASE (RADIAL_UNIFORM)
      DR = (R_MAX - R_MIN) / REAL(N-1, KIND=dp)
      DO I = 1, N
        RGRID%R(I) = R_MIN + REAL(I-1, KIND=dp)*DR
      END DO
    CASE (RADIAL_GEOMETRIC_OUTER)
      IF (.NOT. PRESENT(STRETCH_PARAM)) THEN
        WRITE(*,'(A)') 'BUILD_RADIAL_GRID: STRETCH_PARAM required for RADIAL_GEOMETRIC_OUTER'
        STOP 1
      END IF
      Q = STRETCH_PARAM
      IF (.NOT. (Q > 0.0_dp .AND. Q <= 1.0_dp)) THEN
        WRITE(*,'(A)') 'BUILD_RADIAL_GRID: RADIAL_GEOMETRIC_OUTER requires 0 < STRETCH_PARAM <= 1'
        STOP 1
      END IF
      WIDTH = R_MAX - R_MIN
      ! Sum of a geometric series of N-1 cell widths, ratio Q, equal to
      ! WIDTH: CELL0*(1+Q+Q**2+...+Q**(N-2)) = WIDTH. Q=1 (uniform)
      ! handled separately since the closed form divides by (1-Q).
      IF (ABS(Q - 1.0_dp) < 1.0E-14_dp) THEN
        CELL0 = WIDTH / REAL(N-1, KIND=dp)
      ELSE
        CELL0 = WIDTH * (1.0_dp - Q) / (1.0_dp - Q**(N-1))
      END IF
      RGRID%R(1) = R_MIN
      S = R_MIN
      DO I = 2, N
        S = S + CELL0*Q**(I-2)
        RGRID%R(I) = S
      END DO
      RGRID%R(N) = R_MAX   ! exact, not accumulated-roundoff
    CASE DEFAULT
      WRITE(*,'(A,I0)') 'BUILD_RADIAL_GRID: unknown STRETCH_TYPE=', REQ_STRETCH
      STOP 1
    END SELECT
  END IF
END SUBROUTINE BUILD_RADIAL_GRID

!> Builds a shell RADIAL_GRID_T directly from a caller-supplied, already
!> strictly-increasing radii array -- e.g. a TOV profile's own crust-range
!> subset, reusing whatever physically-motivated point distribution was
!> already computed elsewhere (dense wherever density/pressure change
!> fastest -- which, for this project's crust, is the same near-surface
!> region the Hall dynamics are most sensitive to) instead of fitting a
!> separate synthetic formula (cf. RADIAL_GEOMETRIC_OUTER, which guesses
!> a clustering ratio by hand). R_MIN/R_MAX are taken directly from
!> R_ARRAY's own first/last entries. FULL_SPHERE is not supported this
!> way -- it needs the half-step-staggered origin construction
!> (APPLY_ORIGIN_REGULARITY's own convention), not an arbitrary external
!> array that may not even include r=0.
!> @param RGRID Output grid.
!> @param R_ARRAY Strictly increasing radii, size>=2 (hard STOP
!>   otherwise -- a non-monotonic or degenerate array would silently
!>   produce nonsense finite-difference stencils downstream).
!> @param FD_ORDER Optional nominal interior finite-difference order,
!>   default 4.
SUBROUTINE BUILD_RADIAL_GRID_FROM_ARRAY(RGRID, R_ARRAY, FD_ORDER)
  TYPE(RADIAL_GRID_T), INTENT(OUT) :: RGRID
  REAL(KIND=dp),       INTENT(IN)  :: R_ARRAY(:)
  INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: FD_ORDER
  INTEGER(KIND=i4) :: N, I

  N = SIZE(R_ARRAY)
  IF (N < 2) THEN
    WRITE(*,'(A)') 'BUILD_RADIAL_GRID_FROM_ARRAY: R_ARRAY needs at least 2 points'
    STOP 1
  END IF
  DO I = 1, N-1
    IF (R_ARRAY(I+1) <= R_ARRAY(I)) THEN
      WRITE(*,'(A,I0,A,ES16.8,A,ES16.8)') &
        'BUILD_RADIAL_GRID_FROM_ARRAY: R_ARRAY not strictly increasing at i=', I, &
        ': ', R_ARRAY(I), ' -> ', R_ARRAY(I+1)
      STOP 1
    END IF
  END DO

  RGRID%N             = N
  RGRID%FULL_SPHERE    = .FALSE.
  RGRID%R_MIN          = R_ARRAY(1)
  RGRID%R_MAX          = R_ARRAY(N)
  RGRID%STRETCH_TYPE   = RADIAL_EXPLICIT
  RGRID%STRETCH_PARAM  = 0.0_dp
  RGRID%FD_ORDER       = 4
  IF (PRESENT(FD_ORDER)) RGRID%FD_ORDER = FD_ORDER
  ALLOCATE(RGRID%R(N))
  RGRID%R = R_ARRAY
END SUBROUTINE BUILD_RADIAL_GRID_FROM_ARRAY

!> Reads a radial grid written by app/mhdvsh_tov.f90's equidistribution
!> step (format: a single '#'-prefixed header line, then rows "i
!> r_km dr_km"; dr_km is not re-used here -- BUILD_RADIAL_GRID_FROM_ARRAY
!> recomputes whatever it needs from R(:) itself -- only r_km is read).
!> Delegates straight to BUILD_RADIAL_GRID_FROM_ARRAY once the column
!> is in memory, so it inherits the same strictly-increasing check.
!> @param RGRID Output grid.
!> @param PATH Path to the grid file.
!> @param FD_ORDER Optional nominal interior finite-difference order,
!>   default 4.
SUBROUTINE BUILD_RADIAL_GRID_FROM_FILE(RGRID, PATH, FD_ORDER)
  TYPE(RADIAL_GRID_T), INTENT(OUT) :: RGRID
  CHARACTER(LEN=*),    INTENT(IN)  :: PATH
  INTEGER(KIND=i4), INTENT(IN), OPTIONAL :: FD_ORDER
  INTEGER(KIND=i4), PARAMETER :: FILE_UNIT = 71
  INTEGER(KIND=i4) :: IOS, N, IDX, I
  CHARACTER(LEN=1024) :: LINE
  REAL(KIND=dp) :: R_VAL, DR_VAL
  REAL(KIND=dp), ALLOCATABLE :: R_ARRAY(:)

  OPEN(UNIT=FILE_UNIT, FILE=TRIM(PATH), STATUS='OLD', ACTION='READ', IOSTAT=IOS)
  IF (IOS /= 0) THEN
    WRITE(*,'(A,A)') 'BUILD_RADIAL_GRID_FROM_FILE: cannot open ', TRIM(PATH)
    STOP 1
  END IF

  ! First pass: count data rows (non-comment lines).
  N = 0
  DO
    READ(FILE_UNIT, '(A)', IOSTAT=IOS) LINE
    IF (IOS /= 0) EXIT
    IF (LEN_TRIM(LINE) == 0) CYCLE
    IF (LINE(1:1) == '#') CYCLE
    N = N + 1
  END DO

  IF (N < 2) THEN
    WRITE(*,'(A,A)') 'BUILD_RADIAL_GRID_FROM_FILE: fewer than 2 data rows in ', TRIM(PATH)
    STOP 1
  END IF

  ALLOCATE(R_ARRAY(N))
  REWIND(FILE_UNIT)
  I = 0
  DO
    READ(FILE_UNIT, '(A)', IOSTAT=IOS) LINE
    IF (IOS /= 0) EXIT
    IF (LEN_TRIM(LINE) == 0) CYCLE
    IF (LINE(1:1) == '#') CYCLE
    I = I + 1
    READ(LINE, *) IDX, R_VAL, DR_VAL
    R_ARRAY(I) = R_VAL
  END DO
  CLOSE(FILE_UNIT)

  IF (PRESENT(FD_ORDER)) THEN
    CALL BUILD_RADIAL_GRID_FROM_ARRAY(RGRID, R_ARRAY, FD_ORDER=FD_ORDER)
  ELSE
    CALL BUILD_RADIAL_GRID_FROM_ARRAY(RGRID, R_ARRAY)
  END IF
END SUBROUTINE BUILD_RADIAL_GRID_FROM_FILE

END MODULE GRID_RADIAL
