!> Regression test for GRID_RADIAL, covering the new (2026-10-16)
!> RADIAL_GEOMETRIC_OUTER stretching option and BUILD_RADIAL_GRID_FROM_ARRAY,
!> alongside the existing RADIAL_UNIFORM default. Five checks:
!>  (1) RADIAL_UNIFORM: unchanged from before -- R(1)=R_min, R(N)=R_max,
!>      uniform spacing, confirming the new optional STRETCH_TYPE/
!>      STRETCH_PARAM arguments don't disturb existing callers;
!>  (2) RADIAL_GEOMETRIC_OUTER: R(1)=R_min, R(N)=R_max exactly,
!>      monotonically increasing, and consecutive cell widths follow the
!>      requested ratio Q to tight tolerance (direct check of the
!>      geometric-series construction, not just "looks clustered");
!>  (3) smaller Q clusters more aggressively toward R_max than larger Q
!>      -- the outermost cell is strictly smaller for the smaller-Q grid,
!>      confirming the stretching direction and relative strength are
!>      both correct, not just that SOME non-uniform grid was produced;
!>  (4) BUILD_RADIAL_GRID_FROM_ARRAY: an arbitrary (not formula-generated)
!>      strictly-increasing array round-trips exactly into RGRID%R, with
!>      R_MIN/R_MAX/N taken from the array itself;
!>  (5) BUILD_RADIAL_GRID_FROM_ARRAY feeds RADIAL_OPERATORS::
!>      BUILD_RADIAL_OPERATORS without error and the resulting D1
!>      differentiates a known polynomial correctly on a non-uniform
!>      array -- confirms FD_WEIGHTS_GENERIC's arbitrary-spacing claim is
!>      actually exercised correctly for an externally-supplied (not
!>      formula-clustered) grid, not just asserted.
!>  (6) BUILD_RADIAL_GRID_FROM_FILE: writes the same hand-picked R_CUSTOM
!>      array from (4) to a file in the same "i r_km dr_km" format
!>      app/mhdvsh_tov.f90's equidistribution step produces (including a
!>      leading '#' header line, exercising that the reader skips it),
!>      reads it back, and checks an exact round-trip of R(:) -- the
!>      dr_km column itself is not re-checked here since
!>      BUILD_RADIAL_GRID_FROM_FILE recomputes everything downstream from
!>      R(:) alone.
PROGRAM TEST_GRID_RADIAL
USE KINDS,            ONLY: dp, i4
USE GRID_RADIAL,      ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID, BUILD_RADIAL_GRID_FROM_ARRAY, &
                             BUILD_RADIAL_GRID_FROM_FILE, RADIAL_UNIFORM, RADIAL_GEOMETRIC_OUTER
USE RADIAL_OPERATORS, ONLY: RADIAL_OPERATOR_T, BUILD_RADIAL_OPERATORS
IMPLICIT NONE

REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-10_dp
INTEGER(KIND=i4), PARAMETER :: N_R = 20_i4
REAL(KIND=dp),    PARAMETER :: R_MIN = 10.0_dp, R_MAX = 12.0_dp
INTEGER(KIND=i4), PARAMETER :: FILE_UNIT = 81_i4
CHARACTER(LEN=*), PARAMETER :: TMP_GRID_PATH = 'test_grid_radial_tmp.dat'

TYPE(RADIAL_GRID_T) :: RGRID_UNIFORM, RGRID_Q8, RGRID_Q5, RGRID_ARR, RGRID_FILE
TYPE(RADIAL_OPERATOR_T) :: OPS_ARR
INTEGER(KIND=i4) :: N_FAIL, I
REAL(KIND=dp)    :: DR_UNIFORM, CELL_WIDTH, PREV_WIDTH, RATIO
REAL(KIND=dp), ALLOCATABLE :: R_CUSTOM(:), F(:), DF_NUMERIC(:)

N_FAIL = 0

! --- (1) RADIAL_UNIFORM: existing behavior unchanged ---
CALL BUILD_RADIAL_GRID(RGRID_UNIFORM, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL CHECK_CLOSE("uniform_r_min", RGRID_UNIFORM%R(1), R_MIN, TOL, N_FAIL)
CALL CHECK_CLOSE("uniform_r_max", RGRID_UNIFORM%R(N_R), R_MAX, TOL, N_FAIL)
DR_UNIFORM = (R_MAX-R_MIN)/REAL(N_R-1, KIND=dp)
DO I = 1, N_R
  CALL CHECK_CLOSE("uniform_spacing", RGRID_UNIFORM%R(I), R_MIN+REAL(I-1,KIND=dp)*DR_UNIFORM, TOL, N_FAIL)
END DO

! --- (2) RADIAL_GEOMETRIC_OUTER: boundaries exact, ratio matches Q ---
CALL BUILD_RADIAL_GRID(RGRID_Q8, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE., &
  STRETCH_TYPE=RADIAL_GEOMETRIC_OUTER, STRETCH_PARAM=0.8_dp)
CALL CHECK_CLOSE("geo_r_min", RGRID_Q8%R(1), R_MIN, TOL, N_FAIL)
CALL CHECK_CLOSE("geo_r_max", RGRID_Q8%R(N_R), R_MAX, TOL, N_FAIL)

PREV_WIDTH = RGRID_Q8%R(2) - RGRID_Q8%R(1)
DO I = 2, N_R-1
  CALL CHECK_TRUE("geo_monotonic", RGRID_Q8%R(I+1) > RGRID_Q8%R(I), N_FAIL)
  CELL_WIDTH = RGRID_Q8%R(I+1) - RGRID_Q8%R(I)
  RATIO = CELL_WIDTH / PREV_WIDTH
  CALL CHECK_CLOSE("geo_ratio_matches_q", RATIO, 0.8_dp, TOL, N_FAIL)
  PREV_WIDTH = CELL_WIDTH
END DO

! --- (3) smaller Q clusters more aggressively than larger Q ---
CALL BUILD_RADIAL_GRID(RGRID_Q5, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE., &
  STRETCH_TYPE=RADIAL_GEOMETRIC_OUTER, STRETCH_PARAM=0.5_dp)
CALL CHECK_TRUE("smaller_q_smaller_outer_cell", &
  (RGRID_Q5%R(N_R)-RGRID_Q5%R(N_R-1)) < (RGRID_Q8%R(N_R)-RGRID_Q8%R(N_R-1)), N_FAIL)
! Both must still span the exact same domain despite the different ratio.
CALL CHECK_CLOSE("q5_r_min", RGRID_Q5%R(1), R_MIN, TOL, N_FAIL)
CALL CHECK_CLOSE("q5_r_max", RGRID_Q5%R(N_R), R_MAX, TOL, N_FAIL)

! --- (4) BUILD_RADIAL_GRID_FROM_ARRAY: exact round-trip ---
ALLOCATE(R_CUSTOM(7))
R_CUSTOM = (/ 10.0_dp, 10.1_dp, 10.4_dp, 10.9_dp, 11.5_dp, 11.8_dp, 12.0_dp /)  ! irregular, hand-picked
CALL BUILD_RADIAL_GRID_FROM_ARRAY(RGRID_ARR, R_CUSTOM)
CALL CHECK_TRUE("array_n_matches", RGRID_ARR%N == SIZE(R_CUSTOM), N_FAIL)
CALL CHECK_CLOSE("array_r_min", RGRID_ARR%R_MIN, R_CUSTOM(1), TOL, N_FAIL)
CALL CHECK_CLOSE("array_r_max", RGRID_ARR%R_MAX, R_CUSTOM(SIZE(R_CUSTOM)), TOL, N_FAIL)
DO I = 1, SIZE(R_CUSTOM)
  CALL CHECK_CLOSE("array_roundtrip", RGRID_ARR%R(I), R_CUSTOM(I), TOL, N_FAIL)
END DO

! --- (5) The resulting grid feeds BUILD_RADIAL_OPERATORS and D1
! differentiates f(r)=r**2 correctly (df/dr=2r) at interior points, on
! this genuinely irregular (not formula-clustered) spacing.
CALL BUILD_RADIAL_OPERATORS(OPS_ARR, RGRID_ARR)
ALLOCATE(F(RGRID_ARR%N), DF_NUMERIC(RGRID_ARR%N))
F = RGRID_ARR%R**2
DF_NUMERIC = MATMUL(OPS_ARR%D1, F)
DO I = 2, RGRID_ARR%N-1
  CALL CHECK_CLOSE("array_d1_matches_analytic", DF_NUMERIC(I), 2.0_dp*RGRID_ARR%R(I), TOL, N_FAIL)
END DO

! --- (6) BUILD_RADIAL_GRID_FROM_FILE: round-trip through a file in
! app/mhdvsh_tov.f90's own output format.
OPEN(UNIT=FILE_UNIT, FILE=TMP_GRID_PATH, STATUS='REPLACE', ACTION='WRITE')
WRITE(FILE_UNIT,'(A)') '# i  r_km  dr_km'
DO I = 1, SIZE(R_CUSTOM)
  WRITE(FILE_UNIT,'(I6,2ES16.8)') I, R_CUSTOM(I), 0.0_dp   ! dr_km unused by the reader
END DO
CLOSE(FILE_UNIT)

CALL BUILD_RADIAL_GRID_FROM_FILE(RGRID_FILE, TMP_GRID_PATH)
CALL CHECK_TRUE("file_n_matches", RGRID_FILE%N == SIZE(R_CUSTOM), N_FAIL)
DO I = 1, SIZE(R_CUSTOM)
  CALL CHECK_CLOSE("file_roundtrip", RGRID_FILE%R(I), R_CUSTOM(I), TOL, N_FAIL)
END DO
OPEN(UNIT=FILE_UNIT, FILE=TMP_GRID_PATH, STATUS='OLD')
CLOSE(FILE_UNIT, STATUS='DELETE')

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_grid_radial"
END IF

CONTAINS

SUBROUTINE CHECK_TRUE(NAME, COND, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  LOGICAL,          INTENT(IN)    :: COND
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (.NOT. COND) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A)') "FAIL  ", NAME
  ELSE
    WRITE(*,'(A,A)') "PASS  ", NAME
  END IF
END SUBROUTINE CHECK_TRUE

SUBROUTINE CHECK_CLOSE(NAME, VAL, EXPECTED, TOL, N_FAIL)
  CHARACTER(LEN=*), INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: VAL, EXPECTED, TOL
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ERR
  ERR = ABS(VAL-EXPECTED)/MAX(ABS(EXPECTED), 1.0E-300_dp)
  IF (ERR > TOL) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  relerr=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  relerr=", ERR
  END IF
END SUBROUTINE CHECK_CLOSE

END PROGRAM TEST_GRID_RADIAL
