!> Regression test for INITIAL_CONDITIONS. Three checks:
!>  (1) READ_IC_NAMELIST round-trip -- write a known namelist file
!>      (list-valued &IC_SELECT TAGS, confirmed with user 2026-08-29),
!>      read it back, confirm the parsed IC_TAGS(:)/IC_MODE_T array
!>      match exactly;
!>  (2) BUILD_INITIAL_CONDITION(['bessel_riccati'], ...) reproduces a
!>      hand-computable closed form: at n=0, f_0(x)=x*j_0(x)=sin(x) (the
!>      a_n branch) and f_0(x)=x*y_0(x)=-cos(x) (the b_n branch) --
!>      independent of SPHERICAL_BESSEL's own already-tested recursion,
!>      this checks the mode-accumulation/dispatch wiring itself;
!>  (3) the composability mechanism itself -- calling
!>      BUILD_INITIAL_CONDITION with the SAME tag listed twice produces
!>      exactly double the single-tag result, confirming each active tag
!>      entry accumulates into the same PHI/PSI rather than the last one
!>      overwriting. Only one real category exists so far (see the
!>      module's own header), so this is the closest available check on
!>      the list-of-tags mechanism until a second category exists.
PROGRAM TEST_INITIAL_CONDITIONS
USE KINDS,              ONLY: dp, i4
USE VSH,                ONLY: YLM_INDEX
USE GRID_RADIAL,        ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE FIELD_TYPES,        ONLY: SPECTRAL_SCALAR_T
USE INITIAL_CONDITIONS, ONLY: IC_MODE_T, READ_IC_NAMELIST, BUILD_INITIAL_CONDITION
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: SCRATCH_PATH = 'test_initial_conditions_scratch.nml'
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-10_dp
INTEGER(KIND=i4), PARAMETER :: N_R = 10_i4, LMAX = 3_i4

CHARACTER(LEN=64), ALLOCATABLE :: IC_TAGS(:)
TYPE(IC_MODE_T), ALLOCATABLE   :: MODES(:)
TYPE(RADIAL_GRID_T)            :: RGRID
TYPE(SPECTRAL_SCALAR_T)        :: PHI, PSI, PHI2, PSI2
INTEGER(KIND=i4) :: N_FAIL, UNIT, IR, IDX_10, IDX_20
REAL(KIND=dp)    :: R, EXPECTED_PHI

N_FAIL = 0

! --- (1) NAMELIST round-trip: single-entry tag list, 2 modes ---
OPEN(NEWUNIT=UNIT, FILE=SCRATCH_PATH, STATUS='REPLACE', ACTION='WRITE')
WRITE(UNIT,'(A)') "&IC_SELECT"
WRITE(UNIT,'(A)') "  N_TAGS = 1"
WRITE(UNIT,'(A)') "  TAGS = 'bessel_riccati'"
WRITE(UNIT,'(A)') "/"
WRITE(UNIT,'(A)') "&BESSEL_RICCATI_MODES"
WRITE(UNIT,'(A)') "  N_MODES      = 2"
WRITE(UNIT,'(A)') "  MODE_FIELD   = 'phi', 'psi'"
WRITE(UNIT,'(A)') "  MODE_L       = 1, 2"
WRITE(UNIT,'(A)') "  MODE_M       = 0, 1"
WRITE(UNIT,'(A)') "  MODE_N       = 0, 1"
WRITE(UNIT,'(A)') "  MODE_X_SCALE = 1.0, 2.5"
WRITE(UNIT,'(A)') "  MODE_A       = 1.0, 0.3"
WRITE(UNIT,'(A)') "  MODE_B       = 0.0, 0.7"
WRITE(UNIT,'(A)') "/"
CLOSE(UNIT)

CALL READ_IC_NAMELIST(SCRATCH_PATH, IC_TAGS, MODES)

CALL CHECK_TRUE("n_tags_parsed", SIZE(IC_TAGS) == 1, N_FAIL)
IF (SIZE(IC_TAGS) == 1) CALL CHECK_TRUE("tag1", TRIM(IC_TAGS(1)) == 'bessel_riccati', N_FAIL)
CALL CHECK_TRUE("n_modes_parsed", SIZE(MODES) == 2, N_FAIL)
IF (SIZE(MODES) == 2) THEN
  CALL CHECK_TRUE("mode1_field", TRIM(MODES(1)%FIELD) == 'phi', N_FAIL)
  CALL CHECK_TRUE("mode1_l",     MODES(1)%L == 1, N_FAIL)
  CALL CHECK_TRUE("mode1_m",     MODES(1)%M == 0, N_FAIL)
  CALL CHECK_TRUE("mode1_n",     MODES(1)%N == 0, N_FAIL)
  CALL CHECK_CLOSE("mode1_x_scale", MODES(1)%X_SCALE, 1.0_dp, TOL, N_FAIL)
  CALL CHECK_CLOSE("mode1_a",       MODES(1)%A, 1.0_dp, TOL, N_FAIL)
  CALL CHECK_CLOSE("mode1_b",       MODES(1)%B, 0.0_dp, TOL, N_FAIL)
  CALL CHECK_TRUE("mode2_field", TRIM(MODES(2)%FIELD) == 'psi', N_FAIL)
  CALL CHECK_TRUE("mode2_l",     MODES(2)%L == 2, N_FAIL)
  CALL CHECK_TRUE("mode2_m",     MODES(2)%M == 1, N_FAIL)
  CALL CHECK_TRUE("mode2_n",     MODES(2)%N == 1, N_FAIL)
  CALL CHECK_CLOSE("mode2_x_scale", MODES(2)%X_SCALE, 2.5_dp, TOL, N_FAIL)
  CALL CHECK_CLOSE("mode2_a",       MODES(2)%A, 0.3_dp, TOL, N_FAIL)
  CALL CHECK_CLOSE("mode2_b",       MODES(2)%B, 0.7_dp, TOL, N_FAIL)
END IF

OPEN(NEWUNIT=UNIT, FILE=SCRATCH_PATH, STATUS='OLD')
CLOSE(UNIT, STATUS='DELETE')

! --- (2) BUILD_INITIAL_CONDITION against hand-computable closed forms ---
! Single-shell grid, r in [1,3] km -- x=X_SCALE*r stays well away from 0.
CALL BUILD_RADIAL_GRID(RGRID, N_R, 1.0_dp, 3.0_dp, FULL_SPHERE=.FALSE.)

DEALLOCATE(MODES)
ALLOCATE(MODES(2))
! Mode 1: phi, l=1,m=0,n=0, x_scale=1, a=1,b=0 -> f_0(x)=x*j_0(x)=sin(x)=sin(r)
MODES(1) = IC_MODE_T(FIELD='phi', L=1_i4, M=0_i4, N=0_i4, X_SCALE=1.0_dp, A=1.0_dp, B=0.0_dp)
! Mode 2: phi, l=1,m=0,n=0, x_scale=1, a=0,b=1 -> f_0(x)=x*y_0(x)=-cos(x)=-cos(r)
! (same (l,m) as mode 1 -- tests that modes at the same (field,l,m) accumulate)
MODES(2) = IC_MODE_T(FIELD='phi', L=1_i4, M=0_i4, N=0_i4, X_SCALE=1.0_dp, A=0.0_dp, B=1.0_dp)

CALL BUILD_INITIAL_CONDITION([CHARACTER(LEN=64) :: 'bessel_riccati'], MODES, RGRID, LMAX, PHI, PSI)

IDX_10 = YLM_INDEX(1_i4, 0_i4)
DO IR = 1, N_R
  R = RGRID%R(IR)
  EXPECTED_PHI = SIN(R) - COS(R)   ! sum of both modes' closed forms
  CALL CHECK_CLOSE("phi_l1m0_accumulated", REAL(PHI%COEF(IR,IDX_10),KIND=dp), EXPECTED_PHI, TOL, N_FAIL)
END DO

! Untouched (l,m) slots and all of Psi should still be exactly zero.
IDX_20 = YLM_INDEX(2_i4, 0_i4)
CALL CHECK_TRUE("phi_l2m0_untouched_is_zero", ALL(PHI%COEF(:,IDX_20) == (0.0_dp,0.0_dp)), N_FAIL)
CALL CHECK_TRUE("psi_all_zero", ALL(PSI%COEF == (0.0_dp,0.0_dp)), N_FAIL)

! --- (3) Composability: the SAME tag listed twice must accumulate,
! not overwrite -- the mechanism a future second category will rely on.
CALL BUILD_INITIAL_CONDITION([CHARACTER(LEN=64) :: 'bessel_riccati', 'bessel_riccati'], &
  MODES, RGRID, LMAX, PHI2, PSI2)
DO IR = 1, N_R
  CALL CHECK_CLOSE("phi_l1m0_double_tag_accumulates", &
    REAL(PHI2%COEF(IR,IDX_10),KIND=dp), 2.0_dp*REAL(PHI%COEF(IR,IDX_10),KIND=dp), TOL, N_FAIL)
END DO

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_initial_conditions"
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

END PROGRAM TEST_INITIAL_CONDITIONS
