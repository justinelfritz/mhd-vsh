!> Round-trip test for IO_CHECKPOINT: writes a hand-picked (Phi,Psi)
!> state plus step/t/grid metadata, reads it back, and checks every
!> field matches to full double-precision round-trip precision (the
!> checkpoint format uses ES24.16 -- ~16 significant digits, so the
!> comparison tolerance is tight but not exactly machine epsilon, since
!> decimal<->binary conversion isn't perfectly lossless at any finite
!> number of digits).
PROGRAM TEST_IO_CHECKPOINT
USE KINDS,          ONLY: dp, i4
USE VSH,            ONLY: YLM_INDEX
USE GRID_RADIAL,    ONLY: RADIAL_GRID_T, BUILD_RADIAL_GRID
USE FIELD_TYPES,    ONLY: SPECTRAL_SCALAR_T, ALLOC_SPECTRAL_SCALAR
USE IO_CHECKPOINT,  ONLY: WRITE_CHECKPOINT, READ_CHECKPOINT
IMPLICIT NONE

CHARACTER(LEN=*), PARAMETER :: SCRATCH_PATH = 'test_io_checkpoint_scratch.dat'
INTEGER(KIND=i4), PARAMETER :: N_R   = 12
REAL(KIND=dp),    PARAMETER :: R_MIN = 10.8_dp, R_MAX = 11.7_dp
INTEGER(KIND=i4), PARAMETER :: LMAX  = 3
INTEGER(KIND=i4), PARAMETER :: ISTEP_WRITTEN = 12345
REAL(KIND=dp),    PARAMETER :: T_WRITTEN = 67.891_dp
REAL(KIND=dp),    PARAMETER :: TOL = 1.0E-14_dp

TYPE(RADIAL_GRID_T)     :: RGRID
TYPE(SPECTRAL_SCALAR_T) :: PHI, PSI, PHI2, PSI2
INTEGER(KIND=i4) :: IR, L, M, IDX, N_FAIL
INTEGER(KIND=i4) :: ISTEP_READ, N_R_READ, LMAX_READ
REAL(KIND=dp)    :: T_READ, R_MIN_READ, R_MAX_READ

CALL BUILD_RADIAL_GRID(RGRID, N_R, R_MIN, R_MAX, FULL_SPHERE=.FALSE.)
CALL ALLOC_SPECTRAL_SCALAR(PHI, N_R, LMAX)
CALL ALLOC_SPECTRAL_SCALAR(PSI, N_R, LMAX)

! Fill every mode with a distinct, non-trivial complex value per row so
! a mis-indexed round-trip (wrong (ir,idx) pairing) would be caught, not
! masked by a repeated/symmetric pattern.
DO L = 0, LMAX
  DO M = -L, L
    IDX = YLM_INDEX(L, M)
    DO IR = 1, N_R
      PHI%COEF(IR,IDX) = CMPLX(0.1_dp*IR + REAL(IDX,KIND=dp), -0.2_dp*IR + 0.5_dp*REAL(IDX,KIND=dp), KIND=dp)
      PSI%COEF(IR,IDX) = CMPLX(1.3_dp*IR - REAL(IDX,KIND=dp), 0.7_dp*IR + 0.1_dp*REAL(IDX,KIND=dp), KIND=dp)
    END DO
  END DO
END DO

CALL WRITE_CHECKPOINT(SCRATCH_PATH, ISTEP_WRITTEN, T_WRITTEN, PHI, PSI, RGRID)
CALL READ_CHECKPOINT(SCRATCH_PATH, ISTEP_READ, T_READ, PHI2, PSI2, N_R_READ, LMAX_READ, R_MIN_READ, R_MAX_READ)

N_FAIL = 0
CALL CHECK_INT("n_r",   N_R_READ,   N_R,   N_FAIL)
CALL CHECK_INT("lmax",  LMAX_READ,  LMAX,  N_FAIL)
CALL CHECK_INT("istep", ISTEP_READ, ISTEP_WRITTEN, N_FAIL)
CALL CHECK_REAL("r_min", R_MIN_READ, R_MIN, TOL, N_FAIL)
CALL CHECK_REAL("r_max", R_MAX_READ, R_MAX, TOL, N_FAIL)
CALL CHECK_REAL("t",     T_READ,     T_WRITTEN, TOL, N_FAIL)
CALL CHECK_REAL("phi_coef_max_err", MAXVAL(ABS(PHI2%COEF-PHI%COEF)), 0.0_dp, TOL, N_FAIL)
CALL CHECK_REAL("psi_coef_max_err", MAXVAL(ABS(PSI2%COEF-PSI%COEF)), 0.0_dp, TOL, N_FAIL)

OPEN(UNIT=52, FILE=SCRATCH_PATH, STATUS='OLD')
CLOSE(UNIT=52, STATUS='DELETE')

IF (N_FAIL > 0) THEN
  WRITE(*,'(A,I0,A)') "RESULT: FAILED - ", N_FAIL, " check(s) failed"
  STOP 1
ELSE
  WRITE(*,'(A)') "RESULT: PASSED - test_io_checkpoint"
END IF

CONTAINS

SUBROUTINE CHECK_INT(NAME, GOT, EXPECTED, N_FAIL)
  CHARACTER(*),     INTENT(IN)    :: NAME
  INTEGER(KIND=i4), INTENT(IN)    :: GOT, EXPECTED
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  IF (GOT /= EXPECTED) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,I0,A,I0)') "FAIL  ", NAME, "  got=", GOT, "  expected=", EXPECTED
  ELSE
    WRITE(*,'(A,A,A,I0)') "PASS  ", NAME, "  value=", GOT
  END IF
END SUBROUTINE CHECK_INT

SUBROUTINE CHECK_REAL(NAME, GOT, EXPECTED, TOL_ARG, N_FAIL)
  CHARACTER(*),     INTENT(IN)    :: NAME
  REAL(KIND=dp),    INTENT(IN)    :: GOT, EXPECTED, TOL_ARG
  INTEGER(KIND=i4), INTENT(INOUT) :: N_FAIL
  REAL(KIND=dp) :: ERR
  ERR = ABS(GOT-EXPECTED)
  IF (ERR > TOL_ARG) THEN
    N_FAIL = N_FAIL + 1
    WRITE(*,'(A,A,A,ES10.3)') "FAIL  ", NAME, "  err=", ERR
  ELSE
    WRITE(*,'(A,A,A,ES10.3)') "PASS  ", NAME, "  err=", ERR
  END IF
END SUBROUTINE CHECK_REAL

END PROGRAM TEST_IO_CHECKPOINT
